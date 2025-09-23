'''
Given a query structure (a crystal structure of a protein-ligand complex) and its 
binding site residues (vdms), determine whether any vdGs can correctly place the chemical 
groups within the ligand of the query structure.

Usage:
    >> python rerank_poses.py $YOUR_YAML_FILE
'''

import sys
import os
import yaml
import numpy as np
import prody as pr
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import match_vdgs as match
import dock_utils as dock
import Frags
import utils
import concurrent.futures  # for ProcessPoolExecutor
import signal              # to catch SIGINT/SIGTERM
from functools import partial
import io
from contextlib import redirect_stdout, redirect_stderr

with open(sys.argv[1], 'r') as f:
    config = yaml.safe_load(f)

# Config variables
rmsd_threshold = config['rmsd_threshold']
query_pdbs_dir = config['query_pdbs_dir']
solved_struct_path = config['solved_struct_path']
lig_smiles = config['lig_smiles']
vdg_lib_dir = config['vdg_lib_dir']
outdir = config['outdir']
frags_to_exclude = config['frags_to_exclude']  # list of SM
frags_to_include = config['frags_to_include']  # list of SM
query_pdbs = config.get('query_pdbs')  # value should be either 'all' 
                                       # (i.e. all pdbs in query_pdbs_dir) or a list
overwrite_existing = config.get('overwrite_existing', False)

print('Config variables:')
print(f"rmsd_threshold: {rmsd_threshold}")
print(f"query_pdbs_dir: {query_pdbs_dir}")
print(f"solved_struct_path: {solved_struct_path}")
print(f"lig_smiles: {lig_smiles}")
print(f"vdg_lib_dir: {vdg_lib_dir}")
print(f"outdir: {outdir}")
print(f"frags_to_exclude: {frags_to_exclude}")
print(f"frags_to_include: {frags_to_include}")
print(f"query_pdbs: {query_pdbs}")
print(f"overwrite_existing: {overwrite_existing}")

MAX_WORKERS = 10

# Returns a tuple: (pdbfile, list_of_failed_vdg_files)
def process_one_pdb(pdbfile, solved_struct_path, rmsd_threshold, query_pdbs_dir, 
    lig_smiles, vdg_lib_dir, outdir, frags_to_exclude, frags_to_include,
    overwrite_existing):

    _stdout_buf = io.StringIO()
    _stderr_buf = io.StringIO()
    with redirect_stdout(_stdout_buf), redirect_stderr(_stderr_buf):

        # Load ground truth structure
        solved_struct = pr.parsePDB(solved_struct_path)

        # Get binding site residues from the solved structure
        lig_obj = solved_struct.hetatm.select(
            'not (ion or resname SEP or resname TPO or resname MSE)') 
        # if lig_obj is None, then raise an error and note to the user that they may 
        # have to reformat the resname in the PDB file
        if lig_obj is None:
            raise ValueError(
                f'Prody could not find ligand atoms in {solved_struct_path}. You may need to '
                f'reformat the lig resname in the PDB file.')
        ligname = set(lig_obj.getResnames())
        assert len(ligname) == 1
        ligname = list(ligname)[0]

        # Gather all binding site residue combinations
        all_bsr_combos = dock.get_bsr_combinations(solved_struct, ligname)
        # For each prediction, get the frags and CG coords in the pred, iterate through 
        # binding site res combos, align, and determine if there are vdg matches. 
    
        # Keep track of which frags are in vdg lib and which aren't
        frags_in_lib = {}

        # Keep track of which database vdg files could not be loaded
        failed_vdg_files = []

        # Load predictions
        print('\n' + '='*10, pdbfile, '='*10)
        output_dir = dock.name_outdir(pdbfile, outdir, 2) # hardcoding to 2 is intentional
        utils.set_up_outdir(output_dir, overwrite=overwrite_existing)

        # Get all fragments. For every frag, iterate over the bsr combos and determine 
        # vdg matches.
        orig_filtered_frags, pdb_mol_reassigned = Frags.get_frags_from_pdbfile(
            os.path.join(query_pdbs_dir, pdbfile), lig_smiles)
        # Remove duplicates of filtered_frags 
        deduplicated_filtered_frags = []
        for (sub, sub_smiles) in orig_filtered_frags:
            # Step 0) is this frag named something else in the vdg lib?
            for existingfrag in os.listdir(vdg_lib_dir):
                if utils.smiles_equiv(existingfrag, sub_smiles):
                    sub_smiles = existingfrag
                    break

            # First of all, is it in frags_to_exclude?
            if Frags.check_in_exclude_list(sub_smiles, frags_to_exclude): 
                # means it's equiv to a frag in frags_to_exclude
                continue

            # If sub_smiles already exists, then skip. Otherwise, add
            if sub_smiles in [i[1] for i in deduplicated_filtered_frags]:
                continue
            if sub_smiles not in [i[1] for i in deduplicated_filtered_frags]:
                deduplicated_filtered_frags.append((sub, sub_smiles))
            if sub_smiles in frags_to_exclude: # Exclude specified frags 
                continue
    
            # Determine if this frag is in vdg lib (and the job didn't break)
            if sub_smiles not in frags_in_lib.keys():
                # Is it in the vdg lib dir?
                if sub_smiles in os.listdir(vdg_lib_dir):
                    # Did the job (to create the vdgs) finish w/o issue?
                    if Frags.check_vdg_job_status(sub_smiles, vdg_lib_dir):
                        frags_in_lib[sub_smiles] = True
                    else:
                        frags_in_lib[sub_smiles] = False 
                else: 
                    # If the str isn't in os.listdir(), it's still possible that a 
                    # degenerate smiles is in vdg_lib_dir
                    for existingfrag in os.listdir(vdg_lib_dir):
                        # Is it equivalent to sub_smiles and was its job completed?
                        is_equivalent = utils.smiles_equiv(existingfrag, sub_smiles)
                        if is_equivalent: # rename the smiles to match what's in db
                            sub_smiles = existingfrag
                        if is_equivalent and Frags.check_vdg_job_status(existingfrag, vdg_lib_dir):
                            frags_in_lib[sub_smiles] = True
                            break
                    else:
                        frags_in_lib[sub_smiles] = False
        # Log
        Frags.summarize_frags(deduplicated_filtered_frags, frags_in_lib, frags_to_exclude, 
                              frags_to_include)

        # Iterate over frags that are in the vdg lib
        for sub, sub_smiles in deduplicated_filtered_frags:
            if frags_to_include == 'all':
                pass
            elif isinstance(frags_to_include, list):
                # If frags_to_include is a list, check if sub_smiles is in the list
                if sub_smiles not in frags_to_include:
                    continue
            else:
                raise ValueError("frags_to_include must be 'all' or a list of SMILES.")
            if sub_smiles not in frags_in_lib.keys():
                # Meaning: not being searched for (possibly being excluded in yml file)
                continue
            if not frags_in_lib[sub_smiles]:
                # Meaning: not in vdg lib
                continue
            print('Searching for', sub_smiles)
            # Get the query CG coords for this substructure in the query struct. Should 
            # return all instances and permutations.
            # get the atom indices for the sub molecule

            q_cg_coord_perms = dock.get_query_cg_coords(pdb_mol_reassigned, sub)
            # Iterate over the query structure's binding site residues and superpose vdGs 
            # (incl. the CG) onto it, keeping track of the vdGs that have low rmsd to the 
            # known bb+CG.
            for combo, coords in all_bsr_combos:
                # Retain PAIRS ONLY
                if len(combo) != 2:
                    continue
                bsr_incl_bb_identities = combo
                input_bsr_bb_coords = coords
                # Reorder the bb coords and bsr_AA_identities in alphabetic order
                paired = list(zip(bsr_incl_bb_identities, input_bsr_bb_coords))
                paired_sorted = sorted(paired, key=lambda pair: pair[0])
                # Unzip them back
                bsr_incl_bb_identities, input_bsr_bb_coords = zip(*paired_sorted)
                bsr_incl_bb_identities = '_'.join(sorted(bsr_incl_bb_identities))
                input_bsr_bb_coords = np.array(input_bsr_bb_coords)
                # Flatten input_bsr_bb_coords
                input_bsr_bb_coords = input_bsr_bb_coords.reshape(-1, 3)
                # Dock fragments by aligning vdg backbones to bsr backbones
                vdgs_path = os.path.join(vdg_lib_dir, sub_smiles, 'nr_vdgs', str(
                    len(bsr_incl_bb_identities.split('_'))), bsr_incl_bb_identities)
                if not os.path.exists(vdgs_path):
                    continue
                vdg_paths = [p for p in os.listdir(vdgs_path) if p.endswith(('.pdb', '.pdb.gz'))]
                # Iterate over each vdg
                for vdg_path in vdg_paths:
                    vdg_full_path = os.path.join(vdgs_path, vdg_path)
                    try:
                        vdg_prody_obj = pr.parsePDB(vdg_full_path)
                    except Exception as e:
                        failed_vdg_files.append(f"{vdg_full_path} :: {e.__class__.__name__}: {e}")
                        continue
                    # Superpose the vdg backbone atoms onto the input structure.
                    # Make sure order of query vdg AAs is consistent with order of input_pdb AAs
                    resind_perms = dock.map_aa_identities_to_vdg_resinds(
                        [vdg_prody_obj.select(f'resindex {_r}').getResnames()[0] 
                            for _r in sorted(set(vdg_prody_obj.getResindices()))], 
                        bsr_incl_bb_identities.split('_'), )

                    # Iterate over each resind permutation, get bb coords, and superpose onto 
                    # the query structure (input structure). Multiple resind perms or cg atom 
                    # perms could have rmsd < threshold, so keep track of minimum rmsd 
                    best_rmsd, best_resind_perm, best_transf = [None, None, None]
                    for resind_perm in resind_perms:
                        vdg_perm_bb_coords = []
                        for _r in resind_perm:
                            for atom_name in ['N', 'CA', 'C']: 
                                _coords = utils.get_atom_coords(vdg_prody_obj.select(
                                    f'resindex {_r}'), atom_name)
                                vdg_perm_bb_coords.append(_coords)
                        vdg_perm_bb_coords = np.array(vdg_perm_bb_coords)
                        # Superpose the database vdg bb + cg atoms onto the input structure 
                        # and calculate rmsd.
                        # --- First, get combined bb+cg coords for db vdg
                        vdg_db_cg_coords = np.array(match.get_database_cg_coords(vdg_prody_obj))
                        try:
                            db_bb_and_cg_coords = np.concatenate((vdg_perm_bb_coords, vdg_db_cg_coords), 
                                                             axis=0)
                        except:
                            continue
                        # --- Then, superpose and calc rmsd for all query CG permutations
                        for q_cg_coord_perm in q_cg_coord_perms:
                            q_bb_and_cg = np.concatenate((input_bsr_bb_coords, q_cg_coord_perm), 
                                                         axis=0)
                            if db_bb_and_cg_coords.shape != q_bb_and_cg.shape:
                                raise ValueError(
                                    'The query CG must have the same num. of atoms as the database CG.')
                            moved_db_vdg_bb_and_cg_coords, db_transf = pr.superpose(
                                db_bb_and_cg_coords, q_bb_and_cg)
                            db_to_query_rmsd = pr.calcRMSD(moved_db_vdg_bb_and_cg_coords, 
                                                           q_bb_and_cg)
                            if db_to_query_rmsd > rmsd_threshold:
                                continue
                            # If this is the best rmsd so far, keep track of it
                            if best_rmsd is None or db_to_query_rmsd < best_rmsd:
                                best_rmsd = db_to_query_rmsd
                                best_resind_perm = resind_perm
                                #best_q_cg_coord_perm = q_cg_coord_perm
                                best_transf = db_transf
                    # If we found a best rmsd, write out the vdg with the best rmsd
                    if best_rmsd is not None:
                        resind_perm = best_resind_perm
                        #q_cg_coord_perm = best_q_cg_coord_perm
                        db_transf = best_transf
                        db_to_query_rmsd = best_rmsd
                        # Determine output name
                        bsr_perm_str = ''
                        for _resix in resind_perm: 
                            resix_sel = vdg_prody_obj.select(f'resindex {_resix}')
                            resname = resix_sel.getResnames()[0]
                            reschain = resix_sel.getChids()[0]
                            resnum = resix_sel.getResnums()[0]
                            resseg = resix_sel.getSegnames()[0]
                            bsr_perm_str += f'_{resname}_{resseg}_{reschain}_{resnum}'
                        # Output vdg name 
                        out_subdir = f'{sub_smiles}{bsr_perm_str}'
                        subdir_path = os.path.join(output_dir, sub_smiles, out_subdir)
                        os.makedirs(subdir_path, exist_ok=True)
                        output_vdg_name = (
                            f"{pdbfile.split('/')[-1].removesuffix('.pdb')}_"
                            f"{sub_smiles}_"
                            f"{vdg_path.rstrip('/').split('/')[-1].removesuffix('.pdb.gz')}_"
                            f"{db_to_query_rmsd:.2f}")

                        output_vdg_path = os.path.join(subdir_path, output_vdg_name+'.pdb')
                        # check if a file in output_vdg_path already exists
                        if os.path.exists(output_vdg_path):
                            print('WARNING: overwriting existing file', output_vdg_path)
                        pr.writePDB(output_vdg_path, pr.applyTransformation(db_transf, 
                            vdg_prody_obj.copy()))

        return (pdbfile, failed_vdg_files, _stdout_buf.getvalue(), _stderr_buf.getvalue())


def main():
    # Create a new process group 
    try:
        os.setpgrp()
    except Exception:
        pass

    def _terminate_all(signum, frame): # for graceful shutdown of all workers
        try:
            # Send SIGTERM to the whole process group
            os.killpg(0, signal.SIGTERM)
        except Exception:
            pass
        # Exit the parent
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGINT, _terminate_all)
    signal.signal(signal.SIGTERM, _terminate_all)

    # Build the filtered pdb list
    all_pdbs = sorted(os.listdir(query_pdbs_dir))
    if query_pdbs == 'all':
        pdb_list = all_pdbs
    elif isinstance(query_pdbs, list):
        pdb_list = [p for p in all_pdbs if p in query_pdbs]
    else:
        raise ValueError("query_pdbs must be 'all' or a list of pdb files.")

    if not pdb_list:
        print("No PDBs to process.")
        return

    # Run in parallel
    failed_vdg_files_global = []

    worker_fn = partial(process_one_pdb,
        solved_struct_path=solved_struct_path, rmsd_threshold=rmsd_threshold,
        query_pdbs_dir=query_pdbs_dir, lig_smiles=lig_smiles, vdg_lib_dir=vdg_lib_dir,
        outdir=outdir, frags_to_exclude=frags_to_exclude, 
        frags_to_include=frags_to_include, overwrite_existing=overwrite_existing)

    max_workers = (MAX_WORKERS if isinstance(MAX_WORKERS, int) and MAX_WORKERS > 0 
                   else os.cpu_count())

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(worker_fn, pdbfile): pdbfile for pdbfile in pdb_list}
        try:
            for fut in concurrent.futures.as_completed(futures):
                pdbfile = futures[fut]


                try:
                    # upack captured logs
                    pdbfile, failed, captured_out, captured_err = fut.result()
                    failed_vdg_files_global.extend(failed)

                    # print this PDB's logs as one coherent block
                    print("\n" + "#"*80)
                    if captured_out:
                        sys.stdout.write(captured_out)
                        if not captured_out.endswith("\n"):
                            sys.stdout.write("\n")
                    if captured_err:
                        print(f"[stderr] ---------------- {pdbfile} ----------------")
                        sys.stdout.write(captured_err)
                        if not captured_err.endswith("\n"):
                            sys.stdout.write("\n")
                    print("#"*80 + "\n")

                except Exception as e:
                    print(f"ERROR: {pdbfile} failed: {e.__class__.__name__}: {e}")
                    continue

        except KeyboardInterrupt:
            # If Ctrl-C occurs while waiting, cancel pending futures and re-raise
            for f in futures:
                f.cancel()
            raise
        except SystemExit:
            # If we called _terminate_all from a signal, propagate
            for f in futures:
                f.cancel()
            raise

    # Log
    if failed_vdg_files_global:
        print('\nvdG files that could not be parsed:')
        for failed_file in failed_vdg_files_global:
            print(f" - {failed_file}")

    return

if __name__ == "__main__":
    main()
