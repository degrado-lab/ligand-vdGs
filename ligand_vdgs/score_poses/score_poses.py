'''
Count the number of vdGs that match the interactions between all ligand CGs and 
the residues they interact with in a query structure. 

Usage:
    >> python score_poses.py $YOUR_YAML_FILE

'''

import sys
import os
import yaml
import numpy as np
import prody as pr
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import match_vdgs as match
import dock_utils as dock
import Frags
import utils
import multiprocessing
import signal
import platform

def _init_worker():
    # Makes it so that subprocesses will die when the main process dies
    signal.signal(signal.SIGINT, signal.SIG_IGN)  # ignore SIGINT in workers
    if platform.system() == "Linux":
        try:
            import ctypes
            libc = ctypes.CDLL("libc.so.6")
            PR_SET_PDEATHSIG = 1
            libc.prctl(PR_SET_PDEATHSIG, signal.SIGTERM)
        except Exception:
            pass

def _process_one_sub_smiles(pdbfile, sub, sub_smiles, vdg_lib_dir, output_dir,
                            all_bsr_combos, rmsd_threshold, lig_smiles):

    failed_vdg_files_local = []

    try:
        # Load struct-specific state we need to compute q_cg_coord_perms for this sub
        # (rebuild the reassigned molecule here inside the process)
        # Get all fragments. For every frag, iterate over the bsr combos and determine vdg matches.
        # We only need pdb_mol_reassigned for q_cg_coord_perms.
        _, pdb_mol_reassigned = Frags.get_frags_from_pdbfile(pdbfile, lig_smiles)

        print('Searching for', sub_smiles, flush=True)
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
                except Exception:
                    failed_vdg_files_local.append(vdg_full_path)
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
                        print('-- WARNING: overwriting existing file', output_vdg_path, flush=True)
                    pr.writePDB(output_vdg_path, pr.applyTransformation(db_transf, 
                        vdg_prody_obj.copy()))
    except Exception as e:
        print(f"[worker:{sub_smiles}] ERROR: {e}", flush=True)

    return failed_vdg_files_local
# -----------------------------------------------------------------------------
def main():
    with open(sys.argv[1], 'r') as f:
        config = yaml.safe_load(f)

    # Config variables
    rmsd_threshold = config['rmsd_threshold']
    query_pdb_paths = config['query_pdb_paths']
    lig_smiles = config['lig_smiles']
    vdg_lib_dir = config['vdg_lib_dir']
    outdir = config['outdir']
    frags_to_exclude = config['frags_to_exclude']  # list of SM
    frags_to_include = config['frags_to_include']  # list of SM
    overwrite_existing = config.get('overwrite_existing', False)
    n_procs = int(config.get('n_procs', os.cpu_count() or 1))
    
    print('Config variables:')
    print(f"rmsd_threshold: {rmsd_threshold}")
    print(f"query_pdb_paths: {query_pdb_paths}")
    print(f"lig_smiles: {lig_smiles}")
    print(f"vdg_lib_dir: {vdg_lib_dir}")
    print(f"outdir: {outdir}")
    print(f"frags_to_exclude: {frags_to_exclude}")
    print(f"frags_to_include: {frags_to_include}")
    print(f"overwrite_existing: {overwrite_existing}")
    print(f"n_procs: {n_procs}")
    
    '''
    For each query struct, get the frags and CG coords in the struct, iterate through 
    binding site res combos, align, and determine if there are vdg matches. 
    ''' 
    
    # Keep track of which frags are in vdg lib and which aren't
    frags_in_lib = {}

    # Keep track of which database vdg files could not be loaded
    failed_vdg_files = []

    # Load query structs
    for pdbfile in query_pdb_paths:
        print('\n' + '='*10, pdbfile, '='*10)
        output_dir = dock.name_outdir(pdbfile, outdir, len(query_pdb_paths))
        if len(query_pdb_paths) == 1:
            output_dir += '_frags'
        utils.set_up_outdir(output_dir, overwrite=overwrite_existing)

        # Load struct and get all binding site residue combinations
        query_obj = pr.parsePDB(pdbfile)
        ligname = set(query_obj.hetatm.select('not (ion or resname SEP or resname TPO or resname MSE)').getResnames())
        assert len(ligname) == 1
        ligname = list(ligname)[0]
        all_bsr_combos = dock.get_bsr_combinations(query_obj, ligname)

        # Precompute lib listing once
        vdg_dir_listing = set(os.listdir(vdg_lib_dir))

        def canonicalize_smiles(smiles: str) -> str:
            "Return lib's canonical name for this SMILES if an equivalent exists; else original."
            if smiles in vdg_dir_listing:
                return smiles
            for existing in vdg_dir_listing:
                if utils.smiles_equiv(existing, smiles):
                    return existing
            return smiles

        # Get all fragments. For every frag, iterate over the bsr combos and determine vdg matches.
        orig_filtered_frags, pdb_mol_reassigned = Frags.get_frags_from_pdbfile(pdbfile, lig_smiles)

        # Remove duplicates of filtered_frags (post-exclude, using canonical names)
        deduplicated_filtered_frags = []
        seen = set()
        for (sub, sub_smiles) in orig_filtered_frags:
            # Is this frag named something else in the vdg lib?
            # Map to library naming if an equivalent SMILES exists in the lib (canonicalize once)
            canon = canonicalize_smiles(sub_smiles)

            # First of all, is it in frags_to_exclude?
            if Frags.check_in_exclude_list(canon, frags_to_exclude):
                continue
            
            # Dedup by canonical name
            if canon in seen:
                continue
            seen.add(canon)
            deduplicated_filtered_frags.append((sub, canon))

            # Mark presence & job completion in lib (fill frags_in_lib once)
            if canon not in frags_in_lib:
                in_dir = canon in vdg_dir_listing
                frags_in_lib[canon] = in_dir and Frags.check_vdg_job_status(canon, vdg_lib_dir)

        # Log
        Frags.summarize_frags(deduplicated_filtered_frags, frags_in_lib, frags_to_exclude, 
                              frags_to_include)

        # ---------------------------------------------------------------------
        # Concurrently iterate over frags that are in the vdg lib
        
        futures = []
        # Ensure the executor definitely uses the 'spawn' context
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=n_procs, initializer=_init_worker, mp_context=ctx) as executor:
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

                # Submit this fragment for concurrent processing
                futures.append(executor.submit(_process_one_sub_smiles,
                    pdbfile, sub, sub_smiles, vdg_lib_dir, output_dir,
                    all_bsr_combos, rmsd_threshold, lig_smiles))

            # Gather results (e.g., any vdG files that failed to parse)
            for fut in as_completed(futures):
                try:
                    failed_vdg_files.extend(fut.result() or [])
                except Exception as e:
                    print(f"[main] A worker failed: {e}", flush=True)
        # ---------------------------------------------------------------------

    if failed_vdg_files:
        print('\nvdG files that could not be parsed:')
        for failed_file in failed_vdg_files:
            print(f" - {failed_file}")

    print('Job completed.')
    return


if __name__ == "__main__":
    multiprocessing.set_start_method("spawn", force=True) 
    main()