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
import multiprocessing as mp
from contextlib import redirect_stdout, redirect_stderr
import logging
from rdkit import rdBase
from rdkit import RDLogger
import io

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

_VDG_DIR_CACHE = {}             # path -> listdir(list[str])
_VDG_PARSE_CACHE = {}           # path -> (vdg_prody_obj, vdg_AAs, 
                                # vdg_db_cg_coords(np.ndarray), 
                                # bb_cache{resindex:{'N':arr,'CA':arr,'C':arr}})
_SOLVED_STRUCT = None           # for caching across workers
_ALL_BSR_COMBOS = None          # for caching across workers

def _listdir_cached(path):
    # cache os.listdir to avoid repeated disk scans.
    lst = _VDG_DIR_CACHE.get(path)
    if lst is None:
        lst = os.listdir(path)
        _VDG_DIR_CACHE[path] = lst
    return lst

def _precompute_vdg(vdg_full_path):
    # Precompute resnames, CG coords, and per-residue backbone coords.
    if vdg_full_path in _VDG_PARSE_CACHE:
        return _VDG_PARSE_CACHE[vdg_full_path]
    vdg_prody_obj = pr.parsePDB(vdg_full_path)
    resindices = sorted(set(vdg_prody_obj.getResindices()))
    num_vdms = int(vdg_full_path.split('/')[-3])
    # Not sure why, but some vdg pdbs have duplicate residues. Resindex 0 is always 
    # the lig, and the 1st and 2nd resindices are the bsrs. Select just those indices 
    # and ignore the duplicated residues.
    if num_vdms == 2:
        resindices = resindices[:3] # first 3 residues are always lig + 2 bsrs.
        bb_resinds = [1, 2] 
    elif num_vdms == 1:
        resindices = resindices[:2]
        bb_resinds = [1]
    vdg_AAs = [vdg_prody_obj.select(f'resindex {_r}').getResnames()[0] for _r in resindices]
    # Precompute backbone atom coords per residue once
    bb_cache = {}
    for _r in bb_resinds: # resindex 0 = lig
        sel = vdg_prody_obj.select(f'resindex {_r}')
        bb_cache[_r] = {
            'N':  utils.get_atom_coords(sel, 'N'),
            'CA': utils.get_atom_coords(sel, 'CA'),
            'C':  utils.get_atom_coords(sel, 'C'),
        }
    vdg_db_cg_coords = np.array(match.get_database_cg_coords(vdg_prody_obj))
    _VDG_PARSE_CACHE[vdg_full_path] = (vdg_prody_obj, vdg_AAs, vdg_db_cg_coords, bb_cache)
    return _VDG_PARSE_CACHE[vdg_full_path]

def _init_rdkit_logging(stream_like):
    # Route RDKit messages to the python logger so warnings/errors are grouped with 
    # each pdbfile's print output in order, preventing interleaving of messages b/n 
    # parallel workers.
    rdBase.LogToPythonLogger()
    logger = logging.getLogger('rdkit')
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(stream_like)
    sh.setFormatter(logging.Formatter('%(levelname)s [%(name)s]: %(message)s'))
    logger.addHandler(sh)

def _process_one_pdb(pdbfile):
    # Run vdg matching. Capture stdout/stderr to an in-memory buffer so that messages 
    # won't be interleaved when using multiple processes.
    global _SOLVED_STRUCT, _ALL_BSR_COMBOS
    if _SOLVED_STRUCT is None or _ALL_BSR_COMBOS is None:
        _SOLVED_STRUCT = pr.parsePDB(solved_struct_path)
        lig_obj = _SOLVED_STRUCT.hetatm.select(
            'not (ion or resname SEP or resname TPO or resname MSE)')
        if lig_obj is None:
            raise ValueError(
                f'Prody could not find ligand atoms in {solved_struct_path}. You may need '
                f'to reformat the lig resname in the PDB file.')
        ligname = list(set(lig_obj.getResnames()))[0]
        _ALL_BSR_COMBOS = dock.get_bsr_combinations(_SOLVED_STRUCT, ligname)

    failed_vdg_files = []
    buf = io.StringIO()

    # Capture all output for this pdbfile into a buffer
    with redirect_stdout(buf), redirect_stderr(buf):
        _init_rdkit_logging(buf)
        output_dir = dock.name_outdir(pdbfile, outdir, make_pdb_subfolder=True)
        utils.set_up_outdir(output_dir, overwrite=overwrite_existing)

        # Get all fragments. For every frag, iterate over the bsr combos and determine 
        # vdg matches.
        frags_in_lib = {}
        dbname_smiles_map = {} # map original -> dbname name used in vdg lib
            
        orig_filtered_frags, pdb_mol_reassigned = Frags.get_frags_from_pdbfile(
            os.path.join(query_pdbs_dir, pdbfile), lig_smiles, quiet=True)
        filtered_frags = {}
        for sub_smiles, substruct_perm_grouped_by_site in orig_filtered_frags.items():
            frags_in_lib, db_name = check_frag_in_vdg_lib(sub_smiles, frags_in_lib, vdg_lib_dir)
            dbname_smiles_map[sub_smiles] = db_name
            if frags_in_lib.get(db_name, False):
                # store under db_name key so downstream lookups are consistent
                filtered_frags[db_name] = substruct_perm_grouped_by_site


        # For every frag that's in the vdg lib, iterate over the bsr combos and determine 
        # vdg matches.
        for sub_smiles, substruct_perm_grouped_by_site in filtered_frags.items():
            should_run = should_run_frag(sub_smiles, frags_to_include, frags_in_lib)
            # Determine frags to skip over
            if not should_run:
                print(f'Skipping {sub_smiles}..........................')
                continue
            # Otherwise, proceed with finding vdg matches for this frag
            print(f'Searching for {sub_smiles}..........................')

            # Get the query CG coords for this substructure in the query struct. 
            grouped_q_cg_perms = [] # all instances of CG + perms, grouped by site on lig
            for substruct_site in substruct_perm_grouped_by_site:
                q_cg_perms_single_site = []
                for substruct_perm in substruct_site:
                    sub, perm_inds, orig_mol_inds = substruct_perm
                    single_perm_cg_coords = dock.get_query_cg_coords(sub, sub_smiles)
                    q_cg_perms_single_site.append((single_perm_cg_coords, perm_inds, 
                                                   orig_mol_inds))
                grouped_q_cg_perms.append(q_cg_perms_single_site)

            # Iterate over the query structure's binding site residues and superpose vdGs 
            # (incl. the CG) onto it, keeping track of the vdGs that have low rmsd to the 
            # known bb+CG.
            for combo, bsrlabel, bsrAAs, coords in _ALL_BSR_COMBOS:
                bsr_incl_bb_identities = combo
                input_bsr_bb_coords = coords
                # Reorder the bb coords and bsr_AA_identities in alphabetic order
                paired = list(zip(bsr_incl_bb_identities, input_bsr_bb_coords))
                paired_sorted = sorted(paired, key=lambda pair: pair[0])
                # Unzip them back
                bsr_incl_bb_identities, input_bsr_bb_coords = zip(*paired_sorted)
                bsr_incl_bb_identities = '_'.join(sorted(bsr_incl_bb_identities))
                # Flatten input_bsr_bb_coords
                input_bsr_bb_coords = np.array(input_bsr_bb_coords).reshape(-1, 3)
                # Dock fragments by aligning vdg backbones to bsr backbones
                vdgs_path = os.path.join(vdg_lib_dir, sub_smiles, 'nr_vdgs',
                                         str(len(bsr_incl_bb_identities.split('_'))),
                                         bsr_incl_bb_identities)
                if not os.path.exists(vdgs_path):
                    continue

                # Iterate over each vdg. Each vdg represents a pair of bsrs.
                vdg_paths = [p for p in _listdir_cached(vdgs_path) if p.endswith(
                    ('.pdb', '.pdb.gz'))]
                for vdg_path in vdg_paths:
                    vdg_full_path = os.path.join(vdgs_path, vdg_path)
                    try: # parse and precompute once
                        vdg_prody_obj, vdg_AAs, vdg_db_cg_coords, bb_cache = \
                            _precompute_vdg(vdg_full_path)
                    except Exception:
                        failed_vdg_files.append(vdg_full_path)
                        continue

                    # Superpose the vdg backbone atoms onto the input structure.
                    # Make sure order of query vdg AAs is consistent with order of 
                    # input_pdb AAs. Use cached AA list & cached mapper.
                    resind_perms = dock.map_aa_identities_to_vdg_resinds(vdg_AAs,
                        bsr_incl_bb_identities.split('_'))

                    # Determine if this vdg has bb+cg geometry that matches the query 
                    # struct. Iterate over each CG site and each CG atom perm within that 
                    # site. Keep track of the lowest rmsd result, b/c there are cases where 
                    # >1 resind perms and/or cg perms yield rmsd < threshold. 
                    best_rmsd, best_resind_perm, best_cg_perm, best_transf = [None, None, 
                        None, None] # best_cg_perm will describe both site and atom perm
                    for resind_perm in resind_perms:
                        # build backbone coords from bb_cache to minimize prody overhead
                        vdg_perm_bb_coords = np.array([coord for _r in resind_perm for coord 
                            in (bb_cache[_r]['N'], bb_cache[_r]['CA'], bb_cache[_r]['C'])])
                        # Superpose the database vdg bb + cg atoms onto the input structure 
                        # and calculate rmsd.
                        # --- First, get combined bb+cg coords for db vdg
                        try: # precomputed cg coords
                            db_bb_and_cg_coords = np.concatenate((vdg_perm_bb_coords, 
                                                                  vdg_db_cg_coords), axis=0)
                        except Exception:
                            continue
                        # --- Then, superpose and calc rmsd for all query CG perms
                        for q_cg_site in grouped_q_cg_perms:
                            for q_cg_perm in q_cg_site:
                                q_cg_coord_perm, perm_inds, orig_mol_inds = q_cg_perm
                                q_bb_and_cg = np.concatenate((input_bsr_bb_coords, 
                                                              q_cg_coord_perm), axis=0)
                                if db_bb_and_cg_coords.shape != q_bb_and_cg.shape:
                                    raise ValueError('The query CG must have the same '
                                    'num. of atoms as the database CG.')
                                moved_db_vdg_bb_and_cg_coords, db_transf = pr.superpose(
                                    db_bb_and_cg_coords, q_bb_and_cg)
                                db_to_query_rmsd = pr.calcRMSD(
                                    moved_db_vdg_bb_and_cg_coords, q_bb_and_cg)
                                if db_to_query_rmsd > rmsd_threshold:
                                    continue
                                # Keep track of best rmsd so far
                                if best_rmsd is None or db_to_query_rmsd < best_rmsd:
                                    best_rmsd = db_to_query_rmsd
                                    best_resind_perm = resind_perm
                                    # orig_mol_inds describes cg site; perm_inds describes 
                                    # cg perm
                                    best_cg_perm = (perm_inds, orig_mol_inds)
                                    best_transf = db_transf

                    # If we found a best rmsd, write out the vdg
                    if best_rmsd is not None:
                        resind_perm = best_resind_perm
                        cg_site_and_perm = best_cg_perm
                        db_transf = best_transf
                        db_to_query_rmsd = best_rmsd
                        # Determine output name.
                        q_bsr_perm_str = ''
                        for _bsr, resname in zip(bsrlabel, bsrAAs):
                            resseg, reschain, resnum = _bsr
                            q_bsr_perm_str += f'_{resname}_{resseg}_{reschain}_{resnum}'
                        out_subdir = f'{sub_smiles}{q_bsr_perm_str}'
                        subdir_path = os.path.join(output_dir, sub_smiles, out_subdir)
                        os.makedirs(subdir_path, exist_ok=True)
                        output_vdg_name = (
                            f"{pdbfile.split('/')[-1].removesuffix('.pdb')}_"
                            f"{sub_smiles}_"
                            f"{vdg_path.rstrip('/').split('/')[-1].removesuffix('.pdb.gz')}_"
                            f"{db_to_query_rmsd:.2f}")
                        output_vdg_path = os.path.join(subdir_path, output_vdg_name+'.pdb')
                        if os.path.exists(output_vdg_path):
                            print('WARNING: overwriting existing file', output_vdg_path)
                        pr.writePDB(output_vdg_path, pr.applyTransformation(db_transf, 
                                                                    vdg_prody_obj.copy()))

    log_text = buf.getvalue()
    return (pdbfile, log_text, failed_vdg_files)

def check_frag_in_vdg_lib(sub_smiles, frags_in_lib, vdg_lib_dir):
    '''Determine whether this sub_smiles has been processed and completed in the 
    vdg_lib_dir and mark its status in the frags_in_lib dict.'''
    # Step 0) is this frag named something else in the vdg lib?
    db_name = sub_smiles

    for existingfrag in _listdir_cached(vdg_lib_dir):
        if utils.smiles_equiv(existingfrag, sub_smiles):
            db_name = existingfrag
            break

    # Is it in frags_to_exclude?
    if Frags.check_in_exclude_list(db_name, frags_to_exclude):
        # does internal check for smiles equivalence. 
        return frags_in_lib, db_name

    # Is it in vdg lib and the job didn't break?
    if db_name not in frags_in_lib:
        # Is it in the vdg lib dir?
        if db_name in _listdir_cached(vdg_lib_dir):
            # Did the job (to create the vdgs) finish w/o issue?
            if Frags.check_vdg_job_status(db_name, vdg_lib_dir):
                frags_in_lib[db_name] = True
            else:
                frags_in_lib[db_name] = False
        else:
            # If the str isn't in os.listdir(), it's still possible that a 
            # degenerate smiles is in vdg_lib_dir
            for existingfrag in _listdir_cached(vdg_lib_dir):
                # Is it equivalent to sub_smiles and was its job completed?
                is_equivalent = utils.smiles_equiv(existingfrag, db_name)
                if is_equivalent:  # rename the smiles to match what's in db
                    db_name = existingfrag
                if is_equivalent and Frags.check_vdg_job_status(existingfrag, vdg_lib_dir):
                    frags_in_lib[db_name] = True
                    break
            else:
                frags_in_lib[db_name] = False

    return frags_in_lib, db_name

def should_run_frag(sub_smiles, frags_to_include, frags_in_lib):
    ''' Determine whether to sample this frag or not based on the yml inputs. 
    `frags_to_include` must be 'all' or a list of SMILES.'''
    # First of all, rule out if not allowed by the yml input 
    if frags_to_include == 'all':
        pass
    elif isinstance(frags_to_include, list):
        if not any(utils.smiles_equiv(sub_smiles, inc) for inc in frags_to_include):
            return False
    else:
        raise ValueError("frags_to_include must be 'all' or a list of SMILES.")

    # Now, check if this frag is in the vdg lib and its job is completed
    if sub_smiles not in frags_in_lib:
        # Meaning: not being searched for (possibly being excluded in yml file)
        return False
    if not frags_in_lib[sub_smiles]:
        # Meaning: not in vdg lib
        return False
    return True

def main():
    
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

    # Parallelize per-PDB processing. Expose to workers so they inherit 
    # precomputed objects (on POSIX fork).
    global _SOLVED_STRUCT, _ALL_BSR_COMBOS
    _SOLVED_STRUCT = solved_struct
    _ALL_BSR_COMBOS = all_bsr_combos

    # Build list of pdbfiles to run (respecting query_pdbs)
    all_pdbs = sorted(os.listdir(query_pdbs_dir))
    if query_pdbs == 'all':
        pdb_list = all_pdbs
    elif isinstance(query_pdbs, list):
        pdb_list = [p for p in all_pdbs if p in query_pdbs]
    else:
        raise ValueError("`query_pdbs` must be 'all' or a list of pdb files.")

    if not pdb_list:
        print("No pdb files selected to run.")
        return

    os.makedirs(outdir, exist_ok=True)

    # Print the per-fragment lines once (using the first pdb) and show ONE shared header
    first_pdb = pdb_list[0]
    prelog_buf = io.StringIO()
    with redirect_stdout(prelog_buf), redirect_stderr(prelog_buf):
        # print the per-fragment lines once
        orig_filtered_frags_once, _pdbmol = Frags.get_frags_from_pdbfile(
            os.path.join(query_pdbs_dir, first_pdb), lig_smiles, quiet=False)

        # compute the same frags_in_lib mapping we do per-PDB, then summarize once
        frags_in_lib_once = {}
        dbname_keys_once = {}
        for sub_smiles in orig_filtered_frags_once.keys():
            frags_in_lib_once, db_name = check_frag_in_vdg_lib(sub_smiles, frags_in_lib_once, 
                                                                 vdg_lib_dir)
            dbname_keys_once[sub_smiles] = db_name

        Frags.summarize_frags(frags_in_lib_once, frags_to_exclude, frags_to_include)

    print("Fragment list and summary:")
    header_text = prelog_buf.getvalue()
    if header_text:
        print(header_text, end="" if header_text.endswith("\n") else "\n", flush=True)

    # Multiprocessing
    n_workers = min(len(pdb_list), max(1, mp.cpu_count() - 1))
    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(_process_one_pdb, pdb_list)  # preserves order of pdb_list

    for pdbfile, log_text, _failed in results:
        sep = f"\n==================== {pdbfile} ====================\n"
        sys.stdout.write(sep)
        if log_text:
            sys.stdout.write(log_text)
            if not log_text.endswith("\n"):
                sys.stdout.write("\n")
        sys.stdout.flush()

    # Union failures across workers
    all_failed = set()
    for _pdbfile, _log_text, _failed in results:
        all_failed.update(_failed)

    if all_failed:
        print("\nvdG files that could not be parsed (unique across all jobs):")
        for f in sorted(all_failed):
            print(f" - {f}")
        print("", flush=True)

    return

if __name__ == "__main__":
    main()
