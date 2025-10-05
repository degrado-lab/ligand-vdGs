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

_VDG_DIR_CACHE = {}         # path -> listdir(list[str])
_VDG_PARSE_CACHE = {}       # path -> (vdg_prody_obj, vdg_AAs, 
                                # vdg_db_cg_coords(np.ndarray), 
                                # bb_cache{resindex:{'N':arr,'CA':arr,'C':arr}})

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
    # Not sure why, but some vdg pdbs have duplicate residues. Resindex 0 is always 
    # the lig, and the 1st and 2nd resindices are the bsrs. Select just those indices 
    # and ignore the duplicated residues.
    resindices = resindices[:3] # first 3 residues are always lig + 2 bsrs. 
    vdg_AAs = [vdg_prody_obj.select(f'resindex {_r}').getResnames()[0] for _r in resindices]
    # Precompute backbone atom coords per residue once
    bb_cache = {}
    for _r in [1, 2]: # resindex 0 = lig
        sel = vdg_prody_obj.select(f'resindex {_r}')
        bb_cache[_r] = {
            'N':  utils.get_atom_coords(sel, 'N'),
            'CA': utils.get_atom_coords(sel, 'CA'),
            'C':  utils.get_atom_coords(sel, 'C'),
        }
    vdg_db_cg_coords = np.array(match.get_database_cg_coords(vdg_prody_obj))
    _VDG_PARSE_CACHE[vdg_full_path] = (vdg_prody_obj, vdg_AAs, vdg_db_cg_coords, bb_cache)
    return _VDG_PARSE_CACHE[vdg_full_path]

def check_frag_in_vdg_lib(sub_smiles, frags_in_lib, vdg_lib_dir):
    '''Determine whether this sub_smiles has been processed and completed in the 
    vdg_lib_dir and mark its status in the frags_in_lib dict.'''
    # Step 0) is this frag named something else in the vdg lib?
    for existingfrag in _listdir_cached(vdg_lib_dir):
        if utils.smiles_equiv(existingfrag, sub_smiles):
            sub_smiles = existingfrag
            break

    # Is it in frags_to_exclude?
    if Frags.check_in_exclude_list(sub_smiles, frags_to_exclude): 
        # does internal check for smiles equivalence. 
        return frags_in_lib

    # Is it in vdg lib and the job didn't break?
    if sub_smiles not in frags_in_lib:
        # Is it in the vdg lib dir?
        if sub_smiles in _listdir_cached(vdg_lib_dir): 
            # Did the job (to create the vdgs) finish w/o issue?
            if Frags.check_vdg_job_status(sub_smiles, vdg_lib_dir):
                frags_in_lib[sub_smiles] = True
            else:
                frags_in_lib[sub_smiles] = False 
        else: 
            # If the str isn't in os.listdir(), it's still possible that a 
            # degenerate smiles is in vdg_lib_dir
            for existingfrag in _listdir_cached(vdg_lib_dir):
                # Is it equivalent to sub_smiles and was its job completed?
                is_equivalent = utils.smiles_equiv(existingfrag, sub_smiles)
                if is_equivalent: # rename the smiles to match what's in db
                    sub_smiles = existingfrag
                if is_equivalent and Frags.check_vdg_job_status(existingfrag, vdg_lib_dir):
                    frags_in_lib[sub_smiles] = True
                    break
            else:
                frags_in_lib[sub_smiles] = False

    return frags_in_lib

def should_run_frag(sub_smiles, frags_to_include, frags_in_lib):
    ''' Determine whether to sample this frag or not based on the yml inputs. 
    `frags_to_include` must be 'all' or a list of SMILES.'''
    # First of all, rule out if not allowed by the yml input 
    if frags_to_include == 'all':
        pass
    elif isinstance(frags_to_include, list):
        # If frags_to_include is a list, check if sub_smiles is in the list
        if sub_smiles not in frags_to_include:
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
    
    # Keep track of which frags are in vdg lib and which aren't
    frags_in_lib = {}

    # Keep track of which database vdg files could not be loaded
    failed_vdg_files = []

    # Load predictions. For each pred, get the frags and CG coords in the pred, iterate 
    # through binding site res combos, align, and determine if there are vdg matches. 
    for pdbfile in sorted(os.listdir(query_pdbs_dir)):
        if query_pdbs == 'all':
            pass
        elif isinstance(query_pdbs, list):
            # If query_pdbs is a list, check if pdbfile is in the list
            if pdbfile not in query_pdbs:
                continue
        else:
            raise ValueError("query_pdbs must be 'all' or a list of pdb files.")
        print('\n' + '='*10, pdbfile, '='*10)
        output_dir = dock.name_outdir(pdbfile, outdir, make_pdb_subfolder=True)
        utils.set_up_outdir(output_dir, overwrite=overwrite_existing)

        # Get all fragment smiles, Mol objs (to later get coords), and perms 
        orig_filtered_frags, pdb_mol_reassigned = Frags.get_frags_from_pdbfile(
            os.path.join(query_pdbs_dir, pdbfile), lig_smiles)
        filtered_frags = {}
        for sub_smiles, substruct_perm_grouped_by_site in orig_filtered_frags.items():
            frags_in_lib = check_frag_in_vdg_lib(sub_smiles, frags_in_lib, vdg_lib_dir)
            
            if frags_in_lib[sub_smiles]: 
                filtered_frags[sub_smiles] = substruct_perm_grouped_by_site

        # Log
        Frags.summarize_frags(frags_in_lib, frags_to_exclude, frags_to_include)

        # For every frag that's in the vdg lib, iterate over the bsr combos and determine 
        # if there are vdg matches.
        for sub_smiles, substruct_perm_grouped_by_site in filtered_frags.items():
            should_run = should_run_frag(sub_smiles, frags_to_include, frags_in_lib)
            # Determine frags to skip over
            if not should_run:
                print(f'Skipping {sub_smiles}..........................\n')
                continue
            # Otherwise, proceed with finding vdg matches for this frag
            print(f'Searching for {sub_smiles}..........................\n')
            # Get the query CG coords for this substructure in the query struct. 
            grouped_q_cg_perms = [] # all instances of CG + perms, grouped by site on lig
            for substruct_site in substruct_perm_grouped_by_site:
                q_cg_perms_single_site = []
                for substruct_perm in substruct_site:
                    sub, perm_inds, orig_mol_inds = substruct_perm
                    # Get the query CG coords for this substructure in the query struct. 
                    single_perm_cg_coords = dock.get_query_cg_coords(sub, sub_smiles)
                    q_cg_perms_single_site.append((single_perm_cg_coords, perm_inds, orig_mol_inds))
                grouped_q_cg_perms.append(q_cg_perms_single_site)

            # Iterate over the query structure's binding site residues and superpose vdGs 
            # (incl. the CG) onto it, keeping track of the vdGs that have low rmsd to the 
            # known bb+CG.
            for combo, bsrlabel, bsrAAs, coords in all_bsr_combos:
                # Retain pairs only
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
                # Iterate over each vdg. Each vdg represents a pair of bsrs.
                vdg_paths = [p for p in _listdir_cached(vdgs_path) if p.endswith(
                    ('.pdb', '.pdb.gz'))]
                for vdg_path in vdg_paths:
                    vdg_full_path = os.path.join(vdgs_path, vdg_path)
                    try: # parse and precompute once 
                        vdg_prody_obj, vdg_AAs, vdg_db_cg_coords, bb_cache = _precompute_vdg(
                            vdg_full_path)
                    except Exception as e:
                        failed_vdg_files.append(vdg_full_path)
                        continue
                    
                    # Superpose the vdg backbone atoms onto the input structure. Make sure 
                    # the order of query vdg AAs is consistent with order of input_pdb AAs.
                    # Use cached AA list & cached mapper.
                    
                    resind_perms = dock.map_aa_identities_to_vdg_resinds(vdg_AAs,
                        bsr_incl_bb_identities.split('_'))

                    # Determine if this vdg has bb+cg geometry that matches the query struct.
                    # Iterate over each CG site and each CG atom perm within that site. 
                    # Keep track of the lowest rmsd result, b/c there are cases where >1 
                    # resind perms and/or cg perms yield rmsd < threshold.
                    best_rmsd, best_resind_perm, best_cg_perm, best_transf = [None, None, 
                        None, None] # best_cg_perm will describe both site and atom perm
                    for resind_perm in resind_perms:
                        # build backbone coords from bb_cache to minimize prody overhead
                        vdg_perm_bb_coords = np.array([coord for _r in resind_perm for coord 
                            in (bb_cache[_r]['N'], bb_cache[_r]['CA'], bb_cache[_r]['C'])])
                        
                        # Superpose the database vdg bb + cg atoms onto the input structure 
                        # and calculate rmsd.
                        # --- First, get combined bb+cg coords for db vdg
                        try: # precomputed cg coords once
                            db_bb_and_cg_coords = np.concatenate((vdg_perm_bb_coords, vdg_db_cg_coords), axis=0)
                        except Exception:
                            continue
                        # --- Then, superpose and calc rmsd for all query CG permutations
                        for q_cg_site in grouped_q_cg_perms:
                            for q_cg_perm in q_cg_site:
                                q_cg_coord_perm, perm_inds, orig_mol_inds = q_cg_perm
                                q_bb_and_cg = np.concatenate((input_bsr_bb_coords, q_cg_coord_perm), 
                                                             axis=0)
                                if db_bb_and_cg_coords.shape != q_bb_and_cg.shape:
                                    raise ValueError(
                                        'The query CG must have the same num. of atoms as the '
                                        'database CG.')
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
                                    # orig_mol_inds describes cg site; perm_inds describes cg perm
                                    best_cg_perm = (perm_inds, orig_mol_inds) 
                                    best_transf = db_transf
                    # If we found a best rmsd, write out the vdg with the best rmsd
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
                        # Output vdg name 
                        out_subdir = f'{sub_smiles}{q_bsr_perm_str}'
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

    if failed_vdg_files:
        print('\nvdG files that could not be parsed:')
        for failed_file in failed_vdg_files:
            print(f" - {failed_file}")

    return

if __name__ == "__main__":
    main()
