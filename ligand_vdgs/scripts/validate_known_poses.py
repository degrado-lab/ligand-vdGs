'''
Given a query structure (a crystal structure of a protein-ligand complex) and its binding site 
residues, determine whether any vdGs can correctly place the chemical groups within the 
ligand of the query structure.

Usage:
    >> python validate_known_poses.py $YOUR_YAML_FILE
'''

import sys
import os
import yaml
from itertools import combinations, permutations
import numpy as np
import prody as pr

script, yml_file = sys.argv

def main():
    # Parse inputs from yml file
    vdgs_dir, cg_name, query_path, bindingsite_residues, query_lig_res, \
               query_cg_atoms, rmsd_threshold, out_dir = parse_yaml_input(yml_file)
    
    # Set up output dir
    pdbname = query_path.split('/')[-1].removesuffix('.pdb')
    output_dir = os.path.join(out_dir, yml_file.removesuffix('.yml').split('/')[-1])
    if os.path.exists(output_dir):
        if os.listdir(output_dir):
            raise ValueError(f'The output directory {output_dir} is not empty. '
                             'Please remove its contents or specify a new output '
                             'directory to prevent accidental overwriting.')
    else: 
        os.makedirs(output_dir)

    # Store the coords of the query structure's ligand chemical group
    query_struct = pr.parsePDB(query_path)
    q_lig_seg, q_lig_ch, q_lig_resnum = query_lig_res
    q_cg_coords = get_coords(query_struct, query_cg_atoms, q_lig_seg, q_lig_ch, 
                           q_lig_resnum, query_path)
    
    # Iterate over the query structure's binding site residues and superpose vdGs 
    # (incl. the CG) onto it, keeping track of the vdGs that have low rmsd to the 
    # known bb+CG.
    query_bb_coords = {}
    for bsr in bindingsite_residues:
        bsr_seg, bsr_chain, bsr_resnum = bsr
        bbcoords = get_coords(query_struct, ['N', 'CA', 'C', 'O'], bsr_seg, bsr_chain, bsr_resnum,
                              query_path)
        query_bb_coords[tuple(bsr)] = bbcoords
    
    # Get subsets of vdg residues
    query_res_sets = [[b] for b in bindingsite_residues] + list(
        combinations(bindingsite_residues, 2))  + list(
        combinations(bindingsite_residues, 3)) + list(
        combinations(bindingsite_residues, 4)
        )
    '''
    # If only considering pairs: # note that there may be FGs in solved structures
    #                            # with only 1 vdM.
    query_res_sets = list(combinations(bindingsite_residues, 2))
    '''
    
    # Iterate over subsets of binding site residues and determine
    # whether any of the vdGs in the database match the "ground truth" binding site
    # residues and CG positions of the query structure.
    for q_res_set in query_res_sets:

        if len(q_res_set) == 1:
            subdir = 'single_res_matches'
            _vdgs_dir = os.path.join(rf"{vdgs_dir}", 'nr_vdgs', '1')

        elif len(q_res_set) == 2:
            subdir = 'double_res_matches'
            _vdgs_dir = os.path.join(rf"{vdgs_dir}", 'nr_vdgs', '2')
        elif len(q_res_set) == 3:
            subdir = 'triple_res_matches'
            _vdgs_dir = os.path.join(rf"{vdgs_dir}", 'nr_vdgs', '3')
        elif len(q_res_set) == 4:
            subdir = 'quad_res_matches'
            _vdgs_dir = os.path.join(rf"{vdgs_dir}", 'nr_vdgs', '4')
        else:
            assert False

        database_vdg_names = os.listdir(_vdgs_dir)
        database_vdg_paths = [os.path.join(_vdgs_dir, d) for d in database_vdg_names]
        pdbnames = [i.split('/')[-1] for i in database_vdg_paths]
        assert len(database_vdg_paths) == len(set(pdbnames))
        database_vdgs = [pr.parsePDB(d) for d in database_vdg_paths]

        # Calc rmsd on subsets of bb residues (incl. CG) 
        q_bbcoords_list = [query_bb_coords[tuple(r)] for r in q_res_set]
        q_flattened_bb_list = [item for sublist in q_bbcoords_list for item in sublist]
        q_bb_and_cg = np.array(q_flattened_bb_list + q_cg_coords)
        
        # Determine the AA resnames of the query vdms.
        query_vdm_resnames = []
        for q_res in q_res_set:
            bsr_seg, bsr_chain, bsr_resnum = q_res
            q_res_obj = get_residue_obj(query_struct, bsr_seg, bsr_chain, bsr_resnum, 
                                        query_path)
            query_vdm_resnames.append(q_res_obj.getResnames()[0])
        
        # Sample the database vdGs.
        for database_vdg_index, database_vdg in enumerate(database_vdgs):
            db_vdg_path = database_vdg_paths[database_vdg_index]
            
            # Get CG coordinates
            db_cg_coords = get_database_cg_coords(database_vdg)
            if len(db_cg_coords) != len(q_cg_coords):
                print('Incorrect number of database CG coordinates.')
                continue
                raise ValueError(
                    'The query CG must have the same num. of atoms as the database CG.')
            
            # Get vdM (vdG residues that interact w/ CG) bb coordinates if they're the
            # same AA identity as the binding site residues being queried. Make sure all
            # valid permutations are sampled.
            all_permutations_db_vdm_bb_coords_and_resinds = get_database_vdm_bb_coords(
                database_vdg, query_vdm_resnames, db_vdg_path)
            if all_permutations_db_vdm_bb_coords_and_resinds is None:
                continue
            
            # Iterate over the vdM permutations and superpose/calc RMSD against the
            # known structure (incl. CG).
            for perm_resinds, permutation_db_vdm_bb_coords in \
                                all_permutations_db_vdm_bb_coords_and_resinds.items():
                db_vdg_bb_and_cg = np.array(permutation_db_vdm_bb_coords + db_cg_coords)
                moved_db_vdg_bb_cg_coords, transf = pr.superpose(db_vdg_bb_and_cg, 
                                                                 q_bb_and_cg)
                rmsd = round(pr.calcRMSD(moved_db_vdg_bb_cg_coords, q_bb_and_cg), 2)
                if rmsd < rmsd_threshold:
                    moved_vdg = pr.applyTransformation(transf, database_vdg.copy())
                    
                    # Write out just the vdms in this permutation
                    perm_resinds_str = ' '.join([str(i) for i in perm_resinds])
                    perm_resinds_sele = f'resindex {perm_resinds_str}'
                    moved_vdg_perm_resinds_obj = moved_vdg.copy().select(
                        f'({perm_resinds_sele}) or (occupancy > 2.8)') # occ 3 = CG
                    vdg_pdbname = db_vdg_path.rstrip('/').split('/')[-1].removesuffix('.pdb')
                    output_vdg_name = f'{vdg_pdbname}'
                    for residue in q_res_set:
                        bsr_flat_list = [str(it) for it in residue]
                        bsr_str = '_'.join(bsr_flat_list)
                        output_vdg_name += f'_{bsr_str}'
                    subdir_path = os.path.join(output_dir, subdir)
                    output_vdg_path = os.path.join(output_dir, subdir, output_vdg_name)
                    if not os.path.exists(subdir_path):
                        os.makedirs(subdir_path)
                    pr.writePDB(output_vdg_path, 
                                moved_vdg_perm_resinds_obj.copy().select(
                                    'occupancy > 1.5'))
    
def get_database_vdm_bb_coords(database_vdg, query_vdm_resnames, db_vdg_path):
    '''
    Returns bb coordinates of all valid permutations that match the order of
    query_vdm_resnames
    '''
    db_vdm_resinds = get_vdg_interacting_resindices(database_vdg)
    db_vdm_firstatoms = [
        database_vdg.select(f'resindex {r}')[0] for r in db_vdm_resinds]
    db_vdm_resnames = [i.getResname() for i in db_vdm_firstatoms]
    # Get the indices of (permuted) db_vdm_resnames that are the same AA identities as 
    # the binding site residues being queried, so that the indices can be mapped to
    # the resindices and eventually be used for extracting backbone coordinates.
    list_match_ind_permutations = find_matching_permutations_indices(query_vdm_resnames, 
                                                                     db_vdm_resnames)
    if not list_match_ind_permutations:
        return None

    all_permutations_bb_coords_and_resinds = {}
    for match_inds in list_match_ind_permutations:
        match_resinds = [db_vdm_resinds[ind] for ind in match_inds] # map to the resinds
        bb_coords = [] # if 2 res: [res1 N, res1 CA, res1 C, res2 N, res2 CA, res2C]
        for match_r in match_resinds:
            for bbatom in ['N', 'CA', 'C', 'O']:
                bbatom_coords = get_coords_from_resindex(database_vdg, match_r, bbatom, 
                                                         db_vdg_path)
                bb_coords.append(bbatom_coords)
        all_permutations_bb_coords_and_resinds[tuple(match_resinds)] = bb_coords

    return all_permutations_bb_coords_and_resinds

def get_coords_from_resindex(struct, resindex, atomname, pdbpath):
    atom_obj = struct.select(f'resindex {resindex} and name {atomname}')
    if len(atom_obj) != 1:
        raise ValueError(f'The selection "resindex {resindex} name {atomname}" '
                 'is expected to select one atom, but it corresponds to '
                 f'{len(atom_obj)} atoms in {pdbpath}.')
    
    return list(atom_obj.getCoords()[0])

def find_matching_permutations_indices(list1, list2):
    '''
    Find all permutations of elements in `list2` that match the elements of `list1`
    and record their indices.
    '''
    if len(list1) > len(list2):
        return None
    
    len_list1 = len(list1)
    all_permutations = permutations(range(len(list2)), len_list1)
    matching_indices = []
    
    for perm_indices in all_permutations:
        perm_elements = [list2[i] for i in perm_indices]
        if perm_elements == list1:
            matching_indices.append(list(perm_indices))
    
    return matching_indices

def get_vdg_interacting_resindices(database_vdg):
    '''Returns resindices of residues (vdms) interacting with the CG'''
    vdg_resinds = []
    for v in database_vdg:
        occ = v.getOccupancy()
        # A vdg residue interacting with the CG has occ == 2, but sometimes 
        # there are numerical errors.
        if occ > 1.9 and occ < 2.1:
            resind = v.getResindex()
            if resind not in vdg_resinds:
                vdg_resinds.append(resind)
    return vdg_resinds

def get_residue_obj(struct, seg, chain, resnum, pdbpath):
    if seg:
        sel_str = f'segment {seg} and chain {chain} and resnum {resnum}'
    else:
        sel_str = f'chain {chain} and resnum {resnum}'
    res_obj = struct.select(sel_str)
    num_res_selected = len(set(res_obj.getResindices())) # Make sure only 1 res is selected
    if num_res_selected != 1:
        raise ValueError(f'The residue selection "segment {seg} chain {chain} resnum {resnum}" '
                         'is expected to select one residue, but it corresponds to '
                         f'{num_res_selected} residues in {pdbpath}.')
    return res_obj

def get_coords(struct, atomnames, seg, ch, res, pdbpath):
    coords = []
    for atom in atomnames:
        if seg:
            atom_sel_str = f'segment {seg} and chain {ch} and resnum {res} and name {atom}'
        else:
            atom_sel_str = f'chain {ch} and resnum {res} and name {atom}'

        atom_sel = struct.select(atom_sel_str) 
        if len(atom_sel) != 1 or atom_sel is None:
            raise ValueError(f'segment {seg} chain {ch} resnum {res} name {atom} is '
                             f'expected to select one atom, but it corresponds to '
                             f'{len(atom_sel)} atoms in {pdbpath}')
        coord = atom_sel.getCoords()[0]
        coords.append(list(coord)) # for ease of adding/flattening lists and not arrays
    return coords

def get_database_cg_coords(database_vdg):
    '''
    Returns the coordinates of CG atoms ordered according to the input SMARTS pattern.
    CG atoms have their occupancies set to 3.0 (but select >2.9 because sometimes,
    there are numerical errors), and the order of the value of each CG atom's occupancy
    corresponds to the order it's represented in the SMARTS string. For example, 
    if the occupancies are [3.1, 3.0, 3.2], the order of the atoms should be indices
    [1, 0, 2].
    '''
    cg_atoms = [i for i in database_vdg if i.getOccupancy() > 2.9]
    occs = [i.getOccupancy() for i in cg_atoms]
    sorted_indices = sorted(range(len(occs)), key=lambda i: occs[i])
    ordered_cg_atoms = [cg_atoms[ind] for ind in sorted_indices]
    ordered_cg_atom_coords = [list(a.getCoords()) for a in ordered_cg_atoms]
    return ordered_cg_atom_coords

def parse_yaml_input(yml_file):
    with open(yml_file, 'r') as inF:
        user_data = yaml.load(inF, Loader=yaml.FullLoader)
    vdgs_dir = user_data['vdgs_dir']
    cg_name = user_data['cg_name']
    query_path = user_data['query_path']
    bindingsite_residues = user_data['bindingsite_residues']
    query_lig_res = user_data['query_lig']
    query_cg_atoms = user_data['query_cg_atoms']
    rmsd_threshold = user_data['rmsd_threshold']
    out_dir = user_data['out_dir']
    return vdgs_dir, cg_name, query_path, bindingsite_residues, query_lig_res, \
           query_cg_atoms, rmsd_threshold, out_dir

if __name__ == "__main__":
    main()
