import os
from itertools import combinations, permutations, product
import prody as pr

def subsets_of_query_residues(bindingsite_residues):
    query_res_sets = [[b] for b in bindingsite_residues] + list(
        combinations(bindingsite_residues, 2))  + list(
        combinations(bindingsite_residues, 3)) + list(
        combinations(bindingsite_residues, 4))
    return query_res_sets

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

def get_database_vdm_bb_coords(database_vdg, query_vdm_resnames, db_vdg_path):
    '''
    Returns bb coordinates of all valid AA permutations that match the order of
    query_vdm_resnames
    '''
    db_vdm_resinds = get_vdg_interacting_resindices(database_vdg)
    db_vdm_firstatoms = [
        database_vdg.select(f'resindex {r}')[0] for r in db_vdm_resinds]
    db_vdm_resnames = [i.getResname() for i in db_vdm_firstatoms]
    # Get the indices of (permuted) db_vdm_resnames that are the same AA identities as 
    # the binding site residues being queried, so that the indices can be mapped to
    # the resindices and eventually be used for extracting backbone coordinates.
    list_match_ind_AA_permutations = find_matching_permutations_indices(query_vdm_resnames, 
                                                                     db_vdm_resnames)
    if not list_match_ind_AA_permutations:
        return None

    all_AA_perms_bb_coords_and_resinds = {}
    for match_inds in list_match_ind_AA_permutations:
        match_resinds = [db_vdm_resinds[ind] for ind in match_inds] # map to the resinds
        bb_coords = [] # if 2 res: [res1 N, res1 CA, res1 C, res2 N, res2 CA, res2C]
        for match_r in match_resinds:
            for bbatom in ['N', 'CA', 'C']:
                bbatom_coords = get_coords_from_resindex(database_vdg, match_r, bbatom, 
                                                         db_vdg_path)
                bb_coords.append(bbatom_coords)
        all_AA_perms_bb_coords_and_resinds[tuple(match_resinds)] = bb_coords

    return all_AA_perms_bb_coords_and_resinds

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

def get_coords_from_resindex(struct, resindex, atomname, pdbpath):
    atom_obj = struct.select(f'resindex {resindex} and name {atomname}')
    if len(atom_obj) == 0:
        raise ValueError(f'"resindex {resindex} and name {atomname}" not found.')
    
    if len(atom_obj) > 1: # not sure why, but sometimes each atom is duplicated,
                          # despite having identical coords.
        ambiguous_coords = None
        for atom in atom_obj:
            if ambiguous_coords is None:
                ambiguous_coords = atom.getCoords()
            else:
                dist = pr.calcDistance(ambiguous_coords, atom.getCoords())
                if dist > 0.2:
                    raise ValueError(f'The selection "resindex {resindex} and name '
                        f' {atomname}" is expected to select one atom, but it corresponds '
                        f'to {len(atom_obj)} atoms in {pdbpath}.')
    
    return list(atom_obj.getCoords()[0])

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

def get_in_dir_and_out_dir(q_res_set, vdgs_dir):
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
    return _vdgs_dir, subdir

def get_database_vdgs_for_spec_AAs(_vdgs_dir, AAs):
    AA_str = '_'.join(sorted(AAs))
    print(AAs)
    database_vdg_names = os.listdir(os.path.join(_vdgs_dir, AA_str))
    database_vdg_paths = [os.path.join(_vdgs_dir, AA_str, d) for d in database_vdg_names]
    #pdbnames = [i.split('/')[-1] for i in database_vdg_paths]
    #assert len(database_vdg_paths) == len(set(pdbnames))
    database_vdgs = [pr.parsePDB(d) for d in database_vdg_paths]
    assert len(database_vdg_paths) == len(database_vdgs)
    return database_vdg_paths, database_vdgs

def determine_permutation_inds(symm_atoms, list_atom_names):
    # Return nested list of indices to permute on, based on atom names declared to be 
    # symmetrically equivalent
    permute_indices = []
    for group in symm_atoms:
        permute_indices.append([list_atom_names.index(atom) for atom in group])
    return permute_indices

def permute_cg_coords(coords, cg_atom_names, symm_atoms):
    # `symm_atoms` is a list of lists of atoms to swap
    if len(symm_atoms) == 0:
        return [coords]
    
    # Make sure `symm_atoms` are in `cg_atoms`
    flattened_symm_atoms = [item for sublist in symm_atoms for item in sublist]
    symm_atoms_in_cg_atoms = [i for i in flattened_symm_atoms if i in cg_atom_names]
    if len(symm_atoms_in_cg_atoms) != len(flattened_symm_atoms):
        raise ValueError('All atoms specified in `symm_atoms` must be in `cg_atoms`.')
    if len(flattened_symm_atoms) != len(set(flattened_symm_atoms)):
        raise ValueError('Symmetric atoms should not be specified more than once.')

    # Generate permutations 
    permute_indices = determine_permutation_inds(symm_atoms, cg_atom_names)
    perm_group_names = [] 

    for indices in permute_indices:
        extracted_names = [cg_atom_names[i] for i in indices]
        perm_group_names.append(list(permutations(extracted_names)))
    
    # Combine the permutations of each group
    all_perm_atom_names = []
    for perm_combination in product(*perm_group_names):  
        new_order_atoms = cg_atom_names[:]
        for group_idx, indices in enumerate(permute_indices):
            for i, index in enumerate(indices):
                new_order_atoms[index] = perm_combination[group_idx][i]
        all_perm_atom_names.append(new_order_atoms)
    
    # Extract the coords based on the indices of the permuted atom names
    all_perm_coords = []
    for perm in all_perm_atom_names:
        perm_coords = [coords[cg_atom_names.index(atom)] for atom in perm]
        all_perm_coords.append(perm_coords)
    
    return all_perm_atom_names, all_perm_coords