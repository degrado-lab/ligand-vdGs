from itertools import combinations, product
import numpy as np
import os
import utils

def get_bsr_combinations(solved_struct, ligname):
    bindingsite_residues = get_bindingsite_residues(solved_struct,
        [], ligname, dist_from_lig=4.5, CA_only=False)
    bsr_combos = get_vdg_subsets(bindingsite_residues)
    all_bsr_combos = []
    for bsr_combo in bsr_combos:
        bsr_AA_identities = []
        input_bsr_bb_coords = []
        if len(bsr_combo) != 2:
            continue
        for bsr in bsr_combo: 
            bsr_seg, bsr_chain, bsr_resnum = bsr
            res_obj = select_residue(solved_struct, bsr_seg, bsr_chain, bsr_resnum)
            res_bb_coords = []
            for atom_name in ['N', 'CA', 'C']:
                res_bb_coords.append(utils.get_atom_coords(res_obj, atom_name))
            input_bsr_bb_coords.append(res_bb_coords)
            AA = get_res_AA_identity(res_obj)
            bsr_AA_identities.append(AA)
        bsr_AA_identities = [AA if AA != 'GLY' else 'bb' for AA in bsr_AA_identities]
        # All of the AAs have potential to provide bb contacts, so sample "bb" vdms
        options = [(x,) if x == 'bb' else (x, 'bb') for x in bsr_AA_identities]
        combo_variants = list(product(*options))
        for c in combo_variants:
            if c not in all_bsr_combos:
                all_bsr_combos.append((c, input_bsr_bb_coords.copy())) 
    return all_bsr_combos

def get_bindingsite_residues(prody_obj, addl_residues, ligname, dist_from_lig=8,
                             CA_only=True):
    res = []
    # Use CA_only when you're doing blind docking and don't want to use sc info
    # Use CA_only=False when you want to use sc positions to define interactions
    if CA_only: 
        _selection = 'name CA'
    else: 
        _selection = 'protein'
    CAs = prody_obj.select(
        f'{_selection} within {dist_from_lig} of resname {ligname} and not element CA') # exclude calcium
    for ca in CAs:
        res_tup = (ca.getSegname(), ca.getChid(), ca.getResnum())
        if res_tup not in res:
            res.append(res_tup)
    if addl_residues:
        res += addl_residues 
    for r in res:
        if not isinstance(r, tuple) or len(r) != 3:
            print(f"Invalid residue tuple: {r}")
    print('\nbinding site residues for pymol selection:\n')
    print('select bindingsite, ' + ' or '.join(
    f'(seg {seg} and chain {chain} and resi {resnum})' for seg, chain, resnum in res) + '\n')
    return res

def get_vdg_subsets(input_list):
    # Initialize an empty list to store all subsets
    all_subsets = []
    # Loop through subset sizes 2 to 3, generate combos, then add to list
    for r in range(2, 4): # exclude subset size 1
        subsets = combinations(input_list, r)
        all_subsets.extend(subsets)
    return all_subsets

def select_residue(prody_obj, seg, chain, resnum):
    if seg == '':
        sele = f'chain {chain} and resnum {resnum}'
    else:
        sele = f'segment {seg} and chain {chain} and resnum {resnum}'
    res_obj = prody_obj.select(sele)
    return res_obj


def get_res_AA_identity(res_obj):
    assert len(set(res_obj.getResnames())) == 1
    AA = res_obj.getResnames()[0]
    return AA

def get_query_cg_coords(mol, sub):
    query_matches = mol.GetSubstructMatches(sub, uniquify=False) # tuple of atom inds
    q_cg_coord_perms = [] # each element is a permutation
    for query_match in query_matches: # each permutation
        perm_coords = []
        # Get atom names, etc.
        conf = mol.GetConformer()  # Get the 3D conformer to get coords 
        assert len(query_match) == sub.GetNumAtoms()
        
        for atom_idx in query_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            pos = conf.GetAtomPosition(atom_idx)  # returns an RDKit Point3D object
            xyz = (pos.x, pos.y, pos.z)
            perm_coords.append(xyz)
        q_cg_coord_perms.append(np.array(perm_coords))
    return q_cg_coord_perms

def name_outdir(pdbfile, outdir, num_query_structs):
    # Name the outdir for this pdb query
    pdbname = pdbfile.split('/')[-1].removesuffix('.pdb')
    pdb_id = pdbname[:4]
    if num_query_structs > 1:
        output_dir = os.path.join(outdir, pdb_id, pdbname)
    else:
        output_dir = os.path.join(outdir, pdbname)
    return output_dir

def map_aa_identities_to_vdg_resinds(vdg_AAs, target_list, wildcard='bb'):
    '''
    Given a list of amino acids in the prody obj (`vdg_AAs`) and a list of AAs to 
    select for (`target_list`), find all permutations of indices (i.e., resindices) in the 
    vdg (`vdg_AAs`) that match the target list elements in order. The amino acid `bb` 
    is treated as a wildcard.
    '''

    # Find all indices (resindices)for each unique element (AA) in the target list
    AA_indices = {}

    # Handle regular AAs (non-wildcards, i.e. not `bb`)
    for AA in set(target_list):
        if AA != wildcard:
            AA_indices[AA] = [i for i, x in enumerate(vdg_AAs) if x == AA]
    
    # For wildcard, it can match any AA in the vdg (exclude lig)
    AA_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
                'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                'THR', 'TRP', 'TYR', 'VAL']
    if wildcard in target_list:
        AA_indices[wildcard] = list(range(len(vdg_AAs)))
        AA_indices[wildcard] = [i for i, x in enumerate(vdg_AAs) if x in AA_codes]
    
    # Generate all permutations of indices
    result = []

    def backtrack(current_indices, target_position):
        # If we've matched all target elements, add the current permutation to result
        if target_position == len(target_list):
            result.append(current_indices.copy())
            return
        
        # Get the next element to match
        current_element = target_list[target_position]
        
        # Iterate over each possible index for this element
        for idx in AA_indices.get(current_element, []):
            # Make sure we don't use the same index twice
            if idx not in current_indices:
                current_indices.append(idx)
                backtrack(current_indices, target_position + 1)
                current_indices.pop()  # Backtrack
    backtrack([], 0)
    return result