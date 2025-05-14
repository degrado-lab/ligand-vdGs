import os
import yaml
import sys
import argparse
from itertools import combinations
import time
import numpy as np
import prody as pr

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', "--yml-file", type=str, required=True,
                        help="Yaml file containing binding pocket information.")
    parser.add_argument('-r', "--rmsd", type=float, default=1.3,
                        help="RMSD cutoff of vdg bb alignment to define a match.")
    return parser.parse_args()

def main():
    print('WARNING: ONLY PROCESSING FRAGMENTS STARTING WITH "g" FOR TESTING PURPOSES')
    # Parse inputs from yml file
    args = parse_args()
    cutoff = args.rmsd
    yml_file = args.yml_file
    (vdgs_dir, input_pdb_path, bindingsite_residues) = parse_yaml_input(yml_file)
    # Parse the pdb and iterate through the combinations of binding site residues
    input_pdb = pr.parsePDB(input_pdb_path)
    bsr_combos = get_vdg_subsets(bindingsite_residues)
    # Iterate over binding site residue combinations
    for bsr_combo in bsr_combos:
        bsr_AA_identities = []
        input_bsr_bb_coords = []
        for bsr in bsr_combo: 
            bsr_seg, bsr_chain, bsr_resnum = bsr
            res_obj = select_residue(input_pdb, bsr_seg, bsr_chain, bsr_resnum)
            for atom_name in ['N', 'CA', 'C']:
                input_bsr_bb_coords.append(get_atom_coords(res_obj, atom_name))
            AA = get_res_AA_identity(res_obj)
            bsr_AA_identities.append(AA)
        input_bsr_bb_coords = np.array(input_bsr_bb_coords)
        bsr_AA_identities = [AA if AA != 'GLY' else 'bb' for AA in bsr_AA_identities]
        bsr_AA_identities = '_'.join(sorted(bsr_AA_identities))
        print(bsr_AA_identities)
        # Dock fragments by aligning vdg backbones to bsr backbones
        for frag in os.listdir(vdgs_dir):
            if not frag.startswith('gu'):
                continue
            print(frag)
            vdgs_path = os.path.join(vdgs_dir, frag, 'nr_vdgs', str(
                len(bsr_combo)), bsr_AA_identities)
            if not os.path.exists(vdgs_path):
                print(f"Directory {vdgs_path} does not exist.")
                continue
            vdg_paths = os.listdir(vdgs_path)
            # Iterate over each vdg
            for vdg_path in vdg_paths:
                vdg_full_path = os.path.join(vdgs_path, vdg_path)
                vdg_prody_obj = pr.parsePDB(vdg_full_path)
                # Superpose the vdg backbone atoms onto the input structure.
                # Make sure order of query vdg AAs is consistent with order of input_pdb AAs
                resind_perms = map_aa_identities_to_vdg_resinds(
                    [vdg_prody_obj.select(f'resindex {_r}').getResnames()[0] 
                        for _r in sorted(set(vdg_prody_obj.getResindices()))], 
                    bsr_AA_identities.split('_'), )
                vdg_aminoacids = [vdg_prody_obj.select(f'resindex {_r}').getResnames()[0] 
                    for _r in sorted(set(vdg_prody_obj.getResindices()))]
                
                # Iterate over each resind permutation, get bb coords, and superpose onto 
                # the query structure (input structure)
                for resind_perm in resind_perms:
                    vdg_perm_bb_coords = []
                    for _r in resind_perm:
                        for atom_name in ['N', 'CA', 'C']: 
                            _coords = get_atom_coords(vdg_prody_obj.select(
                                f'resindex {_r}'), atom_name)
                            vdg_perm_bb_coords.append(_coords)
                    vdg_perm_bb_coords = np.array(vdg_perm_bb_coords)
                    # Superpose the vdg bb atoms onto the input structure and calculate rmsd
                    moved_coords, transf = pr.superpose(vdg_perm_bb_coords, 
                                                        input_bsr_bb_coords)
                    rmsd = pr.calcRMSD(moved_coords, input_bsr_bb_coords)
                    if rmsd < cutoff:
                        print(rmsd)
                        label = f'output/{bsr_combo}_{frag}_{vdg_path}_{resind_perm}.pdb'
                        pr.writePDB(label, pr.applyTransformation(transf, vdg_prody_obj))

def get_atom_coords(prody_obj, atom_name):
    atom = prody_obj.select(f'name {atom_name}')
    coords = atom.getCoords()[0]
    return coords

def get_vdg_subsets(input_list):
    # Initialize an empty list to store all subsets
    all_subsets = []
    # Loop through subset sizes 1 to 3, generate combos, then add to list
    for r in range(1, 4):
        subsets = combinations(input_list, r)
        all_subsets.extend(subsets)
    return all_subsets

def parse_yaml_input(yml_file):
    with open(yml_file, 'r') as inF:
        user_data = yaml.load(inF, Loader=yaml.FullLoader)
    vdgs_dir = user_data['vdgs_dir']
    input_pdb_path = user_data['input_pdb']
    bindingsite_residues = user_data['bindingsite_residues']
    return [vdgs_dir, input_pdb_path, bindingsite_residues]

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

if __name__ == "__main__":
    main()
