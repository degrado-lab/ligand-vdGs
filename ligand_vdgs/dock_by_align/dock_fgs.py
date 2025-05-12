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
    return parser.parse_args()

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

def main():
    print('WARNING: ONLY PROCESSING FRAGMENTS STARTING WITH "g"')
    # Parse inputs from yml file
    args = parse_args()
    yml_file = args.yml_file
    (vdgs_dir, input_pdb_path, bindingsite_residues) = parse_yaml_input(yml_file)
    # Parse the pdb and iterate through the combinations of binding site residues
    input_pdb = pr.parsePDB(input_pdb_path)
    bsr_combos = get_vdg_subsets(bindingsite_residues)
    for bsr_combo in bsr_combos:
        bsr_AA_identities = []
        input_bsr_bb_coords = []
        for bsr in bsr_combo: 
            bsr_seg, bsr_chain, bsr_resnum = bsr
            res_obj = select_residue(input_pdb, bsr_seg, bsr_chain, bsr_resnum)
            bb_coords = [get_atom_coords(res_obj, atom_name) for atom_name in ['N', 'CA', 'C']]
            input_bsr_bb_coords.append(bb_coords)
            AA = get_res_AA_identity(res_obj)
            bsr_AA_identities.append(AA)
        bsr_AA_identities = [AA if AA != 'GLY' else 'bb' for AA in bsr_AA_identities]
        bsr_AA_identities = '_'.join(sorted(bsr_AA_identities))
        print(bsr_AA_identities)
        # Dock fragments by aligning vdg backbones to bsr backbones
        for frag in os.listdir(vdgs_dir):
            if not frag.startswith('ge'):
                continue
            print(frag)
            vdgs_path = os.path.join(vdgs_dir, frag, 'nr_vdgs', str(
                len(bsr_combo)), bsr_AA_identities)
            if not os.path.exists(vdgs_path):
                print(f"Directory {vdgs_path} does not exist.")
                continue
            vdg_paths = os.listdir(vdgs_path)
            for vdg_path in vdg_paths:
                vdg_full_path = os.path.join(vdgs_path, vdg_path)
                vdg_prody_obj = pr.parsePDB(vdg_full_path)
                # Superpose the vdg backbone atoms onto the input structure
                vdg_bb_coords  = []
                # Make sure order of vdg AAs is consistent with order of input_pdb AAs
                map_aa_identities_to_vdg_aas(bsr_AA_identities.split('_'), vdg_prody_obj)
                #vdg_res_objs = []
                #equivalent = [] # to note when vdgs contain >1 of the same AA
                #for aa_identity in bsr_AA_identities.split('_'):
                #    vdg_res_obj = vdg_prody_obj.select(f'resname {aa_identity}')
                #    if len(set(vdg_res_obj.getResindices())) > 1:
                #        print(aa_identity)

def map_aa_identities_to_vdg_aas(AAs, vdg_prody_obj):
    # challenge is that if AA is 'bb', then it can be any of the AAs, so define that last
    vdg_residues = []
    not_bb = [a for a in AAs if a != 'bb']
    visited_resinds = []
    # first, take care of the non-bb AAs
    for AA in not_bb:
        # find the vdg residues that match this AA
        vdg_AA_resinds = set(vdg_prody_obj.select(f'resname {AA}').getResindices()) - set(visited_resinds)
        # add the first resind if there are >1 vdg residues that map to this AA. they'll 
        # be permuted at a later step.
        vdg_residues.append(tuple([vdg_prody_obj.select(f'resindex {list(vdg_AA_resinds)[0]}'), 
                                   AA]))
        visited_resinds.append(list(vdg_AA_resinds)[0])
    print(visited_resinds)
    for a, b in zip(vdg_prody_obj.getResindices(), vdg_prody_obj.getResnames()):
        print(a, b)
    # this doesn't make sense bc when subset is size 1 and just "THR", needs to sample both resinds 0 and 1



    #print(visited_resinds, set(vdg_prody_obj.getResindices()))
    ## then, the remaining AAs are 'bb'
    #remaining_resinds = set(vdg_prody_obj.getResindices()) - set(visited_resinds)






def get_atom_coords(prody_obj, atom_name):
    atom = prody_obj.select(f'name {atom_name}')
    assert len(atom) == 1
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

if __name__ == "__main__":
    main()
