import os
import yaml
import sys
import argparse
from itertools import product
import numpy as np
import prody as pr
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import utils
import dock_utils as dock

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', "--yml-file", type=str, required=True,
                        help="Yaml file containing binding pocket information.")
    return parser.parse_args()

def main():
    # Parse inputs from yml file
    args = parse_args()
    yml_file = args.yml_file
    (vdgs_dir, input_pdb_path, addl_bsr_res, out_dir, cutoff, frags_to_dock,
     ligname) = parse_yaml_input(yml_file)
    utils.set_up_outdir(out_dir, overwrite=False)
    # Parse the pdb and iterate through the combinations of binding site residues
    input_pdb = pr.parsePDB(input_pdb_path)
    bindingsite_residues = dock.get_bindingsite_residues(input_pdb,
        set([tuple(b) for b in addl_bsr_res]), ligname)
    bsr_combos = dock.get_vdg_subsets(bindingsite_residues)
    print('rrr')
    print(bsr_combos)
    # Gather all binding site residue combinations
    all_bsr_combos = []
    for bsr_combo in bsr_combos:
        print(bsr_combo, 'zz')
        bsr_AA_identities = []
        input_bsr_bb_coords = []
        for bsr in bsr_combo: 
            bsr_seg, bsr_chain, bsr_resnum = bsr
            res_obj = dock.select_residue(input_pdb, bsr_seg, bsr_chain, bsr_resnum)
            res_bb_coords = []
            for atom_name in ['N', 'CA', 'C']:
                res_bb_coords.append(dock.get_atom_coords(res_obj, atom_name))
            input_bsr_bb_coords.append(res_bb_coords)
            AA = dock.get_res_AA_identity(res_obj)
            bsr_AA_identities.append(AA)
        bsr_AA_identities = [AA if AA != 'GLY' else 'bb' for AA in bsr_AA_identities]
        # All of the AAs have potential to provide bb contacts, so sample "bb" vdms
        options = [(x,) if x == 'bb' else (x, 'bb') for x in bsr_AA_identities]
        combo_variants = list(product(*options))
        for c in combo_variants:
            if c not in all_bsr_combos:
                all_bsr_combos.append((c, input_bsr_bb_coords.copy())) 

    for combo, coords in all_bsr_combos:
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
        for frag in os.listdir(vdgs_dir):
            if not frag in frags_to_dock: 
                continue
            vdgs_path = os.path.join(vdgs_dir, frag, 'nr_vdgs', str(
                len(bsr_incl_bb_identities.split('_'))), bsr_incl_bb_identities)
            if not os.path.exists(vdgs_path):
                continue
            vdg_paths = os.listdir(vdgs_path)
            # Iterate over each vdg
            for vdg_path in vdg_paths:
                input_pdbname = os.path.basename(input_pdb_path)[:4]
                vdg_pdbname = vdg_path[:4]
                if input_pdbname == vdg_pdbname:
                    continue
                vdg_full_path = os.path.join(vdgs_path, vdg_path)
                vdg_prody_obj = pr.parsePDB(vdg_full_path)
                # Superpose the vdg backbone atoms onto the input structure.
                # Make sure order of query vdg AAs is consistent with order of input_pdb AAs
                resind_perms = map_aa_identities_to_vdg_resinds(
                    [vdg_prody_obj.select(f'resindex {_r}').getResnames()[0] 
                        for _r in sorted(set(vdg_prody_obj.getResindices()))], 
                    bsr_incl_bb_identities.split('_'), )
                
                # Iterate over each resind permutation, get bb coords, and superpose onto 
                # the query structure (input structure)
                for resind_perm in resind_perms:
                    vdg_perm_bb_coords = []
                    for _r in resind_perm:
                        for atom_name in ['N', 'CA', 'C']: 
                            _coords = dock.get_atom_coords(vdg_prody_obj.select(
                                f'resindex {_r}'), atom_name)
                            vdg_perm_bb_coords.append(_coords)
                    vdg_perm_bb_coords = np.array(vdg_perm_bb_coords)
                    # Superpose the vdg bb atoms onto the input structure and calculate rmsd
                    moved_coords, transf = pr.superpose(vdg_perm_bb_coords, 
                                                        input_bsr_bb_coords)
                    rmsd = pr.calcRMSD(moved_coords, input_bsr_bb_coords)
                    if rmsd > cutoff: 
                        continue

                    aligned_vdg_obj = pr.applyTransformation(transf, vdg_prody_obj)
                    # Check for clashes w/ protein backbone
                    assert len(set(aligned_vdg_obj.hetatm.getResindices())) == 1 # make sure 
                                                                    # selecting only 1 lig
                    # If there are atoms within 3.5A, they may be suspicious, so check them
                    neighbs = pr.findNeighbors(aligned_vdg_obj.select('occupancy > 2.9'),
                        3.5, input_pdb.select('name N CA C O and not element CA')) 
                    found_clash = False
                    for n in neighbs:
                        elem1, elem2, dist = n[0].getElement(), n[1].getElement(), n[2]
                        if dist < utils.min_clash_dist(elem1, elem2):
                            found_clash = True
                            break
                    if found_clash:
                        continue # continue to next vdg res permutation
                    # Write out pdb
                    vdgpdb = vdg_path.removesuffix('.pdb.gz')
                    label = f'{out_dir}/{frag}/{frag}_{vdgpdb}_' \
                            f'{"_".join([str(r) for r in resind_perm])}_{rmsd:.2f}'
                    # Make base dirs if they don't exist
                    if not os.path.exists(os.path.join(out_dir, frag)):
                        os.makedirs(os.path.join(out_dir, frag))
                    pr.writePDB(label, aligned_vdg_obj)
    print('DONE')

def parse_yaml_input(yml_file):
    with open(yml_file, 'r') as inF:
        user_data = yaml.load(inF, Loader=yaml.FullLoader)
    vdgs_dir = user_data['vdgs_dir']
    input_pdb_path = user_data['input_pdb']
    addl_bindingsite_residues = user_data['additional_bindingsite_residues']
    if not addl_bindingsite_residues:
        addl_bindingsite_residues = []
    out_dir = user_data['out_dir']
    cutoff = user_data['rmsd_threshold']
    frags_to_dock = user_data['frags_to_dock']
    ligname = user_data['ligname']
    return [vdgs_dir, input_pdb_path, addl_bindingsite_residues, out_dir, cutoff, 
            frags_to_dock, ligname]

if __name__ == "__main__":
    main()
