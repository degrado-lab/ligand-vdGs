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
from rdkit import Chem
from rdkit.Chem import AllChem
import re
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import match_vdgs as match
import utils
import Frags
import dock
import tempfile
from itertools import product

# TODO: DETECT AROMATICITY INSTEAD OF SUPPLYING LIG SMILES

with open(sys.argv[1], 'r') as f:
    config = yaml.safe_load(f)

# Config variables
rmsd_threshold = config['rmsd_threshold']
query_pdbs_dir = config['query_pdbs_dir']
solved_struct_path = config['solved_struct_path']
lig_smiles = config['lig_smiles']
vdg_lib_dir = config['vdg_lib_dir']
outdir = config['outdir']
query_lig_res = config['query_lig_res']  # (segment, chain, resnum
frags_to_exclude = config['frags_to_exclude']  # list of SM


def get_bsr_combinations(solved_struct, ligname):
    bindingsite_residues = dock.get_bindingsite_residues(solved_struct,
        [], ligname, dist_from_lig=4.5, CA_only=False)
    bsr_combos = dock.get_vdg_subsets(bindingsite_residues)
    all_bsr_combos = []
    for bsr_combo in bsr_combos:
        bsr_AA_identities = []
        input_bsr_bb_coords = []
        if len(bsr_combo) != 2:
            continue
        for bsr in bsr_combo: 
            bsr_seg, bsr_chain, bsr_resnum = bsr
            res_obj = dock.select_residue(solved_struct, bsr_seg, bsr_chain, bsr_resnum)
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
    return all_bsr_combos

def get_frags_from_pdbfile(query_pdbs_dir, pdbfile):
    query_struct = pr.parsePDB(os.path.join(query_pdbs_dir, pdbfile))
    # Identify ligand and fragment it
    hetatms = query_struct.hetatm
    assert len(set(hetatms.getResindices())) == 1
    # Convert from prody obj to rdkit Mol obj
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_pdb:
        pr.writePDB(tmp_pdb.name, hetatms)
    
    pdb_mol = Chem.MolFromPDBFile(tmp_pdb.name, removeHs=True)
    lig_template = Chem.MolFromSmiles(lig_smiles) 
    # Assign bond orders (b/c rdkit doesn't calculate this from PDB coords) to detect 
    # correct valence and aromaticity 
    try:
        pdb_mol_assigned_bonds = AllChem.AssignBondOrdersFromTemplate(lig_template, pdb_mol)
        pdb_mol_assigned_bonds_no_H, pdb_mol_smiles_no_H = Frags.manually_remove_Hs(pdb_mol_assigned_bonds, 
                                                                                    pdb_mol_assigned_bonds)
        filtered_frags = Frags.get_fragments(2, pdb_mol_assigned_bonds_no_H, 4, 5)
        return filtered_frags, pdb_mol_assigned_bonds_no_H, query_struct
    except ValueError as e:
        print(f"Error processing {pdbfile}: {e}")
        return [], None, query_struct

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
            #print(atom.GetSymbol())
            pos = conf.GetAtomPosition(atom_idx)  # returns an RDKit Point3D object
            info = atom.GetPDBResidueInfo()
            xyz = (pos.x, pos.y, pos.z)
            perm_coords.append(xyz)
        q_cg_coord_perms.append(np.array(perm_coords))
    return q_cg_coord_perms

def make_outdir(pdbfile, output_dir):
    # Name the outdir for this pdb query
    pdbname = pdbfile.removesuffix('.pdb')
    pdb_id = pdbname[:4]
    output_dir = os.path.join(outdir, pdb_id, pdbname)
    return output_dir

def smiles_equiv(existingfrag, sub_smiles):
    mol1 = Chem.MolFromSmarts(existingfrag)
    mol2 = Chem.MolFromSmarts(sub_smiles)
    is_equivalent = (mol1.HasSubstructMatch(mol2) and 
                 mol2.HasSubstructMatch(mol1) and 
                 mol1.GetNumAtoms() == mol2.GetNumAtoms())
    return is_equivalent

def main():
    
    # Load ground truth structure
    solved_struct = pr.parsePDB(solved_struct_path)

    # Get binding site residues from the solved structure
    ligname = set(solved_struct.hetatm.getResnames())
    assert len(ligname) == 1
    ligname = list(ligname)[0]

    # Gather all binding site residue combinations
    all_bsr_combos = get_bsr_combinations(solved_struct, ligname)
    # For each prediction, get the frags and CG coords in the pred, iterate through 
    # binding site res combos, align, and determine if there are vdg matches. 
    
    # Keep track of which frags are in vdg lib and which aren't
    frags_in_lib = {}
    
    # Load predictions
    for pdbfile in sorted(os.listdir(query_pdbs_dir)):
        print('\n' + '='*10, pdbfile, '='*10)
        output_dir = make_outdir(pdbfile, outdir)

        # Get all fragments. For every frag, iterate over the bsr combos and determine 
        # vdg matches.
        orig_filtered_frags, pdb_mol_reassigned, query_struct = get_frags_from_pdbfile(
            query_pdbs_dir, pdbfile)
        # Remove duplicates of filtered_frags 
        deduplicated_filtered_frags = []
        for (sub, sub_smiles) in orig_filtered_frags:
            # If sub_smiles already exists, then skip
            if sub_smiles in [i[1] for i in deduplicated_filtered_frags]:
                continue
            # Even if not, it may be represented with a diff smiles str
            degenerate = False
            for _s in [i[1] for i in deduplicated_filtered_frags]:
                if smiles_equiv(_s, sub_smiles):
                    degenerate = True
                    break
            if degenerate:
                continue
            # Otherwise, add it to deduplicated_filtered_frags
            if sub_smiles not in [i[1] for i in deduplicated_filtered_frags]:
                deduplicated_filtered_frags.append((sub, sub_smiles))
            # Determine if this frag is in vdg lib
            if sub_smiles in frags_to_exclude: # Exclude specified frags 
                continue
            _degenerate = False
            for _f in frags_to_exclude:
                if smiles_equiv(_f, sub_smiles):
                    degenerate = True
                    break
            if _degenerate:
                continue
            # Check if the sub_smiles is in the vdg lib
            if sub_smiles not in frags_in_lib.keys():
                if sub_smiles in os.listdir(vdg_lib_dir):
                    frags_in_lib[sub_smiles] = True
                else: 
                    # If the str isn't in os.listdir(), it's still possible that a 
                    # degernate smiles is in vdg_lib_dir
                    for existingfrag in os.listdir(vdg_lib_dir):
                        # Is it equivalent to sub_smiles?
                        is_equivalent = smiles_equiv(existingfrag, sub_smiles)
                        if is_equivalent:
                            frags_in_lib[sub_smiles] = True
                            break
                    else:
                        frags_in_lib[sub_smiles] = False
        # Iterate over frags that are in the vdg lib
        for sub, sub_smiles in deduplicated_filtered_frags:
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

            q_cg_coord_perms = get_query_cg_coords(pdb_mol_reassigned, sub)
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
                vdg_paths = os.listdir(vdgs_path)
                # Iterate over each vdg
                for vdg_path in vdg_paths:
                    vdg_full_path = os.path.join(vdgs_path, vdg_path)
                    try:
                        vdg_prody_obj = pr.parsePDB(vdg_full_path)
                    except Exception as e:
                        print(f"WARNING: could not parse, {vdg_full_path}")
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
                                _coords = dock.get_atom_coords(vdg_prody_obj.select(
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
                            db_to_query_rmsd = round(pr.calcRMSD(moved_db_vdg_bb_and_cg_coords, 
                                                           q_bb_and_cg), 2)
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
                        subdir_path = os.path.join(output_dir, out_subdir)
                        os.makedirs(subdir_path, exist_ok=True)
                        if not os.path.exists(subdir_path):
                            os.makedirs(subdir_path)
                        output_vdg_name = f'{sub_smiles}_' + vdg_path.rstrip('/').split(
                            '/')[-1].removesuffix('.pdb.gz') + f'_{db_to_query_rmsd}'
                        output_vdg_path = os.path.join(subdir_path, output_vdg_name+'.pdb')
                        # check if a file in output_vdg_path already exists
                        if os.path.exists(output_vdg_path):
                            print('WARNING: overwriting existing file', output_vdg_path)
                        pr.writePDB(output_vdg_path, pr.applyTransformation(db_transf, 
                            vdg_prody_obj.copy()))

    print('\nFragment SMILES found in the vdg lib:', 
                [smiles for smiles, in_lib in frags_in_lib.items() if in_lib])
    print('\nFragment SMILES not found in the vdg lib:', 
          [smiles for smiles, in_lib in frags_in_lib.items() if not in_lib])
    return

    


if __name__ == "__main__":
    main()
