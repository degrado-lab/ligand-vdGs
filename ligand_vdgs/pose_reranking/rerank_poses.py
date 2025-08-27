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
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import match_vdgs as match
import Frags
import dock
import tempfile
from itertools import product

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
    hetatms = query_struct.hetatm.select('not (ion or resname SEP or resname TPO or resname MSE)')
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
        pdb_mol_assigned_bonds_no_H, pdb_mol_smiles_no_H = Frags.manually_remove_Hs(pdb_mol_assigned_bonds)
        filtered_frags = Frags.get_fragments(2, pdb_mol_assigned_bonds_no_H, 4, 5)
        return filtered_frags, pdb_mol_assigned_bonds_no_H, query_struct
    except ValueError as e:
        print(f"Error processing {pdbfile}: {e}", flush=True)
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
            pos = conf.GetAtomPosition(atom_idx)  # returns an RDKit Point3D object
            xyz = (pos.x, pos.y, pos.z)
            perm_coords.append(xyz)
        q_cg_coord_perms.append(np.array(perm_coords))
    return q_cg_coord_perms

def make_outdir(pdbfile, outdir):
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

def check_vdg_job_status(sub_smiles):
    # Check if the vdg generation job finished without issues
    vdg_log_file = os.path.join(vdg_lib_dir, sub_smiles, 'logs', f'{sub_smiles}_log')
    if not os.path.exists(vdg_log_file):
        return False
    with open(vdg_log_file, 'r') as f:
        log_contents = f.read()
    return 'Job completed.' in log_contents

def summarize_frags(deduplicated_filtered_frags, frags_in_lib, frags_to_exclude, 
                    frags_to_include):
    groups = {
        "Excluded": [],
        "Not in include list": [],
        "Not searched": [],
        "Not in vdg lib or incomplete": [],
        "In vdg lib": []}

    for sub, sub_smiles in deduplicated_filtered_frags:
        if any(smiles_equiv(f, sub_smiles) for f in frags_to_exclude):
            groups["Excluded"].append(sub_smiles)
        elif frags_to_include != 'all' and isinstance(frags_to_include, list) and sub_smiles not in frags_to_include:
            groups["Not in include list"].append(sub_smiles)
        elif sub_smiles not in frags_in_lib:
            groups["Not searched"].append(sub_smiles)
        elif not frags_in_lib[sub_smiles]:
            groups["Not in vdg lib or incomplete"].append(sub_smiles)
        else:
            groups["In vdg lib"].append(sub_smiles)

    print("\n=== Fragment summary ===", flush=True)
    for category, frags in groups.items():
        if frags:
            print(f"{category}:", flush=True)
            print("   ", ", ".join(frags), flush=True)

def check_in_exclude_list(sub_smiles, frags_to_exclude):
    for frag in frags_to_exclude:
        if smiles_equiv(frag, sub_smiles):
            return True
    return False

def main():
    
    # Load ground truth structure
    solved_struct = pr.parsePDB(solved_struct_path)

    # Get binding site residues from the solved structure
    ligname = set(solved_struct.hetatm.select('not (ion or resname SEP or resname TPO or resname MSE)').getResnames())
    assert len(ligname) == 1
    ligname = list(ligname)[0]

    # Gather all binding site residue combinations
    all_bsr_combos = get_bsr_combinations(solved_struct, ligname)
    # For each prediction, get the frags and CG coords in the pred, iterate through 
    # binding site res combos, align, and determine if there are vdg matches. 
    
    # Keep track of which frags are in vdg lib and which aren't
    frags_in_lib = {}

    # Keep track of which database vdg files could not be loaded
    failed_vdg_files = []

    # Load predictions
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
        output_dir = make_outdir(pdbfile, outdir)

        # Get all fragments. For every frag, iterate over the bsr combos and determine 
        # vdg matches.
        orig_filtered_frags, pdb_mol_reassigned, query_struct = get_frags_from_pdbfile(
            query_pdbs_dir, pdbfile)
        # Remove duplicates of filtered_frags 
        deduplicated_filtered_frags = []
        for (sub, sub_smiles) in orig_filtered_frags:
            # Step 0) is this frag named something else in the vdg lib?
            for existingfrag in os.listdir(vdg_lib_dir):
                if smiles_equiv(existingfrag, sub_smiles):
                    sub_smiles = existingfrag
                    break

            # First of all, is it in frags_to_exclude?
            if check_in_exclude_list(sub_smiles, frags_to_exclude): 
                # means it's equiv to a frag in frags_to_exclude
                continue

            # If sub_smiles already exists, then skip. Otherwise, add
            if sub_smiles in [i[1] for i in deduplicated_filtered_frags]:
                continue
            if sub_smiles not in [i[1] for i in deduplicated_filtered_frags]:
                deduplicated_filtered_frags.append((sub, sub_smiles))
            if sub_smiles in frags_to_exclude: # Exclude specified frags 
                continue
            
            # Skip if it's in frags_to_exclude. Use a separate var (_degenerate) 
            # to avoid confusion.
            if check_in_exclude_list(sub_smiles, frags_to_exclude):
                continue
        
            # Determine if this frag is in vdg lib (and the job didn't break)
            if sub_smiles not in frags_in_lib.keys():
                # Is it in the vdg lib dir?
                if sub_smiles in os.listdir(vdg_lib_dir):
                    # Did the job (to create the vdgs) finish w/o issue?
                    if check_vdg_job_status(sub_smiles):
                        frags_in_lib[sub_smiles] = True
                    else:
                        frags_in_lib[sub_smiles] = False 
                else: 
                    # If the str isn't in os.listdir(), it's still possible that a 
                    # degenerate smiles is in vdg_lib_dir
                    for existingfrag in os.listdir(vdg_lib_dir):
                        # Is it equivalent to sub_smiles and was its job completed?
                        is_equivalent = smiles_equiv(existingfrag, sub_smiles)
                        if is_equivalent: # rename the smiles to match what's in db
                            sub_smiles = existingfrag
                        if is_equivalent and check_vdg_job_status(existingfrag):
                            frags_in_lib[sub_smiles] = True
                            break
                    else:
                        frags_in_lib[sub_smiles] = False
        # Log
        summarize_frags(deduplicated_filtered_frags, frags_in_lib, frags_to_exclude, frags_to_include)

        # Iterate over frags that are in the vdg lib
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
                vdg_paths = [p for p in os.listdir(vdgs_path) if p.endswith(('.pdb', '.pdb.gz'))]
                # Iterate over each vdg
                for vdg_path in vdg_paths:
                    vdg_full_path = os.path.join(vdgs_path, vdg_path)
                    try:
                        vdg_prody_obj = pr.parsePDB(vdg_full_path)
                    except Exception as e:
                        failed_vdg_files.append(vdg_full_path)
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
                        if not os.path.exists(subdir_path):
                            os.makedirs(subdir_path)
                        output_vdg_name = f'{pdbfile.removesuffix(".pdb")}_{sub_smiles}_' + vdg_path.rstrip('/').split(
                            '/')[-1].removesuffix('.pdb.gz') + f'_{round(db_to_query_rmsd, 2):0.2f}'
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
