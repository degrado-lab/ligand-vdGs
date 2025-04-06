'''
Generates fragments from all PDB ligands. Fragments are defined as all atoms that are 
within 3 bonds of any ligand atom. This script takes a few minutes to run.

Inputs: 
    CCD_smiles: tab-delimited file containing the smiles of molecules in the Chemical 
        Component Dictionary (CCD), downloaded from: https://www.wwpdb.org/data/ccd. This 
        is provided in ligand-vdGs under resources/Components-smiles-stereo-oe.smi 
    
Usage: 
    >> python ligand_vdgs/scripts/fragment_pdb_ligs.py
'''

import os
import sys
import argparse
import re
import csv
import pickle as pkl
from rdkit import Chem, RDLogger
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import fragment as frag

RDLogger.DisableLog('rdApp.*') 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ccd', default='resources/Components-smiles-stereo-oe.smi',
                        help="Path to CCD file containing smiles of molecules. Defaults "
                             "to resources/Components-smiles-stereo-oe.smi.")
    parser.add_argument('--outdir', default='resources', 
                        help="Path to output a pickled dict of database fragments and sdf "
                        "file of frags. Defaults to ligand-vdGs/resources.")
    parser.add_argument('--min-counts-per-frag', default=50, type=int, 
                        help="Minimum number of unique ccd ligands each fragment must be a "
                        "substruct of to qualify for being stored in the output pkl and sdf "
                        "files; defaults to 50.")
    parser.add_argument('--logfile', default='fragment_pdb_ligs.log',
                        help="Path to output log file.")
    parser.add_argument('--bond-radius', default=2, type=int,
                        help="Number of bonds away from an atom to cut during "
                        "fragmentation. Bonds to H's count. Defaults to 2.") 
    return parser.parse_args()

def main():
    args = parse_args()
    ccd = args.ccd
    outdir = args.outdir
    min_counts_per_frag = args.min_counts_per_frag
    logfile = args.logfile
    bond_radius = args.bond_radius

    # ------------------------------------------------------------------------------------------
    # Ensure that the input ccd file is valid
    if not os.path.isfile(ccd):
        raise FileNotFoundError(f"CCD file {ccd} does not exist.")
    
    # Set up output dir
    out_dict_path = os.path.join(args.outdir, 'database_frags_dict.pkl')
    out_sdf_path = os.path.join(args.outdir, 'database_frags.sdf')

    if not os.path.exists(outdir):
        os.makedirs(outdir) 

    for output_file in [out_dict_path, out_sdf_path]:
        if os.path.exists(out_dict_path):
            raise FileExistsError(f'Output file {output_file} already exists. Exiting.')

    '''
    Use frag_dict to store fragments and the molecules from which they come from.
        Dict Key: An alphabetically-sorted list of elements in the substruct (fragment)
                  to group isomers.
        Dict Value: Subdicts where: 
            Subdict Key: smiles of fragments with the same composition as the Dict Key 
            Subdict Value: num of unique ligands in the ccd that contain that fragment.
    
    The Dict key is essentially a representation of its formula, and is used because when 
    each new molecule frag is being added, it can be expensive to check whether the frag is 
    already represented in the dict, as opposed to checking only its isomers.
    '''

    frag_dict = {} 

    # Load file of all ligands in the PDB
    num_ligs_skipped = 0
    num_ligs_passed = 0
    num_ligs_total = 0
    with open(ccd, mode="r", newline="") as file:
        reader = csv.reader(file, delimiter="\t")
        # Process each ligand
        for line_ind, row in enumerate(reader):
            num_ligs_total += 1
            assert len(row) == 3, f"Row {row} expected to have exactly 3 cols."
            smiles, lig_resname, chem_name = row

            if undesired_elements(smiles): # filter out non-druglike elements
                num_ligs_skipped += 1
                continue

            # Convert SMILES to RDKit molecule and remove H's
            orig_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            mol, smiles_no_Hs = manually_remove_Hs(orig_mol) # rdkit's remove H methods 
                                                             # don't work 

            # Make sure there's at least a C, O, or N (rule out Fe-S clusters, ions, etc.) 
            elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
            if not any(elem in elements for elem in ['C', 'O', 'N']):
                num_ligs_skipped += 1
                continue

            # Decompose the ligand into fragments and store the fragment SMILES. Use SMILES 
            # instead of SMARTS b/c only SMILES (from rdkit) differentiates aliphatic and 
            # aryl.
            substructs = frag.fragment_on_bond_d(mol, bond_radius)
            # Update frag_dict by adding the new fragments; skip if all-carbon
            for orig_sub in substructs: # contains H's that need to be scrubbed
                sub, sub_smiles = manually_remove_Hs(orig_sub) 
                frag_dict = record_frag(sub, frag_dict, sub_smiles)

            num_ligs_passed += 1

    output_results(frag_dict, out_sdf_path, out_dict_path, min_counts_per_frag)
    report_stats(frag_dict, num_ligs_skipped, num_ligs_total, num_ligs_passed)

def output_results(frag_dict, out_sdf, out_pkl, min_counts_per_frag):
    # Determine which fragments to write out to an SDF file (based on the 
    # min_counts_per_frag) and order by size
    sorted_mols_to_output, sorted_labels = sort_frag_dict(frag_dict, min_counts_per_frag)

    # Write out fragments to an SDF file .
    frag.output_sdf(sorted_mols_to_output, sorted_labels, out_sdf) 
    
    # Write out frag_dict as a pickle file. 
    # ### should be database_frags_dict.pkl 
    pkl.dump(frag_dict, open(out_pkl, 'wb'))


def sort_frag_dict(frag_dict, min_counts_per_frag):
    # Create a list of (mol_object, smiles, count) tuples
    fragment_data = []

    for fragments in frag_dict.values():  # Iterate over elements
        for smiles, count in fragments.items():
            if count < min_counts_per_frag:
                continue
            mol_from_smiles = Chem.MolFromSmarts(smiles) # treat as smarts
            fragment_data.append((mol_from_smiles, f'{smiles}_{count}', count))

    # Sort by count (most to least)
    fragment_data.sort(key=lambda x: x[2], reverse=True)

    # Unpack sorted data
    mols_to_output, labels, _ = zip(*fragment_data) if fragment_data else ([], [], [])

    return mols_to_output, labels
    

def report_stats(frag_dict, num_ligs_skipped, num_ligs_total, num_ligs_passed):
    # Report statistics    
    final_num_frags  = 0 # number of unique frags
    for elements, frags in frag_dict.items():
        final_num_frags += len(frags)

    print('='*30, ' Report ', '='*30)
    print(f'Number of unique fragments recorded in frag_dict: {final_num_frags}')
    print(f'Number of mols passed: {num_ligs_passed}')
    print(f'Number of mols unsuccessfully parsed: {num_ligs_skipped} out of {num_ligs_total}') 

def manually_remove_Hs(orig_mol):
    # Remove hydrogens. Docs say that Chem.RemoveHs() implicit and explicit are removed, 
    # but this isn't true, so need to manually remove H's. 
    mol = Chem.RemoveHs(orig_mol, sanitize=False)
    # Convert back to SMILES, without showing explicit hydrogens
    smiles = Chem.MolToSmiles(mol, allHsExplicit=False, isomericSmiles=False) # still has H's
    # Remove standalone [H] completely
    smiles_no_Hs = re.sub(r'\[H\]', '', smiles)
    # Remove all 'H' except when part of '[Hg]'
    smiles_no_Hs = re.sub(r'H(?!g\])', '', smiles_no_Hs)

    return mol, smiles_no_Hs

def record_frag(substruct, frag_dict, smiles):

    '''
    Update frag_dict to store these substructs and the molecule from which they are from.
    Dict Key: An alphabetically-sorted list of elements in the substruct (representation of a 
              formula to group substructs into isomers)
    Dict Value: Subdicts where: 
        Subdict Key: smiles of fragments with the same composition as the Dict Key 
        Subdict Value: num of unique ligands in the ccd that contain that fragment.
    
    The Dict key is essentially a representation of its formula, and is used because when 
    each new molecule frag is being added, it can be expensive to check whether the frag is 
    already represented in the dict, as opposed to checking only its isomers.
    
    Note: Skip fragment if all its atoms are carbons. 
    '''

    if 'Mg' in smiles:
        return frag_dict
    
    elements = "".join(sorted([a.GetSymbol() for a in substruct.GetAtoms() if 
                               a.GetSymbol() != 'H'])) # b/c rdkit will add implicit H's
    
    # Skip if the elements are all carbons
    num_carbons = len([i for i in elements if i == 'C' or i == 'c'])
    if len(elements) == num_carbons:
        return frag_dict

    
    if elements not in frag_dict.keys():
        frag_dict[elements] = {}
        frag_dict[elements][smiles] = 1 
        return frag_dict

    if smiles in frag_dict[elements]: # Is the smiles already recorded?
        frag_dict[elements][smiles] += 1 
        return frag_dict # automatically duplicate

    else: 
        # Even if the smiles isn't already recorded, it still might be represented b/c the 
        # smiles could be degenerate, even w/ canonical. Use HasSubstructMatch() to be sure. 
        for existing_smiles in frag_dict[elements].keys():
            try:
                existing_sub_mol = Chem.MolFromSmarts(existing_smiles) # treat as smarts
                if existing_sub_mol.HasSubstructMatch(substruct) and \
                    substruct.HasSubstructMatch(existing_sub_mol): # then it's a duplicate
                    frag_dict[elements][existing_smiles] += 1 
                    return frag_dict
            except:
                return frag_dict

    # Passed all checks and deemed nonredundant.
    frag_dict[elements][smiles] = 1
    return frag_dict

def undesired_elements(smiles): 
    if re.search(r'(Be|Pt|Ru|Ir|Fe|Zn|Mg|Cu|Sn|Zr|Pb|Ga|Hg|Pd|Si|Mo|W|Se|Mn|Ti|Y|V|Ni|Rh|Te|Au|Ag|Co|Sb|Tb|Re|As|Cd|Hf|Lu|Na|Ca|Os|Cr|In|Al)', smiles):
        return True
    if re.search(r'B(?!r)', smiles): # bypass if it's Br, but return as undesired if it's B w/o r
        return True
    return False

if __name__ == "__main__":
    main()
