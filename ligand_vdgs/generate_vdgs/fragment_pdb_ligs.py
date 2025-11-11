'''
Generates fragments from all PDB ligands. Fragments are defined as all atoms that are 
within 3 bonds of any ligand atom. This script takes 5-10 mins to run on 48k smiles.

Inputs: 
    CCD_smiles: tab-delimited file containing the smiles of molecules in the Chemical 
        Component Dictionary (CCD), downloaded from: https://www.wwpdb.org/data/ccd. This 
        is provided in ligand-vdGs under resources/Components-smiles-cactvs.smi 
    
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
import Frags 

RDLogger.DisableLog('rdApp.*') 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ccd', default='resources/Components-smiles-cactvs.smi',
                        help="Path to CCD file containing smiles of molecules. Defaults "
                             "to resources/CComponents-smiles-cactvs.smi.")
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
                        "fragmentation. Defaults to 2.") 
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
    num_ligs_failed = 0
    num_ligs_not_druglike = 0
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
                num_ligs_not_druglike += 1
                continue

            # Convert SMILES to RDKit molecule and remove H's
            orig_mol = Chem.MolFromSmiles(smiles, sanitize=False)

            results = Frags.manually_remove_Hs(orig_mol, 'single') # rdkit's remove 
                                                    # H method isn't good enough.

            if results is None:
                num_ligs_failed += 1
                continue
            mol_info, _ = results
            mol, _ = mol_info
            # Make sure there's at least a C, O, or N (rule out Fe-S clusters, ions, etc.)
            if not Frags.is_organic(mol):
                num_ligs_not_druglike += 1
                continue

            # Decompose the ligand into fragments and store the fragment SMILES. Use SMILES 
            # instead of SMARTS b/c only SMILES (from rdkit) differentiates aliphatic and 
            # aryl (C,c vs. [#6]). Fragment on bond radii `bond_radius` AND the postive 
            # integers less than `bond_radius`, because for example, drugs containing 
            # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O. 
            filtered_frags = Frags.get_fragments(bond_radius, mol, 4, 5)
            for sub_smiles, substruct_site_groups in filtered_frags.items():
                substruct_site = substruct_site_groups[0] # only need one site group 
                substruct_perm = substruct_site [0] # only need one perm. 
                sub = substruct_perm[0] # select Mol obj (Mol obj, perm inds, orig mol inds)
                frag_dict = record_frag(sub, frag_dict, sub_smiles, lig_resname)

            num_ligs_passed += 1

    output_results(frag_dict, out_sdf_path, out_dict_path, min_counts_per_frag)
    report_stats(frag_dict, num_ligs_failed, num_ligs_total, num_ligs_passed, 
                 num_ligs_not_druglike, logfile)

def output_results(frag_dict, out_sdf, out_pkl, min_counts_per_frag):
    # Determine which fragments to write out to an SDF file (based on the 
    # min_counts_per_frag) and order by size
    '''
    sorted_mols_to_output, sorted_labels = sort_frag_dict(frag_dict, min_counts_per_frag)
    # Write out fragments to an SDF file .
    frag.output_sdf(sorted_mols_to_output, sorted_labels, out_sdf) 
    '''
    
    # Write out frag_dict as a pickle file. 
    pkl.dump(frag_dict, open(out_pkl, 'wb'))


def sort_frag_dict(frag_dict, min_counts_per_frag):
    # Create a list of (mol_object, smiles, count) tuples
    fragment_data = []

    for fragments in frag_dict.values():  # Iterate over elements
        for smiles, lignames in fragments.items():
            count = len(set(lignames))
            if count < min_counts_per_frag:
                continue
            mol_from_smiles = Chem.MolFromSmarts(smiles) # treat as smarts
            fragment_data.append((mol_from_smiles, f'{smiles}_{count}', count))

    # Sort by count (most to least)
    fragment_data.sort(key=lambda x: x[2], reverse=True)

    # Unpack sorted data
    mols_to_output, labels, _ = zip(*fragment_data) if fragment_data else ([], [], [])

    return mols_to_output, labels

def report_stats(frag_dict, num_ligs_failed, num_ligs_total, num_ligs_passed, 
                 num_ligs_not_druglike, logfile):
    # Report statistics    
    final_num_frags  = 0 # number of unique frags
    for elements, frags in frag_dict.items():
        final_num_frags += len(frags)
    # Print to logfile
    with open(logfile, 'w') as log:
        log.write('='*30 + ' Report ' + '='*30 + '\n')
        log.write(f'Number of mols in input CCD: {num_ligs_total}\n')
        log.write(f'Number of unique fragments recorded in frag_dict: {final_num_frags}\n')
        log.write(f'Number of mols parsed: {num_ligs_passed}\n')
        log.write(f'Number of mols unsuccessfully parsed: {num_ligs_failed}\n')
        log.write(f'Number of mols skipped b/c not druglike: {num_ligs_not_druglike}\n')

def record_frag(substruct, frag_dict, smiles, lig_resname):

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

    elements = "".join(sorted([a.GetSymbol() for a in substruct.GetAtoms() if 
                               a.GetSymbol() != 'H'])) # b/c rdkit will add implicit H's

    # Determine whether to add fragment to the dict
    if elements not in frag_dict.keys(): # automatically a new entry 
        frag_dict[elements] = {}
        frag_dict[elements][smiles] = [lig_resname] 
        return frag_dict

    if smiles in frag_dict[elements].keys(): # Is the smiles already recorded?
        if lig_resname not in frag_dict[elements][smiles]:
            frag_dict[elements][smiles].append(lig_resname)
        return frag_dict 

    else: 
        # Even if the smiles isn't already recorded, it still might be represented b/c the 
        # smiles could be degenerate, even w/ canonical. Use HasSubstructMatch() to be sure. 
        for existing_smiles in frag_dict[elements].keys():
            try:
                existing_sub_mol = Chem.MolFromSmarts(existing_smiles) # treat as smarts
                # The smiles are equivalent only if they are both substructs of each other, 
                # b/c rdkit is often wrong. 
                mol1_has_mol2 = existing_sub_mol.HasSubstructMatch(substruct)
                mol2_has_mol1 = substruct.HasSubstructMatch(existing_sub_mol)

                if mol1_has_mol2 and mol2_has_mol1: # then it's a duplicate, 
                                                           # but still add bc diff lig name
                    if lig_resname not in frag_dict[elements][existing_smiles]:
                        frag_dict[elements][existing_smiles].append(lig_resname) 
                    return frag_dict
            except:
                return frag_dict

    # Passed all checks and deemed nonredundant.
    if elements not in frag_dict.keys():
        frag_dict[elements] = {}
    if smiles not in frag_dict[elements].keys():
        frag_dict[elements][smiles] = []
    if lig_resname not in frag_dict[elements][smiles]:
        frag_dict[elements][smiles].append(lig_resname)
    return frag_dict

def undesired_elements(smiles): 
    if re.search(r'(Be|Pt|Ru|Ir|Fe|Zn|Mg|Cu|Sn|Zr|Pb|Ga|Hg|Pd|Si|Mo|W|Se|Mn|Ti|Y|V|Ni|Rh|Te|Au|Ag|Co|Sb|Tb|Re|As|Cd|Hf|Lu|Na|Ca|Os|Cr|In|Al|Pr|se|U|K|Ar)', smiles):
    #if re.search(r'(Be|Pt|Ru|Ir|Fe|Zn|Mg|Cu|Sn|Zr|Pb|Ga|Hg|Pd|Si|Mo|W|Se|Mn|Ti|Y|V|Ni|Rh|Te|Au|Ag|Co|Sb|Tb|Re|As|Cd|Hf|Lu|Na|Ca|Os|Cr|In|Al|Pr|se|U|K)', smiles):
        return True
    if re.search(r'B(?!r)', smiles): # bypass if it's Br, but return as undesired if it's B w/o r
        return True
    return False

if __name__ == "__main__":
    main()
