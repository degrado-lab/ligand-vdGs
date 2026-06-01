'''
Generates fragments from all database ligands. Fragments are defined as all atoms that are within
`--bond-radius` bonds (default: 2) of any ligand atom. This script takes 5-10 mins to run on
48k smiles.

Inputs:
    CCD_smiles: tab-delimited file with 2-3 columns:
                    col 1 (required): SMILES string
                    col 2 (required): short ligand identifier (e.g. CCD 3-letter code or
                                      internal compound ID)
                    col 3 (optional): ligand name; not used by this script

                Defaults are set up for PDB ligands: the CCD file is provided in
                ligand-vdGs under resources/Components-smiles-cactvs.smi, downloaded
                from https://www.wwpdb.org/data/ccd.

Usage (PDB defaults):
    >> python ligand_vdgs/generate_vdgs/fragment_database_ligs.py

Usage (custom ligand list):
    Specify --ccd and --outdir; --logfile defaults to logs/fragment_database_ligs.log.
    >> python ligand_vdgs/generate_vdgs/fragment_database_ligs.py \
           --ccd <path/to/your_ligands> \
           --outdir <output_dir>

Output (<outdir>/database_frags_dict.pkl):
    A pickled nested dict with structure:
        {
            elements_str (str): {
                smiles (str): [lig_resname (str), ...]
            }
        }
    where:
        elements_str  -- alphabetically sorted concatenation of non-H element symbols
                         for all atoms in the fragment (e.g. "CCNO"); used to group
                         structural isomers so duplicate-checking only scans isomers.
        smiles        -- canonical SMILES of the fragment.
        lig_resname   -- identifier of each database ligand (e.g. CCD 3-letter code for PDB ligands).
    Within each elements_str bucket, fragments are sorted by number of unique ligands
    (descending).
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
                             "to resources/Components-smiles-cactvs.smi.")
    parser.add_argument('--outdir', default='resources',
                        help="Path to output a pickled dict of database fragments. "
                        "Defaults to ligand-vdGs/resources.")
    parser.add_argument('--logfile', default='logs/fragment_database_ligs.log',
                        help="Path to output log file. Defaults to logs/fragment_database_ligs.log.")
    parser.add_argument('--bond-radius', default=2, type=int,
                        help="Number of bonds away from an atom to cut during "
                        "fragmentation. Defaults to 2.") 
    return parser.parse_args()

def main():
    args = parse_args()
    ccd = args.ccd
    outdir = args.outdir
    logfile = args.logfile
    bond_radius = args.bond_radius

    # ------------------------------------------------------------------------------------------
    # Ensure that the input ccd file is valid
    if not os.path.isfile(ccd):
        raise FileNotFoundError(f"CCD file {ccd} does not exist.")
    
    # Set up output dir
    out_dict_path = os.path.join(outdir, 'database_frags_dict.pkl')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    log_dir = os.path.dirname(logfile)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    if os.path.exists(out_dict_path):
        raise FileExistsError(f'Output file {out_dict_path} already exists. Exiting.')

    # frag_dict: {elements_str: {smiles: [lig_resnames]}}
    # Keyed by alphabetically-sorted element string to group isomers together; duplicate-
    # checking via HasSubstructMatch() is expensive for large numbers of fragments, so we 
    # only run it within isomer groups.
    frag_dict = {}

    # Load file of ligands
    num_ligs_failed = 0
    num_ligs_not_druglike = 0
    num_ligs_passed = 0
    num_ligs_total = 0
    failed_smiles = []
    with open(ccd, mode="r", newline="") as file:
        reader = csv.reader(file, delimiter="\t")
        # Process each ligand
        for line_num, row in enumerate(reader, start=1):
            num_ligs_total += 1
            assert 2 <= len(row) <= 3, f"Row {row} expected to have 2 or 3 cols."
            smiles, lig_resname = row[0], row[1]

            if undesired_elements(smiles): # filter out non-druglike elements
                num_ligs_not_druglike += 1
                continue

            # Convert SMILES to RDKit molecule and remove H's
            orig_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if orig_mol is None:
                print(f'[WARNING] line {line_num} ({lig_resname}): could not parse SMILES: {smiles}')
                failed_smiles.append((line_num, lig_resname, smiles))
                num_ligs_failed += 1
                continue

            results = Frags.manually_remove_Hs(orig_mol, 'single') # rdkit's remove
                                                    # H method isn't good enough.

            if results is None:
                print(f'[WARNING] line {line_num} ({lig_resname}): could not parse SMILES: {smiles}')
                failed_smiles.append((line_num, lig_resname, smiles))
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
            # aryl (C,c vs. [#6]). Fragment on bond radii `bond_radius` AND the positive
            # integers less than `bond_radius`, because for example, drugs containing 
            # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O. 
            filtered_frags = Frags.get_fragments(bond_radius, mol, 4, 5)
            for sub_smiles, substruct_site_groups in filtered_frags.items():
                substruct_site = substruct_site_groups[0] # only need one site group 
                substruct_perm = substruct_site[0] # only need one perm.
                sub = substruct_perm[0] # select Mol obj (Mol obj, perm inds, orig mol inds)
                frag_dict = record_frag(sub, frag_dict, sub_smiles, lig_resname)

            num_ligs_passed += 1

    output_results(frag_dict, out_dict_path)
    report_stats(frag_dict, num_ligs_failed, num_ligs_total, num_ligs_passed,
                 num_ligs_not_druglike, failed_smiles, logfile)

def output_results(frag_dict, out_pkl):
    pkl.dump(sort_frag_dict(frag_dict), open(out_pkl, 'wb'))


def sort_frag_dict(frag_dict):
    sorted_dict = {}
    for elements, fragments in frag_dict.items():
        sorted_dict[elements] = dict(
            sorted(fragments.items(), key=lambda x: len(set(x[1])), reverse=True)
        )
    return sorted_dict

def report_stats(frag_dict, num_ligs_failed, num_ligs_total, num_ligs_passed,
                 num_ligs_not_druglike, failed_smiles, logfile):
    final_num_frags = 0
    for _, frags in frag_dict.items():
        final_num_frags += len(frags)
    with open(logfile, 'w') as log:
        log.write('='*30 + ' Report ' + '='*30 + '\n')
        log.write(f'Number of mols in input file: {num_ligs_total}\n')
        log.write(f'Number of unique fragments recorded: {final_num_frags}\n')
        log.write(f'Number of mols parsed: {num_ligs_passed}\n')
        log.write(f'Number of mols unsuccessfully parsed: {num_ligs_failed}\n')
        log.write(f'Number of mols skipped b/c not druglike: {num_ligs_not_druglike}\n')
        if failed_smiles:
            log.write('\nFailed SMILES:\n')
            for line_num, lig_resname, smi in failed_smiles:
                log.write(f'  line {line_num} ({lig_resname}): {smi}\n')

def record_frag(substruct, frag_dict, smiles, lig_resname):

    '''
    Update frag_dict to store these substructs and the molecule from which they are from.
    Dict Key: Alphabetically-sorted concatenation of element symbols in the substruct
              (elements_str; used to group structural isomers together)
    Dict Value: Subdicts where: 
        Subdict Key: smiles of fragments with the same composition as the Dict Key
        Subdict Value: list of CCD resnames of ligands that contain that fragment.

    The Dict key is essentially a representation of its formula, and is used because when
    each new molecule frag is being added, it can be expensive to check whether the frag is
    already represented in the dict, as opposed to checking only its isomers.

    Note: All-carbon fragments are excluded upstream by Frags.get_fragments() before
          reaching this function.
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
            except Exception:
                return frag_dict

    # Passed all checks and deemed nonredundant.
    frag_dict[elements][smiles] = [lig_resname]
    return frag_dict

def undesired_elements(smiles):
    metals     = r'Be|Pt|Ru|Ir|Fe|Zn|Mg|Cu|Sn|Zr|Pb|Ga|Hg|Pd|Si|Mo|W|Se|Mn|Ti|Y|V|Ni|Rh|Te|Au|Ag|Co|Sb|Re|As|Cd|Hf|Na|Ca|Os|Cr|In|Al|U|K|se'
    lanthanides = r'La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu'
    noble_gases = r'He|Ne|Ar|Kr|Xe|Rn'
    if re.search(f'({metals}|{lanthanides}|{noble_gases})', smiles):
        return True
    if re.search(r'B(?!r)', smiles): # bypass Br, flag bare B
        return True
    return False

if __name__ == "__main__":
    main()
