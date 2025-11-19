'''
For a list of ligands, determine if a vdg library for each of their fragments exist, 
and identify fragments that don't currently have a vdg library.
'''

import sys
import os
from rdkit import Chem, RDLogger
sys.path.append(os.path.join(os.path.dirname(__file__), 'ligand_vdgs', 'functions'))
import Frags 
import utils

# ------- Example settings --------------------------------
list_smiles = ['CNCc1ccc(cc1)c2[nH]c3cc(F)cc4C(=O)NCCc2c34', 
               'O=C1OC2(C3=C(OC4=C2C=C(Br)C(O)=C4Br)C(Br)=C(O)C(Br)=C3)C5=C1C=CC=C5', 
               'C1CNC[C@@H]1OC2=C(C=C(C=C2NC(=O)C3=CC(=NC=N3)C(=O)NC4=CC(=CC(=C4O[C@@H]5CCNC5)NC(=O)CCCCN=C(N)N)C(F)(F)F)C(F)(F)F)NC(=O)CCCCN=C(N)N', 
               'C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O']

fraglibs = ['/wynton/home/degradolab/skt/docking/frag_lib', 
    '/wynton/group/degradolab/skt/docking/databases/frag_vdg_lib']
# ------------------------------------------------------------
frags_to_query = {}
frags_to_query['already_have'] = [] # list of tup(smiles, fraglib_dir)
frags_to_query['started_but_not_completed'] = []
frags_to_query['need'] = []

# Iterate over smiles
for smiles in list_smiles:
    # Convert SMILES to RDKit molecule and remove H's
    orig_mol = Chem.MolFromSmiles(smiles, sanitize=False)
    
    results = Frags.manually_remove_Hs(orig_mol, 'single') # rdkit's remove 
                                            # H method isn't good enough.
    
    mol_info, _ = results
    mol, _ = mol_info
    
    # Decompose the ligand into fragments and store the fragment SMILES. Use SMILES 
    # instead of SMARTS b/c only SMILES (from rdkit) differentiates aliphatic and 
    # aryl (C,c vs. [#6]). Fragment on bond radii `bond_radius` AND the postive 
    # integers less than `bond_radius`, because for example, drugs containing 
    # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O. 
    filtered_frags = Frags.get_fragments(2, mol, 4, 5)
    for sub_smiles, substruct_site_groups in filtered_frags.items():
        # Is this smiles represented in any fragment library?
        found_match = False
        for fraglib in fraglibs:
            for frag_smiles in os.listdir(fraglib):
                if utils.smiles_equiv(frag_smiles, sub_smiles, check_atom_order=False):
                    # Has the job finished yet? Check the log file in 
                    # fraglib/frag_smiles/logs/<frag_smiles>_log to find the line 
                    # "Job completed."
                    logfile = os.path.join(fraglib, frag_smiles, 'logs', 
                        f'{frag_smiles}_log')
                    if os.path.exists(logfile):
                        with open(logfile, 'r') as f:
                            loglines = f.readlines()
                        if any('Job completed.' in line for line in loglines):
                            frags_to_query['already_have'].append((sub_smiles, fraglib))
                        else:
                            frags_to_query['started_but_not_completed'].append(
                                (sub_smiles, fraglib))
                    else:
                        frags_to_query['started_but_not_completed'].append(
                            (sub_smiles, fraglib))
                    
                    found_match = True
        if not found_match:
            frags_to_query['need'].append(sub_smiles)
    
# pretty print the results
for key, val in frags_to_query.items():
    print(f'{key}:')
    for item in val:
        print(f'  {item}')
