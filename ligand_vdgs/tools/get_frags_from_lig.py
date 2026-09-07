'''
For a list of ligands, determine if a vdg library for each of their fragments exist, 
and identify fragments that don't currently have a vdg library.
'''

import os
from rdkit import Chem
from ligand_vdgs.functions import Frags, utils

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
    if orig_mol is None:
        print(f'[WARNING] could not parse SMILES: {smiles}')
        continue
    # Perceive aromaticity before fragmenting: unsanitized, a Kekule-written ring
    # keys as [C;r6]=[C;r6](...) instead of cc(...), which silently matches nothing
    # in the library. fragment_database_ligs sanitizes at the same point.
    try:
        Chem.SanitizeMol(orig_mol)
    except Exception as e:
        print(f'[WARNING] sanitization failed ({e}): {smiles}')
        continue

    results = Frags.manually_remove_Hs(orig_mol, 'single') # rdkit's remove
                                            # H method isn't good enough.
    if results is None:
        # H-removal/canonicalization failed (see the warning it printed); skip
        # this ligand rather than aborting the whole scan.
        print(f'[WARNING] could not process SMILES: {smiles}')
        continue
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
                if utils.smiles_equiv(frag_smiles, sub_smiles):
                    if Frags.check_vdg_job_status(frag_smiles, fraglib):
                        frags_to_query['already_have'].append((sub_smiles, fraglib))
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
