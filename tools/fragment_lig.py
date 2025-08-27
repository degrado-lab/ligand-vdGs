import os
import sys
from rdkit import Chem            
import pickle as pkl
sys.path.append("/wynton/home/degradolab/skt/docking/ligand-vdGs")
from ligand_vdgs.functions import Frags

#smiles = 'C[C@]1([C@H]2C[C@H]2OC(=N1)N)c3cc(ccc3F)NC(=O)c4cnc(cn4)OC' # EJ7 (6c2i)
#smiles = 'CCCCn1c(cn2c1nc3c2C(=O)NC(=O)N3C)c4ccccc4C' # L87 (4gk3)
#smiles = 'CC(C)c1ccc2c(c1)C(=O)c3cc(c(nc3O2)N)c4[nH]nnn4' # E0M (6bny)
#smiles = 'c1ccc2c(c1)C(=O)c3ccc(cc3S2(=O)=O)c4[nH]nnn4' # 91X (5yva)
#smiles = 'N=C(N)N[C@H]1C[C@H](O[C@H]([C@@H]1NC(=O)C)[C@@H]([C@@H](CO)O)O)S(=O)(=O)O' # GYG (6br5) 
#smiles = 'CC(=O)N1C[C@@H](C[C@H]1C(=O)N[C@@H](CCCNC(=N)N)C(=O)c2nc3ccccc3s2)O' # AAB (1NC6) 
###### Metalloproteins w/ hydroxamic acids: #####
#smiles = 'CC(C)C[C@H]([C@H](CSc1cccs1)C(=O)NO)C(=O)N[C@@H](Cc2ccccc2)C(=O)NC' # BAT (1rm8)
#smiles = 'CC(C)C[C@H]([C@@H](C(=O)NO)O)C(=O)N[C@H](C(=O)NC)C(C)(C)C ' # 097 (3hy7)
#smiles = 'CC(C)C[C@H]([C@H](CNC(=O)c1nccs1)C(=O)NO)C(=O)N[C@H](C(=O)NC)C(C)(C)C' # WR2 (2w12, 2w15)
###### Kinases w/ ureas: #####
# * smiles = 'COc1cc(c(cc1NC(=O)Nc2cnc(cn2)C#N)Cl)OC' # A42 (2ywp)
#smiles = 'CC(C)(C)c1cc(no1)NC(=O)Nc2ccc(cc2)Oc3ccncc3' # SR8 (2oh4)
# * smiles = 'COC(=O)Nc1[nH]c2ccc(cc2n1)Oc3ccc(cc3)NC(=O)Nc4cc(ccc4F)C(F)(F)F' # GIG (2oh4)
#smiles = 'Cc1cc(ccn1)c2c3cnc(cc3[nH]n2)NC(=O)NCc4ccccc4' # S69 (5ke0)
# * smiles = 'CC1=CN(C(=O)NC1=O)[C@@H]2C[C@@H]([C@H](O2)CNC(=O)Nc3ccc(cc3)[N+](=O)[O-])O' # WMJ (2yoh)
#smiles = 'C[C@H](c1c(cncc1Cl)Cl)Oc2cc(ncc2C#N)NC(=O)Nc3c(cc(c(n3)C=O)CN4CCN(CC4=O)C)OCCN5CCOCC5' # VVW (8kh7)
# * smiles = 'CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CNC(=S)Nc3ccc(c(c3)C(F)(F)F)Cl)O' # 74W (2yof)
# * smiles = 'CCN(CC)S(=O)(=O)c1cc(ccc1Cl)Nc2nccc(n2)c3ccnc(c3)c4ccc(cc4)NC(=O)NC' # 88C (7pw7)
# * smiles = 'c1cc2c(cc1CN3CCOCC3)nc([nH]2)c4c(c[nH]n4)NC(=O)NC5CC5' # 35R (5n23)
###### Factor Xa or thrombin inhibitors: #####
# * smiles = r'[H]/N=C(/c1ccc(cc1)NCc2nc3cc(ccc3n2C)Cn4c(nc5c4cccc5)C)\N' # R11 (1g2m)
# * smiles = 'C[C@@H](C(=O)N1CCOCC1)N2CC[C@@H](C2=O)NS(=O)(=O)\C=C(/C)\c3ccc(s3)Cl' # 701 (2uwo)
# smiles = 'CCc1c(cccc1S(=O)(=O)N[C@@H](Cc2cc(on2)c3ccc(s3)Cl)C(=O)N4CCC(CC4)OC)N5CCOCC5=O' # VYR (4btt) # there's an insertion code res 61 that prody processes incorrectly
#smiles = '[H]/N=C(\C)/N1CC[C@@H](C1)Oc2ccc(cc2)[C@H](Cc3ccc4ccc(cc4c3)/C(=N/[H])/N)C(=O)O' # DX9 (1fax)
#smiles = 'CC(=N)N1CCC(CC1)Oc2ccc(cc2)[C@H](Cc3ccc4ccc(cc4c3)C(=N)N)C(=O)O' # BX3 (1mtu)
#smiles = 'c1ccc(c(c1)Cc2nc3c(o2)ccnc3NCC([C@H]4CCCCN4)(F)F)n5cncn5' # 382 (1zgi) # insertion code that prody does wrong
#smiles = 'CN(C)[C@H]1CCN(C1)C(=O)[C@H](CNC(=O)c2ccc(s2)Cl)NS(=O)(=O)c3cccc(c3OC(F)F)N4CCCCC4=O' # 7R9 (4lxb) # insertion code
#smiles = 'CC(C)N1CCC(CC1)NC(=O)c2cc3ccccc3n2Cc4cc(on4)c5ccc(s5)Cl' # IIA (2boh) # insertion code
#smiles = r'[H]/N=C(\N)/N1CC=C(C1)CCNC(=O)[C@@H]2C[C@@H]3C[C@@H]([C@H](C[C@@H]3N2C(=O)[C@@H]([C@@H](C(C)C)Cl)NC(=O)[C@@H](COS(=O)(=O)O)OC)O)O' # SN3 (2gde)
###### Muscarinic receptors: #####
# * smiles = r'C[N+]1([C@@H]2CC(C[C@H]1[C@H]3[C@@H]2O3)OC(=O)Nc4ccc(cc4c5cccs5)F)C' # 9EC (5zhp)
#smiles = r'C[N+](C)(C)CC#CCOC1=NOCC1' # IXO (7trk, 6oij, 6oik)
#smiles = r'c1ccc(cc1)C(c2ccccc2)(C(=O)O[C@H]3CN4CCC3CC4)O' # QNB (3uon)
#smiles = r'CN1[C@H]2CC[C@@H]1CC(C2)OC(=O)[C@H](CO)c3ccccc3' # OIN (6wjc)
#smiles = r'C[N+]1([C@@H]2CC(C[C@H]1[C@H]3[C@@H]2O3)OC(=O)[C@H](CO)c4ccccc4)C' # 3C0 (5zkc)
#smiles = r'C[N+]1([C@@H]2CC(C[C@H]1[C@H]3[C@@H]2O3)OC(=O)C(c4cccs4)(c5cccs5)O)C' # 0HK (4daj)
###### B-adrenergic GPCRs w/ ??anilines??: #####
#smiles = r'CN[C@@H]1CCc2c(ccc(c2O)O)[C@H]1O' # G1I (8uo4, 8ge2
#smiles = r'Cc1ccccc1CC(C)(C)NC[C@@H](c2ccc(c3c2OCC(=O)N3)O)O' # P0G (3sn6)
#smiles = r'CC(C)NC[C@@H](COc1ccccc1CC=C)O' # JTZ (8jjo)
#smiles = r'CCCc1ccccc1OC[C@H](CNCCCc2cn(nn2)CCCCc3ccc(cc3)OC[C@H](CNC(C)C)O)O' # A1AE2 (8w1v)
#smiles = r'CCOC(=O)c1c(c2c(o1)cccc2OC[C@H](CNC(C)C)O)C' # JSZ (3ny9)
#smiles = r'CC(C)NC[C@@H](COc1cccc2c1c3ccccc3[nH]2)O' # CAU (2rh1)
###### aIIbB3 carboxylates: #####
#smiles = r'c1cc(cnc1)S(=O)(=O)N[C@H](CNC(=O)c2cc3cc(sc3s2)CCC4CCNCC4)C(=O)O' # 180 (2vc2)
#smiles = r'CCCCS(=O)(=O)N[C@@H](Cc1ccc(cc1)OCCCCC2CCNCC2)C(=O)O' # AGG (2vdm)
####### avB3 guanidines: #####
#smiles = r'c1ccc2c(c1)nc(s2)C(=O)N[C@@H](CCNC(=O)CCCCc3ccc4cccnc4n3)C(=O)O' # JUY (6mk0)
####### Factor Xa proline-like/pyrrolidine: #####
smiles = r'COc1ccc(cc1)n2c3c(c(n2)C(=O)N)CCN(C3=O)c4ccc(cc4)N5CCCCC5=O' # GG2 (6w70_chainA, 6w70_chainB, 2p16)
#smiles = r'c1cc(ccc1N2CCOCC2=O)N3C[C@@H](OC3=O)CNC(=O)c4ccc(s4)Cl' # RIV (2w26)

# idk what mols these are: 
#smiles = r'CN(C1=CC=C(C(O2)=C1)C(C3=CC=CC=C3C([N-]S(=O)(N(C)C)=O)=O)=C(C=C/4)C2=CC4=[N+](C)/C)C'
#smiles = r'CN(C1=CC=C(C(O2)=C1)C(C3=CC=CC=C3C([N]S(=O)(N(C)C)=O)=O)=C(C=C/4)C2=CC4=[N](C)/C)C'

database_frags_path = 'resources/database_frags_dict.pkl'
bond_radius = 2
min_frag_size = 4
max_frag_size = 5

def main():

    frags = []

    # Convert SMILES to RDKit molecule and remove H's
    orig_mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol, smiles_no_Hs = Frags.manually_remove_Hs(orig_mol) # rdkit's remove H methods 
                                                     # don't work 

    # Make sure there's at least a C, O, or N (rule out Fe-S clusters, ions, etc.)
    if not Frags.is_organic(mol):
        print('compound is not organic:', smiles) 
        sys.exit(1)
    
    # Decompose the ligand into fragments and store the fragment SMILES. Use SMILES 
    # instead of SMARTS b/c only SMILES (from rdkit) differentiates aliphatic and 
    # aryl (C,c vs. [#6]). Fragment on bond radii `bond_radius` AND the postive 
    # integers less than `bond_radius`, because for example, drugs containing 
    # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O.
    filtered_frags = Frags.get_fragments(bond_radius, mol, min_frag_size, max_frag_size)
    for (sub, sub_smiles) in filtered_frags:
        if sub_smiles not in frags:
            frags.append(sub_smiles)

    is_in_database_frags(frags)

def is_in_database_frags(lig_frags):
    # Determine whether these smiles are in database_frags_dict.pkl
    db_frags = pkl.load(open('resources/database_frags_dict.pkl', 'rb'))
    dict_smiles = []
    for elements, subdict in db_frags.items():
        for smiles, lignames in subdict.items():
            dict_smiles.append(smiles)

    was_not_run = []
    ran_not_finished = []
    ran_and_finished = []
    not_database =[]

    for lig_smiles in lig_frags: 
        lig_smiles_mol = Chem.MolFromSmarts(lig_smiles)
        # check against all database frags
        found_dict_match = False
        for _dict_smiles in dict_smiles:
            if 'K' in _dict_smiles:
                continue
            dict_smiles_mol = Chem.MolFromSmarts(_dict_smiles)

            if _dict_smiles == lig_smiles:
                found_dict_match = _dict_smiles
                break
            elif lig_smiles_mol.HasSubstructMatch(dict_smiles_mol
                ) and dict_smiles_mol.HasSubstructMatch(lig_smiles_mol): 
                found_dict_match = _dict_smiles 
                break
        if found_dict_match:
            # Was this fragment run/completed?
            vdglib_dir = '/wynton/group/degradolab/skt/docking/databases/frag_vdg_lib/'
            vdg_dir = os.path.join(vdglib_dir, found_dict_match)
            if not os.path.exists(vdg_dir):
                was_not_run.append(found_dict_match)
                continue
            else:
                logfile = f'/wynton/group/degradolab/skt/docking/databases/frag_vdg_lib/{found_dict_match}/logs/{found_dict_match}_log'
                with open(logfile, 'r') as inF:
                    lines = inF.readlines()
                if 'Job completed.' in lines[-1]:
                    ran_and_finished.append(found_dict_match)
                else: 
                    ran_not_finished.append(found_dict_match)

        else:
            not_database.append(lig_smiles)
    if not_database:
        print('SMILES NOT IN DATABASE:', not_database)
    
    print('SMILES IN DATABASE:')
    print_counts(was_not_run, db_frags, '> Smiles that I did not run yet:')
    print_counts(ran_not_finished, db_frags, '> Smiles that I ran but did not finish yet:')
    print_counts(ran_and_finished, db_frags, '> Smiles that I ran and finished:')

def print_counts(_list, db_frags, str_to_print):
    print()
    print(str_to_print)
    smiles_counts = {}
    for db_smiles in _list:
        for elements, subdict in db_frags.items():
            if db_smiles in subdict.keys():
                smiles_counts[db_smiles] = len(subdict[db_smiles])
    # Sort by number of ligands
    for k, v in sorted(smiles_counts.items(), key=lambda item: item[1], reverse=True):
        print(f'{k}: {v} ligands')
    



main()
