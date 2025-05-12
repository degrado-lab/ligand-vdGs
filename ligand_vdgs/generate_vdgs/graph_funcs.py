import os
import sys
import time
start_time = time.time()
import torch
from torch_geometric.data import Data
from rdkit import Chem
#from rdkit.Chem import AllChem
import numpy as np
from scipy.spatial import cKDTree

#import matplotlib.pyplot as plt
import prody as pr
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import graph_utils as utils

'''
If graph-building is too slow, considering turning off hydrogen counting.
'''

def main():
    vdgs_dir = '/wynton/group/degradolab/skt/docking/databases/vdg_lib/'
    cg = 'carbamate'
    # -----------------------------------------------------------------------
    nr_vdgs_path = os.path.join(vdgs_dir, cg, 'nr_vdgs')
    for size_subset in os.listdir(nr_vdgs_path):
        size_subset_path = os.path.join(nr_vdgs_path, size_subset)
        for subset_AAs in os.listdir(size_subset_path):
            # TODO: ADD DESCRIPTION
            process_subset(size_subset_path, subset_AAs, vdgs_dir, cg) 

def parse_protein_atoms(pdbpath, ligand_resname='LIG'):
    try:
        par = pr.parsePDB(pdbpath).protein
    except:
        print('COULD NOT OPEN', pdbpath)
        return
    protein_atoms = []
    residue_info = []
    for res in set(par.getResindices()):
        res_obj = par.select(f'resindex {res} not element H D')
        b = res_obj[0]
        seg, chain, resnum, prot_resname = b.getSegindex(), b.getChid(), b.getResnum(), b.getResname()
        for atom_obj in res_obj: 
            xyz, charge = atom_obj.getCoords(), atom_obj.getCharge() # formal charge, not partial
            z = Chem.GetPeriodicTable().GetAtomicNumber(atom_obj.getElement())
            atom_name = atom_obj.getName()
            num_Hs = utils.prody_num_hyds(atom_obj, par) 
            arom = utils.prody_is_arom(atom_name, atom_obj.getResname())
            hyb = 0 # hybridization state not in PDB
            prot_flag = 1  # if ligand, prot_flag = 0 
            features = torch.tensor([z, charge, num_Hs, hyb, arom, prot_flag], dtype=torch.float32)  
            # Add to lists
            protein_atoms.append((features, xyz))
            residue_info.append(f"{prot_resname}{resnum} {atom_name}")
    return protein_atoms, residue_info

def build_graph(protein_atoms_data, residue_info, bindingsite_data, cg, cg_rdmol):
    pos_list = [] # coords
    node_feats = []
    node_residues = []
    node_ligand_ids = []

    # Add ligand atoms
    lig_key = [i for i in bindingsite_data.keys() if i[4] == 'LIG']
    assert len(lig_key) == 1
    lig_key = lig_key[0]
    lig_dict = bindingsite_data[lig_key]
    for ligatom_name, ligatom_data in lig_dict.items():
        coords = ligatom_data['coords']
        pos_list.append(coords)
        features = torch.tensor([
            ligatom_data['atomic_num'],
            ligatom_data['partial_charge'],
            ligatom_data['num_Hs'],
            ligatom_data['hybridization'],
            ligatom_data['arom'],
            0], # prot_flag = 0 for ligand
            dtype=torch.float32) 
        node_feats.append(features) 
        node_residues.append(ligatom_name)
        node_ligand_ids.append(cg) # use cg name instead of lig name

    # Add protein atoms
    for i, (feat, xyz) in enumerate(protein_atoms_data):
        pos_list.append(xyz)
        node_feats.append(feat)
        node_residues.append(residue_info[i]) # resname
        node_ligand_ids.append("PROT")

    pos = torch.tensor(np.array(pos_list), dtype=torch.float32)
    x = torch.stack(node_feats) # node feature matrix

    # Add edges
    edge_index = []
    edge_attr = []

    # Covalent bonds within ligand
    for bond in cg_rdmol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        feat = bond_features(bond)
        edge_index.append([i, j])
        edge_index.append([j, i])
        edge_attr.append(feat)
        edge_attr.append(feat)

    # Add protein atom bonds
    for i, atom1 in enumerate(protein_atoms_data):
        for j, atom2 in enumerate(protein_atoms_data):
            if i < j:
                distance = np.linalg.norm(atom1[1] - atom2[1])  
                if distance < 1.8:  # Then covalent bond 
                    bond = torch.tensor([1, 1], dtype=torch.float32)  # Example: using single bond
                    edge_index.append([i + len(lig_dict), j + len(lig_dict)])  # Adjust indices
                    edge_index.append([j + len(lig_dict), i + len(lig_dict)])
                    edge_attr.append(bond)
                    edge_attr.append(bond)
    
    
    edge_index = torch.tensor(edge_index).t().contiguous()
    edge_attr = torch.stack(edge_attr)

    # Encode residue info and ligand ID (as metadata)
    data = Data(x=x, pos=pos, edge_index=edge_index, edge_attr=edge_attr)
    data.residue_info = node_residues
    data.ligand_ids = node_ligand_ids
    return data

def bond_features(bond):
    # map bond type to numerical value
    bond_type = bond.GetBondType()
    bond_map = {
        Chem.rdchem.BondType.SINGLE: 1,
        Chem.rdchem.BondType.DOUBLE: 2,
        Chem.rdchem.BondType.TRIPLE: 3,
        Chem.rdchem.BondType.AROMATIC: 4}
    
    bond_val = bond_map.get(bond_type, 0)
    return torch.tensor([bond_val, 1], dtype=torch.float32)  # is_covalent = 1

def get_features_from_bindingsite(bindingsite_path, _pdbname):
    features_dict = {} 
    '''
    features_dict = {residue: residue subdict} 
            residue_subdict = {atom name: atom_subdict}
                    atom_subdict = {feat name: feat value}
    '''
    
    # SHOULD BE ITERATING OVER LIG RESIDUE AND PROT RESIDUES
    # key = atom name, val = residue subdict, 
    par = pr.parsePDB(bindingsite_path)
    '''
    Select lig residue
    '''
    lig_scr = _pdbname.split('_')[1:4]
    seg, chain, resnum = lig_scr
    if seg: 
        sele = f'segname {seg} and chain {chain} and resnum {resnum}' 
    else:
        sele = f'chain {chain} and resnum {resnum}'

    lig_pr_obj = par.select(sele)
    frag_atoms = lig_pr_obj.select('occupancy > 2.9') # occ=3.0 is CG. look at CG only, 
                                                   # not the entire lig.
    frag_names = frag_atoms.getNames()
    assert len(frag_names) == len(set(frag_names)) 
    res_descript = (bindingsite_path, seg, chain, resnum, 'LIG') # !!!!!!! LIG!!!
    assert res_descript not in features_dict.keys()
    features_dict[res_descript] = {} 
    '''
    TODO: can do protein residues at some point later 
    '''
    res_dict, cg_rdmol = utils.get_feats_for_res(lig_pr_obj, bindingsite_path, sele, 
                                                 frag_names, lig_scr)
    if type(res_dict) == tuple or res_dict is None: # failed
        return res_dict, None

    features_dict[res_descript] = res_dict
    return features_dict, cg_rdmol


def process_subset(size_subset_path, subset_AAs, vdgs_dir, cg):
    # TODO: ADD DESCRIPTION
    failed = []
    passed = []
    for num_pdb, vdg_pdbfile in enumerate(os.listdir(os.path.join(size_subset_path, 
                                                                  subset_AAs))):
        # Get PDB of just the CG + vdMs 
        vdg_pdbpath = os.path.join(size_subset_path, subset_AAs, vdg_pdbfile)
        # Get PDB of CG + all binding site residues
        b_site_path, _pdbname = utils.get_binding_site_path(vdg_pdbpath, 
                                                            vdgs_dir, cg)
        bindingsite_results = get_features_from_bindingsite(b_site_path, _pdbname)
        if bindingsite_results is None:
            continue
        if len(bindingsite_results) == 2:
            bindingsite, cg_rdmol = bindingsite_results
        else:
            print('WTF')
            print(type(bindingsite_results))
            print(bindingsite_results)
        if type(bindingsite) == tuple:
            failed.append(bindingsite)
            continue # skip over
        passed.append(bindingsite)
        protein_atoms_data, residue_info = parse_protein_atoms(vdg_pdbpath)
        graph_data = build_graph(protein_atoms_data, residue_info, bindingsite, cg, cg_rdmol)
    #print('='* 60)
    #print('SUBSET', subset_AAs)
    #print('PASSED:', len(passed))
    #print('FAILED:', len(failed))
    #assert len(passed) + len(failed) == len(os.listdir(os.path.join(size_subset_path, subset_AAs)))
    #print('='* 60)

    # Visualize graph
        utils.visualize_graph(graph_data, node_label_attr='residue_info', 
                              title_label=str(subset_AAs)+'_'+_pdbname)  

main()
end_time = time.time()

# Calculate the elapsed time
elapsed_time = end_time - start_time

minutes = elapsed_time // 60  # Whole minutes
seconds = elapsed_time % 60   # Remaining seconds

print(f"Time taken: {int(minutes)} minutes and {seconds:.2f} seconds")
