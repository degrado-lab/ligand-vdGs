import os
from io import StringIO
from rdkit import Chem
import prody as pr
import matplotlib.pyplot as plt
from openff.toolkit import topology
import espaloma as esp
import networkx as nx
from torch_geometric.utils import to_networkx
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
GLOBAL_TOOLKIT_REGISTRY._toolkits.pop(0)

def get_binding_site_path(vdg_pdbpath, vdgs_dir, cg):
    # Using pdb name from nr_vdgs (which only contains CG + vdMs), get the path for the 
    # pdb containing the binding site (CG + all binding site residues) for calculating 
    # partial charges, SASA (if applicable), etc.
    _pdbname = vdg_pdbpath.split('/')[-1]
    biol_assem_cg = '_'.join(_pdbname.split('_')[:5]) + '.pdb.gz'
    bindingsite_path = os.path.join(vdgs_dir, cg, 'vdg_pdbs', biol_assem_cg)
    return bindingsite_path, _pdbname 

def calc_partial_charges(molecule, frag_names):
    # Use espaloma to calculate partial charges, given an openff molecule obj
    molecule_graph = esp.Graph(molecule)
    espaloma_model = esp.get_model("latest") # load pre-trained model
    espaloma_model(molecule_graph.heterograph) # apply the model
    # create an OpenMM System and calculate the partial charges
    try:
        openmm_system = esp.graphs.deploy.openmm_system_from_graph(molecule_graph)
        # retrieve partial charges from molecule obj
        dict_from_esp = {} # key = atom name, val = partial_charge 
        for a in molecule.atoms:
            name, partial_charge = (a.name, a.partial_charge.magnitude)
            #if name not in frag_names:
            #    continue
            assert name not in dict_from_esp.keys()
            dict_from_esp[name] = partial_charge
            import numpy as np
            partial_charge = np.round(partial_charge, 2)
        return dict_from_esp
    except Exception as e:
        print(f"Error: {e}")
        return None
    
def prody_num_hyds(atom_obj, par):
    neighbs = pr.findNeighbors(atom_obj, 1.1, par.select('element H'))
    return len(neighbs)

def prody_is_arom(atom_name, resname):
    if resname not in ['PHE', 'TYR', 'TRP']:
        return False
    if atom_name in ['CG', 'CD2', 'CD1', 'CD2', 'CE1', 'CE2', 'CH2', 'NE1', 'CE3', 'CZ3', 'CZ2', 'CZ']:
        return True
    else:
        return False

def get_rdkit_feats(atom):
    # Get atomic number, charge, aromaticity, hybridization state, and number of Hs
    atomic_num = atom.GetAtomicNum()
    hybridization = atom.GetHybridization()
    is_aromatic = atom.GetIsAromatic()
    is_in_ring = atom.IsInRing()
    total_degree = atom.GetTotalDegree() # incl. H's
    degree_noH = atom.GetDegree() # no H's
    num_Hs = atom.GetTotalNumHs()
    return (atomic_num, hybridization, is_aromatic, is_in_ring, total_degree, 
            degree_noH, num_Hs)

def add_rdkit_feats_to_dict(rdmol_atom, features_dict, CG_atom_name):
    (atomic_num, hybridization, is_aromatic, is_in_ring, total_degree, 
        degree_noH, num_Hs) = get_rdkit_feats(rdmol_atom)
    features_dict[CG_atom_name]['atomic_num'] = atomic_num
    features_dict[CG_atom_name]['hybridization'] = hybridization
    features_dict[CG_atom_name]['arom'] = is_aromatic
    features_dict[CG_atom_name]['in_ring'] = is_in_ring
    features_dict[CG_atom_name]['total_degree'] = total_degree
    features_dict[CG_atom_name]['degree_noH'] = degree_noH
    features_dict[CG_atom_name]['num_Hs'] = num_Hs
    return features_dict

def prody_to_rdmol_and_openff_mol(prody_obj, bindingsitepdb, prody_sel):
    # Convert prody obj to an rdkit mol and openff Topology obj to run espaloma. 
    # To avoid writing and reading from disk, write to pdb stream.

    with StringIO() as pdb_stream: # write to memory buffer
        pr.writePDBStream(pdb_stream, prody_obj)
        pdb_data = pdb_stream.getvalue()

    # Convert to RDKit Mol
    rdmol = Chem.MolFromPDBBlock(pdb_data, sanitize=False, removeHs=False)
    num_atoms_in_rdmol = rdmol.GetNumAtoms()
    
    # Convert to openFF Molecule
    try: 
        ff_mol = topology.Molecule.from_rdkit(rdmol, hydrogens_are_explicit=True,
    
                                              allow_undefined_stereo=True)
    
        num_atoms_in_ff_mol = len(ff_mol.atoms) 
    
        if num_atoms_in_rdmol != num_atoms_in_ff_mol:
    
            print('WARNING: this lig was processed incorrectly. The # of atoms in the rdkit '
    
                  'mol must be the same as the # of atoms in the openff mol.')
    
            print(bindingsitepdb) 
    
            print('Ligand selection:', prody_sel)
    
        return rdmol, ff_mol
    except:
        return None, None

def get_feats_for_res(pr_obj, bindingsite_path, sele, list_atom_names, scr):
    res_dict = {}
    # If it's a lig residue, only get features for the CG atoms, not the entire lig
    assert len(set(pr_obj.getResindices())) == 1 # select only 1 res
    seg, chain, resnum = scr
    # Convert prody  obj to rdkit mol and openff mol objs
    rdmol, openff_mol = prody_to_rdmol_and_openff_mol(pr_obj, bindingsite_path, sele)
    if rdmol is None:
        return None, None
    Chem.SanitizeMol(rdmol) # sanitize the mol to get hybridization
    '''
    err = Chem.SanitizeMol(rdmol, catchErrors=True) # sanitize to get hybridization 
    if err != 0: # if it returned 0, then everything was successful
        print(err)
    '''
    dict_from_esp = calc_partial_charges(openff_mol, list_atom_names)
    if dict_from_esp is None:
        failed = (bindingsite_path, seg, chain, resnum)
        return failed, None
    # Gather features from mol objects
    cg_atom_idxs = []
    for CG_atom_name in list_atom_names:
        assert CG_atom_name not in res_dict.keys()
        res_dict[CG_atom_name] = {}
        rdmol_atom = [i for i in rdmol.GetAtoms() if 
                      i.GetPDBResidueInfo().GetName().strip() == CG_atom_name] 
        assert len(rdmol_atom) == 1
        rdmol_atom = rdmol_atom[0]
        # Get charge from openff / espaloma
        res_dict[CG_atom_name]['partial_charge'] = dict_from_esp[CG_atom_name] 
        # Get coords from prody
        coords = pr_obj.select(f'name {CG_atom_name}')[0].getCoords()
        res_dict[CG_atom_name]['coords'] = coords
        # Get features from RDKit mol
        res_dict = add_rdkit_feats_to_dict(rdmol_atom, res_dict, CG_atom_name)
        # Add atom idx to the list of CG atom idxs
        cg_atom_idxs.append(rdmol_atom.GetIdx())

    # Make a new rdmol obj that's just the CG. DO NOT use PathToSubmol because despite 
    # specifying the exact atom indices to put in the substruct, it adds other atoms
    # automatically, and the behavior can't be turned off. Instead, create a new obj.

    cg_rdmol = Chem.RWMol() # will be referred to as submol
    submol_map_dict = {}
    # Add atoms
    for _submol_idx, _orig_idx in enumerate(cg_atom_idxs):
        atom = rdmol.GetAtomWithIdx(_orig_idx)
        cg_rdmol.AddAtom(atom)
        submol_map_dict[_orig_idx] = _submol_idx

    # Add bonds
    for original_idx1 in cg_atom_idxs: # bond start
        for original_idx2 in cg_atom_idxs: # bond end
            if not original_idx1 < original_idx2: # don't add reciprocal bonds
                continue
            bond = rdmol.GetBondBetweenAtoms(original_idx1, original_idx2)
            if not bond:
                continue
            cg_rdmol.AddBond(submol_map_dict[original_idx1], 
                                submol_map_dict[original_idx2], bond.GetBondType())


    return res_dict, cg_rdmol

#def visualize_graph(data, title_label, node_label_attr='residue_info', color_by='ligand_ids', figsize=(10, 8)):
#    """
#    Visualize graph using networkx. Color by ligand/protein and show residue or atom labels.
#    """
#    G = to_networkx(data, to_undirected=True)
#    pos_3d = data.pos.cpu().numpy()
#    pos_2d = {i: (p[0], p[1]) for i, p in enumerate(pos_3d)}
#
#    if color_by == 'ligand_ids':
#        colors = ['skyblue' if lid == 'PROT' else 'pink' for lid in data.ligand_ids]
#    else:
#        colors = ['gray'] * len(data.x)
#
#    labels = None
#    if node_label_attr == 'residue_info':
#        labels = {i: data.residue_info[i] for i in range(len(data.residue_info))}
#    #elif node_label_attr == 'atom_type':
#    #    labels = {i: int(data.x[i][0].item()) for i in range(data.x.size(0))}
#
#    plt.figure(figsize=figsize)
#    nx.draw(G, pos=pos_2d, with_labels=True if labels else False,
#            node_color=colors, labels=labels,
#            node_size=200, font_size=12, edge_color='gray', alpha=0.8)
#    #plt.title('Protein-Ligand Binding Site Graph')
#    plt.axis('off')
#    plt.savefig(f'graph_visualization_{title_label}.png', dpi=300, bbox_inches='tight')
