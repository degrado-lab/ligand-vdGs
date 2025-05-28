from rdkit import Chem
import re

def get_fragments(bond_radius, mol, min_frag_size=4, max_frag_size=7):
    # Decompose the ligand into fragments and store the fragment SMILES. Use SMILES 
    # instead of SMARTS b/c only SMILES (from rdkit) differentiates aliphatic and 
    # aryl (C,c vs. [#6]). Fragment on bond radii `bond_radius` AND the postive 
    # integers less than `bond_radius`, because for example, drugs containing 
    # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O. 
    filtered_frags = []
    substructs = []
    for rad in list(range(1, bond_radius + 1)):
        substructs += fragment_on_bond_d(mol, rad)
    # Update frag_dict by adding the new fragments; skip if all-carbon
    for orig_sub in substructs: # contains H's that need to be scrubbed.
        # will have duplicates b/c we didn't check for duplicates yet.
        sub, sub_smiles = manually_remove_Hs(orig_sub) 
        # apply size threshold.
        frag_size = len([i for i in sub.GetAtoms() if i.GetSymbol() != 'H'])
        if frag_size < min_frag_size or frag_size > max_frag_size:
            continue
        # skip if the elements are all carbons
        frag_elements = "".join(sorted([a.GetSymbol() for a in sub.GetAtoms() if 
                    a.GetSymbol() != 'H'])) # b/c rdkit will add implicit H's
        frag_num_carbons = len([i for i in frag_elements if i == 'C' or i == 'c'])
        if len(frag_elements) == frag_num_carbons:
            continue
        filtered_frags.append((sub, sub_smiles))

    return filtered_frags

def manually_remove_Hs(orig_mol):
    # Remove hydrogens. Docs say that Chem.RemoveHs() implicit and explicit are removed, 
    # but this isn't true for [nH], [OH], [Ho], etc. so need to manually remove H's. 
    mol = Chem.RemoveHs(orig_mol, sanitize=False)
    # Convert back to SMILES, without showing explicit hydrogens
    smiles = Chem.MolToSmiles(mol, allHsExplicit=False, isomericSmiles=False) # still has H's
    # Remove standalone [H] completely
    smiles_no_Hs = re.sub(r'\[H\]', '', smiles)
    # Remove all 'H' except when part of '[Hg]'
    smiles_no_Hs = re.sub(r'H(?!g\])', '', smiles_no_Hs)
    # Is there exactly one standalone character inside brackets? (i.e. the result of H 
    # stripping for [nH], [OH], etc.). If so, remove the brackets. (e.g. [n] -> n) 
    smiles_no_Hs = re.sub(r'\[([^\[\]])\]', r'\1', smiles_no_Hs)

    return mol, smiles_no_Hs

def fragment_on_bond_d(mol, radius):
    # Code from https://iwatobipen.wordpress.com/2020/08/12/get-and-draw-molecular-fragment-with-user-defined-path-rdkit-memo/
    atoms = mol.GetAtoms()
    submols = []
    for atom in atoms:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom.GetIdx(), 
            enforceSize=False) # don't enforce size to also get frags that are < radius away
        amap = {}
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
        submols.append(submol)
    return submols

def is_organic(mol):
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if not any(elem in elements for elem in ['C', 'O', 'N']):
        False
    else:
        return True
