import re
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
import prody as pr
import os
import utils

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

def _remove_bracket_hydrogens_in_smiles(smiles: str) -> str:
    """
    Strip explicit hydrogens from bracket atoms in a SMILES string (string-only).
    - [CH3] -> C
    - [NH2] -> N
    - [NH3+] -> [N+]
    - [H]   -> '' (remove the token)
    Brackets are dropped only if the remaining content is a plain element/atom symbol
    (letters only, e.g., C, N, n). If charge/chirality/isotope/class remains,
    brackets are kept.
    """
    # 1) Remove standalone [H]
    s = re.sub(r'\[H\]', '', smiles)

    # 2) Rewrite every bracket atom content
    def _rewrite(m):
        inner = m.group(1)

        # Remove explicit H or H<number> that is not part of an element symbol like Hg/He
        cleaned = re.sub(r'H(?![a-z])\d*', '', inner)

        # If anything special remains (charge '+/-', chirality '@', isotope/class digits,
        # colon, or any non-letter), we must keep the brackets (e.g., [N+], [C@H], [13C])
        if re.search(r'[\+\-\@\:\d]', cleaned) or re.search(r'[^A-Za-z]', cleaned):
            return f'[{cleaned}]'

        # Otherwise if it's only letters (e.g., C, N, n), drop brackets:
        # [CH3] -> C, [NH2] -> N, [nH] -> n
        return cleaned

    s = re.sub(r'\[([^\]]+)\]', _rewrite, s)
    return s


def manually_remove_Hs(orig_substruct):
    # Remove hydrogens. Docs say that Chem.RemoveHs() implicit and explicit are removed, 
    # but this isn't true for [nH], [OH], [Ho], etc. so need to manually remove H's. 
    # Convert back to SMILES, without showing explicit hydrogens

    # First, export SMILES *with explicit Hs shown* so that patterns like [NH3+] are present, 
    # then regex-strip bracket hydrogens while keeping charge and other annotations.
    smiles_w_H = Chem.MolToSmiles(orig_substruct, allHsExplicit=True, 
                              isomericSmiles=False) # shows [NH3+], [CH3], etc.

    # Remove standalone [H], remove H/Hn within brackets, keep charges (e.g., [NH3+] -> [N+])
    smiles_no_Hs = _remove_bracket_hydrogens_in_smiles(smiles_w_H)

    # Is there exactly one standalone character inside brackets? (i.e. the result of H 
    # stripping for [nH], [OH], etc.). If so, remove the brackets. (e.g. [n] -> n) 
    # NOTE: This may already be handled by the helper; this line is safe as a final pass.
    smiles_no_Hs = re.sub(r'\[([A-Za-z])\]', r'\1', smiles_no_Hs)

    # Complication: after converting the Mol obj to smiles, the Mol obj won't have the same 
    # atom order as the original Mol obj, which is important when extract CG coords. 
    # We need to rearrange the atom order of orig_substruct by getting the atom indices in 
    # the orig full molecule.
    # Parse the SMILES back into a new molecule: return the original substructure (Mol 
    # unchanged) and the H-stripped SMILES string
    return orig_substruct, smiles_no_Hs

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

def get_frags_from_pdbfile(pdbfile, lig_smiles):
    if pdbfile.endswith('.pdb') or pdbfile.endswith('.pdb.gz'):
        query_struct = pr.parsePDB(pdbfile)
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
        pdb_mol_assigned_bonds_no_H, pdb_mol_smiles_no_H = manually_remove_Hs(pdb_mol_assigned_bonds)
        filtered_frags = get_fragments(2, pdb_mol_assigned_bonds_no_H, 4, 5)
        return filtered_frags, pdb_mol_assigned_bonds_no_H
    except ValueError as e:
        print(f"Error processing {pdbfile}: {e}", flush=True)
        return [], None

def check_vdg_job_status(sub_smiles, vdg_lib_dir):
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
        "In vdg lib and in include list": []}

    for sub, sub_smiles in deduplicated_filtered_frags:
        if any(utils.smiles_equiv(f, sub_smiles) for f in frags_to_exclude):
            groups["Excluded"].append(sub_smiles)
        elif frags_to_include != 'all' and isinstance(frags_to_include, list) and sub_smiles not in frags_to_include:
            groups["Not in include list"].append(sub_smiles)
        elif sub_smiles not in frags_in_lib:
            groups["Not searched"].append(sub_smiles)
        elif not frags_in_lib[sub_smiles]:
            groups["Not in vdg lib or incomplete"].append(sub_smiles)
        else:
            groups["In vdg lib and in include list"].append(sub_smiles)

    print("\n=== Fragment summary ===", flush=True)
    for category, frags in groups.items():
        if frags:
            print(f"{category}:", flush=True)
            print("   ", ", ".join(frags), flush=True)
    print("\n========================", flush=True)

def check_in_exclude_list(sub_smiles, frags_to_exclude):
    for frag in frags_to_exclude:
        if utils.smiles_equiv(frag, sub_smiles):
            return True
    return False