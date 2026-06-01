import os
import re
import io
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import defaultdict
import prody as pr
import utils

def get_fragments(bond_radius, mol, min_frag_size=4, max_frag_size=5, quiet=True): 
    # Decompose the ligand into fragments and store the fragment SMILES. Use SMILES 
    # instead of SMARTS b/c only SMILES (from rdkit) differentiates aliphatic and 
    # aryl (C,c vs. [#6]). Fragment on bond radii `bond_radius` AND the postive 
    # integers less than `bond_radius`, because for example, drugs containing 
    # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O. 
    filtered_frags = {} # key=sub_smiles, value = list of Mol objs
    substructs = []
    for rad in list(range(1, bond_radius + 1)):
        substructs += fragment_on_bond_d(mol, rad)
    # Add substructs to `filtered_frags`.
    for orig_sub, orig_mol_inds in substructs: # contains H's that need to be scrubbed.
        # apply size threshold.
        frag_size = orig_sub.GetNumHeavyAtoms()
        if frag_size < min_frag_size or frag_size > max_frag_size:
            continue
        _results = manually_remove_Hs(orig_sub, return_single_mol_or_perms='perms')
        if _results is None:
            continue

        sub_perms_mols_inds, sub_smiles = _results
        for sub, perm_inds in sub_perms_mols_inds: 
            # skip if the elements are all carbons
            frag_elements = "".join(sorted([a.GetSymbol() for a in sub.GetAtoms() if 
                        a.GetSymbol() != 'H'])) # b/c rdkit will add implicit H's
            frag_num_carbons = len([i for i in frag_elements if i == 'C' or i == 'c'])
            if len(frag_elements) == frag_num_carbons:
                continue
            # add to dict
            substruct_data = (sub, perm_inds, orig_mol_inds)
            # get elements list based on perm_inds
            if sub_smiles not in filtered_frags:
                filtered_frags[sub_smiles] = [substruct_data]
            else:
                # check if this exact substructure (w/ same atom order) is already present.
                # not sure why, but there are duplicates within the same lig obj.
                already_present = False
                for existing_sub, existing_sub_inds, existing_orig_mol_inds in filtered_frags[sub_smiles]:
                    if (perm_inds, orig_mol_inds) == (existing_sub_inds, existing_orig_mol_inds):
                        already_present = True
                        break
                if not already_present:
                    filtered_frags[sub_smiles].append(substruct_data)

    # Group substructs by instances (sites) in the ligand. For example, if a substruct has
    # orig_mol_inds [14, 15, 16, 17] and another has [16, 17, 18, 19], they are the same 
    # site. The purpose of grouping is to prevent overcounting of matches when there are diff 
    # permutations of the same atoms or slight overlap of sites.
    grouped_frags = {}
    for sub_smiles, substruct_data in filtered_frags.items():
        groups = group_lig_sites_by_overlap(substruct_data)
        grouped_frags[sub_smiles] = groups
        if not quiet: 
            print(f"Fragment: {sub_smiles}, # perms: {len(substruct_data)}, "
                  f"# sites: {len(groups)}", flush=True)
    return grouped_frags

def group_lig_sites_by_overlap(data, key_index=2, threshold=0.5):
    """
    Input: list of tuples (substruct Mol obj, perm_inds, orig_mol_inds) that describe 
        instances of a frag in a lig. Determine whether each site has >1 instances, as 
        determined by sharing (overlapping) >=1/2 of the atoms b/n one frag instance and 
        another. For example, if one instance has orig_mol_inds [14, 15, 16, 17] and 
        another has [16, 17, 18, 19], they are the considered the same site on the lig.
        This avoids overcounting of matches when there are diff permutations of the CG.
    Output: list of groups.
    Method: Two frag instances are related if |A ∩ B| > threshold * max(|A|, |B|).
    """
    
    sets = [set(item[key_index]) for item in data]
    n = len(data)
    # Union–find
    parent = list(range(n))
    rank = [0]*n
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb: 
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1
    # Connect pairs that overlap 
    for i in range(n):
        for j in range(i+1, n):
            inter = len(sets[i] & sets[j])
            denom = max(len(sets[i]), len(sets[j])) # >= half
            if inter >= threshold * denom:
                union(i, j)
    # Collect groups
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(data[i])

    return list(groups.values())

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

def reorder_substruct_based_on_smiles(orig_substruct):
    # Convert the Mol obj (`sub`) to SMILES, and then use the SMILES to reorder 
    # the atoms in the Mol obj. This seems roundabout, but MolToSmiles() doesn't 
    # cleanly map atom orders, and we want this sane SMILES that correctly 
    # describes branching.

    # 0. Get initial SMILES
    smiles_w_H = Chem.MolToSmiles(orig_substruct, allHsExplicit=True, 
                              isomericSmiles=False) # shows [NH3+], [CH3], etc.
    # 1. Parse the SMILES into a Mol
    mol = Chem.MolFromSmiles(smiles_w_H)
    # 2. Annotate atoms with their SMILES position
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(i)
    # 3. Re-parse the SMILES so that GetAtoms() order = SMILES order
    mol2 = Chem.MolFromSmiles(
        Chem.MolToSmiles(mol, canonical=False, allHsExplicit=True, isomericSmiles=False))
    # 4. Build permutation: new index -> original index
    order = [a.GetAtomMapNum() for a in mol2.GetAtoms()]
    # 5. Renumber the molecule to that order
    reordered_sub = Chem.RenumberAtoms(mol, order)
    return smiles_w_H, reordered_sub

def manually_remove_Hs(orig_substruct, return_single_mol_or_perms):
    '''return_single_mol_or_perms must be 'single' (for processing a whole ligand) or 
    'perms' (when processing frags). '''
    # Remove hydrogens. RDKit docs say that Chem.RemoveHs() implicit and explicit are removed, 
    # but this isn't true for [nH], [OH], [Ho], etc. so need to manually remove H's. 
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
    # -- Parse the SMILES back into a new molecule
    mol_from_smarts = Chem.MolFromSmarts(smiles_no_Hs)
    # -- Map original atoms to new atom order using substructure matching
    try:
        matches_to_map_to_sub = orig_substruct.GetSubstructMatches(mol_from_smarts, 
            uniquify=False) # returns all permutations of order of atom inds
    except:
        return None
    # -- Reorder atoms in the original mol. matches_to_map_to_sub is a list of perms, so 
    #    we just need to use the first one
    cg_atom_perms = []
    for perm_inds in list(matches_to_map_to_sub):
        mol_copy = Chem.Mol(orig_substruct) # for immutable deep copy; don't use copy.deepcopy()
        try:
            renumbered_substruct = Chem.RenumberAtoms(mol_copy, perm_inds)
        except:
            continue
        cg_atom_perms.append((renumbered_substruct, perm_inds))
    if len(cg_atom_perms) == 0:
        return None
    if return_single_mol_or_perms == 'single':
        return cg_atom_perms[0], smiles_no_Hs # return the first permutation only
    elif return_single_mol_or_perms == 'perms':
        return cg_atom_perms, smiles_no_Hs # return all permutations. `cg_atom_perms` is a 
                                           # list of (substructure Mol objs, perm_inds)

def fragment_on_bond_d(mol, radius):
    # Code from https://iwatobipen.wordpress.com/2020/08/12/get-and-draw-molecular-fragment-with-user-defined-path-rdkit-memo/
    atoms = mol.GetAtoms()
    submols = []
    for atom in atoms:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom.GetIdx(), 
            enforceSize=False) # don't enforce size to also get frags that are < radius away
        amap = {}
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
        # Store the submol and its atom indices in the original mol to ensure that if there 
        # are >1 instances of a frag in a single ligand, they won't get skipped
        orig_inds = sorted(amap.keys())
        submols.append((submol, orig_inds))
    return submols

def is_organic(mol):
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if not any(elem in elements for elem in ['C', 'O', 'N']):
        return False
    else:
        return True

def get_frags_from_structure(struct_or_path, lig_smiles, quiet=False):
    """
    Fragment a ligand from either a prody obj or a path to a pdb/pdb.gz/cif file

    Returns:
        filtered_frags: dict[sub_smiles -> list of groups (sites)]
        pdb_mol_assigned_bonds_no_H: RDKit Mol (ligand, no H, single permutation)
    """
    
    if isinstance(struct_or_path, str):
        if struct_or_path.endswith('.pdb') or struct_or_path.endswith('.pdb.gz'):
            query_struct = pr.parsePDB(struct_or_path)
        elif struct_or_path.endswith('.cif'):
            query_struct = pr.parseCIF(struct_or_path)
    else:
        # assume prody obj
        query_struct = struct_or_path

    # Identify ligand and fragment it
    hetatms = query_struct.hetatm.select(
        'not (ion or resname SEP or resname TPO or resname MSE)')
    
    if hetatms is None:
        raise ValueError("get_frags_from_structure: could not find ligand HETATM atoms.")
    assert len(set(hetatms.getResnames())) == 1, "Expected exactly one ligand residue."

    # Convert from ProDy object to RDKit Mol via an in-memory PDB block
    buf = io.StringIO()
    pr.writePDBStream(buf, hetatms)
    pdb_block = buf.getvalue()
    pdb_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True)
    lig_template = Chem.MolFromSmiles(lig_smiles)

    # Assign bond orders (valence + aromaticity)
    try:
        pdb_mol_assigned_bonds = AllChem.AssignBondOrdersFromTemplate(lig_template, pdb_mol)
        pdb_mol_assigned_bonds_no_H_perm_inds, pdb_mol_smiles_no_H = manually_remove_Hs(
            pdb_mol_assigned_bonds, return_single_mol_or_perms='single')
        pdb_mol_assigned_bonds_no_H, pdb_mol_perm_inds = pdb_mol_assigned_bonds_no_H_perm_inds
        # Fragment the H-stripped ligand
        filtered_frags = get_fragments(2, pdb_mol_assigned_bonds_no_H, 4, 5, quiet=quiet)
        return filtered_frags, pdb_mol_assigned_bonds_no_H
    except ValueError as e:
        print(f"Error processing {struct_or_path}: {e}", flush=True)
        return {}, None

def check_vdg_job_status(sub_smiles, vdg_lib_dir):
    # Check if the vdg generation job finished without issues
    vdg_log_file = os.path.join(vdg_lib_dir, sub_smiles, f'{sub_smiles}_log')
    if not os.path.exists(vdg_log_file):
        return False
    with open(vdg_log_file, 'r') as f:
        log_contents = f.read()
    return 'Job completed.' in log_contents

def summarize_frags(frags_in_lib, frags_to_exclude, frags_to_include, logfile_fh):
    groups = {
        "Excluded": [],
        "Not in include list": [],
        "Not in vdg lib or incomplete": [],
        "In vdg lib and in include list": []}

    for sub_smiles in frags_in_lib:
        if sub_smiles in frags_to_exclude:
            groups["Excluded"].append(sub_smiles)
        elif (frags_to_include != 'all' and isinstance(frags_to_include, list) and 
              sub_smiles not in frags_to_include):
            groups["Not in include list"].append(sub_smiles)
        elif not frags_in_lib[sub_smiles]:
            groups["Not in vdg lib or incomplete"].append(sub_smiles)
        else:
            groups["In vdg lib and in include list"].append(sub_smiles)

    # Log summary
    if groups:
        print("\n--- Fragment vdG Library Summary ---", file=logfile_fh)
        for category, frags in groups.items():
            if frags:
                print(f"{category}:", file=logfile_fh)
                print("   ", ", ".join(frags), file=logfile_fh)
        print("\n", file=logfile_fh)
