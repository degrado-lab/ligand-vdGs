import os
import shutil
from collections import defaultdict
from rdkit import Chem

def handle_existing_files(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    if len(os.listdir(out_dir)) > 0:
        raise ValueError(f'The output dir {out_dir} is not empty. Remove files or define a new '
             'output dir name to prevent accidental overwriting.')

def set_up_outdir(outdir, overwrite=False):
    '''Create outdir if it does not exist, and require overwrite if it does, so that there are '
    no stale files.'''
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise ValueError(f'The filename you designated as the output directory, {outdir}, '
                             'already exists and is not a directory.')
        if overwrite:
            print(f'\nWarning: overwriting existing output directory {outdir} because '
                  'overwrite_existing was set to True.\n')
            # Another process might remove/create concurrently; guard with try/except.
            try:
                shutil.rmtree(outdir)
            except FileNotFoundError:
                pass
            os.makedirs(outdir, exist_ok=True)
        else:
            # Allow empty dirs; error only if files are present.
            if any(os.scandir(outdir)):
                raise ValueError(f'The output directory {outdir} is not empty. Remove files or '
                                 'set overwrite_existing to True to prevent accidental overwriting.')
    else:
        parent_dir = os.path.dirname(outdir)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        os.makedirs(outdir, exist_ok=True)

def valid_database_subdir_format(input_dir):
    ''' Ensure that the pdb database dir has subdirs formatted similarly to the RCSB
    mirror format (see docs/database_generation.txt). '''
    
    def error_message():
        print('Database structure must be similar to the RCSB mirror format; '
              'see docs/database_generation.txt')
    
    # Input path must be a dir
    if not os.path.exists(input_dir):
        error_message()
        return False
    # Input dir must not be empty
    if not len(os.listdir(input_dir)) > 0: 
        error_message()
        return False
    # Ensure that there is at least one correctly formatted subdir
    for p in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, p)
        if os.path.isdir(subdir_path):
            # Subdirs should be the inner 2 characters of a 4-char pdb file name
            if len(p) == 2: 
                for pdb in os.listdir(subdir_path):
                    if pdb[1:3] == p:
                        return True
    # No valid subdir was found
    error_message()
    return False

vdw_radii = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 
            'Br': 1.85,'I': 1.98,}

def min_clash_dist(elem1, elem2, hbond_tol=0.6, general_tol=0.4):
    # Params: tolerance (Å) for decreasing the minimum allowable distance for atoms 
    # that may participate in hbonds (hbond_tol) or non-hbond interactions (general_tol).
    # Returns the minimum distance at which atoms are considered to clash.

    can_hbond = {frozenset(['N', 'O']), frozenset(['O', 'O']), frozenset(['N', 'N'])}

    if elem1 not in vdw_radii:
        print('Undefined vdw radius for element:', elem1)
        vdw1 = max(vdw_radii.values())
    else:
        vdw1 = vdw_radii[elem1]

    if elem2 not in vdw_radii:
        print('Undefined vdw radius for element:', elem2)
        vdw2 = max(vdw_radii.values())
    else:
        vdw2 = vdw_radii[elem2]

    min_dist = vdw1 + vdw2

    # Apply tolerance
    if frozenset([elem1, elem2]) in can_hbond:
        min_dist -= hbond_tol 
    else:
        min_dist -= general_tol 

    return min_dist

def get_atom_coords(prody_obj, atom_name):
    atom = prody_obj.select(f'name {atom_name}')
    coords = atom.getCoords()[0]
    return coords

def smiles_equiv(existingfrag, sub_smiles, check_atom_order):
    # check_atom_order requires that the order of elements is identical between the SMILES
    mol1 = Chem.MolFromSmarts(existingfrag)
    mol2 = Chem.MolFromSmarts(sub_smiles)
    mol1_has_mol2 = mol1.HasSubstructMatch(mol2) 
    mol2_has_mol1 = mol2.HasSubstructMatch(mol1)

    mol1_charge = Chem.GetFormalCharge(mol1)
    mol2_charge = Chem.GetFormalCharge(mol2)

    # is the element order correct?
    mol1_elements = [atom.GetSymbol() for atom in mol1.GetAtoms()]
    mol2_elements = [atom.GetSymbol() for atom in mol2.GetAtoms()]

    is_equivalent = (mol1.GetNumAtoms() == mol2.GetNumAtoms() and 
                        mol1_has_mol2 and mol2_has_mol1 and 
                        mol1_charge == mol2_charge)

    if check_atom_order:
        is_equivalent = (is_equivalent and mol1_elements == mol2_elements)

    return is_equivalent

def _find_resonance_O_groups(mol):
    """
    Find resonance X–O groups where:
      - A central atom X has at least one multiple bond to O, and
      - Either:
          * X has >= 3 O neighbors, or
          * X has >= 2 O neighbors AND (any O is anionic OR X is charged).

    All X–O bonds in such a group will be treated as resonance-equivalent
    (ignore single vs double).

    Returns: dict {bond_idx: (center_atomic_num, group_id)}
    """
    resonance_groups = {}
    next_group_id = 0

    # X candidates where multi-O resonance is common
    X_CANDIDATES = {6, 7, 8, 15, 16}  # C, N, O, P, S

    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z not in X_CANDIDATES:
            continue

        center_charge = atom.GetFormalCharge()
        oxy_bonds = []
        has_multiple_to_O = False
        has_anionic_O = False

        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 8:  # oxygen neighbor
                oxy_bonds.append(bond.GetIdx())
                if bond.GetBondTypeAsDouble() > 1.01:  # > single
                    has_multiple_to_O = True
                if other.GetFormalCharge() < 0:
                    has_anionic_O = True

        num_O = len(oxy_bonds)
        if not has_multiple_to_O or num_O < 2:
            continue

        # Decide if this X center should be treated as a resonance X–O center
        is_resonance_center = (num_O >= 3 or (num_O >= 2 and 
                              (has_anionic_O or center_charge != 0)))

        if not is_resonance_center:
            continue

        for bidx in oxy_bonds:
            resonance_groups[bidx] = (z, next_group_id)
        next_group_id += 1

    return resonance_groups

def identify_mol_symmetry(mol):
    """
    Identify atom symmetry using Weisfeiler–Lehman (WL) color refinement.
    - Atom labels: (atomic number, is_aromatic)  [charges ignored to allow
      resonance-equivalent atoms (e.g. carboxylate, nitrate) to become equivalent]
    - Bond labels:
        * Normal bonds: (0, bond_order, is_aromatic, is_conjugated)
        * Resonance X–O bonds: (1, center_atomic_num, group_id)

    Returns a list of lists, where atom indices in each sublist are equivalent.
    """
    atoms = list(mol.GetAtoms())
    num_atoms = mol.GetNumAtoms()

    labels = [(atom.GetAtomicNum(), atom.GetIsAromatic()) for atom in atoms]

    # Precompute resonance X–O bond groups
    resonance_groups = _find_resonance_O_groups(mol)

    def bond_featurizer(bond):
        bidx = bond.GetIdx()
        if bidx in resonance_groups:
            center_z, group_id = resonance_groups[bidx]
            # 1 = resonance X–O type
            return (1, center_z, group_id)

        # 0 = normal bond type
        return (0, float(bond.GetBondTypeAsDouble()), bond.GetIsAromatic(),
            bond.GetIsConjugated())

    # Precompute neighbors and bond labels
    nbrs = [[] for _ in range(num_atoms)]
    bond_labels = [[] for _ in range(num_atoms)]

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        b_label = bond_featurizer(bond)

        nbrs[i].append(j)
        nbrs[j].append(i)
        bond_labels[i].append(b_label)
        bond_labels[j].append(b_label)

    # WL refinement
    changed = True
    while changed:
        changed = False
        sigs = []

        for i in range(num_atoms):
            neigh_sig = []
            for j, b_lbl in zip(nbrs[i], bond_labels[i]):
                neigh_sig.append((b_lbl, labels[j]))

            neigh_sig.sort()
            sigs.append((labels[i], tuple(neigh_sig)))

        uniq = {}
        refined_labels = []
        next_id = 0
        for s in sigs:
            if s not in uniq:
                uniq[s] = next_id
                next_id += 1
            refined_labels.append(uniq[s])

        # Repeat until labels stop changing
        if refined_labels != labels:
            labels = refined_labels
            changed = True

    # Collect symmetry classes: atoms with the same final label are grouped
    groups = defaultdict(list)
    for i, lbl in enumerate(labels):
        groups[lbl].append(i)

    # Sort classes by size, then lexicographically by indices (for stability)
    return sorted(groups.values(), key=lambda g: (len(g), g))

def define_symmetry(smiles):
    # returns None if no symmetry, or symm_classes if symmetry
    mol_obj = Chem.MolFromSmarts(smiles)
    symm_groups = identify_mol_symmetry(mol_obj)
    has_symm = len(symm_groups) < mol_obj.GetNumAtoms()
    if not has_symm:
        return None

    # Step 1: assign raw labels based on group index
    max_index = max(max(sub) for sub in symm_groups)
    raw = [None] * (max_index + 1)
    for label, sub in enumerate(symm_groups):
        for idx in sub:
            raw[idx] = label

    # Step 2: renumber labels in order of first appearance (so sequence starts with 0)
    mapping = {}
    next_label = 0
    normalized = []
    for lab in raw:
        if lab not in mapping:
            mapping[lab] = next_label
            next_label += 1
        normalized.append(mapping[lab])

    return " ".join(map(str, normalized))