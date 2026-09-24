import itertools
import os
import re
import shutil
from rdkit import Chem
from rdkit.Chem import rdMolAlign as MA
import numpy as np

import hashlib

def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def _int_or_none(v):
    if v.lower() == 'none':
        return None
    return int(v)

def set_up_outdir(outdir, overwrite=False):
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise ValueError(f'[ERROR] The filename you designated as the output directory, {outdir}, '
                             'already exists and is not a directory.')
        if overwrite:
            print(f'[WARNING] Overwriting existing output directory {outdir} because '
                  'overwrite_existing was set to True.')
            try:
                shutil.rmtree(outdir)
            except FileNotFoundError:
                pass
            os.makedirs(outdir, exist_ok=True)
        else:
            with os.scandir(outdir) as entries:
                not_empty = any(entries)
            if not_empty:
                raise ValueError(f'[ERROR] The output directory {outdir} is not empty. Remove files or '
                                 'set overwrite_existing to True to prevent accidental overwriting.')
    else:
        parent_dir = os.path.dirname(outdir)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        os.makedirs(outdir, exist_ok=True)

def smiles_equiv(existingfrag, sub_smiles):
    return fragment_keys_equivalent(existingfrag, sub_smiles)

_TERMINAL_RESONANCE_ELEMENTS = frozenset({7, 8, 16})
_TERMINAL_RESONANCE_CENTERS = frozenset({5, 6, 7, 8, 15, 16, 17, 35, 53})

_QUERY_DEGREE_LINE = re.compile(r'^\s*AtomExplicitDegree (\d+) (!?)= val\s*$', re.M)
_QUERY_HCOUNT_LINE = re.compile(r'^\s*AtomHCount (\d+) (!?)= val\s*$', re.M)

def _query_primitive(atom, pattern):
    if not atom.HasQuery():
        return None
    description = atom.DescribeQuery()
    if 'RecursiveStructure' in description:
        return None
    found = pattern.findall(description)
    if not found:
        return None
    return tuple(sorted((int(value), bool(negated)) for value, negated in found))

def _query_allows_terminal(atom):
    degree = _query_primitive(atom, _QUERY_DEGREE_LINE)
    return degree is None or degree == ((1, False),)

def _has_deliberate_charge_assignment(mol, bond_indices):
    if len(bond_indices) < 4:
        return False
    if any(mol.GetBondWithIdx(i).GetBondTypeAsDouble() != 1.0 for i in bond_indices):
        return False

    charges = set()
    for bond_idx in bond_indices:
        bond = mol.GetBondWithIdx(bond_idx)
        for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
            if atom.GetDegree() == 1:
                charges.add(atom.GetFormalCharge())
    return len(charges) > 1

def _find_resonance_terminal_groups(mol):
    resonance_groups = {}
    next_group_id = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _TERMINAL_RESONANCE_CENTERS:
            continue

        terminal_bonds_by_element = {}
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            atomic_num = other.GetAtomicNum()
            if atomic_num not in _TERMINAL_RESONANCE_ELEMENTS:
                continue
            if other.GetIsAromatic() or other.GetDegree() != 1:
                continue
            if not _query_allows_terminal(other):
                continue
            terminal_bonds_by_element.setdefault(atomic_num, []).append(bond.GetIdx())

        for atomic_num, bond_indices in sorted(terminal_bonds_by_element.items()):
            if len(bond_indices) < 2:
                continue
            if _has_deliberate_charge_assignment(mol, bond_indices):
                continue
            for bond_idx in bond_indices:
                resonance_groups[bond_idx] = (atomic_num, next_group_id)
            next_group_id += 1

    return resonance_groups

def _find_resonance_center_groups(mol, center_z, neighbor_z, charge_rule):
    if charge_rule not in (None, 'positive', 'nonzero'):
        raise ValueError(f"unknown charge_rule {charge_rule!r}")

    resonance_groups = {}
    next_group_id = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != center_z or atom.GetIsAromatic():
            continue

        neighbor_bonds = []
        has_double_bond = False
        for bond in atom.GetBonds():
            if bond.GetOtherAtom(atom).GetAtomicNum() != neighbor_z:
                continue
            neighbor_bonds.append(bond.GetIdx())
            if bond.GetBondTypeAsDouble() == 2.0:
                has_double_bond = True

        if len(neighbor_bonds) < 2:
            continue
        if not has_double_bond:
            charge = atom.GetFormalCharge()
            if charge_rule is None:
                continue
            if charge_rule == 'positive' and charge <= 0:
                continue
            if charge_rule == 'nonzero' and charge == 0:
                continue

        for bond_idx in neighbor_bonds:
            resonance_groups[bond_idx] = (center_z, next_group_id)
        next_group_id += 1

    return resonance_groups

def _find_resonance_CN_groups(mol):
    return _find_resonance_center_groups(mol, 6, 7, 'positive')

def _find_resonance_NO_groups(mol):
    return _find_resonance_center_groups(mol, 7, 8, 'nonzero')

def _find_resonance_SN_groups(mol):
    return _find_resonance_center_groups(mol, 16, 7, None)

def _find_resonance_NN_groups(mol):
    nitrogen_indices = {
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic()
    }
    visited = set()
    resonance_groups = {}
    next_group_id = 0

    for start in sorted(nitrogen_indices):
        if start in visited:
            continue

        component = set()
        component_bonds = set()
        stack = [start]
        visited.add(start)
        while stack:
            atom_idx = stack.pop()
            component.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                neighbor_idx = bond.GetOtherAtomIdx(atom_idx)
                if neighbor_idx not in nitrogen_indices:
                    continue
                component_bonds.add(bond.GetIdx())
                if neighbor_idx not in visited:
                    visited.add(neighbor_idx)
                    stack.append(neighbor_idx)

        has_double_bond = any(
            mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() == 2.0
            for bond_idx in component_bonds
        )
        if len(component) < 3 or not has_double_bond:
            continue

        for bond_idx in component_bonds:
            resonance_groups[bond_idx] = (7, next_group_id)
        next_group_id += 1

    return resonance_groups

def _find_resonance_aromatic_N_atoms(mol):
    aromatic_atom_indices = {
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()
    }
    visited = set()
    resonance_atoms = set()

    for start in sorted(aromatic_atom_indices):
        if start in visited:
            continue

        component = set()
        stack = [start]
        visited.add(start)
        while stack:
            atom_idx = stack.pop()
            component.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                if not bond.GetIsAromatic():
                    continue
                neighbor_idx = bond.GetOtherAtomIdx(atom_idx)
                if neighbor_idx not in aromatic_atom_indices:
                    continue
                if neighbor_idx in visited:
                    continue
                visited.add(neighbor_idx)
                stack.append(neighbor_idx)

        nitrogen_atoms = {
            atom_idx for atom_idx in component
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7
        }
        if len(nitrogen_atoms) >= 2:
            resonance_atoms.update(nitrogen_atoms)

    return resonance_atoms

def mol_from_fragment(fragment):
    return init_query_ring_info(Chem.MolFromSmarts(fragment))

_BRACKET_ATOM_INNER = re.compile(r'\[([^\]]+)\]')
_RING_TOKEN = re.compile(r'\A(?:!R|R|r\d+)\Z')
_DEGREE_TOKEN = re.compile(r'\AD(\d+)\Z')
_H_TOKEN = re.compile(r'\A!?H0\Z')
_RING_ISOTOPE = {'!R': 1, 'R': 2}
_H_ISOTOPE = {'H0': 1, '!H0': 2}
_MAX_ISOTOPE = 65535

def split_bracket_annotations(inner, fragment=None):
    atom_part, *tokens = inner.split(';')
    ring_token = degree = h_token = None
    for token in tokens:
        if _RING_TOKEN.match(token):
            seen, name = ring_token, 'ring context'
        elif _DEGREE_TOKEN.match(token):
            seen, name = degree, 'heavy degree'
        elif _H_TOKEN.match(token):
            seen, name = h_token, 'hydrogen flag'
        else:
            raise ValueError(
                f'unrecognised annotation {token!r} in {inner!r} of '
                f'{fragment!r}; key comparison cannot encode it')
        if seen is not None:
            raise ValueError(
                f'conflicting annotations in {inner!r} of {fragment!r}: '
                f'a second {name} token {token!r}')
        if _RING_TOKEN.match(token):
            ring_token = token
        elif _DEGREE_TOKEN.match(token):
            degree = int(_DEGREE_TOKEN.match(token).group(1))
        else:
            h_token = token
    if degree is not None and h_token is not None:
        raise ValueError(
            f'{inner!r} of {fragment!r} carries both a degree and a hydrogen '
            'flag; they share one isotope field and are exclusive by construction')
    return atom_part, ring_token, degree, h_token

def _annotations_as_isotope(fragment):
    def _encode(match):
        inner = match.group(1)
        atom_part, ring_token, degree, h_token = split_bracket_annotations(
            inner, fragment)
        ring_code = 0 if ring_token is None else _RING_ISOTOPE.get(
            ring_token, 10 + int(ring_token[1:]) if ring_token[0] == 'r' else 0)
        dh_code = (_H_ISOTOPE[h_token] if h_token is not None
                   else 10 + degree if degree is not None else 0)
        isotope = ring_code + 100 * dh_code
        if isotope > _MAX_ISOTOPE:
            raise ValueError(f'annotation isotope {isotope} for {inner!r} of '
                             f'{fragment!r} exceeds RDKit\'s 16-bit isotope field')
        return f'[{isotope or ""}{atom_part}]'
    return _BRACKET_ATOM_INNER.sub(_encode, fragment)

def fragment_keys_equivalent(key_a, key_b):
    return fragment_query_mols_equivalent(fragment_key_query_mol(key_a),
                                          fragment_key_query_mol(key_b))

def fragment_key_query_mol(key):
    return init_query_ring_info(Chem.MolFromSmarts(_annotations_as_isotope(key)))

def fragment_query_mols_equivalent(mol_a, mol_b):
    if mol_a is None or mol_b is None:
        return False
    if mol_a.GetNumAtoms() != mol_b.GetNumAtoms():
        return False
    if Chem.GetFormalCharge(mol_a) != Chem.GetFormalCharge(mol_b):
        return False
    return bool(mol_a.HasSubstructMatch(mol_b) and mol_b.HasSubstructMatch(mol_a))

def init_query_ring_info(mol):
    if mol is None:
        return mol
    Chem.FastFindRings(mol)
    return mol

_AROMATIC_H_ATOM = re.compile(r'\[([a-z][a-z]?)(H\d?)([^\]]*)\]')

def aromatic_h_constrained_atoms(fragment):
    return [match.group(0) for match in _AROMATIC_H_ATOM.finditer(fragment or '')]

def _require_hydrogen_free(mol, context):
    if any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        raise ValueError(
            f"{context} requires a hydrogen-free molecular graph; strip explicit "
            "hydrogen atom nodes upstream"
        )

def _automorphism_graph(mol):
    num_atoms = mol.GetNumAtoms()
    resonance_bonds = _find_resonance_terminal_groups(mol)
    resonance_bonds.update(_find_resonance_CN_groups(mol))
    resonance_bonds.update(_find_resonance_NO_groups(mol))
    resonance_bonds.update(_find_resonance_SN_groups(mol))
    resonance_bonds.update(_find_resonance_NN_groups(mol))
    aromatic_resonance_atoms = _find_resonance_aromatic_N_atoms(mol)
    resonance_atoms = set(aromatic_resonance_atoms)
    for bond_idx in resonance_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        resonance_atoms.add(bond.GetBeginAtomIdx())
        resonance_atoms.add(bond.GetEndAtomIdx())

    atom_labels = []
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        if atom_idx in aromatic_resonance_atoms:
            query_label = ("aromatic_N_resonance",)
        elif atom.HasQuery() and atom_idx not in resonance_atoms:
            query_label = atom.GetSmarts()
        else:
            query_label = None
        degree_label = _query_primitive(atom, _QUERY_DEGREE_LINE)
        h_label = (_query_primitive(atom, _QUERY_HCOUNT_LINE)
                   if atom.GetAtomicNum() == 6 else None)
        atom_labels.append((
            atom.GetAtomicNum(),
            atom.GetIsAromatic(),
            atom.GetIsotope(),
            0 if atom.GetIdx() in resonance_atoms else atom.GetFormalCharge(),
            atom.GetNumRadicalElectrons(),
            int(atom.GetChiralTag()),
            query_label,
            degree_label,
            h_label,))

    adjacency = [[None] * num_atoms for _ in range(num_atoms)]
    for bond in mol.GetBonds():
        if bond.GetIdx() in resonance_bonds:
            bond_label = (1,)
        else:
            bond_label = (
                0,
                float(bond.GetBondTypeAsDouble()),
                bond.GetIsAromatic(),
                bond.GetIsConjugated(),
                int(bond.GetStereo()),
                int(bond.GetBondDir()),
                bond.GetSmarts() if bond.HasQuery() else None,
            )
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        adjacency[i][j] = bond_label
        adjacency[j][i] = bond_label

    return atom_labels, adjacency

def _wl_colors(atom_labels, adjacency):
    num_atoms = len(atom_labels)
    colors = list(atom_labels)

    for _ in range(max(1, num_atoms)):
        signatures = []
        for i in range(num_atoms):
            neighbors = [
                (adjacency[i][j], colors[j])
                for j in range(num_atoms)
                if adjacency[i][j] is not None
            ]
            neighbors.sort(key=repr)
            signatures.append((colors[i], tuple(neighbors)))

        signature_ids = {}
        refined = []
        for signature in signatures:
            if signature not in signature_ids:
                signature_ids[signature] = len(signature_ids)
            refined.append(signature_ids[signature])

        if len(set(refined)) == len(set(colors)):
            return refined
        colors = refined

    return colors

def validate_atom_permutations(permutations, num_atoms=None):
    if permutations is None:
        return None

    normalized = []
    seen = set()
    inferred_num_atoms = num_atoms
    for permutation in permutations:
        try:
            perm = tuple(int(i) for i in permutation)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid atom permutation: {permutation!r}") from exc
        if inferred_num_atoms is None:
            inferred_num_atoms = len(perm)
        if len(perm) != inferred_num_atoms:
            raise ValueError(
                f"Atom permutation has {len(perm)} entries; expected {inferred_num_atoms}"
            )
        if set(perm) != set(range(inferred_num_atoms)):
            raise ValueError(f"Not a permutation of 0..{inferred_num_atoms - 1}: {perm}")
        if perm not in seen:
            seen.add(perm)
            normalized.append(perm)

    if not normalized:
        raise ValueError("At least one atom permutation is required")

    identity = tuple(range(inferred_num_atoms))
    if identity not in seen:
        raise ValueError("Atom permutations must include the identity mapping")

    normalized.sort()
    normalized.remove(identity)
    normalized.insert(0, identity)
    return tuple(normalized)

def group_preserving_permutations(labels):
    labels = list(labels)
    n = len(labels)

    groups = {}
    for i, label in enumerate(labels):
        groups.setdefault(label, []).append(i)
    groups_list = list(groups.values())
    group_iperms = [list(itertools.permutations(g)) for g in groups_list]

    result = []
    for combo in itertools.product(*group_iperms):
        perm = [None] * n
        for pos_list, perm_of_pos in zip(groups_list, combo):
            for pos, orig in zip(pos_list, perm_of_pos):
                perm[pos] = orig
        result.append(perm)
    return result

def identify_mol_automorphisms(mol, max_automorphisms=10000):
    if mol is None:
        raise ValueError("Cannot identify automorphisms of a null molecule")
    _require_hydrogen_free(mol, "Automorphism enumeration")
    if not isinstance(max_automorphisms, int) or isinstance(max_automorphisms, bool) \
            or max_automorphisms < 1:
        raise ValueError("max_automorphisms must be a positive integer")

    atom_labels, adjacency = _automorphism_graph(mol)
    num_atoms = len(atom_labels)
    if num_atoms == 0:
        return ((),)

    colors = _wl_colors(atom_labels, adjacency)
    color_members = {}
    for atom_idx, color in enumerate(colors):
        color_members.setdefault(color, []).append(atom_idx)

    order = sorted(
        range(num_atoms),
        key=lambda i: (len(color_members[colors[i]]),
                       -sum(edge is not None for edge in adjacency[i]), i),
    )
    mapping = [-1] * num_atoms
    used_targets = [False] * num_atoms
    automorphisms = []

    def search(depth):
        if len(automorphisms) > max_automorphisms:
            return
        if depth == num_atoms:
            automorphisms.append(tuple(mapping))
            return

        source = order[depth]
        candidates = sorted(
            color_members[colors[source]], key=lambda target: target != source
        )
        for target in candidates:
            if used_targets[target]:
                continue

            preserves_graph = True
            for other_source, other_target in enumerate(mapping):
                if other_target < 0:
                    continue
                if adjacency[source][other_source] != adjacency[target][other_target]:
                    preserves_graph = False
                    break
            if not preserves_graph:
                continue

            mapping[source] = target
            used_targets[target] = True
            search(depth + 1)
            used_targets[target] = False
            mapping[source] = -1
            if len(automorphisms) > max_automorphisms:
                return

    search(0)
    if len(automorphisms) > max_automorphisms:
        raise ValueError(
            f"Molecule has more than max_automorphisms={max_automorphisms}; "
            "refusing to use an incomplete symmetry mapping set"
        )
    return validate_atom_permutations(automorphisms, num_atoms)

def _proper_kabsch_rotations(H):
    if not np.isfinite(H).all():
        raise ValueError("Non-finite Kabsch cross-covariance matrix")

    U, _, Vt = np.linalg.svd(H, full_matrices=False)
    UVt = np.matmul(U, Vt)
    det_uvt = _det3(UVt)
    if not np.isfinite(det_uvt).all():
        raise ValueError("Non-finite determinant from Kabsch SVD")

    d = np.where(det_uvt < 0.0, -1.0, 1.0).astype(H.dtype, copy=False)
    D = np.zeros_like(H)
    D[..., 0, 0] = 1.0
    D[..., 1, 1] = 1.0
    D[..., 2, 2] = d
    return np.matmul(U, np.matmul(D, Vt))

def _det3(H):
    a, b, c = H[..., 0, 0], H[..., 0, 1], H[..., 0, 2]
    d, e, f = H[..., 1, 0], H[..., 1, 1], H[..., 1, 2]
    g, h, i = H[..., 2, 0], H[..., 2, 1], H[..., 2, 2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

def _singular_values_3x3(H):
    """Singular values of batched 3x3 matrices without per-matrix LAPACK calls."""
    A = np.matmul(np.transpose(H, (0, 2, 1)), H)
    q = (A[:, 0, 0] + A[:, 1, 1] + A[:, 2, 2]) / 3.0
    A[:, 0, 0] -= q
    A[:, 1, 1] -= q
    A[:, 2, 2] -= q
    p = np.sqrt(np.maximum(np.einsum('nij,nij->n', A, A) / 6.0, 0.0))
    A /= np.maximum(p, np.finfo(H.dtype).tiny)[:, None, None]
    phi = np.arccos(np.clip(_det3(A) * 0.5, -1.0, 1.0)) / 3.0
    eig1 = q + 2.0 * p * np.cos(phi)
    cross01 = np.cross(H[:, 0], H[:, 1])
    cross02 = np.cross(H[:, 0], H[:, 2])
    cross12 = np.cross(H[:, 1], H[:, 2])
    eig23_sum = (np.einsum('ni,ni->n', cross01, cross01)
                 + np.einsum('ni,ni->n', cross02, cross02)
                 + np.einsum('ni,ni->n', cross12, cross12))
    eig1 = np.maximum(eig1, 0.0)
    eig23_product = _det3(H) ** 2 / np.maximum(eig1, np.finfo(H.dtype).tiny)
    eig23_sum = np.maximum((eig23_sum - eig23_product)
                           / np.maximum(eig1, np.finfo(H.dtype).tiny), 0.0)
    disc = np.sqrt(np.maximum(eig23_sum ** 2 - 4.0 * eig23_product, 0.0))
    eig2 = 0.5 * (eig23_sum + disc)
    eig3 = eig23_product / np.maximum(eig2, np.finfo(H.dtype).tiny)
    s1 = np.sqrt(eig1)
    s2 = np.sqrt(eig2)
    s3 = np.minimum(np.sqrt(eig3), s2)
    return np.stack((s1, s2, s3), axis=1)

def kabsch_ssd(X, Y, chunk_size=30000):
    X = np.asarray(X, dtype=np.float32)
    Y = np.asarray(Y, dtype=np.float32)

    if Y.ndim != 3:
        raise ValueError(f"kabsch_ssd: Y must have ndim=3, got Y.ndim={Y.ndim}, Y.shape={Y.shape}")
    if X.ndim not in (2, 3):
        raise ValueError(f"kabsch_ssd: X must have ndim 2 or 3, got X.ndim={X.ndim}")
    if X.shape[-2] != Y.shape[1]:
        raise ValueError(
            f"kabsch_ssd: point-count mismatch, X has {X.shape[-2]} points per structure, "
            f"Y has {Y.shape[1]}")
    if X.ndim == 3 and X.shape[0] != Y.shape[0]:
        raise ValueError(
            f"kabsch_ssd: batch mismatch, X has {X.shape[0]} structures, Y has {Y.shape[0]}")
    if not np.isfinite(X).all() or not np.isfinite(Y).all():
        raise ValueError(
            "Non-finite values detected in kabsch_ssd input. "
            f"X.shape={X.shape}; Y.shape={Y.shape}")

    M = Y.shape[0]
    if M == 0:
        return np.empty((0,), dtype=np.float64)

    n_pts = Y.shape[1]
    inv_n = 1.0 / n_pts
    fixed_X = X.ndim == 2
    if fixed_X:
        Xc = (X - np.add.reduce(X, axis=0) * inv_n).astype(np.float64)
        XcT = np.ascontiguousarray(Xc.T)
        x_norm = float(np.add.reduce(np.add.reduce(Xc * Xc)))

    chunks = []
    for start in range(0, M, chunk_size):
        Yc = Y[start:start + chunk_size].astype(np.float64)
        Yc -= np.add.reduce(Yc, axis=1)[:, None, :] * inv_n

        if fixed_X:
            H = np.matmul(XcT[None, :, :], Yc)
            xn = x_norm
        else:
            Xc_b = X[start:start + chunk_size].astype(np.float64)
            Xc_b -= np.add.reduce(Xc_b, axis=1)[:, None, :] * inv_n
            H = np.matmul(np.transpose(Xc_b, (0, 2, 1)), Yc)
            xn = np.add.reduce(np.add.reduce(Xc_b * Xc_b, axis=2), axis=1)

        sv = _singular_values_3x3(H)
        d = np.where(_det3(H) < 0.0, -1.0, 1.0)
        y_norm = np.add.reduce(np.add.reduce(Yc * Yc, axis=2), axis=1)
        trace = sv[:, 0] + sv[:, 1] + d * sv[:, 2]
        chunks.append(np.maximum(xn + y_norm - 2.0 * trace, 0.0))

    return np.concatenate(chunks, axis=0) if len(chunks) > 1 else chunks[0]

def kabsch(X, Y, chunk_size=30000):
    X = np.asarray(X, dtype=np.float32)
    Y = np.asarray(Y, dtype=np.float32)

    if Y.ndim != 3:
        raise ValueError(f"kabsch: Y must have ndim=3, got Y.ndim={Y.ndim}, Y.shape={Y.shape}")

    if X.ndim not in (2, 3):
        raise ValueError(f"kabsch: X must have ndim 2 or 3, got X.ndim={X.ndim}")
    if X.shape[-2] != Y.shape[1]:
        raise ValueError(
            f"kabsch: point-count mismatch, X has {X.shape[-2]} points per structure, "
            f"Y has {Y.shape[1]}")
    if X.ndim == 3 and X.shape[0] != Y.shape[0]:
        raise ValueError(
            f"kabsch: batch mismatch, X has {X.shape[0]} structures, Y has {Y.shape[0]}")

    x_finite = np.isfinite(X)
    y_finite = np.isfinite(Y)
    if not x_finite.all() or not y_finite.all():
        x_invalid = np.argwhere(~x_finite)
        y_invalid = np.argwhere(~y_finite)
        x_msg = ("none" if x_invalid.size == 0
                 else f"first invalid X value at {tuple(x_invalid[0])}")
        y_msg = ("none" if y_invalid.size == 0
                 else f"first invalid Y value at {tuple(y_invalid[0])}")
        raise ValueError(
            "Non-finite values detected in kabsch input. "
            f"{x_msg}; {y_msg}; "
            f"X.shape={X.shape}; Y.shape={Y.shape}"
        )

    M = Y.shape[0]
    R_chunks = []
    t_chunks = []
    ssd_chunks = []

    if X.ndim == 2:
        Xbar = X.mean(axis=0, keepdims=True)
        Xc = X - Xbar
        XcT = Xc.T

        for start in range(0, M, chunk_size):
            stop = min(start + chunk_size, M)
            Y_chunk = Y[start:stop]

            Ybar = Y_chunk.mean(axis=1, keepdims=True)
            Yc = Y_chunk - Ybar

            H = np.matmul(XcT[None, :, :], Yc)

            R = _proper_kabsch_rotations(H).astype(np.float32, copy=False)

            XRbar = np.matmul(Xbar[None, :, :], R)[:, 0, :]
            t = (Ybar[:, 0, :] - XRbar).astype(np.float32, copy=False)

            XR = np.matmul(Xc[None, :, :], R)
            diff = XR - Yc
            ssd = np.sum(diff * diff, axis=(1, 2), dtype=np.float64)

            R_chunks.append(R)
            t_chunks.append(t)
            ssd_chunks.append(ssd)

    else:
        for start in range(0, M, chunk_size):
            stop = min(start + chunk_size, M)
            X_chunk = X[start:stop]
            Y_chunk = Y[start:stop]

            Xbar = X_chunk.mean(axis=1, keepdims=True)
            Ybar = Y_chunk.mean(axis=1, keepdims=True)

            Xc = X_chunk - Xbar
            Yc = Y_chunk - Ybar

            H = np.matmul(np.transpose(Xc, (0, 2, 1)), Yc)

            R = _proper_kabsch_rotations(H).astype(np.float32, copy=False)
            t = (Ybar - np.matmul(Xbar, R)).reshape(-1, 3).astype(np.float32, copy=False)

            diff = np.matmul(Xc, R) - Yc
            ssd = np.sum(diff * diff, axis=(1, 2), dtype=np.float64)

            R_chunks.append(R)
            t_chunks.append(t)
            ssd_chunks.append(ssd)

    if R_chunks:
        R = np.concatenate(R_chunks, axis=0)
        t = np.concatenate(t_chunks, axis=0)
        ssd = np.concatenate(ssd_chunks, axis=0)
    else:
        R = np.empty((0, 3, 3), dtype=np.float32)
        t = np.empty((0, 3), dtype=np.float32)
        ssd = np.empty((0,), dtype=np.float64)

    return R, t, ssd

def convert_time_elapsed(seconds):
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = round(seconds % 60, 2)
    return h, m, s

def best_inplace_symmetry_rmsd(ref_mol, query_mol):
    if ref_mol.GetNumAtoms() != query_mol.GetNumAtoms():
        raise ValueError(f"[ERROR] Atom count mismatch: ref={ref_mol.GetNumAtoms()}, "
            f"query={query_mol.GetNumAtoms()}")
    _require_hydrogen_free(ref_mol, "In-place symmetry RMSD")
    _require_hydrogen_free(query_mol, "In-place symmetry RMSD")
    return MA.CalcRMS(
        query_mol,
        ref_mol,
        maxMatches=0,
        symmetrizeConjugatedTerminalGroups=True,
    )

_SMILES_TO_FILENAME = {'/': '_fs_', '\\': '_bs_'}
_SMILES_TO_JOB_NAME = {'#': '_tp_', '/': '_fs_', '\\': '_bs_'}

def normalize_rmsd(num_atoms, atoms):
    if atoms == 'flankbb':
        max_threshold, min_threshold = 1.5, 0.5
    elif atoms == 'cgvdmbb':
        max_threshold, min_threshold = 1.0, 0.5
    else:
        raise ValueError(f"Unknown atom set for normalize_rmsd: {atoms}")
    min_atoms, max_atoms = 8, 15
    if num_atoms < min_atoms:
        return min_threshold
    if num_atoms > max_atoms:
        return max_threshold
    return min_threshold + (num_atoms - min_atoms) / (max_atoms - min_atoms) * (max_threshold - min_threshold)

_SMILES_ATOM_RE = re.compile(
    r"\[\d*(?P<bracket>[A-Za-z][a-z]?|\*)"
    r"|(?P<organic>Br|Cl|[BCNOFPSI]|[bcnosp])")

def extract_elements(smiles: str):
    return [m.group('bracket') or m.group('organic')
            for m in _SMILES_ATOM_RE.finditer(smiles)]

def smiles_to_filename(smiles):
    return ''.join(_SMILES_TO_FILENAME.get(c, c) for c in smiles)

def filename_to_smiles(name):
    for ch, enc in _SMILES_TO_FILENAME.items():
        name = name.replace(enc, ch)
    return name

def smiles_to_job_name(smiles):
    return ''.join(_SMILES_TO_JOB_NAME.get(c, c) for c in smiles)

def job_name_to_filename(job_name):
    for ch, enc in _SMILES_TO_JOB_NAME.items():
        job_name = job_name.replace(enc, ch)
    return smiles_to_filename(job_name)
