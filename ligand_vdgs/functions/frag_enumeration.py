import itertools
import json
import re
from rdkit import Chem

ANNOTATION_PROP = '_vdgAtomAnnotation'

def _no_heteroatoms(sub):
    heavy = [a for a in sub.GetAtoms() if a.GetAtomicNum() != 1]
    return not heavy or all(a.GetAtomicNum() == 6 for a in heavy)

def induced_submol(mol, atom_indices, annotations):
    editable = Chem.RWMol()
    new_index = {}
    for new_idx, orig_idx in enumerate(atom_indices):
        editable.AddAtom(Chem.Atom(mol.GetAtomWithIdx(orig_idx)))
        new_index[orig_idx] = new_idx
    for orig_i, orig_j in itertools.combinations(atom_indices, 2):
        bond = mol.GetBondBetweenAtoms(orig_i, orig_j)
        if bond is not None:
            editable.AddBond(new_index[orig_i], new_index[orig_j], bond.GetBondType())
    sub = editable.GetMol()
    for orig_idx in atom_indices:
        sub.GetAtomWithIdx(new_index[orig_idx]).SetProp(ANNOTATION_PROP, annotations[orig_idx])
    return sub

def connected_atom_subsets(mol, min_frag_size, max_frag_size):
    if max_frag_size < 2:
        return []
    subsets = set()
    for bond_subgraphs in Chem.FindAllSubgraphsOfLengthMToN(mol, 1, max_frag_size - 1):
        for bond_indices in bond_subgraphs:
            atoms = set()
            for bond_idx in bond_indices:
                bond = mol.GetBondWithIdx(bond_idx)
                atoms.add(bond.GetBeginAtomIdx())
                atoms.add(bond.GetEndAtomIdx())
            if min_frag_size <= len(atoms) <= max_frag_size:
                subsets.add(tuple(sorted(atoms)))
    return sorted(subsets)

def enumerate_induced_fragments(mol, min_frag_size=4, max_frag_size=5, annotations=None):
    if annotations is None:
        annotations = parent_atom_annotations(mol)
    fragments = {}
    for atom_indices in connected_atom_subsets(mol, min_frag_size, max_frag_size):
        sub = induced_submol(mol, atom_indices, annotations)
        if _no_heteroatoms(sub):
            continue
        result = manually_remove_Hs(sub, 'single')
        if result is None:
            continue
        fragments.setdefault(result[1], []).append(atom_indices)
    return fragments

def get_fragments(bond_radius, mol, min_frag_size=4, max_frag_size=5, quiet=True):
    filtered_frags = {}
    for orig_sub, orig_mol_inds in (s for rad in range(1, bond_radius + 1)
            for s in fragment_on_bond_d(mol, rad, min_frag_size, max_frag_size)):
        results = manually_remove_Hs(orig_sub, return_single_mol_or_perms='perms')
        if results is None:
            continue
        sub_perms_mols_inds, sub_smiles = results
        orig_mol_inds_key = tuple(orig_mol_inds)
        seen = {(existing_inds, tuple(existing_orig))
                for existing_sub, existing_inds, existing_orig in filtered_frags.get(sub_smiles, [])}
        for sub, perm_inds in sub_perms_mols_inds:
            if _no_heteroatoms(sub) or (perm_inds, orig_mol_inds_key) in seen:
                continue
            seen.add((perm_inds, orig_mol_inds_key))
            filtered_frags.setdefault(sub_smiles, []).append((sub, perm_inds, orig_mol_inds))

    grouped_frags = {}
    for sub_smiles, substruct_data in filtered_frags.items():
        groups = group_lig_sites_by_overlap(substruct_data)
        grouped_frags[sub_smiles] = groups
        if not quiet:
            print(f"Fragment: {sub_smiles}, # perms: {len(substruct_data)}, "
                  f"# sites: {len(groups)}", flush=True)
    return grouped_frags

def group_lig_sites_by_overlap(data, key_index=2, threshold=0.5):
    sets = [set(item[key_index]) for item in data]
    n = len(data)
    parent = list(range(n))
    rank = [0] * n
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
    for i in range(n):
        for j in range(i + 1, n):
            if len(sets[i] & sets[j]) >= threshold * max(len(sets[i]), len(sets[j])):
                union(i, j)
    roots = [find(i) for i in range(n)]
    return [[d for d, r in zip(data, roots) if r == root] for root in dict.fromkeys(roots)]

def ring_query_for_atom(atom):
    if atom.GetIsAromatic():
        return None
    if not atom.IsInRing():
        return '!R'
    return f'r{min(atom.GetOwningMol().GetRingInfo().AtomRingSizes(atom.GetIdx()))}'

def atom_annotations(atom):
    if atom.GetAtomicNum() == 1:
        return ''
    parts = []
    ring = ring_query_for_atom(atom)
    if ring is not None:
        parts.append(ring)
    if atom.GetAtomicNum() == 6:
        parts.append('H0' if atom.GetTotalNumHs(includeNeighbors=True) == 0 else '!H0')
    else:
        parts.append(f'D{sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1)}')
    return ';'.join(parts)

def parent_atom_annotations(mol):
    return [atom_annotations(atom) for atom in mol.GetAtoms()]

_BRACKET_ATOM = re.compile(r'\[([^\]]+)\]')
_BRACKET_ATOM_LOOSE = re.compile(r'\[[^\]]*\]')
_EXPLICIT_H = re.compile(r'H(?![a-z])\d*')
_NON_LETTER = re.compile(r'[^A-Za-z]')

def _remove_bracket_hydrogens_in_smiles(smiles, annotations=None):
    out, pos = [], 0
    for idx, m in enumerate(_BRACKET_ATOM.finditer(smiles)):
        out.append(smiles[pos:m.start()])
        pos = m.end()
        inner = m.group(1)
        if inner == 'H':
            continue
        cleaned = _EXPLICIT_H.sub('', inner)
        annotation = None if annotations is None else annotations[idx]
        if annotation:
            out.append(f'[{cleaned};{annotation}]')
        elif _NON_LETTER.search(cleaned) or len(cleaned) > 1:
            out.append(f'[{cleaned}]')
        else:
            out.append(cleaned)
    out.append(smiles[pos:])
    return ''.join(out)

class _AnnotationUnmappable(Exception):
    pass

def _annotations_in_written_order(mol, smiles_w_H):
    if not mol.GetNumAtoms() or not mol.GetAtomWithIdx(0).HasProp(ANNOTATION_PROP):
        return None
    order = mol.GetPropsAsDict(includePrivate=True, includeComputed=True).get('_smilesAtomOutputOrder')
    if order is None:
        raise _AnnotationUnmappable('MolToSmiles did not record _smilesAtomOutputOrder')
    order = list(json.loads(order)) if isinstance(order, str) else list(order)
    if len(order) != mol.GetNumAtoms():
        raise _AnnotationUnmappable(
            f'output order covers {len(order)} of {mol.GetNumAtoms()} atoms')
    n_brackets = len(_BRACKET_ATOM_LOOSE.findall(smiles_w_H))
    if n_brackets != len(order):
        raise _AnnotationUnmappable(
            f'{n_brackets} bracket atoms in {smiles_w_H!r} for {len(order)} atoms')
    queries = []
    for atom_idx in order:
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.HasProp(ANNOTATION_PROP):
            raise _AnnotationUnmappable('fragment is only partially annotated')
        queries.append(atom.GetProp(ANNOTATION_PROP) or None)
    return queries

_MAX_CG_ATOM_PERMS = 1000

def manually_remove_Hs(orig_substruct, return_single_mol_or_perms):
    if return_single_mol_or_perms not in ('single', 'perms'):
        raise ValueError("return_single_mol_or_perms must be 'single' or 'perms', "
                         f"got {return_single_mol_or_perms!r}")
    if orig_substruct is None:
        return None

    try:
        substruct = Chem.RemoveAllHs(orig_substruct, sanitize=False)
        substruct.UpdatePropertyCache(strict=False)
    except Exception as e:
        print(f'[WARNING] manually_remove_Hs: RemoveAllHs failed ({e}); '
              'dropping fragment.', flush=True)
        return None
    if substruct.GetNumAtoms() == 0:
        return None

    try:
        smiles_w_H = Chem.MolToSmiles(substruct, allHsExplicit=True, isomericSmiles=False)
    except Exception as e:
        print(f'[WARNING] manually_remove_Hs: MolToSmiles failed on an unsanitized '
              f'fragment ({e}); dropping it.', flush=True)
        return None

    smiles_no_Hs = _remove_bracket_hydrogens_in_smiles(smiles_w_H)
    try:
        annotations = _annotations_in_written_order(substruct, smiles_w_H)
    except _AnnotationUnmappable as e:
        print(f'[WARNING] manually_remove_Hs: ring annotation for {smiles_no_Hs!r} '
              f'could not be mapped onto the written SMILES ({e}); dropping '
              'fragment.', flush=True)
        return None
    annotated_key = (smiles_no_Hs if annotations is None
                     else _remove_bracket_hydrogens_in_smiles(smiles_w_H, annotations))

    mol_from_smarts = Chem.MolFromSmarts(smiles_no_Hs)
    if mol_from_smarts is None or mol_from_smarts.GetNumAtoms() == 0:
        print(f'[WARNING] manually_remove_Hs: H-stripped SMARTS {smiles_no_Hs!r} is '
              'unusable; dropping fragment.', flush=True)
        return None
    try:
        matches_to_map_to_sub = substruct.GetSubstructMatches(
            mol_from_smarts, uniquify=False,
            maxMatches=1 if return_single_mol_or_perms == 'single' else _MAX_CG_ATOM_PERMS)
    except Exception as e:
        print(f'[WARNING] manually_remove_Hs: substructure match failed for '
              f'{smiles_no_Hs!r} ({e}); dropping fragment.', flush=True)
        return None
    if (return_single_mol_or_perms == 'perms'
            and len(matches_to_map_to_sub) >= _MAX_CG_ATOM_PERMS):
        print(f'[WARNING] manually_remove_Hs: {smiles_no_Hs!r} hit the '
              f'{_MAX_CG_ATOM_PERMS}-match cap, so its permutation set would be '
              'incomplete; dropping fragment.', flush=True)
        return None

    cg_atom_perms = []
    for perm_inds in matches_to_map_to_sub:
        if len(perm_inds) != substruct.GetNumAtoms():
            continue
        try:
            cg_atom_perms.append((Chem.RenumberAtoms(Chem.Mol(substruct), perm_inds), perm_inds))
        except Exception:
            continue
    if len(cg_atom_perms) == 0:
        print(f'[WARNING] manually_remove_Hs: {smiles_no_Hs!r} produced no '
              'whole-molecule atom mapping back onto the fragment; dropping it.',
              flush=True)
        return None
    if return_single_mol_or_perms == 'single':
        return cg_atom_perms[0], annotated_key
    return cg_atom_perms, annotated_key

def fragment_on_bond_d(mol, radius, min_frag_size=None, max_frag_size=None, annotations=None):
    if annotations is None:
        annotations = parent_atom_annotations(mol)
    submols = []
    for atom in mol.GetAtoms():
        amap = {}
        submol = Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(
            mol, radius, atom.GetIdx(), enforceSize=False), atomMap=amap)
        frag_size = submol.GetNumHeavyAtoms()
        if ((min_frag_size is not None and frag_size < min_frag_size) or
                (max_frag_size is not None and frag_size > max_frag_size)):
            continue
        for orig_idx, sub_idx in amap.items():
            submol.GetAtomWithIdx(sub_idx).SetProp(ANNOTATION_PROP, annotations[orig_idx])
        submols.append((submol, sorted(amap.keys())))
    return submols

def is_organic(mol):
    return any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
