import os
import itertools
import json
import re
import io
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import Counter, defaultdict
import prody as pr

KEY_SCHEMA = 'annot-1'

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
        sub.GetAtomWithIdx(new_index[orig_idx]).SetProp(ANNOTATION_PROP,
                                                        annotations[orig_idx])
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

def enumerate_induced_fragments(mol, min_frag_size=4, max_frag_size=5,
                                annotations=None):
    if annotations is None:
        annotations = parent_atom_annotations(mol)
    fragments = {}
    for atom_indices in connected_atom_subsets(mol, min_frag_size, max_frag_size):
        sub = induced_submol(mol, atom_indices, annotations)
        heavy = [a for a in sub.GetAtoms() if a.GetAtomicNum() != 1]
        if not heavy or all(a.GetAtomicNum() == 6 for a in heavy):
            continue
        result = manually_remove_Hs(sub, 'single')
        if result is None:
            continue
        fragments.setdefault(result[1], []).append(atom_indices)
    return fragments

def get_fragments(bond_radius, mol, min_frag_size=4, max_frag_size=5, quiet=True):
    filtered_frags = {}
    substructs = []
    annotations = parent_atom_annotations(mol)
    for rad in list(range(1, bond_radius + 1)):
        substructs += fragment_on_bond_d(mol, rad, min_frag_size, max_frag_size,
                                         annotations)
    for orig_sub, orig_mol_inds in substructs:
        frag_size = orig_sub.GetNumHeavyAtoms()
        if frag_size < min_frag_size or frag_size > max_frag_size:
            continue
        _results = manually_remove_Hs(orig_sub, return_single_mol_or_perms='perms')
        if _results is None:
            continue

        sub_perms_mols_inds, sub_smiles = _results
        for sub, perm_inds in sub_perms_mols_inds: 
            heavy_atoms = [a for a in sub.GetAtoms() if a.GetAtomicNum() != 1]
            if not heavy_atoms or all(a.GetAtomicNum() == 6 for a in heavy_atoms):
                continue
            substruct_data = (sub, perm_inds, orig_mol_inds)
            if sub_smiles not in filtered_frags:
                filtered_frags[sub_smiles] = [substruct_data]
            else:
                already_present = False
                for existing_sub, existing_sub_inds, existing_orig_mol_inds in filtered_frags[sub_smiles]:
                    if (perm_inds, orig_mol_inds) == (existing_sub_inds, existing_orig_mol_inds):
                        already_present = True
                        break
                if not already_present:
                    filtered_frags[sub_smiles].append(substruct_data)

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
    for i in range(n):
        for j in range(i+1, n):
            inter = len(sets[i] & sets[j])
            denom = max(len(sets[i]), len(sets[j]))
            if inter >= threshold * denom:
                union(i, j)
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(data[i])

    return list(groups.values())

ANNOTATION_PROP = '_vdgAtomAnnotation'

def atom_annotations(atom):
    if atom.GetAtomicNum() == 1:
        return ''
    parts = []
    ring = ring_query_for_atom(atom)
    if ring is not None:
        parts.append(ring)
    if atom.GetAtomicNum() == 6:
        parts.append('H0' if atom.GetTotalNumHs(includeNeighbors=True) == 0
                     else '!H0')
    else:
        heavy_degree = sum(1 for nbr in atom.GetNeighbors()
                           if nbr.GetAtomicNum() != 1)
        parts.append(f'D{heavy_degree}')
    return ';'.join(parts)

def ring_query_for_atom(atom):
    if atom.GetIsAromatic():
        return None
    if not atom.IsInRing():
        return '!R'
    ring_info = atom.GetOwningMol().GetRingInfo()
    return f'r{min(ring_info.AtomRingSizes(atom.GetIdx()))}'

_BRACKET_ATOM = re.compile(r'\[([^\]]+)\]')
_BRACKET_ATOM_LOOSE = re.compile(r'\[[^\]]*\]')
_EXPLICIT_H = re.compile(r'H(?![a-z])\d*')
_NON_LETTER = re.compile(r'[^A-Za-z]')

def _remove_bracket_hydrogens_in_smiles(smiles: str, annotations=None) -> str:
    counter = itertools.count()

    def _rewrite(m):
        inner = m.group(1)
        idx = next(counter)

        if inner == 'H':
            return ''

        cleaned = _EXPLICIT_H.sub('', inner)

        annotation = None if annotations is None else annotations[idx]
        if annotation:
            return f'[{cleaned};{annotation}]'

        if _NON_LETTER.search(cleaned):
            return f'[{cleaned}]'

        if len(cleaned) > 1:
            return f'[{cleaned}]'
        return cleaned

    s = _BRACKET_ATOM.sub(_rewrite, smiles)
    return s

class _AnnotationUnmappable(Exception):
    pass

def _annotations_in_written_order(mol, smiles_w_H):
    if not mol.GetNumAtoms() or not mol.GetAtomWithIdx(0).HasProp(ANNOTATION_PROP):
        return None
    props = mol.GetPropsAsDict(includePrivate=True, includeComputed=True)
    order = props.get('_smilesAtomOutputOrder')
    if order is None:
        raise _AnnotationUnmappable('MolToSmiles did not record '
                                        '_smilesAtomOutputOrder')
    if isinstance(order, str):
        order = json.loads(order)
    order = list(order)
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
        smiles_w_H = Chem.MolToSmiles(substruct, allHsExplicit=True, 
                                  isomericSmiles=False)
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
    max_matches = 1 if return_single_mol_or_perms == 'single' else _MAX_CG_ATOM_PERMS
    try:
        matches_to_map_to_sub = substruct.GetSubstructMatches(mol_from_smarts,
            uniquify=False,
            maxMatches=max_matches)
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
    n_atoms = substruct.GetNumAtoms()
    cg_atom_perms = []
    for perm_inds in list(matches_to_map_to_sub):
        if len(perm_inds) != n_atoms:
            continue
        mol_copy = Chem.Mol(substruct)
        try:
            renumbered_substruct = Chem.RenumberAtoms(mol_copy, perm_inds)
        except Exception:
            continue
        cg_atom_perms.append((renumbered_substruct, perm_inds))
    if len(cg_atom_perms) == 0:
        print(f'[WARNING] manually_remove_Hs: {smiles_no_Hs!r} produced no '
              'whole-molecule atom mapping back onto the fragment; dropping it.',
              flush=True)
        return None
    if return_single_mol_or_perms == 'single':
        return cg_atom_perms[0], annotated_key
    return cg_atom_perms, annotated_key
def parent_atom_annotations(mol):
    return [atom_annotations(atom) for atom in mol.GetAtoms()]

def fragment_on_bond_d(mol, radius, min_frag_size=None, max_frag_size=None,
                       annotations=None):
    atoms = mol.GetAtoms()
    if annotations is None:
        annotations = parent_atom_annotations(mol)
    submols = []
    for atom in atoms:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom.GetIdx(), 
            enforceSize=False)
        amap = {}
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
        if min_frag_size is not None or max_frag_size is not None:
            frag_size = submol.GetNumHeavyAtoms()
            if ((min_frag_size is not None and frag_size < min_frag_size) or
                    (max_frag_size is not None and frag_size > max_frag_size)):
                continue
        for orig_idx, sub_idx in amap.items():
            submol.GetAtomWithIdx(sub_idx).SetProp(ANNOTATION_PROP,
                                                   annotations[orig_idx])
        orig_inds = sorted(amap.keys())
        submols.append((submol, orig_inds))
    return submols

def is_organic(mol):
    return any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())

_LIGAND_EXCLUDE_CLAUSE = 'ion or water or resname SEP or resname TPO or resname MSE'

def _heavy_element_counts(sel):
    heavy = sel.select("not element H D") if sel is not None else None
    return Counter(el.upper() for el in heavy.getElements()) if heavy is not None else Counter()

def _template_element_counts(template, n=1):
    counts = Counter(atom.GetSymbol().upper() for atom in template.GetAtoms())
    return counts if n == 1 else Counter({el: c * n for el, c in counts.items()})

def _disambiguate_ligand_pool(pool, query_lig_smiles):
    resnames = sorted(set(pool.getResnames()))
    if len(resnames) == 1:
        return pool, resnames
    template = Chem.MolFromSmiles(query_lig_smiles)
    if template is None:
        return pool, resnames
    matches = []
    for name in resnames:
        sel = pool.select(f"resname {name}")
        if sel is None:
            continue
        n_residues = len(set(sel.getResindices().tolist()))
        if n_residues > 0 and _heavy_element_counts(sel) \
                == _template_element_counts(template, n_residues):
            matches.append(name)
    if len(matches) == 1:
        return pool.select(f"resname {matches[0]}"), matches
    return pool, resnames

def identify_ligand_selection(query_struct, query_lig_smiles):
    hetatms = query_struct.select(f'hetatm and not ({_LIGAND_EXCLUDE_CLAUSE})')
    if hetatms is not None:
        sel, resnames = _disambiguate_ligand_pool(hetatms, query_lig_smiles)
        if len(resnames) == 1:
            return sel, resnames[0]

    all_pool = query_struct.select(
        f'not (protein or nucleic or {_LIGAND_EXCLUDE_CLAUSE})')
    if all_pool is None:
        raise ValueError("get_query_ligand_mol: could not find ligand atoms.")
    sel, resnames = _disambiguate_ligand_pool(all_pool, query_lig_smiles)
    if len(resnames) != 1:
        raise ValueError(
            f"get_query_ligand_mol: expected exactly one ligand "
            f"residue, got: {resnames}")
    return sel, resnames[0]

def classify_ligand_instances(query_struct, query_lig_smiles):
    sel, resname = identify_ligand_selection(query_struct, query_lig_smiles)
    template = Chem.MolFromSmiles(query_lig_smiles)
    if template is None:
        return resname, []
    resindices = sorted(set(sel.getResindices().tolist()))
    instances = []
    for resindex in resindices:
        res = sel.select(f"resindex {resindex}")
        if res is not None and _heavy_element_counts(res) == _template_element_counts(template):
            instances.append((resindex,
                f"{res.getChids()[0]}{int(res.getResnums()[0])}{res.getIcodes()[0]}".strip()))
    if len(instances) != len(resindices):
        return resname, []
    return resname, instances

def select_struct_for_ligand_instance(query_struct, resname, other_resindices):
    if not other_resindices:
        return query_struct
    return query_struct.select(
        f"not (resname {resname} and resindex {' '.join(str(i) for i in other_resindices)})")

def get_query_ligand_mol(query_struct_or_path, query_lig_smiles):
    if isinstance(query_struct_or_path, str):
        if query_struct_or_path.endswith('.pdb') or query_struct_or_path.endswith('.pdb.gz'):
            query_struct = pr.parsePDB(query_struct_or_path)
        elif query_struct_or_path.endswith('.cif') or query_struct_or_path.endswith('.cif.gz'):
            query_struct = pr.parseCIF(query_struct_or_path)
        else:
            raise ValueError(f"Unsupported file format: {query_struct_or_path}")
    else:
        query_struct = query_struct_or_path

    query_hetatms, ligname = identify_ligand_selection(query_struct, query_lig_smiles)

    n_query_ligand_residues = len(set(query_hetatms.getResindices().tolist()))
    if n_query_ligand_residues != 1:
        raise ValueError(
            f"get_query_ligand_mol: expected exactly one ligand "
            f"residue instance of {ligname!r}, got {n_query_ligand_residues} "
            "residues.")

    buf = io.StringIO()
    pr.writePDBStream(buf, query_hetatms)
    pdb_block = buf.getvalue()
    query_pdb_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True)
    query_lig_template = Chem.MolFromSmiles(query_lig_smiles)
    if query_pdb_mol is None or query_lig_template is None:
        print(f"[ERROR] Processing {query_struct_or_path}: RDKit failed to parse molecule",
              flush=True)
        return None

    if query_pdb_mol.GetNumAtoms() != query_lig_template.GetNumAtoms():
        print(f"[WARNING] Processing {query_struct_or_path}: ligand selection has "
              f"{query_pdb_mol.GetNumAtoms()} heavy atoms but the query SMILES has "
              f"{query_lig_template.GetNumAtoms()}; bond orders are assigned only to the "
              "atoms the template matches, and the rest stay single-bonded.",
              flush=True)

    try:
        query_pdb_mol_assigned_bonds = AllChem.AssignBondOrdersFromTemplate(
            query_lig_template, query_pdb_mol)
        stripped = manually_remove_Hs(query_pdb_mol_assigned_bonds,
                                      return_single_mol_or_perms='single')
        if stripped is None:
            print(f"[ERROR] Processing {query_struct_or_path}: could not strip hydrogens "
                  "from the ligand (see the warning above)", flush=True)
            return None
        (query_pdb_mol_assigned_bonds_no_H, _pdb_mol_perm_inds), _smiles_no_H = stripped
        Chem.GetSymmSSSR(query_pdb_mol_assigned_bonds_no_H)
        return query_pdb_mol_assigned_bonds_no_H
    except Exception as e:
        print(
            f"[ERROR] Processing {query_struct_or_path}: "
            f"ligand hydrogen removal failed: {type(e).__name__}: {e}",
            flush=True,
        )
        return None

def submol_from_match(mol, match):
    conf_src = mol.GetConformer()
    em = Chem.RWMol()
    conf = Chem.Conformer(len(match))
    for new_idx, orig_idx in enumerate(match):
        em.AddAtom(Chem.Atom(mol.GetAtomWithIdx(orig_idx)))
        conf.SetAtomPosition(new_idx, conf_src.GetAtomPosition(orig_idx))
    positions = {orig_idx: new_idx for new_idx, orig_idx in enumerate(match)}
    for i, orig_i in enumerate(match):
        for bond in mol.GetAtomWithIdx(orig_i).GetBonds():
            other = bond.GetOtherAtomIdx(orig_i)
            j = positions.get(other)
            if j is not None and j > i:
                em.AddBond(i, j, bond.GetBondType())
    sub = em.GetMol()
    sub.AddConformer(conf, assignId=True)
    return sub

def check_vdg_job_status(sub_smiles, vdg_lib_dir):
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
        "In vdg lib but incomplete": [],
        "In vdg lib and in include list": []}

    for sub_smiles in frags_in_lib:
        if sub_smiles in frags_to_exclude:
            groups["Excluded"].append(sub_smiles)
        elif (frags_to_include != 'all' and isinstance(frags_to_include, list) and 
              sub_smiles not in frags_to_include):
            groups["Not in include list"].append(sub_smiles)
        elif not frags_in_lib[sub_smiles]:
            groups["In vdg lib but incomplete"].append(sub_smiles)
        else:
            groups["In vdg lib and in include list"].append(sub_smiles)

    if groups:
        print("\n--- Fragment vdG Library Summary ---", file=logfile_fh)
        for category, frags in groups.items():
            if frags:
                print(f"{category}:", file=logfile_fh)
                print("   ", ", ".join(frags), file=logfile_fh)
        print("\n", file=logfile_fh)
