import os
import re
import hashlib
import argparse
import pickle as pkl
from collections import defaultdict
from itertools import groupby

from ligand_vdgs.functions import Frags

from ligand_vdgs.functions import utils as utils_module
from ligand_vdgs.functions.utils import (
    fragment_key_query_mol, fragment_query_mols_equivalent, mol_from_fragment)
from ligand_vdgs.functions.vdg_npz_utils import FRAGMENT_ALIASES_FILENAME

_SOLVENT_ARTIFACT_HALOGENS = frozenset({17, 35, 53})

def is_solvent_artifact(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _SOLVENT_ARTIFACT_HALOGENS:
            continue
        nbrs = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
        if nbrs.count(8) >= 2 and 6 not in nbrs:
            return True
    return False

def load_frags_dict(path):
    with open(path, 'rb') as handle:
        payload = pkl.load(handle)
    schema = payload.get('key_schema')
    if schema != Frags.KEY_SCHEMA:
        raise ValueError(
            f'{path} was written under key_schema {schema!r}, but this code emits '
            f'{Frags.KEY_SCHEMA!r}. Regenerate it with fragment_database_ligs.py.')
    return payload['frags'], payload['support_pooled'], {
        k: v for k, v in payload.items() if k not in ('frags', 'support', 'support_pooled')}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract fragment SMILES as a one-column work list.")
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl',
                        help="Path to database_frags_dict.pkl. "
                             "Default: resources/database_frags_dict.pkl.")
    parser.add_argument('--output', default='resources/fragment_work_list.txt',
                        help="Path to write the one-column work list. "
                             "Default: resources/fragment_work_list.txt.")
    parser.add_argument('--min-support', default=100, type=int,
                        help="Minimum support -- distinct parent biounits containing "
                             "the fragment with every atom observed -- for a fragment "
                             "to qualify. Read from the dict, which counted it over "
                             "the whole parent database. The threshold is monotone "
                             "and cheap to lower afterwards via "
                             "make_sge_scripts_for_frags --include-only, so start "
                             "high. Default: 100.")
    parser.add_argument('--max-size', default=5, type=int,
                        help="Max fragment heavy-atom count. Default: 5.")
    parser.add_argument('--vdg-lib-dir', default=None,
                        help="Library dir to write fragment_aliases.tsv into, so "
                             "load_fragment_aliases finds it. Omit only for scratch runs.")

    return parser.parse_args()

_BARE_ATOMS = frozenset({'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I',
                         'b', 'c', 'n', 'o', 'p', 's'})
_BRACKET_ATOM = re.compile(r'\[([^\]]*)\]')
_CHARGE = re.compile(r'(?:\+\d+|-\d+|\++|-+)$')

def charge_normalized_fragment(fragment):
    def _strip(match):
        head, sep, tail = match.group(1).partition(';')
        stripped_head = _CHARGE.sub('', head)
        if stripped_head == head:
            return match.group(0)
        stripped = stripped_head + sep + tail
        return stripped if stripped in _BARE_ATOMS else f'[{stripped}]'

    normalized = _BRACKET_ATOM.sub(_strip, fragment)
    return normalized if normalized != fragment else None

def group_protonation_variants(qualifying):
    query_mols = {smiles: fragment_key_query_mol(smiles) for smiles in qualifying}
    _by_size = {n: [s for _, s in grp] for n, grp in groupby(
        sorted((mol.GetNumAtoms(), smiles) for smiles, mol in query_mols.items()
               if mol is not None), key=lambda t: t[0])}

    representatives, aliases = {}, {}
    charge_groups = []
    for smiles in sorted(qualifying):
        normalized = charge_normalized_fragment(smiles)
        mol = fragment_key_query_mol(normalized) if normalized is not None else None
        if mol is None:
            representatives[smiles] = qualifying[smiles]
            continue
        for group in charge_groups:
            if (group[1].GetNumAtoms() == mol.GetNumAtoms()
                    and fragment_query_mols_equivalent(group[1], mol)):
                group[2].append(smiles)
                break
        else:
            charge_groups.append([normalized, mol, [smiles]])

    for normalized, norm_mol, members in charge_groups:
        # norm_mol is always resolvable (it seeded this charge_group), so no None guard.
        twin = normalized if normalized in qualifying else next(
            (cand for cand in _by_size.get(norm_mol.GetNumAtoms(), ())
             if fragment_query_mols_equivalent(norm_mol, query_mols[cand])), None)
        if twin is not None:
            rep, pooled = twin, [twin] + members
        elif len(members) > 1:
            rep, pooled = normalized, members
        else:
            representatives[members[0]] = qualifying[members[0]]
            continue
        representatives[rep] = set().union(*(qualifying[m] for m in pooled))
        aliases |= {member: rep for member in members}
    return representatives, aliases

def prepare_fragments(frags_dict, max_size):
    qualifying = defaultdict(set)
    unparseable = []
    oversized = solvent = 0
    for smiles_dict in frags_dict.values():
        for smiles, lignames in smiles_dict.items():
            mol = mol_from_fragment(smiles)
            if mol is None:
                print(f'[WARNING] skipping unparseable fragment: {smiles}')
                unparseable.append(smiles)
                continue
            if mol.GetNumHeavyAtoms() > max_size:
                oversized += 1
                continue
            if is_solvent_artifact(mol):
                solvent += 1
                continue
            qualifying[smiles].update(lignames)

    for count, desc in ((len(unparseable), 'failed to parse'),
                        (oversized, f'exceeded --max-size ({max_size})'),
                        (solvent, 'were solvent artifacts')):
        if count:
            print(f'[WARNING] {count} fragment(s) {desc} and were excluded.')

    return group_protonation_variants(dict(qualifying))

def _prepared_cache_key(frags_dict_path, max_size):
    "sha256 over the dict file and the sources that define the grouping, plus max_size."
    digest = hashlib.sha256(f'max_size={max_size}\n'.encode())
    for path in (frags_dict_path, __file__, utils_module.__file__):
        digest.update(f'\x00{os.path.basename(path)}\x00'.encode())
        with open(path, 'rb') as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b''):
                digest.update(chunk)
    return digest.hexdigest()[:24]

def prepared_cache_path(frags_dict_path, max_size):
    return os.path.join(
        os.path.expanduser(os.environ.get('VDG_SCRATCH', '~/docking/scratch')),
        'prepared_fragments', _prepared_cache_key(frags_dict_path, max_size) + '.pkl')

def prepare_fragments_cached(frags_dict, max_size, frags_dict_path, quiet=False):
    """Memoized prepare_fragments under $VDG_SCRATCH, keyed on the sha256 of
    frags_dict_path and the grouping sources (not mtime), plus max_size. Saves
    minutes over the production dict on repeated dry runs/top-ups. No staleness
    fallback -- a miss recomputes."""
    path = prepared_cache_path(frags_dict_path, max_size)
    if os.path.exists(path):
        with open(path, 'rb') as handle:
            prepared = pkl.load(handle)
        if not quiet:
            print(f'Reusing prepared fragments from {path}')
        return prepared
    prepared = prepare_fragments(frags_dict, max_size)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f'{path}.{os.getpid()}.tmp'
    with open(tmp, 'wb') as handle:
        pkl.dump(prepared, handle, protocol=pkl.HIGHEST_PROTOCOL)
    os.replace(tmp, path)
    if not quiet:
        print(f'Cached prepared fragments at {path}')
    return prepared

def _unresolved_fragment_error(unresolved, max_size):
    too_big, solventy, absent = [], [], []
    for frag in unresolved:
        mol = mol_from_fragment(frag)
        if mol is None:
            absent.append(frag)
        elif mol.GetNumHeavyAtoms() > max_size:
            too_big.append((frag, mol.GetNumHeavyAtoms()))
        elif is_solvent_artifact(mol):
            solventy.append(frag)
        else:
            absent.append(frag)
    problems = []
    if too_big:
        problems.append(
            "excluded by --max-size {}: {}".format(
                max_size, ", ".join(f"{f} ({n} heavy atoms)" for f, n in too_big)))
    if solventy:
        problems.append(f"excluded as crystallographic solvent artifacts: {solventy}")
    if absent:
        problems.append(
            f"absent from the fragment dict even up to equivalent atom "
            f"ordering, so they have no vdG sites to mine -- check the "
            f"annotated-SMARTS spelling, ring primitives included, with "
            f"scripts/lookup_fragment_key.py: {absent}")
    return (f"{len(unresolved)} requested fragment(s) could not be selected. "
            + "; ".join(problems))

def select_fragments(frags_dict, min_support, max_size, support=None,
                     return_aliases=False, prepared=None, include=None,
                     include_only=False):
    prepared = (prepared if prepared is not None
                else prepare_fragments(frags_dict, max_size))
    representatives, aliases = prepared

    included, unresolved = resolve_include_fragments(include, prepared)
    if unresolved:
        raise ValueError(_unresolved_fragment_error(unresolved, max_size))

    if include_only:
        selected = set(included.values())
    elif min_support <= 0:
        selected = set(representatives)
    else:
        if support is None:
            raise ValueError(
                f"select_fragments needs `support` (from load_frags_dict) to apply "
                f"a threshold of {min_support}.")
        missing = [s for s in representatives if s not in support]
        if missing:
            raise ValueError(
                f"support is missing {len(missing)} fragment(s), e.g. {missing[:3]}. "
                "It must cover every fragment surviving the size/solvent filters -- "
                "check for a --max-size mismatch against the fragment dict.")
        selected = {s for s in representatives if support[s] >= min_support}

    selected = sorted(selected | set(included.values()))
    _kept = set(selected)
    aliases = {alias: rep for alias, rep in aliases.items() if rep in _kept}
    return (selected, aliases) if return_aliases else selected

def resolve_include_fragments(include, prepared):
    representatives, aliases = prepared
    resolved, unresolved = {}, []
    if not include:
        return resolved, unresolved

    sized = sorted((m.GetNumAtoms(), key) for key in list(representatives) + list(aliases)
                   if (m := fragment_key_query_mol(key)) is not None)
    _by_size = {n: [k for _, k in grp] for n, grp in groupby(sized, key=lambda t: t[0])}

    for requested in include:
        if requested in representatives:
            resolved[requested] = requested
            continue
        if requested in aliases:
            resolved[requested] = aliases[requested]
            continue
        mol = fragment_key_query_mol(requested)
        match = None if mol is None else next(
            (cand for cand in _by_size.get(mol.GetNumAtoms(), ())
             if fragment_query_mols_equivalent(mol, fragment_key_query_mol(cand))), None)
        if match is None:
            unresolved.append(requested)
        else:
            resolved[requested] = aliases.get(match, match)
    return resolved, unresolved

def fragment_dict_keys(frags_dict):
    return {smiles for smiles_dict in frags_dict.values() for smiles in smiles_dict}

def alias_kind(representative, dict_keys):
    return 'neutral_twin' if representative in dict_keys else 'promoted'

def write_fragment_aliases(path, aliases, dict_keys):
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, 'w') as f:
        f.write('# alias\trepresentative\tkind\n')
        f.write('# kind=promoted: the representative is the charge-stripped SMARTS of '
                'its aliases and is NOT a key of the fragment dict (no CCD ligand is '
                'drawn that way); use scripts/lookup_fragment_key.py to find the '
                'dict keys it covers.\n')
        for alias in sorted(aliases):
            rep = aliases[alias]
            f.write(f'{alias}\t{rep}\t{alias_kind(rep, dict_keys)}\n')

def resolve_fragment_key(key, candidates):
    query = fragment_key_query_mol(key)
    loose = fragment_key_query_mol(charge_normalized_fragment(key) or key)
    found = []
    for cand in candidates:
        if cand == key:
            found.append((cand, 'exact'))
            continue
        cand_mol = fragment_key_query_mol(cand)
        if cand_mol is None:
            continue
        if query is not None and fragment_query_mols_equivalent(query, cand_mol):
            found.append((cand, 'equivalent'))
            continue
        if loose is None:
            continue
        cand_loose_key = charge_normalized_fragment(cand)
        if fragment_query_mols_equivalent(
                loose, fragment_key_query_mol(cand_loose_key)
                if cand_loose_key is not None else cand_mol):
            found.append((cand, 'charge_variant'))
    order = {'exact': 0, 'equivalent': 1, 'charge_variant': 2}
    return sorted(found, key=lambda pair: (order[pair[1]], pair[0]))

def main():
    args = parse_args()

    if not os.path.isfile(args.frags_dict):
        raise FileNotFoundError(f"--frags-dict path does not exist: {args.frags_dict}")

    frags_dict, support, meta = load_frags_dict(args.frags_dict)
    print(f'{args.frags_dict}: key_schema {meta.get("key_schema")}, support counted '
          f'in {meta.get("params", {}).get("support_unit", "?")} over '
          f'{meta.get("db_identity", {}).get("sha256", "?")[:16]}...')

    smiles_to_run, aliases = select_fragments(
        frags_dict, args.min_support, args.max_size, support=support,
        return_aliases=True,
        prepared=prepare_fragments_cached(frags_dict, args.max_size, args.frags_dict))

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, 'w') as f:
        for smiles in smiles_to_run:
            f.write(f"{smiles}\n")

    if args.vdg_lib_dir:
        alias_path = os.path.join(args.vdg_lib_dir, FRAGMENT_ALIASES_FILENAME)
    else:
        alias_path = os.path.splitext(args.output)[0] + '_aliases.tsv'
        print(f'[WARNING] --vdg-lib-dir not given; writing aliases to {alias_path}, '
              'which load_fragment_aliases will not find.')
    dict_keys = fragment_dict_keys(frags_dict)
    write_fragment_aliases(alias_path, aliases, dict_keys)
    if aliases:
        reps = set(aliases.values())
        promoted = sorted(r for r in reps if alias_kind(r, dict_keys) == 'promoted')
        print(f'Collapsed {len(aliases)} protonation variant(s) into '
              f'{len(reps)} representative(s); wrote {alias_path}.')
        if promoted:
            print(f'{len(promoted)} representative(s) are promoted charge-stripped keys '
                  f'absent from the fragment dict (see {alias_path}): {promoted}')

    print(f'Wrote {len(smiles_to_run)} fragments to {args.output}.')

if __name__ == '__main__':
    main()
