'''
Builds the fragment vocabulary: every connected subgraph of every database ligand,
scored by biounit support.

Support definition:
    support(f) = number of biounits with a fully-observed match of f in some ligand
    instance. Set semantics over biounits: NCS copies don't multiply support, and
    a biounit holding both ATP and ADP contributes 1 (union), not 2, to an adenine
    fragment's support.
'''

import argparse
import concurrent.futures
import multiprocessing as mp
import os
import pickle as pkl
import re
import time
from collections import defaultdict
from concurrent.futures.process import BrokenProcessPool

import numpy as np
from rdkit import Chem, RDLogger

from ligand_vdgs.functions import Frags, ccd_templates, ligand_perception
from ligand_vdgs.functions.utils import fragment_key_query_mol, fragment_query_mols_equivalent

RDLogger.DisableLog('rdApp.*')

def parse_args():
    parser = argparse.ArgumentParser(
        description='Enumerate all connected subgraphs of the database ligands and '
                    'count each one\'s support in distinct parent biounits.')
    parser.add_argument('--roster', default='resources/ligand_roster.pkl',
        help='Ligand roster from build_ligand_roster.py.')
    parser.add_argument('--outdir', default='resources', help='Where to write database_frags_dict.pkl.')
    parser.add_argument('--logfile', default='logs/fragment_database_ligs.log')
    parser.add_argument('--min-frag-size', default=4, type=int, help='Minimum heavy atoms per fragment.')
    parser.add_argument('--max-frag-size', default=5, type=int,
        help='Maximum heavy atoms per fragment. Settled at 5: a 6-ring is only seen via its '
             '5-atom heteroatom-bearing arcs, by design -- the whole ring is too large/diffuse '
             'a signal. 5-rings still close fully.')
    parser.add_argument('--record-min-support', default=1, type=int,
        help='File-size floor on the output, NOT the selection threshold (that is --min-support '
             'in extract_fragment_smiles.py, applied at read time).')
    parser.add_argument('--num-procs', default=1, type=int,
        help='Worker processes; result is independent of this value.')
    parser.add_argument('--force', action='store_true',
        help='Overwrite an existing output file (default: error, so a qsub cannot silently '
             'clobber an in-progress build).')
    return parser.parse_args()

_FRAG_PARAMS = None

def _init_worker(min_frag_size, max_frag_size):
    global _FRAG_PARAMS
    _FRAG_PARAMS = (min_frag_size, max_frag_size)
    RDLogger.DisableLog('rdApp.*')

def enumerate_one_resname(resname):
    min_frag_size, max_frag_size = _FRAG_PARAMS
    out = {'resname': resname, 'status': None, 'frags': [], 'warnings': []}

    perceived = ligand_perception.perceive_ligand_graph(resname)
    mol = perceived.mol
    if mol is None:
        out['status'] = 'no_template'
        out['warnings'].append(f'[WARNING] {resname}: in the roster but has no usable CCD template now; the roster and the template store disagree.')
        return out

    names = _template_heavy_atom_names(resname)
    if names is None or len(names) != mol.GetNumAtoms():
        out['status'] = 'name_map_failed'
        out['warnings'].append(f'[WARNING] {resname}: could not map template atom names onto its Mol ({"None" if names is None else len(names)} names for {mol.GetNumAtoms()} atoms); skipping.')
        return out

    if undesired_elements(mol):
        out['status'] = 'not_druglike'
        return out

    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        out['status'] = 'sanitize_failed'
        out['warnings'].append(f'[WARNING] {resname}: sanitization failed ({e}).')
        return out

    if not Frags.is_organic(mol):
        out['status'] = 'not_druglike'
        return out

    for idx, name in enumerate(names):
        mol.GetAtomWithIdx(idx).SetProp(_NAME_PROP, name)

    fragments = Frags.enumerate_induced_fragments(mol, min_frag_size, max_frag_size)
    for key, atom_index_sets in fragments.items():
        name_sets = []
        elements = None
        for atom_indices in atom_index_sets:
            name_sets.append(tuple(sorted(mol.GetAtomWithIdx(i).GetProp(_NAME_PROP) for i in atom_indices)))
            if elements is None:
                elements = fragment_elements(mol, atom_indices)
        out['frags'].append((elements, key, name_sets))
    out['status'] = 'ok' if fragments else 'no_frags'
    return out

_NAME_PROP = '_vdgPdbAtomName'

def _template_heavy_atom_names(resname):
    from ligand_vdgs.functions import ccd_templates
    try:
        template = ccd_templates.get_template(resname)
    except Exception:
        return None
    if template is None:
        return None
    return [atom.name for atom in template.atoms if atom.element not in ('H', 'D')]

def support_from_roster(roster, keys_by_resname):
    types = roster['types']
    biounits = roster['biounits']

    key_ids = {}
    keyset_ids = {}
    keysets = []
    type_keyset = np.empty(len(types), dtype=np.int64)

    stats = defaultdict(int)
    for type_idx, (resname, observed) in enumerate(types):
        observed_set = set(observed)
        present = []
        for key, name_sets in keys_by_resname.get(resname, {}).items():
            if any(observed_set.issuperset(name_set) for name_set in name_sets):
                key_id = key_ids.get(key)
                if key_id is None:
                    key_id = len(key_ids)
                    key_ids[key] = key_id
                present.append(key_id)
        arr = np.array(sorted(present), dtype=np.int32)
        token = arr.tobytes()
        keyset_id = keyset_ids.get(token)
        if keyset_id is None:
            keyset_id = len(keysets)
            keyset_ids[token] = keyset_id
            keysets.append(arr)
        type_keyset[type_idx] = keyset_id
    stats['num_roster_types'] = len(types)
    stats['num_distinct_keysets'] = len(keysets)

    counts = np.zeros(len(key_ids), dtype=np.int64)
    for type_ids in biounits.values():
        parts = [keysets[type_keyset[i]] for i in type_ids]
        if not parts:
            continue
        union = np.unique(np.concatenate(parts)) if len(parts) > 1 else parts[0]
        if union.size:
            counts[union] += 1
    stats['num_biounits'] = len(biounits)
    return {key: int(counts[key_id]) for key, key_id in key_ids.items()}, dict(stats)

def fragment_elements(mol, atom_indices):
    return ''.join(sorted(mol.GetAtomWithIdx(i).GetSymbol() for i in atom_indices
                          if mol.GetAtomWithIdx(i).GetSymbol() != 'H'))

def collapse_keys(elements_of_key, comparison_errors, log=print):
    groups = defaultdict(list)
    for key in sorted(elements_of_key):
        groups[(elements_of_key[key], _key_token_signature(key))].append(key)

    alias = {}
    smarts_cache = {}
    n_pairs = 0
    for group in groups.values():
        if len(group) == 1:
            alias[group[0]] = group[0]
            continue
        representatives = []
        for key in group:
            if key not in smarts_cache:
                smarts_cache[key] = fragment_key_query_mol(key)
            query = smarts_cache[key]
            found = None
            for existing in ([] if query is None else representatives):
                existing_query = smarts_cache[existing]
                if existing_query is None:
                    continue
                n_pairs += 1
                try:
                    if fragment_query_mols_equivalent(existing_query, query):
                        found = existing
                        break
                except Exception as e:
                    entry = comparison_errors.setdefault(str(e), [0, (key, existing)])
                    entry[0] += 1
                    continue
            if found is None:
                representatives.append(key)
                alias[key] = key
            else:
                alias[key] = found
    n_collapsed = sum(1 for k, v in alias.items() if k != v)
    log(f'  collapsed {n_collapsed} degenerate spellings out of {len(alias)} raw '
        f'keys ({n_pairs} equivalence checks)')
    return alias, n_collapsed

_BRACKET_ATOM = re.compile(r'\[([^\]]*)\]')

def _key_token_signature(key):
    tokens = sorted(_BRACKET_ATOM.findall(key))
    bare = sorted(ch for ch in _BRACKET_ATOM.sub('', key) if ch.isalpha())
    return (tuple(tokens), tuple(bare))

def sort_frag_dict(frag_dict):
    return {elements: dict(sorted(fragments.items(), key=lambda x: len(set(x[1])), reverse=True))
            for elements, fragments in frag_dict.items()}

UNDESIRED_ELEMENTS = {
    'Be', 'Ba', 'Bi', 'Pt', 'Ru', 'Ir', 'Fe', 'Zn', 'Mg', 'Cu', 'Sn', 'Zr', 'Pb', 'Ga',
    'Hg', 'Pd', 'Si', 'Mo', 'W', 'Se', 'Mn', 'Ti', 'Y', 'V', 'Ni', 'Rh', 'Te', 'Au', 'Ag',
    'Co', 'Sb', 'Re', 'As', 'Cd', 'Hf', 'Na', 'Li', 'Ca', 'Os', 'Cr', 'In', 'Al', 'U', 'K',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn',
}
def undesired_elements(mol):
    return any(a.GetSymbol() in UNDESIRED_ELEMENTS for a in mol.GetAtoms())

def load_roster(path):
    with open(path, 'rb') as f:
        roster = pkl.load(f)
    for field in ('types', 'biounits', 'db_identity', 'ccd_identity'):
        if field not in roster:
            raise SystemExit(f'[ERROR]: {path} is not a ligand roster (no {field!r} key). Build it with build_ligand_roster.py.')
    if roster.get('partial'):
        raise SystemExit(f'[ERROR]: {path} was built with --limit, so it covers only part of the parent database and every support value from it is an undercount. Rebuild it without --limit.')
    current = ccd_templates.store_identity()
    if roster['ccd_identity'] != current:
        raise SystemExit(
            f'[ERROR]: {path} was built from a different CCD template store. Recorded '
            f'{roster["ccd_identity"]}; current {current}. Restore that store or rebuild '
            'the roster before generating fragments.')
    return roster

def enumerate_resnames(resnames, min_frag_size, max_frag_size, num_procs, log=print):
    keys_by_resname = {}
    elements_of_key = {}
    stats = defaultdict(int)
    warnings = []
    started = time.time()

    def tally(out):
        stats['num_ligs_total'] += 1
        stats[f'num_ligs_{out["status"]}'] += 1
        warnings.extend(out['warnings'])
        if not out['frags']:
            return
        per_key = {}
        for elements, key, name_sets in out['frags']:
            per_key[key] = name_sets
            elements_of_key.setdefault(key, elements)
        keys_by_resname[out['resname']] = per_key
        stats['num_ligs_with_frags'] += 1

    _init_worker(min_frag_size, max_frag_size)
    if num_procs > 1:
        ctx = mp.get_context('spawn')
        pool = concurrent.futures.ProcessPoolExecutor(
            max_workers=num_procs, mp_context=ctx, initializer=_init_worker,
            initargs=(min_frag_size, max_frag_size))
        try:
            futures = {pool.submit(enumerate_one_resname, r): r for r in resnames}
            for done, fut in enumerate(concurrent.futures.as_completed(futures), 1):
                try:
                    tally(fut.result())
                except BrokenProcessPool:
                    raise
                except Exception as exc:
                    stats['num_ligs_errored'] += 1
                    warnings.append(f'[WARNING] {futures[fut]}: {exc!r}')
                if done % 2000 == 0:
                    log(f'  enumerated {done}/{len(resnames)} ligands, {len(elements_of_key)} distinct keys, {done / max(time.time() - started, 1e-9):.1f}/s')
        except BrokenProcessPool as exc:
            raise RuntimeError(
                f'An enumeration worker died without returning a result ({exc}). The '
                f'usual cause is the OOM killer; re-run with fewer --num-procs or a '
                f'larger -l mem_free (which is PER SLOT under -pe smp).') from exc
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
    else:
        for done, resname in enumerate(resnames, 1):
            try:
                tally(enumerate_one_resname(resname))
            except Exception as exc:
                stats['num_ligs_errored'] += 1
                warnings.append(f'[WARNING] {resname}: {exc!r}')
            if done % 2000 == 0:
                log(f'  enumerated {done}/{len(resnames)} ligands')
    stats['enumeration_s'] = round(time.time() - started, 1)
    return keys_by_resname, elements_of_key, dict(stats), warnings

def support_frontier(support):
    if not support:
        return []
    values = np.fromiter(support.values(), dtype=np.int64, count=len(support))
    return [(t, int((values >= t).sum())) for t in (1, 2, 5, 10, 25, 50, 100, 250, 500, 1000)]

def apply_aliases(keys_by_resname, elements_of_key, alias):
    aliased = {}
    elements_of_rep = {}
    resnames_of_rep = defaultdict(list)
    for resname, per_key in keys_by_resname.items():
        merged = defaultdict(set)
        for key, name_sets in per_key.items():
            rep = alias.get(key, key)
            merged[rep].update(name_sets)
            elements_of_rep.setdefault(rep, elements_of_key[key])
        aliased[resname] = {rep: sorted(names) for rep, names in merged.items()}
        for rep in merged:
            resnames_of_rep[rep].append(resname)
    return aliased, elements_of_rep, dict(resnames_of_rep)

def pool_protonation_variants(elements_of_rep, log=print):
    from ligand_vdgs.generate_vdgs.extract_fragment_smiles import group_protonation_variants
    started = time.time()
    representatives, aliases = group_protonation_variants(
        {key: set() for key in elements_of_rep})
    pooled_of_key = {key: aliases.get(key, key) for key in elements_of_rep}
    elements_of_pool = {}
    for key, pooled in pooled_of_key.items():
        elements_of_pool.setdefault(pooled, elements_of_rep[key])
    log(f'  pooled {len(aliases)} protonation variant(s) into '
        f'{len(representatives)} representative(s) in {time.time() - started:.1f}s')
    return pooled_of_key, elements_of_pool

def map_keys(keys_by_resname, mapping):
    out = {}
    for resname, per_key in keys_by_resname.items():
        merged = defaultdict(set)
        for key, name_sets in per_key.items():
            merged[mapping.get(key, key)].update(name_sets)
        out[resname] = {target: sorted(names) for target, names in merged.items()}
    return out

def drop_below_floor(keys_by_resname, elements_of_rep, resnames_of_rep, support, record_min_support):
    keep = {key for key, value in support.items() if value >= record_min_support}
    n_dropped = len(elements_of_rep) - len(keep & set(elements_of_rep))
    elements_of_rep = {k: v for k, v in elements_of_rep.items() if k in keep}
    resnames_of_rep = {k: v for k, v in resnames_of_rep.items() if k in keep}
    keys_by_resname = {resname: {k: v for k, v in per_key.items() if k in keep} for resname, per_key in keys_by_resname.items()}
    support = {k: v for k, v in support.items() if k in keep}
    return keys_by_resname, elements_of_rep, resnames_of_rep, support, n_dropped

def build_vocabulary(elements_of_rep, resnames_of_rep):
    frag_dict = defaultdict(dict)
    stats = defaultdict(int)
    for key, resnames in resnames_of_rep.items():
        frag_dict[elements_of_rep[key]][key] = sorted(set(resnames))
        stats['num_frags_recorded'] += 1
    return dict(frag_dict), dict(stats)

def stamp_log(logfile, roster, params, roster_path, n_resnames):
    log_dir = os.path.dirname(logfile)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
    with open(logfile, 'w') as log:
        log.write(f'{"=" * 25} Fragment vocabulary {"=" * 25}\n')
        log.write(f'IN PROGRESS -- started {time.strftime("%Y-%m-%dT%H:%M:%S")}. If this line survives, the run did not finish.\n')
        log.write(f'roster:             {roster_path}\n')
        log.write(f'roster db_identity: {roster.get("db_identity")}\n')
        log.write(f'roster ccd_identity: {roster.get("ccd_identity")}\n')
        log.write(f'ccd codes:          {n_resnames}\n')
        log.write(f'params:             {params}\n')

def report_stats(logfile, stats, warnings, comparison_errors, frontier, roster,
                 params):
    log_dir = os.path.dirname(logfile)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
    with open(logfile, 'w') as log:
        log.write(f'{"=" * 25} Fragment vocabulary {"=" * 25}\n')
        log.write(f'roster db_identity: {roster.get("db_identity")}\n')
        log.write(f'roster pdb_dir:     {roster.get("pdb_dir")}\n')
        log.write(f'roster ccd_identity: {roster.get("ccd_identity")}\n')
        log.write(f'key_schema:         {Frags.KEY_SCHEMA}\n')
        for key in sorted(params):
            log.write(f'{key}: {params[key]}\n')
        for key in sorted(stats):
            log.write(f'{key}: {stats[key]}\n')
        if frontier:
            log.write('\nSupport frontier (min biounits -> protonation-pooled representatives at or above; one representative is one mining job):\n')
            for threshold, count in frontier:
                log.write(f'  >= {threshold:5d}: {count}\n')
        if comparison_errors:
            log.write('\nSubstructure comparison errors (count, example pair):\n')
            for reason, (count, (smi, existing)) in sorted(comparison_errors.items(), key=lambda kv: kv[1][0], reverse=True):
                log.write(f'  {count}x {reason}\n      e.g. {smi} vs {existing}\n')
        if warnings:
            log.write(f'\nWarnings ({len(warnings)}):\n')
            for warning in warnings[:500]:
                log.write(f'  {warning}\n')
            if len(warnings) > 500:
                log.write(f'  ... and {len(warnings) - 500} more\n')

def output_results(payload, out_pkl):
    tmp = out_pkl + '.tmp'
    with open(tmp, 'wb') as f:
        pkl.dump(payload, f, protocol=pkl.HIGHEST_PROTOCOL)
    os.replace(tmp, out_pkl)

def main():
    args = parse_args()
    out_dict_path = os.path.join(args.outdir, 'database_frags_dict.pkl')
    if os.path.exists(out_dict_path) and not args.force:
        raise SystemExit(
            f'[ERROR]: {out_dict_path} already exists. Pass --force to overwrite it, '
            f'or move it aside. Nothing was run, so {args.logfile} (if it exists) '
            f'still describes an earlier run -- do not read it as this one.')
    if args.min_frag_size < 1 or args.max_frag_size < args.min_frag_size:
        raise SystemExit(f'[ERROR]: need 1 <= --min-frag-size <= --max-frag-size, got {args.min_frag_size} and {args.max_frag_size}.')

    ligand_perception.require_template_store()

    roster = load_roster(args.roster)
    resnames = sorted({resname for resname, _names in roster['types']})
    params = {'min_frag_size': args.min_frag_size,
              'max_frag_size': args.max_frag_size,
              'record_min_support': args.record_min_support,
              'enumeration': 'connected-induced-subgraphs',
              'support_unit': 'distinct parent biounits'}
    print(f'{len(resnames)} CCD codes in the roster, {len(roster["biounits"])} '
          f'biounits.', flush=True)
    stamp_log(args.logfile, roster, params, args.roster, len(resnames))

    stats = {}
    warnings = []
    comparison_errors = {}
    frontier = []
    try:
        keys_by_resname, elements_of_key, enum_stats, warnings = enumerate_resnames(
            resnames, args.min_frag_size, args.max_frag_size, args.num_procs)
        stats.update(enum_stats)
        stats['num_raw_keys'] = len(elements_of_key)

        alias, n_collapsed = collapse_keys(elements_of_key, comparison_errors)
        stats['num_keys_collapsed'] = n_collapsed
        keys_by_resname, elements_of_rep, resnames_of_rep = apply_aliases(keys_by_resname, elements_of_key, alias)
        stats['num_keys'] = len(elements_of_rep)

        support, join_stats = support_from_roster(roster, keys_by_resname)
        stats.update(join_stats)

        keys_by_resname, elements_of_rep, resnames_of_rep, support, n_dropped = drop_below_floor(keys_by_resname, elements_of_rep, resnames_of_rep, support, args.record_min_support)
        stats['num_frags_below_record_floor'] = n_dropped
        stats['num_keys_recorded'] = len(elements_of_rep)

        pooled_of_key, elements_of_pool = pool_protonation_variants(elements_of_rep)
        support_pooled, _ = support_from_roster(roster, map_keys(keys_by_resname, pooled_of_key))
        for pooled in elements_of_pool:
            support_pooled.setdefault(pooled, 0)
        stats['num_pooled_keys'] = len(support_pooled)
        frontier = support_frontier(support_pooled)

        frag_dict, vocab_stats = build_vocabulary(elements_of_rep, resnames_of_rep)
        stats.update(vocab_stats)
        stats['num_final_frags'] = sum(len(f) for f in frag_dict.values())
    finally:
        report_stats(args.logfile, stats, warnings, comparison_errors, frontier,
                     roster, params)

    if not stats.get('num_final_frags'):
        raise SystemExit(f'[ERROR]: recorded 0 fragments from {len(resnames)} CCD codes; refusing to write an empty {out_dict_path}. See {args.logfile}.')

    os.makedirs(args.outdir, exist_ok=True)
    output_results({'frags': sort_frag_dict(frag_dict), 'support': support, 'support_pooled': support_pooled,
                    'key_schema': Frags.KEY_SCHEMA, 'db_identity': roster['db_identity'], 'roster_pdb_dir': roster.get('pdb_dir'),
                    'ccd_identity': roster['ccd_identity'],
                    'params': params, 'stats': stats}, out_dict_path)
    print(f'Wrote {out_dict_path}: {stats["num_final_frags"]} fragments. Log: {args.logfile}', flush=True)
    print('Support frontier, in representatives -- one representative is one mining job:', flush=True)
    for threshold, count in frontier:
        print(f'  support >= {threshold:5d}: {count} representatives', flush=True)

if __name__ == '__main__':
    main()
