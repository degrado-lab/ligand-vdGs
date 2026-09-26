"""Rank tetrazole/carboxylate partners and audit charge and nested-key effects.

The primary test is the pooled-backbone single-AA Pearson rank in both directions.
Charge partitions and nested keys are separate diagnostics, not ranking substitutes.
"""
import argparse
import csv
import hashlib
import itertools
import os
from collections import defaultdict

import numpy as np

from ligand_vdgs.functions import ligand_structure, utils
from ligand_vdgs.identify_bioisosteres.common import (
    calc_single_aa_propensities, fragment_keys_nested)
from ligand_vdgs.identify_bioisosteres.compare_aa_profiles import aligned_vectors, _pearson
from ligand_vdgs.identify_bioisosteres.fig1_known_partner_rank import load_profiles

MAIN_TET = '[c;H0]1[n;D2][n;D2][n;D2][n;D2]1'
MAIN_CARB = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D1]'
LONG_CARB = '[C;!R;!H0][C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D1]'
AROM_CARB = '[c;H0][C;!R;H0](=[O;!R;D1])[O;!R;D1]'
CARB_MOTIF = '[C;!R;H0](=[O;!R;D1])[O;!R;D1]'
SIGNS = ('neut', 'neg')

def _profile(key, counts):
    enrichment, _ = calc_single_aa_propensities('', counts=counts, bb_mode='pooled')
    labels = list(enrichment)
    return key, np.array([enrichment[k] for k in labels]), labels, np.array([counts[k] for k in labels]), sum(counts.values())

def _similarity(a, b):
    va, vb, shared = aligned_vectors(a[2], a[1], a[3], b[2], b[1], b[3], min_coverage=0.5)
    return _pearson(va, vb, 5)[0], len(shared)

def _rank(pool, partner):
    matches = [score for label, score in pool if label == partner]
    if len(matches) != 1:
        raise AssertionError(f'[ERROR] partner absent from eligible pool: {partner}')
    return sum(score >= matches[0] for _, score in pool), len(pool)

def _is_tetrazole(key):
    mol = utils.fragment_key_query_mol(key)
    return (mol.GetNumAtoms() == 5 and sorted(a.GetAtomicNum() for a in mol.GetAtoms()) == [6, 7, 7, 7, 7]
            and key.count('1') == 2)

def _sign_counts(vdg_lib_dir, key, sign):
    base = os.path.join(vdg_lib_dir, key, 'nr_vdgs', '1', sign)
    if not os.path.isdir(base): return {}
    counts = {}
    for name in os.listdir(base):
        if not name.endswith('.npz') or 'X' in name[:-4].split('_'): continue
        with np.load(os.path.join(base, name)) as data:
            counts[name[:-4]] = int(data['cluster_num_parents'].sum())
    return counts

def _entry(stem):
    return str(stem).split('_', 1)[0]

def _fold(stem):
    return int(hashlib.md5(_entry(stem).encode()).hexdigest(), 16) % 2

def _fold_counts(vdg_lib_dir, key, fold):
    """Counts from disjoint PDB entries, including nr and member observations."""
    counts, used = defaultdict(int), set()
    for sign in ('pos', 'neut', 'neg', 'unreadable'):
        base = os.path.join(vdg_lib_dir, key, 'nr_vdgs', '1', sign)
        if not os.path.isdir(base): continue
        for name in os.listdir(base):
            if not name.endswith('.npz') or 'X' in name[:-4].split('_'): continue
            with np.load(os.path.join(base, name)) as data:
                by_cluster = defaultdict(set)
                for ids, stems in ((data['cluster_id'], data['nr_parent_biounit']),
                                   (data['mem_cluster_id'], data['mem_parent_biounit'])):
                    for cid, stem in zip(ids, stems):
                        entry = _entry(stem)
                        if _fold(entry) == fold:
                            by_cluster[int(cid)].add(entry)
                            used.add(entry)
                counts[name[:-4]] += sum(len(entries) for entries in by_cluster.values())
    return dict(counts), used

def _holdout_partner_rows(vdg_lib_dir, profiles, min_support):
    """Reciprocal ranks with target profiles built from opposite PDB-entry folds."""
    rows = []
    for partner in (MAIN_CARB, AROM_CARB):
        for tet_fold in (0, 1):
            tet_counts, tet_parents = _fold_counts(vdg_lib_dir, MAIN_TET, tet_fold)
            carb_counts, carb_parents = _fold_counts(vdg_lib_dir, partner, 1 - tet_fold)
            if tet_parents & carb_parents:
                raise AssertionError('[ERROR] parent-entry folds overlap')
            tet, carb = _profile(MAIN_TET, tet_counts), _profile(partner, carb_counts)
            if min(tet[4], carb[4]) < min_support:
                raise ValueError('[ERROR] held-out target falls below min-support')
            similarity, shared = _similarity(tet, carb)
            row = dict(tetrazole=MAIN_TET, carboxylate=partner, tet_fold=tet_fold,
                       support_tet=tet[4], support_carb=carb[4], parent_entries_tet=len(tet_parents),
                       parent_entries_carb=len(carb_parents), overlapping_entries=0,
                       pearson_r=similarity, shared_categories=shared)
            for anchor, other, anchor_profile, other_profile, suffix in (
                    (MAIN_TET, partner, tet, carb, 'tet_to_carb'),
                    (partner, MAIN_TET, carb, tet, 'carb_to_tet')):
                pool = []
                for key, profile in profiles.items():
                    if key == anchor or fragment_keys_nested(anchor, key): continue
                    score, _ = _similarity(anchor_profile, other_profile if key == other else profile)
                    if np.isfinite(score): pool.append((key, score))
                rank, size = _rank(pool, other)
                row[f'rank_{suffix}'] = rank
                row[f'pool_{suffix}'] = size
                row[f'pct_{suffix}'] = rank / size
            row['top_5pct_both'] = row['pct_tet_to_carb'] <= .05 and row['pct_carb_to_tet'] <= .05
            rows.append(row)
    return rows

def _write_tsv(path, rows, fields):
    with open(path, 'w') as file:
        writer = csv.DictWriter(file, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vdg-lib-dir', required=True)
    parser.add_argument('--profiles-dir', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--min-support', type=int, default=50)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    by_key = {p[0]: p for p in load_profiles(args.profiles_dir) if p[4] >= args.min_support}
    lib_keys = [key for key in os.listdir(args.vdg_lib_dir)
                if os.path.isdir(os.path.join(args.vdg_lib_dir, key, 'nr_vdgs'))
                and ligand_structure.check_vdg_job_status(key, args.vdg_lib_dir)]
    tet_keys = sorted(key for key in lib_keys if _is_tetrazole(key))
    carb_keys = sorted(key for key in lib_keys if CARB_MOTIF in key)
    for key in tet_keys + carb_keys:
        if key not in by_key:
            counts = defaultdict(int)
            for sign in ('pos', 'neut', 'neg', 'unreadable'):
                for bucket, support in _sign_counts(args.vdg_lib_dir, key, sign).items():
                    counts[bucket] += support
            if sum(counts.values()) >= args.min_support:
                by_key[key] = _profile(key, counts)
    for key in (MAIN_TET, MAIN_CARB, LONG_CARB, AROM_CARB):
        if key not in by_key: raise ValueError(f'[ERROR] required profile missing or below support gate: {key}')

    keys = sorted(by_key)
    anchors = sorted(set(tet_keys + carb_keys) & set(by_key))
    ranked = {}
    for key in anchors:
        eligible = []
        for other in keys:
            if other == key or fragment_keys_nested(key, other): continue
            r, _ = _similarity(by_key[key], by_key[other])
            if np.isfinite(r):
                eligible.append((other, r))
        ranked[key] = eligible

    pairs = [(a, b, 'tetrazole_carboxylate') for a in tet_keys for b in carb_keys]
    pairs += [(a, b, 'carboxylate_carboxylate') for a, b in itertools.combinations(carb_keys, 2)]
    rows = []
    for a, b, kind in pairs:
        if a not in by_key or b not in by_key: continue
        r, shared = _similarity(by_key[a], by_key[b])
        nested = fragment_keys_nested(a, b)
        row = dict(kind=kind, frag_A=a, frag_B=b, support_A=by_key[a][4], support_B=by_key[b][4],
                   pearson_r=r, shared_categories=shared, nested=nested)
        for anchor, partner, suffix in ((a, b, 'A_to_B'), (b, a, 'B_to_A')):
            pool = ranked[anchor]
            rank = _rank(pool, partner)[0] if not nested and np.isfinite(r) else 0
            row[f'rank_{suffix}'] = rank or ''
            row[f'pool_{suffix}'] = len(pool)
            row[f'pct_{suffix}'] = rank / len(pool) if rank and pool else ''
        row['top_1pct_both'] = bool(not nested and row['pct_A_to_B'] != '' and
                                    row['pct_A_to_B'] <= .01 and row['pct_B_to_A'] <= .01)
        row['top_5pct_both'] = bool(not nested and row['pct_A_to_B'] != '' and
                                    row['pct_A_to_B'] <= .05 and row['pct_B_to_A'] <= .05)
        rows.append(row)
    _write_tsv(os.path.join(args.outdir, 'pair_rankings.tsv'), rows, list(rows[0]))
    holdout_rows = _holdout_partner_rows(args.vdg_lib_dir, by_key, args.min_support)
    _write_tsv(os.path.join(args.outdir, 'parent_holdout_rankings.tsv'), holdout_rows, list(holdout_rows[0]))

    state, state_support = {}, {}
    for key in lib_keys:
        for sign in SIGNS:
            counts = _sign_counts(args.vdg_lib_dir, key, sign)
            state_support[(key, sign)] = sum(counts.values())
            if state_support[(key, sign)] >= args.min_support:
                state[(key, sign)] = _profile(key, counts)
    state_pairs = [((key, 'neut'), (key, 'neg')) for key in tet_keys + carb_keys]
    state_pairs += [((MAIN_TET, ts), (key, cs)) for key in (MAIN_CARB, LONG_CARB, AROM_CARB)
                    for ts in SIGNS for cs in SIGNS]
    state_rows = []
    state_anchors = {x for pair in state_pairs for x in pair if x in state}
    state_ranked = {}
    for anchor in state_anchors:
        eligible = []
        for other, profile in state.items():
            if other == anchor or (other[0] != anchor[0] and fragment_keys_nested(anchor[0], other[0])):
                continue
            r, _ = _similarity(state[anchor], profile)
            if np.isfinite(r): eligible.append((other, r))
        state_ranked[anchor] = eligible
    for left, right in state_pairs:
        if left not in state or right not in state:
            state_rows.append(dict(frag_A=left[0], sign_A=left[1], support_A=state_support.get(left, 0),
                                   frag_B=right[0], sign_B=right[1], support_B=state_support.get(right, 0),
                                   pearson_r='', shared_categories='',
                                   eligible=False, rank_A_to_B='', pool_A_to_B='', rank_B_to_A='',
                                   pool_B_to_A='', pct_A_to_B='', pct_B_to_A=''))
            continue
        a, b = state[left], state[right]
        r, shared = _similarity(a, b)
        row = dict(frag_A=left[0], sign_A=left[1], support_A=a[4], frag_B=right[0],
                   sign_B=right[1], support_B=b[4], pearson_r=r, shared_categories=shared,
                   eligible=np.isfinite(r))
        for anchor, partner, suffix in ((left, right, 'A_to_B'), (right, left, 'B_to_A')):
            rank, pool_size = _rank(state_ranked[anchor], partner) if np.isfinite(r) else (0, len(state_ranked[anchor]))
            row[f'rank_{suffix}'] = rank or ''
            row[f'pool_{suffix}'] = pool_size
            row[f'pct_{suffix}'] = rank / pool_size if rank and pool_size else ''
        state_rows.append(row)
    _write_tsv(os.path.join(args.outdir, 'charge_comparisons.tsv'), state_rows, list(state_rows[0]))

    fold_rows = []
    for short, long in ((MAIN_CARB, LONG_CARB),):
        for left_fold in (0, 1):
            ca, pa = _fold_counts(args.vdg_lib_dir, short, left_fold)
            cb, pb = _fold_counts(args.vdg_lib_dir, long, 1 - left_fold)
            if pa & pb: raise AssertionError('[ERROR] parent-entry folds overlap')
            a, b = _profile(short, ca), _profile(long, cb)
            r, shared = _similarity(a, b)
            fold_rows.append(dict(short=short, long=long, short_fold=left_fold,
                                  short_entries=len(pa), long_entries=len(pb), support_A=a[4],
                                  support_B=b[4], pearson_r=r, shared_categories=shared,
                                  overlapping_entries=len(pa & pb)))
    _write_tsv(os.path.join(args.outdir, 'nested_parent_folds.tsv'), fold_rows, list(fold_rows[0]))
    primary = next(r for r in rows if r['frag_A'] == MAIN_TET and r['frag_B'] == MAIN_CARB)
    print(f"PRIMARY: r={primary['pearson_r']:.4f} A->B={primary['rank_A_to_B']}/{primary['pool_A_to_B']} "
          f"B->A={primary['rank_B_to_A']}/{primary['pool_B_to_A']} "
          f"top5_both={primary['top_5pct_both']} top1_both={primary['top_1pct_both']}")
    print(f'WROTE: {len(rows)} pairs, {len(state_rows)} state comparisons, '
          f'{len(fold_rows)} nested folds, {len(holdout_rows)} parent holdouts')

if __name__ == '__main__':
    main()
