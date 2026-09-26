'''
Compare AA interaction profiles across CGs to identify bioisosteres (see
docs/bioisostere_identification.md for methodology and output reference).
Input: .npz files from compute_aa_profiles.py, pair mode (default) or
single-AA mode (--single-aa). Metrics: --spearman (default), --pearson.
Writes similarity_matrix_{metric}.npz + top_pairs_{metric}.tsv.
'''
import os
import argparse
import numpy as np
from scipy.stats import spearmanr, pearsonr

from ligand_vdgs.functions.utils import filename_to_smiles
from ligand_vdgs.identify_bioisosteres.common import is_bb_category

# Split points compute_aa_profiles.py puts between CG name and run suffixes
# (norm method, --bb-mode), so CG naming survives new suffixes being added.
_PAIR_MARKER = '_aa_pair_freq'
_SINGLE_MARKER = '_single_aa_freq'

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--npz', nargs='+', metavar='FILE', help='.npz profile files from compute_aa_profiles.py.')
    p.add_argument('--profiles-dir', metavar='DIR', default=None,
                   help='Dir to search recursively for profile .npz files (default: '
                        'outputs/bioisosteres/aa_profiles, used only when --npz is not given). '
                        'Can be combined with --npz by passing --profiles-dir explicitly.')
    p.add_argument('--single-aa', action='store_true',
                   help='Compare 1-D single-AA enrichment vectors instead of 2-D AA-pair '
                        'matrices. Output filenames get an "_single" suffix.')
    p.add_argument('--spearman', action='store_true', help='Compute Spearman ρ similarity (default if no metric given).')
    p.add_argument('--pearson', action='store_true', help='Compute Pearson r similarity.')
    p.add_argument('--min-shared', type=int, default=5, metavar='N',
                   help='Min shared non-NaN entries to report a score; below this the pair is NaN (default: 5).')
    p.add_argument('--min-coverage', type=float, default=0.5, metavar='F',
                   help='[Single-AA] Min fraction of the larger profile\'s AA set that must be in the shared '
                        'intersection (default: 0.5); prevents a sparse CG from defining the comparison space.')
    p.add_argument('--weight-by-count', action='store_true',
                   help='[Single-AA] Scale each AA\'s enrichment by min(1, count / --weight-threshold) '
                        'before comparing, pulling low-confidence entries toward zero.')
    p.add_argument('--weight-threshold', type=int, default=30, metavar='N',
                   help='[Single-AA] Cluster count at which confidence weight reaches 1.0 (default: 30).')
    p.add_argument('--top-n', type=int, default=20, metavar='N', help='Top pairs to print (default: 20); TSV always has all.')
    p.add_argument('--exclude-bb', action='store_true',
                   help='Drop backbone contact categories (bb, <AA>-bb) before comparing, so similarity reflects '
                        'sidechain chemistry alone. Enrichments are not recomputed, so their denominators still '
                        'include the backbone categories.')
    p.add_argument('--outdir', default=os.path.join('outputs', 'bioisosteres', 'similarity'), help='Output directory.')
    return p.parse_args()

def find_npz_files(profiles_dir, stem_filter='_aa_pair_freq'):
    return [os.path.join(root, name) for root, _, files in os.walk(profiles_dir)
            for name in sorted(files) if name.endswith('.npz') and stem_filter in name]

def cg_label(path):
    '''CG name from npz path: cut at the freq marker, reverse the filename encoding.'''
    stem = os.path.splitext(os.path.basename(path))[0]
    for marker in (_PAIR_MARKER, _SINGLE_MARKER):
        head, sep, _ = stem.partition(marker)
        if sep:
            return filename_to_smiles(head)
    return filename_to_smiles(stem)

def load_profile(path):
    d = np.load(path)
    return d['matrix'], list(d['aa_labels'])

def load_single_profile(path):
    '''Returns (enrichments, aa_labels, counts, total_count).'''
    d = np.load(path)
    return d['enrichments'].astype(float), list(d['aa_labels']), d['counts'].astype(int), int(d['total_count'])

def profile_provenance(path):
    return str(np.load(path)['bb_mode'])

def drop_bb_categories(labels, *arrays):
    keep = [i for i, l in enumerate(labels) if not is_bb_category(l)]
    return ([labels[i] for i in keep],) + tuple(a[keep] for a in arrays)

def aligned_upper_tri(matrix_a, labels_a, matrix_b, labels_b):
    '''Align finite shared-AA upper triangles, including the diagonal.'''
    idx_b = {aa: i for i, aa in enumerate(labels_b)}
    shared = [aa for aa in labels_a if aa in idx_b]
    if len(shared) < 2:
        return np.array([]), np.array([])
    ia, ib = [labels_a.index(aa) for aa in shared], [idx_b[aa] for aa in shared]
    rows, cols = np.triu_indices(len(shared), k=0)
    va, vb = matrix_a[np.ix_(ia, ia)][rows, cols], matrix_b[np.ix_(ib, ib)][rows, cols]
    mask = np.isfinite(va) & np.isfinite(vb)
    return va[mask], vb[mask]

def aligned_vectors(labels_a, enrichments_a, counts_a, labels_b, enrichments_b, counts_b,
                    min_coverage=0.0, weight_threshold=0):
    '''Align finite shared AAs; optionally filter by coverage and weight by count.'''
    idx_b = {aa: i for i, aa in enumerate(labels_b)}
    shared = [aa for aa in labels_a if aa in idx_b]
    if len(shared) < 2 or len(shared) / max(len(labels_a), len(labels_b)) < min_coverage:
        return np.array([]), np.array([]), []
    ia, ib = np.array([labels_a.index(aa) for aa in shared]), np.array([idx_b[aa] for aa in shared])
    va, vb = enrichments_a[ia].astype(float).copy(), enrichments_b[ib].astype(float).copy()
    mask = np.isfinite(va) & np.isfinite(vb)
    if not mask.any():
        return np.array([]), np.array([]), []
    va, vb = va[mask], vb[mask]
    if weight_threshold > 0:
        va *= np.minimum(1.0, counts_a[ia[mask]] / weight_threshold)
        vb *= np.minimum(1.0, counts_b[ib[mask]] / weight_threshold)
    return va, vb, [aa for aa, keep in zip(shared, mask) if keep]

def _spearman(va, vb, min_shared):
    if len(va) < min_shared:
        return np.nan, np.nan
    r, p = spearmanr(va, vb)
    return float(r), float(p)

def _pearson(va, vb, min_shared):
    if len(va) < min_shared:
        return np.nan, np.nan
    r, p = pearsonr(va, vb)
    return float(r), float(p)

def build_sim_matrix(profiles, metric, min_shared, mode='pair', min_coverage=0.0, weight_threshold=0):
    '''Returns (sim_matrix, pval_matrix) of Spearman/Pearson correlations.'''
    n = len(profiles)
    sim, pval = np.full((n, n), np.nan), np.full((n, n), np.nan)
    np.fill_diagonal(sim, 1.0)
    np.fill_diagonal(pval, 0.0)
    corr = _spearman if metric == 'spearman' else _pearson
    for i in range(n):
        for j in range(i + 1, n):
            if mode == 'single':
                _, enr_i, lbl_i, cnt_i, _ = profiles[i]
                _, enr_j, lbl_j, cnt_j, _ = profiles[j]
                va, vb, _ = aligned_vectors(lbl_i, enr_i, cnt_i, lbl_j, enr_j, cnt_j, min_coverage, weight_threshold)
            else:
                va, vb = aligned_upper_tri(profiles[i][1], profiles[i][2], profiles[j][1], profiles[j][2])
            sim[i, j], pval[i, j] = sim[j, i], pval[j, i] = corr(va, vb, min_shared)
    return sim, pval

def fast_single_aa_pearson(profiles, min_shared=5, min_coverage=0.5):
    '''Vectorized pairwise-complete Pearson for the single-AA Figure 1/2 matrices.'''
    labels = sorted({label for _, _, aa, _, _ in profiles for label in aa})
    values = np.full((len(profiles), len(labels)), np.nan)
    values[np.repeat(np.arange(len(profiles)), [len(aa) for _, _, aa, _, _ in profiles]),
           np.concatenate([[labels.index(label) for label in aa] for _, _, aa, _, _ in profiles])] = (
        np.concatenate([enrichment for _, enrichment, _, _, _ in profiles]))
    finite = np.isfinite(values)
    mask = finite.astype(float)
    x = np.where(finite, values, 0.0)
    n = mask @ mask.T
    count = mask.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        sx = x @ mask.T
        var = (x * x) @ mask.T - sx * sx / n
        sim = (x @ x.T - sx * sx.T / n) / np.sqrt(var * var.T)
    sim[(n < min_shared) | (n / np.maximum(count[:, None], count[None, :]) < min_coverage) | (var <= 0) | (var.T <= 0)] = np.nan
    np.clip(sim, -1.0, 1.0, out=sim)
    np.fill_diagonal(sim, 1.0)
    return sim

def write_pairs_tsv(sim, labels, outpath, profiles=None, mode='pair', pval=None):
    '''Write finite pairs in descending score order and return them.'''
    pairs = sorted(((sim[i, j], i, j) for i in range(len(labels)) for j in range(i + 1, len(labels))
                   if np.isfinite(sim[i, j])), reverse=True)
    single = mode == 'single' and profiles is not None
    total_counts = {p[0]: int(p[4]) for p in profiles} if single else {}
    with open(outpath, 'w') as f:
        f.write('\t'.join(['rank', 'cg1', 'cg2']
                          + (['n_clusters_cg1', 'n_clusters_cg2'] if single else [])
                          + ['score', 'pvalue']) + '\n')
        for rank, (score, i, j) in enumerate(pairs, 1):
            f.write('\t'.join([str(rank), labels[i], labels[j]]
                              + ([str(total_counts.get(labels[i], 0)), str(total_counts.get(labels[j], 0))] if single else [])
                              + [f'{score:.4f}', f'{pval[i,j]:.3e}']) + '\n')
    return pairs

def run_metric(metric, profiles, labels, args, mode='pair'):
    print(f'\n--- {metric} ---')
    sim, pval = build_sim_matrix(
        profiles, metric, args.min_shared, mode=mode,
        min_coverage=args.min_coverage if mode == 'single' else 0.0,
        weight_threshold=args.weight_threshold if mode == 'single' and args.weight_by_count else 0)

    tag = f'{metric}_single' if mode == 'single' else metric
    np.savez(os.path.join(args.outdir, f'similarity_matrix_{tag}.npz'),
             matrix=sim, cg_labels=np.array(labels, dtype=str), pval=pval)

    pairs = write_pairs_tsv(sim, labels, os.path.join(args.outdir, f'top_pairs_{tag}.tsv'),
                            profiles=profiles if mode == 'single' else None, mode=mode, pval=pval)

    top = min(args.top_n, len(pairs))
    if pairs:
        w = max(len(lbl) for lbl in labels)
        print(f'Top {top} of {len(pairs)} pairs by {"Spearman ρ" if metric == "spearman" else "Pearson r"}:')
        for rank, (score, i, j) in enumerate(pairs[:top], 1):
            parts = [f'  {rank:3d}. {labels[i]:{w}s}  <->  {labels[j]:{w}s}  score = {score:.3f}']
            if np.isfinite(pval[i, j]):
                parts.append(f'  p={pval[i,j]:.2e}')
            if mode == 'single':
                parts.append(f'  [N={int(profiles[i][4])}/{int(profiles[j][4])}]')
            print(''.join(parts))
    else:
        print(f'[WARNING] No pairs with sufficient shared data for {metric}. Try lowering --min-shared or --min-coverage.')

    print(f'  similarity_matrix_{tag}.npz')
    print(f'  top_pairs_{tag}.tsv')

def main():
    args = parse_args()
    mode = 'single' if args.single_aa else 'pair'
    stem_filter = '_single_aa_freq' if mode == 'single' else '_aa_pair_freq'

    npz_files = list(args.npz or [])
    # Fall back to the default profiles-dir tree only when --npz was not given;
    # an explicit --profiles-dir always searches, per the "combine with --npz" doc.
    profiles_dir = args.profiles_dir
    if profiles_dir is None and not args.npz:
        profiles_dir = os.path.join('outputs', 'bioisosteres', 'aa_profiles')
    if profiles_dir:
        npz_files.extend(find_npz_files(profiles_dir, stem_filter))
    npz_files = sorted(set(npz_files))

    if len(npz_files) < 2:
        print(f'[ERROR] Need at least 2 .npz profile files. Use --npz or --profiles-dir. '
              f'(Mode: {mode}, searching for files containing "{stem_filter}")')
        return

    print(f'Loading {len(npz_files)} profiles ({mode} mode)...')
    profiles, provenances = [], set()
    for path in npz_files:
        label = cg_label(path)
        provenances.add(profile_provenance(path))
        if mode == 'single':
            enr, aa_labels, counts, total_count = load_single_profile(path)
            if args.exclude_bb:
                aa_labels, enr, counts = drop_bb_categories(aa_labels, enr, counts)
            profiles.append((label, enr, aa_labels, counts, total_count))
            print(f'  {label}  ({len(aa_labels)} categories, N={total_count})')
        else:
            matrix, aa_labels = load_profile(path)
            if args.exclude_bb:
                keep = [i for i, l in enumerate(aa_labels) if not is_bb_category(l)]
                aa_labels = [aa_labels[i] for i in keep]
                matrix = matrix[np.ix_(keep, keep)]
            profiles.append((label, matrix, aa_labels))
            print(f'  {label}  ({len(aa_labels)} categories with pair data)')

    if len(provenances) > 1:
        print(f'[WARNING] Profiles were built under more than one category set '
              f'{sorted(str(p) for p in provenances)}. Their labels align but their denominators differ, so the '
              f'similarities are not meaningful. Recompute them with one --bb-mode.')

    labels = [p[0] for p in profiles]
    seen = {}
    for i, lbl in enumerate(labels):
        if lbl in seen:
            print(f'[WARNING] Duplicate label "{lbl}" at indices {seen[lbl]} and {i}. Consider using files from a single norm method only.')
        seen[lbl] = i

    metrics = [metric for metric in ('spearman', 'pearson') if getattr(args, metric)] or ['spearman']

    os.makedirs(args.outdir, exist_ok=True)
    print(f'\nMode:    {mode}')
    print(f'Metrics: {metrics}')
    print(f'outdir:  {args.outdir}')
    if mode == 'single':
        print(f'min-coverage:     {args.min_coverage}')
        print(f'weight-by-count:  {args.weight_by_count}' + (f'  (threshold={args.weight_threshold})' if args.weight_by_count else ''))

    for metric in metrics:
        run_metric(metric, profiles, labels, args, mode=mode)

    print(f'\nDone. Outputs in {args.outdir}/')

if __name__ == '__main__':
    main()
