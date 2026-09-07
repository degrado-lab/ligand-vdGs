'''
Compare AA interaction profiles across CGs to identify bioisosteres. See
docs/bioisostere_identification.md for the full methodology (modes, metrics,
quality filters) and output reference. visualize_bioisosteres.py handles all
plotting downstream of this script's output.

Input: .npz files from compute_aa_profiles.py, in pair mode (default,
*_aa_pair_freq*.npz) or single-AA mode (--single-aa, *_single_aa_freq*.npz;
recommended starting point).

Metrics: --spearman (default if none given), --pearson, --jaccard.

Outputs (per active metric m, in --outdir):
  similarity_matrix_{m}[_single].npz  — similarity matrix + CG labels (+ p-values)
  top_pairs_{m}[_single].tsv          — all finite pairs ranked descending

Examples:
    python compare_aa_profiles.py --profiles-dir outputs/bioisosteres/aa_profiles/ --single-aa
    python compare_aa_profiles.py --npz cg1.npz cg2.npz cg3.npz --pearson --jaccard
'''
import os
import argparse
import numpy as np
from scipy.stats import spearmanr, pearsonr

from ligand_vdgs.functions.utils import filename_to_smiles
from ligand_vdgs.identify_bioisosteres.common import is_bb_category

# Markers compute_aa_profiles.py puts between the CG name and the run's suffixes.
# Matched as a split point rather than a fixed suffix list, so adding a suffix
# there (norm method, --bb-mode) does not silently break CG naming.
_PAIR_MARKER = '_aa_pair_freq'
_SINGLE_MARKER = '_single_aa_freq'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--npz', nargs='+', metavar='FILE',
                   help='One or more .npz profile files (from compute_aa_profiles.py).')
    p.add_argument('--profiles-dir', metavar='DIR', default=None,
                   help='Directory to search recursively for profile .npz files '
                        '(default: outputs/bioisosteres/aa_profiles, used only when '
                        '--npz is not given). '
                        'In pair mode looks for *_aa_pair_freq*.npz; '
                        'in single-AA mode looks for *_single_aa_freq*.npz. '
                        'Can be combined with --npz by passing --profiles-dir explicitly.')
    p.add_argument('--single-aa', action='store_true',
                   help='Compare 1-D single-AA enrichment vectors instead of 2-D '
                        'AA-pair matrices. Simpler and more interpretable. Output '
                        'filenames are suffixed with "_single".')
    p.add_argument('--spearman', action='store_true',
                   help='Compute Spearman ρ similarity (default when no metric flag given).')
    p.add_argument('--pearson', action='store_true',
                   help='Compute Pearson r similarity.')
    p.add_argument('--jaccard', action='store_true',
                   help='Compute top-k Jaccard similarity.')
    p.add_argument('--jaccard-k', type=int, default=10, metavar='K',
                   help='Number of top-enriched AA (pairs) to use for Jaccard (default: 10).')
    p.add_argument('--min-shared', type=int, default=5, metavar='N',
                   help='Minimum shared non-NaN AA (pair) entries required to report a '
                        'similarity score; pairs below this are NaN (default: 5).')
    p.add_argument('--min-coverage', type=float, default=0.5, metavar='F',
                   help='[Single-AA mode] Minimum fraction of the larger profile\'s AA set '
                        'that must appear in the shared intersection (default: 0.5). '
                        'Prevents a sparse CG from defining the comparison space.')
    p.add_argument('--weight-by-count', action='store_true',
                   help='[Single-AA mode] Scale each AA\'s enrichment by '
                        'min(1, count / --weight-threshold) before computing similarity. '
                        'Entries backed by few clusters are pulled toward zero.')
    p.add_argument('--weight-threshold', type=int, default=30, metavar='N',
                   help='[Single-AA mode] Cluster count at which confidence weight '
                        'reaches 1.0 (default: 30). Used only with --weight-by-count.')
    p.add_argument('--top-n', type=int, default=20, metavar='N',
                   help='Number of top pairs to print to console (default: 20). '
                        'All pairs are always written to the TSV.')
    p.add_argument('--jaccard-positive-only', action='store_true',
                   help='For Jaccard, restrict top-k selection to AAs with enrichment > 0 '
                        'before selecting top-k. Prevents universally hydrophobic/aromatic '
                        'AAs from dominating the intersection across all CGs.')
    p.add_argument('--score-cutoff', type=float, default=None, metavar='F',
                   help='Only write pairs with score >= this value to the TSV and console '
                        '(default: no cutoff). More principled than --top-n.')
    p.add_argument('--exclude-bb', action='store_true',
                   help='Drop backbone contact categories (bb, <AA>-bb) before '
                        'comparing, so similarity is decided by sidechain chemistry '
                        'alone. The enrichments themselves are not recomputed, so the '
                        'values still carry the backbone categories in their '
                        'denominator.')
    p.add_argument('--outdir', default=os.path.join('outputs', 'bioisosteres', 'similarity'),
                   help='Output directory (default: outputs/bioisosteres/similarity/).')
    return p.parse_args()


def find_npz_files(profiles_dir, stem_filter='_aa_pair_freq'):
    found = []
    for root, _, files in os.walk(profiles_dir):
        for f in sorted(files):
            if f.endswith('.npz') and stem_filter in f:
                found.append(os.path.join(root, f))
    return found


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
    '''Load 1-D single-AA profile.
    Returns (enrichments, aa_labels, counts, total_count).'''
    d = np.load(path)
    return (d['enrichments'].astype(float),
            list(d['aa_labels']),
            d['counts'].astype(int),
            int(d['total_count']))


def profile_provenance(path):
    '''bb_mode a profile was built under.

    Two profiles built under different backbone modes have different label sets
    and different denominators, so aligning them on shared labels quietly compares
    incomparable numbers -- the alignment succeeds, which is what makes it worth
    checking rather than trusting.
    '''
    d = np.load(path)
    return str(d['bb_mode'])


def drop_bb_categories(labels, *arrays):
    '''Filter labels (and any parallel arrays) down to non-backbone categories.'''
    keep = [i for i, l in enumerate(labels) if not is_bb_category(l)]
    return ([labels[i] for i in keep],) + tuple(a[keep] for a in arrays)


def aligned_upper_tri(matrix_a, labels_a, matrix_b, labels_b):
    '''
    Return (va, vb, pair_labels) for the shared-AA upper-triangle entries
    (including diagonal) where both values are finite.
    pair_labels[i] is e.g. "ASP_HIS" for the i-th entry.
    '''
    idx_a = {aa: i for i, aa in enumerate(labels_a)}
    idx_b = {aa: i for i, aa in enumerate(labels_b)}
    shared = [aa for aa in labels_a if aa in idx_b]
    if len(shared) < 2:
        return np.array([]), np.array([]), []
    ia = [idx_a[aa] for aa in shared]
    ib = [idx_b[aa] for aa in shared]
    sub_a = matrix_a[np.ix_(ia, ia)]
    sub_b = matrix_b[np.ix_(ib, ib)]
    n = len(shared)
    rows, cols = np.triu_indices(n, k=0)
    all_labels = [f'{shared[r]}_{shared[c]}' for r, c in zip(rows, cols)]
    va = sub_a[rows, cols]
    vb = sub_b[rows, cols]
    mask = np.isfinite(va) & np.isfinite(vb)
    keep = np.where(mask)[0]
    return va[mask], vb[mask], [all_labels[i] for i in keep]


def aligned_vectors(labels_a, enrichments_a, counts_a,
                    labels_b, enrichments_b, counts_b,
                    min_coverage=0.0, weight_threshold=0):
    '''
    Return (va, vb, aa_labels) aligned on shared AAs where both are finite.

    Coverage filter: if the shared set covers less than min_coverage of the
    larger profile's AA set, returns empty arrays (pair treated as incomparable).

    Confidence weighting: if weight_threshold > 0, each AA's enrichment is
    scaled by min(1, count / weight_threshold) independently for each CG.
    Entries with fewer supporting clusters are pulled toward zero and
    contribute less to downstream correlation.
    '''
    idx_a = {aa: i for i, aa in enumerate(labels_a)}
    idx_b = {aa: i for i, aa in enumerate(labels_b)}
    shared = [aa for aa in labels_a if aa in idx_b]

    if len(shared) < 2:
        return np.array([]), np.array([]), []

    if min_coverage > 0:
        max_n = max(len(labels_a), len(labels_b))
        if len(shared) / max_n < min_coverage:
            return np.array([]), np.array([]), []

    ia = np.array([idx_a[aa] for aa in shared])
    ib = np.array([idx_b[aa] for aa in shared])
    va = enrichments_a[ia].astype(float).copy()
    vb = enrichments_b[ib].astype(float).copy()

    mask = np.isfinite(va) & np.isfinite(vb)
    if not mask.any():
        return np.array([]), np.array([]), []

    ia_m = ia[mask]
    ib_m = ib[mask]
    va, vb = va[mask], vb[mask]
    shared_out = [shared[k] for k in np.where(mask)[0]]

    if weight_threshold > 0:
        ca = counts_a[ia_m].astype(float)
        cb = counts_b[ib_m].astype(float)
        va = va * np.minimum(1.0, ca / weight_threshold)
        vb = vb * np.minimum(1.0, cb / weight_threshold)

    return va, vb, shared_out


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


def _jaccard(va, vb, pair_labels, k, min_shared, positive_only=False):
    # np.argsort(v)[-k:] is the whole array at k=0, which would score 1.0.
    if k < 1:
        raise ValueError(f'--jaccard-k must be >= 1, got {k}')
    if positive_only:
        mask = (va > 0) & (vb > 0)
        va = va[mask]
        vb = vb[mask]
        pair_labels = [pair_labels[i] for i in np.where(mask)[0]]
    if len(va) < min_shared or len(va) < k:
        return np.nan, 0, []
    top_a = set(np.argsort(va)[-k:])
    top_b = set(np.argsort(vb)[-k:])
    inter = top_a & top_b
    union = top_a | top_b
    score = len(inter) / len(union) if union else np.nan
    shared = [pair_labels[i] for i in sorted(inter, key=lambda i: -va[i])]
    return float(score), len(inter), shared


def build_sim_matrix(profiles, metric, min_shared, mode='pair',
                     min_coverage=0.0, weight_threshold=0, jaccard_k=10,
                     jaccard_positive_only=False):
    '''
    Returns (sim_matrix, pval_matrix, extras).
    pval_matrix: p-values for spearman/pearson; NaN for jaccard.
    extras: dict (i,j) -> (n_shared, shared_labels) for jaccard; empty otherwise.
    '''
    n = len(profiles)
    sim = np.full((n, n), np.nan)
    pval = np.full((n, n), np.nan)
    np.fill_diagonal(sim, 1.0)
    np.fill_diagonal(pval, 0.0)
    extras = {}

    for i in range(n):
        for j in range(i + 1, n):
            if mode == 'single':
                _, enr_i, lbl_i, cnt_i, _ = profiles[i]
                _, enr_j, lbl_j, cnt_j, _ = profiles[j]
                va, vb, plabels = aligned_vectors(
                    lbl_i, enr_i, cnt_i, lbl_j, enr_j, cnt_j,
                    min_coverage, weight_threshold)
            else:
                va, vb, plabels = aligned_upper_tri(
                    profiles[i][1], profiles[i][2],
                    profiles[j][1], profiles[j][2])

            if metric == 'spearman':
                r, p = _spearman(va, vb, min_shared)
                sim[i, j] = sim[j, i] = r
                pval[i, j] = pval[j, i] = p
            elif metric == 'pearson':
                r, p = _pearson(va, vb, min_shared)
                sim[i, j] = sim[j, i] = r
                pval[i, j] = pval[j, i] = p
            elif metric == 'jaccard':
                score, n_sh, shared = _jaccard(va, vb, plabels, jaccard_k, min_shared,
                                               positive_only=jaccard_positive_only)
                sim[i, j] = sim[j, i] = score
                extras[(i, j)] = extras[(j, i)] = (n_sh, shared)

    return sim, pval, extras


def write_pairs_tsv(sim, labels, extras, outpath, metric,
                    profiles=None, mode='pair', pval=None, score_cutoff=None):
    '''Write all finite pairs ranked by score (descending).
    In single-AA mode, includes n_clusters columns for each CG.
    For spearman/pearson, includes a pvalue column when pval is provided.
    Pairs below score_cutoff are excluded when score_cutoff is set.
    '''
    n = len(labels)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            score = sim[i, j]
            if np.isfinite(score):
                if score_cutoff is None or score >= score_cutoff:
                    pairs.append((score, i, j))
    pairs.sort(reverse=True)

    has_pval = pval is not None and metric in ('spearman', 'pearson')

    with open(outpath, 'w') as f:
        if mode == 'single' and profiles is not None:
            total_counts = {p[0]: int(p[4]) for p in profiles}
            if metric == 'jaccard':
                f.write('rank\tcg1\tcg2\tn_clusters_cg1\tn_clusters_cg2\t'
                        'score\tn_shared_aas\tshared_aas\n')
                for rank, (score, i, j) in enumerate(pairs, 1):
                    n_sh, shared = extras.get((i, j), (0, []))
                    nc1 = total_counts.get(labels[i], 0)
                    nc2 = total_counts.get(labels[j], 0)
                    f.write(f'{rank}\t{labels[i]}\t{labels[j]}\t{nc1}\t{nc2}\t'
                            f'{score:.4f}\t{n_sh}\t{",".join(shared)}\n')
            else:
                pval_col = '\tpvalue' if has_pval else ''
                f.write(f'rank\tcg1\tcg2\tn_clusters_cg1\tn_clusters_cg2\tscore{pval_col}\n')
                for rank, (score, i, j) in enumerate(pairs, 1):
                    nc1 = total_counts.get(labels[i], 0)
                    nc2 = total_counts.get(labels[j], 0)
                    pval_str = f'\t{pval[i,j]:.3e}' if has_pval else ''
                    f.write(f'{rank}\t{labels[i]}\t{labels[j]}\t{nc1}\t{nc2}\t'
                            f'{score:.4f}{pval_str}\n')
        else:
            if metric == 'jaccard':
                f.write('rank\tcg1\tcg2\tscore\tn_shared\tshared_pairs\n')
                for rank, (score, i, j) in enumerate(pairs, 1):
                    n_shared, shared = extras.get((i, j), (0, []))
                    f.write(f'{rank}\t{labels[i]}\t{labels[j]}\t{score:.4f}'
                            f'\t{n_shared}\t{",".join(shared)}\n')
            else:
                pval_col = '\tpvalue' if has_pval else ''
                f.write(f'rank\tcg1\tcg2\tscore{pval_col}\n')
                for rank, (score, i, j) in enumerate(pairs, 1):
                    pval_str = f'\t{pval[i,j]:.3e}' if has_pval else ''
                    f.write(f'{rank}\t{labels[i]}\t{labels[j]}\t{score:.4f}{pval_str}\n')

    return pairs


def run_metric(metric, profiles, labels, args, mode='pair'):
    print(f'\n--- {metric} ---')

    weight_threshold = args.weight_threshold if (mode == 'single' and args.weight_by_count) else 0
    min_coverage = args.min_coverage if mode == 'single' else 0.0

    sim, pval, extras = build_sim_matrix(
        profiles, metric, args.min_shared,
        mode=mode,
        min_coverage=min_coverage,
        weight_threshold=weight_threshold,
        jaccard_k=args.jaccard_k,
        jaccard_positive_only=args.jaccard_positive_only)

    tag = f'{metric}_single' if mode == 'single' else metric

    npz_data = dict(matrix=sim, cg_labels=np.array(labels, dtype=str))
    if metric in ('spearman', 'pearson'):
        npz_data['pval'] = pval
    np.savez(os.path.join(args.outdir, f'similarity_matrix_{tag}.npz'), **npz_data)

    tsv_path = os.path.join(args.outdir, f'top_pairs_{tag}.tsv')
    pairs = write_pairs_tsv(sim, labels, extras, tsv_path, metric,
                            profiles=profiles if mode == 'single' else None,
                            mode=mode,
                            pval=pval if metric in ('spearman', 'pearson') else None,
                            score_cutoff=args.score_cutoff)

    top = min(args.top_n, len(pairs))
    if pairs:
        w = max(len(lbl) for lbl in labels)
        metric_label = ('Spearman ρ' if metric == 'spearman'
                        else 'Pearson r' if metric == 'pearson'
                        else f'Jaccard (k={args.jaccard_k})')
        print(f'Top {top} of {len(pairs)} pairs by {metric_label}:')
        for rank, (score, i, j) in enumerate(pairs[:top], 1):
            line = f'  {rank:3d}. {labels[i]:{w}s}  <->  {labels[j]:{w}s}  score = {score:.3f}'
            if metric in ('spearman', 'pearson') and np.isfinite(pval[i, j]):
                line += f'  p={pval[i,j]:.2e}'
            if mode == 'single':
                nc_i = int(profiles[i][4])
                nc_j = int(profiles[j][4])
                line += f'  [N={nc_i}/{nc_j}]'
            if metric == 'jaccard':
                n_sh, shared = extras.get((i, j), (0, []))
                label_word = 'AAs' if mode == 'single' else 'AA pairs'
                line += (f'  [{n_sh} shared {label_word}: '
                         f'{", ".join(shared[:5])}{"..." if len(shared) > 5 else ""}]')
            print(line)
    else:
        print(f'[WARNING] No pairs with sufficient shared data for {metric}. '
              'Try lowering --min-shared or --min-coverage.')

    print(f'  similarity_matrix_{tag}.npz')
    print(f'  top_pairs_{tag}.tsv')


def main():
    args = parse_args()

    mode = 'single' if args.single_aa else 'pair'
    stem_filter = '_single_aa_freq' if mode == 'single' else '_aa_pair_freq'

    npz_files = list(args.npz or [])
    # Only fall back to the default profiles-dir tree when the user didn't
    # explicitly point us at files with --npz; an explicit --profiles-dir
    # always searches, matching the documented "combine with --npz" behavior.
    profiles_dir = args.profiles_dir
    if profiles_dir is None and not args.npz:
        profiles_dir = os.path.join('outputs', 'bioisosteres', 'aa_profiles')
    if profiles_dir:
        npz_files.extend(find_npz_files(profiles_dir, stem_filter))
    npz_files = sorted(set(npz_files))

    if len(npz_files) < 2:
        print(f'[ERROR] Need at least 2 .npz profile files. '
              f'Use --npz or --profiles-dir. '
              f'(Mode: {mode}, searching for files containing "{stem_filter}")')
        return

    print(f'Loading {len(npz_files)} profiles ({mode} mode)...')
    profiles = []
    provenances = set()
    for path in npz_files:
        label = cg_label(path)
        if mode == 'single':
            provenances.add(profile_provenance(path))
            enr, aa_labels, counts, total_count = load_single_profile(path)
            if args.exclude_bb:
                aa_labels, enr, counts = drop_bb_categories(aa_labels, enr, counts)
            profiles.append((label, enr, aa_labels, counts, total_count))
            print(f'  {label}  ({len(aa_labels)} categories, N={total_count})')
        else:
            provenances.add(profile_provenance(path))
            matrix, aa_labels = load_profile(path)
            if args.exclude_bb:
                keep = [i for i, l in enumerate(aa_labels) if not is_bb_category(l)]
                aa_labels = [aa_labels[i] for i in keep]
                matrix = matrix[np.ix_(keep, keep)]
            profiles.append((label, matrix, aa_labels))
            print(f'  {label}  ({len(aa_labels)} categories with pair data)')

    if len(provenances) > 1:
        print(f'[WARNING] Profiles were built under more than one category set '
              f'{sorted(str(p) for p in provenances)}. Their labels align but their '
              f'denominators differ, so the similarities are not meaningful. '
              f'Recompute them with one --bb-mode.')

    labels = [p[0] for p in profiles]

    seen = {}
    for i, lbl in enumerate(labels):
        if lbl in seen:
            print(f'[WARNING] Duplicate label "{lbl}" at indices {seen[lbl]} and {i}. '
                  'Consider using files from a single norm method only.')
        seen[lbl] = i

    any_explicit = args.spearman or args.pearson or args.jaccard
    metrics = []
    if args.spearman or not any_explicit:
        metrics.append('spearman')
    if args.pearson:
        metrics.append('pearson')
    if args.jaccard:
        metrics.append('jaccard')

    os.makedirs(args.outdir, exist_ok=True)
    print(f'\nMode:    {mode}')
    print(f'Metrics: {metrics}')
    print(f'outdir:  {args.outdir}')
    if mode == 'single':
        print(f'min-coverage:     {args.min_coverage}')
        print(f'weight-by-count:  {args.weight_by_count}'
              + (f'  (threshold={args.weight_threshold})' if args.weight_by_count else ''))

    for metric in metrics:
        run_metric(metric, profiles, labels, args, mode=mode)

    print(f'\nDone. Outputs in {args.outdir}/')


if __name__ == '__main__':
    main()
