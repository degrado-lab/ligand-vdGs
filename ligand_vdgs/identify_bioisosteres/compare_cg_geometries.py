"""
Compare fragment libraries for bioisostere identification (AA enrichment +
CG geometry), unifying compute_aa_profiles.py + compare_aa_profiles.py.

For each pair of fragment libraries (A, B) and each shared AA bucket:
  enrichment_A/B  — log(observed / expected) cluster counts for this AA
                    type, size-normalised by background frequency. Positive
                    = interacts with this AA more than expected by chance.
  min_dist        — min pairwise CG centroid distance (Å) after Kabsch-
                    aligning backbone + CG centroid of each vdG pair. The CG
                    centroid is included in the fit because a 3-atom (N, CA,
                    C) backbone superposition, while well determined, is very
                    sensitive to coordinate noise; adding the centroid ties
                    the alignment to the CG's own position. For AA buckets
                    with repeated residue types, all permutations of
                    identical residues are tried and the best kept. No
                    atom-to-atom correspondence between fragments is needed.
  frac_A/B_matched — fraction of each library's nr vdGs with a match
                    (< threshold) in the other library.

Usage:
    python compare_cg_geometries.py --vdg-lib-dir /path/to/frag_lib \\
        --frags "cnnnn" "CC(=O)[O-]" --outdir outputs/bioisosteres

    # All fragment pairs; AA enrichment only (faster):
    python compare_cg_geometries.py --vdg-lib-dir /path/to/frag_lib \\
        --outdir outputs/bioisosteres --skip-geometry
"""
import os
import sys
import math
import argparse
import itertools

import numpy as np

from ligand_vdgs.functions.vdg_struct_utils import BB_LABELS
from ligand_vdgs.functions import utils, Frags
from ligand_vdgs.functions.vdg_npz_utils import (load_vdg_bucket, aa_perm_indices,
                                                 resolve_fragment_alias)

from ligand_vdgs.identify_bioisosteres.common import calc_single_aa_propensities


# ---------------------------------------------------------------------------
# AA enrichment
# ---------------------------------------------------------------------------

def compute_single_aa_enrichments(nr_vdgs_dir):
    """Log-enrichment per contact category (bg_weighted norm). {} if no single-AA vdGs.

    bb_mode is pinned explicitly rather than left to the default: the category set
    is part of what the numbers mean, and this branch reports them next to
    geometry, where a silently-changing label set would be hard to notice.
    """
    try:
        enrich, _ = calc_single_aa_propensities(nr_vdgs_dir, 'bg_weighted',
                                                bb_mode='per-residue')
    except FileNotFoundError:
        return {}
    return enrich


# ---------------------------------------------------------------------------
# Geometric comparison
# ---------------------------------------------------------------------------

def compare_bucket_geometry(cg_A, bb_A, cg_B, bb_B, aa_bucket_parts,
                             threshold, max_per_lib=1000):
    """
    Compare CG centroid distances between two libraries for one AA bucket.

    cg_A/B [N,n_cg,3], bb_A/B [N,num_vdms,3,3]; aa_bucket_parts is the bucket
    label split on '_' and defines which residue slots may be permuted. Each library is
    subsampled to at most max_per_lib nr vdGs before the O(N_A x N_B)
    Kabsch comparison. threshold (Å) sets the "matched" cutoff.

    Returns (min_dist, frac_A_matched, frac_B_matched).
    """
    cg_A = np.asarray(cg_A, np.float32)
    bb_A = np.asarray(bb_A, np.float32)
    cg_B = np.asarray(cg_B, np.float32)
    bb_B = np.asarray(bb_B, np.float32)

    N_A, N_B    = cg_A.shape[0], cg_B.shape[0]
    num_vdms    = bb_A.shape[1]
    aa_perms    = aa_perm_indices(list(aa_bucket_parts))

    # Subsample to keep runtime feasible
    rng = np.random.default_rng(0)
    if N_A > max_per_lib:
        idx = rng.choice(N_A, max_per_lib, replace=False)
        cg_A, bb_A = cg_A[idx], bb_A[idx]
        N_A = max_per_lib
    if N_B > max_per_lib:
        idx = rng.choice(N_B, max_per_lib, replace=False)
        cg_B, bb_B = cg_B[idx], bb_B[idx]
        N_B = max_per_lib

    cg_A_cents = cg_A.mean(axis=1)     # [N_A, 3]
    cg_B_cents = cg_B.mean(axis=1)     # [N_B, 3]

    # Build reference vectors for B (fixed)
    ref_B = np.empty((N_B, num_vdms * 3 + 1, 3), dtype=np.float32)
    ref_B[:, :num_vdms * 3] = bb_B.reshape(N_B, -1, 3)
    ref_B[:, -1]             = cg_B_cents

    # Pairwise minimum distances [N_A, N_B] — minimised over AA permutations
    dists_min = np.full((N_A, N_B), np.inf, dtype=np.float32)

    for perm in aa_perms:
        # Build A's references with this residue permutation
        bb_A_perm = bb_A[:, perm, :, :]        # [N_A, num_vdms, 3, 3]
        ref_A = np.empty((N_A, num_vdms * 3 + 1, 3), dtype=np.float32)
        ref_A[:, :num_vdms * 3] = bb_A_perm.reshape(N_A, -1, 3)
        ref_A[:, -1]             = cg_A_cents

        for i in range(N_A):
            # Kabsch fast path: one X [N_ref, 3] against many Y [N_B, N_ref, 3]
            # R, t satisfy: ref_A[i] @ R[j] + t[j] ≈ ref_B[j] (right-multiply)
            R, t, _ = utils.kabsch(ref_A[i], ref_B)    # R [N_B,3,3], t [N_B,3]

            # Transform CG centroid of A[i] into each B[j]'s frame
            cg_A_aligned = cg_A_cents[i] @ R + t        # [N_B, 3]
            dists = np.linalg.norm(cg_A_aligned - cg_B_cents, axis=1)  # [N_B]
            np.minimum(dists_min[i], dists, out=dists_min[i])

    min_dist       = float(dists_min.min())
    frac_A_matched = float((dists_min.min(axis=1) < threshold).mean())
    frac_B_matched = float((dists_min.min(axis=0) < threshold).mean())

    return min_dist, frac_A_matched, frac_B_matched


# ---------------------------------------------------------------------------
# Library discovery helpers
# ---------------------------------------------------------------------------

def _resolve_requested_frags(vdg_lib_dir, requested):
    """Map user-given `--frags` names onto the library directories that hold them.

    Two ways a name that looks right has no directory: it is a protonation
    variant that was collapsed onto a neutral representative at build time
    (`fragment_aliases.tsv`), and a plain typo. Resolve the first, and fail on
    what is left -- silently comparing fewer fragments than were asked for is
    indistinguishable from a fragment with no vdGs.
    """
    resolved, missing = [], []
    for name in requested:
        if os.path.isdir(os.path.join(vdg_lib_dir, name, 'nr_vdgs')):
            resolved.append(name)
            continue
        alias = utils.smiles_to_filename(
            resolve_fragment_alias(vdg_lib_dir, utils.filename_to_smiles(name)))
        if alias != name and os.path.isdir(os.path.join(vdg_lib_dir, alias, 'nr_vdgs')):
            print(f'[INFO] {name!r} was collapsed onto {alias!r} at build time; '
                  f'using that directory.')
            resolved.append(alias)
        else:
            missing.append(name)
    if missing:
        raise SystemExit(
            f'ERROR: no fragment directory in {vdg_lib_dir} for {missing} (checked '
            f'fragment_aliases.tsv for collapsed protonation variants).')
    return resolved


def _find_completed_frags(vdg_lib_dir, frag_filter):
    """Return fragment dir names with nr_vdgs/ present AND a completed job log.

    A directory can have nr_vdgs/ populated from a run that later crashed
    (e.g. mid-clustering); only the '<frag>_log' completion line is authoritative
    (Frags.check_vdg_job_status). Incomplete fragments are dropped, with a warning.
    """
    candidates = (_resolve_requested_frags(vdg_lib_dir, frag_filter) if frag_filter
                  else sorted(os.listdir(vdg_lib_dir)))
    completed, incomplete = [], []
    for d in candidates:
        if not os.path.isdir(os.path.join(vdg_lib_dir, d, 'nr_vdgs')):
            continue
        if Frags.check_vdg_job_status(d, vdg_lib_dir):
            completed.append(d)
        else:
            incomplete.append(d)
    if incomplete:
        print(f'[WARNING] skipping {len(incomplete)} fragment(s) with incomplete '
              f'vdG-generation jobs (no "Job completed." in <frag>_log): {incomplete}')
    return completed


def _bucket_counts(vdg_lib_dir, frag, subset_size):
    """
    bucket -> cluster count for every npz in one fragment/subset-size dir.

    Deliberately not common.load_bucket_counts: that one drops X and, under
    bb_mode='off', every bb-containing bucket (it serves the propensity path).
    Here bb buckets are ordinary geometry rows -- and
    the size-1 `bb` bucket is usually the library's largest.
    """
    d = os.path.join(vdg_lib_dir, frag, 'nr_vdgs', str(subset_size))
    counts = {}
    if not os.path.isdir(d):
        return counts
    for fname in os.listdir(d):
        if not fname.endswith('.npz'):
            continue
        with np.load(os.path.join(d, fname)) as npz:
            counts[fname[:-4]] = len(npz['cluster_id'])
    return counts


# Per-fragment data (enrichments, bucket counts) is reused across every pair the
# fragment appears in; without these caches a full sweep re-reads it O(n^2) times.
def _cached_counts(cache, vdg_lib_dir, frag, subset_size):
    key = (frag, subset_size)
    if key not in cache:
        cache[key] = _bucket_counts(vdg_lib_dir, frag, subset_size)
    return cache[key]


def _cached_enrichments(cache, vdg_lib_dir, frag):
    if frag not in cache:
        cache[frag] = compute_single_aa_enrichments(
            os.path.join(vdg_lib_dir, frag, 'nr_vdgs'))
    return cache[frag]


def _shared_aa_buckets(vdg_lib_dir, frag_A, frag_B, subset_sizes):
    """Return set of (subset_size, aa_bucket) tuples present in both libraries."""
    def _buckets(frag, size):
        d = os.path.join(vdg_lib_dir, frag, 'nr_vdgs', str(size))
        if not os.path.isdir(d):
            return set()
        return {f[:-4] for f in os.listdir(d) if f.endswith('.npz')}

    shared = set()
    for size in subset_sizes:
        for bucket in _buckets(frag_A, size) & _buckets(frag_B, size):
            shared.add((int(size), bucket))
    return shared


# ---------------------------------------------------------------------------
# Main comparison loop
# ---------------------------------------------------------------------------

def _fmt(x):
    """Format a float for TSV output, or 'nan' for NaN."""
    return 'nan' if math.isnan(x) else f'{x:.4f}'


def compare_pair(vdg_lib_dir, frag_A, frag_B, subset_sizes, threshold,
                 max_per_lib, skip_geometry, enrich_cache=None, counts_cache=None):
    """
    Compare two fragment libraries.  Returns list of row dicts for TSV output.

    enrich_cache / counts_cache are optional cross-pair caches (see
    _cached_enrichments / _cached_counts); omit them to run standalone.
    """
    if enrich_cache is None:
        enrich_cache = {}
    if counts_cache is None:
        counts_cache = {}

    # AA enrichments (single-AA, subset_size=1 only)
    enrich_A = _cached_enrichments(enrich_cache, vdg_lib_dir, frag_A)
    enrich_B = _cached_enrichments(enrich_cache, vdg_lib_dir, frag_B)

    rows = []
    shared = _shared_aa_buckets(vdg_lib_dir, frag_A, frag_B, subset_sizes)

    for subset_size, aa_bucket in sorted(shared):
        cnt_A = _cached_counts(counts_cache, vdg_lib_dir, frag_A, subset_size)
        cnt_B = _cached_counts(counts_cache, vdg_lib_dir, frag_B, subset_size)
        N_A   = cnt_A.get(aa_bucket, 0)
        N_B   = cnt_B.get(aa_bucket, 0)

        # Single-AA enrichment is only defined per residue type, so for
        # multi-residue buckets these columns describe aa_parts[0] alone.
        aa_parts = [p for p in aa_bucket.split('_') if p not in BB_LABELS]
        enrich_a = enrich_A.get(aa_parts[0], float('nan')) if aa_parts else float('nan')
        enrich_b = enrich_B.get(aa_parts[0], float('nan')) if aa_parts else float('nan')

        row = dict(
            frag_A=frag_A, frag_B=frag_B,
            subset_size=subset_size, aa_bucket=aa_bucket,
            N_A=N_A, N_B=N_B,
            enrichment_A=_fmt(enrich_a), enrichment_B=_fmt(enrich_b),
            min_dist='nan', frac_A_matched='nan', frac_B_matched='nan',
        )

        if not skip_geometry and N_A > 0 and N_B > 0:
            bucket_A = load_vdg_bucket(vdg_lib_dir, frag_A, subset_size, aa_bucket)
            bucket_B = load_vdg_bucket(vdg_lib_dir, frag_B, subset_size, aa_bucket)
            if bucket_A is not None and bucket_B is not None:
                try:
                    # Interchangeable slots come from aa_bucket_parts, the array
                    # contractually aligned to the slot axis of 'bb' that perm
                    # indexes -- not from an nr vdG's resnames (which vary
                    # within bb-containing buckets).
                    aa_bucket_parts = bucket_A['aa_bucket_parts']
                    if list(bucket_B['aa_bucket_parts']) != list(aa_bucket_parts):
                        raise ValueError(
                            f'aa_bucket_parts disagree between libraries: '
                            f'{list(aa_bucket_parts)} vs {list(bucket_B["aa_bucket_parts"])}')
                    min_d, frac_a, frac_b = compare_bucket_geometry(
                        bucket_A['cg'], bucket_A['bb'],
                        bucket_B['cg'], bucket_B['bb'],
                        aa_bucket_parts, threshold, max_per_lib)
                    row['min_dist']        = _fmt(min_d)
                    row['frac_A_matched']  = _fmt(frac_a)
                    row['frac_B_matched']  = _fmt(frac_b)
                except Exception as e:
                    print(f'    [WARNING] geometry comparison failed for '
                          f'{aa_bucket}: {e}')

        rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_min_dist_heatmap(rows, frag_A, frag_B, outpath):
    """
    Plot min_dist as a heatmap: rows = aa_bucket, cols = subset_size.
    Only buckets with a finite min_dist are shown.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('[WARNING] matplotlib not available; skipping heatmap.')
        return

    finite_rows = [r for r in rows if r['min_dist'] != 'nan']
    if not finite_rows:
        return

    # Group by subset_size; sort buckets by min_dist ascending
    sizes   = sorted({r['subset_size'] for r in finite_rows})
    buckets = sorted({r['aa_bucket'] for r in finite_rows},
                     key=lambda b: min(float(r['min_dist'])
                                       for r in finite_rows if r['aa_bucket'] == b))

    matrix = np.full((len(buckets), len(sizes)), np.nan)
    buck_idx = {b: i for i, b in enumerate(buckets)}
    size_idx = {s: j for j, s in enumerate(sizes)}
    for r in finite_rows:
        matrix[buck_idx[r['aa_bucket']], size_idx[r['subset_size']]] = float(r['min_dist'])

    fig, ax = plt.subplots(figsize=(max(2, len(sizes)), max(3, len(buckets) * 0.35)))
    cmap = matplotlib.colormaps['RdYlGn_r'].copy()
    cmap.set_bad(color='lightgrey')
    im = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=3.0, aspect='auto',
                   interpolation='nearest')
    ax.set_xticks(range(len(sizes)))
    ax.set_xticklabels([f'size {s}' for s in sizes], fontsize=7)
    ax.set_yticks(range(len(buckets)))
    ax.set_yticklabels(buckets, fontsize=6)
    ax.set_title(f'CG centroid distance (Å)\n{frag_A}  vs  {frag_B}', fontsize=8, pad=8)
    plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02).ax.tick_params(labelsize=6)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

TSV_FIELDS = ['frag_A', 'frag_B', 'subset_size', 'aa_bucket',
              'N_A', 'N_B', 'min_dist', 'frac_A_matched', 'frac_B_matched',
              'enrichment_A', 'enrichment_B']


def _safe_label(s):
    """Make a filesystem-safe filename label; prefix uppercase with '_' so
    case-insensitive filesystems (macOS HFS+) don't collide.

    Input is a library directory name, i.e. already `utils.smiles_to_filename`
    output, so '/' and '\\' cannot reach here. Brackets and parens are kept:
    they are legal in filenames, and dropping them made the map non-injective
    ('CC(=O)[O-]' and 'CC=OO-' collided onto one PNG)."""
    out = []
    for ch in s:
        if ch.isupper():
            out.append('_' + ch)
        elif ch in r'/\:*?"<>|':
            out.append('_')
        else:
            out.append(ch)
    return ''.join(out)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vdg-lib-dir', required=True,
                        help='Path to the fragment library root (frag_lib/).')
    parser.add_argument('--frags', nargs='+', default=None,
                        help='Fragment SMILES (library dir names) to compare. '
                             'Default: all fragments found in vdg-lib-dir. '
                             'If 2 are given, compare only that pair.')
    parser.add_argument('--outdir', default=os.path.join('outputs', 'bioisosteres',
                                                          'cg_geometry'),
                        help='Output directory for TSV and PNG files.')
    parser.add_argument('--match-threshold', type=float, default=1.5,
                        dest='threshold',
                        help='Å cutoff to count an nr vdG as "matched" '
                             '(default: 1.5).')
    parser.add_argument('--subset-sizes', nargs='+', type=int, default=[1, 2],
                        dest='subset_sizes', choices=[1, 2],
                        help='Subset sizes to include (default: 1 2).')
    parser.add_argument('--max-per-lib', type=int, default=1000,
                        dest='max_per_lib',
                        help='Maximum nr vdGs per library per bucket for '
                             'geometry comparison (default: 1000). Larger '
                             'libraries are randomly subsampled.')
    parser.add_argument('--skip-geometry', action='store_true',
                        help='Skip geometric comparison; output AA enrichment only.')
    parser.add_argument('--no-plot', action='store_true',
                        help='Skip heatmap generation.')
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isdir(args.vdg_lib_dir):
        sys.exit(f'[ERROR] vdg-lib-dir not found: {args.vdg_lib_dir}')

    frags = _find_completed_frags(args.vdg_lib_dir, args.frags)
    if len(frags) < 2:
        sys.exit(f'[ERROR] Need at least 2 completed fragment libraries; '
                 f'found: {frags}')

    pairs = list(itertools.combinations(frags, 2))
    print(f'Comparing {len(pairs)} fragment pair(s) across '
          f'subset sizes {args.subset_sizes}.')

    os.makedirs(args.outdir, exist_ok=True)
    tsv_path = os.path.join(args.outdir, 'cg_geometry_comparison.tsv')

    enrich_cache, counts_cache = {}, {}

    with open(tsv_path, 'w') as tsv:
        tsv.write('\t'.join(TSV_FIELDS) + '\n')

        for frag_A, frag_B in pairs:
            print(f'\n{frag_A}  vs  {frag_B}')
            rows = compare_pair(
                args.vdg_lib_dir, frag_A, frag_B,
                args.subset_sizes, args.threshold, args.max_per_lib,
                args.skip_geometry, enrich_cache, counts_cache)

            for row in rows:
                tsv.write('\t'.join(str(row[f]) for f in TSV_FIELDS) + '\n')
            tsv.flush()

            if rows:
                print(f'  {len(rows)} AA buckets compared.')
                finite = [r for r in rows if r['min_dist'] != 'nan']
                if finite:
                    best = min(finite, key=lambda r: float(r['min_dist']))
                    print(f'  Best match: {best["aa_bucket"]} '
                          f'(size {best["subset_size"]}) '
                          f'min_dist={best["min_dist"]} Å')

            if not args.no_plot:
                png = os.path.join(
                    args.outdir,
                    f'min_dist_{_safe_label(frag_A)}_vs_{_safe_label(frag_B)}.png')
                plot_min_dist_heatmap(rows, frag_A, frag_B, png)

    print(f'\nResults written to: {tsv_path}')


if __name__ == '__main__':
    main()
