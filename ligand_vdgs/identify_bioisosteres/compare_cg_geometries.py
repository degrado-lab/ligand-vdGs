"""
Compare fragment libraries for bioisostere identification (AA enrichment + CG geometry),
unifying compute_aa_profiles.py + compare_aa_profiles.py. See
docs/bioisostere_identification.md ("Independent branch") for the mmd2 formula, output
TSV reference, and flags.

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
import hashlib
import argparse
import itertools

import numpy as np

from ligand_vdgs.functions.vdg_struct_utils import BB_LABELS
from ligand_vdgs.functions import utils, ligand_structure
from ligand_vdgs.functions.vdg_npz_utils import (resolve_fragment_alias,
    load_vdg_bucket_all_signs, iter_bucket_files)

from ligand_vdgs.identify_bioisosteres.common import calc_single_aa_propensities

# Fields compare_pair needs from each bucket (subset of vdg_npz_utils.BUCKET_ROW_FIELDS).
BUCKET_FIELDS = ('cg', 'bb', 'cluster_id', 'cluster_num_parents', 'resnames',
                 'slot_flags', 'parent_biounit')

def compute_single_aa_enrichments(nr_vdgs_dir):
    """Return per-residue, background-weighted enrichments or {} if absent."""
    try:
        return calc_single_aa_propensities(nr_vdgs_dir, 'bg_weighted',
                                           bb_mode='per-residue')[0]
    except FileNotFoundError:
        return {}

# Idealized N-CA-C: N-CA 1.458 A, CA-C 1.525 A, N-CA-C 111.2 deg, CA at the origin. A
# constant, data-independent target -- never another vdG -- so the fit cannot be biased
# by whichever library it is being compared against.
_ang = np.radians(111.2)
IDEAL_NCAC = np.array([[1.458, 0.0, 0.0], [0.0, 0.0, 0.0],
                       [1.525 * np.cos(_ang), 1.525 * np.sin(_ang), 0.0]], dtype=np.float64)

def _ideal_frame_points(bb_slot0, cg):
    """CG centroids [N,3] after Kabsch-fitting each vdG's own slot-0 N/CA/C (bb_slot0,
    [N,3,3]) onto the fixed IDEAL_NCAC. The CG never enters the fit."""
    R, t, _ = utils.kabsch(bb_slot0, np.broadcast_to(IDEAL_NCAC, (bb_slot0.shape[0], 3, 3)))
    return np.einsum('ni,nij->nj', cg.mean(axis=1), R) + t

def mmd2(pts_A, w_A, pts_B, w_B):
    """Unbiased weighted squared MMD (distance-induced kernel): the U-statistic
    2E|X-Y| - E|X-X'| - E|Y-Y'|, excluding i==j in the self terms. Each library needs
    >= 2 points; returns nan otherwise."""
    def offdiag_mean(pts, w):
        if len(w) < 2:
            return float('nan')
        weights = w[:, None] * w[None, :]
        np.fill_diagonal(weights, 0.0)
        return float((np.linalg.norm(pts[:, None] - pts[None, :], axis=-1) * weights).sum()
                     / weights.sum())
    def cross_mean(pts_A, w_A, pts_B, w_B):
        weights = w_A[:, None] * w_B[None, :]
        return float((np.linalg.norm(pts_A[:, None] - pts_B[None, :], axis=-1) * weights).sum()
                     / weights.sum())
    return (2 * cross_mean(pts_A, w_A, pts_B, w_B)
            - offdiag_mean(pts_A, w_A) - offdiag_mean(pts_B, w_B))

def _self_split_mmd2(pts, w, parent_biounit):
    """mmd2 between two halves of one library, split by parent biounit (a deterministic
    hash, not Python's randomized str hash) so a structure's clusters land in one half.
    nan if fewer than 2 clusters end up in either half."""
    in_half_a = np.array([int(hashlib.md5(p.encode()).hexdigest(), 16) % 2 == 0
                          for p in parent_biounit])
    if in_half_a.sum() < 2 or (~in_half_a).sum() < 2:
        return float('nan')
    return mmd2(pts[in_half_a], w[in_half_a], pts[~in_half_a], w[~in_half_a])

def compare_bucket_geometry(cg_A, bb_A, w_A, parent_A, cg_B, bb_B, w_B, parent_B,
                            max_per_lib=1000):
    """
    Compare CG centroid geometry between two libraries for one AA bucket via mmd2.

    cg_A/B [N,n_cg,3], bb_A/B [N,num_vdms,3,3] (only slot 0 anchors the fit -- see
    module docstring), w_A/B = cluster_num_parents, parent_A/B = nr_parent_biounit.
    Each library is subsampled to at most max_per_lib nr vdGs first (a cost cap; unlike
    the old per-pair Kabsch fit, subsampling here cannot bias the result).

    Returns dict(mmd2=..., mmd2_self_A=..., mmd2_self_B=...).
    """
    rng = np.random.default_rng(0)
    def subsample(cg, bb, w, parent):
        if len(w) <= max_per_lib:
            return cg, bb, w, parent
        idx = rng.choice(len(w), max_per_lib, replace=False)
        return cg[idx], bb[idx], w[idx], parent[idx]

    cg_A, bb_A, w_A, parent_A = subsample(cg_A, bb_A, w_A, parent_A)
    cg_B, bb_B, w_B, parent_B = subsample(cg_B, bb_B, w_B, parent_B)

    pts_A = _ideal_frame_points(bb_A[:, 0].astype(np.float64), cg_A.astype(np.float64))
    pts_B = _ideal_frame_points(bb_B[:, 0].astype(np.float64), cg_B.astype(np.float64))
    w_A, w_B = w_A.astype(np.float64), w_B.astype(np.float64)

    return dict(mmd2=mmd2(pts_A, w_A, pts_B, w_B),
                mmd2_self_A=_self_split_mmd2(pts_A, w_A, parent_A),
                mmd2_self_B=_self_split_mmd2(pts_B, w_B, parent_B))

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
            f'[ERROR] no fragment directory in {vdg_lib_dir} for {missing} (checked '
            f'fragment_aliases.tsv for collapsed protonation variants).')
    return resolved

def _find_completed_frags(vdg_lib_dir, frag_filter):
    """Return fragment dir names with nr_vdgs/ present AND a completed job log.

    A directory can have nr_vdgs/ populated from a run that later crashed
    (e.g. mid-clustering); only the '<frag>_log' completion line is authoritative
    (ligand_structure.check_vdg_job_status). Incomplete fragments are dropped, with a warning.
    """
    candidates = (_resolve_requested_frags(vdg_lib_dir, frag_filter)
                  if frag_filter else sorted(os.listdir(vdg_lib_dir)))
    completed, incomplete = [], []
    for d in candidates:
        if not os.path.isdir(os.path.join(vdg_lib_dir, d, 'nr_vdgs')):
            continue
        if ligand_structure.check_vdg_job_status(d, vdg_lib_dir):
            completed.append(d)
        else:
            incomplete.append(d)
    if incomplete:
        print(f'[WARNING] skipping {len(incomplete)} fragment(s) with incomplete '
              f'vdG-generation jobs (no "Job completed." in <frag>_log): {incomplete}')
    return completed

def _bucket_counts(vdg_lib_dir, frag, subset_size):
    """bucket -> cluster count for every npz in one fragment/subset-size dir.

    Deliberately not common.load_bucket_counts: that one drops X and, under bb_mode='off',
    every bb-containing bucket (it serves the propensity path). Here bb buckets are ordinary
    geometry rows -- and the size-1 `bb` bucket is usually the library's largest.
    """
    counts = {}
    for _subset, _sign, bucket, path in iter_bucket_files(vdg_lib_dir, frag, subset_size):
        with np.load(path) as npz:
            counts[bucket] = counts.get(bucket, 0) + len(npz['cluster_id'])
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

def _fmt(x):
    """Format a float for TSV output, or 'nan' for NaN."""
    return 'nan' if math.isnan(x) else f'{x:.4f}'

def compare_pair(vdg_lib_dir, frag_A, frag_B, subset_sizes,
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
    shared = {(size, bucket)
              for size in subset_sizes
              for bucket in (_cached_counts(counts_cache, vdg_lib_dir, frag_A, size).keys()
                            & _cached_counts(counts_cache, vdg_lib_dir, frag_B, size).keys())}

    for subset_size, aa_bucket in sorted(shared):
        N_A = _cached_counts(counts_cache, vdg_lib_dir, frag_A, subset_size).get(aa_bucket, 0)
        N_B = _cached_counts(counts_cache, vdg_lib_dir, frag_B, subset_size).get(aa_bucket, 0)

        # Single-AA enrichment is only defined per residue type, so for
        # multi-residue buckets these columns describe aa_parts[0] alone.
        aa_parts = [p for p in aa_bucket.split('_') if p not in BB_LABELS]
        aa0 = aa_parts[0] if aa_parts else None

        row = dict(frag_A=frag_A, frag_B=frag_B, subset_size=subset_size, aa_bucket=aa_bucket,
                   N_A=N_A, N_B=N_B,
                   enrichment_A=_fmt(enrich_A.get(aa0, float('nan')) if aa0 else float('nan')),
                   enrichment_B=_fmt(enrich_B.get(aa0, float('nan')) if aa0 else float('nan')),
                   mmd2='nan', mmd2_self_A='nan', mmd2_self_B='nan')

        if not skip_geometry and N_A > 0 and N_B > 0:
            bucket_A = load_vdg_bucket_all_signs(vdg_lib_dir, frag_A, subset_size, aa_bucket, BUCKET_FIELDS)
            bucket_B = load_vdg_bucket_all_signs(vdg_lib_dir, frag_B, subset_size, aa_bucket, BUCKET_FIELDS)
            if bucket_A is not None and bucket_B is not None:
                try:
                    # aa_bucket_parts is the array contractually aligned to bb's slot
                    # axis (unlike an nr vdG's resnames, which vary within bb-containing
                    # buckets); a mismatch means the same bucket label means something
                    # different in each library.
                    if list(bucket_B['aa_bucket_parts']) != list(bucket_A['aa_bucket_parts']):
                        raise ValueError(
                            f'aa_bucket_parts disagree between libraries: '
                            f'{list(bucket_A["aa_bucket_parts"])} vs '
                            f'{list(bucket_B["aa_bucket_parts"])}')
                    geom = compare_bucket_geometry(
                        bucket_A['cg'], bucket_A['bb'], bucket_A['cluster_num_parents'],
                        bucket_A['parent_biounit'],
                        bucket_B['cg'], bucket_B['bb'], bucket_B['cluster_num_parents'],
                        bucket_B['parent_biounit'], max_per_lib)
                    row['mmd2']        = _fmt(geom['mmd2'])
                    row['mmd2_self_A'] = _fmt(geom['mmd2_self_A'])
                    row['mmd2_self_B'] = _fmt(geom['mmd2_self_B'])
                except Exception as e:
                    print(f'    [WARNING] geometry comparison failed for '
                          f'{aa_bucket}: {e}')

        rows.append(row)

    return rows

def plot_mmd2_heatmap(rows, frag_A, frag_B, outpath):
    """
    Plot mmd2 as a heatmap: rows = aa_bucket, cols = subset_size.
    Only buckets with a finite mmd2 are shown. Diverging colormap centered at 0, since
    mmd2 is an unbiased estimator and can be slightly negative near identity.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('[WARNING] matplotlib not available; skipping heatmap.')
        return

    finite_rows = [r for r in rows if r['mmd2'] != 'nan']
    if not finite_rows:
        return

    # Group by subset_size; sort buckets by mmd2 ascending
    sizes   = sorted({r['subset_size'] for r in finite_rows})
    buckets = sorted({r['aa_bucket'] for r in finite_rows},
                     key=lambda b: min(float(r['mmd2'])
                                       for r in finite_rows if r['aa_bucket'] == b))

    matrix = np.full((len(buckets), len(sizes)), np.nan)
    buck_idx = {b: i for i, b in enumerate(buckets)}
    size_idx = {s: j for j, s in enumerate(sizes)}
    for r in finite_rows:
        matrix[buck_idx[r['aa_bucket']], size_idx[r['subset_size']]] = float(r['mmd2'])

    fig, ax = plt.subplots(figsize=(max(2, len(sizes)), max(3, len(buckets) * 0.35)))
    cmap = matplotlib.colormaps['RdYlGn_r'].copy()
    cmap.set_bad(color='lightgrey')
    bound = np.nanmax(np.abs(matrix)) if np.isfinite(matrix).any() else 1.0
    ax.set_xticks(range(len(sizes)))
    ax.set_xticklabels([f'size {s}' for s in sizes], fontsize=7)
    ax.set_yticks(range(len(buckets)))
    ax.set_yticklabels(buckets, fontsize=6)
    ax.set_title(f'mmd2 (weighted, unbiased)\n{frag_A}  vs  {frag_B}', fontsize=8, pad=8)
    plt.colorbar(ax.imshow(matrix, cmap=cmap, vmin=-bound, vmax=bound, aspect='auto',
                          interpolation='nearest'),
                ax=ax, fraction=0.05, pad=0.02).ax.tick_params(labelsize=6)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()

TSV_FIELDS = ['frag_A', 'frag_B', 'subset_size', 'aa_bucket',
              'N_A', 'N_B', 'mmd2', 'mmd2_self_A', 'mmd2_self_B',
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
                args.subset_sizes, args.max_per_lib,
                args.skip_geometry, enrich_cache, counts_cache)

            for row in rows:
                tsv.write('\t'.join(str(row[f]) for f in TSV_FIELDS) + '\n')
            tsv.flush()

            if rows:
                print(f'  {len(rows)} AA buckets compared.')
                finite = [r for r in rows if r['mmd2'] != 'nan']
                if finite:
                    best = min(finite, key=lambda r: float(r['mmd2']))
                    print(f'  Closest bucket: {best["aa_bucket"]} '
                          f'(size {best["subset_size"]}) '
                          f'mmd2={best["mmd2"]} (self_A={best["mmd2_self_A"]})')

            if not args.no_plot:
                png = os.path.join(
                    args.outdir,
                    f'mmd2_{_safe_label(frag_A)}_vs_{_safe_label(frag_B)}.png')
                plot_mmd2_heatmap(rows, frag_A, frag_B, png)

    print(f'\nResults written to: {tsv_path}')

if __name__ == '__main__':
    main()
