'''
Calculate the likelihood of each AA to interact with a given chemical group (CG).

1. Single-AA propensities: bar chart of log-enrichment per AA.
2. AA-pair propensities: symmetric heatmap of log-enrichment per AA pair.

Propensity = log(observed / expected), where:
  `observed` = number of non-redundant vdG clusters for that AA bucket,
               optionally divided by the AA's heavy-atom count.
  `expected` = total cluster count × background AA frequency,
               optionally divided by an average AA (pair) heavy-atom count.

Size normalization (--norm-AAsize-method, default: bg_weighted):
  bg_weighted  Single-AA: observed / AA size; expected / background-frequency-
               weighted mean size. AA-pair: analogously, expected is divided by
               2 × the frequency-weighted mean single-AA size (= the expected
               combined size of two independently frequency-sampled AAs).
               Enrichment = 0 when an AA's count is fully explained by both
               its background frequency and its size; positive values indicate
               contacts beyond what size and frequency together predict.
               This is the theoretically correct normalization.
  unweighted   Single-AA: observed / AA size; expected / unweighted mean size
               across all 20 AAs. AA-pair: same logic applied to combined pair
               sizes and their unweighted mean over all 210 unique pairs.
               Less principled than bg_weighted; provided for comparison only.
  none         Raw counts vs. background frequency only (no size correction).

For AA-pair buckets, `expected` is multiplied by 2 for heterogeneous pairs
(e.g. ASP_HIS) because the canonical sorted bucket absorbs both orderings
(ASP+HIS and HIS+ASP); homogeneous pairs (ASP_ASP) have only one ordering,
so no factor of 2 applies.

Using cluster count (not raw PDB occurrence) corrects for database bias
toward overrepresented protein families (e.g. kinases).

Backbone-containing buckets (e.g. "bb", "ALA_bb") are excluded: backbone
geometry is conserved across all AAs and does not contribute discriminating
power to sidechain-based propensity analysis.

The heatmap is also saved as a .npz (matrix + AA labels) for downstream
bioisostere comparison across CGs.

Example usage (default: bg_weighted normalization, AA-pair heatmap only):
    python visualize_aa_profiles.py \
        --smiles "CNC1CC1" \
        --vdglib-dir /path/to/frag_lib

    # Can also produce the single-AA bar chart (pass both flags to get both outputs):
    python visualize_aa_profiles.py \
        --smiles "CNC1CC1" \
        --vdglib-dir /path/to/frag_lib \
        --plot-single-aa-freqs --plot-aa-pair-freqs
'''
import os
import sys
import argparse
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.append(os.path.join(os.path.dirname(__file__), '../functions'))
from utils import resolve_cg

sys.path.append(os.path.join(os.path.dirname(__file__), '../tools'))
from reference_values import (AA_size, AA_background_freq,
                               avg_aa_size_unweighted, avg_aa_size_bg_weighted,
                               avg_aa_pair_size_unweighted, avg_aa_pair_size_bg_weighted)

NORM_SUFFIX = {
    'none':        '',
    'unweighted':  '_norm_AA_size',
    'bg_weighted': '_norm_AA_size_bg_weighted',
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--smiles', type=str, required=True,
                        help='CG SMILES string. Resolved to the matching frag_lib '
                             'directory via canonical equivalence.')
    parser.add_argument('--vdglib-dir', type=str, required=True,
                        help='Path to fragment library root (frag_lib/).')
    parser.add_argument('--norm-AAsize-method',
                        choices=['none', 'unweighted', 'bg_weighted'],
                        action='append',
                        dest='norm_methods',
                        default=None,
                        metavar='METHOD',
                        help='Size-normalization method(s) to run. Choices: none, '
                             'unweighted, bg_weighted. Can be repeated to produce '
                             'outputs for multiple methods in one call. '
                             'Default when omitted: bg_weighted.')
    parser.add_argument('--plot-single-aa-freqs', action='store_true',
                        help='Plot bar chart of single-AA interaction propensities. '
                             'Suppresses the AA-pair heatmap unless --plot-aa-pair-freqs '
                             'is also passed.')
    parser.add_argument('--plot-aa-pair-freqs', action='store_true',
                        help='Plot heatmap of AA-pair interaction propensities '
                             '(on by default; only needed when combining with '
                             '--plot-single-aa-freqs).')
    parser.add_argument('--outdir', type=str, default=None,
                        help='Output directory for plots and data '
                             '(default: outputs/aa_profiles/<cg>/).')
    return parser.parse_args()


def load_bucket_counts(npz_dir):
    '''
    Read all .npz files in npz_dir and return a dict mapping bucket name
    (e.g. "LYS" or "ASP_HIS") -> number of non-redundant vdG clusters.

    Cluster count (len of cluster_id array) is used rather than raw PDB
    occurrence (sum of cluster_size) to avoid bias from overrepresented
    protein families in the database.

    Backbone-containing buckets (any part == "bb") are excluded.
    '''
    counts = {}
    for fname in os.listdir(npz_dir):
        if not fname.endswith('.npz'):
            continue
        bucket = fname[:-4]  # strip .npz
        if 'bb' in bucket.split('_'):
            continue
        with np.load(os.path.join(npz_dir, fname), allow_pickle=True) as npz:
            counts[bucket] = len(npz['cluster_id'])
    return counts


def calc_single_aa_propensities(nr_vdgs_dir, norm_method):
    '''Returns dict of AA -> log-enrichment, sorted ascending.'''
    single_dir = os.path.join(nr_vdgs_dir, '1')
    if not os.path.isdir(single_dir):
        raise FileNotFoundError(f'No single-AA vdGs found at: {single_dir}')
    counts = load_bucket_counts(single_dir)
    total = sum(counts.values())  # total single-AA clusters; independent of pair analysis

    avg_size = None
    if norm_method == 'unweighted':
        avg_size = avg_aa_size_unweighted()
    elif norm_method == 'bg_weighted':
        avg_size = avg_aa_size_bg_weighted()
    elif norm_method != 'none':
        raise ValueError(f'Unknown norm_method: {norm_method!r}')

    propensities = {}
    for aa, count in counts.items():
        bg_freq = AA_background_freq.get(aa, 0)
        if bg_freq <= 0:
            continue
        if norm_method == 'none':
            observed = count
            expected = total * bg_freq
        else:
            aa_size = AA_size.get(aa, 0)
            if aa_size == 0:
                continue
            observed = count / aa_size
            expected = total * bg_freq / avg_size
        if observed > 0 and expected > 0:
            propensities[aa] = math.log(observed / expected)

    return dict(sorted(propensities.items(), key=lambda x: x[1]))


def calc_aa_pair_propensities(nr_vdgs_dir, norm_method):
    '''Returns dict of "AA1_AA2" -> log-enrichment (unsorted).'''
    pair_dir = os.path.join(nr_vdgs_dir, '2')
    if not os.path.isdir(pair_dir):
        raise FileNotFoundError(f'No AA-pair vdGs found at: {pair_dir}')
    counts = load_bucket_counts(pair_dir)
    total = sum(counts.values())  # total pair-AA clusters; independent of single-AA analysis

    avg_pair_size = None
    if norm_method == 'unweighted':
        avg_pair_size = avg_aa_pair_size_unweighted()
    elif norm_method == 'bg_weighted':
        avg_pair_size = avg_aa_pair_size_bg_weighted()
    elif norm_method != 'none':
        raise ValueError(f'Unknown norm_method: {norm_method!r}')

    propensities = {}
    for bucket, count in counts.items():
        parts = bucket.split('_')
        if len(parts) != 2:
            continue
        aa1, aa2 = parts
        bg1 = AA_background_freq.get(aa1, 0)
        bg2 = AA_background_freq.get(aa2, 0)
        if bg1 <= 0 or bg2 <= 0:
            continue
        # Buckets are canonical sorted pairs (alphabetical), so a heterogeneous
        # bucket like ASP_HIS absorbs both orderings (ASP+HIS and HIS+ASP).
        # The null-model expected count must account for both orderings too:
        #   expected(ASP_HIS) = total × (bg_ASP×bg_HIS + bg_HIS×bg_ASP)
        #                     = total × 2 × bg_ASP × bg_HIS
        # Homogeneous pairs (ASP_ASP) have only one ordering, so no factor of 2.
        symmetry_factor = 1 if aa1 == aa2 else 2
        if norm_method == 'none':
            observed = count
            expected = total * symmetry_factor * bg1 * bg2
        else:
            pair_size = AA_size.get(aa1, 0) + AA_size.get(aa2, 0)
            if pair_size == 0:
                continue
            observed = count / pair_size
            expected = total * symmetry_factor * bg1 * bg2 / avg_pair_size
        if observed > 0 and expected > 0:
            propensities[bucket] = math.log(observed / expected)

    return propensities


def build_heatmap_matrix(propensities):
    '''
    Build a symmetric NxN enrichment matrix from AA-pair propensities.
    AAs are sorted by summed pair enrichment (descending), so the AA with the
    highest total pair enrichment score appears first (row/column 0).
    Missing pairs are NaN (rendered as white in the heatmap).

    Returns (matrix, aa_list). Save aa_list alongside matrix so two CGs
    can be aligned by label for downstream comparison.
    '''
    aa_scores = {}
    for pair, val in propensities.items():
        aa1, aa2 = pair.split('_')
        aa_scores[aa1] = aa_scores.get(aa1, 0) + val
        if aa2 != aa1:
            aa_scores[aa2] = aa_scores.get(aa2, 0) + val
    aas = sorted(aa_scores, key=aa_scores.get, reverse=True)

    idx = {aa: i for i, aa in enumerate(aas)}
    matrix = np.full((len(aas), len(aas)), np.nan)
    for pair, val in propensities.items():
        aa1, aa2 = pair.split('_')
        i, j = idx[aa1], idx[aa2]
        matrix[i, j] = val
        matrix[j, i] = val  # symmetric
    return matrix, aas


def plot_heatmap(propensities, cg, norm_method, outdir):
    matrix, aas = build_heatmap_matrix(propensities)

    cmap = matplotlib.colormaps['RdBu_r'].copy()
    cmap.set_bad(color='white')  # NaN (missing pairs) → white

    fig, ax = plt.subplots(figsize=(4, 3))
    vmax = np.nanmax(np.abs(matrix))
    img = ax.imshow(matrix, cmap=cmap, vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    cbar = plt.colorbar(img, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=5)
    cbar.set_label('Enrichment', fontsize=7)
    ax.set_xticks(np.arange(len(aas)))
    ax.set_yticks(np.arange(len(aas)))
    ax.set_xticklabels(aas, rotation=90, fontsize=6)
    ax.set_yticklabels(aas, fontsize=6)
    ax.set_title(f'AA Pair Propensities: {cg}', pad=10, fontsize=6)

    suffix = NORM_SUFFIX[norm_method]
    base = os.path.join(outdir, f'{cg}_aa_pair_freq{suffix}')
    plt.savefig(base + '.png', dpi=400, transparent=True, bbox_inches='tight')
    plt.close()

    # Save matrix + labels for downstream bioisostere comparison.
    # Load with: d = np.load(...npz, allow_pickle=True)
    #            matrix, aas = d['matrix'], list(d['aa_labels'])
    np.savez(base + '.npz', matrix=matrix, aa_labels=np.array(aas, dtype=str))


def plot_bar_chart(propensities, cg, norm_method, outdir):
    aas = list(propensities.keys())
    freqs = list(propensities.values())

    norm = plt.Normalize(min(freqs) - 0.2, max(freqs))  # shift min to avoid dark purples
    colors = cm.viridis(norm(freqs))

    fig, ax = plt.subplots(figsize=(4, 1.5))
    ax.bar(np.arange(len(aas)), freqs, color=colors, edgecolor='black', width=0.7)
    ax.set_ylabel('Enrichment', fontsize=5)
    ax.tick_params(axis='y', labelsize=4)
    ax.set_xticks(np.arange(len(aas)))
    ax.set_xticklabels(aas, rotation=90, fontsize=7)
    ax.tick_params(axis='x', length=0)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)

    ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5, label='Null (0)')
    ax.legend(loc='upper left', fontsize=5)

    suffix = NORM_SUFFIX[norm_method]
    plt.savefig(os.path.join(outdir, f'{cg}_single_aa_freq{suffix}.png'),
                dpi=400, transparent=True, bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()
    norm_methods = args.norm_methods if args.norm_methods else ['bg_weighted']

    try:
        cg = resolve_cg(args.vdglib_dir, args.smiles)
    except FileNotFoundError:
        print(f'[ERROR] No fragment library entry matches SMILES "{args.smiles}" in {args.vdglib_dir}.\n'
              f'        Check that the SMILES is correct and that the library is built for this fragment.')
        return

    outdir = args.outdir if args.outdir is not None else os.path.join('outputs', 'aa_profiles', cg)
    os.makedirs(outdir, exist_ok=True)
    nr_vdgs_dir = os.path.join(args.vdglib_dir, cg, 'nr_vdgs')

    print(f'CG (frag_lib dir): {cg}')
    print(f'norm methods:      {norm_methods}')
    print(f'outdir:            {outdir}')

    # Default: run the AA-pair heatmap unless the user explicitly asked only for
    # single-AA plots (in which case they would not pass --plot-aa-pair-freqs).
    run_single = args.plot_single_aa_freqs
    run_pairs  = args.plot_aa_pair_freqs or not args.plot_single_aa_freqs

    for norm_method in norm_methods:
        if run_single:
            propensities = calc_single_aa_propensities(nr_vdgs_dir, norm_method)
            if not propensities:
                print(f'[WARN] No qualifying single-AA vdGs for {cg} ({norm_method}); skipping bar chart.')
            else:
                plot_bar_chart(propensities, cg, norm_method, outdir)

        if run_pairs:
            propensities = calc_aa_pair_propensities(nr_vdgs_dir, norm_method)
            if not propensities:
                print(f'[WARN] No qualifying AA-pair vdGs for {cg} ({norm_method}); skipping heatmap.')
            else:
                plot_heatmap(propensities, cg, norm_method, outdir)


if __name__ == '__main__':
    main()
