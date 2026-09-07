'''
Compute AA interaction profiles (log-enrichment per AA / AA-pair, relative to background
frequency) for all completed CGs in a vdG library. See
docs/bioisostere_identification.md for the full methodology (formula, size-normalization
methods, CG filtering) and output reference.

Example:
    python compute_aa_profiles.py --vdglib-dir /path/to/frag_lib --skip-pairs
    python compute_aa_profiles.py --vdglib-dir /path/to/frag_lib --plot-single-aa
'''
import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from ligand_vdgs.functions.Frags import check_vdg_job_status

from ligand_vdgs.identify_bioisosteres.common import (BB_MODES, load_bucket_counts,
    calc_single_aa_propensities, calc_aa_pair_propensities, is_bb_category)

NORM_SUFFIX = {
    'none':        '',
    'bg_weighted': '_norm_AA_size_bg_weighted',
}

# The category set is part of what a profile *is*, so it has to be part of the
# filename: two --bb-mode runs into one --outdir would otherwise overwrite each
# other silently, and the flags exist precisely to invite that comparison.
BB_MODE_TAG = {'per-residue': '_bbPR', 'pooled': '_bbPOOL', 'off': '_bbOFF'}


def profile_suffix(norm_method, bb_mode):
    return NORM_SUFFIX[norm_method] + BB_MODE_TAG[bb_mode]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vdglib-dir', type=str, required=True,
                        help='Path to fragment library root (frag_lib/); all completed CGs '
                             'found here are processed.')
    parser.add_argument('--norm-AAsize-method',
                        choices=['none', 'bg_weighted'],
                        action='append',
                        dest='norm_methods',
                        default=None,
                        metavar='METHOD',
                        help='Size-normalization method(s) to run; repeatable to produce '
                             'outputs for several methods in one call. Default: bg_weighted.')
    parser.add_argument('--skip-pairs', action='store_true',
                        help='Skip AA-pair propensity computation (single-AA .npz still '
                             'saved). Use for a faster run when only single-AA comparison '
                             'is needed.')
    parser.add_argument('--plot-single-aa', action='store_true',
                        help='Plot single-AA enrichment bar charts. Single-AA .npz files '
                             'are always saved regardless of this flag.')
    parser.add_argument('--outdir', type=str,
                        default=os.path.join('outputs', 'bioisosteres', 'aa_profiles'),
                        help='Root output directory; each CG gets a subdirectory '
                             '(default: outputs/bioisosteres/aa_profiles/<cg>/).')
    parser.add_argument('--skip-log', type=str, default=None, metavar='FILE',
                        help='Write skipped (incomplete/low-count) CGs to this file '
                             'instead of stdout.')
    parser.add_argument('--bb-mode', choices=list(BB_MODES), default='per-residue',
                        help="How backbone-mediated contacts enter the profile. "
                             "'per-residue' (default) resolves each backbone slot to "
                             "the residue that donated it (GLY-bb, ALA-bb, ...); "
                             "'pooled' keeps one 'bb' category; 'off' drops backbone "
                             "buckets entirely (required with "
                             "--norm-AAsize-method none). Backbone categories are "
                             "shown hatched.")
    parser.add_argument('--pair-bb-mode', choices=list(BB_MODES), default='pooled',
                        help="Same, for the AA-pair heatmap. Defaults to 'pooled' "
                             "because per-residue would make a 39x39 grid whose "
                             "rarer cells hold single-digit counts.")
    parser.add_argument('--min-clusters', type=int, default=50, metavar='N',
                        help='Minimum total single-AA cluster count required to process a '
                             'CG (default: 50); below this, statistics are too unreliable.')
    return parser.parse_args()


def build_heatmap_matrix(propensities):
    '''
    Build a symmetric NxN enrichment matrix from AA-pair propensities. AAs are sorted by
    summed pair enrichment (descending), so the highest-scoring AA appears first
    (row/column 0). Missing pairs are NaN (rendered white in the heatmap).

    Returns (matrix, aa_list); save aa_list alongside matrix so two CGs can be aligned
    by label for downstream comparison.
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


def plot_heatmap(propensities, cg, norm_method, outdir, counts=None, extras=None,
                 bb_mode='pooled'):
    matrix, aas = build_heatmap_matrix(propensities)

    cmap = matplotlib.colormaps['RdBu_r'].copy()
    cmap.set_bad(color='white')  # NaN (missing pairs) → white

    # Scale with the category count: per-residue bb makes this 39 wide, and a
    # fixed 4in canvas overlaps the tick labels well before that.
    side = max(4.0, 0.22 * len(aas) + 1.5)
    fig, ax = plt.subplots(figsize=(side, side * 0.78))
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
    # Backbone categories are a different kind of claim from a sidechain one
    # ("contacts via backbone", any residue), so mark their ticks rather than
    # letting them read as a 21st amino acid.
    for axis in (ax.get_xticklabels(), ax.get_yticklabels()):
        for tick, aa in zip(axis, aas):
            if is_bb_category(aa):
                tick.set_color('#8c510a')
                tick.set_fontweight('bold')
    ax.set_title(f'AA Pair Propensities: {cg}', pad=10, fontsize=6)

    suffix = profile_suffix(norm_method, bb_mode)
    base = os.path.join(outdir, f'{cg}_aa_pair_freq{suffix}')
    plt.savefig(base + '.png', dpi=400, transparent=True, bbox_inches='tight')
    plt.close()

    # Save matrix + labels for downstream bioisostere comparison, with the same
    # provenance the single-AA npz carries: pair values moved twice (additive-size
    # null -> product null, and bb pair buckets from always-excluded to included),
    # so a stale file aligns on labels against a fresh one and correlates silently.
    # Load with: d = np.load(...npz); matrix, aas = d['matrix'], list(d['aa_labels'])
    extras = extras or {}
    np.savez(base + '.npz', matrix=matrix, aa_labels=np.array(aas, dtype=str),
             bb_mode=np.array(bb_mode),
             total_count=np.array(sum((counts or {}).values()), dtype=int),
             excluded_noncanonical=np.array(extras.get('X', 0), dtype=int),
             excluded_noncanonical_bb=np.array(extras.get('noncanonical_bb', 0), dtype=int))


def plot_bar_chart(propensities, cg, norm_method, outdir, counts=None, extras=None,
                   bb_mode='per-residue'):
    aas = list(propensities.keys())
    freqs = list(propensities.values())

    non_bb_freqs = [f for aa, f in zip(aas, freqs) if not is_bb_category(aa)]
    vmin = min(non_bb_freqs) - 0.2 if non_bb_freqs else -0.2
    vmax = max(non_bb_freqs) if non_bb_freqs else 1.0
    norm = plt.Normalize(vmin, vmax)

    # Per-residue backbone deconvolution makes this up to 39 bars, each with a
    # two-line label; a fixed 4in canvas is unreadable past ~20.
    fig, ax = plt.subplots(figsize=(max(4.0, 0.20 * len(aas)), 2))
    for i, (aa, freq) in enumerate(zip(aas, freqs)):
        if is_bb_category(aa):
            ax.bar(i, freq, color='#cccccc', edgecolor='black', width=0.7,
                   hatch='//', label='backbone')
        else:
            ax.bar(i, freq, color=cm.viridis(norm(freq)), edgecolor='black', width=0.7)

    ax.set_ylabel('Enrichment', fontsize=5)
    ax.tick_params(axis='y', labelsize=4)
    ax.set_xticks(np.arange(len(aas)))
    # One line, not two: the labels are rotated 90 deg, so a second line sits
    # *beside* the first and the two collide once there are 39 of them.
    if counts:
        xlabels = [f'{aa} (N={counts.get(aa, 0)})' for aa in aas]
    else:
        xlabels = aas
    ax.set_xticklabels(xlabels, rotation=90, fontsize=5)
    ax.tick_params(axis='x', length=0)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)

    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)

    # Title: SMILES bold red; N= smaller gray below it
    total = sum(counts.values()) if counts else 0
    ax.set_title(cg, fontsize=8, fontweight='bold', color='darkred', pad=14)
    if total:
        # Report the X clusters that were excluded rather than letting them be
        # invisible: they are real geometry with no defensible AA prior.
        n_x = (extras or {}).get('X', 0)
        sub = f'N={total:,}'
        if n_x:
            sub += f'   (+{n_x:,} X excluded)'
        ax.text(0.5, 1.0, sub, transform=ax.transAxes,
                ha='center', va='bottom', fontsize=5, color='#666666')

    suffix = profile_suffix(norm_method, bb_mode)
    plt.savefig(os.path.join(outdir, f'{cg}_single_aa_freq{suffix}.png'),
                dpi=400, transparent=True, bbox_inches='tight')
    plt.close()


def save_single_aa_npz(propensities, counts, cg, norm_method, outdir, extras=None,
                       bb_mode='per-residue'):
    '''Save single-AA enrichment profile as .npz for downstream comparison.

    Arrays saved:
      aa_labels   — contact-category names, in the same order as enrichments/counts
      enrichments — log-enrichment value per category
      counts      — vdG cluster count per category (for confidence weighting)
      total_count — sum of all counted single-AA clusters (scalar)
      bb_mode     — how backbone contacts were categorized; profiles built under
                    different modes have different label sets and must not be
                    compared cell-for-cell
      excluded_*  — clusters left out of both numerator and denominator
    '''
    aas = list(propensities.keys())
    suffix = profile_suffix(norm_method, bb_mode)
    extras = extras or {}
    np.savez(os.path.join(outdir, f'{cg}_single_aa_freq{suffix}.npz'),
             aa_labels=np.array(aas, dtype=str),
             enrichments=np.array([propensities[aa] for aa in aas]),
             counts=np.array([counts.get(aa, 0) for aa in aas], dtype=int),
             total_count=np.array(sum(counts.values()), dtype=int),
             bb_mode=np.array(bb_mode),
             excluded_noncanonical=np.array(extras.get('X', 0), dtype=int),
             excluded_noncanonical_bb=np.array(extras.get('noncanonical_bb', 0), dtype=int))


def _write_log(text, log_path, label):
    if log_path:
        os.makedirs(os.path.dirname(log_path) or '.', exist_ok=True)
        open_mode = 'a' if os.path.exists(log_path) else 'w'
        with open(log_path, open_mode) as f:
            f.write(text + '\n')
        print(f'{label} written to: {log_path}')
    else:
        print(text)


def _report_list(items, header, skip_log, label):
    if not items:
        return
    _write_log('\n'.join([header] + [f'  {i}' for i in items]), skip_log, label)


def report_skipped(skipped, skip_log):
    _report_list([f'{s}' for s in skipped], f'Skipped {len(skipped)} incomplete CG(s):',
                 skip_log, 'Skipped CG list')


def main():
    args = parse_args()
    norm_methods = args.norm_methods if args.norm_methods else ['bg_weighted']

    # 'none' has no size term, so a 4-atom backbone category has no defensible
    # prior against a whole sidechain. Fail loudly instead of inventing one.
    if 'none' in norm_methods and (args.bb_mode != 'off' or args.pair_bb_mode != 'off'):
        print("[ERROR] --norm-AAsize-method none requires --bb-mode off and "
              "--pair-bb-mode off (see common.category_priors).")
        return

    if not os.path.isdir(args.vdglib_dir):
        print(f'[ERROR] vdglib-dir not found: {args.vdglib_dir}')
        return

    plot_single_aa_bars = args.plot_single_aa
    run_pairs = not args.skip_pairs

    print(f'vdglib-dir:   {args.vdglib_dir}')
    print(f'norm methods: {norm_methods}')
    print(f'bb mode:      {args.bb_mode} (single) / {args.pair_bb_mode} (pairs)')
    print(f'min-clusters: {args.min_clusters}')
    print(f'outdir root:  {args.outdir}')
    print()

    skipped = []
    skipped_low_count = []
    processed = []

    for dirname in sorted(os.listdir(args.vdglib_dir)):
        if not os.path.isdir(os.path.join(args.vdglib_dir, dirname)):
            continue
        if not check_vdg_job_status(dirname, args.vdglib_dir):
            skipped.append(dirname)
            continue

        cg = dirname
        nr_vdgs_dir = os.path.join(args.vdglib_dir, cg, 'nr_vdgs')

        # Preload single-AA bucket counts once per CG (reused across all norm_methods)
        single_dir = os.path.join(nr_vdgs_dir, '1')
        if not os.path.isdir(single_dir):
            skipped_low_count.append((cg, 0))
            continue
        preloaded_single, single_extras = load_bucket_counts(single_dir,
                                                             bb_mode=args.bb_mode)
        # Keep the filtering denominator consistent with the profile and NPZ
        # output: whatever --bb-mode admits is part of the total cluster count.
        total_count = sum(preloaded_single.values())
        if total_count < args.min_clusters:
            skipped_low_count.append((cg, total_count))
            continue

        # Preload pair bucket counts once per CG
        preloaded_pair, pair_extras = None, {}
        if run_pairs:
            pair_dir = os.path.join(nr_vdgs_dir, '2')
            if os.path.isdir(pair_dir):
                preloaded_pair, pair_extras = load_bucket_counts(
                    pair_dir, bb_mode=args.pair_bb_mode)
            else:
                print(f'  [WARNING] No AA-pair vdGs found at: {pair_dir}; skipping pair heatmaps for {cg}.')

        outdir = os.path.join(args.outdir, cg)
        os.makedirs(outdir, exist_ok=True)

        print(f'Processing: {cg}  (N={total_count})')
        for norm_method in norm_methods:
            # Single-AA: always compute and save NPZ; optionally plot bar chart.
            try:
                single_props, single_counts = calc_single_aa_propensities(
                    nr_vdgs_dir, norm_method, preloaded_single,
                    bb_mode=args.bb_mode)
            except FileNotFoundError as e:
                print(f'  [WARNING] {e}; skipping single-AA.')
                single_props, single_counts = {}, {}

            if single_props:
                save_single_aa_npz(single_props, single_counts, cg, norm_method, outdir,
                                   single_extras, args.bb_mode)
                if plot_single_aa_bars:
                    plot_bar_chart(single_props, cg, norm_method, outdir, single_counts,
                                   single_extras, args.bb_mode)
            elif plot_single_aa_bars:
                print(f'  [WARNING] No qualifying single-AA vdGs ({norm_method}); skipping bar chart.')

            if run_pairs and preloaded_pair is not None:
                pair_props = calc_aa_pair_propensities(
                    nr_vdgs_dir, norm_method, preloaded_pair,
                    bb_mode=args.pair_bb_mode)
                if not pair_props:
                    print(f'  [WARNING] No qualifying AA-pair vdGs ({norm_method}); skipping heatmap.')
                else:
                    plot_heatmap(pair_props, cg, norm_method, outdir, preloaded_pair,
                                 pair_extras, args.pair_bb_mode)

        processed.append(cg)

    print(f'\nDone. Processed {len(processed)} CG(s).')

    _report_list([f'{cg}  (N={n})' for cg, n in skipped_low_count],
                 f'Skipped {len(skipped_low_count)} low-count CG(s) '
                 f'(< {args.min_clusters} total clusters):',
                 args.skip_log, 'Low-count CG list')

    report_skipped(skipped, args.skip_log)


if __name__ == '__main__':
    main()
