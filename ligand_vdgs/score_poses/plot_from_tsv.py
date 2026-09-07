"""
Correlate vdG hit-finder scores against ligand RMSD using vdg_hit_finder.py's
TSV output. Per structure dir under --matches-root: ligand_rmsd_summary.tsv
(pdbfile, ligand_rmsd) and vdg_hits.tsv (pdbfile, frag, subset_size, vdg_rmsd).

One panel per condition (RMSD cutoff x subset-size/weight config) in --out-fig.
Scores are only comparable within a structure, so each title leads with the mean
per-structure Spearman; the pooled value follows, for spread.

--detail-figures adds a per-condition 4-panel figure with external-rank panels.
Needs sample names ending in digits (e.g. "<label>_3"); samples containing
"ground_truth" are starred/outlined, and the figure is skipped if no sample
matches that convention.

Usage:
  python plot_from_tsv.py --matches-root <dir> [--detail-figures] [--per-struct-detail]
"""

import os, re, csv, math, argparse
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

GT_MARK = "ground_truth"
IGNORE_FRAGS = []  # e.g. ['CC(C)N', 'CN(C)C']
RMSD_CUTS = [0.6, 0.7]
SIZE_CONFIGS = [
    ("sz1_only_w1.0",   {1: 1.0}),
    ("sz2_only_w1.0",   {2: 1.0}),
    ("sz12_unweighted", {1: 1.0, 2: 1.0}),
    ("sz12_weighted",   {1: 0.5, 2: 1.0}),
]


# ------------------------------- stats -------------------------------------

def _nan_safe_corr(fn):
    def wrapped(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        if x.size < 2 or np.allclose(x.std(), 0) or np.allclose(y.std(), 0):
            return float("nan")
        return float(fn(x, y).statistic)
    return wrapped

pearson = _nan_safe_corr(stats.pearsonr)
spearman = _nan_safe_corr(stats.spearmanr)


def linreg(x, y):
    """Intercept, slope for y ~ a + b*x; NaN-safe."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.size < 2:
        return float("nan"), float("nan")
    if np.allclose(x.var(), 0):
        return float(y.mean()), float("nan")
    r = stats.linregress(x, y)
    return float(r.intercept), float(r.slope)


def mean_or_nan(vals):
    vals = [v for v in vals if not math.isnan(v)]
    return sum(vals) / len(vals) if vals else float("nan")


# ------------------------------- I/O ----------------------------------------

def sample_name(pdbfile):
    """Strip compound extensions (.pdb.gz etc.) to get the sample id."""
    for ext in (".pdb.gz", ".pdb", ".cif.gz", ".cif"):
        if pdbfile.endswith(ext):
            return pdbfile[: -len(ext)]
    return os.path.splitext(pdbfile)[0]


def read_tsv(path, parse):
    """Rows of `path` mapped through `parse`; unparseable rows skipped."""
    if not os.path.isfile(path):
        print(f"  [WARNING] missing {path}")
        return []
    out = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out.append(parse(row))
            except (KeyError, ValueError):
                continue
    return out


def score_samples(hit_rows, rmsd_cut, weights, normalize, all_samples):
    """
    {sample: sum(weights[subset_size]) over qualifying hits}. If `normalize`,
    each fragment's count is first divided by its mean count across samples (so
    ubiquitous fragments don't dominate), then averaged over fragments.

    `all_samples` is the full roster: zero-hit samples score 0 and must also sit
    in the denominator of every fragment's mean, else every score in the
    condition is scaled by n_all/n_hit -- harmless within a structure, but it
    distorts anything pooled across structures.
    """
    counts = {}  # sample -> frag -> weighted count
    for sample, frag, subset_size, vdg_rmsd in hit_rows:
        w = weights.get(subset_size, 0.0)
        if w <= 0 or vdg_rmsd > rmsd_cut or frag in IGNORE_FRAGS:
            continue
        per_frag = counts.setdefault(sample, {})
        per_frag[frag] = per_frag.get(frag, 0.0) + w

    # Union, not intersect: a sample with hits but no ligand RMSD is still
    # scored; condition_points drops it later for want of an x.
    samples = sorted(set(counts) | set(all_samples))
    get = lambda s, f: counts.get(s, {}).get(f, 0.0)
    if not normalize:
        return {s: sum(counts.get(s, {}).values()) for s in samples}
    frags = sorted({f for fc in counts.values() for f in fc})
    avg = {f: sum(get(s, f) for s in samples) / max(1, len(samples)) for f in frags}
    return {s: sum(get(s, f) / avg[f] for f in frags if avg[f] > 0) / max(1, len(frags))
            for s in samples}


def condition_data(matches_root, labels, rmsd_cut, weights, normalize):
    """{label: (scores, ligand_rmsd_map)} for one condition."""
    data = {}
    for label in labels:
        d = os.path.join(matches_root, label)
        rmsd_map = dict(read_tsv(os.path.join(d, "ligand_rmsd_summary.tsv"),
                                 lambda r: (sample_name(r["pdbfile"]), float(r["ligand_rmsd"]))))
        hits = read_tsv(os.path.join(d, "vdg_hits.tsv"),
                        lambda r: (sample_name(r["pdbfile"]), r["frag"],
                                   int(r["subset_size"]), float(r["vdg_rmsd"])))
        data[label] = (score_samples(hits, rmsd_cut, weights, normalize, rmsd_map), rmsd_map)
    return data


def condition_points(data, lig_rmsd_max=5.0):
    """-> pooled [(rmsd, score, label)] and per-structure correlations.

    A raw score is a weighted hit count over one structure's ligand and pocket
    (and, when normalized, over its own per-fragment mean), so the pooled y axis
    mixes scales and its correlation partly measures which structures sit high
    or low. Putting every structure on a common scale first (per-structure ranks
    or z-scores, in condition_data's output) would fix it; until then the mean
    per-structure Spearman is the defined statistic and the pooled one is spread.
    """
    pooled, per_struct = [], []
    for label, (scores, rmsd_map) in data.items():
        pts = [(rmsd_map[s], sc) for s, sc in scores.items()
               if rmsd_map.get(s, math.inf) < lig_rmsd_max]
        if not pts:
            continue
        xs, ys = zip(*pts)
        pooled += [(x, y, label) for x, y in pts]
        per_struct.append(dict(label=label, pearson=pearson(xs, ys),
                               spearman=spearman(xs, ys), n_pairs=len(pts)))
    return pooled, per_struct


def summarize(cfg_name, weights, per_struct, detail=False):
    print(f"\n~~~ Condition: {cfg_name} (weights: {weights}) ~~~")
    if not per_struct:
        print("  [INFO] No usable structures for this condition.")
        return
    print(f"  structures used: {len(per_struct)}")
    for key in ("spearman", "pearson"):
        m = mean_or_nan([s[key] for s in per_struct])
        print(f"  mean {key.capitalize():8s}(ligRMSD, score): "
              + ("    NA" if math.isnan(m) else f"{m:6.3f}"))
    if detail:
        for m in per_struct:
            print(f"    {m['label']}: n={m['n_pairs']:3d}, "
                  f"Spearman={m['spearman']:6.3f}, Pearson={m['pearson']:6.3f}")


# ------------------------------- plotting -----------------------------------

def label_colors(labels):
    cmap = plt.colormaps.get_cmap("nipy_spectral")
    return {lbl: cmap(0.05 + 0.90 * i / max(1, len(labels) - 1))
            for i, lbl in enumerate(labels)}


def despine(*axes):
    for ax in axes:
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)


def plot_conditions(cond_results, out_fig):
    """One scatter panel per condition; cond_results[cfg] holds points + stats."""
    if not cond_results:
        print("[INFO] No conditions to plot.")
        return

    all_labels = sorted({p[2] for r in cond_results.values() for p in r["points"]})
    color = label_colors(all_labels)
    ncols = min(3, len(cond_results))
    nrows = (len(cond_results) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows), squeeze=False)
    flat = [ax for row in axes for ax in row]

    for ax, (cfg, res) in zip(flat, cond_results.items()):
        pts = res["points"]
        if not pts:
            ax.set_title(cfg + " (no points)"); ax.axis("off"); continue

        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        for lbl in all_labels:
            sub = [(p[0], p[1]) for p in pts if p[2] == lbl]
            if sub:
                ax.scatter(*zip(*sub), s=12, alpha=0.75, edgecolors="none",
                           label=lbl, color=color[lbl])
        a, b = linreg(xs, ys)
        if not (math.isnan(a) or math.isnan(b)):
            xl = [min(xs), max(xs)]
            ax.plot(xl, [a + b * x for x in xl], "--", lw=2.0, color="black")

        ax.set_xlabel("Ligand RMSD to crystal (Å)")
        ax.set_ylabel("vdG score (weighted hit count)")
        ax.set_xlim(0.0, max(5.0, max(xs) * 1.05))
        ax.set_title(f"{cfg}\nper-struct ρ={res['mean_spearman']:5.2f} "
                     f"({res['n_structs']} structs)  |  pooled ρ={res['spearman']:5.2f}, "
                     f"r={res['pearson']:5.2f}, n={len(xs)}", fontsize=9)
        despine(ax)

    for ax in flat[len(cond_results):]:
        ax.axis("off")

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_fig, dpi=300)
    print(f"\n[INFO] Saved figure: {out_fig}")


def _cluster_by_y(named_y, tol):
    """[(name, y)] -> [(y_center, [names])], grouping points within tol."""
    clusters, cur, lo, hi = [], [], None, None
    for name, y in sorted(named_y, key=lambda t: t[1]):
        if cur and y - hi > tol:
            clusters.append(((lo + hi) / 2.0, cur)); cur, lo = [], y
        cur.append(name); lo, hi = (y if lo is None else lo), y
    return clusters + [((lo + hi) / 2.0, cur)] if cur else clusters


def _star_ground_truth(ax, named_y, tol, x_star, name_color):
    """Star + label each cluster of ground-truth samples at the panel's right edge."""
    for y_center, names in _cluster_by_y(named_y, tol):
        names = sorted(names)
        text = "\n".join(", ".join(names[i:i + 4]) for i in range(0, len(names), 4))
        color = name_color.get(names[0], "black") if len(names) == 1 else "black"
        ax.scatter([x_star], [y_center], s=90, marker="*", color=color,
                   edgecolors="none", zorder=5, clip_on=False)
        ax.annotate(text, (x_star, y_center), textcoords="offset points", xytext=(5, 5),
                    ha="left", va="center", fontsize=7, linespacing=1.0, clip_on=False,
                    path_effects=[pe.withStroke(linewidth=1, foreground="white")])


def make_detail_figure(cfg_name, data, out_dir):
    """4 panels: external rank vs score, external rank vs vdG rank, ligand RMSD
    vs score (ground truth outlined), legend. No-op without ranked sample names."""
    ranked, ground_truth = {}, []   # label -> [(rank, score)]; (label, score, sample)
    for label, (scores, _) in data.items():
        for sample, score in scores.items():
            m = re.search(r"(\d+)$", sample)
            if m:
                ranked.setdefault(label, []).append((int(m.group(1)), score))
            if GT_MARK in sample:
                ground_truth.append((label, score, sample))
    if not ranked:
        return

    color = label_colors(sorted(ranked))
    fig, (ax1, ax2, ax3, axL) = plt.subplots(
        1, 4, figsize=(18, 4.8), gridspec_kw={"width_ratios": [1, 1, 1, 0.12]})

    gt_ranks = {}
    for label, pts in ranked.items():
        ranks, scores_ = zip(*sorted(pts))
        ax1.plot(ranks, scores_, lw=1, color=color[label], label=label)
        gt_scores = [s for (lbl, s, _) in ground_truth if lbl == label]
        # Rank ground truth among the ranked samples without plotting it inline.
        all_ranks = stats.rankdata(-np.array(list(scores_) + gt_scores, float), method="average")
        ax2.plot(ranks, all_ranks[:len(scores_)], lw=1, color=color[label])
        if gt_scores:
            gt_ranks[label] = all_ranks[len(scores_):]

    for label, (scores, rmsd_map) in data.items():
        pts = [(rmsd_map[s], sc, GT_MARK in s) for s, sc in scores.items()
               if rmsd_map.get(s, math.inf) <= 5.0]
        for is_gt, edge, z in ((False, "none", 4), (True, "black", 5)):
            sub = [(r, sc) for r, sc, g in pts if g == is_gt]
            if sub:
                ax3.scatter(*zip(*sub), s=14, color=color.get(label, "black"),
                            edgecolors=edge, linewidths=0.6, zorder=z)

    x_star = max(r for pts in ranked.values() for r, _ in pts) + 0.5
    for ax, xlab, ylab in ((ax1, "external rank", "vdG score"),
                           (ax2, "external rank", "vdG rank"),
                           (ax3, "Ligand RMSD to crystal (Å)", "vdG score")):
        ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    ax1.set_xlim(0, x_star + 0.5); ax2.set_xlim(0, x_star + 0.5); ax3.set_xlim(0, 5.1)
    despine(ax1, ax2, ax3)

    name_color = {s[:4].upper(): color.get(lbl, "black") for lbl, _, s in ground_truth}
    left = [(s[:4].upper(), sc) for (_, sc, s) in ground_truth]
    right = [(s[:4].upper(), float(np.min(gt_ranks[l])))
             for (l, _, s) in ground_truth if l in gt_ranks]
    tol = 0.06 * (max(y for _, y in left) - min(y for _, y in left)) if left else 1
    _star_ground_truth(ax1, left, tol, x_star, name_color)
    _star_ground_truth(ax2, right, 0.45, x_star, name_color)

    axL.axis("off")
    axL.legend(*ax1.get_legend_handles_labels(), loc="center", fontsize=8, frameon=False)
    fig.suptitle(cfg_name)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = os.path.join(
        out_dir, "detail_" + "".join(c if (c.isalnum() or c in "-._") else "_" for c in cfg_name) + ".png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"[INFO] Saved detail figure: {out_path}")


# ---------------------------------- main -------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--matches-root", required=True,
                   help="Root with one subdir per structure (each has ligand_rmsd_summary.tsv + vdg_hits.tsv).")
    p.add_argument("--labels", nargs="+", default=None,
                   help="Structure dir names under matches-root (default: all).")
    p.add_argument("--per-struct-detail", action="store_true",
                   help="Print per-structure Pearson/Spearman (default: condition-level summary only).")
    p.add_argument("--normalize-by-frag-avg", action="store_true",
                   help="Divide each fragment's per-sample count by its mean count across samples before summing.")
    p.add_argument("--detail-figures", action="store_true",
                   help="Also save a per-condition 4-panel figure with external-rank correlation panels "
                        "(needs sample names ending in digits; see make_detail_figure docstring).")
    p.add_argument("--out-fig", default="vdg_score_vs_ligrmsd_conditions.png",
                   help="Output PNG for the overview grid (one panel per condition).")
    return p.parse_args()


def main():
    args = parse_args()
    root = os.path.abspath(args.matches_root)
    labels = args.labels or sorted(
        d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    print("Structures / labels:\n" + "\n".join(f"  - {l}" for l in labels))

    out_dir = os.path.dirname(os.path.abspath(args.out_fig)) or "."
    cond_results = {}
    for base, weights in SIZE_CONFIGS:
        for rmsd_cut in RMSD_CUTS:
            cfg = f"{base}_rmsd{rmsd_cut}"
            data = condition_data(root, labels, rmsd_cut, weights, args.normalize_by_frag_avg)
            pts, per_struct = condition_points(data)
            summarize(cfg, weights, per_struct, detail=args.per_struct_detail)
            if not pts:
                print("  [INFO] No valid (rmsd<5Å) samples for this condition.")

            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            spr = [m["spearman"] for m in per_struct if not math.isnan(m["spearman"])]
            cond_results[cfg] = dict(
                points=pts, pearson=pearson(xs, ys) if xs else float("nan"),
                spearman=spearman(xs, ys) if xs else float("nan"),
                mean_spearman=mean_or_nan(spr), n_structs=len(spr))
            if args.detail_figures:
                make_detail_figure(cfg, data, out_dir)

    plot_conditions(cond_results, args.out_fig)


if __name__ == "__main__":
    main()
