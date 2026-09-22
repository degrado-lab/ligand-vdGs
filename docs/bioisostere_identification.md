# Bioisostere Identification via AA Interaction Profiles

CGs (chemical groups) are candidate bioisosteres if they make similar interactions with the
binding site. Compute each CG's AA interaction profile from its vdGs (log-enrichment vs. PDB
background), then compare profiles across CGs.

---

## Step 1 — `compute_aa_profiles.py`

```bash
python ligand_vdgs/identify_bioisosteres/compute_aa_profiles.py \
    --vdglib-dir /path/to/frag_lib --skip-pairs   # fast, single-AA only
    # omit --skip-pairs for AA-pair heatmaps too (~10x slower)
```

Outputs → `outputs/bioisosteres/aa_profiles/<cg>/` (`--outdir` to change).

**Enrichment = log(observed / expected).** Observed = summed `cluster_num_parents` for the
bucket's non-redundant clusters. Expected = total loaded parent-entry support × `p(category)`,
where categories are *(residue, moiety)* pairs (e.g. `ASP` sidechain vs. `GLY-bb` backbone) and
`p(category) = bg_freq × size / Σ(bg_freq × size)` (`size` = sidechain heavy-atom count, or 4 for
any backbone category; glycine has no sidechain category). Pair prior: `p(cat1,cat2) = sf·p1·p2`
(`sf=2` heterogeneous, 1 homogeneous). `total` excludes `X` (non-canonical, ~0.24% of slots) and
is not comparable across `--bb-mode` settings.

Key flags: `--min-clusters N` (default 50, skip low-support CGs) · `--bb-mode
per-residue|pooled|off` (backbone attribution; default `per-residue`, attributes via the
cluster's nr vdG `nr_scrr_resname`) · `--pair-bb-mode` (same, for pairs; default `pooled`) ·
`--norm-AAsize-method` (default `bg_weighted`; `none` requires `--bb-mode off`) · `--plot-single-aa`
· `--skip-log FILE`.

**Caveats:**
- `<AA>-bb` attribution is per-cluster from the nr vdG's resname (backbone geometry alone can't
  distinguish donors); measured 87% single-resname clusters, 95% avg member coverage. Understates
  glycine (55.8% of clusters vs. 67.7% of members).
- Only `GLY-bb` is reliably enriched (median +0.45, positive in 94% of 303 CGs); other `<AA>-bb`
  rows are typically negative — compare backbone rows to each other, not to zero.
- Zero-count AAs are omitted, not recorded as depleted (no pseudocount) — shrinks the shared-AA
  set used in Step 2.
- `X` is never a contact category (no background frequency); tracked separately as
  `excluded_noncanonical(_bb)`.

Outputs per CG: `<cg>_single_aa_freq{suffix}.npz`, `<cg>_aa_pair_freq{suffix}.npz` (unless
`--skip-pairs`), and matching `.png` with `--plot-single-aa`. `{suffix}` encodes norm method +
bb-mode (e.g. `_norm_AA_size_bg_weighted_bbPR`). Skipped CGs (incomplete run or low support) are
reported at the end or to `--skip-log`.

---

## Step 2 — `compare_aa_profiles.py`

```bash
# Single-AA mode (recommended — simpler, more interpretable):
python ligand_vdgs/identify_bioisosteres/compare_aa_profiles.py \
    --profiles-dir outputs/bioisosteres/aa_profiles/ \
    --single-aa --spearman --pearson --jaccard --jaccard-positive-only --score-cutoff 0.99

# AA-pair mode:
python ligand_vdgs/identify_bioisosteres/compare_aa_profiles.py \
    --profiles-dir outputs/bioisosteres/aa_profiles/ --spearman --pearson --jaccard
```

`--single-aa` compares `*_single_aa_freq*.npz` (~20-D vectors); default compares
`*_aa_pair_freq*.npz` (upper-triangle, ~210 values). No metric flag → Spearman only.

Metrics: **Pearson** (`--pearson`, recommended — most discriminating), Spearman (`--spearman`,
robust but compresses near 1.0), top-k Jaccard (`--jaccard`,
`|top_k(A)∩top_k(B)|/|top_k(A)∪top_k(B)|`). P-values (scipy) for Spearman/Pearson are
**anticonservative in pair mode** (upper-triangle entries aren't independent) with no
multiple-testing correction — treat as a ranking aid only.

Quality filters: `--min-shared N` (5) · `--min-coverage F` (0.5, single-AA: shared AAs must cover
≥F of the larger set) · `--weight-by-count`/`--weight-threshold` (30, single-AA only, scales
enrichment by `min(1, count/threshold)`) · `--jaccard-positive-only` (restrict top-k to
both-positive AAs — avoids hydrophobic dominance but inflates scores) · `--score-cutoff F` ·
`--top-n` (20) · `--jaccard-k` (10) · `--exclude-bb` (drops backbone categories post hoc;
denominators still include them) · `--npz FILE ...` (compare specific files instead of
directory discovery).

> `--profiles-dir` discovery loads **every** norm-method/`--bb-mode` variant on disk — if Step 1
> produced multiple, the same CG enters multiple times and its own variants show up as near-1.0
> "pairs". Use `--npz` or a single-variant directory. `bb_mode` mismatches are warned, not blocked.

Outputs (per metric `m`, `_single` suffix in single-AA mode): `similarity_matrix_{m}.npz`
(matrix + labels + p-values), `top_pairs_{m}.tsv` (ranked pairs). No plotting in this step.

---

## Plotting

Not shipped — see `~/docking/scratch/recovered_visualize_bioisosteres/recovered_spec.txt`.
Consider a network graph, clustermap, or UMAP over the similarity matrix.

---

## Known limitations

- Enrichment reflects PDB deposition bias, not true chemical affinity.
- Heavy-atom size normalization is a proxy for binding surface; SASA would be more principled.
- Single-AA profiles lose cooperativity (e.g. ASP+HIS vs. ASP+LYS) — pair mode captures it, less
  interpretably.
- Jaccard without `--jaccard-positive-only` is hydrophobic-biased (MET/CYS/PHE/TRP/TYR dominate),
  and saturates to 1.0 when a filtered profile has exactly `--jaccard-k` entries left.
- Zero-count AAs are dropped, not penalized — depletion is invisible to every metric.

---

## Independent branch — `compare_cg_geometries.py`

Compares where the CG sits relative to the interacting backbone between two fragment libraries:
for each shared AA bucket, Kabsch-aligns backbone + CG centroid across every nr vdG pair and
reports closest approach.

```bash
python ligand_vdgs/identify_bioisosteres/compare_cg_geometries.py \
    --vdg-lib-dir <path/to/frag_lib> --frags "cnnnn" "CC(=O)[O-]" \
    --outdir outputs/bioisosteres/cg_geometry

# All fragment pairs, AA enrichment only (faster):
python ligand_vdgs/identify_bioisosteres/compare_cg_geometries.py \
    --vdg-lib-dir <path/to/frag_lib> --skip-geometry
```

Flags: `--frags A B ...` (default all) · `--match-threshold Å` (1.5, distance below which an nr
vdG counts as matched) · `--subset-sizes N ...` (`1 2`) · `--max-per-lib N` (1000, centroids
sampled per library per bucket) · `--skip-geometry` · `--no-plot`.

Output TSV caveats: `min_dist` is a global minimum over N_A×N_B pairs — shrinks with sample size,
not comparable across buckets with different N (CG centroid is part of the Kabsch fit, deflating
it further); use comparatively within a bucket only. `frac_A/B_matched` are computed on the
`--max-per-lib`-capped subsample. `enrichment_A/B` describe only the first residue of a
multi-residue bucket (`nan` if that residue has no prior). Backbone buckets are **not** excluded
(this script reads its own bucket counts, not the propensity path).
