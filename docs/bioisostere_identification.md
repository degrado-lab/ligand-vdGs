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

Outputs per CG: `<cg>_single_aa_freq{suffix}.npz` and, unless `--skip-pairs`,
`<cg>_aa_pair_freq{suffix}.npz` plus its heatmap `.png`. `--plot-single-aa` additionally writes
the single-AA chart. `{suffix}` encodes norm method + bb-mode (e.g.
`_norm_AA_size_bg_weighted_bbPR`). Skipped CGs (incomplete run or low support) are
reported at the end or to `--skip-log`.

---

## Step 2 — `compare_aa_profiles.py`

```bash
# Single-AA mode (recommended — simpler, more interpretable):
python ligand_vdgs/identify_bioisosteres/compare_aa_profiles.py \
    --profiles-dir outputs/bioisosteres/aa_profiles/ --single-aa --spearman --pearson

# AA-pair mode:
python ligand_vdgs/identify_bioisosteres/compare_aa_profiles.py \
    --profiles-dir outputs/bioisosteres/aa_profiles/ --spearman --pearson
```

`--single-aa` compares `*_single_aa_freq*.npz` (~20-D vectors); default compares
`*_aa_pair_freq*.npz` (upper-triangle, ~210 values). No metric flag → Spearman only.

Metrics: **Pearson** (`--pearson`, recommended — most discriminating), Spearman (`--spearman`,
robust but compresses near 1.0). P-values (scipy) for both are **anticonservative in pair mode**
(upper-triangle entries aren't independent) with no multiple-testing correction — treat as a
ranking aid only.

Quality filters: `--min-shared N` (5) · `--min-coverage F` (0.5, single-AA: shared AAs must cover
≥F of the larger set) · `--weight-by-count`/`--weight-threshold` (30, single-AA only, scales
enrichment by `min(1, count/threshold)`) · `--top-n` (20) · `--exclude-bb` (drops backbone
categories post hoc; denominators still include them) · `--npz FILE ...` (compare specific files
instead of directory discovery).

> `--profiles-dir` discovery loads **every** norm-method/`--bb-mode` variant on disk — if Step 1
> produced multiple, the same CG enters multiple times and its own variants show up as near-1.0
> "pairs". Use `--npz` or a single-variant directory. `bb_mode` mismatches are warned, not blocked.

Outputs (per metric `m`, `_single` suffix in single-AA mode): `similarity_matrix_{m}.npz`
(matrix + labels + p-values), `top_pairs_{m}.tsv` (ranked pairs). No plotting in this step.

---

## W1 figures and fragment-library benchmark

```bash
python ligand_vdgs/identify_bioisosteres/fig1_known_partner_rank.py \
    --profiles-dir <pooled_single_aa_profiles> --outdir <fig1_dir>
python ligand_vdgs/identify_bioisosteres/fig2_clustermap.py \
    --profiles-dir <pooled_single_aa_profiles> \
    --sim-npz <fig1_dir>/sim_matrix_pearson_single_bbPOOL.npz --outdir <fig2_dir>
python ligand_vdgs/identify_bioisosteres/fig3_geometry_overlay.py \
    --vdg-lib-dir <frag_lib> --tsv <geometry_comparison.tsv> --outdir <fig3_dir>
python ligand_vdgs/identify_bioisosteres/bioisostere_benchmark.py \
    --vdg-lib-dir <frag_lib> --profiles-dir <pooled_single_aa_profiles> \
    --outdir <benchmark_dir>
```

Figure 1 reports each preselected partner's reciprocal percentile among finite Pearson
comparisons, excluding graph-nested fragment keys; Figure 2 reuses its matrix. The benchmark
reports tetrazole/carboxylate and carboxylate/carboxylate pairs in `pair_rankings.tsv`, separate
neutral/negative charge partitions in `charge_comparisons.tsv`, parent-disjoint reciprocal ranks
in `parent_holdout_rankings.tsv`, and a parent-disjoint diagnostic for nested aliphatic
carboxylates in `nested_parent_folds.tsv`. Its primary finding criterion is
top 5% in both directions; top 1% and exact ranks are also reported. Charge-partition ranks
require at least 50 parent-support counts in each profile. The partition label comes from the
library's chemistry perception and does not independently prove where a proton sits.

AA profiles describe contact-category preference, not directional hydrogen bonding. Figure 3
and `compare_cg_geometries.py` compare CG centroids in a protein-backbone frame; they do not
establish matching polar atom positions or electrostatics. A graph-nested pair shares source
observations, so its ordinary profile correlation is descriptive; the benchmark's disjoint
parent-entry folds remove direct overlap for that diagnostic.

**Figure 1 caption (N3 amide↔phenol, coordinator ruling 2026-09-25):** AA profiles track shared
donor/acceptor chemistry, so amide↔phenol isn't a clean negative for profiles alone (pct 7.1%,
stronger-ranked than every known pair). Reported as-is, not smoothed over. Geometry
(`mmd2`, same spec as P1) does discriminate: amide/phenol's two highest-support shared single-AA
buckets are PHE (N_A=1084, N_B=691, mmd2=0.137 vs. self floors -0.065/0.002) and, at P1's own
LEU pick, mmd2=0.182 (vs. self floors -0.010/0.072) — both well above their noise floors, unlike
P1 amide/ester at LEU (mmd2=0.047, near its -0.010 floor). Geometry separates P1 from N3 even
though the AA-profile metric alone does not.

---

## Known limitations

- Enrichment reflects PDB deposition bias, not true chemical affinity.
- Heavy-atom size normalization is a proxy for binding surface; SASA would be more principled.
- Single-AA profiles lose cooperativity (e.g. ASP+HIS vs. ASP+LYS) — pair mode captures it, less
  interpretably.
- Zero-count AAs are dropped, not penalized — depletion is invisible to every metric.

---

## Independent branch — `compare_cg_geometries.py`

Compares where the CG sits relative to the interacting backbone between two fragment libraries:
for each shared AA bucket, an unbiased weighted MMD (`mmd2`) between the two libraries' CG
centroids in a fixed ideal-backbone frame. Each vdG's own vdM slot-0 N/CA/C is Kabsch-fit once
onto a constant ideal N-CA-C triangle (never onto another vdG, never including the CG), so the
alignment cannot be biased by the library it's compared against. For 2-residue buckets only slot
0 (`aa_bucket_parts` order) anchors the fit; slot 1 and the CG ride along rigidly.

```bash
python ligand_vdgs/identify_bioisosteres/compare_cg_geometries.py \
    --vdg-lib-dir <path/to/frag_lib> --frags "cnnnn" "CC(=O)[O-]" \
    --outdir outputs/bioisosteres/cg_geometry

# All fragment pairs, AA enrichment only (faster):
python ligand_vdgs/identify_bioisosteres/compare_cg_geometries.py \
    --vdg-lib-dir <path/to/frag_lib> --skip-geometry
```

Flags: `--frags A B ...` (default all) · `--subset-sizes N ...` (`1 2`) · `--max-per-lib N`
(1000, centroids sampled per library per bucket) · `--skip-geometry` · `--no-plot`.

Output TSV caveats: `mmd2` is an unbiased U-statistic estimator of 0 for two samples of the same
population, so it can be slightly negative — it is not itself required to be `>= 0`.
`mmd2_self_A/B` (same library, split by parent biounit) is that bucket's noise floor; calibrate
`mmd2` against it rather than against zero. Both are computed on the `--max-per-lib`-capped
subsample. `enrichment_A/B` describe only the first residue of a multi-residue bucket (`nan` if
that residue has no prior). Backbone buckets are **not** excluded (this script reads its own
bucket counts, not the propensity path).
