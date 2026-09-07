# Bioisostere Identification via AA Interaction Profiles

## Concept

Chemical groups (CGs) are **bioisosteres** of each other if one can potentially replace the other in a drug molecule because they make similar interactions with the binding site. This workflow identifies candidate bioisosteres by comparing the amino acid (AA) interaction profiles of different CGs.

Each CG in the fragment library has a set of vdGs — structural exemplars of how that CG interacts with nearby protein residues. By computing which AAs are enriched or depleted in those vdGs relative to their background frequency in the PDB, we get an AA interaction profile for each CG. CGs with similar profiles likely interact with similar binding environments, making them candidates for bioisosteric substitution.

---

## Step 1 — Compute AA profiles (`compute_aa_profiles.py`)

For every completed CG in the library, compute log-enrichment of each AA relative to background. Outputs a single-AA `.npz` per CG and an AA-pair heatmap (see the outputs table for exactly when each is written).

```bash
# Fast run — single-AA profiles only:
python ligand_vdgs/identify_bioisosteres/compute_aa_profiles.py \
    --vdglib-dir /wynton/home/degradolab/skt/docking/frag_lib \
    --skip-pairs

# Full run including AA-pair heatmaps:
python ligand_vdgs/identify_bioisosteres/compute_aa_profiles.py \
    --vdglib-dir /wynton/home/degradolab/skt/docking/frag_lib
```

Outputs go to `outputs/bioisosteres/aa_profiles/<cg>/` by default.

### Key flags

| Flag | Default | Effect |
|---|---|---|
| `--vdglib-dir DIR` | *(required)* | Path to fragment library root (`frag_lib/`) |
| `--skip-pairs` | off | Skip AA-pair computation (~10× faster; single-AA NPZ still saved) |
| `--plot-single-aa` | off | Also plot single-AA enrichment bar charts |
| `--min-clusters N` | 50 | Skip CGs with fewer than N total vdG clusters (too sparse for reliable statistics) |
| `--norm-AAsize-method METHOD` | `bg_weighted` | Size-normalization method. Can be repeated for multiple outputs |
| `--outdir DIR` | `outputs/bioisosteres/aa_profiles/` | Root output directory |
| `--skip-log FILE` | stdout | Write skipped CG names to FILE. Appends if the file exists, so use a fresh path per run. |
| `--bb-mode MODE` | `per-residue` | How backbone-mediated contacts enter the single-AA profile. `per-residue` attributes every backbone slot to a residue via its cluster's nr vdG `nr_scrr_resname`, giving categories `GLY-bb`, `ALA-bb`, … (see the accuracy caveat below); `pooled` keeps one `bb` category; `off` drops backbone buckets entirely. Backbone categories are hatched in the bar chart and their heatmap ticks are brown/bold. |
| `--pair-bb-mode MODE` | `pooled` | Same, for the AA-pair heatmap. Pooled by default because per-residue would make a 39×39 grid whose rarer cells hold single-digit counts. Mixed buckets (`ALA_bb`, `bb_bb`) are included. |

A sidechain category is sized by its sidechain heavy-atom count only, which is
what the bucket label claims and makes the categories partition the residue.

### Enrichment formula

**Enrichment = log(observed / expected)**, where:

- **Observed** = number of non-redundant vdG *clusters* for that AA bucket. Cluster count is used instead of raw PDB occurrence to avoid inflating counts from overrepresented protein families (e.g. kinases).
- **Expected** = *total clusters × p(category)*, where `p` is the prior below.
- **`total`** is the sum of the counts actually loaded — every `.npz` bucket in
  `nr_vdgs/<size>/` except `X` (always dropped) and, under `--bb-mode off`, every
  bucket with a backbone label in it. It is *not* the fragment's total vdG count,
  and it shifts when `--bb-mode` changes which buckets are read, so enrichments
  computed under different modes are not directly comparable.

Positive values = the category interacts with this CG more than chance; negative = less.

### Contact categories and the null model

A bucket label is a claim about **which moiety** of a residue contacts the CG, not
just which residue: `ASP` means "ASP's sidechain", `bb` means "a backbone". So the
null partitions contacts by *(residue, moiety)*:

| category | `bg_freq` — fraction of residues presenting it | `size` — heavy atoms it exposes |
|---|---|---|
| `ASP` (sidechain) | PDB frequency of ASP | sidechain heavy atoms (`AA_size − 4`) |
| `GLY-bb` (per-residue backbone) | PDB frequency of GLY | 4 (N, CA, C, O) |
| `bb` (pooled backbone) | 1.0 — *every* residue has one | 4 |

    p(category) = bg_freq × size / Σ(bg_freq × size)

The denominator is the background-weighted mean residue size, so **the priors sum
to exactly 1 in every backbone mode**, including the pooled one, because
`Σ_r freq(r)·n_sc(r) + 1.0·4 = Σ_r freq(r)·AA_size(r)`. Glycine has no sidechain
category by construction (0 sidechain heavy atoms — which is exactly why no
`GLY.npz` is ever generated), and under `per-residue` its share of the prior is
claimed by `GLY-bb`.

**Pairs** use the independent-draw product, `p(cat1, cat2) = sf × p1 × p2`, with
`sf = 2` for a heterogeneous bucket (it absorbs both orderings) and 1 for a
homogeneous one. Each `p` already carries its size term, so size normalization is
applied once, and `Σ sf·p1·p2 = (Σp)² = 1`.

**How exact is `<AA>-bb`?** Attribution is **per cluster, via its nr vdG**: backbone clustering is on N/CA/C
only, so geometry cannot separate a glycine donor from an alanine one and a `bb`
cluster may mix donors. Measured on one large `bb.npz` (1,883 clusters, 24,697
members): 87% of clusters are single-resname and the nr vdG's resname covers 95%
of its cluster's members on average. The residual is conservative for the
headline result -- glycine is 55.8% of clusters by nr vdG against 67.7% of
members, so nr-vdG attribution *understates* it.

**Reading the backbone rows.** **`GLY-bb` is the only backbone category that is typically enriched.** Measured
over 303 CGs: `GLY-bb` has median enrichment **+0.45** and is positive for **94%**
of CGs, while every one of the other 19 `<AA>-bb` categories has a *negative*
median (`CYS-bb` −1.10 down to `ILE-bb` −2.33) and is positive for at most 6% of
CGs. That contrast is mostly a statement about glycine's backbone being the most
exposed one there is, not about any particular CG. It is not, however, a
dominant category: its median rank is ~13 of 39 (67th percentile), and it lands
in a profile's top 5 only 3% of the time. **Compare backbone rows against each
other and across CGs, not against zero.**

> An AA with an observed count of exactly zero is **omitted** from the profile
> rather than recorded as strongly depleted, and no pseudocount is applied. A
> fully-depleted AA is therefore indistinguishable from one that was never
> sampled, which shrinks the shared-AA set used in Step 2.

### Size normalization

Larger AAs have more heavy atoms and therefore more geometric opportunities to interact, which inflates their raw counts. Size normalization corrects for this.

| Method | Filename suffix | Description |
|---|---|---|
| `bg_weighted` *(default)* | `_norm_AA_size_bg_weighted` | `p = bg×size / Σ(bg×size)`. Priors sum to exactly 1; enrichment = 0 when a count is fully explained by frequency and size. Theoretically correct. |
| `none` | *(no suffix)* | `p = bg / Σbg` — no size term. **Requires `--bb-mode off` and `--pair-bb-mode off`**: without a size term there is no defensible way to weigh a 4-atom backbone against a whole sidechain, and inventing one is the approximation this design removed. |

### Outputs per CG

| File | Produced when |
|---|---|
| `<cg>_single_aa_freq{suffix}.npz` | When `nr_vdgs/1` exists and at least one propensity is computable. Contains `aa_labels`, `enrichments`, `counts`, `total_count`, plus the provenance needed to tell two profiles apart: `bb_mode`, `excluded_noncanonical` (the `X` clusters), `excluded_noncanonical_bb` |
| `<cg>_aa_pair_freq{suffix}.npz` | Unless `--skip-pairs`, **and** `nr_vdgs/2` exists with non-empty pair propensities. Written as a side effect of drawing the heatmap. Carries `matrix`, `aa_labels`, and the same provenance keys as the single-AA file. |
| `<cg>_single_aa_freq{suffix}.png` | With `--plot-single-aa` |
| `<cg>_aa_pair_freq{suffix}.png` | Same condition as the pair `.npz` |

`{suffix}` is the norm-method suffix from the table above, then the `--bb-mode`
tag (`_bbPR`, `_bbPOOL`, `_bbOFF`) — e.g.
`cn(c)C_single_aa_freq_norm_AA_size_bg_weighted_bbPR.npz`. The mode is in the
name because it is part of what the numbers mean: without it, two `--bb-mode`
runs into one `--outdir` would overwrite each other, and the flags exist to
invite exactly that comparison.

Plot titles use the on-disk directory name, which is the filename-escaped SMILES
(`/` → `_fs_`, `\` → `_bs_`), not the original SMILES.

### Non-canonical residues (`X`)

`X` is a slot label meaning *non-canonical atoms are what contact the CG* — the
real residue name is still in `nr_scrr_resname`. `X` is **never** a contact
category, at any `--bb-mode`: it has no background frequency to divide by, and hit
finding can never produce it, so counting one as an observation of the residue
whose name it wears would be wrong. It is dropped from numerator and denominator
alike. The count is written to the NPZ as
`excluded_noncanonical` and annotated on the bar chart (`N=13,952 (+20 X
excluded)`). Library-wide it is 0.24% of slots, about half of which is
selenomethionine.

Under `--bb-mode per-residue`, a backbone slot whose residue is not canonical is
likewise dropped and counted in `excluded_noncanonical_bb`; under `pooled` it stays
in `bb`.

### CG filtering

CGs are skipped for two reasons, reported at the end (or to `--skip-log`):
- **Incomplete**: no `"Job completed."` line in the log file
- **Low-count**: total single-AA cluster count < `--min-clusters` (default 50). These have too few examples for reliable enrichment statistics.

---

## Step 2 — Compare profiles (`compare_aa_profiles.py`)

Compute pairwise similarity between all CG profiles and rank candidate bioisostere pairs.

```bash
# Single-AA mode (recommended starting point — simpler and more interpretable):
python ligand_vdgs/identify_bioisosteres/compare_aa_profiles.py \
    --profiles-dir outputs/bioisosteres/aa_profiles/ \
    --single-aa --spearman --pearson --jaccard \
    --jaccard-positive-only --score-cutoff 0.99

# AA-pair mode (2-vdM profiles; more detail but harder to interpret):
python ligand_vdgs/identify_bioisosteres/compare_aa_profiles.py \
    --profiles-dir outputs/bioisosteres/aa_profiles/ \
    --spearman --pearson --jaccard
```

### Modes

| Mode | Flag | Input files | Description |
|---|---|---|---|
| Single-AA *(recommended)* | `--single-aa` | `*_single_aa_freq*.npz` | Compares 1-D enrichment vectors (~20 values). Simpler and more interpretable. Output filenames are suffixed `_single`. |
| AA-pair | *(default)* | `*_aa_pair_freq*.npz` | Compares upper-triangle of the NxN AA-pair enrichment matrix (~210 values). More detail, harder to interpret. |

### Metrics

Three metrics are available. No metric flags → Spearman only. Any metric flag(s) → only those run.

| Metric | Flag | Range | Formula | Notes |
|---|---|---|---|---|
| Spearman ρ | `--spearman` | −1 to 1 | Rank correlation of enrichment vectors | Robust to outliers; asks if CGs rank AAs in the same order. With ~20 AAs, scores compress near 1.0 — limited discriminatory power. |
| Pearson r | `--pearson` | −1 to 1 | Linear correlation of enrichment values | Sensitive to magnitude differences, not just rank. Most discriminating for this application. **Recommended.** |
| Top-k Jaccard | `--jaccard` | 0 to 1 | \|top_k(A) ∩ top_k(B)\| / \|top_k(A) ∪ top_k(B)\| | Directly interpretable. See `--jaccard-positive-only` below. |

P-values are computed for Spearman and Pearson (via `scipy`) and included in the TSV and console output.

> **Caveat on pair-mode p-values.** In pair mode the correlated vectors are the
> upper triangle of a symmetric AA-pair matrix, so each AA recurs in ~N entries
> and the ~210 entries are far from independent. `scipy` assumes independence, so
> the reported p-values are anticonservative. No multiple-testing correction is
> applied across the O(n²) CG pairs either. Treat them as a ranking aid, not
> as inference.

### Quality filters

| Flag | Default | Effect |
|---|---|---|
| `--min-shared N` | 5 | Min shared non-NaN entries required to report a score |
| `--min-coverage F` | 0.5 | *(single-AA mode)* Shared AAs must cover ≥ F of the larger profile's AA set |
| `--weight-by-count` | off | *(single-AA mode)* Scale each AA's enrichment by min(1, count/threshold) before correlation. Entries with few supporting clusters are pulled toward zero. |
| `--weight-threshold N` | 30 | Count threshold for confidence weighting (used with `--weight-by-count`) |
| `--jaccard-positive-only` | off | Restrict Jaccard top-k to AAs enriched (> 0) in **both** profiles. Prevents universally hydrophobic AAs from dominating the intersection, but because it requires both, AAs enriched in only one CG are dropped from the union — which inflates the score. |
| `--score-cutoff F` | none | Only report pairs with score ≥ F in TSV and console. More principled than `--top-n`. |
| `--top-n N` | 20 | Number of pairs printed to the console |
| `--jaccard-k N` | 10 | Size of the top-k set used by the Jaccard metric |
| `--npz FILE ...` | none | Compare specific profile NPZ files instead of discovering them under `--profiles-dir` |
| `--exclude-bb` | off | Drop backbone categories (`bb`, `<AA>-bb`) before comparing, so similarity is decided by sidechain chemistry alone. The enrichments themselves are not recomputed, so the values still carry the backbone categories in their denominator. |

> **`--profiles-dir` discovery does not filter by norm method.** If Step 1 was run
> with several `--norm-AAsize-method` values, every variant of the same CG is on
> disk, and discovery matches all of them — the same CG then enters the similarity
> matrix two or three times and its variants surface as near-1.0 "pairs" at the top
> of `top_pairs_*.tsv`. The same applies to `--bb-mode` variants.
> Step 3 has `--norm-suffix` for this; **Step 2 has no equivalent.** Use `--npz`
> to name the one variant you want, or point `--profiles-dir` at a directory
> holding a single variant's output.
>
> Step 2 *does* warn when the loaded profiles disagree on the `bb_mode`
> recorded in their NPZs. The two profiles align on labels and correlate
> silently while their denominators differ, so recompute rather than reading the
> numbers.

> `--weight-by-count`, `--weight-threshold` and `--min-coverage` apply to
> single-AA mode only and are silently ignored in pair mode.

### Recommended score cutoffs (single-AA Pearson)

Pair counts below are from one particular run over the acetate/tetrazole-family
library; they are illustrative, not general. Re-derive them for your own library.

| Cutoff | Pairs returned | Interpretation |
|---|---|---|
| ≥ 0.999 | ~2 | Very high confidence only |
| ≥ 0.995 | ~18 | High confidence, compact set |
| ≥ 0.990 | ~37 | Broader high-confidence set |

### Outputs (per active metric `m`)

| File | Description |
|---|---|
| `similarity_matrix_{m}.npz` | Raw similarity matrix + CG labels + p-values (Spearman/Pearson) |
| `top_pairs_{m}.tsv` | All pairs ranked by score; includes N columns (single-AA mode) and p-values |

In single-AA mode both filenames carry an extra `_single` suffix, e.g.
`similarity_matrix_pearson_single.npz`.

Plotting (heatmaps, network graph, UMAP) happens entirely in Step 3
(`visualize_bioisosteres.py`), which reads the `.npz` this step produces.

---

## Step 3 — Visualize results (`visualize_bioisosteres.py`)

Produces three publication-quality figures from the similarity matrix and (optionally) the raw profiles.

```bash
# Network + clustermap only:
python ligand_vdgs/identify_bioisosteres/visualize_bioisosteres.py \
    --sim-npz outputs/bioisosteres/similarity/similarity_matrix_pearson_single.npz

# All three including UMAP:
python ligand_vdgs/identify_bioisosteres/visualize_bioisosteres.py \
    --sim-npz outputs/bioisosteres/similarity/similarity_matrix_pearson_single.npz \
    --profiles-dir outputs/bioisosteres/aa_profiles/ \
    --umap
```

### Outputs

| File | Description |
|---|---|
| `network_{metric}.png` | Network graph: nodes = CGs, edges = pairs above `--network-threshold`. Connected components are bioisostere clusters. Edge thickness and color encode score. |
| `clustermap_{metric}.png` | Seaborn clustermap of top `--top-cgs` CGs with explicit dendrograms. Easier to read than a plain heatmap. |
| `umap_profiles_{metric}.png` | 2-D UMAP of raw enrichment profiles (one point per CG). Proximity = similar AA preference. Points colored by log₁₀(total clusters). |

### Key flags

| Flag | Default | Effect |
|---|---|---|
| `--sim-npz FILE` | *(required)* | Similarity matrix NPZ from `compare_aa_profiles.py` |
| `--profiles-dir DIR` | none | Directory of per-CG single-AA NPZ files (required for `--umap`) |
| `--umap` | off | Produce UMAP embedding (requires `umap-learn`) |
| `--network-threshold F` | 0.995 | Minimum score to draw a network edge |
| `--top-cgs N` | 20 | Number of top CGs shown in the clustermap |
| `--min-clusters N` | 50 | Min cluster count for a CG to appear in the UMAP |
| `--min-aa-freq F` | 0.5 | For UMAP: only include AAs present in ≥ F of profiles |
| `--norm-suffix S` | `_norm_AA_size_bg_weighted_bbPR` | Which profile variant to load — everything after `_single_aa_freq`, so norm method *and* `--bb-mode` tag. Must match Step 1's suffix exactly, or no profiles load. |
| `--umap-neighbors N` | 10 | UMAP `n_neighbors` |
| `--umap-label-top N` | 40 | How many points to label in the UMAP (0 = all) |
| `--metric NAME` | inferred | Metric name used in output filenames |
| `--outdir DIR` | `outputs/bioisosteres/similarity/` | Output directory |

`--min-clusters` here compares against the `total_count` stored in the NPZ, which
includes the `bb` bucket unless Step 1 was run with `--bb-mode off`. The effective
threshold is then not identical to Step 1's.

### Choosing `--network-threshold`

Based on observed Pearson score distributions:

Same caveat as above — these counts came from one library and are illustrative.

| Threshold | Typical edges | Use case |
|---|---|---|
| 0.999 | ~2 pairs | Very high confidence only |
| 0.995 | ~18 pairs | Recommended starting point |
| 0.990 | ~37 pairs | Broader view |

---

## Interpreting results

**Network graph** is the easiest to explain to others: nodes that are connected are predicted bioisosteres. Connected components (islands) are bioisostere clusters — e.g. all aliphatic alcohols in one cluster, all aromatic N-heterocycles in another.

**Clustermap** shows the same information as a matrix, with dendrograms revealing the hierarchical grouping. Useful for methods figures.

**UMAP** shows the full CG landscape: nearby points have similar AA preferences. Points far from others have unique interaction profiles with few bioisostere candidates. Color by total cluster count to identify which CGs have reliable statistics.

**P-values** for Spearman/Pearson are in the TSV (`pvalue` column). With ~20 AAs, ρ > 0.45 is already p < 0.05, so top pairs will all have p ≈ 0 — but p-values are useful for filtering out weakly-similar pairs lower in the ranking. In pair mode they are anticonservative (see the caveat in Step 2).

### Known limitations

- **PDB bias**: enrichment is relative to PDB background AA frequencies, which reflect what proteins have been crystallized, not true chemical affinities.
- **Size normalization**: heavy-atom count is a proxy for binding surface; SASA would be more principled but requires per-structure computation.
- **Single-AA profiles lose cooperativity**: two CGs may prefer ASP individually but differ in whether they see ASP+HIS vs. ASP+LYS together. Pair-mode profiles capture this at the cost of interpretability.
- **Score compression**: with ~20 AAs, Spearman scores compress near 1.0. Pearson is more discriminating.
- **Jaccard hydrophobic bias**: without `--jaccard-positive-only`, top-k is dominated by universally hydrophobic/aromatic AAs (MET, CYS, PHE, TRP, TYR) across many CG types.
- **Zero-count AAs are dropped, not penalized** (see Step 1), so strong depletion is invisible to every metric.
- **Jaccard saturates**: if a profile has exactly `--jaccard-k` entries left after filtering, the top-k set is the whole set and the score is trivially 1.0.

---

## Independent branch — CG geometry (`compare_cg_geometries.py`)

Rather than comparing AA profiles, this compares where the CG itself sits
relative to the interacting backbone, between two fragment libraries. For each
shared AA bucket it Kabsch-aligns backbone + CG centroid across every pair of
nr vdGs and reports the closest approach.

```bash
python ligand_vdgs/identify_bioisosteres/compare_cg_geometries.py \
    --vdg-lib-dir <path/to/frag_lib> \
    --frags "cnnnn" "CC(=O)[O-]" \
    --outdir outputs/bioisosteres/cg_geometry

# All fragment pairs, AA enrichment only (much faster):
python ligand_vdgs/identify_bioisosteres/compare_cg_geometries.py \
    --vdg-lib-dir <path/to/frag_lib> --skip-geometry
```

| Flag | Default | Effect |
|---|---|---|
| `--vdg-lib-dir DIR` | *(required)* | Fragment library root |
| `--frags A B ...` | all | Which fragment libraries to compare |
| `--match-threshold Å` | 1.5 | Distance below which an nr vdG counts as "matched" |
| `--subset-sizes N ...` | `1 2` | Subset sizes to include |
| `--max-per-lib N` | 1000 | Centroids sampled per library per bucket |
| `--skip-geometry` | off | AA enrichment only |
| `--no-plot` | off | Skip the heatmap |

Reading the output TSV:

- `min_dist` is the **global minimum** over all N_A × N_B nr vdG pairs, so it
  decreases as sample size grows and is not comparable across buckets with very
  different N. The CG centroid is part of the Kabsch fit, which further deflates
  it. Use it comparatively within a bucket, not as an absolute.
- `frac_A_matched` / `frac_B_matched` are computed on the subsample capped by
  `--max-per-lib`, not the full library.
- `enrichment_A` / `enrichment_B` describe only the **first** residue of a
  multi-residue bucket; single-AA enrichment is not defined for a bucket as a whole.
  They come from the single-AA propensity path, so they carry that path's scope:
  a bucket whose only non-backbone part has no prior reports `nan` here even
  though `N_A`/`N_B` and the geometry columns are populated.
- Backbone buckets (`bb`, `ASP_bb`, …) are **not** excluded: this script reads its
  own bucket counts rather than the propensity path's, so `N_A`/`N_B`, `min_dist`
  and `frac_*_matched` are computed for them like any other bucket.
