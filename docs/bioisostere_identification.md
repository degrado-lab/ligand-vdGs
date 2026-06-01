# Bioisostere Identification via AA Interaction Profiles

## Concept


Chemical groups are **bioisosteres** of each other if one can potentially replace the other in a drug molecule because they make similar interactions with the binding site. This workflow identifies candidate bioisosteres by comparing the amino acid interaction profiles of different chemical groups (CGs).


Each CG in the fragment library has a set of vdGs (van der Graphs) — structural exemplars of how that CG interacts with nearby protein residues. By computing which AAs are enriched or depleted in those vdGs relative to their background frequency in the PDB, we get an **AA interaction profile** for each CG. CGs with similar profiles likely interact with similar binding environments, making them candidates for bioisosteric substitution.

## Running the analysis

```bash
python ligand_vdgs/identify_bioisosteres/visualize_aa_profiles.py \
    --smiles "CNC1CC1" \
    --vdglib-dir /wynton/home/degradolab/skt/docking/frag_lib
```

Outputs go to `outputs/aa_profiles/<cg>/` by default, where `<cg>` is the canonical name of the fragment in the library (matching the `frag_lib/` subdirectory name). Pass `--outdir` to override.

By default, only the AA-pair heatmap is produced using `bg_weighted` normalization. The two plot flags are not additive by themselves: passing only `--plot-single-aa-freqs` produces the single-AA bar chart and suppresses the pair heatmap. To get both outputs, pass both flags.

| Flag | Effect |
|---|---|
| `--plot-aa-pair-freqs` | Heatmap: enrichment of each AA pair (subset_size=2 vdGs) — **on by default** |
| `--plot-single-aa-freqs` | Bar chart: enrichment of each individual AA — suppresses pair heatmap unless `--plot-aa-pair-freqs` is also passed |
| `--norm-AAsize-method METHOD` | Size-normalization method(s). Default: `bg_weighted`. Can be repeated, e.g. `--norm-AAsize-method bg_weighted --norm-AAsize-method none` |
| `--outdir DIR` | Output directory for plots and data. Default: `outputs/aa_profiles/<cg>/` |

## Interpreting the output

**Enrichment = log(observed / expected)**, where:

- **Observed** = number of non-redundant vdG *clusters* for that AA bucket. Cluster count is used instead of raw PDB occurrence to avoid inflating counts from overrepresented protein families (e.g. kinases).
- **Expected** = total clusters × background frequency of that AA in the PDB.

For AA-pair buckets, expected is multiplied by 2 for heterogeneous pairs (e.g. ASP_HIS) because pairs are stored alphabetically, so the ASP_HIS bucket captures both orderings (ASP+HIS and HIS+ASP). Homogeneous pairs (ASP_ASP) have only one ordering, so no factor of 2.

**Size normalization** corrects for the fact that larger AAs have more heavy atoms and therefore more geometric opportunities to interact, inflating their raw counts. When normalization is on, observed is divided by the AA's heavy-atom count and expected is divided by an average heavy-atom count, so enrichment reflects preference *beyond* what size and frequency predict.

The default normalization is `bg_weighted`. The other methods are available via `--norm-AAsize-method` and can be repeated to produce multiple outputs in one call.

| Method | Filename suffix | Null model |
|---|---|---|
| `bg_weighted` *(default)* | `_norm_AA_size_bg_weighted` | Observed / AA size; expected / background-frequency-weighted mean size. Enrichment = 0 when count is fully explained by both frequency and size. The theoretically correct normalization. |
| `unweighted` | `_norm_AA_size` | Observed / AA size; expected / unweighted mean size across all 20 AAs. Less principled; provided for comparison. |
| `none` | *(no suffix)* | Raw counts vs. background frequency only (no size correction) |

Backbone-containing buckets (`bb`) are excluded — backbone geometry is conserved across all AAs and does not distinguish which sidechain is interacting.

## Comparing profiles across CGs (bioisostere scoring)

> **Future feature.** Cross-CG profile comparison is not yet implemented. The outputs below are saved to support it when ready.

Each run saves a `.npz` alongside each pair heatmap PNG (single-AA bar charts produce only a PNG, no `.npz`). The filename matches the PNG (suffix depends on norm method):
- `<cg>_aa_pair_freq_norm_AA_size_bg_weighted.npz` (bg_weighted — default)
- `<cg>_aa_pair_freq_norm_AA_size.npz` (unweighted)
- `<cg>_aa_pair_freq.npz` (none)

Example load:

```python
import numpy as np
d = np.load('outputs/aa_profiles/<cg>/<cg>_aa_pair_freq_norm_AA_size_bg_weighted.npz', allow_pickle=True)
matrix = d['matrix']   # NxN enrichment values (NaN = pair not observed)
aas    = list(d['aa_labels'])
```

To compare two CGs, align their matrices by `aa_labels` (row/column order differs per CG) and compute a similarity score (e.g. Pearson correlation on non-NaN entries). CGs whose profiles correlate strongly are candidates for bioisosteric substitution.
