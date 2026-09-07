# Calibration evidence for the backbone-slot and CG-contact thresholds

This file is evidence, not narrative: the numbers behind three constants whose
rationale is already stated in code, kept here so the thresholds can be checked
rather than taken on faith. All measurements were made once, pre-rebuild, on a
193-fragment library (`delete_frag_lib`, 534,460 size-1 clusters). None of the
three constants depend on which library was mined — they're calibrated against
real amino-acid geometry (CB/CD position, H-bond donor distance) — so they don't
need re-measuring after a rebuild unless one starts looking wrong in practice.

---

## 1. Backbone-slot clash thresholds

`hit_finder_core.BB_SLOT_SIDECHAIN_CLASH` (3.4 Å) and `PRO_NH_DONOR_CUTOFF`
(3.5 Å) let a single `bb` bucket label stay safe without splitting by donor
residue: at read time, the query residue's real sidechain/CD/backbone-N is
checked against the CG directly (`backbone_slot_blockers`).

**Sidechain clash (stands in for a glycine-specific check).** Reconstructed CB
position (standard virtual-CB from N/CA/C) to nearest CG atom, over all
`bb`-labeled nr vdGs:

| population | p1 | p5 | p10 | median | n |
|---|---|---|---|---|---|
| GLY-derived | 1.90 | 2.21 | 2.41 | 3.56 | 27,506 |
| non-GLY-derived | 4.47 | 4.56 | 4.63 | 5.21 | 28,072 |

The non-GLY row is real backbone contacts at CB-bearing residues, so its
clearance is by construction non-clashing (p1 = 4.47 Å validates the virtual-CB
reconstruction). At 3.4 Å the clash test misclassifies 0.08% of these
known-good geometries while correctly passing 54% of the glycine population —
geometry a per-donor-residue split would have blocked outright regardless of
whether it actually clashes.

**Donor cutoff (stands in for a proline-specific check).** Of 1,531
proline-derived `bb` nr vdGs — whose real backbone N demonstrably did not
donate — zero have a CG acceptor (N/O/S) within 3.5 Å of N. `PRO_NH_DONOR_CUTOFF`
is validated at that value with no observed false rejections.

Proline's CD (which sits where the donated H would go) is not a separate
constant: `split_residue_heavy_atoms` includes CD in proline's canonical
sidechain set, so `BB_SLOT_SIDECHAIN_CLASH` screens it as an ordinary sidechain
clash. A dedicated virtual-CD reconstruction was used only to confirm this
during calibration (fitted on 156 real prolines, mean deviation 0.25 Å from the
rigid pyrrolidine geometry) — it found no case the general clash test misses.

---

## 2. Cluster-size caveat: long sidechains resolve more loosely

Clustering uses only CG + vdM N/CA/C, so the sidechain itself never enters the
metric — including for buckets whose entire meaning is "this sidechain makes
the contact." Measured on `CC(=O)[O-]` by re-deriving each member's real
sidechain and fitting it onto its nr vdG under the clustering frame:

| bucket | tip atom | rotatable χ | tip RMSF @ 0.5 Å tol. |
|---|---|---|---|
| ASP | OD1/OD2 | 2 | 0.62 Å |
| ARG | CZ | 4 | 0.66 Å |
| LYS | NZ | 4 | 1.14 Å |

Arginine, despite four rotatable angles, is as tight as aspartate — fixing the
CG relative to the backbone constrains a bidentate guanidinium almost fully.
Lysine is the outlier at roughly 2x: an ammonium's single-point, non-directional
contact lets its CH2 chain reach the same position by many paths. Net: a
residue-dependent caveat for cross-residue `cluster_size` comparisons in
`identify_bioisosteres`, not a reason to change the clustering metric.

---

## 3. CG-to-vdM contact cutoff

`clus_and_deduplicate_vdgs.CG_VDM_CONTACT_CUTOFF` (4.5 Å) rejects a vdM whose
nearest heavy atom is farther than that from every CG atom. Real contacts sit
well inside it — over 2,503 sampled nr vdGs, CG-to-nearest-heavy-atom distance:

| p50 | p75 | p90 | p95 | p99 | p99.9 | max |
|---|---|---|---|---|---|---|
| 3.50 | 3.76 | 3.98 | 4.12 | 4.30 | 5.87 | 9.05 |

The distribution has a hard edge around 4.4 Å and then a sparse tail; nothing
in 4.5–10 Å trades away real data. This guard backstops two independent
mechanisms that used to slip vdMs 12+ Å from the CG into the library: (1) a vdM
selected by proximity to the whole ligand rather than to the CG specifically
(large cofactors — FAD, HEM, GSH — where a contact near one end got attributed
to a CG matched at the other), and (2) duplicate atom names in prepwizard
output, the mechanism `docs/pitfalls.md` documents under "Prepwizard relabels
atoms it cannot build as part of a ligand." The guard runs per vdM slot and
fails the whole environment on a miss, so at subset size 2 one non-contacting
slot can discard a genuine partner — measured at 0.133% of size-2 nr vdGs, and
not worth a partial-slot rewrite at that rate (the contacting slot is not lost:
size-1 vdGs are enumerated independently). Revisit if `MAX_SUBSET_SIZE` rises,
since the odds of at least one bad slot grow with subset size.
