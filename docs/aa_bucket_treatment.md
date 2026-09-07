# How residue slots are labelled, and what that means downstream

Explanation of the treatment of backbone, glycine, proline, and
non-canonical residues. Full detail: `docs/representation_assessment.md` (calibration 
evidence for the backbone-slot clash thresholds), 
`docs/bioisostere_identification.md` (the enrichment null model), 
`docs/pitfalls.md` (interpretation caveats).

## The one idea

A vdG bucket label answers **"which moiety of a residue is contacting the
chemical group?"** — not "which residue is nearby". So the label space is
*(residue, moiety)*, and there are exactly three kinds of label:

| label | meaning |
|---|---|
| `ASP`, `TRP`, … | that residue's **sidechain** is what contacts the CG |
| `bb` | the residue's **backbone** is what contacts the CG |
| `X` | **non-canonical atoms** are what contact the CG |

The true residue name is always kept alongside, in `nr_scrr_resname`. So
`bb` and `X` lose nothing — they say what kind of contact it is, and the resname
says whose.

## Backbone: one label, and why not three

There is a **single** `bb` label, covering every residue including glycine and
proline. This was briefly split into `bb` / `bbGLY` / `bbPRO`, on the reasoning
that glycine has no CB to get in the way and proline has no amide N–H to donate.

The chemistry is real; encoding it as a *partition of the library* was not. A
partition is per-bucket, but the property is per-vdG:

- 54% of glycine-derived backbone geometries are perfectly hostable by a
  CB-bearing residue. `bbGLY` made **100%** of them unreachable — and glycine is
  about half the entire backbone pool.
- ~81% of backbone geometries are hostable by a proline, against the 2.75% that
  `bbPRO` exposed.

Both properties are now tested **on the read path, against the query residue's own
atoms**, which is strictly better than a label describing the *library* residue:

- `BB_SLOT_SIDECHAIN_CLASH` = 3.4 Å — reject if the query residue's sidechain
  (or any non-canonical heavy atom it carries) would sit on top of the CG. Rejects 0.08% of known-good geometries (the 1st
  percentile of real CB-to-CG distance is 4.47 Å), while excluding 46% of
  glycine-derived ones, which is exactly the CB collision `bbGLY` was guarding.
- `PRO_NH_DONOR_CUTOFF` = 3.5 Å — reject a proline query if the geometry needs a
  backbone N–H donation. Validated on 1,531 proline-derived backbone vdGs: none
  donate.

A **sidechain**-labelled slot is deliberately *not* screened — there the sidechain
contact is the claim being tested, not an obstacle to it.

## Glycine

Glycine has no sidechain, so it produces **no `GLY` sidechain bucket, ever**
(verified: no `GLY.npz` in any of 193 fragments). Every glycine observation lives
under `bb`. That is not a gap in the data — it is the chemistry.

The consequence for anything counting residues: **if you look for `GLY` you will
find nothing, and if you pool `bb` you will hide the largest single contributor
to it.** Glycine is ~49% of all backbone clusters library-wide.

## Non-canonical residues (`X`)

`X` is a label, not a resname. It means the residue carries heavy atoms its name
does not have, *and those atoms are what contact the CG*. Detection is by
comparing atom **names** against a canonical table, not by trusting the resname —
that is what catches e.g. a GFP chromophore, which ProDy hands back as an
ordinary sidechain.

One deliberate exception: **selenomethionine** renamed to `MET` with `SE` for `SD`
is a phasing reagent and a methionine isostere, not a functional modification, so
it stays `MET`. It was ~51% of every `X` slot before that carve-out. The rest of
`X` is real: oxidized/adducted CYS ~20%, alkylated LYS ~18%, scattered ncAAs.

`X` buckets are **unreachable by hit finding** — a query structure never produces
an `X` label. That is intentional: the geometry is real and worth keeping, but its
chemical environment is not reproducible, so it must not be counted as an
observation of the residue whose name it wears. Any analysis that aggregates over
buckets has to decide explicitly whether `X` is in scope.

This costs less coverage than it sounds like, because `X` needs the *modification
itself* to be the contacting group (a non-canonical atom within 4.5 Å of the CG,
`align_and_cluster.reorder_vdg_subset`). A modified residue whose canonical part
makes the contact keeps its resname or `bb` label and stays fully reachable, with
`SLOT_MODIFIED` recording the modification separately. On the query side
`query_slot_labels` returns `(resname, bb)`, so a query's modified residue still
matches every backbone bucket; only its non-canonical sidechain identity is
unmatched, which is the same statement.

**Do not read `X` as a census of non-canonical chemistry in the PDB.** Prepwizard
renames some modified residues into standard amino acids (`LLP` → `LYS` keeping its
PLP atoms, chromophores → `GLY`/`ASN`); these keep a full backbone, so they do become
slots, and their foreign atoms land them in `X`. Unreachability is therefore doing
double duty as a quarantine for preprocessing artifacts. `nr_scrr_resname` is
still the true residue, but the *composition* of `X` is partly an artifact of s02.

## Scale (194-fragment library, 3.64 M slots)

| | share |
|---|---|
| sidechain contact | 86.9% |
| backbone contact | 12.7% *(of which glycine-derived 54.7%, proline-derived 2.3%)* |
| `X` | 0.24% |
| modified residue, but contacting canonically | 0.11% |

## What this changes for the bioisostere profiles

Because a label is a *moiety* claim, the enrichment null model has to be a
partition over moieties, not over residues:

    p(category) = bg_freq × size / Σ(bg_freq × size)

with `bg_freq` = "fraction of residues presenting this category" (a sidechain:
that residue's PDB frequency; a backbone: 1.0, since every residue has one) and
`size` = heavy atoms that moiety exposes (sidechain = `AA_size − 4`; backbone = 4).

This sums to exactly 1 and the categories partition each residue's heavy atoms.

Backbone contacts now appear on the profiles in one of two ways:

- **`--bb-mode per-residue`** (default for single-AA): each backbone slot is
  attributed to a residue through its cluster's nr vdG `nr_scrr_resname`,
  giving `GLY-bb`, `ALA-bb`, … This is what puts glycine on the plot.

  Attribution is **per cluster, via its nr vdG**: backbone clustering is on N/CA/C
  only, so geometry cannot separate a glycine donor from an alanine one and a `bb`
  cluster may mix donors. Measured on one large `bb.npz` (1,883 clusters, 24,697
  members): 87% of clusters are single-resname and the nr vdG's resname covers 95%
  of its cluster's members on average. The residual is conservative for the
  headline result -- glycine is 55.8% of clusters by nr vdG against 67.7% of
  members, so nr-vdG attribution *understates* it.
- **`--bb-mode pooled`** (default for the pair heatmap): one `bb` category, so
  the grid stays 20×20 instead of 39×39 with single-digit cells.

`X` is never a category, in either mode — it has no background frequency to divide
by. Its count is reported (`excluded_noncanonical` in the NPZ, annotated on the
bar chart) rather than silently dropped.

**Caveat to carry:** measured over 303 CGs, `GLY-bb` is the *only* backbone
category that is typically enriched — median **+0.45**, positive for **94%** of
CGs — while every other `<AA>-bb` has a negative median (`CYS-bb` −1.10 down to
`ILE-bb` −2.33) and is positive for at most 6%. That mostly says glycine's
backbone is the most exposed one there is, not something about the CG. It is not
a dominant category, though: median rank ~13 of 39, top-5 in only 3% of profiles.
Read backbone rows across CGs, not against zero.
