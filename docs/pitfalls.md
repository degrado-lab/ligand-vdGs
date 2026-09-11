# Pitfalls

Things that will mislead you if you don't already know them. Keyed on
function/file names, not line numbers — re-verify after refactors.

- **Assumptions** — premises about inputs and scope. Violating one gives wrong or
  missing results, usually without an error.
- **Caveats** — works as designed, but means something narrower than it looks, or
  carries an invariant you must respect. Not bugs.
- **Known issues** — real limitations, **settled**: alternatives were measured and
  were worse. Don't re-open.
- **Open bugs** — code does something other than intended.

---

## Assumptions

### Exactly one ligand resname per query structure

Copies are fine. With more than one resname, `score_one_model`'s
`n_heavy * n_residues` tie-break picks one; failing that, hit finding prints
`[ERROR]` and returns empty **for that model only**, without aborting the run.

### A ligand is a single HETATM residue

Nothing treats a covalently linked chain of HETATM residues (a glycan like
NAG-NAG-BMA) as one ligand — no LINK parsing, no cross-residue connectivity
merge. **Mining** (`interactions.py`) enumerates one PDB residue at a time, so
each monomer is silently mined as its own ligand. **Hit finding** used to do the
same for homopolymers sharing a resname (two `NAG` fed to RDKit as one molecule);
it now raises `ValueError` from a residue-instance count check in
`Frags.get_query_ligand_mol`. Fail-hard is the only fix applied — real support
means reworking both the mining enumeration and the single-residue assumption.

### Ligand–protein covalent bonds are not represented

`find_cg_matches` reads one HETATM residue at a time and OpenBabel perceives
bonds only within it. A CG spanning the link can't be mined anyway (smallest CG
is 4 connected heavy atoms); the protein-side modified CYS is labelled `X` if its
extra atoms contact the CG, excluding it from hit finding. Covalent-warhead
terminal methyls (`CF0`/`0QE`) are dropped in preprocessing. Don't read the
library as saying anything about covalent binders (11/65,600 structures
discarded).

### The query ligand's tautomer is whatever the CCD template says

Query-side bond orders come from the template SMILES in
`AssignBondOrdersFromTemplate`, so the tautomer is the CCD's choice, not
density-derived. A different tautomer changes the fragment key at the tautomeric
atom: 2-hydroxypyridine gives `cc(n)O`, 2-pyridone `cc(n)=O`. (Ring-only
fragments are shared; `[nH]` vs `n` never splits — the bracket-H regex strips
them.)

**Settled: the split is correct, tautomers aren't enumerated.** A hydroxyl donor
vs. a carbonyl acceptor is exactly the distinction the library exists to make;
pooling merges different H-bonding roles, worse than a miss. Measured with
`rdMolStandardize.TautomerEnumerator` (keys from deposited form vs. all forms):
2-pyridone 3→7, guanine 5→10, imatinib 10→12, warfarin 5→12. Querying the union
roughly doubles fragments per ligand with forms the structure doesn't have, and
some tautomers dearomatize a ring. A tautomer decision belongs at the parent,
once, before fragmentation — not a fan-out at match time.

**H-stripping does pool donor and acceptor aromatic N (resolved by the rebuild's `D<n>` + carbon `H0`/`!H0` keys; evidence below).** The
bracket-H regex writes pyrrole-type `[nH]` and pyridine-type `n` as the same
bare `n`, and keys record no degree, so an N-substituted ring N (`n(R)`, neither
role) matches too. Unlike the pyridone split above these are not tautomers of
one another — pyridine vs. pyrrole is a fixed difference — so it violates the
rule in the previous paragraph. Measured 2026-09-09 over the 29
production-threshold keys (≥250 ligands) containing aromatic `n`, classifying
each matched N site by the BioLiP CCD SMILES (79k sites; ~8% of ligands had no
parsable SMILES): 15% `[nH]` donor, 54% acceptor, 31% N-substituted. Per key it
runs from acceptor-dominated (`cc(n)[N;!R]`: 2/92/6%) to donor-dominated
(`cc(n)=[O;!R]`: 63/5/32%); `cc(c)n` (9,896 ligands) is 21/56/22%. The cost
falls mostly on statistics, not geometry — the protein side sorts the roles into
different buckets and clusters, but `cluster_num_parents` and per-key totals are
normalized over the mixed pool. Options, safest first: (1) degree in the key
(`[n;D2]` vs `[n;D3]`) is graph-only and needs no H trust, but touches every `n`
key (rebuild + `fragment_keys_equivalent`); (2) record each CG atom's H count
per observation from the protonated mirror at mining time (cf. `nr_slot_flag`)
and filter at read time — reversible, no library split; (3) `[nH]`/`[n;H0]` in
the key — splits genuine tautomers (imidazole N1/N3) and inherits prepwizard's
and the CCD's H placement, which the pipeline otherwise refuses to trust.

The substituted share is one instance of a general gap: keys record no degree,
so any heteroatom on a fragment's boundary has unknown substitution. Same
measurement over all 178 production keys: 47 have a boundary heteroatom whose
substitution is genuinely mixed (5–95% of sites) — phenol vs. aryl ether
`cc(c)[O;!R]` (63% substituted, 11,129 ligands), aniline `cc(c)[N;!R]` (90%),
carboxylic acid vs. ester `[C;!R][C;!R](=[O;!R])[O;!R]` (21%), primary vs.
higher amine `[C;!R][C;!R]([C;!R])[N;!R]` (70%), sulfonamide N
`c[S;!R]([N;!R])(=[O;!R])=[O;!R]` (68%). Option (1) generalizes to `D<n>` on
every heteroatom. Caveat: SMARTS `D` counts *explicit* connections, so both
matchers must see H-free graphs — the miner reads the protonated mirror through
OpenBabel, where H's are real atoms (cf. the `r<n>` disagreement below). Unlike `r<n>`, the two toolkits agree on what
`D` means; only the input differs, so `DeleteHydrogens()` before `Match` is a
complete fix (it renumbers atoms — check anything mapping match indices back to
block order). Vocabulary cost of `D<n>`, same 178 keys, splitting each key's
ligands by degree signature: non-carbon atoms → median 2 variants/key, 196
variant keys keep ≥250 ligands, 28 original keys lose every variant; all atoms →
median 9 (max 68), 263 keep ≥250, 70 lose every variant. Non-carbon preserves
the vocabulary; all-atom fragments it.
Carbon measured separately (2026-09-09, `uncapped_profile`,
`[C;!R][C;!R](=[O;!R])[O;!R]`, bb/ASP/GLU, 4,425 obs): H-bearing boundary carbons
sit ≤3.5 Å from an acceptor O in 15–24% of observations vs 9–11% for no-H carbons
(bb: 12–25% vs 3%), with the O in the H hemisphere 55–68% vs 11–25% — a real
C–H···O signal, and 18% of clusters mix the classes. Heavy degree is the wrong
carbon primitive (D3 conflates sp3 C–H, 24%, with sp2 no-H, 11%; 1/2/3 H don't
order); a binary `H0`/`!H0` captures it at 50 of 178 keys losing every ≥250
variant (vs 28 for non-carbon `D` alone, 74 for full carbon H count). Confirmed
on `frag_lib` (16 keys chosen for balanced H/no-H classes across alkyl-N/O/S,
aromatic c–n/s/o, substituted aromatic, carbonyl-α, ring sp3; 34k obs): on `bb`,
H-bearing > no-H with the bootstrap CI clear of 0 at 40 of 45 slots where both
classes have ≥30 obs; ASP/GLU 28 of 44 (the 5 `bb` misses: three same-sign CIs
grazing 0, one underpowered, one acyl-phosphate vs CH2–O–P class confound; none
is a confident reversal). Pooled over `bb`: CH2 14%, sp3 CH 9%, sp2 no-H 5%,
quaternary sp3 4% within 3.5 Å of backbone O, so the no-H signal is not just C=O. The
earlier alkyl C–O null was underpowered (now +0.07 [0.01, 0.12] on both C–O
carbons). Ligand Cα-like positions are ≥94% H-bearing, so the flag is
single-class there and never hides a Cα–H···O contact. Scripts and tables:
`~/docking/scratch/carbon_degree_test/` (`keys_fraglib_ranked.txt`,
`obs_fraglib_summary.tsv`, `report/carbon_flag_report.pdf` with forest plots).

Formal charge comes from the same place but is pinned by sanitization, not by the
template's drawing: `MolFromPDBBlock` and the template's `MolFromSmiles` both
sanitize (only the final H-strip is `sanitize=False`, which preserves charges).
Verified through `Frags.get_query_ligand_mol`: nitrobenzene reaches the matcher
as `[-1, +1]` from a hypervalent template (`N(=O)=O`), a charge-separated one,
and from no template, matching `c[N+;!R](=[O;!R])[O-;!R]` every time. A new way
of building a query mol must preserve that — library keys were normalized the
same way and nothing checks the two agree.

### Fragment cost estimates are counted from `--pdb-dir`, not from the CCD

`estimate_frag_cost.py` samples the **parent PDB mirror** given as `--pdb-dir`
and matches fragment keys against the ligands found there. It never reads the CCD
or a previous build. A fragment whose ligands are absent from the mirror
estimates 0 structures / 0 occurrences, falls below `--min-instances` in
`make_sge_scripts_for_frags.py`, and gets no job — silently, since 0 is also what
a genuinely rare fragment gets.

Check this when **adding novel ligands** beyond the CCD: their structures must be
in the mirror these scripts are pointed at, and the estimate regenerated, or
their fragments count as if the ligands didn't exist. A fragment known to be
wanted but scoring 0 goes in via `--include`/`--include-file`, which bypasses the
occurrence threshold and applies a resource-tier floor instead.

### Other scope limits

- RMSD alignment blows up combinatorially for large/symmetric CGs (phosphates);
  prefilters help, nothing bounds it.
- BioLiP2 contains pseudo-ligands (`ALA1`, `ALAA`, …) needing cleanup — why
  preprocessing is slow.

---

## Caveats

### Data model and labels

#### CG atom correspondence comes from the automorphism group, not from stored order

Two standing prohibitions

1. **Never replace a consumer's enumeration with a reliance on the stored CG atom
   order.** Every atom-correspondence comparison must minimize over the group in
   `nr_vdgs/cg_symmetry.npz` (`vdg_npz_utils.load_cg_symmetry`); the stored order
   is one arbitrary representative of an orbit, and per-record perception failures
   can even contradict the SMARTS topology.
2. **Never drop "improper" permutations from the group** — ones that map a
   structure onto its mirror image rather than onto itself, e.g. the swap on
   2,3-dimethyloxirane (`CC1OC1C`) that merges cis and trans. They look like
   spurious extra permutations and are not: dropping them makes RMSDs
   orientation-dependent and splits clusters that should merge.

Where the derived group is known to be too coarse, and why it isn't narrowed: see
"Ring positions are merged with exocyclic ones in automorphism orbits" under Known
issues, and `docs/symmetry_edge_cases.md`.

#### A backbone label is a role, not a residue — read `aa_bucket_parts`, never resnames

`bb` means "this slot contacts the CG via backbone"
(`align_and_cluster.reorder_vdg_subset`), covering every donor residue incl.
Gly/Pro. The same LEU is `LEU` in one vdG, `bb` in another.

1. **Slot-swap eligibility comes from `aa_bucket_parts` only** —
   `vdg_npz_utils.aa_perm_indices(bucket["aa_bucket_parts"])`. Resnames fail
   silently both ways: a `bb` slot holding a second ASP invents a spurious swap;
   one holding a non-canonical name (MSE, SEP, TPO, PTR) matches nothing. Fixed
   in `hit_finder_core.py`; `compare_cg_geometries.py` was already correct.
2. A backbone label says nothing about *that slot's own* identity — bucket
   `PRO_bb` names the sidechain-contacting slot; the `bb` slot may itself be PRO.
3. Not a wildcard a resname can be promoted into — `query_slot_labels` decides
   which labels a query may be looked up under, deliberately not onto.
   Unreachable buckets (incl. `X`) are by design; don't add fallbacks.
4. Bucket aggregations disagree on scope: `common.load_bucket_counts` always
   drops `X`, plus backbone buckets under `bb_mode='off'`;
   `compare_cg_geometries._bucket_counts` keeps everything. New aggregations must
   choose explicitly.
5. A backbone label conflates "sidechain unresolved" with "backbone genuinely
   closer" — check `nr_slot_flag`/`slot_reason()` (`SLOT_*` in
   `functions/vdg_struct_utils.py`). `SLOT_NO_SC` doesn't distinguish
   GLY-by-chemistry from unresolved (only `nr_scrr_resname == GLY` does). It's
   the **nr vdG's** flag only — clustering ignores it, so a cluster can pool
   mixed reasons; don't weight by `cluster_size`.

`nr_scrr_resname` always gives the true residue — fine to report/filter, never to
derive permutations from.

#### GLY/PRO backbones are pooled into `bb`, screened at read time

GLY has no CB (54.7% of the backbone pool — a glycine contact may sit where
another residue's CB would clash); PRO has no amide N-H (2.3%, also restricts
phi). Neither is in the label; both are tested per vdG at hit-finding time
against the *query* residue's atoms
(`hit_finder_core.backbone_slot_blockers`/`backbone_slots_can_host`,
`BB_SLOT_SIDECHAIN_CLASH`=3.4 Å, `PRO_NH_DONOR_CUTOFF`=3.5 Å). Sidechain-labeled
slots are deliberately not screened.

Caveats: uses the query's CG coordinates (can differ from the matched nr vdG's by
~1 Å, so a borderline clash can pass); near-no-op on crystal queries, so it only
earns its keep — and can only be validated — on docked/predicted poses (covered
by `tests/test_backbone_slot_hosting.py` only). Calibration:
`docs/representation_assessment.md`.

#### `cg_atoms_dict` is a resname-keyed union, not an indexable match list

In `VDG.structure_contacts` (`vdg.py`), the non-proteinaceous branch unions the
match lists of **every copy** of a resname.

- **Positions are meaningless** — entries from different copies are concatenated,
  so index *k* isn't `cg_idx` *k*. Never use it to resolve a CG site; the
  per-copy lookup keyed on `(struct, seg, chain, resnum, resname)` is in the
  `cg_idxs` loop below. The proteinaceous branch builds a genuinely indexable
  dict — same variable, two shapes.
- **The union is required for correctness** — it feeds `preprocess_lines`'s probe
  prefilter, which ORs over `(resname, atomname)`. A last-wins dict drops probe
  lines for atoms the last copy lacks (disorder/obabel perception differ by
  chain). In 4a05, `BGC`'s last-wins copy covered 11 of 12 atoms across 5 copies,
  masking every probe line for `O1`. Across 15 fragments, 2.4% of multi-copy
  groups had an under-covering copy.

Over-inclusion is safe: the prefilter only decides which probe lines are
*considered*; the per-copy lookup still assigns each contact to a site.

#### The occupancy protocol is duplicated in the vdG-miner submodule

#Atom roles in a vdG PDB are encoded in occupancy: CG atom *i* at
#`3.0 + 0.01 * i`, vdM residues `2.0`, non-CG ligand atoms `1.0`. `3.00-3.99` is
#closed, `4.0`+ reserved (constants in `functions/vdg_struct_utils.py`).
#
#`external/vdG-miner` isn't installable, so its reader `fingerprint_helpers.py`
#hardcodes its own copy of these constants (the writer,
#`clus_and_deduplicate_vdgs.py`, imports them) — changing the encoding means
#editing the submodule in the same commit. The 100-slot cap is a real CG-size
#limit (`cg_slot_occupancy` raises past it; ~20x headroom vs. the largest fragment
#at 5 atoms). Readers use only the *order* and distinctness of CG occupancies,
#never the value.

#### Chain IDs are one column; two-character chains are remapped, never read as-is

Every reader in the pipeline takes the chain from PDB column 22 alone (ProDy,
probe, vdG-miner's `line[12:26]` hash and its HETATM ligand scan). Source files
for large assemblies write two-character chain IDs across columns 21-22
(`ASPA4  66` = chain `A4`), so chain `A4` reads as chain `4` and residue 66 of
both chains merges into one residue. ProDy 2.6 gives the merged atoms a single
resindex, and `writePDB` round-trips the collision with column 21 blanked, so
parsing before fixing hides it. In the 2026-09 database 670 of 65,600 structures
had such chains and vdG-miner dropped all of them whole (probe and PDB lines
hash differently). s01 fixes the text before parsing
(`preprocessing/_chain_ids.remap_chain_ids`: first free character from
`CHAIN_POOL`, mapping recorded as `REMARK 900 CHAIN ID REMAPPED` lines) and
skips structures with more chains than characters. The same source segments are
written as HETATM, so `embedded_amino_acid_hetatm_to_atom` runs in the same
pass: otherwise the remap alone would let the miner scan protein residues as
ligands. Don't compare library chain IDs for these structures against the
deposited entry without the REMARK mapping.

#### Output PDB names parse from the right, not the left

`materialize_vdg_pdbs.py`/`write_vdg_hit_pdbs.py` join name fields with a single
`_` (`functions/vdg_pdb_io.py`). The *head* has a variable field count (`<frag>`
is a sanitized SMILES; `<AA_BUCKET>` joins labels with `_`, so `ASP` is one field
and `ASP_bb` two; `<source>` may carry an assembly suffix like `1f8s_1`) — don't
index from the front.

Only the tail is fixed: exactly `1 + subset_size` residue tags of 4 fields each
(`seg_chain_resnum_resname`), then the ligand. Take the last
`4 × (subset_size + 1)` fields after stripping the extension and any `~<n>`
duplicate suffix. Hit files from `write_vdg_hit_pdbs.py` differ twice: no ligand
tag, and the name ends `..._<vdM tags>_<rmsd>` — drop one field before counting
from the right.

The library's `.npz` bucket file names are unrelated to this grammar.

#### Fragment directory names are glob patterns — `glob.glob` silently returns nothing

An annotated fragment key is also the directory name under `frag_lib/`, and
`[C;!R]`/`[N;r5]` are valid glob bracket expressions. So
`glob.glob("frag_lib/[C;!R][C;!R][C;!R][N;r5]/nr_vdgs/1/*.npz")` returns `[]`
while `os.listdir` on the same path returns 20 files — the pattern is
interpreted, matches no 4-character name, and the empty result is
indistinguishable from an unbuilt fragment. Aromatic keys (`ccnnn`) have no
metacharacters, so such a script appears to work on part of the library.

Use `os.listdir`/`os.scandir`, or `glob.escape(path)`. Nothing in the repo globs
fragment names today, so this is a rule for new scripts, not a known bug. In the
shell, always quote: `ls "frag_lib/$frag"` — unquoted works only because bash
leaves non-matching patterns literal and no fragment name happens to match
another (checked across the 81 present).

### Geometry and alignment

#### Degenerate Kabsch fits

`utils.kabsch` always returns a proper rotation (`det(R)=+1`). Rank ≥2 (≥3
non-collinear points) → unique rotation; rank 0/1 → SSD still meaningful but
rotation isn't. Fine if you only use SSD; not if you apply `R` to points outside
the fitted set.

`vdg_npz_utils._append_full_ligand_from_parent` skips full-ligand reconstruction
when either CG point set's secondary/primary singular-value ratio ≤ `1e-3`, to
avoid rotating unmatched atoms arbitrarily about a near-collinear axis.

#### Chain breaks are masked, not penalized

CA–CA > 4.5 Å sets the flanking residue and everything past it to
`FLANK_CHAIN_BREAK` (`'!'`, NaN); a missing/non-protein/ambiguous/no-CA residue
is `FLANK_MISSING` (`'-'`) — only positions *past* a break become
`FLANK_CHAIN_BREAK`. Positive/negative flanks are scanned independently. Neither
marker is `'X'` (`NONCANONICAL_AA_LABEL`, a vdM slot label).

Stage-2 RMSD clustering drops non-finite rows rather than penalizing them
(denominator is the shared-finite-row count, but the cutoff uses the full
expected array size), so a terminal vdG with missing flanks can cluster with an
interior one on core geometry alone. Sequence-coverage similarity gives both
markers the configured threshold as a neutral prior, so a missing residue and a
chain break score the same even though the data distinguishes them.

#### `EXCELLENT_MATCH_CUTOFF` reports the first good match, not the best one

`hit_finder_core.EXCELLENT_MATCH_CUTOFF = 0.3` Å — once a permutation hits RMSD ≤
cutoff the loop `break`s, so `best_*` is the first sub-0.3 Å result, not the
minimum over permutations.

**Intended, not a correctness problem**: the break fires only after a hit is
established, so no hit is lost; only the reported value can be slightly
pessimistic, and 0.3 Å over a full backbone+CG motif is inside crystallographic
error. Matters for small upward bias when ranking/thresholding on `vdg_rmsd`, and
`write_vdg_hit_pdbs.py` may write the second-best superposition. For the true
minimum, set the cutoff to `0.0` rather than deleting the branch.

### Counting and statistics

#### Overlapping CG sites credit a shared contact to every site

A ligand can match one CG's SMARTS several overlapping ways — SAH's ribose
matches `OCCO` three ways; UPA in 11ba, eight. 83% of multi-match CG sites
(745,207/894,974) overlap, worst on `CCOC`, `OCCO`, `OCCCO`, `CC(O)CO`, `COCCO`.

`structure_contacts` emits one neighbor row per site, so the same physical
contact appears in more than one vdG: observations across overlapping sites are
**correlated, not independent** (inflating propensities/cluster sizes), and the
inflation is **non-uniform**, scaling with a fragment's SMARTS self-overlap.

Deliberate: crediting only one site truncated the other sites' environments and,
when a shared atom carried a site's only contact, dropped that vdG entirely (SAH
`O4' C4' C3' O3'`: two residues genuinely contacting it, no vdG). Incomplete
geometry can't be recovered downstream; correlated counts can be normalized
later, since parent biounit and per-record CG atom names are retained
(`nr_parent_biounit`/`nr_cg_names`, `mem_*`).

#### Every chain copy of a CG is mined, so NCS copies inflate cluster sizes

`update_sc_info` builds one `sc_info` entry per `(struct, segi, chain)` carrying
a CG match, and `mine_environments` iterates all of them (it used to keep only
the best-contacted copy, discarding ~31% of chain-level CG copies without
checking geometric identity — redundancy is now measured downstream, in
clustering).

Cost: a homo-octamer contributes up to 8 near-identical members to one cluster,
so `cluster_size` reflects crystallographic multiplicity as much as interaction
frequency, non-uniformly across fragments. Count **distinct parent structures**
instead (`nr_parent_biounit`/`mem_parent_biounit`). Compounds with the
overlapping-site double-counting above. Partly counteracted upstream: s01
trimming drops intra-PDB copies with matching contact fingerprints, so NCS-copy
count depends on whether the library was built from a trimmed database.

Related fix: a per-chain "no selectable CG residues" case used to `return` and
abandon the whole structure; `mine_environments` now skips just that chain.

#### Hit deduplication is broader than "cross-library"

`deduplicate_hits` keys on (BSR residue *set*, ≥`min_shared_atoms` shared ligand
atoms), ignoring `frag`, `aa_bucket`, and nr vdG — so it also collapses multiple
nr vdGs within one bucket and multiple `aa_bucket` variants of one fragment,
leaving at most one hit per (residue set, ligand site). Counts in
`results_summary.txt` are post-dedup; `--no-dedup` gives raw counts.

The key is the residue *set*, not the `bsr_combo` string: that string is ordered
to match `aa_bucket`, and backbone labels sort after every uppercase resname, so
one physical residue pair emits several orderings across its `aa_bucket`
variants. Grouping on the raw string — including `--no-dedup` output — splits
those variants apart.

#### `min_dist` is not comparable across buckets

In `compare_cg_geometries.py`, `min_dist` is the smallest CG-nr vdG distance over
all N_A × N_B pairs — a minimum over a growing sample can only decrease, so it
tracks bucket population as much as similarity. On `CC(=O)O` vs `CC(=O)[O-]`,
bucket `ARG`: 0.074 Å at 25 sampled nr vdGs, 0.051 Å at 50, 0.0 Å at 100+.
Compare only between buckets of similar N, or use a sample-size-stable statistic
(e.g. 5th percentile of pair distances). The CG nr vdG is also inside the Kabsch
fit whose residual is then measured on it, deflating the value further.

#### Other `compare_cg_geometries.py` columns

- `frac_A_matched`/`frac_B_matched` are computed on the `--max-per-lib`
  subsample (default 1000) but read as full-library fractions.
- `enrichment_A`/`enrichment_B` describe only the **first** residue of a
  multi-residue bucket — for `ASP_HIS` you get ASP; HIS is dropped.
- Backbone buckets are counted like any other here — this script uses its own
  `_bucket_counts`, not `common.load_bucket_counts`, so the propensity path's
  `X`/backbone exclusions don't apply.
- The heatmap's `vmax=3.0` is hardcoded, unrelated to `--match-threshold`.

#### Enrichment statistics

- AAs with an observed count of exactly zero are **omitted** rather than recorded
  as depleted, with no pseudocount (`identify_bioisosteres/common.py`) — a
  fully-depleted AA is indistinguishable from an unsampled one.
- Pair mode (`compare_aa_profiles.py`) correlates the upper triangle of a
  symmetric matrix; entries are non-independent, so p-values are
  anticonservative, with no multiple-testing correction across O(n²) CG pairs.
- `--jaccard-positive-only` masks entries positive in **both** profiles,
  inflating the score for AAs enriched in only one CG.
- Enrichment is relative to PDB background frequency (crystallized, not intrinsic
  affinity); size normalization uses heavy-atom count, cruder than SASA.
- **Profiles built under different `--bb-mode` are not comparable and align
  silently** — category sets overlap, so `compare_aa_profiles.py` correlates them
  with only denominators differing. It warns when the NPZs' recorded `bb_mode`
  disagree — heed it.
- **`GLY-bb` is the only backbone category typically enriched** (median +0.45
  over 303 CGs, positive for 94%; every other `<AA>-bb` has negative median,
  positive ≤6%), mostly reflecting glycine's backbone being the most exposed
  (~49% of all backbone clusters). Not dominant (median rank ~13/39) — compare
  backbone rows across CGs, not against zero.
- **Backbone categories thin out fast under `--bb-mode per-residue`** (~278
  backbone clusters per CG split 20 ways; `GLY-bb` ~135, `TRP-bb` ~3) — rare
  `<AA>-bb` categories are noise and a zero count disappears entirely. Why
  `--pair-bb-mode` defaults to `pooled`.
- **Supplying `counts=` yourself changes the denominator** — both propensity
  functions sum `counts.values()` but skip categories without a usable prior, so
  mixing `bb_mode` between the count dict and the scoring call drops categories
  from the numerator while keeping them in the denominator. Pass the same
  `bb_mode` to both.

#### `cg_rmsd` in `benchmark_pose_recovery.py` is not a prediction

`score_one_model` fits `R, t` to a target containing backbone **and** the query
CG (`Y[idx, N_bb:] = q_cg_coords`), then measures `cg_rmsd` against the crystal
CG using that same transform — part of what was fitted. It's the CG-only
component of `vdg_rmsd`, bounded by it, not an independent recovery measurement.
Correct for scoring a supplied pose; wrong for docking/prediction, where the fit
must use backbone only — nothing here measures that version.

#### `plot_from_tsv.py` pools scores that are not on a common scale

`condition_grid_points` pools `(ligand RMSD, score)` across structures for the
panel's `pooled ρ`/`r`. A raw score isn't comparable between structures — it's a
weighted hit count that grows with ligand size and pocket richness, and
`--normalize-by-frag-avg` divides by *that structure's own* per-fragment mean, so
`1.0` means different things per structure. Pooled correlation thus measures
partly the real score/RMSD relationship and partly which structures sit high or
low.

The fix — per-structure ranks or z-scores before pooling — is not implemented.
Until then, **mean per-structure Spearman** (what `summarize_condition` prints
and each panel title leads with) is the well-defined statistic; the pooled value
shows spread, not the headline.

Related: `score_samples` uses the structure's full sample roster, so zero-hit
samples sit in `frag_avg`'s denominator — omitting them would scale every score
by a different per-structure constant, invisible to per-structure correlation but
landing on the pooled value.

### Database preparation

#### Database trimming counts protein neighbours only

`interactions.get_nr_res_interactions_with_ligs_in_pdb` builds each ligand's
fingerprint from `protein and exwithin` (disjoint from the ligand selection),
replacing a plain `within` that also swept in the ligand's own residue and nearby
waters/hetero atoms. Consequences: networks are smaller, so
`redundancy.check_networks`'s fixed tolerances (intra-PDB dedup `tol=2`) are
relatively more permissive; a ligand whose only neighbours are non-protein is
dropped outright (nucleic acids and ProDy-missed modified residues land in the
*ligand* selection themselves). Expect the trimmed database to shrink,
concentrated in nucleic-acid structures.

Smaller fingerprints also make intra-PDB dedup catch more chain copies — the
upstream half of "Every chain copy of a CG is mined"; read both before drawing
conclusions from `cluster_size`.

#### Prepwizard relabels atoms it cannot build as part of a ligand

An amino acid modeled with only its backbone N can't be built as an amino acid by
prepwizard (s02); it reclassifies the orphan N as a free molecule, protonates it
to ammonium, and writes it under the **resname, chain, resnum of a ligand** —
coordinates unchanged, so the ligand gets a stray atom tens of Å away (`2y1x` SAH
A:1001's second `N` is really THR D:478's, 54 Å off).

`_prep_filters.drop_prepwizard_hazard_residues` removes the cause (runs from s01
and s02), so a database prepared with current code shouldn't show this. Two
things to know anyway: the name collision is incidental and is what the guard
keys on (599/65,600 prepared PDBs have a ligand atom >5 Å from every other heavy
atom in its residue; only a handful collide with a real CG atom name, so the
duplicate-name warning in `_get_atomgroup_for_env` undercounts the defect);
single-atom HET groups are affected the same way and dropped too (589 lone
backbone `N`, 11 lone `C1` from one-atom HET residues like `CF0`/`0QE`, one
lutetium ion that hosted a lone N — `single_carbon_het_resindices` is restricted
to carbon so monatomic ions survive).

Prepwizard also **renames some HET residues into `ATOM` records under a standard
amino acid resname**; `find_cg_matches` reads only `HETATM`, so an affected
ligand disappears from mining with no warning (`CYT`→`CYS` in `5buv`/`5epu`,
`HSE`→`HIS` in `6a0s`, `SRO`→`SER` in `7bs2`; 37 per 921,702 residues). s02 undoes
this via `snapshot_restorable_ligands`/`restore_renamed_ligands`, matched by the
literal 24-character coordinate field (`_prep_filters._coord_key`, not
chain/resnum, so it survives renumbering) — relying on `-noimpref` not moving
heavy atoms. A mismatch **skips a restore rather than corrupting anything**. If
you change prepwizard flags, re-check `restore_renamed_ligands` still reports
non-zero on `5buv`.

Only *free* ligands are restored — a residue covalently bonded to a standard
amino acid is left as prepwizard wrote it (`CR8` in `3tmr`, `MDO` in `2o7d`),
since reverting a chain-embedded chromophore to HETATM would make it mineable as
a ligand (a composition change, not a fix). Renamed *atoms* (`6a0s` `HSE`:
`N`→`NA`) aren't restored — matching and per-name selection read the same file,
so internal consistency suffices.

#### Modified residues are protein records, decided by CCD type, not by atom names

Without prepwizard, a modified residue stays `HETATM` under its own resname
(`KCX`, `SEP`, `LLP`, `MSE`). vdG-miner then rejects it as a slot (20-resname
whitelist in `vdg.py`) *and* scans it as a ligand (`cg.py` reads `HETATM` only),
so a phospho-serine would be mined as a phosphate CG. s01 therefore runs
`_prep_filters.modified_residues_to_protein` on the text before parsing:

- `MSE`→`MET`, `SEC`→`CYS` (`MODIFIED_RESIDUE_RENAMES`), atom names unchanged;
  `SE` is admitted by `_ALTERNATE_HEAVY_ATOMS`, so these are full slots. In the
  2026-09 database MSE is 4881 residues, 1148 near a ligand.
- Every other `HETATM` residue whose CCD `_chem_comp.type` ends in
  `PEPTIDE LINKING` and has a heavy atom within 1.8 Å of a standard amino acid
  becomes an `ATOM` record with its own resname: neither slot nor ligand.
- `NON-POLYMER` residues bonded to the chain (PLP, HEM, HEC, FAD, SAM) stay
  `HETATM` and are mined as ligands, as today.

**Backbone atom names are not the test.** 353 chain-bonded cofactor residues in
the audited database (SAM, SAH, NXL, 0G6) carry `N`/`CA`/`C`, and 23
peptide-linking residues do not. The CCD type table is
`resources/ccd_polymer_types.tsv` (`scripts/fetch_ccd_polymer_types.py`, login
node only: compute nodes have no network). Resnames absent from the table are
left as written and listed once at the end of s01; a chain-bonded modified amino
acid among them would be mined as a ligand, so refresh the table if the list is
not empty. Known small loss: a peptide-linking ncAA bonded to a *standard*
residue of a peptide ligand is converted too (85 chain-bonded peptide-linking
residues in 65,600 structures, 15 near a ligand, 10 of them `SEC`).

Prepwizard-era databases differ: there the modification sits on a parent-named
residue (`LYS` with PLP atoms), which is what the `X` slot label was built for.
`X` stays as the safety net for anything either path misses.

### Job generation

#### `make_sge_scripts_for_hit_finder.py` treats `--ref-pdb` as global

`resources/validation_set.csv` is `(num_procs, query_dir, smiles)`; any
`--ref-pdb` passed to the generator is copied into every generated job. Fine for
one reference across the batch; needs a 4th CSV column and generator support if
different query directories need different references.

---

## Known issues

Settled. Read to avoid re-opening a closed question; no action expected.

### vdG automorphisms and hit-finder `CalcRMS` use different policies, deliberately

vdG generation calls `utils.identify_mol_automorphisms`, applying this project's
resonance/aromatic normalization to unsanitized **fragment** SMARTS. The hit
finder's `utils.best_inplace_symmetry_rmsd` instead uses RDKit's
`CalcRMS(maxMatches=0, symmetrizeConjugatedTerminalGroups=True)`, which doesn't
consume those automorphisms. **Making them agree would be wrong** — the fragment
normalization exists because fragments are truncated (a degree-1 atom may have an
omitted parent bond; see next entry), while `best_inplace_symmetry_rmsd` runs on
**full ligand graphs**, where a terminal oxygen really is terminal. Nothing keyed
on it enters the library — the only caller is `hit_finder_core.score_one_model`,
computing a reported pose-vs-reference metric (`lig_rmsd_value`) used only in
`rmsd_records` and a `plot_from_tsv.py` threshold.

### Fragment symmetry lacks parent-boundary information

Automorphisms for vdG generation come from the fragment SMARTS alone. An atom
with degree one in the truncated graph can still have an omitted bond in the
parent — `CCP(=O)O` doesn't say whether both oxygens are genuinely terminal. A
dictionary-wide P-O permutation can be right for one occurrence and wrong for
another with a different omitted attachment. The terminal-atom rule applies
uniformly to terminal N, O, S on B/C/N/O/P/S/Cl/Br/I centers.

### Ring positions are merged with exocyclic ones in automorphism orbits

`utils._automorphism_graph` carries no ring-membership label, so in keys like
`cn(c)c` and `nc(n)n` a ring nitrogen and an exocyclic one land in the same
orbit — the derived group permutes positions that are not chemically
interchangeable.

**Settled: no fix.** Adding the label was considered and rejected: each of these
keys mixes a bridgehead-nitrogen chemotype (no exocyclic atom at all) with an
N-substituted one, so any group restriction that helps the minority breaks the
majority. There is no speedup on offer either way — the batched Kabsch cost is
sublinear in permutation count.

### Fragments inherit the parent's perception and are never sanitized

The parent ligand is sanitized (`fragment_database_ligs.main`), but cut fragments
are not: `Frags.fragment_on_bond_d` returns raw `Chem.PathToSubmol` output, so a
fragment's canonical SMILES — its identity everywhere in the pipeline — encodes
the *parent's* perception. Since 2026-09-04 the key also carries each atom's ring
context *in the parent* as a SMARTS primitive (`[C;!R]`, `[C;r5]`, `[C;r6]`, …,
never plain `R`), recorded before `PathToSubmol` cuts the ring
(`Frags.ring_query_for_atom`).

**Settled: sanitizing fragments was measured and is strictly worse.** Over 25,424
fragments from 1,500 CCD ligands, `SanitizeMol` fails outright on 52.7%
(ring-opening leaves aromatic-flagged atoms with no ring to justify the flag);
the 47.3% that survive get re-described as a free-standing species — bracket-form
cut atoms gain radical electrons (`CC(N)=O`→`[C]C(N)=O`), non-bracket ones
silently gain implicit hydrogens. Either way the fragment stops describing a
substructure. Inheriting the parent's flags is wanted (aromatic `c` in a key
means "aromatic in the ligand it came from"); partial sanitization isn't worth
trying either, since the dominant failure is the ring-opening cut itself.

**Downstream cost: measured at zero.** Over the full 48,709-ligand CCD
(2026-08-27): 0 SMARTS re-parse failures, 0 substructure-comparison errors, 0
empty fragment-site groups. A future non-zero count is this known issue, not
corruption — failed fragments are dropped by `manually_remove_Hs`; a failed
comparison falls back to "not a duplicate" (can double-record under two keys,
never loses data). Run totals: 47,226/48,709 ligands parsed, 173 failed, 146
failed sanitization, 1,164 not druglike, 2,196 unique fragments.

**The one real action: parse fragment keys with `utils.mol_from_fragment`, never
`MolFromSmiles`** — keys are stable only under the unsanitized path that produced
them. Already the convention everywhere (`clus_and_deduplicate_vdgs`,
`hit_finder_core`, `dock_utils`, `extract_fragment_smiles`,
`vdg_generation_wrapper`, `scripts/` audits).

Sanitizing the *parent* (added 2026-08-27) shifts a minority of keys: on a
3,000-ligand sample, 765/~810 unchanged, 44 disappeared (un-perceived parsing
artifacts), 49 appeared (their correctly aromatized forms).

Related — see "Ligand bond perception can silently swap atom identities."

### Fragment keys match exactly on aromaticity, and nothing bridges a mismatch

A fragment key is a SMARTS, so aromatic `c` matches only an aromatic query atom,
and hit finding matches those keys against the query ligand directly
(`match_library_frags_to_query`). A key whose aromaticity disagrees with the
query mol's simply finds no match. (Before 2026-09-05 the same was true via an
isomorphism fallback on two `MolFromSmarts` mols, which covered canonical-*string*
differences only. Key resolution by name, and that fallback, are gone: the key
set is now read off the library.)

**Settled: nothing to do, the two paths agree.** Kekulé- and aromatic-written
parents give identical fragment key sets; on the query path
(`MolFromPDBBlock`→`AssignBondOrdersFromTemplate`→`manually_remove_Hs`) a
two-heterocycle amide comes out fully aromatic with a key set identical to the
library path. If you introduce a third way of building a query mol, verify its
key set against the library path first — a mismatch shows as fragments silently
absent, not an error.

### Ligand bond perception can silently swap atom identities

`Frags.get_query_ligand_mol` writes the query ligand to a PDB block parsed with
`Chem.MolFromPDBBlock` — no bond records, so RDKit infers connectivity from
distances. `AllChem.AssignBondOrdersFromTemplate` raises only when the inferred
graph can't match the template *at all*; a wrong-but-isomorphic graph is
accepted, attaching template chemistry to the wrong physical atoms.

**Scope is query-side only** (`vdg_hit_finder.py`/`hit_finder_core.py`). Library
generation fragments CCD SMILES directly, where connectivity comes from the
SMILES — a mis-bonded query can produce a wrong hit but can't corrupt the
library. A swap between automorphism-related atoms is harmless; the damaging case
is a coincidental isomorphism mapping chemistry onto genuinely different atoms.

**A second, more common mechanism: ambiguous template matches.** With no bond
orders, like-element atoms with the same connectivity are interchangeable; RDKit
prints `More than one matching pattern found - picking one` and picks by atom
index, ignoring geometry. Demonstrated on acetic acid (hydroxyl oxygen written
first): the double bond goes to the atom at 1.34 Å, the single bond to the one at
1.22 Å. Real cases: carboxylic acids, esters, carbamates, amidines, phosphate
mono/diesters, N-oxides — anywhere a like-element pair member carries an H the
file doesn't have. Genuinely-automorphic pairs (carboxylate, nitro, sulfonamide)
make the swap harmless.

**Settled: not being fixed.** Full detection needs ground-truth connectivity
keyed by PDB atom name (CCD `_chem_comp_bond` — a ~500 MB dependency plus
plumbing), not worth it against a query-side-only failure that can't reach the
library; cheaper options (bond-length checks, disambiguating by distance) flag
suspects without resolving them. **Unquantified by choice** — the ceiling is a
wrong hit on one query pose, a limit on trust rather than queued work.

### Protonation variants are collapsed, so charge resolves one way only

Fragment selection pools protonation-state variants onto one representative
before thresholding (`group_protonation_variants` in
`extract_fragment_smiles.py`), and only that representative is mined — the
charge-loose SMARTS `cC(=O)O` matches both benzoate and benzoic acid in the
prepwizard-protonated mirror, so the anion's library would be a strict subset.
The alias map goes to `fragment_aliases.tsv` (library root) or
`<work list stem>_aliases.tsv`.

**A representative need not be a frags-dict key.** When two or more charged
spellings share a charge-stripped form and the dict holds no neutral twin
(aromatic nitro: RDKit writes `c[N+](=O)[O-]` or `c[N+](=O)O`, never neutral),
the charge-stripped SMARTS itself is *promoted* to representative and names the
library directory. Looking that key up in `database_frags_dict.pkl` by string
finds nothing, and no CCD ligand is drawn that way; that is expected, not a
missing fragment. The alias file marks these rows `kind=promoted`, and
`scripts/lookup_fragment_key.py KEY` lists the dict keys a key covers (exact,
reordered, or charge-variant). On the 2026-09 dict there are 11 promoted groups
(nitro, azides, isonitrile, nitrone, N-oxide); only nitro clears the threshold.
A lone charged spelling with no twin is *not* promoted: nothing would pool, and
the key that exists in the dict is the more useful directory name. The
2026-09-06 library predates this rule and holds both nitro spellings as two
nested libraries.

**Settled: the query side reaches the collapsed library without the table.** Hit
finding no longer derives a query key to look up — it matches the library's own
keys against the ligand, and the neutral key's `[O;!R]` matches an anionic
oxygen, so a carboxylate ligand reaches the acid directory. (The earlier
2026-09-04 fix resolved the query key through the cached alias table instead;
either was needed because the isomorphism fallback couldn't bridge the gap —
`[O;!R]` and `[O-;!R]` don't mutually match, and `GetFormalCharge` on a SMARTS
mol reports 0 for both, so the charge test passed vacuously.)

Charge therefore resolves one way only: a key that is itself charged doesn't
match a query drawn neutral. Currently empty in practice — at the production
threshold (`--min-instances 250`, `--max-size 5`) 1 of 177 representatives
carries a formal charge, at the looser 100 cut 2 of 324, and no neutralized key
occurs in `database_frags_dict.pkl`. Neither has a realizable neutral form: a
4-substituent N is always +1, and sanitization normalizes nitro to the
charge-separated form, so `N(=O)=O` and `[N+](=O)[O-]` both match. It becomes
real only if a future key set admits a charged key with a realizable neutral
twin, which `hit_finder_core._warn_charge_only_miss` reports at read time — it
fires only when the neutralized key matches the ligand (i.e. a hit was actually
lost) and says whether a usable neutral twin exists.

Do **not** extend pooling to chemically distinct forms — a donor and an acceptor
should stay split (see the tautomer assumption above). The general form, a
directed "query may draw on" relation covering protonation/H-count but not ring
context, is a TODO item.

Related: the `min_dist` caveat's worked example (`CC(=O)O` vs `CC(=O)[O-]`)
involves nested rather than independent libraries — shared observations are a
second cause of `min_dist → 0`, alongside sample size.

### RDKit and OpenBabel disagree on what `r<n>` means

RDKit reads SMARTS `r<n>` as "the atom's *smallest* ring has size n"; OpenBabel
reads it as "the atom is a member of *some* n-ring". Fragment keys are written by
RDKit (`Frags.ring_query_for_atom`) and matched against query ligands by RDKit,
but matched against database ligands by OpenBabel in the miner
(`external/vdG-miner/vdg_miner/vdg/cg.py`). Uncorrected, the miner accepts CG
instances hit finding can never reproduce.

On the shipped fragment dict, 4292 of 6144 keys carry `r<n>` and the two toolkits
differ for 3.2% of them, in both directions. The divergence needs an atom in two
rings of *different* sizes — spiro/fused/bridged saturated systems, e.g. on
1-oxaspiro[4.5]decane the spiro carbon is `r5` to RDKit and both `r5` and `r6` to
OpenBabel.

Settled: `cg.find_cg_matches` and `estimate_frag_cost._count_one` re-filter
OpenBabel's matches through `cg.match_satisfies_ring_sizes`, which requires
`OBAtom::MemberOfRingSize()` (the smallest ring, matching RDKit) to equal the
annotated size at every constrained position. The constraint list is parsed by
`cg.ring_size_constraints` and its length checked against
`OBSmartsPattern.NumAtoms()`, so an unparsed SMARTS construct raises rather than
silently disabling the filter. Don't "simplify" this back to a raw
`GetUMapList()`.

### OpenBabel does not perceive large macrocycles, so large-`r<n>` keys mine nothing

The `r<n>` annotation is uncapped (`r12`, `r72` occur), but OpenBabel's ring
perception doesn't report membership in rings that large, so such a key matches
nowhere on the mining side even where RDKit matches it on the CCD SMILES it came
from. ~3% of the 1338 keys with n>8 are in this state — the fragment is selected,
a job is scheduled, it runs, and it produces an empty library, indistinguishable
from a genuinely rare fragment.

Not fixed. The `r<n>` filter above doesn't help: it only removes matches, and
there are none. Options if it becomes a problem: cap the annotation above some
ring size, or drop large-ring fragments at selection time. A plain-`R` fallback
is deliberately excluded (it would nest libraries), so "just use `R` for big
rings" isn't available.

### ProDy caps some string fields regardless of the dtype written

`vdg_npz_utils` no longer narrows `name`/`resname`/`chid` widths on the way into
an AtomGroup, but ProDy's own field widths still apply on write-back (segnames
`U6`, elements `U2`), so a segname >6 chars silently truncates regardless of
dtype. **Settled: not fixable at the npz layer** — truncation happens inside
ProDy's field definitions.

---

## Open bugs

### `group_lig_sites_by_overlap` merges disjoint sites transitively

`Frags.group_lig_sites_by_overlap` groups matches with union-find, so A-B and B-C
overlapping puts A and C in one group even when A and C are disjoint
(single-linkage). Verified on `OCCOCCOCCOCCO`: six matches of
`[C;!R][O;!R][C;!R][C;!R][O;!R]`, spanning atoms 0-4 through 8-12, collapse into
one group. `hit_finder_core.py:560` uses the group's list index as `q_site_idx`
and keeps one best hit per site, so two genuinely distinct binding sites on a
polyether or polysaccharide are reported as one.

Fix is complete-linkage (require overlap with every current member, not any one);
the loop is already O(n²) so it costs nothing. Not applied because it changes
`q_site_idx` semantics for hit finding.

Otherwise, code doing something other than intended goes here; a settled
limitation belongs under Known issues above.
