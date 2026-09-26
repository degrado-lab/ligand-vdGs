# Pitfalls

Standing contracts and known silent-failure modes. Key entries by function/file, not
line number; re-check after refactors.

## Standing contracts

- Contact-strength arrays (`buried_area`, `shared_area`, `n_atom_pairs`,
  `min_heavy_dist`) are residue-aligned with `env[1:]`; re-key by
  `(seg, chain, resnum)`, never by position. `shared_area` is a 1/k credit,
  including the ligand residue, so buried plus shared area counts each occluded
  surface point once. `buried_area == 0` does not mean “no contact.”
- `smarts_to_cgs.py` always writes the sibling `<cg>_matches.annot.pkl`; staging or
  cleaning `*_matches.pkl` must carry it. Templated OBMols can contain phantom atoms;
  discard matches touching indices absent from `pdb_atom_names(mol)`.
- `_annotations_as_isotope` and `_query_primitive` intentionally use different
  grammars; unknown annotation tokens must raise, while query primitives may return
  `None`. Keep `tests/test_key_primitive_readers_agree.py` meaningful.
- Ligand atom names live in two vocabularies; never mix them.
  - Deposited PDB names (`pdb_atom_names`, `cg_match_dict`, clusterer selections):
    use only against the same file.
  - CCD template names (roster `types`, fragment-dict name sets): use only against
    each other.
  - Cross only via `ccd_templates.resolve_names`. An alt (legacy) name can be ANOTHER
    atom's current name (ABU legacy `CD` = current `C`), so a hand-rolled lookup
    silently picks the wrong atom.
  - Residual: legacy names that permute current names resolve wrongly by name alone
    (B3M/B3K, 6 roster instances); only geometry can catch it.
- Compare time-varying npz schemas (legacy term) after parsing JSON and removing `build_date`.
  Identify parent databases by `functions/db_identity.identity_of`, never by path.
- `cg_formal_charge == -1` is ambiguous. Treat `cg_heavy_degree == -1` as the
  unreadable sentinel; unreadable CGs belong in `unreadable`, not `neg`.
- `load_vdg_bucket`/`vdg_npz_path` require the fourth positional `sign` argument.
  Flat (unpartitioned) schema-v2 libraries must fail with `BucketSchemaMismatch`, not look
  like empty buckets.

## Assumptions

- A query structure has exactly one ligand resname. Multiple resnames select one by
  `n_heavy * n_residues`, or return an error for that model. A ligand is one HETATM
  residue: cross-residue polymers and covalent protein links are unsupported and are
  not represented by the library.
- Query tautomer and bond orders come from the CCD template, not density. Tautomers
  are not enumerated. H-stripping intentionally merges aromatic N roles; N/O degree
  stays out of the key. H status is the per-observation `cg_placed_h` label.
  Formal charge must survive any new query-molecule construction path.
- Fragment cost estimates sample `--pdb-dir`, not the CCD. Novel ligands must be in
  that parent database and the estimate must be regenerated; use `--include` for a
  wanted fragment with zero observed occurrences.
- Large/symmetric CGs can make RMSD alignment combinatorial. BioLiP2 includes
  pseudo-ligands requiring preprocessing cleanup.

## Caveats

### Representation and parsing

- `prep_benchmark_set.py` checks CCD codes against the current RCSB CCD before treating a
  missing roster code as novel. Obsolete source codes (for example PLINDER `HSR` and
  `ZN2`) can be absent from current coordinate files; map a code to a current ligand
  and verify the coordinates before considering it for a benchmark.
- `prep_benchmark_set.py` uses `altloc="all"` only to discover ligand altloc labels.
  It retains separate coordinates for the same atom in each altloc, so do not write
  that combined selection as one ligand. Some PanDDA entries have the bound ligand
  only in altloc B; prep reparses a single altloc before writing the complex.
- Atom correspondence must minimize over `nr_vdgs/cg_symmetry.npz`, including improper
  permutations; stored atom order is not authoritative. See
  `docs/symmetry_edge_cases.md` for known over-broad ring/exocyclic orbits.
- `aa_bucket_parts`, not residue names, defines slot permutability. `bb` is a role,
  not a residue or wildcard; use `nr_slot_flag`/`slot_reason()` to distinguish why it
  was assigned. `nr_scrr_resname` is for reporting/filtering only. GLY/PRO backbone
  compatibility is screened at hit-finding time against the query.
- `cg_atoms_dict` is a resname-keyed union across copies, not an indexable match list;
  resolve sites through the per-copy key. Occupancy encodes vdG atom roles; readers
  use order/distinctness, and the 100-slot CG cap is real.
- PDB chain IDs occupy one column. Preprocessing remaps two-character chains and
  records `REMARK 900 CHAIN ID REMAPPED`; never compare those IDs directly with the
  deposited entry. Output PDB names parse from the right: the tail contains fixed
  residue-tag fields, while fragment/source heads are variable.
- **Hydrogen-free keys pool H states.** Only carbon atoms carry an H flag (`H0`/`!H0`);
  heteroatoms carry heavy degree (`D`) only (`frag_enumeration.atom_annotations`). `D` does not
  fix H count when H moves (tautomers) or protonation changes, so one slot pools donor and
  acceptor forms. A query
  acceptor then inherits contacts made to the donor form, and vice versa. For example, a
  pyridine-type `n` scores an Asp/Glu or backbone C=O as favorable because pooled `[nH]`
  observations donated to it. Riskiest cases:
  - Aromatic N tautomers: imidazole, pyrazole, triazoles, tetrazole, purines.
  - Lactam/lactim: 2-pyridone vs 2-hydroxypyridine, uracil-like rings. The ring `n` pools;
    the exocyclic O may key separately by bond order.
  - Ionizable groups: carboxylic acid, phosphate/phosphonate, tetrazole, acyl-sulfonamide and
    other acidic N–H, amidine, guanidine. Hit finding concatenates all charge-sign partitions
    (`hit_finder_core.py`), and `fragment_aliases.tsv` maps charged variants to neutral keys.
  - Phenol/phenolate and hydroxamic acid O–H.
  In this build, H states come from perception without waters/ions, so the pooled labels are
  unreliable too. Decision (2026-09-24): keep keys hydrogen-free. Evidence: 93% of slot tests
  lacked support, and the median H − noH contact-rate difference was −1.7 pp. Per-observation
  `*_cg_num_h` fields stay in buckets, so slot-level H filtering is possible (not implemented).
  Archive:
  `~/docking/frag_lib_validation/`.
- **Key reachability (W2 smoke test, 2026-09-24; 23 DR-38 xtals + 4 PiB).** Basis: keys that
  build-time CCD-template perception enumerates for each ligand (26 of 27; 9psz has no RCSB
  entry). 0 library keys unmatched by the hit finder, 0 charge-only misses. Of the 2,816 keys,
  957 are in the library; 1,678 are below `min_support` and 181 never occur in the parent DB.
  The hit finder also matches 72 extra library keys (11 structures): a key's implicit bond
  between aromatic atoms means "single or aromatic" in SMARTS, so all-aromatic keys match biaryl
  single bonds, while keys written with explicit `-` do not match fused rings (RDKit verified;
  the OpenBabel miner is assumed to follow the same Daylight semantics). Script and logs:
  `~/docking/scratch/W2/key_reachability.py`, `key_reach2_200031.log`.
- `~/docking/query_structs/*_xtal.pdb` ligand resnames are truncated 5-char CCD IDs
  (`A1JC2` -> `A1J`) that collide with real 3-char CCD codes. Resolve these ligands by SMILES
  (as `get_query_ligand_mol` does), never by resname -> CCD template.
- Fragment directories are literal names, not glob patterns: use `os.listdir`/
  `os.scandir` or `glob.escape`. Parse fragment keys with `utils.mol_from_fragment`,
  not sanitized `MolFromSmiles`.

### Geometry and counting

- Mining requires every SMARTS edge to satisfy the provisional broad envelope
  `0.70 <= d/(r_cov,i + r_cov,j) <= 1.25`. Rejection is per row, never per fragment.
  SMARTS bond class/order is preserved for diagnostics and later calibration; ambiguous
  bonds use the same broad fallback. This catches gross failures but deliberately accepts
  residual borderline contamination until element-pair/class bounds are calibrated.
- `utils.kabsch` returns `R, t, ssd` for `Y ≈ X @ R + t`; rank-0/1 fits have a
  meaningful SSD but no unique rotation. Chain breaks become NaN and are dropped by
  clustering, not penalized.
- `EXCELLENT_MATCH_CUTOFF` stops at the first permutation under 0.3 Å, so reported
  `best_*` is not necessarily the minimum; use cutoff 0 for the true minimum.
- `compare_cg_geometries.py`'s old `min_dist`/`frac_matched` were vacuous: the CG
  centroid was one of the Kabsch fit points, and the O(N_A x N_B) per-pair fit let
  dense libraries always "match" (2026-09-24, frac=1.00 for both a real pair and its
  negative control). Replaced by `mmd2`: each vdG's own vdM slot-0 N/CA/C is fit once
  onto a fixed ideal N-CA-C (never onto another vdG, never including the CG), then
  weighted unbiased squared MMD is computed between the resulting CG points. `mmd2` is
  an unbiased estimator of 0 for two samples of the same population, so it can be
  slightly negative -- `mmd2_self_A/B` (same-library, parent-split) is the noise floor,
  not required to be >= 0. 2-residue buckets anchor on slot 0 only (aa_bucket_parts
  order); slot 1 and the CG ride along rigidly -- no ideal 2-residue geometry exists.
- Enrichment uses crystallized PDB background and omits zero-observed AAs. Pairwise
  p-values are non-independent and uncorrected; positive-only Jaccard is optimistic.
  Profiles with different `bb_mode` are not comparable. Supply the same mode to both
  propensity counts and scoring.
- Placement (backbone-only match) at subset size 1: N/CA/C of one residue fits every vdG
  in the bucket (9rlo: median `bb_rmsd` 0.02 Å), so ranking by `bb_rmsd` is noise there.
  CG position depends on the unknown rotamer. Pose scoring (CG known) is unaffected.
- `benchmark_pose_recovery.py`'s `cg_rmsd` is part of a fit containing the CG, not an
  independent prediction metric. `plot_from_tsv.py` pools raw, structure-dependent
  scores; use mean per-structure Spearman as the headline statistic.

### Database and jobs

- Database trimming fingerprints use protein neighbors only; ligands, waters, and
  other hetero atoms do not count. Prepwizard can create stray atoms or relabel
  ligands as `ATOM`; current s01/s02 filters and restores known cases. If prepwizard
  flags change, verify `restore_renamed_ligands` on a known affected structure.
- Modified residues are classified by CCD polymer type, not atom names. Unknown CCD
  types are reported by s01; `X` remains the safety net.
- `make_sge_scripts_for_hit_finder.py` treats `--ref-pdb` as global for the batch.

## Known issues (settled)

- `apply_template_to_obmol` destroys deposited atom IDs. Read names from the PDB block;
  roster code therefore uses raw HETATM columns. CCD alternate names can also template
  an instance under a spelling enumeration does not emit; normalize through
  `index_by_name(template)` first, and remove placed hydrogens before lookup.
- Apply the record floor before protonation pooling so write-time and read-time
  representatives agree.
- Fragment-generation automorphisms and hit-finder `CalcRMS` intentionally differ;
  fragment symmetry lacks parent-boundary information. Ring/exocyclic over-merging is
  accepted, and fragments inherit unsanitized parent perception.
- Fragment keys match aromaticity exactly. Query bond perception can attach template
  chemistry to the wrong physical atom in an isomorphic but incorrectly perceived
  graph; this is query-side only and is not currently fixed.
- Protonation variants collapse to one representative. Promoted charge-stripped
  representatives need not be fragment-dict keys; resolve them through
  `fragment_aliases.tsv` or `scripts/lookup_fragment_key.py`.
- RDKit and OpenBabel interpret `r<n>` differently; the miner must retain its
  `match_satisfies_ring_sizes` re-filter. OpenBabel may miss very large macrocycles.
- ProDy can truncate overlong string fields (for example segnames) on write-back;
  npz dtype changes cannot prevent it.

## Open bugs

- Library side only: `frag_enumeration.group_lig_sites_by_overlap` (used by `get_fragments`) still
  merges sites by single linkage at 50% overlap. Query sites no longer use it: each distinct
  atom set is its own `q_site_idx` (2026-09-25).
- Some library directories use promoted SMARTS names absent from the fragment dict.
  This is expected; resolve through `fragment_aliases.tsv` or
  `scripts/lookup_fragment_key.py`, never by comparing directory names with dict keys.
