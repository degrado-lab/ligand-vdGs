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
- Compare time-varying npz schemas after parsing JSON and removing `build_date`.
  Identify parent databases by `functions/db_identity.identity_of`, never by path.
- `cg_formal_charge == -1` is ambiguous. Treat `cg_heavy_degree == -1` as the
  unreadable sentinel; unreadable CGs belong in `unreadable`, not `neg`.
- `load_vdg_bucket`/`vdg_npz_path` require the fourth positional `sign` argument.
  Pre-DR-61 flat/schema-v2 libraries must fail with `BucketSchemaMismatch`, not look
  like empty buckets.

## Assumptions

- A query structure has exactly one ligand resname. Multiple resnames select one by
  `n_heavy * n_residues`, or return an error for that model. A ligand is one HETATM
  residue: cross-residue polymers and covalent protein links are unsupported and are
  not represented by the library.
- Query tautomer and bond orders come from the CCD template, not density. Tautomers
  are not enumerated. H-stripping intentionally merges aromatic N roles; N/O degree
  stays out of the key. H status is the per-observation `cg_placed_h` label (DR-62).
  Formal charge must survive any new query-molecule construction path.
- Fragment cost estimates sample `--pdb-dir`, not the CCD. Novel ligands must be in
  that parent database and the estimate must be regenerated; use `--include` for a
  wanted fragment with zero observed occurrences.
- Large/symmetric CGs can make RMSD alignment combinatorial. BioLiP2 includes
  pseudo-ligands requiring preprocessing cleanup.

## Caveats

### Representation and parsing

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
- `min_dist` depends strongly on sample size and is further deflated because the nr
  vdG is fitted. Compare similarly sized buckets or use a stable percentile.
- `compare_cg_geometries.py` reports subsample-based match fractions, first-residue
  enrichment, and its own backbone/X counts; its heatmap `vmax` is fixed at 3.0.
- Enrichment uses crystallized PDB background and omits zero-observed AAs. Pairwise
  p-values are non-independent and uncorrected; positive-only Jaccard is optimistic.
  Profiles with different `bb_mode` are not comparable. Supply the same mode to both
  propensity counts and scoring.
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

- `group_lig_sites_by_overlap` uses single-linkage union-find, so disjoint sites can
  merge transitively and share one `q_site_idx`. Complete-linkage would fix it, but
  would change hit-finder semantics.
- Some library directories use promoted SMARTS names absent from the fragment dict.
  This is expected; resolve through `fragment_aliases.tsv` or
  `scripts/lookup_fragment_key.py`, never by comparing directory names with dict keys.
