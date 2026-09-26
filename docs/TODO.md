# TODO, decision records, and cross-session inbox

Single tracker for open work, open decisions (DR-*), and session-to-session messages.
Delete an item and its citations once resolved.

## Current library status

Known defects baked into this build (fixes deferred to next parent DB / rebuild, see §3):
- Preprocessed DB has zero waters/monatomic ions — corrupts protonation everywhere;
  metal-coordination vdGs can't exist in this library; defer to the next rebuild.
- 535 roster instances carry a stray orphan ammonium atom glued onto the ligand, defeating
  the CCD template and dropping the whole instance; defer to the next rebuild.
- Perception changed after the roster was built (H handling, legacy atom names, no more
  element-symbol renaming) — roster/fragment-dict/library are stale relative to current
  code, and no gate (`load_roster`'s `ccd_identity`/`db_identity` checks) catches it.

## 2. Analysis (bioisosteres, hit finding, scoring)

- **Rebuild coupling (analysis code written against the DR-63 library).** Revisit on rebuild:
  `vdg_npz.load_vdg_bucket` (`hit_finder_core.py:486-487`, `compare_cg_geometries.py:211`),
  `check_vdg_job_status` (hit finder, `identify_bioisosteres/{common,compare_cg_geometries,
  compute_aa_profiles}.py`), `nr_*`/`mem_*` field names, `BUCKET_SCHEMA_VERSION`, fragment
  keys (perception drift, §"Known defects"). Workers append coupling points here as found.
  - W3 (`tools/self_recovery.py:78-90,180-195`): reads buckets directly: `{nr,mem}_parent_biounit`,
    `{nr,mem}_cg_{chain,resnum,resname,names}`, `{nr,mem}_scrr_*`, `cluster_id`/`mem_cluster_id`,
    `cluster_num_parents`. Hit-record keys (`:51-53`): `bsr_combo` format `seg:chain:resnum;...`
    (`:133`), `lig_instance`, `vdg_cluster_*`, `vdg_rmsd`, `q_atom_indices`. Also
    `score_one_model_multi_instance` 5-tuple return (`:215`), `dock.get_bsr_combinations` 4-tuples
    (`:236`).
  - W4 (`tools/score_eval.py`, planned): M5 enrichment reads per-(fragment, subset size, `aa_bucket`)
    totals of `cluster_num_parents` across the whole library; `hits.npz` contract fields (plan_W4).
- **DR-36 (DECIDED 2026-09-24, user)** — (a) CG-free matching for recovery/placement.
  - Problem: hits are matched on joint bb+CG RMSD (`hit_finder_core.py:697`); recovery
    sampled that way is a lower bound.
  - Decision: add a placement mode (backbone-only match, per-hit held-out CG error). Joint
    matching stays for pose scoring (does this pose match a vdG).
  - Cost accepted: lever-arm CG error grows with distance from backbone — measure it.
  - Match cutoff: `bb_rmsd <= tau * sqrt(n_atoms / N_bb)` (joint hits ⊆ placement hits).
  - Revisit fired 2026-09-24 (9rlo pilot: 3.1M hits). Subset 1 is degenerate: one residue's
    N/CA/C matches the whole bucket, and CG position depends on the rotamer.
    - Revised ruling (user skeptical, n=1 structure): subset 1 is NOT demoted. W4 evaluates
      subset 1, subset >= 2, and pooled as separate strata, with every method. Consensus
      (M4) and cross-fragment consistency (M7) may recover the rotamer mode.
    - Decide only on dev results.
  - Never use side-chain coordinates in placement (user: rotamers are unknown when docking
    or generating poses). Library side-chain storage is ruled out.
  - Known limitation (accepted for this benchmark): BSR selection uses the current code's
    definition (true-ligand contacts, may use query side chains), so placement results are
    "conditioned on the true contacting residues". A deployable placement mode must choose
    BSRs without the ligand.
- Decoys: synthetic decoys (perturbed CG, non-contacting residue combos) are enough for now.
  Eventually: add a real, carefully curated decoy set (user).
- **DR-38 (DECIDED 2026-09-24, user)** — cohort: 23 hand-picked PDB entries deposited after
  the library build, ligands absent from the library: 9r30 9rci 9gs6 9r8y 9y76 9d5q 9gx3
  9fqy 9d0s 9foc 9rlo 9jj5 9rfe 8yqe 9ddg 9g0h 7hc7 9ei7 9fs0 9psz 9r9k 9qff 9d9i, at
  `~/docking/query_structs/<id>/*_xtal.pdb` only (Boltz `*_{1..10}.pdb`, `*_preds/` unused). De novo
  set: `query_structs/PiB_{mef,nir,ruc,vel}` (crystal of designed complex). External sets
  (W3) must state date cutoff and ligand-novelty filter the same way. Parent-DB cutoff
  (RCSB, 2026-09-24): max deposit 2023-12-14, max release 2024-01-03, over 65,598 IDs
  (48 obsolete); `~/docking/benchmarks/parent_db_dates.tsv`. The BioLiP2 download is from
  Aug 2024 (user); contents end at the Jan 2024 release. External-set filter: released after
  2024-01-03, plus exclusion by PDB ID and roster ligand.
- Concern (user): fragments with large libraries / nonspecific vdW contacts may dominate
  small-library / polar-directional ones in scoring and bioisostere comparison. W3
  self-recovery measures this before any normalization is chosen.
- Protonation-matching strictness: charge-loose today (`charge_normalized_fragment`,
  `group_protonation_variants`, `extract_fragment_smiles.py:68,80`). Decide charge-strict vs.
  keep pooling + rely on `hit_finder_core._warn_charge_only_miss:237`.
- Any library-wide statistic vs. a null must permute at the parent-biounit level (structure-
  level permutation null failed where a multinomial null passed, 2.7x too permissive —
  2026-09-10 contact-gate analysis).
- Nested fragments (e.g. `cccn` inside `cccnc`) share observations — correlated counts across
  fragments. Buckets themselves are independent (clustering unaffected); normalize only when
  aggregating across fragments (bioisosteres, hit scoring).
- Calibrate Stage-1 RMSD cutoff `T` against recall on real mined environments (not stored nr
  vdGs) — blocks the Stage-2 flanking-context question (whether splitting pose clusters by
  flanking context is a better use of representative budget than a tighter `T`).
- Query/library key reachability: build an explicit relaxation map (protonation, H-count)
  instead of the ad hoc `resolve_fragment_alias` + charge-only-miss warning. Ring-size
  differences stay hard constraints (`frag_enumeration.ring_query_for_atom:116`).
- `_has_deliberate_charge_assignment` (`utils.py:72`, called at `:109`): only fragment that
  would exercise it is below the selection threshold — delete, or lower the threshold.
- Overlapping fragments. Decision 2026-09-25 (evidence in
  `private/project_planning/W5_overlap_proposal.md`):
  - Keep raw hits. Site grouping merges only identical atom sets, not single linkage at 50%
    (W2). `deduplicate_hits` is deleted (W2); it dropped the largest site in 30.5% of DR-38 BSR
    sets and keyed on truth in bb mode.
  - Aggregation across fragments happens in scoring, over an overlap graph built from topology or
    placed geometry (W5, `functions/fragment_overlap.py`). The default is containment-aware.
  - Open: nested fragments share ~93% of parent observations, so summing support double-counts.
    Per-observation support is possible without a rebuild: the library stores mem-level parent
    identity (biounit, ligand seg/chain/resnum, vdM residues, CG atom names; W5). Bioisostere counts also need the containment-aware rule (W5).
- Icode handling (de facto policy already in code): residues with >1 distinct icode are
  skipped on both sides (`dock_utils.py:47-49`, `align_and_cluster.py:761-762`); a lone icode
  is fine.

## 3. If re-running vdG generation (new `~/docking/frag_lib`)

- Rebuild roster, fragment dict, and library together (perception drift since last build —
  see "Known defects" above).
- Drop write-only bucket fields `{nr,mem}_cg_num_h`, `_cg_placed_h`, `_cg_heavy_degree`
  (only reader: `post_build_buckets.contract_problems`) and bump `BUCKET_SCHEMA_VERSION`.
  Keep the in-memory `heavy_degree` annotation (`clus_and_deduplicate_vdgs._charge_sign`).
  Only in a full rebuild: the changed schema string breaks `include-only` top-ups into an old
  library. Skip if H-filtered scoring (pitfalls, "Hydrogen-free keys") is planned.
- Fix orphan-atom-in-ligand bug upstream in prep (535 roster instances affected). Verify with
  `~/docking/scratch/roster_fallback_diag/diag.py` — `N/7` rows should vanish.
- Re-decide `sasa.MIN_CONTACT_AREA` (currently 0.0) against geometric recurrence once a
  library exists: per contact-area band (0-0.25, 0.25-0.5, 0.5-1, 1-2, 2-5, >5 A^2), check
  `cluster_num_parents` histogram (distinct PDB entries, never `cluster_size`) and
  `cluster_pose_radius`. Tight clusters recurring across unrelated biounits = real geometry;
  diffuse singletons = noise. Cross-check pi classes separately. Re-run the CG-scoped sweep
  too (smaller fragments have less buriable surface).
- Cation-pi / S-pi to protein rings: 96%/94% of missed pairs are inside the 6.5 A prefilter
  with zero credited area (ligand-side moiety is a single atom, self-buried). Widening the
  prefilter buys ~4%. Decide whether these need a non-SASA term.
- Calibrate element-pair/bond-class-specific bond-length bounds against trusted structures
  (current broad Cordero envelope admits some suspect short C-S/aromatic C-N and long C-O/C-C
  bonds).
- Re-run s01* and preprocessing before frag_lib generation to fix parsing warnings (many
  environments discarded due to duplicates).
- Quantify current errors and fix: phosphate fragment cost is extreme (47G, many hours).
- Convert all "centroid" terminology to "medoid" throughout.
- Linear 5- and 6-atom CGs from both linear and ring structures collapse into the same
  cluster, losing information — needs an actual fix.
- Decide `include_water=True` (or equivalent) for vdG-miner, or defer further downstream.
- Add RSCC/RSR/RSRZ to per-vdG quality records — plumbing was removed, not left as a no-op
  (`tests/test_vdg_miner_environments.py:140-142` asserts they're not entry-point params), so
  this is a from-scratch write.
- CG element-pinning gap: `dock_utils.cg_element_symbols:157-177` returns None for any
  unpinned (OR/negation) SMARTS atom (atomic number 0); callers must treat those slots as
  matching any element (`clus_and_deduplicate_vdgs.py:1488-1499`, `dock_utils.py:191-199`;
  `hit_finder_core.py` never calls this). Believed unreachable today — that's a claim about
  the fragment dict, not a code guarantee. Add a pitfalls.md entry if the dict ever gains such
  a key.
- Stage-1 clustering merges mirror-image (face-inverted) CG poses: the minimizing automorphism
  can superimpose a near-planar CG's substituents while its central atom sits on the opposite
  face. Add a signed-volume (face) consistency check to stage-1 edges, generalized beyond P,
  and measure over-splitting before adopting. Decided: face check, not a tighter cutoff.
  Details: `private/project_planning/cg_face_check.md`.
- Library sites come from one Kekule form per parent (`frag_enumeration.get_fragments` bonds one
  mol), so a bond-order key (`P(=[O;D1])[O;D1]`, `C=[O;D1]`) records only O-sets containing that
  parent's double-bonded O (2j9y/FOO: never {O1P,O3P}). The query side now unions all
  terminal-resonance forms (`hit_finder_core.query_resonance_forms`). A rebuild should do the same.
- Scrub irrelevant PDBs (PSI/PSII, ribosomes, heavy Fe-S cluster complexes). The >62-chain
  skip in s01 incidentally catches some but is a chain-ID-encoding guard, not a content scrub.

## 4. When updating the parent database

- Zero waters/monatomic ions in the preprocessed DB (HOH/ZN/MG/CA/NA/CL/MN/FE, checked all
  65,600 files) — corrupts protonation and blocks metal-coordination vdGs. Raising trim
  radius does not fix this.
- Exercise orphan-backbone-N (`_prep_filters.incomplete_backbone_resindices:22`) and
  modified-residue handling (`_prep_filters.modified_residues_to_protein:314-395`) on a fresh
  BioLiP input — both are fixed and applied by s01 but untested since the fix landed.
- `UNDESIRED_ELEMENTS` (`fragment_database_ligs.py:217-222`) is a denylist that already leaks
  (12 CCD elements are neither organic nor denylisted). Switch to an allowlist
  (C,N,O,S,P,F,Cl,Br,I,B) before adding non-CCD ligands. Se stays denylisted (:219); B stays
  allowed (real covalent warheads in boronic acids/esters).
- Remove cryoprotectants in `fragment_database_ligs.py` — no resname skip list exists there
  (only `extract_fragment_smiles.is_solvent_artifact:18`, which covers halogen oxyanions
  only).
- H handling for non-CCD ligands: `ligand_perception.py:182` deletes H after
  `PerceiveBondOrders` (needs positions first). CCD-template path doesn't care (H counts come
  from template). Open question is narrower now: a new non-CCD ligand has no template, lands
  on OpenBabel fallback, and does need H — decide how to build them.
- Consider a `keys_only` fast path in `fragment_database_ligs.py` skipping the
  permutation/site-grouping machinery (unneeded when only the element string is used) —
  verified no-op on CCD output, doesn't exist today.
- Ligand-inclusion floor: require >=1 C, >=4 atoms (threshold undecided, carried over).
- Strip salts before calculating formal charge — a co-crystallized counterion isn't part of
  the compound.
- Apply quality filters to fragment atoms, not the whole ligand (some ligands are huge; a
  fragment-level filter is the right granularity).
- Re-evaluate PLIP vs. Probe for contact detection — likely superseded by the SASA refactor;
  confirm whether this is resolved. PLIP auto-adds H (trust unclear, want an off-flag) and by
  default keeps only one altloc (flag exists to change that).
- BioLiP2 includes NMR/EM structures, not just crystal.
- Can't reconstruct the biological complex from annotations alone — BioLiP2 may flag
  "chain A binds lig A" while the ligand also makes tertiary/quaternary contacts. Proposed
  workflow: use the BioLiP2 flag as a signal, then build the protein from that.
- H-bond caution: placed H's are an artifact of the modeling program, not density (most
  crystal structures have no H density). An MM forcefield used in refinement may
  under-model CH H-bonds even where real. Lean toward storing heavy-atom positions plus the
  H-placement method, not trusting placed H's outright. Open question whether QM/MM scoring
  downstream still needs explicit H's.
- Store altloc info instead of collapsing to one, when an alt conformation has an interesting
  interaction.
- Altloc resolution before PrepWizard (needed for planned PanDDA depositions). The current
  `prepwizard_BioLiP2_repaired` has 0 altloc labels in 65,598 files (scan 2026-09-25):
  PrepWizard collapsed them, and which label it kept is unaudited.
  - Move `prep_benchmark_set.resolve_altlocs` (highest-occupancy ligand label, near-tie
    flag, per-residue protein fallback, never drop atoms) into `functions/`. Share it across
    preprocessing, mining, and benchmark prep. It replaces `vdg_struct_utils._pick_best_altloc`,
    which picks per backbone atom (occupancy, then 'A') and can mix N/CA/C from different altlocs
    within one residue.
  - Deferred audit: for about 500 sampled raw RCSB mmCIFs, measure ligand-altloc prevalence
    (PanDDA vs not) and whether PrepWizard kept the highest-occupancy label.
- Store ASU, biological assembly, and symmetry-mate information. Build crystallographic
  symmetry mates: ligands at a lattice interface are otherwise invisible. Biounit copies and
  symmetry mates are not independent evidence (W5 redundancy metric).
- Joint X-ray/neutron structures (good heavy-atom coords + good H/proton positions in the
  same structure) are a good validation subset for unusual H-bonding.
- BioLiP(3) accessed 2026-08-30; `ligand.tsv.gz` downloaded.

## 5. De novo readiness (score_poses/, for scoring designed/unsolved models)

Full detail: `private/project_planning/de_novo_readiness_audit.md`. Audit only, no fix
implemented; not why benchmark numbers are embargoed (that's DR-36/DR-38, §2).

- B2: multi-instance ligand fallback (mixed composition / SMILES parse failure) still
  hard-fails instead of degrading gracefully.
- B3: undiagnosed OOM on large ligands (1zap/A70, 124 atoms); repro qsub written, not run.
- B4: per-structure scoring cost — backbone-RMSD batching fix applied and verified, but
  wall-time win unmeasured and a second `kabsch_ssd` site untouched. Combo-necessity
  attribution (4ij8, 1zap) shows the real waste is 99%+ unnecessary frag x combo tasks,
  not batching — a third structure (6ift) is queued to confirm, blocked on library
  post-build validation (§1).
- A1-A11: eleven open assumptions (modified-residue resnames, occupancy-as-data-channel,
  crystal-tuned geometric cutoffs, etc.) that could silently degrade de novo scoring —
  none fixed, see file for specifics.

## 6. Test-suite vacuity audit — open items

Full detail: `private/project_planning/test_vacuity_audit.md`. Not a build gate unless
that behavior is exercised or changed again.

- **S2 (HIGH, OPEN)** — `utils.kabsch` forces a proper rotation even for degenerate
  (collinear) fits, so `hit_finder_core.held_out_placement` can return a bogus distance instead of
  `None`. Needs a rank/degeneracy check; numerical threshold needs spec confirmation.
- `tests/test_resource_tiers.py:87` fails today (pre-existing, not vacuity):
  `resources_for()` no longer takes `short_queue`.
- Partial audit (23/~25 files, fan-out killed early): two CONFIRMED-vacuous tests
  (`test_butina_clustering.py:73,87`, `test_env_assembly_guards.py:161`) plus V3-V17
  weaker findings and a clean list — see file.

## Deferred design decisions (no timeline)

- `private/project_planning/literature_min_support_redundancy.md` — phosphate/redundancy
  follow-up.
