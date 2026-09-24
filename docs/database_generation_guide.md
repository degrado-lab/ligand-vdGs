# Database Setup Guide

Builds a current-schema vdG library from a protonated PDB database. Run from repo root;
editable install exposes `ligand_vdgs.*`.

## Install

```bash
conda env create -f environment.yml
conda activate lig_vdgs
```

Installs the repo (`pip install -e .`). Requires OpenBabel Python bindings
(`from openbabel import openbabel`); the CLI executable alone is insufficient.

## 1. Prepare the parent database

Layout: `<pdb-dir>/<lowercase inner two chars>/XXXX.pdb` (biounits `XXXX_1.pdb` also read).
`scripts/format_parent_database.py` copies/validates standard filenames; biounit mirrors are
laid out by hand. Generation discovers structures through `functions.parent_db`; preserve the
inner-two-character directory layout and filename stems.

**Optional trim/repair:** `preprocessing/s01_trim_database.py` (module-level settings, no CLI)
deduplicates structures and extracts 20 Å binding sites. Use `skip_to_output_pdbs = False` on a
fresh database. It remaps two-column chain IDs to unique one-column IDs (as `REMARK 900`) and
converts peptide-linked HETATM residues to ATOM; repair a database trimmed before this fix with
`scripts/remap_chain_ids.py`. Modified residues are classified via
`resources/ccd_polymer_types.tsv` (MSE/SEC→MET/CYS, other PEPTIDE LINKING→ATOM,
NON-POLYMER→HETATM ligand); regenerate with `scripts/fetch_ccd_polymer_types.py` for newer CCD
entries. Both s01 and s02 drop amino acids missing N/CA/C and one-heavy-atom carbon HET residues.

**Required protonation:** every input PDB needs hydrogens (Reduce2 or Schrödinger PrepWizard via
`preprocessing/s02_run_prepwizard.sh`). Its array task count must match `num_batches_total = 15`:

```bash
qsub -t 1-15 ligand_vdgs/preprocessing/s02_run_prepwizard.sh 15
```

s02 reapplies the orphan filter and restores free ligands PrepWizard relabeled as protein. Its
output is the `--pdb-dir` used below.

## 2. Build the fragment dictionary

The roster pass reads the parent PDB database. It records CCD-templated ligand types and the
canonical heavy-atom names observed in each biounit; OpenBabel-fallback and unreadable ligand
instances do not enter the roster. NCS copies count once per biounit. The fragment pass reads
that roster and the CCD templates, enumerates connected induced fragments, then counts support
as distinct biounits with every fragment atom observed. A SMILES-only ligand cannot qualify.

```bash
python ligand_vdgs/generate_vdgs/build_ligand_roster.py \
  --pdb-dir <protonated-pdb-dir> --out resources/ligand_roster.pkl --num-procs 8
python ligand_vdgs/generate_vdgs/fragment_database_ligs.py \
  --roster resources/ligand_roster.pkl --outdir <output-dir> --num-procs 8
```

- Enumerates connected induced subgraphs (4–5 heavy atoms by default) and writes
  `<output-dir>/database_frags_dict.pkl`.
- Chemistry comes from CCD templates. To add a non-CCD ligand, add a CCD-style template under
  the exact PDB residue name, then rebuild the roster and dictionary.

## 3. Select and generate vdGs

A fragment qualifies if it passes `Frags.is_organic`, has no `UNDESIRED_ELEMENTS` (metals,
lanthanides, noble gases, Si/Se/As/Te; boron allowed), ≤5 heavy atoms, is not a halogen
oxyanion, and has support ≥ `--min-support` (distinct parent biounit stems, pooled across
protonation variants, computed over the full dictionary).
`include-only` does not lower a threshold: it builds only explicitly named fragments 
regardless of support. Raising the threshold changes the selected vocabulary and requires 
a fresh library build.

```bash
python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py \
  -s "<SMARTS>" -p <protonated-pdb-dir> -o <vdg-library> --num-procs <n>
```

`-s` is parsed with RDKit `MolFromSmarts` (no sanitization/valence checking/H inference).

| Option | Meaning |
|---|---|
| -s, --smarts | Required fragment SMARTS |
| -c, --cg | Output label (default -s; encoded via `utils.smiles_to_filename`) |
| -p, --pdb-dir | Required protonated parent database |
| -o, --out-dir | Required library root; fragment subdir must be empty/absent |
| --num-procs | Processes (default 10) |
| --subset-sizes 1 2 | Only 1/2 accepted; default both |
| --no-profile-compute | Skip `<cg>_compute_profile.json` sidecar |
| -m, --max-num-vdgs-to-clus | Debug cap on PDB IDs/bucket; never production |

Outputs cluster under `nr_vdgs/<subset-size>/<pos|neut|neg|unreadable>/`. No overwrite/resume
flag.

**Geometry and symmetry:**

- Stage 1 uses deterministic sphere-exclusion (Butina/GROMOS) clustering on CG atoms and each
  vdM's N/CA/C. Stage 2 subdivides by flanking sequence and CA similarity.
- Pose distance minimizes over the CG automorphism group and interchangeable same-label vdM
  slots. Stored atom order is not authoritative.
- Fragment-generation symmetry differs from hit-finder `CalcRMS` symmetry.
- Each SMARTS bond must satisfy the Cordero envelope
  `0.70 <= d/(r_cov,i+r_cov,j) <= 1.25`; otherwise the row is rejected
  (`stream_skipped_cg_bond_geometry`). Per-element calibration is deferred.
- `cg_symmetry.npz` records automorphism normalization. Terminal N/O/S atoms on one
  B/C/N/O/P/S/Cl/Br/I center may exchange within element, except substituted, bridging, or
  aromatic atoms. With four or more terminal atoms, mixed charges retain their distinctions.
- Environments are reconstructed once per fragment for all requested subset sizes.

### SGE (Wynton)

```bash
python ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py \
  --vdg-lib-dir <library> --pdb-dir <pdb-dir> \
  --log-dir <logs> --sge-out-dir <empty-script-dir> \
  --frags-dict resources/database_frags_dict.pkl --min-support <n>
python ligand_vdgs/generate_vdgs/submit_frag_jobs_by_size.py <empty-script-dir>
```

Cost is estimated inline by sampling `--pdb-dir` (`--sample-size`/`--estimate-procs` control the
sample) and sizes jobs into tiers; `--h-rt`/`--num-procs` fix a single value for every fragment
instead. Fragment prep is memoized under `$VDG_SCRATCH/prepared_fragments/`; the generated
`submission_order.tsv` controls submission, descending by slots and then estimated tier count.

```bash
MIN_SUPPORT=<n> ./run_production_frags.sh --mode threshold-plus-include --no-submit
MIN_SUPPORT=<n> ./run_production_frags.sh --mode threshold-plus-include
```

`threshold-plus-include` requires no existing fragment directories and builds the threshold set
plus INCLUDE. `include-only` requires an existing library and builds only explicitly named
INCLUDE fragments, without applying the threshold. `--no-submit` generates scripts and still
updates `fragment_aliases.tsv`; it does not submit jobs. Job tiers use SMARTS passes over
`--pdb-dir`. Force a wanted fragment with INCLUDE/`--include`.

For direct script generation, `--resume` requires matching library provenance and skips
fragments whose completion marker is present. It rewrites scripts for unfinished fragments, but
does not clean partial output directories; remove a leftover non-empty fragment directory before
resubmitting it because the wrapper refuses that output path. The production shell wrapper
rebuilds its script directory each run, and top-up mode filters completed fragments and reports
partial directories before submission.

### SLURM or another scheduler

```bash
python ligand_vdgs/generate_vdgs/extract_fragment_smiles.py \
  --frags-dict resources/database_frags_dict.pkl --min-support <n> \
  --output <fragment-work-list.txt>
```

Pass `--vdg-lib-dir` so aliases land at `<vdg-lib-dir>/fragment_aliases.tsv` (otherwise a
`<stem>_aliases.tsv` is written beside the list and library readers won't find it). Submit one
wrapper invocation per line:

```bash
while IFS= read -r smarts; do
  your-submit-command python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py \
    -s "$smarts" -p "$PDB_DIR" -o "$OUT_DIR" --num-procs "$NPROCS" --subset-sizes 1 2
done < fragment-work-list.txt
```

## 4. Validate the library

```text
<library>/fragment_aliases.tsv
<library>/<cg_label>/<cg_label>_log
<library>/<cg_label>/nr_vdgs/cg_symmetry.npz
<library>/<cg_label>/nr_vdgs/{1,2}/{pos,neut,neg,unreadable}/<aa_bucket>.npz
```

- `nr_*` rows are coordinate-bearing cluster representatives; `mem_*` rows are identity-only
  members linked by `mem_cluster_id`. Observations are `nr+mem`; clusters are `nr`.
- `cluster_size` counts observations. `cluster_num_parents` counts distinct PDB
  entries/depositions, not support or cluster size. Selection support counts distinct biounit
  stems.
- Parent records store a biounit stem and `parent_pdb_dir`. For copied libraries, use
  `vdg_npz_utils.resolve_parent_pdb_path`; use `vdg_npz_utils.rederive_member_coords` for member
  coordinates.
- H-class fields use negative values for unreadable atoms; zero is a real count. Carbonyl-O
  coordinates are separate and excluded from RMSD.

After the fleet drains:

```bash
python scripts/check_library_after_build.py --vdg-lib-dir <library>
python ligand_vdgs/tools/h_class_diagnostic.py --lib <library> --out h_class.tsv
```

The first checks aliases, charged-forward resolution, completion markers, and H-class fields. A
directory's existence doesn't prove completion — direct walkers must call
`Frags.check_vdg_job_status` (as `load_vdg_bucket` does). Run the H-class diagnostic after every
full build/rebuild and whenever fragment keys change.

Optional, for troubleshooting: the H-class check above only runs `h_class_diagnostic.py
--dry-run` internally and reports pass/fail. If it fails, or you want to see which buckets are
affected, run the same check standalone for detail:

```bash
python ligand_vdgs/tools/h_class_diagnostic.py --lib <library> --dry-run
```

## 5. Inspect selected vdGs

`materialize_vdg_pdbs.py` writes PDBs for PyMOL. `-c` must be a fragment's `nr_vdgs/` directory;
`-o` must be empty or absent:

```bash
python ligand_vdgs/generate_vdgs/materialize_vdg_pdbs.py \
  -c <library>/<cg_label>/nr_vdgs/ -o <empty-output-dir> --top-clusters 10
```

Use `--aa-buckets`/`--subset-sizes` to filter. Default mode writes nr vdGs; `--members` adds
individual members and requires `--top-clusters` or `--clusters`; `--reps` includes the nr vdG in
its cap.

Output names have variable heads — parse from the right: strip `.pdb.gz` and any `~<n>` suffix,
then take the final `4 * (subset_size + 1)` underscore-separated fields as four-field residue
tags (infer subset_size from AA-bucket tokens in nr mode, from the path in members mode). Hit
files use the same residue-tag rule but have no ligand tag and end with an RMSD.
