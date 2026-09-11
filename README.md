# ligand-vdGs

A pipeline for scoring docked ligand poses using van der Graphs (vdGs) — recurring interaction geometries between small-molecule chemical groups and protein residues, extracted from the PDB.

Builds on [vdG-miner](https://github.com/degrado-lab/vdG-miner) by Rian Kormos.

## Installation

From the `ligand-vdGs/` directory:

```bash
conda env create -f environment.yml
conda activate lig_vdgs
```

> **Note:** vdG-miner requires OpenBabel *with Python bindings* (`from openbabel import openbabel`), not just the CLI tool. The `openbabel-wheel` package in `environment.yml` provides this.

`environment.yml` installs the repo itself in editable mode (`pip install -e .`),
which is what lets `ligand_vdgs.*` imports resolve. If you set the environment up
some other way, run `pip install -e .` from the repository root once.

Everything is still invoked by path from the repository root, e.g.
`python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py ...`.

## Build the vdG Library

For full instructions (prerequisites, database preparation, SGE/SLURM submission), see the [Database Generation Guide](docs/database_generation_guide.md).

The core command, run once per fragment SMILES:

```bash
python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py \
    -s "<fragment SMILES, matched as SMARTS>" -c "<CG label; defaults to the SMARTS if omitted>" \
    -p <path/to/pdb_database/> \
    -b <path/to/probe_output/> \
    -o <path/to/vdg_library/> \
    --num-procs <n> \
    --subset-sizes 1 2
```

All requested subset sizes share one environment-reconstruction pass and are
written independently under `nr_vdgs/<subset_size>/`.

Each bucket npz holds one `nr_*` row per non-redundant vdG and one `mem_*` row per
clustered vdG that is not the nr one, so `cluster_size == 1 + (mem_ rows)`. Stage-1
clustering is sphere exclusion (Butina 1999 / GROMOS), which guarantees every member
lies within the cutoff of the vdG that represents it; the achieved radius is stored
per cluster. See
[docs/database_generation_guide.md](docs/database_generation_guide.md) for the full
array schema. A per-phase timing sidecar is written by default; `--no-profile-compute` skips it.

`-b` is the Probe output directory produced by preprocessing step 3.
`-c` (or its default, `-s`) is encoded with `utils.smiles_to_filename` before use,
so a SMILES containing `/` or `\` can be passed as-is.
Exact graph automorphisms are derived automatically from `-s` in
`clus_and_deduplicate_vdgs.py` and persisted to `nr_vdgs/cg_symmetry.npz`, which is
what consumers read (`vdg_npz_utils.load_cg_symmetry`) — independent of the
scheduler. Stored CG coordinates retain the SMARTS-slot order, while comparisons
that depend on atom correspondence enumerate the recorded automorphisms. Fragment
strings are parsed only with `MolFromSmarts`; they are not sanitized as standalone
molecules and no
hydrogens are inferred. Ordinary graph symmetry applies to every atom type. A
terminal-atom rule additionally ignores untrusted resonance, charge, and proton
placement: terminal N, O and S on B/C/N/O/P/S/Cl/Br/I centers are merged, plus
C-N, N-O, S-N and conjugated N-N normalization and an aromatic-N charge/H
normalization. Substituted, bridging and aromatic terminal atoms remain distinct, as
does a saturated center whose four or more single-bonded terminal atoms carry
non-uniform charges (the drawn charges are taken as deliberate).
See [docs/database_generation_guide.md](docs/database_generation_guide.md) for the
exact rule, and `docs/symmetry_edge_cases.md` for where it misfires.

On Wynton (SGE), generate one job per fragment and submit them:

```bash
# generate and submit in one step, longest job first (--mode is required)
MAX_H_RT=36:00:00 ./run_production_frags.sh --mode threshold-plus-include
```

The generated scripts use `#$ -cwd` with repo-relative paths, so submit from the
repository root.

- After a build or rebuild (and after any change to the fragment key scheme), run
  `python ligand_vdgs/tools/h_class_diagnostic.py --lib <vdg_library> --out h_class.tsv`
  before deciding which `H0`/`!H0` key variants to pool at read time; see
  ["When to run the H-class diagnostic"](docs/database_generation_guide.md#when-to-run-the-h-class-diagnostic).

---

## Bioisostere Identification

Identify candidate bioisosteres by comparing amino acid interaction profiles across chemical groups in the fragment library. See [docs/bioisostere_identification.md](docs/bioisostere_identification.md) for full details.

Steps 1 → 2 → 3 form a chain; `compare_cg_geometries.py` is an independent
branch that reads the vdG library directly rather than the profiles.

```bash
# 1. AA enrichment profiles per chemical group
python ligand_vdgs/identify_bioisosteres/compute_aa_profiles.py \
    --vdglib-dir <path/to/vdg_library/>

# 2. Pairwise profile similarity  (reads step 1's outputs from --profiles-dir)
python ligand_vdgs/identify_bioisosteres/compare_aa_profiles.py --single-aa

# 3. Visualize  (--sim-npz is the similarity matrix written by step 2)
python ligand_vdgs/identify_bioisosteres/visualize_bioisosteres.py \
    --sim-npz <path/to/similarity_matrix_*.npz>

# Independent: compare CG geometries directly between libraries
python ligand_vdgs/identify_bioisosteres/compare_cg_geometries.py \
    --vdg-lib-dir <path/to/vdg_library/>
```

---

## Hit Finding and Visualization

**Find vdG hits in query structures:**

```bash
python ligand_vdgs/score_poses/vdg_hit_finder.py \
    --smiles "<SMILES>" \
    --query-dir <dir_of_pdb_files> \
    --vdg-lib-dir <path/to/vdg_library/> \
    --nprocs 4 \
    --outdir <output_dir>
```

`--query-dir` takes a *directory* of structures, not individual files.

Scoring-relevant options, none of which appear above: `--rmsd-threshold`
(default: derived per combo from `normalize_rmsd`, matching the window used
during library clustering), `--contact-cutoff` (off by default; 3.8 Å is
reasonable), `--ref-pdb`, `--min-shared-atoms`, and `--no-dedup`.

**Hits are deduplicated by default.** Records sharing a BSR combo and
overlapping ligand atoms collapse to the single lowest-RMSD hit, so
`results_summary.txt` counts are post-dedup. Pass `--no-dedup` for raw counts.

On Wynton (SGE), `python ligand_vdgs/score_poses/make_sge_scripts_for_hit_finder.py`
generates one job per query set.

**Write matched geometries as PDB files for PyMOL:**

```bash
python ligand_vdgs/tools/write_vdg_hit_pdbs.py \
    --hits-tsv <path/to/vdg_hits.tsv> \
    --vdg-lib-dir <path/to/vdg_library/> \
    --outdir <output_dir>
```

Output is gzipped (`.pdb.gz`), which PyMOL opens directly, and
`--max-per-frag-bsr` (default 10) caps how many are written per leaf. File and
directory names follow the same grammar as the library materializer — fields joined
with a single `_`, residue tags of exactly four fields (`seg_chain_resnum_resname`),
`~<n>` for duplicates — described under
[Output file names](docs/database_generation_guide.md#output-file-names).

Each `<query>/vdg_matches/<frag>/<bsr>/` leaf must be empty when the run reaches
it: files are never overwritten or skipped on write, so leftovers from an earlier
run can't be mistaken for this run's output. There is no `--overwrite` and no
resume after a partial run — clear `--outdir` by hand and rerun.
