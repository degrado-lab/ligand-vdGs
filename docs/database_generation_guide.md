# Database Setup Guide

This guide explains how to generate a van der Graph (vdG) database for small molecule functional groups. If you prefer to use a pre-generated vdG database, you may skip the instructions below and download our database [here](insert-url-here).

> **All commands in this guide assume you are running from the root directory of the `ligand-vdGs` package.**

## Installation

1. **Set up the Python environment.** A conda environment file is provided at `environment.yml`:

   ```bash
   conda env create -f environment.yml
   conda activate vdg_env
   ```

2. **Install this package** with your environment active:

   Pick one:
   ```bash
   pip install -e .   # editable — source changes take effect immediately (recommended during development)
   pip install .      # standard install for end users / distribution
   ```

> **Note on OpenBabel:** The vdG-miner component requires OpenBabel *with Python bindings* (`from openbabel import openbabel`), not just the command-line tool.

## Prerequisites: Set Up a Parent Database

Before running any pipeline steps, you need a preprocessed PDB database to extract vdGs from. You can use any collection of PDB structures — your own custom set or a mirror of the RCSB PDB.

**Directory layout.** Structures must be organized in RCSB mirror format: each file is named `XXXX.pdb` (4-character code) and placed in a subdirectory named after the inner 2 characters of the code. For example, `1ABC.pdb` → `AB/1ABC.pdb`. Use [`ligand-vdGs/scripts/format_parent_database.py`](../scripts/format_parent_database.py) to reformat an existing directory.

Once your source PDBs are in the right layout, complete Steps 1–3 in order before moving on to Steps 4–5.

- **Step 1 (`s01_trim_database.py`) — Prune the database (optional, recommended).** Run [`ligand-vdGs/ligand_vdgs/preprocessing/s01_trim_database.py`](../ligand_vdgs/preprocessing/s01_trim_database.py) to filter by ligand b-factor and extract 10 Å binding-site regions, reducing database size and removing redundant structures. The pipeline works without this step, but skipping it means running on full PDB files (slower, larger, and noisier). If you run it, use the output as `-p` in Step 5 below; otherwise use the formatted database from above.

- **Step 2 (`s02_run_prepwizard.sh`) — Add hydrogens.** The pipeline requires hydrogens to be present in every PDB. Use [Reduce2](https://github.com/cctbx/cctbx_project/tree/master/mmtbx/reduce) (open-source) or Schrödinger's PrepWizard (more accurate, requires a license). To run PrepWizard jobs in parallel on an SGE cluster, see [`ligand-vdGs/ligand_vdgs/preprocessing/s02_run_prepwizard.sh`](../ligand_vdgs/preprocessing/s02_run_prepwizard.sh). Run this on the output of Step 1, or the formatted database above if you skipped trimming.

- **Step 3 (`s03_run_probe.sh`) — Run Probe.** Probe must be installed separately; it is available from the [Richardson lab / MolProbity project](https://github.com/rlabduke/probe). Run [`ligand-vdGs/ligand_vdgs/preprocessing/_run_probe.py`](../ligand_vdgs/preprocessing/_run_probe.py) on each PDB to compute atomic contacts. To submit all PDBs as an SGE array job, use [`ligand-vdGs/ligand_vdgs/preprocessing/s03_run_probe.sh`](../ligand_vdgs/preprocessing/s03_run_probe.sh). The output directory of Probe files is used as `-b` in Step 5 below.

## Step 4: Build the Fragment Dictionary

Build a fragment dictionary from your database ligands by running [`ligand_vdgs/generate_vdgs/fragment_database_ligs.py`](../ligand_vdgs/generate_vdgs/fragment_database_ligs.py). This script outputs `database_frags_dict.pkl`, which enumerates all qualifying chemical fragments across your ligand set.

### Option A: RCSB PDB ligands (default)

The wwPDB Chemical Component Dictionary (CCD) file is already provided at `resources/Components-smiles-cactvs.smi` and is the default input. Run with no arguments:

```bash
python ligand_vdgs/generate_vdgs/fragment_database_ligs.py
```

Output: `resources/database_frags_dict.pkl`

### Option B: Custom database ligands

Prepare a **tab-delimited** file (`.smi` or `.tsv`) with 2–3 columns:

| Column | Required | Content |
|--------|----------|---------|
| 1 | yes | SMILES string |
| 2 | yes | Short ligand identifier (e.g. your internal compound ID) |
| 3 | no | Ligand name — not used by the script |

**The identifier in column 2 must match the residue name used in your PDB files**, as it is propagated into the fragment library and used during hit-finding.

Then run:

```bash
python ligand_vdgs/generate_vdgs/fragment_database_ligs.py \
    --ccd <path/to/your_ligands.tsv> \
    --outdir <output_dir>
```

The output `<output_dir>/database_frags_dict.pkl` is used in the next step.

## Step 5: Generate the vdG Database

This step requires `database_frags_dict.pkl` from Step 4, your protonated PDB database from Step 2 (use the trimmed version from Step 1 if you ran it), and the Probe output directory from Step 3.

A qualifying fragment is one present in at least 80 database ligands and containing at most 5 heavy atoms. The core operation is running [`ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py`](../ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py) once per qualifying fragment SMILES. This calls the `vdG-miner` package to extract and cluster vdGs for that fragment; full usage is in the script header.

```bash
python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py \
    -s "<SMILES>" \
    -p <path/to/pdb_database/> \
    -b <path/to/probe_output/> \
    -o <path/to/vdg_library/> \
    --num-procs <n>
```

`-c` is optional — use it to customize the output subdirectory name (see flag table below). Defaults to `-s`.

| Flag | Description |
|------|-------------|
| `-s` | Fragment SMILES, interpreted as a SMARTS pattern for substructure matching |
| `-c` | Chemical group label used as the output subdirectory name under `-o`. By convention, pass the same SMILES as `-s`. Defaults to `-s` if omitted. |
| `-p` | Path to your protonated PDB database (Step 2). Use the trimmed version from Step 1 if you ran that optional step. |
| `-b` | Path to the Probe output directory (output of Step 3) |
| `-o` | Root output directory for the vdG library |
| `--num-procs` | Number of parallel processes per job |

Since the wrapper must be run once per qualifying fragment, use your cluster's scheduler to parallelize — see below.

### SGE (e.g., Wynton)

[`make_sge_scripts_for_frags.py`](../ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py) extracts qualifying fragments from `database_frags_dict.pkl` and writes one ready-to-submit SGE script per fragment. Default paths are set for Wynton; non-Wynton users must pass `--vdg-lib-dir`, `--pdb-dir`, `--probe-dir`, and `--log-dir` explicitly:

```bash
python ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py
```

Then submit all generated scripts:

```bash
for script in ligand_vdgs/generate_vdgs/frag_submit_scripts/*.sh; do
    qsub "$script"
done
```

> **Note:** generated scripts use `conda activate py3.10`. If your conda environment is named differently or you use `module load`, edit [`resources/frag_sge_template.sh`](../resources/frag_sge_template.sh) before running `make_sge_scripts_for_frags.py`.

### SLURM and other schedulers

First extract the qualifying fragment SMILES to a text file (one SMILES per line) using [`extract_fragment_smiles.py`](../ligand_vdgs/generate_vdgs/extract_fragment_smiles.py). Both flags default to paths under `resources/`, so if you used Option A in Step 4, you can run with no arguments:

```bash
python ligand_vdgs/generate_vdgs/extract_fragment_smiles.py \
    --frags-dict <path/to/database_frags_dict.pkl> \
    --output     <path/to/fragment_smiles.txt>
```

Then submit one job per line, calling `vdg_generation_wrapper.py` with the flags shown above. [`resources/frag_sge_template.sh`](../resources/frag_sge_template.sh) shows how the wrapper is invoked inside a job — use it as a reference when writing your own submission script.

### Expected output

Once all jobs complete, the vdG library directory (`-o`) contains one subdirectory per fragment SMILES, each with non-redundant vdGs stored as NPZ files:

```
<vdg_library>/
  <SMILES>/
    nr_vdgs/
      1/          # single-residue vdGs
        <aa_bucket>.npz
      2/          # two-residue vdGs
        <aa_bucket>.npz
```

`aa_bucket` encodes the amino acid composition of the interacting residues (e.g., `ASP_bb.npz`). These files are the input to the hit-finding step.

TODO: summarize in top-level README
