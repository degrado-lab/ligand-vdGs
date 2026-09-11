# Database Setup Guide

This guide explains how to generate a van der Graph (vdG) database for small molecule functional groups. If you prefer to use a pre-generated vdG database, you may skip the instructions below and download our database [here](insert-url-here).

> **All commands in this guide assume you are running from the root directory of the `ligand-vdGs` package.**

## Installation

1. **Set up the Python environment.** A conda environment file is provided at `environment.yml`:

   ```bash
   conda env create -f environment.yml
   conda activate lig_vdgs
   ```

2. **Run scripts by path, from the repository root.**

   `environment.yml` also installs the repo itself in editable mode
   (`pip install -e .`), which is what makes `ligand_vdgs.*` imports resolve; run
   that command manually from the repository root if you set the environment up
   some other way. `environment.yml` remains the single source of truth for
   dependencies.

> **Note on OpenBabel:** The vdG-miner component requires OpenBabel *with Python bindings* (`from openbabel import openbabel`), not just the command-line tool.

## Prerequisites: Set Up a Parent Database

Before running any pipeline steps, you need a preprocessed PDB database to extract vdGs from. You can use any collection of PDB structures — your own custom set or a mirror of the RCSB PDB.

**Directory layout.** Structures must be organized in RCSB mirror format: each file is named `XXXX.pdb` (4-character code) or `XXXX_N.pdb` (biounit assembly N) and placed in a subdirectory named after the inner 2 characters of the code, **lowercased**: `1ABC.pdb` → `ab/1ABC.pdb`. Every reader (mining and Probe output alike) lowercases those two characters when it looks a structure up. Use [`ligand-vdGs/scripts/format_parent_database.py`](../scripts/format_parent_database.py) to reformat an existing directory; it copies (the source is left intact) and validates every file name before touching anything.

> **Biounits:** the pipeline reads `XXXX_N.pdb` fine — `_pdb_id_from_path` splits the assembly suffix off, the npz stores the biounit stem, and Probe output is named by the full stem. Only `format_parent_database.py` is stricter (4 characters exactly), so lay biounit mirrors out by hand.

Once your source PDBs are in the right layout, complete Steps 1–3 in order before moving on to Steps 4–5.

- **Step 1 (`s01_trim_database.py`) — Prune the database (optional, recommended).** Run [`ligand-vdGs/ligand_vdgs/preprocessing/s01_trim_database.py`](../ligand_vdgs/preprocessing/s01_trim_database.py) to filter by ligand b-factor and extract 20 Å binding-site regions (the `radius` constant; 20 Å is the smallest sphere that keeps a ligand's *third* coordination shell, which is what Step 2's protonation needs — 10 Å kept only 8% of it), reducing database size and removing redundant structures. It has no CLI: the input/output directories, b-factor cutoff and `skip_to_output_pdbs` switch are module-level constants at the top of the script. The dedup pass writes a database JSON that the trimming pass then reads, so a fresh database needs `skip_to_output_pdbs = False` on the first run. The pipeline works without this step, but skipping it means running on full PDB files (slower, larger, and noisier). If you run it, use the output as `-p` in Step 5 below; otherwise use the formatted database from above.

  Before parsing each file, this step rewrites two-character chain IDs (columns 21-22) to unique single characters and turns peptide-bonded HETATM amino acids into ATOM records (`_chain_ids.remap_chain_ids`, `_prep_filters.embedded_amino_acid_hetatm_to_atom`); the chain mapping is kept as `REMARK 900` lines in the output. Both must happen on the text, before ProDy sees it -- see "Chain IDs are one column" in [`pitfalls.md`](pitfalls.md). A database trimmed before this fix can be repaired with `scripts/remap_chain_ids.py` (into a new directory; re-run Step 3 on the repaired files).

  It also handles modified residues, which prepwizard used to rename to their parent amino acid (`_prep_filters.modified_residues_to_protein`): MSE/SEC become MET/CYS (keeping `SE`), and every other chain-bonded HETATM residue whose CCD type is `*PEPTIDE LINKING` (KCX, SEP, LLP, ...) becomes an ATOM record under its own resname, so it is neither mined as a ligand nor accepted as a slot; covalent cofactors (`NON-POLYMER`: PLP, HEM, SAM) stay HETATM ligands. This needs `resources/ccd_polymer_types.tsv`; regenerate it on a login node with `python scripts/fetch_ccd_polymer_types.py` if your database uses CCD entries newer than the table (s01 prints the resnames it did not find). Without the table s01 warns and only the renames happen. See "Modified residues" in [`pitfalls.md`](pitfalls.md).

  This step also drops the residues prepwizard mangles (`_prep_filters.drop_prepwizard_hazard_residues`): amino acid residues missing any of N, CA or C, and one-atom HET residues whose single heavy atom is carbon. Step 2 applies the same filter, so it holds even if you skip Step 1. Prepwizard cannot build an amino acid from a lone backbone N, and re-emits that orphan atom under a *ligand's* resname, chain and resnum, giving that ligand a stray atom tens of Å from the rest of the residue. See "Prepwizard relabels atoms it cannot build as part of a ligand" in [`pitfalls.md`](pitfalls.md).

- **Step 2 (`s02_run_prepwizard.sh`) — Add hydrogens.** The pipeline requires hydrogens to be present in every PDB. This step re-applies the Step 1 filter (so it holds if you skipped Step 1) and, afterwards, restores any free ligand PrepWizard renamed into a protein residue — see [`pitfalls.md`](pitfalls.md). Use [Reduce2](https://github.com/cctbx/cctbx_project/tree/master/mmtbx/reduce) (open-source) or Schrödinger's PrepWizard (more accurate, requires a license). [`s02_run_prepwizard.sh`](../ligand_vdgs/preprocessing/s02_run_prepwizard.sh) is an array job that runs `_protonate_pdbs.py` in batches; the input/output directories and prepwizard path are module-level constants in `_protonate_pdbs.py` (its only CLI flag is `--batch_index`). **The batch count lives in two independent places and must be kept in sync by hand:** the array task count is positional `$1` to the shell script, while `num_batches_total = 15` is a constant in `_protonate_pdbs.py`. If they disagree the database is silently mis-partitioned — tasks either overlap or leave subdirectories unprotonated — so pass the same number:

```bash
qsub -t 1-15 ligand_vdgs/preprocessing/s02_run_prepwizard.sh 15
```
 Run this on the output of Step 1, or the formatted database above if you skipped trimming.

- **Step 3 (`s03_run_probe.sh`) — Run Probe.** Probe must be installed separately; it is available from the [Richardson lab / MolProbity project](https://github.com/rlabduke/probe). Submit all PDBs as an SGE array job with [`ligand-vdGs/ligand_vdgs/preprocessing/s03_run_probe.sh`](../ligand_vdgs/preprocessing/s03_run_probe.sh), which is the entry point; it invokes the private helper `_run_probe.py` per structure, resolved next to the script, so submit from anywhere. Output lands in `<probe_dir>/<xx>/<stem>.probe.gz` in the same lowercased mirror layout as the PDBs, and that directory is `-b` in Step 5 below. It takes four positionals and one task per PDB, so size `-t` to the file count (`$3`, the package path, is accepted but unused):

```bash
N=$(find <path/to/pdb_database/> -type f -name '*.pdb' | wc -l)
qsub -t 1-$N -o <path/to/logs/> -e <path/to/logs/> \
    ligand_vdgs/preprocessing/s03_run_probe.sh \
    <path/to/pdb_database/> <path/to/probe_output/> \
    <path/to/SMARTS-vdg_package/> <path/to/probe_binary>
```


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

This step requires your protonated PDB database from Step 2 (use the trimmed version from Step 1 if you ran it) and the Probe output directory from Step 3. `database_frags_dict.pkl` from Step 4 is read only by the job schedulers below, which decide *which* fragments to run.

A qualifying fragment is one that:

- carries no element from `UNDESIRED_ELEMENTS` (`fragment_database_ligs.py`: metals, lanthanides, noble gases, plus Si/Se/As/Te) and passes `Frags.is_organic` — this is applied at Step 4, so a disqualified element never reaches the fragment dict. **Boron is deliberately not excluded**: boronic acids and boronate esters are real covalent warheads and worth mining.
- has at most 5 heavy atoms,
- is not a halogen oxyanion (a crystallization salt, not a binding moiety), and
- has at least `--min-instances` (default 250) estimated CG occurrences in the parent database — SMARTS matches summed over every ligand copy in the mirror, i.e. candidate vdG sites, not CCD ligand counts and not structure counts — after protonation-state variants are pooled onto one representative (`select_fragments` in `extract_fragment_smiles.py`).

The occurrence count comes from the same sampling pass as the per-fragment cost estimate (`estimate_frag_cost.estimate_fragment_counts`), which also returns the structure count that sizes each job's SGE slots. The core operation is running [`ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py`](../ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py) once per qualifying fragment SMILES. This calls the `vdG-miner` package to extract and cluster vdGs for that fragment; full usage is in the script header.

```bash
python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py \
    -s "<SMILES>" \
    -p <path/to/pdb_database/> \
    -b <path/to/probe_output/> \
    -o <path/to/vdg_library/> \
    --num-procs <n>
```

`-c` is optional — use it to customize the output subdirectory name (see flag table below). Defaults to `-s`. The wrapper encodes it with `utils.smiles_to_filename` (`/` → `_fs_`, `\` → `_bs_`) and passes the *encoded* label to every subprocess, so a SMILES containing `/` (e.g. `C/C=C/O`) works whether given explicitly or reached via the `-s` default.

| Flag | Description |
|------|-------------|
| `-s` | Fragment SMILES, interpreted as a SMARTS pattern for substructure matching |
| `-c` | Chemical group label. Encoded with `utils.smiles_to_filename` and used as the output subdirectory name under `-o`; `/` and `\` are handled for you. Defaults to `-s` if omitted. |
| `-p` | Path to your protonated PDB database (Step 2). Use the trimmed version from Step 1 if you ran that optional step. |
| `-b` | Path to the Probe output directory (output of Step 3) |
| `-o` | Root output directory for the vdG library |
| `--num-procs` | Number of parallel processes per job |
| `--subset-sizes` | Which vdG subset sizes to generate; only `1` and `2` are accepted (default: `1 2`) |
| `--no-profile-compute` | Skip the per-phase timing sidecar (`<cg_label>_compute_profile.json`), which is written by default |
| `-m`, `--max-num-vdgs-to-clus` | **Debugging only.** Cap on distinct PDB IDs per AA bucket (one PDB can still contribute several vdGs), to bound the RMSD step for a quick test run. Never set it for a production build: a capped run silently yields an incomplete library. |

The requested subset sizes share one environment-reconstruction pass and are
clustered independently into `nr_vdgs/<subset_size>/`.

Clustering runs in two stages. **Stage 1** groups vdGs on pose — CG atoms plus
each vdM's N/CA/C — by sphere-exclusion clustering (Butina 1999; equivalently the
GROMOS algorithm of Daura et al. 1999): the exact within-cutoff neighbour graph is
built, then the unassigned vdG with the most neighbours is repeatedly taken as a
cluster seed along with its unassigned neighbours. Every member is therefore within
the cutoff of the vdG that represents it, and the partition does not depend on
input order. The distance is minimised over the fragment's CG automorphisms *and*
over orderings of interchangeable same-label vdM slots, so one physical environment
is one record. **Stage 2** subdivides each pose cluster by flanking-sequence and
flanking-CA similarity; each subgroup stores its own pose-minimax member, and the
exact pose radius of that choice is recorded in `cluster_pose_radius`.

Fragment strings are substructure queries, not standalone molecules. They are
parsed directly with RDKit `MolFromSmarts`, without SMILES sanitization, valence
checking, or hydrogen inference. Keep fragment definitions free of explicit H
atom nodes; hydrogens in the prepared PDB database for Probe are a separate
concern.

Exact automorphisms cover ordinary graph symmetry for every atom type. On top of
that, resonance groups are normalized so that their drawn bond order and charge do
not split symmetric atoms: amidine/guanidine C-N, N-O, S-N and conjugated N-N
groups, plus charge/H on nitrogens within one aromatic component. The main rule is
the terminal-atom one, which ignores untrusted bond order, formal charge, and
proton placement: two or more terminal N, O, or S atoms on a shared
B/C/N/O/P/S/Cl/Br/I center are treated as one set of interchangeable positions.
Thus the two carboxylate oxygens, the two neutral carboxylic-acid oxygens, terminal
phosphate `P=O`/`O-`/`OH` positions, sulfinic-acid oxygens, and dithioacid
sulfurs can exchange. Each element forms its own set, so a terminal O is never
mapped onto a terminal S (in `OP(O)(=S)[S-]` the O pair and S pair permute
independently). Substituted and bridging atoms remain distinguished by
heavy-atom connectivity, and aromatic terminal atoms are excluded because a
truncated ring atom is a real ring position rather than a resonance form. One
carve-out trusts the drawing: a center with four or more terminal atoms drawn
with all-single bonds but mixed charges keeps its charges, so
`[O-]S([O-])([O-])O` permutes its three anionic oxygens and leaves the hydroxyl
fixed.

CG coordinates and atom names are stored in the fragment SMARTS's slot order;
generation does not relabel each record from its geometry. That order is an
indexing convention, not a unique correspondence for symmetric atoms. Clustering
and every atom-correspondence-dependent downstream comparison must minimize over
the complete recorded automorphism group.

Fragment selection also discards halogen oxyanions (perchlorate, chlorate,
periodate, bromate). These are crystallization and cryoprotectant salts, not
ligand chemistry. A generated dictionary never contains the free ions, because
carbon-free ligands are dropped before fragmentation; the filter catches their
esters and any hand-written work list.

Since the wrapper must be run once per qualifying fragment, use your cluster's scheduler to parallelize — see below.

### SGE (e.g., Wynton)

[`make_sge_scripts_for_frags.py`](../ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py) extracts qualifying fragments from `database_frags_dict.pkl` and writes one ready-to-submit SGE script per fragment. It samples `--pdb-dir` to estimate each fragment's cost (or reads `--frag-cost-estimate`) and tiers slots and `h_rt` per fragment, clamped by `--max-h-rt`, which is required. It also writes `<--vdg-lib-dir>/fragment_aliases.tsv`, the map from each collapsed protonation variant to the representative that was mined. The cheapest tier (< 750 estimated structures) requests `h_rt 0:29:00` so the job qualifies for Wynton's short queue and starts far sooner; `--no-short-queue` turns that off. A short-tier job that overruns is killed with a partial fragment directory, so after a build re-run any fragment whose `<cg_label>/<cg_label>_log` lacks `Job completed.` — `rm -rf` that fragment directory first, since the wrapper refuses a non-empty one and has no resume flag (see "Expected output" below). It refuses to run if `--sge-out-dir` already has files. Default paths are set for Wynton. Non-Wynton users must pass `--vdg-lib-dir`, `--pdb-dir`, `--probe-dir`, and `--log-dir` explicitly:

```bash
python ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py \
    --max-h-rt <HH:MM:SS> \
    --vdg-lib-dir  <path/to/vdg_library/> \
    --pdb-dir      <path/to/pdb_database/> \
    --probe-dir    <path/to/probe_output/> \
    --log-dir      <path/to/logs/> \
    --sge-out-dir  <path/to/empty/script_dir/>
```

Then submit all generated scripts:

```bash
for script in <path/to/empty/script_dir/>/*.sh; do
    qsub "$script"
done
```

[`run_production_frags.sh`](../run_production_frags.sh) at the repository root is the driver that does all of the above in one step, and is the only path that gets the submission order right: it clears `--sge-out-dir`, regenerates the fleet, and qsubs **longest estimated job first**, because SGE will not schedule a job whose `h_rt` reaches past the maintenance boundary, so the runway for a long request shrinks by an hour every hour. `--mode` is required (`threshold-plus-include` for a full build, `include-only` for a top-up); `--dry-run` generates and prints the order without submitting. Because the generated scripts use `#$ -cwd` with repo-relative paths, submit from the repository root either way.

```bash
./run_production_frags.sh --mode threshold-plus-include --dry-run   # inspect first
MAX_H_RT=36:00:00 ./run_production_frags.sh --mode threshold-plus-include
```

It requires `resources/frag_cost_estimate.tsv`; regenerate it with [`estimate_frag_cost.py`](../ligand_vdgs/generate_vdgs/estimate_frag_cost.py) whenever `database_frags_dict.pkl` changes (the TSV records the dict's path and SHA-256, and `make_sge_scripts_for_frags.py` warns — it does not refuse — when that hash no longer matches `--frags-dict`):

```bash
python ligand_vdgs/generate_vdgs/estimate_frag_cost.py \
    --pdb-dir <path/to/pdb_database/> \
    --output  resources/frag_cost_estimate.tsv
```

> **Counts come from `--pdb-dir`, not from the CCD.** The estimate is a SMARTS
> pass over a sample of the parent PDB mirror, so a fragment whose ligands are
> not in that mirror estimates 0 occurrences and is dropped by
> `--min-instances` with no job written. When adding novel ligands beyond the
> CCD, get their structures into the mirror you pass as `--pdb-dir` and
> regenerate the TSV before generating the fleet; force in a wanted fragment
> that still scores low with `--include`/`--include-file`.

### SLURM and other schedulers

Use [`extract_fragment_smiles.py`](../ligand_vdgs/generate_vdgs/extract_fragment_smiles.py)
to write a scheduler-agnostic one-column work list. It applies the same selection
rule as the SGE path—in fact both call the same `select_fragments()`—and emits one
fragment SMARTS per line, plus the alias map next to it as
`<output stem>_aliases.tsv` (columns `alias`, `representative`, `kind`; `kind=promoted`
marks a representative that is the charge-stripped SMARTS of its aliases and is not
itself a frags-dict key — trace it with `scripts/lookup_fragment_key.py`; see
docs/pitfalls.md "Protonation variants are collapsed"):

```bash
python ligand_vdgs/generate_vdgs/extract_fragment_smiles.py \
    --pdb-dir    <path/to/pdb_database/> \
    --frags-dict <path/to/database_frags_dict.pkl> \
    --output     <path/to/fragment_work_list.txt>
```

Only the SMARTS varies per job. Paths (`-p`, `-b`, `-o`) and run settings belong
in the submission script. SGE, Slurm, and bare-shell runs all invoke the same
wrapper, which derives exact automorphisms from that SMARTS and writes them to
`nr_vdgs/cg_symmetry.npz` — no scheduler-provided classes or serialized mappings.
`-c` can be omitted: the wrapper defaults it to the SMARTS and applies the same
`utils.smiles_to_filename` encoding `make_sge_scripts_for_frags.py` does, so both
paths land in the same library directory even for the work list's SMARTS
containing `/` (e.g. `C/C=C/O`):

```bash
while IFS= read -r smarts; do
    args=(-s "$smarts" \
          -p "$PDB_DIR" -b "$PROBE_DIR" -o "$OUT_DIR" \
          --num-procs "$NPROCS" --subset-sizes 1 2)
    # replace with sbatch/qsub/srun or a bare call
    your-submit-command python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py "${args[@]}"
done < fragment_work_list.txt
```

[`resources/frag_sge_template.sh`](../resources/frag_sge_template.sh) shows how the wrapper is invoked inside a real job, if you want a reference for the non-fragment-specific arguments.

### Expected output

Once all jobs complete, the vdG library directory (`-o`) contains one subdirectory per fragment SMILES, each with non-redundant vdGs stored as NPZ files:

```
<vdg_library>/
  fragment_aliases.tsv        # SGE path only: collapsed variant -> representative
  <cg_label>/                 # utils.smiles_to_filename(-c), not the raw SMILES
    <cg_label>_log
    <cg_label>_matches.pkl    # SMARTS match lists, left by smarts_to_cgs.py
    <cg_label>_ligands.sdf
    <cg_label>_compute_profile.json   # unless --no-profile-compute
    nr_vdgs/
      cg_symmetry.npz         # cg_smarts + cg_automorphisms; load_cg_symmetry reads this
      1/                      # single-residue vdGs
        <aa_bucket>.npz
      2/                      # two-residue vdGs
        <aa_bucket>.npz
```

`cg_symmetry.npz` is not optional: it is the only record of which SMARTS the
automorphism group was derived from, and consumers must read the group from there
rather than re-deriving it from the directory name. It also defines the slot order
used by stored CG coordinate and atom-name columns.

`aa_bucket` encodes the slot labels of the interacting residues, sorted and joined with `_` (e.g., `ASP_bb.npz`, `SER_bb.npz`). A label is a resname, the backbone label `bb`, or `X`; none contains `_`, so the file name splits back apart unambiguously. Each fragment directory also contains a `<cg_label>_log` file; downstream tools treat a fragment whose log is missing or incomplete as absent from the library. A fragment directory cannot be re-run in place: the wrapper refuses a non-empty `<-o>/<cg_label>`, and there is no overwrite or resume flag.

#### What a bucket npz contains

Each bucket holds **two disjoint blocks of arrays**, plus per-cluster columns.
Every vdG appears in exactly one block.

| prefix | one row per | holds |
|---|---|---|
| `nr_*` | non-redundant vdG (one per cluster) | CG and vdM-backbone coordinates, CG and vdM identity, slot flag, measured quality |
| `mem_*` | clustered vdG that is *not* the nr vdG | identity only, no coordinates; tied back by `mem_cluster_id` |

A *member* is an observation that is not the nr vdG, so

```
cluster_size == 1 + (number of mem_ rows for that cluster)
```

Per-cluster columns: `cluster_id`, `cluster_size` (raw observations),
`cluster_num_parents` (distinct parent structures — the figure to do statistics
on, since NCS copies and homologous entries inflate `cluster_size`),
`cluster_pose_radius` (greatest symmetry-aware RMSD from any member to the stored
row), and `first_stage_cluster_id` / `second_stage_cluster_id`.

Bucket-level columns: `aa_bucket_parts` (the label of each vdM slot — this, not
`nr_scrr_resname`, is what says which slots may be permuted), `cg_elements`, and
`parent_pdb_dir`.

Member rows carry no coordinates by design: they are re-derived from the parent
PDB on demand (`vdg_npz_utils.rederive_member_coords`). Parent structures are
stored as a biounit stem plus one directory scalar, and
`resolve_parent_pdb_path` rebuilds the path — pass `pdb_dir=` to point a copied
library at a different PDB mirror.

Quality is recorded per row rather than only filtered on: `nr_cg_max_b`,
`nr_cg_min_occ`, `nr_vdm_max_b`, `nr_vdm_min_occ` (and the `mem_` equivalents)
are measured over the atoms that actually enter the vdG — the CG's own atoms and
the contacting residues' heavy atoms. Mining applies only a loose floor, so a
stricter cut is a read-path decision rather than a reason to re-mine.

These files are the input to the hit-finding step.

#### When to run the H-class diagnostic

Run [`h_class_diagnostic.py`](../ligand_vdgs/tools/h_class_diagnostic.py) after every
full library build or rebuild, and again whenever the fragment key scheme changes,
before deciding which `H0`/`!H0` key variants downstream code should pool at read
time. It reads the per-observation H-class fields each bucket stores next to its
row arrays -- `nr_cg_heavy_degree`, `nr_cg_num_h`, `mem_cg_heavy_degree`,
`mem_cg_num_h` (int8, one column per CG atom); `--dry-run` only checks that a
library carries them. Writer contract: an atom that could not be read gets a
**negative** value in both fields; 0 is a real count and must never stand for
"unreadable" (the diagnostic errors on any `heavy_degree == 0`). The writer also
stores `nr_vdm_o_coords` (n_nr, num_vdms, 3), the backbone carbonyl O per vdM slot
in the same frame as `nr_vdm_bb_coords` (NaN where absent), kept out of every RMSD; with it,
`--contact-atoms N,CA,C,O` measures real C···O contacts.

```bash
python ligand_vdgs/tools/h_class_diagnostic.py --lib <vdg_library> --out h_class.tsv
```

## Inspect the Library

[`materialize_vdg_pdbs.py`](../ligand_vdgs/generate_vdgs/materialize_vdg_pdbs.py)
writes library vdGs out as PDB files for visual inspection in PyMOL. Two flags are
required: `-c/--cg-nr-vdgs-root` (a fragment's `nr_vdgs/` directory, not the library
root) and `-o/--out-dir`, which must be empty or absent — `fresh_dir` refuses a
populated one, so a partial earlier run has to be cleared by hand.

```bash
python ligand_vdgs/generate_vdgs/materialize_vdg_pdbs.py \
    -c <path/to/vdg_library/><cg_label>/nr_vdgs/ \
    -o <path/to/empty/output_dir/> \
    --top-clusters 10
```

It has two modes, and one selection vocabulary shared by both: `--aa-buckets` and
`--subset-sizes` say where to look, `--top-clusters N` or `--clusters ID...` say which
clusters. The default writes **nr vdGs** — one frame per selected cluster's stored row,
with `--min-cluster-size` dropping the sparse tail. Adding `--members` switches to
**members** mode, which writes each selected cluster's nr vdG *and* the individual
vdGs inside it, one directory per cluster, and requires `--top-clusters` or
`--clusters`; it is the only consumer of the `mem_*` arrays in the npz;
`--reps` caps how many structures are written per cluster **including the nr vdG**,
so `--reps 20` is the nr vdG plus 19 randomly sampled members.

### Output file names

Names are built by joining fields with a **single `_`**, so they can be parsed
programmatically. Both this script and
[`write_vdg_hit_pdbs.py`](../ligand_vdgs/tools/write_vdg_hit_pdbs.py) share the naming
code in `functions/vdg_pdb_io.py`.

```
nr vdGs (default) — every file directly in <out_dir>
    O=S(=O)(O)O_ASP_bb_clus1_size4_2ldb__A_46_ASP__A_43_GLY__A_4_SO4.pdb.gz
    └─ frag ──┘ └bucket┘ └cluster┘ └src┘ └vdM tags──────────┘└lig tag─┘

members (--members)
  <out_dir>/<subset_size>/<AA_BUCKET>/clus<id>_size<size>/
    O=S(=O)(O)O_NR_clus1_size6_1yp2__D_370_ASP__D_2003_SO4.pdb.gz
    O=S(=O)(O)O_1yp3__C_370_ASP__C_1002_SO4.pdb.gz
```

A **residue tag** is always exactly four `_`-separated fields,
`seg_chain_resnum_resname`. The segment field is kept even when the entry has no
segment, which is the usual case — hence the leading empty field that makes
`_D_370_ASP` look like it starts with a doubled underscore. Keeping it is what fixes
the field count; a field that would otherwise contain `_` (a stray `/`, `\`, or space)
becomes `-` for the same reason.

Every name carries exactly **`1 + subset_size` tags — one per vdM slot in slot
order, then the ligand last**:

| subset size | tags |
|---|---|
| 1 | `<res1>_<lig>` |
| 2 | `<res1>_<res2>_<lig>` |

A vdG is one CG cut out of **one** ligand residue plus its vdM residues, so there is
always exactly one ligand tag. The npz stores the CG's `seg`/`chain`/`resnum`/
`resname` once per record for exactly that reason, and the invariant is enforced at
mining time: `vdg_struct_utils.get_cg_atoms` reports and skips a CG whose atoms span
more than one residue, since that would mean fragment matching crossed a residue
boundary. Measured before the schema change, it never happened: 0 of 5,973,427
nr vdG records across a full 349-fragment library.

**Parse from the right.** The head of the name is not positionally parseable: `<frag>`
is a sanitized SMILES (`smiles_to_filename` encodes `/` as `_fs_`) and `<AA_BUCKET>`
joins its labels with `_`, so both contribute a variable number of fields. Before
splitting, strip:

| strip | when |
|---|---|
| `.pdb.gz` | always |
| `~<n>` | when two files would otherwise share a name — see below |

Then take the last `4 × (subset_size + 1)` fields and cut them into groups of 4:

```python
stem   = re.sub(r"~\d+$", "", basename.removesuffix(".pdb.gz"))
fields = stem.split("_")[-4 * (subset_size + 1):]
tags   = ["_".join(fields[i:i + 4]) for i in range(0, len(fields), 4)]
vdms, lig = tags[:-1], tags[-1]
```

**Getting `subset_size`.** Members mode puts it in the directory path. A nr vdGs
run has no subset-size directory — that read as a second "size" next to the
cluster's `_size<n>` — so recover it from the AA bucket label, which sits between
the fragment and `_clus<id>` and has one token per vdM slot (`ASP` → 1, `ASP_bb` → 2).
You know `<frag>`, since you ran the script against one fragment's `nr_vdgs`:

```python
rest   = stem[len(frag) + 1:]                       # drop "<frag>_"
bucket = re.match(r"(.+?)_clus\d+_size\d+_", rest).group(1)
subset_size = len(bucket.split("_"))
```

Don't try to infer it by counting fields from the right instead: `<source>` is a
biounit stem that may itself contain `_` (`1f8s` vs `1f8s_1`), so the field count
after the cluster anchor is not a fixed function of subset size.

Hit files from `write_vdg_hit_pdbs.py` follow the same rule with `subset_size` tags
and no ligand tag — the CG there is the query's own, and the query is already named
at the front of the file. Their tags are **not** the last fields, though: the name
ends `..._<vdM tags>_<rmsd>`, so drop one trailing field before counting from the
right. If two files in one run would get the same name, they are numbered `~1`, `~2`, ….
