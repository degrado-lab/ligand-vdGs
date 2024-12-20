# Database Setup Guide

This guide explains how to generate a van der Graph (vdG) database for small molecule functional groups. If you prefer to use a pre-generated vdG database, you may skip the instructions below and download our database [here](insert-url-here).

## Setting Up a Parent Database
To generate your own vdG database, you must first set up a parent database from which to extract vdGs. You have two options:

### Option A: Use a Custom Parent Database
If you have your own custom PDB database that you want to use as the parent database, follow these instructions to set it up for vdG extraction.
    
1. Format your database so that it is in the same database structure as the RCSB PDB mirror, where each structure is represented by 4 characters (XXXX.pdb) and sorted into subdirectories based on the inner 2 characters. For example, a PDB file named `1ABC.pdb` will be placed into a subdirectory named `AB`, resulting in the file path `AB/1ABC.pdb`. Use [`ligand-vdGs/ligand_vdgs/tools/format_parent_database.py`](../ligand_vdgs/tools/format_parent_database.py).

2. Prune your parent database by running [`ligand-vdGs/ligand_vdgs/preprocessing/s01_trim_database.py`](../ligand_vdgs/preprocessing/s01_trim_database.py), which filters ligand b-factors and isolates binding sites for downstream processing.

3. Ensure that your PDBs contain hydrogens. If they don't, we recommend you use [Reduce2](https://github.com/cctbx/cctbx_project/tree/master/mmtbx/reduce) (open-source) or Schrodinger's prepwizard program (more accurate, but requires a Schrodinger license). If you would like to run prepwizard jobs in parallel on a slurm cluster, you can use [`ligand-vdGs/ligand_vdgs/preprocessing/s02_run_prepwizard.sh`](../ligand_vdgs/preprocessing/s02_run_prepwizard.sh).

4. Run Probe on your database PDBs using [`ligand-vdGs/ligand_vdgs/preprocessing/run_probe.py`](../ligand_vdgs/preprocessing/run_probe.py) to define all contacts made by a ligand. To submit all PDBs in parallel on an SGE cluster, you can use [`ligand-vdGs/ligand_vdgs/preprocessing/s03_submit_probe.sh`](../ligand_vdgs/preprocessing/s03_submit_probe.sh).

### Option B: Use the RCSB PDB as a Parent Database for vdG Extraction
If you want to use the RCSB PDB as your parent database, follow these instructions to extract vdGs from the RCSB PDB.
- [placeholder]

## Generating a vdG Database from the Parent Database

Now that you've prepared your parent database, you can generate vdGs of SMARTS-defined chemical groups (CGs) using [`ligand-vdGs/ligand_vdgs/scripts/vdg_generation_wrapper.py`](../ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py), which calls on scripts within the `vdG-miner` package. Usage instructions are in the header, and all submission scripts for the CGs in the distributed vdG library can be found in [`ligand-vdGs/ligand_vdgs/generate_vdgs/lib_gen_submit_scripts`](../ligand_vdgs/generate_vdgs/lib_gen_submit_scripts/).

[placeholder for deduplication step]