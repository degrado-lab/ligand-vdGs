# Database Setup Guide

This guide explains how to generate a van der Graph (vdG) database for small molecule functional groups. If you prefer to use a pre-generated vdG database, you may skip the instructions below and download our database [here](insert-url-here).

## Setting Up a Parent Database
To generate your own vdG database, you must first set up a parent database from which to extract vdGs. You have two options:

### Option 1: Setting Up a Custom Parent Database for vdG Extraction
If you have your own custom PDB database that you want to use as the parent database, follow these instructions to set it up for vdG extraction.
- Format your database so that it is in the same database structure as the RCSB PDB mirror, where each structure is represented by 4 characters (XXXX.pdb) and sorted into subdirectories based on the inner 2 characters. For example, a PDB file named `1ABC.pdb` will be placed into a subdirectory named `AB`, resulting in the file path `AB/1ABC.pdb`. Use `ligand_vdgs/tools/format_parent_database.py`.
- Ensure that your PDBs contain hydrogens. If they don't, we recommend you use [Reduce2](https://github.com/cctbx/cctbx_project/tree/master/mmtbx/reduce) (open-source) or Schrodinger's prepwizard program (more accurate, but requires a Schrodinger license). If you would like to run prepwizard on multiple structures in parallel, you can use `ligand_vdgs/preprocessing/protonate_pdbs.py`.
- ... placeholder for Probe instructions ...

### Option 2: Setting Up the RCSB PDB as a Parent Database for vdG Extraction
If you want to use the RCSB PDB as your parent database, follow these instructions to extract vdGs from the RCSB PDB.
- [placeholder]

## Generating a vdG Database from the Parent Database

`ligand-vdGs` can be used to generate vdGs of SMARTS-defined chemical groups (CGs). Currently, the pipeline is to run (1) `ligand-vdGs/ligand_vdgs/programs/vdG-miner/vdg_miner/programs/smarts_to_cgs.py` to identify the ligands and PDBs that contain the CGs represented by the specified SMARTS strings, (2) `ligand-vdGs/ligand_vdgs/programs/vdG-miner/vdg_miner/programs/generate_fingerprints.py` to generate fingerprints that describe the environment of the CGs, and finally (3) `ligand-vdGs/ligand_vdgs/programs/vdG-miner/vdg_miner/programs/fingerprints_to_pdbs.py` to output PDB files of the resulting vdGs. A current streamlined pipeline is being developed.