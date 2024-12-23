# ligand-vdGs

## Description
This package builds upon https://github.com/degrado-lab/vdG-miner by Rian Kormos for generating a [vdG](hyperlink_to_eventual_preprint) library. The vdG library is then used for docking ligands into known binding sites.

## Database Generation
For prerequisites and detailed instructions on obtaining a vdG database, please refer to the [Database Generation Guide](docs/database_generation_guide.md) located in `vdG-ligands/docs/`.


## Quick Start
### Set up the conda environment for this package.
From the `vdG-ligands` top-level directory, run 
```bash
conda env create -f environment.yml
conda activate vdg_env
```
You only need to create the environment once, but the environment needs to be active before running the scripts contained in this package. In most cases, that requires running the `conda activate` command whenever you open a new terminal/window.


[placeholder for instructions on installing ligand-vdGs as pip package]
