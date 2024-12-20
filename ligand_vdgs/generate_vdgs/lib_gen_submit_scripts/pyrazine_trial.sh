#!/bin/bash

#$ -S /bin/bash
#$ -o $YOUR_LOG_DIR
#$ -N pyrazine_trial
#$ -cwd
#$ -j y               # tells system STDERR and STDOUT should be joined
#$ -l h_rt=5:00:00    #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes             #-- SGE host reservation

date
hostname

conda activate smart_vdms_env

python ligand_vdgs/scripts/vdg_generation_wrapper.py -s "[nX2]1cc[nX2]cc1" -c pyrazine -p $YOUR_PDB_DIR/ -b $YOUR_PROBE_FILES -o $YOUR_OUTPUT_DIR

date
