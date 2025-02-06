#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/sge_logs/
#$ -N sulfonamide_nh2
#$ -cwd
#$ -j y         # tells system STDERR and STDOUT should be joined
#$ -l h_rt=24:00:00   #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes               #-- SGE host reservation

date
hostname

conda activate smart_vdms_env

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "[NX3&H2][SD4](C)(=O)(=O)" -c sulfonamide_nh2 -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/vdg_lib --symmetry-classes 0 1 2 3 3

date
