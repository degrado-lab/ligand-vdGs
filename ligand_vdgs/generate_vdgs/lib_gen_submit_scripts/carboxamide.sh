#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/sge_logs
#$ -N carboxamide 
#$ -cwd
#$ -j y                 # tells system STDERR and STDOUT should be joined
#$ -l h_rt=24:00:00     #-- runtime limit - max 2 weeks == 336 hours
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -R yes               #-- SGE host reservation

date
hostname

conda activate smart_vdms_env

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "CC(=O)[NH2]" -c carboxamide -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/vdg_lib

date
