#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/sge_logs
#$ -N sat_morpholine_NH
#$ -cwd
#$ -j y               # tells system STDERR and STDOUT should be joined
#$ -l h_rt=240:00:00  #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes             #-- SGE host reservation
#$ -l h=!qb3-as4
#$ -l mem_free=50G

date
hostname

conda activate smart_vdms_env

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "[OX2]1[CX4][CX4][NH][CX4][CX4]1" -c sat_morpholine_NH -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/vdg_lib --symmetry-classes 0 1 2 3 2 1

date
