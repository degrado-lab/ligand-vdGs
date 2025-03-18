#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/sge_logs
#$ -N pyridone_NH
#$ -cwd
#$ -j y                  # tells system STDERR and STDOUT should be joined
#$ -l h_rt=300:00:00     #-- runtime limit - max 2 weeks == 336 hours
#$ -l mem_free=100G
#$ -R yes                #-- SGE host reservation
#$ -l h=!qb3-as4

date
hostname

conda activate smart_vdms_env

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "[NH,nH]~1~[CX3,cX3]~[CX3,cX3]~[CX3,cX3]~[CX3,cX3]~[CX3,cX3]1(=O)" -c pyridone_NH -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/vdg_lib 

date
