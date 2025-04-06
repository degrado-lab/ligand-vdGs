#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/sge_logs/
#$ -N aryl_hydroxyl
#$ -cwd
#$ -j y                 # tells system STDERR and STDOUT should be joined
#$ -l h_rt=300:00:00     #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes               #-- SGE host reservation
#$ -l hostname='!(qb3-as4|qb3-id188|qb3-id340|qb3-id225)'
#$ -l mem_free=30G
#$ -pe smp 10          # Request 10 slots in the SMP parallel environment
            

date
hostname

conda activate vdgs

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "c[OX2H1]" -c aryl_hydroxyl -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/vdg_lib 

qstat -j $JOB_ID

date
