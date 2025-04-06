#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/sge_logs
#$ -N CF3
#$ -cwd
#$ -j y               # tells system STDERR and STDOUT should be joined
#$ -l h_rt=100:00:00   #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes             #-- SGE host reservation
#$ -l hostname='!(qb3-as4|qb3-id188|qb3-id340|qb3-id225)'
#$ -l mem_free=10G
#$ -pe smp 10          # Request 10 slots in the SMP parallel environment

date
hostname

conda activate vdgs

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "[C,c,O,o,N,n,S,s,P,p]C(F)(F)(F)" -c CF3 -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/vdg_lib --symmetry-classes 0 1 2 2 2

date
