#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/sge_logs/
#$ -N pyrazine800
#$ -cwd
#$ -j y                 # tells system STDERR and STDOUT should be joined
#$ -l h_rt=300:00:00    #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes               #-- SGE host reservation
#$ -l mem_free=100G
#$ -pe smp 5            # Request 5 slots in the SMP parallel environment
            

date
hostname

conda activate vdgs

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "[nX2,NX2,NX3H0,nX3H0]~1~c~c~[nX2,NX2,NX3H0,nX3H0]~c~c1" -c pyrazine800 -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/vdg_lib --symmetry-classes 0 1 1 0 1 1 -m 800

date
echo "DONE"
