#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/frag_sge_logs 
#$ -N "CNC(=N)N" 
#$ -cwd
#$ -j y                #-- tells system STDERR and STDOUT should be joined
#$ -l h_rt=300:00:00   #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes              #-- SGE host reservation
#$ -l mem_free=20G
#$ -pe smp 10          #-- Request 10 slots in the SMP parallel environment
            
date # start time

hostname

conda activate vdgs

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s "CNC(=N)N" -c "CNC(=N)N" -p /wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/ -b /wynton/group/degradolab/skt/docking/databases/probe_output/ -o /wynton/group/degradolab/skt/docking/databases/frag_vdg_lib/ 

date # end time

echo "DONE"
