# Template for submitting jobs to SGE. This template is used by 
# ligand_vdgs/score_poses/make_sge_scripts_for_hit_finder.py 

#!/bin/bash

#$ -S /bin/bash
#$ -o $LOG_DIR 
#$ -N $JOB_NAME 
#$ -cwd
#$ -j y                 #-- join STDERR and STDOUT
#$ -l h_rt=5:00:00      #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes               #-- SGE host reservation
#$ -pe smp $NUM_PROCS           #-- Request # of slots in SMP parallel environment
date # start time

hostname

source ~/.bashrc
conda activate lig_vdgs

python ligand_vdgs/score_poses/vdg_hit_finder.py --smiles $SMILES --query-dir $QUERY_DIR --vdg-lib-dir $VDG_LIB_DIR --nprocs $NUM_PROCS --outdir $OUTDIR $OPTIONAL_ARGS

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

date # end time

echo "DONE"
