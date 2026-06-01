# Template for submitting jobs to SGE. This template is used by 
# ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py 

#!/bin/bash

#$ -S /bin/bash
#$ -o $LOG_DIR 
#$ -N $JOB_NAME 
#$ -cwd
#$ -j y                 #-- join STDERR and STDOUT
#$ -l h_rt=$RUN_TIME    #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes               #-- SGE host reservation
##$ -l mem_free=20G
##$ -l scratch=50G
#$ -pe smp $NUM_PROCS           #-- Request # of slots in SMP parallel environment
date # start time

hostname

conda activate py3.10

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s $SMILES -c $CG -p $PDB_DIR -b $PROBE_DIR -o $OUTPUT_DIR -m $MAX_NUM_CLUS --num-procs $NUM_PROCS 

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

date # end time

echo "DONE"
