# Template for submitting jobs to SGE. This template is used by 
# ligand_vdgs/generate_vdgs/make_sge_scripts.py 

#!/bin/bash

#$ -S /bin/bash
#$ -o $LOG_DIR 
#$ -N $JOB_NAME 
#$ -cwd
#$ -j y                #-- tells system STDERR and STDOUT should be joined
#$ -l h_rt=$RUN_TIME   #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes              #-- SGE host reservation
#$ -l mem_free=5G
#$ -pe smp 10          #-- Request 10 slots in the SMP parallel environment
            
date # start time

hostname

conda activate py3.10

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s $SMILES -c $CG -p $PDB_DIR -b $PROBE_DIR -o $OUTPUT_DIR -m $MAX_NUM_CLUS 

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

date # end time

echo "DONE"
