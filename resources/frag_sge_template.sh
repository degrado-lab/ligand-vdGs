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
#$ -l mem_free=$MEM_FREE   #-- PER SLOT under -pe smp, so total = this x $NUM_PROCS
#$ -l scratch=$SCRATCH     #-- node-local scratch; the wrapper writes shards to $TMPDIR
#$ -pe smp $NUM_PROCS           #-- Request # of slots in SMP parallel environment
date # start time

hostname

source ~/miniconda3/etc/profile.d/conda.sh
conda activate lig_vdgs

python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py -s $SMILES -c $CG -p $PDB_DIR -b $PROBE_DIR -o $OUTPUT_DIR -m $MAX_NUM_CLUS --subset-sizes $SUBSET_SIZES --num-procs $NUM_PROCS
# Captured before anything else runs, and used as this script's exit status
# below. Without it the trailing commands make the job exit 0 whatever the
# wrapper did, so a crashed fragment is indistinguishable from a finished one in
# `qacct` and only the 'Job completed.' log line reveals the difference.
STATUS=$?

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

date # end time

if [ "$STATUS" -eq 0 ]; then
    echo "DONE"
else
    echo "FAILED (exit $STATUS)" >&2
fi
exit "$STATUS"
