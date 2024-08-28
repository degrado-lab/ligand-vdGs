#!/bin/bash
#SBATCH --job-name=prepwizard
#SBATCH --partition=gpueight
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --output=array_%A-%a.log
#SBATCH --array=1-15
pwd; hostname; date


if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  batchInd="0"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
  batchInd="1"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 3 ]; then
  batchInd="2"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 4 ]; then
  batchInd="3"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 5 ]; then
  batchInd="4"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 6 ]; then
  batchInd="5"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 7 ]; then
  batchInd="6"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
  batchInd="7"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 9 ]; then
  batchInd="8"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 10 ]; then
  batchInd="9"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 11 ]; then
  batchInd="10"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 12 ]; then
  batchInd="11"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 13 ]; then
  batchInd="12"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 14 ]; then
  batchInd="13"

elif [ "$SLURM_ARRAY_TASK_ID" -eq 15 ]; then
  batchInd="14"

else
  echo "SLURM_ARRAY_TASK_ID does not match any expected values."
fi

python -m smart_vdms.tools.protonate_pdbs --batch_index $batchInd

done
