#!/bin/bash
#SBATCH --job-name=prepwizard
#SBATCH --partition=gpueight
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --output=array_%A-%a.log

pwd; hostname; date

# USAGE INSTRUCTIONS: provide the number of tasks to submit concurrently, e.g.
# >> bash s02_run_prepwizard.sh <NUMBER_OF_TASKS>

# Ensure that the user has provided a number of tasks as the first arg
if [ -z "$1" ]; then
  echo "Usage: $0 <number_of_tasks>"
  echo "  <number_of_tasks>: The number of parallel tasks to submit."
  echo "Example: $0 15  # This will submit 15 tasks."
  exit 1
fi

# Set the number of tasks
NUM_TASKS=$1

# Adjust the SLURM array size (dynamically) based on the number of tasks
#SBATCH --array=1-$NUM_TASKS

# Determine the batch index for each task
batchInd=$((SLURM_ARRAY_TASK_ID - 1))

# Print the batch index 
echo "Running task with batch index: $batchInd"

# Call protonation script with the batch index
python -m ligand_vdgs.preprocessing._protonate_pdbs --batch_index $batchInd
