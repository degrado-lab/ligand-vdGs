#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -R yes
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -l h_rt=24:00:00

# Description:
#   SGE array job: runs Probe on one PDB file per task.
#   The PDB database must be in RCSB mirror layout (subdirs by inner 2 chars of PDB code).
#
# Usage:
#   qsub \
#     -t 1-<number_of_pdbs> \        # one task per PDB; SGE_TASK_ID is 1-based
#     -o <path_to_log_directory> \   # SGE may require this dir to already exist
#     -e <path_to_log_directory> \
#     <path_to_this_script> \
#     <path_to_pdb_dir> \            # $1: PDB database directory
#     <path_to_probe_output_dir> \   # $2: directory to write Probe output files
#     <path_to_SMARTS-vdg_package> \ # $3: path to the SMARTS-vdg package
#     <path_to_probe_executable>     # $4: probe binary

date
hostname

# Get command line arguments
pdb_dir=$1
output_probe_dir=$2
SMARTS_vdg_package=$3
probe_path=$4

# Find all PDB files and store them in an array
pdb_files=($(find "$pdb_dir" -type f -name "*.pdb"))

# Define a specific PDB file for this task using the SGE_TASK_ID
pdb_file=${pdb_files[$((SGE_TASK_ID - 1))]} # SGE_TASK_ID is 1-based

# Run Probe
echo "Processing $pdb_file" 
python _run_probe.py \
  --input-pdb $pdb_file \
  --outdir $output_probe_dir \
  --probe-path $probe_path

# Check if the command was successful
if [ $? -ne 0 ]; then
    echo "Error processing $pdb_file" 
else
    echo "Successfully processed $pdb_file" 
fi

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
