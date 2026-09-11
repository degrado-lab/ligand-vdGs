'''
This script automates batch execution of the Schrodinger prepwizard program to optimize
hydrogens for all pdb structures in the parent databse. This script isn't meant for 
direct execution in most cases. Use s02_run_prepwizard.sh as the main entry point instead. 

Prepwizard flags: 
-disulfides (add?)
-glycosylation (add?)
-rehtreat
-noimpref

There shouldn't be waters in the BioliP2 database, so you shouldn't need flags for
treatment of waters, but double check if there are waters in any of the PDBs.

TODO: Describe how to input the batch numbers.
'''

# Script settings. TODO: convert to command line args

input_pdb_database_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_trimmed'
output_dir = '/home/sophia/DockDesign/databases/prepwizard_BioLiP2'
overwrite = False
prepwizard_bin_path = '/home/rkormos/schrodinger2024-1/utilities/prepwizard'
num_batches_total = 15 # TODO: SHOULD NOT BE HARDCODED
logdir = './prepwizard_logs'

import os
import sys
import subprocess
import time
import argparse
from ligand_vdgs.functions import parent_db
from ligand_vdgs.functions.utils import set_up_outdir
from ligand_vdgs.preprocessing._prep_filters import (
    clean_pdb_file, snapshot_restorable_ligands, restore_renamed_ligands)

def determine_subdirs_in_batch(all_subdirs, num_batches, batch_index):
    # Ensure num_batches is at least 1 to avoid division by zero
    if num_batches < 1:
        raise ValueError('num_batches must be at least 1')
    # Ensure batch_index is < num_batches
    if batch_index >= num_batches:
        raise ValueError(
            'The specified batch_index cannot be greater than the number of batches requested.')
    # Calculate the size of each batch
    batch_size = len(all_subdirs) // num_batches
    remainder = len(all_subdirs) % num_batches

    batches = []
    start = 0
    
    for i in range(num_batches):
        # Calculate the end index for the current batch
        end = start + batch_size + (1 if i < remainder else 0)
        batches.append(all_subdirs[start:end])
        start = end

    return batches[batch_index]

def main():
    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch_index', type=int, help='Index to define batch.')
    args = parser.parse_args()
    batch_index = args.batch_index

    
    # Ensure that the pdb database dir has subdirs formatted similarly to the RCSB
    # mirror format (see docs/database_generation.txt).
    if not parent_db.is_mirror(input_pdb_database_dir):
        sys.exit('Database structure must be similar to the RCSB mirror format; '
                 'see docs/database_generation_guide.md')

    # Create output dir and handle whether the user specified to overwrite or not
    set_up_outdir(output_dir, overwrite)

    # Determine which subdirs are part of this batch
    all_subdirs = sorted(os.listdir(input_pdb_database_dir))
    subdirs_in_batch = determine_subdirs_in_batch(all_subdirs, num_batches_total, batch_index)
    
    # Set up log file
    logpath = os.path.join(logdir, f'batchIndex{batch_index}.log')
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    if os.path.exists(logpath):
        os.remove(logpath)
    print('Writing out log to', logpath)
    
    # Use Schrodinger's prepwizard to optimize hydrogens in each PDB structure
    for subdir in  os.listdir(input_pdb_database_dir):
        if subdir not in subdirs_in_batch:
            continue
        subdir_path = os.path.join(input_pdb_database_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue
        # Iterate over PDB files
        for pdbname in os.listdir(subdir_path):
            output_path = os.path.join(
                output_dir, parent_db.shard(parent_db.stem_of(pdbname)), pdbname)
            # Skip if prepped pdb exists and overwrite == False
            if os.path.exists(output_path) and not overwrite:
                continue
            # Otherwise, run prepwizard on the pdb
            input_pdbpath = os.path.join(subdir_path, pdbname)
            run_prepwizard(pdbname, input_pdbpath, prepwizard_bin_path, logpath, output_path)
            prep_log = f'{parent_db.stem_of(pdbname)}.log'
            if os.path.exists(prep_log):
                os.remove(prep_log)
    with open(logpath, 'a') as file:
        file.write('Script successfully completed.') 

def run_prepwizard(pdbname, input_pdbpath, prepwizard_bin_path, logpath, output_path):
    # Create a tmp pdb path because prepwizard will only output if the output path 
    # is within the dir that this script is submitted. After the output pdb is 
    # generated, you may move the tmp pdb to its intended path.
    cwd = os.getcwd()
    tmpdir = os.path.join(cwd, 'tmp')
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    tmp_pdbpath = os.path.join(tmpdir, 'temp'+pdbname) 

    # Strip residues prepwizard mangles into free molecules relabelled onto a ligand
    # (see _prep_filters). s01 does this too, but s01 is optional, and this is the last
    # point before the damage.
    clean_pdbpath = os.path.join(tmpdir, 'clean'+pdbname)
    try:
        n_dropped = clean_pdb_file(input_pdbpath, clean_pdbpath)
    except Exception as e:
        n_dropped = None
        with open(logpath, 'a') as file:
            file.write(f'[WARNING] could not pre-filter {input_pdbpath}: {e}\n')
    if n_dropped:
        with open(logpath, 'a') as file:
            file.write(f'Dropped {n_dropped} hazard residue(s) from '
                       f'{pdbname} before prepwizard.\n')
        input_pdbpath = clean_pdbpath

    # Prepwizard also rewrites some HET ligands into ATOM records under a standard amino
    # acid resname, which hides them from find_cg_matches (HETATM-only). Record what went
    # in so the rename can be undone afterwards.
    try:
        ligand_snapshot = snapshot_restorable_ligands(input_pdbpath)
    except Exception as e:
        ligand_snapshot = {}
        with open(logpath, 'a') as file:
            file.write(f'[WARNING] could not snapshot ligands in {input_pdbpath}: {e}\n')

    cmd =  [prepwizard_bin_path, input_pdbpath, tmp_pdbpath, '-rehtreat', '-noimpref']
    cmd = ' '.join(cmd)
    with open(logpath, 'a') as file: 
        file.write(f'Full command: {cmd} \n')
    
    # Execute the command and redirect stdout and stderr to the log file
    with open(logpath, 'a') as log_file:
        try:
            # There's no good way to handle waiting until the prepwizard-generated child process is
            # completed, so make this script wait until the expected pdb file is created.
            subprocess.run(cmd, stdout=log_file, stderr=log_file, shell=True, 
                           check=True)
            
            if not wait_for_completion(tmp_pdbpath):
                raise TimeoutError('Prepwizard did not complete within timeout period.')

            try:
                n_restored = restore_renamed_ligands(tmp_pdbpath, ligand_snapshot)
            except Exception as e:
                n_restored = 0
                log_file.write(f'[WARNING] could not restore renamed ligands in '
                               f'{tmp_pdbpath}: {e}\n')
            if n_restored:
                log_file.write(f'Restored {n_restored} ligand residue(s) renamed by '
                               f'prepwizard in {pdbname}.\n')

            move_temp_pdb(tmp_pdbpath, output_path)

        except subprocess.CalledProcessError as e:
            log_file.write('-'*40+'\n')
            log_file.write(f'PREPWIZARD FAILED ON {input_pdbpath} \n')
            log_file.write(str(e)+'\n')
            log_file.write('-'*40+'\n')

        except TimeoutError as e:
            log_file.write('-'*40+'\n')
            log_file.write(f'PREPWIZARD FAILED ON {input_pdbpath} \n')
            log_file.write(str(e)+'\n')
            log_file.write('-'*40+'\n')
    
    # After the child process is completed, move the pdb file in the temp path to the intended path.

def move_temp_pdb(tmp_pdbpath, output_path):
    # Ensure that the parent dirs of the output path exist, then output the pdb.
    parent_dir = os.path.dirname(output_path)
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

    os.system(f'mv {tmp_pdbpath} {output_path}')

def wait_for_completion(expected_output, timeout=1800):
    start_time = time.time()
    while time.time() - start_time < timeout: 
        if os.path.exists(expected_output):
            time.sleep(60) # wait another minute to make sure all output is fully written.
            return True
        time.sleep(10) # wait 10 seconds
    return False

if __name__ == "__main__":
    main()
