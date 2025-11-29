'''
Create SGE submission scripts to run vdg generation pipeline on fragments defined in 
the database_frags_dict.pkl file.

NOTE: the `symmetry_classes` flag of vdg_generation_wrapper.py has to be manually added to 
the submission script (when applicable) after it's generated. 
TODO: automate symmetry detection
'''

import os
import sys
import re
import pickle as pkl
from rdkit import Chem
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import utils

frags_dict_path = 'resources/database_frags_dict.pkl'
template = 'resources/frag_sge_template.sh' 
counts_threshold = 80 # min number of ligs containing this frag in order to run vdg 
                        # vdg generation on
max_size_frag = 5 # 5 atoms. 
out_dir_for_sge_scripts = 'ligand_vdgs/generate_vdgs/frag_submit_scripts/'
out_dir_for_vdg_lib = '/wynton/home/degradolab/skt/docking/frag_lib'
log_dir = '/wynton/home/degradolab/skt/docking/frag_sge_logs'
pdb_parent_dir = '/wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/'
probe_dir = '/wynton/group/degradolab/skt/docking/databases/probe_output/'
max_num_clus = '5000' # max number of vdgs to cluster for each vdg subset
h_rt = '300:00:00 '
num_procs = '20'
# Provide a file with smiles to run (one per line), or set to None
smiles_to_run_file = '/wynton/home/degradolab/skt/logs/needed_frags.txt'

# ------------------------------------------------------------------------------------------

def smiles_in_file(smiles_to_run_file):
    smiles_to_run = set()
    with open(smiles_to_run_file, 'r') as f:
        for line in f:
            smiles_to_run.add(line.strip())
    return smiles_to_run

replace = {'$LOG_DIR': log_dir, 
           '$PDB_DIR': pdb_parent_dir,
           '$PROBE_DIR': probe_dir,
           '$OUTPUT_DIR': out_dir_for_vdg_lib,
           '$MAX_NUM_CLUS': max_num_clus, 
           '$RUN_TIME': h_rt,
           '$NUM_PROCS': num_procs}

def main():
    if not os.path.exists(out_dir_for_sge_scripts):
        os.makedirs(out_dir_for_sge_scripts)
    if os.listdir(out_dir_for_sge_scripts):
        raise FileExistsError(f"Output directory {out_dir_for_sge_scripts} already has "
                              "files. Terminating to prevent overwriting.")
    # Determine whether to run smiles for smiles_to_run_file or from frags_dict.
    if smiles_to_run_file:
        if os.path.exists(smiles_to_run_file):
            smiles_to_run = smiles_in_file(smiles_to_run_file)

    else:
        smiles_to_run = set()
        frags_dict = pkl.load(open(frags_dict_path, 'rb'))
        for elements, smiles in frags_dict.items():
            for smiles, lignames in smiles.items():
                # apply counts threshold
                counts = len(set(lignames))
                if counts < counts_threshold:
                    continue
                # apply size threshold
                copy_smiles = smiles
                copy_smiles = re.sub(r'[(){}\[\]=#@+\-1234]', '', copy_smiles)
                copy_smiles = re.sub(r'(Cl|Br)', 'X', copy_smiles)
                elems = copy_smiles
                size = len(elems)
                #if size != frag_size:
                if size > max_size_frag:
                    continue
                smiles_to_run.add(smiles)

    # iterate over smiles and create scripts
    for smiles in smiles_to_run:
        # determine symmetry
        symm_classes = utils.define_symmetry(smiles)
        # output the job script for this frag
        output_script(template, smiles, symm_classes)

    print(f'RUNNING {len(smiles_to_run)} FRAGMENTS')
    return

def output_script(template, smiles, symm_classes):
    # per-fragment copy of the base mapping
    local_replace = dict(replace)

    script_name = os.path.join(out_dir_for_sge_scripts, smiles + '.sh')
    local_replace['$SMILES'] = f'"{smiles}"'
    local_replace['$CG'] = f'"{smiles}"'
    local_replace['$JOB_NAME'] = f'"{smiles}"'.replace("#", "hash") # safe name

    # always define how to handle --num-procs
    if symm_classes:
        local_replace['--num-procs'] = f'--symmetry-classes {symm_classes} --num-procs'
    else:
        local_replace['--num-procs'] = '--num-procs'

    with open(template, 'r') as f:
        lines = f.readlines()

    with open(script_name, 'w') as f:
        start_copy = False
        for line in lines:
            # skip header. wait until #!/bin/bash
            if line.startswith('#!/bin/bash'):
                start_copy = True
            if start_copy:
                # replace placeholders based on the `replace` dict
                for key, value in local_replace.items():
                    line = line.replace(key, value)
                f.write(line)

main()
