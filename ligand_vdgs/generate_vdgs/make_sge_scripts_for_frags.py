'''
Create SGE submission scripts to run vdg generation pipeline on fragments defined in 
the database_frags_dict.pkl file.

NOTE: the `symmetry_classes` flag of vdg_generation_wrapper.py has to be manually added to 
the submission script (when applicable) after it's generated. 
TODO: automate symmetry detection
'''

import os
import re
import pickle as pkl

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
h_rt = '150:00:00 '
num_procs = '20'

# ------------------------------------------------------------------------------------------

replace = {'$LOG_DIR': log_dir, 
           '$PDB_DIR': pdb_parent_dir,
           '$PROBE_DIR': probe_dir,
           '$OUTPUT_DIR': out_dir_for_vdg_lib,
           '$MAX_NUM_CLUS': max_num_clus, 
           '$RUN_TIME': h_rt,
           '$NUM_PROCS': num_procs}

def main():
    print("NOTE TO USER: the `symmetry_classes` flag of vdg_generation_wrapper.py has to be "
          "manually added to the submission script (when applicable) after it's generated.")
    if not os.path.exists(out_dir_for_sge_scripts):
        os.makedirs(out_dir_for_sge_scripts)
    if os.listdir(out_dir_for_sge_scripts):
        raise FileExistsError(f"Output directory {out_dir_for_sge_scripts} already has "
                              "files. Terminating to prevent overwriting.")
    frags_dict = pkl.load(open(frags_dict_path, 'rb'))
    for elements, smiles in frags_dict.items():
        for smiles, lignames in smiles.items():
            # apply counts threshold
            counts = len(set(lignames))
            if counts < counts_threshold:
                continue
            # apply size threshold
            copy_smiles = smiles
            copy_smiles = re.sub(r'[(){}\[\]=#@+\-1234]', '', copy_smiles)  # Remove unwanted characters
            copy_smiles = re.sub(r'(Cl|Br)', 'X', copy_smiles)  
            elems = copy_smiles
            size = len(elems)
            #if size != frag_size:
            if size > max_size_frag:
                continue
            # create SGE submission script
            output_script(template, smiles)

def output_script(template, smiles):
    # create output script name
    script_name = os.path.join(out_dir_for_sge_scripts, smiles + '.sh')
    # replace $SMILES in template with `smiles`. dict $SMILES will be overwritten 
    # each time output_script is called, but it's ok.
    replace['$SMILES'] = f'"{smiles}"'
    replace['$CG'] = f'"{smiles}"' 
    replace['$JOB_NAME'] = f'"{smiles}"'.replace("#", "hash") # safe name
    # read template
    with open(template, 'r') as f:
        lines = f.readlines()
    # write new script
    with open(script_name, 'w') as f:
        start_copy = False
        for line in lines:
            # skip header. wait until #!/bin/bash
            if line.startswith('#!/bin/bash'):
                start_copy = True
            if start_copy:
                # replace placeholders based on the `replace` dict
                for key, value in replace.items():
                    line = line.replace(key, value) 
                # write out line
                f.write(line)

main()