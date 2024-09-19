'''
Reduces size of the 50G consolidated_BioLiP2_split/ database by 
- ...placeholder...

Usage: (TODO: verify) 
cd $YOUR_SMART-VDMS_DIR
pip install -e . # for debugging and developing
pip install .    # for users
python -m smart_vdms.tools.redo_untrimmed_pdbs
'''

import os
import traceback
import numpy as np
import pickle as pkl
import json
import prody as pr

from smart_vdms.functions.interactions import add_pdb_to_nr_db_dict

# TODO: convert to command-line args
origin_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_split'
target_dir = 'redo_pdbs'
lig_avg_bfactor_cutoff = 40 # default is 40
output_database_dict_name = '20240810_database_dict_remainingpdbs.json'
overwrite_pdbs = False # TODO: If set to true NEED TO ADD code for checking whether a pdbfile already exists

def main():
                
    print('Starting...', flush=True)
    checkpoint_name = 'checkpoint'+output_database_dict_name

    # Iterate through pdb files
    pdbfile_ix = 0

    # Initialize nested dict to keep track of pdbs, segs, chains, resnums, resnames, etc.
    # Nested database_ligands_dict format: lig resname:
    #                                           pdb: 
    #                                                lig res: 
    #                                                      list of interacting residues 
    #                                                      (segment, chain, rensum, resname)
    # See get_lig_interacting_chains() for more details. 
    database_dict = {} 

    for idx, subdir in enumerate(os.listdir(origin_dir)):
        if idx < 746:
            continue
        subdir_path = os.path.join(origin_dir, subdir)

        pdbfiles = os.listdir(subdir_path)

        for pdbfile in pdbfiles:
            pdbfile_ix += 1
            print(pdbfile, flush=True)
            try:

                pdbpath = os.path.join(subdir_path, pdbfile)

                # Add the ligand(s) and interacting residues to the ligand dict
                # if nonredundant. 
                database_dict = add_pdb_to_nr_db_dict(
                    database_dict, pdbpath, lig_bfactor_cutoff=lig_avg_bfactor_cutoff,
                    unrefined=True) # unrefined is a quick and dirty way to check
                                    # for redundancy (computationally cheaper. Option
                                    # to refine later in the vdg creation process.)
            except Exception as e:
                print('PDB FAILED: ', pdbfile)
                print(e)
                traceback.print_exc()

            if pdbfile_ix % 100 == 0:
                print(pdbfile_ix)
                print(f"Dumping on pdbfile_ix: {pdbfile_ix}", flush=True)
                json.dump(database_dict, open(checkpoint_name, "w"))

    print('Ending...')

def convert_nested_dict_keys_and_values(d):
    """Recursively converts the keys (tuples to underscore-separated strings)
    and values (tuples to hyphen-separated strings) of a nested dictionary in place."""
    
    for outer_key in list(d.keys()):  # Iterate over the outer dictionary's keys
        nested_dict = d[outer_key]  # Get the nested dictionary
        
        # Create a new dictionary to hold the updated keys and values
        updated_nested_dict = {}
        
        for inner_key, inner_value in nested_dict.items():
            # Convert the inner key tuple to an underscore-separated string
            new_inner_key = f"{inner_key[0]} {' '.join([str(i) for i in inner_key[1]])}"
            # Convert the inner value tuple to a hyphen-separated string
            new_inner_value = []
            for n in inner_value:
                n_sublist = []
                for sub_n in n:
                    if type(sub_n)==np.int64:
                        sub_n = int(sub_n) # json does not recognize numpy int
                    else:
                        sub_n = str(sub_n) # json does not recognize numpy str either
                    
                    
                    n_sublist.append(sub_n)
                new_inner_value.append(n_sublist)

            # Add the converted key-value pair to the updated dictionary
            updated_nested_dict[new_inner_key] = new_inner_value
        
        # Update the outer dictionary with the updated nested dictionary
        d[outer_key] = updated_nested_dict

'''
    Select each ligand and determine what chains it interacts with. 
    
    Downstream use: 
        -- If it only interacts with 1 chain, then take that monomer and make it
           a separate pdb to isolate it from irrelevant chains to reduce the database size.
        -- If it interacts with >1 chain, then need to keep those interacting chains. 
    
    Store the ligs (seg/chain/resnum) and interacting residues in a dict to further
    reduce the database size by making a guess at whether monomers within a pdb are redundant
    by looking up the lig and vdm resnums and seeing if they're the same across the different 
    monomers (intra-pdb redundancy). You can additionally make a guess about whether 2 pdbs are
    redundant by seeing if their acc. codes are similar (i.e. in the same series), and their ligs 
    and vdms are on the same chains/resnums (inter-pdb redundancy).
'''

if __name__ == "__main__":
    main()
