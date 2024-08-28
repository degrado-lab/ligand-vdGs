'''
Reduces size of the 50G consolidated_BioLiP2_split/ database by 
- ...

Usage: 
cd $YOUR_SMART-VDMS_DIR
pip install -e .
??????? pip install without -e
python -m smart_vdms.scripts.trim_database


'''




import os
import traceback
import numpy as np
import pickle as pkl
import json
import prody as pr

from smart_vdms.functions.interactions import add_pdb_to_nr_db_dict
from smart_vdms.functions.redundancy import print_lig_shell

origin_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_split'
target_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_trimmed'
lig_avg_bfactor_cutoff = 40 # default is 40
output_database_dict_name = '20240809_database_dict.json'
overwrite_pdbs = True # If set to true NEED TO ADD code for checking whether a pdbfile already exists (TODO)


restart = True
pdbs_already_done = 'log_processed_pdbs.txt'
previous_checkpoint_dict = 'checkpoint.pkl'


def main():
                
    print('Starting...', flush=True)
    checkpoint_name = output_database_dict_name.rstrip('.pkl') 

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

    if restart:
        with open(pdbs_already_done) as inF:
            pdbs_done = [i.strip() for i in inF.readlines()]
        print("Starting pickle load from previous run", flush=True)
        database_dict = pkl.load(open(previous_checkpoint_dict, 'rb'))
        #database_dict = json.load(open(previous_checkpoint_dict, 'r'))
        print("Done with pickle load")
        print('Num of PDBs processed in a previous run: ', 
            len([i for i in pdbs_done if i.endswith('.pdb')]))

        print('Converting dict')
        # Converting keys to strings
        convert_nested_dict_keys_and_values(database_dict)

    for subdir in os.listdir(origin_dir):
        subdir_path = os.path.join(origin_dir, subdir)

        pdbfiles = os.listdir(subdir_path)

        for pdbfile in pdbfiles:
            pdbfile_ix += 1
            if restart:
                if pdbfile in pdbs_done:
                    continue
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

            if pdbfile_ix % 1000 == 0:
                print(pdbfile_ix)
                print(f"Dumping on pdbfile_ix: {pdbfile_ix}", flush=True)
                json.dump(database_dict, open(checkpoint_name, "w"))

        #if pkl_dumped: # TODO: get rid
        #    break

    
    print("Starting final JSON dump", flush=True)
    json.dump(database_dict, open(output_database_dict_name, "w"))
    print("Done with final JSON dump", flush=True)

    # Now that the ligand instances have been deduplicated, print out only the binding sites
    # (10A around the lig)
    print_lig_shell(database_dict, origin_dir, target_dir, overwrite_pdbs)
    #print_only_nr_chains_interacting_with_lig(database_dict, origin_dir, target_dir, overwrite_pdbs)
    

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
