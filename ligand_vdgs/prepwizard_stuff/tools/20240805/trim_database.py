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
import pickle as pkl
import prody as pr

from smart_vdms.functions.interactions import add_pdb_to_nr_db_dict
from smart_vdms.functions.redundancy import print_only_nr_chains_interacting_with_lig

origin_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_split'
target_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_trimmed'
lig_avg_bfactor_cutoff = 40 # default is 40
output_database_dict_name = '20240805_database_dict.pkl'
overwrite = True # If set to true NEED TO ADD code for checking whether a pdbfile already exists (TODO)



def main():

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

    for subdir in os.listdir(origin_dir):
        subdir_path = os.path.join(origin_dir, subdir)

        pdbfiles = os.listdir(subdir_path)

        #pkl_dumped = False # TODO: TAKE THIS OUT, FOR TESTING ONLY
        for pdbfile in pdbfiles:
            print(pdbfile)
            try:
                #if pdbfile_ix == 60:
                #    pkl_dumped = True
                #    break

                pdbpath = os.path.join(subdir_path, pdbfile)

                # Add the ligand(s) and interacting residues to the ligand dict
                # if nonredundant. 
                database_dict = add_pdb_to_nr_db_dict(
                    database_dict, pdbpath, lig_bfactor_cutoff=lig_avg_bfactor_cutoff,
                    unrefined=True) # unrefined is a quick and dirty way to check
                                    # for redundancy (computationally cheaper. Option
                                    # to refine later in the vdg creation process.)
            except:
                print('PDB FAILED: ', pdbfile)
                pass

            pdbfile_ix += 1
            if pdbfile_ix % 1000 == 0:
                print(pdbfile_ix)
                checkpoint_name = output_database_dict_name.rstrip('.pkl') + '_checkpoint.pkl'
                pkl.dump(database_dict, open(checkpoint_name, 'wb'))

        #if pkl_dumped: # TODO: get rid
        #    break

    pkl.dump(database_dict, open(output_database_dict_name, 'wb'))


    database_dict = pkl.load(open('database_dict.pkl', 'rb'))
    # Now that the ligand instances have been deduplicated, print out only the deduplicated chains
    # in the deduplicated pdbs.
    # Print out only the pdbs and chains that are in the deduplciated database_dict
    print_only_nr_chains_interacting_with_lig(database_dict, origin_dir, target_dir, overwrite)
    



    

        



    

























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
