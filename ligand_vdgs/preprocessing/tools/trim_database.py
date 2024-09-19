'''
Stand-alone script.
Reduces size of the 50G consolidated_BioLiP2_split/ database by 
- filtering by b-factor
- crudely determining if a binding site is equivalent by collecting all the vdm
  resnums and resnames of a lig, and determining if that same set of vdm resnums,
  resnames, and lig resnums/resnames are repeated in that same pdb or a diff pdb.
  (there's some logic about taking a guess at whether the pdbs are within the same
  series or not - double check that and describe here.)

Usage: (TODO: verify)
cd $YOUR_SMART-VDMS_DIR
pip install -e . # for debugging and developing
pip install .    # for users
python -m smart_vdms.scripts.trim_database

NOTE the reason this script takes forever is because there are a lot of 
"ligands" that are like ALAA, ALAB, ALA1, ALA2, etc. and are not actually
ligands; they're amino acids, but the pdb files are so large that there are not
enough chain columns that they run into each other.
'''

import os
import traceback
import numpy as np
import pickle as pkl
import json
import prody as pr

from smart_vdms.functions.interactions import add_pdb_to_nr_db_dict
from smart_vdms.functions.utils import set_up_outdir

origin_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_split'
target_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_trimmed'
lig_avg_bfactor_cutoff = 40 # default is 40
output_database_dict_name = '20240809_database_dict.json'
overwrite_pdbs = False # If set to true NEED TO ADD code for checking whether a pdbfile already exists (TODO)
skip_to_output_pdbs = True
radius = 10

pdbs_already_done = 'log_processed_pdbs.txt'
previous_checkpoint_dict = 'checkpoint.pkl'


def main():
                
    checkpoint_name = output_database_dict_name + '.checkpoint'

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

    #if restart: # the name restart could be confusng bc now it should be reload
    #    with open(pdbs_already_done) as inF:
    #        pdbs_done = [i.strip() for i in inF.readlines()]
    #    print("Starting json load from previous run", flush=True) # TODO: problem; there might be >1 prev runs
    #    database_dict = json.load(open(previous_checkpoint_dict, 'r'))
    #    print('Num of PDBs processed in a previous run: ', 
    #        len([i for i in pdbs_done if i.endswith('.pdb')]))

    if not skip_to_output_pdbs:
        for subdir in os.listdir(origin_dir):
            subdir_path = os.path.join(origin_dir, subdir)

            pdbfiles = os.listdir(subdir_path)

            for pdbfile in pdbfiles:
                pdbfile_ix += 1
                #if restart:
                #    if pdbfile in pdbs_done:
                #        continue
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
    
    AAs = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 
            'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    
    database_dict = json.load(open(output_database_dict_name, 'r'))
    
    set_up_outdir(target_dir, overwrite_pdbs)
    # First, get all the lig and vdm residues in each pdb, so you only have to process each 
    # pdb once.
    bindingsite_dict = {} # key=pdb, values=list of ligs, where lig=[segment, chain, resnum]
    for ligresn, lig_interactions in database_dict.items():
        # SKIP OVER LIGAND IF IT'S ACTUALLY AN AA
        if len(ligresn) == 4 and ligresn[:3] in AAs:
            continue
        for lig_res in lig_interactions.keys():
            lig_res_ = lig_res.split(' ')
            pdbname = lig_res_[0]
            seg_ch_res = lig_res_[1:4]
            if pdbname not in bindingsite_dict.keys():
                bindingsite_dict[pdbname] = [seg_ch_res]
            else:
                bindingsite_dict[pdbname].append(seg_ch_res)

    # Then, select *radius* around the lig residues, and write out the pdb.
    for pdb_name, list_ligs in bindingsite_dict.items():
        # First, determine if this pdb was already written out in a previous run (if this script is
        # resuming from a previous incomplete run).
        original_pdb_subdir = pdb_name[1:3]
        output_subdir = os.path.join(target_dir, original_pdb_subdir)
        output_path = os.path.join(output_subdir, pdb_name)
        if not overwrite_pdbs and os.path.exists(output_path):
            continue

        # Load pdb
        original_pdbpath = os.path.join(origin_dir, original_pdb_subdir, pdb_name)
        parsed = pr.parsePDB(original_pdbpath)
        # Select all ligand residues
        all_ligs_sel = ''
        for lig_residue in list_ligs:
            _seg, _ch, _resnum = lig_residue 
            if _seg == '':
                lig_sel = f'(chain {_ch} and resnum {_resnum})'
            else:
                lig_sel = f'(segment {_seg} chain {_ch} and resnum {_resnum})'
            if all_ligs_sel == '':
                all_ligs_sel = lig_sel
            else:
                all_ligs_sel += f' or {lig_sel}'
        try:
            # Select sphere around all ligand residues (entire residues)
            sele_around_lig_str = f'within {radius} of ({all_ligs_sel})'
            resinds_around_lig = set(parsed.select(sele_around_lig_str).getResindices())
            resindices_around_lig = ' '.join([str(i) for i in resinds_around_lig])
            around_lig = parsed.select(f'resindex {resindices_around_lig}')

            # Write out PDB
            if not os.path.isdir(output_subdir):
                os.makedirs(output_subdir)
            pr.writePDB(output_path, around_lig)
            print(output_path)
        except Exception as e:
            print('--------------------------------------')
            print('FAILED TO OUTPUT', output_path)
            print(e)
            traceback.print_exc()
            print('--------------------------------------')
    # Check number of PDBs that were actually output
    num_output_pdbs = 0
    for pdb_output_subdir in os.listdir(target_dir):
        for pdb_file in os.listdir(os.path.join(target_dir, pdb_output_subdir)):
            if pdb_file.endswith('.pdb'):
                num_output_pdbs += 1
    print('Number of deduplicated PDBs processed:', len(bindingsite_dict.keys()))
    print('Number of PDBs successfully output:', num_output_pdbs)





    print('Script successfully completed.')

    


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
