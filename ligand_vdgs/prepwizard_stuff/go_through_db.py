'''
Reduces size of the 50G consolidated_BioLiP2_split/ database by: 
- ...


'''

import os
import prody as pr

origin_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_split'
target_dir = ''

# Define the size limit in bytes
size_limit = 2 * 1024 * 1024


def main():

    # Iterate through pdb files
    num_large_files = 0
    total_num_pdbfiles = 0

    # Initialize nested dict to keep track of pdbs, segs, chains, resnums, resnames, etc.
    # Nested pdbs_dict format: {pdb: 
    #                            {lig res: [[nested list of interacting vdms 
    #                            (each element has seg, chain, resnum, resname)]] }
    # See get_lig_interacting_chains() for more details. 
    pdbs_dict = {} 



    # Get the file size 
    for subdir in os.listdir(origin_dir):
        subdir_path = os.path.join(origin_dir, subdir)

        pdbfiles = os.listdir(subdir_path)

        for pdbfile in pdbfiles:
            total_num_pdbfiles += 1
            pdbpath = os.path.join(subdir_path, pdbfile)
            pdbs_dict = get_lig_interacting_chains(pdbs_dict, pdbpath)
            if total_num_pdbfiles % 1000 == 0:
                print(total_num_pdbfiles)
            print(total_num_pdbfiles)
 
             
        

def get_lig_interacting_chains(pdbs_dict, pdbpath):
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
    
    Nested pdbs_dict format: {pdb: 
                                   {lig residue: [[nested list of interacting vdms 
                                   (each element has seg, chain, resnum, resname)]] }
    '''

    pdbfile = pdbpath.split('/')[-1]
    if not pdbfile.endswith('pdb'):
        print('NOT A PDB:', pdbpath)
        return pdbs_dict
    pdbs_dict[pdbpath] = {}
    
    # Identify ligand(s)
    complex = pr.parsePDB(pdbpath)
    ligand = complex.select('not (water or ion or protein)') # prody "hetero" selection not suitable
    if ligand is None: # could be None if the "ligand" is actually a noncanonical AA
        return pdbs_dict
    ligand_resinds = set(ligand.getResindices())
    for lig_resind in ligand_resinds: # but don't store the lig resind b/c it will get re-indexes
                                      # in the trimmed pdb
        interacting_residues = []
        # Determine which segments, chains, and residues each lig interacts with
        neighbs = complex.select(f'within 5 of resindex {lig_resind}')
        for neighb in neighbs:
            seg, chain, resnum, resname = neighb.getSegname(), neighb.getChid(), neighb.getResnum(), \
                neighb.getResname()
            vdm_id = [seg, chain, resnum, resname]
            if vdm_id not in interacting_residues:
                interacting_residues.append(vdm_id)

    return pdbs_dict


main()