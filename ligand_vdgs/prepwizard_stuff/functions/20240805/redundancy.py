import os
import prody as pr
from smart_vdms.functions.utils import set_up_outdir

def check_pdbnames(pdbfile1, pdbfile2):
    # Returns True if the pdb names share at least 3 out of 4 characters 
    # (i.e. likely deposited as part of the same project)
    num_matches = 0
    pdbacc1 = pdbfile1[:4]
    pdbacc2 = pdbfile2[:4]
    for char1, char2 in zip(pdbacc1, pdbacc2):
        if char1 == char2:
            num_matches += 1
    if num_matches >= 3:
        return True
    else:
        return False


def check_networks(candidate_interacting_residues, already_deduplicated_network, tol=0):
    # Quick and dirty method of checking redundancy by checking if the sets of interacting
    # residues share the same resnums and resnames. Returns True if redundant.
    # General rule is tolerate 2 missing residues for checking intra-pdb, and tolerate 0 for
    # checking inter-pdb (b/c it may be missing one interaction with a whole other ligand, or it
    # maybe a mutant, like a Lys mutant may show up in the network but not Ala b/c not close enough.)

    a = [(i[2], i[3]) for i in candidate_interacting_residues]
    b = [(i[2], i[3]) for i in already_deduplicated_network]

    num_matches = len(set(a) & set(b)) # intersection of the 2 sets
    num_vdms_in_vdg_between_both_sets = len(set(a) | set(b))

    # Determine to be redundant if the vdms in set a are in set b, tolerating up to 
    # 2 "missing" residues.
    if num_matches >= num_vdms_in_vdg_between_both_sets - tol:
        return True
    return False


def get_interacting_chains_from_each_pdb(database_dict):
    # Return dict where key = pdbfile, value = set of chains 
    
    chains_to_keep = {} # i.e. chains that were found to be interacting

    for lig_resname, interactions_dict in database_dict.items():
        for lig_instance, interactions in interactions_dict.items():
            pdbfile, lig_res = lig_instance
            if pdbfile not in chains_to_keep.keys():
                chains_to_keep[pdbfile] = set()
            # Add ligand's segment and chain
            lig_res_seg_ch = (lig_res[0], lig_res[1])
            chains_to_keep[pdbfile].add(lig_res_seg_ch)
            # Add each "vdm"'s segment and chain 
            for vdm in interactions:
                vdm_seg_ch = (vdm[0], vdm[1])
                chains_to_keep[pdbfile].add(vdm_seg_ch)

    return chains_to_keep

def print_only_nr_chains_interacting_with_lig(database_dict, origin_dir, target_dir, overwrite):
    set_up_outdir(target_dir, overwrite)
    chains_to_keep = get_interacting_chains_from_each_pdb(database_dict)
    num_pdbs_processed = 0
    for pdbfile, chains in chains_to_keep.items():
        original_pdb_subdir = pdbfile[1:3]
        original_pdbpath = os.path.join(origin_dir, original_pdb_subdir, pdbfile)
        parsed = pr.parsePDB(original_pdbpath)
        selection_str = ''
        chains = sorted(chains)
        for seg_ch in chains:
            _seg, _ch = seg_ch
            if _seg == '':
                chain_sel = f'(chain {_ch})'
            else:
                chain_sel = f'(segment {_seg} and chain {_ch})'
            if selection_str == '':
                selection_str = chain_sel
            else:
                selection_str += f' or {chain_sel}'

        if _seg != '':
            print('sel str: ', selection_str)

        output_subdir = os.path.join(target_dir, original_pdb_subdir)
        if not os.path.isdir(output_subdir):
            os.makedirs(output_subdir)
        output_path = os.path.join(output_subdir, pdbfile)
        prody_sel_chains = parsed.select(selection_str)
        pr.writePDB(output_path, prody_sel_chains)
        num_pdbs_processed += 1
    # Check number of PDBs that were actually output
    num_output_pdbs = 0
    for pdb_output_subdir in os.listdir(target_dir):
        for pdb_file in os.listdir(os.path.join(target_dir, pdb_output_subdir)):
            if pdb_file.endswith('.pdb'):
                num_output_pdbs += 1
    print('Number of deduplicated PDBs processed:', num_pdbs_processed)
    print('Number of PDBs successfully output:', num_output_pdbs)