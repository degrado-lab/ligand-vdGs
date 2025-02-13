'''
Given a query structure (a crystal structure of a protein-ligand complex) and its 
binding site residues (vdms), determine whether any vdGs can correctly place the chemical 
groups within the ligand of the query structure.

Usage:
    >> python validate_known_poses.py $YOUR_YAML_FILE
'''

import sys
import os
import yaml
import numpy as np
import prody as pr
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import match_vdgs as match
import utils

script, yml_file = sys.argv

def main():
    # Parse inputs from yml file
    (vdgs_dir, cg_name, query_path, bindingsite_residues, query_lig_res, 
        query_cg_atoms, symm_query_cg_atoms, rmsd_threshold, 
        out_dir) = parse_yaml_input(yml_file)
    
    # Set up output dir
    pdbname = query_path.split('/')[-1].removesuffix('.pdb')
    output_dir = os.path.join(out_dir, yml_file.removesuffix('.yml').split('/')[-1])
    utils.handle_existing_files(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Get the coords of the query structure's ligand chemical group. Permute any
    # symmetrically equivalent atoms.
    query_struct = pr.parsePDB(query_path)
    q_lig_seg, q_lig_ch, q_lig_resnum = query_lig_res
    q_cg_coords = match.get_coords(query_struct, query_cg_atoms, q_lig_seg, q_lig_ch, 
                                   q_lig_resnum, query_path)
    q_cg_name_perms, q_cg_coord_perms = match.permute_cg_coords(q_cg_coords, query_cg_atoms, 
                                                                symm_query_cg_atoms)
     
    # Iterate over the query structure's binding site residues and superpose vdGs 
    # (incl. the CG) onto it, keeping track of the vdGs that have low rmsd to the 
    # known bb+CG.
    query_bb_coords = {}
    for bsr in bindingsite_residues:
        bsr_seg, bsr_chain, bsr_resnum = bsr
        bbcoords = match.get_coords(query_struct, ['N', 'CA', 'C'], bsr_seg, 
            bsr_chain, bsr_resnum, query_path)
        query_bb_coords[tuple(bsr)] = bbcoords
    
    # Get subsets of vdg residues in query struct.
    query_res_sets = match.subsets_of_query_residues(bindingsite_residues)

    # Iterate over subsets of binding site residues and determine
    # whether any of the vdGs in the database match the "ground truth" binding 
    # site residues and CG positions of the query structure.
    for q_res_set in query_res_sets:

        # Determine the AA resnames of the subset of query vdms.
        query_vdm_resnames = []
        for q_res in q_res_set:
            bsr_seg, bsr_chain, bsr_resnum = q_res
            q_res_obj = match.get_residue_obj(query_struct, bsr_seg, bsr_chain, 
                                              bsr_resnum, query_path)
            query_vdm_resnames.append(q_res_obj.getResnames()[0])
        
        # Retrieve database vdGs
        _vdgs_dir, out_subdir = match.get_in_dir_and_out_dir(q_res_set, vdgs_dir)
        database_vdg_paths, database_vdgs = match.get_database_vdgs_for_spec_AAs(
            _vdgs_dir, query_vdm_resnames)


        # Calc rmsd on subsets of bb residues (incl. CG) on CG atom permutations. 
        q_bbcoords_list = [query_bb_coords[tuple(r)] for r in q_res_set]
        q_flattened_bb_list = [item for sub in q_bbcoords_list for item in sub]
        # Iterate over the database vdGs.
        for database_vdg_index, database_vdg in enumerate(database_vdgs):
            db_vdg_path = database_vdg_paths[database_vdg_index]
            db_cg_coords = match.get_database_cg_coords(database_vdg)
            # Iterate over the query CG permutations and try to match them to the 
            # database vdGs.
            for q_cg_coord_perm in q_cg_coord_perms:
                q_bb_and_cg = np.array(q_flattened_bb_list + q_cg_coord_perm)
                if len(db_cg_coords) != len(q_cg_coord_perm) : 
                    raise ValueError(
                        'The query CG must have the same num. of atoms as the database CG.')

                # Get vdm bb coordinates if they're the same AA identity as the binding site 
                # residues being queried. Make sure all valid AA permutations are sampled.
                all_AA_perms_db_vdm_bb_coords_and_resinds = match.get_database_vdm_bb_coords(
                    database_vdg, query_vdm_resnames, db_vdg_path)
                if all_AA_perms_db_vdm_bb_coords_and_resinds is None:
                    continue
                
                # Iterate over the vdM permutations and superpose/calc RMSD against the
                # known structure (incl. CG).
                for perm_resinds, AA_perm_db_vdm_bb_coords in \
                                    all_AA_perms_db_vdm_bb_coords_and_resinds.items():
                    db_vdg_bb_and_cg = np.array(AA_perm_db_vdm_bb_coords + db_cg_coords)
                    
                    moved_db_vdg_bb_cg_coords, transf = pr.superpose(db_vdg_bb_and_cg, 
                                                                     q_bb_and_cg)
                    rmsd = round(pr.calcRMSD(moved_db_vdg_bb_cg_coords, q_bb_and_cg), 2)
                    if rmsd < rmsd_threshold:
                        moved_vdg = pr.applyTransformation(transf, database_vdg.copy())

                        # Write out just the vdms in this AA permutation
                        AA_perm_resinds_str = ' '.join([str(i) for i in perm_resinds])
                        AA_perm_resinds_sele = f'resindex {AA_perm_resinds_str}'
                        moved_vdg_AA_perm_resinds_obj = moved_vdg.copy().select(
                            f'({AA_perm_resinds_sele}) or (occupancy > 2.8)') # occ 3 = CG
                        vdg_pdbname = db_vdg_path.rstrip('/').split('/')[-1].removesuffix('.pdb')
                        output_vdg_name = f'{vdg_pdbname}'
                        for residue in q_res_set:
                            bsr_flat_list = [str(it) for it in residue]
                            bsr_str = '_'.join(bsr_flat_list)
                            output_vdg_name += f'_{bsr_str}'
                        subdir_path = os.path.join(output_dir, out_subdir)
                        output_vdg_path = os.path.join(output_dir, out_subdir, output_vdg_name)
                        if not os.path.exists(subdir_path):
                            os.makedirs(subdir_path)
                        pr.writePDB(output_vdg_path, 
                                    moved_vdg_AA_perm_resinds_obj.copy().select(
                                        'occupancy > 1.5'))

def parse_yaml_input(yml_file):
    with open(yml_file, 'r') as inF:
        user_data = yaml.load(inF, Loader=yaml.FullLoader)
    vdgs_dir = user_data['vdgs_dir']
    cg_name = user_data['cg_name']
    query_path = user_data['query_path']
    bindingsite_residues = user_data['bindingsite_residues']
    query_lig_res = user_data['query_lig']
    query_cg_atoms = user_data['query_cg_atoms']
    try:
        symm_query_cg_atoms = user_data['symmetric_query_cg_atoms']
    except:
        symm_query_cg_atoms = [[]]
    if type(symm_query_cg_atoms[0]) is not list:
        symm_query_cg_atoms = [symm_query_cg_atoms] # User forgot nested list
    rmsd_threshold = user_data['rmsd_threshold']
    out_dir = user_data['out_dir']
    return [vdgs_dir, cg_name, query_path, bindingsite_residues, query_lig_res, 
           query_cg_atoms, symm_query_cg_atoms, rmsd_threshold, out_dir]

if __name__ == "__main__":
    main()
