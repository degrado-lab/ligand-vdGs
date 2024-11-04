'''
Removes redundant vdGs from the vdG library. vdGs are redundant if they meet all criteria: 
    1. CG and vdM AA identities/positions (backbones within a specified RMSD threshold)
    2. sequence similarity of residues flanking the vdMs (specified by seq. similarity 
       threshold)
    3. similar positions of the flanking residues (CAs within a specified RMSD 
       threshold)
'''

import os
import sys
from itertools import combinations
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import prody as pr

# TODO: convert the following arguments to argparse-compatible 

script, CG = sys.argv

vdg_pdbs_dir = '/wynton/home/degradolab/skt/docking/redun_trial/sulfonamide_tert/vdg_pdbs'
#vdg_pdbs_dir = f'/wynton/group/degradolab/skt/docking/databases/vdg_lib/{CG}/vdg_pdbs' # TODO: remove hardcoding

rmsd_thresh = 0.75
seq_sim_thresh = 0.75
num_flanking = 5

'''
Characterize all subsets of vdMs within each vdG, up to 4 vdMs (i.e., all quadruples, 
triplets, pairs, and singles).

Each key of the `vdm_combos` dict denotes how many residues (vdMs) within the vdG are 
being considered, and each value is a subdict whose keys are the combinations of vdM 
AAs and whose values are a featurized list of the vdG subsets defined by those 
combinations of vdM AAs. The features of the vdG subsets include the CG coordinates, 
the backbone (N, CA, C) coordinates of the vdMs, the sequences flanking those vdMs 
(+/- 5 residues), the PDB from which the vdG is derived, and the PDB segments/chains/
resnums of the CG and vdMs.

Example:

vdm_combos = { 
              4: {
                    (Ala, Ala, Cys, Phe): [
                                             [  [CG coords], 
                                                [bb N-Ca-C coords of Ala1, Ala2, Cys3, Phe4], 
                                                [seq. +/- 5 of Ala1, Ala2, Cys3, Phe4], 
                                                [CA   +/- 5 of Ala1, Ala2, Cys3, Phe4], 
                                                [PDB path],
                                                [PDB seg/chain/resnum of Ala1, Ala2, Cys3, Phe4]
                                             ], 
                                             ...
                                          ], 
                    (Ala, Ala, Cys, Trp): ...
                 }, 
              3: {
                    (Ala, Ala, Cys): ...,
                 } 
             }

Note that with there should be an order-preserving one-to-one mapping of each vdM in the 
lists that enumerate the vdM characteristics. When dealing with multiple vdMs of the same 
AA identity (such as having 2 Ala's in this example), all permutations must be sampled 
when checking for redundancy.
'''

def main():
   # Initialize vdm_combos dict to store the subsets of vdMs within a vdG
   vdm_combos = {}

   # Iterate over the PDBs and CGs that were identified as containing the SMARTS group
   for pdbname in os.listdir(vdg_pdbs_dir):
      pdbpath = os.path.join(vdg_pdbs_dir, pdbname)
      prody_obj = pr.parsePDB(pdbpath)
      cg_coords = get_cg_coords(prody_obj)
      # Characterize the vdM residues (bb coords, flanking residues, pdb paths, etc.)
      vdms_dict = get_vdm_res_features(prody_obj, pdbpath)
      # Determine the vdM combinations, up to 4 residues
      vdm_resinds = list(vdms_dict.keys())
      vdm_subsets = get_vdm_subsets(vdm_resinds)
      # Iterate over subsets
      for vdm_subset in vdm_subsets:
         # Record these features in the same order as in vdm_subset. Then, sort all based on
         # alphabetical order of the vdm AAs.
         re_ordered_aas, re_ordered_bbcoords, re_ordered_flankingseqs, re_ordered_CAs, \
            re_ordered_scrr = reorder_vdg_subset(vdm_subset, vdms_dict)
         # Add to `vdm_combos` dict
         vdm_combos = add_vdgs_to_dict(vdm_combos, vdm_subset, re_ordered_aas, 
            re_ordered_bbcoords, re_ordered_flankingseqs, re_ordered_CAs, re_ordered_scrr,
            cg_coords, pdbpath)
   # Evaluate the complete collection of vdGs and determine redundancy 
   for num_vdms_in_subset, _subsets in vdm_combos.items():
      for _reordered_AAs, _vdgs in _subsets.items():
         # vdG subsets that have identical vdm AA compositions may be redudant 
         if len(_vdgs) <= 1:
            continue
         print()
         print(_reordered_AAs)
         for b in [v[5] for v in _vdgs]:
            print(b)
         compare_vdgs_of_same_AA_comp(_vdgs, rmsd_thresh)

def compare_vdgs_of_same_AA_comp(_vdgs, rmsd_thresh):
   # vdG subsets that have identical vdm AA compositions may be redudant 
   all_cg_coords =    [v[0] for v in _vdgs]
   all_vdm_bbcoords = [v[1] for v in _vdgs]
   all_flankingseqs = [v[2] for v in _vdgs] 
   all_flankingCAs =  [v[3] for v in _vdgs]
   all_pdbpaths =     [v[4] for v in _vdgs]
   all_vdm_scrr =     [v[5] for v in _vdgs]
   all_cg_and_vdmbb_coords = []
   for _cg, _vdmbb in zip(all_cg_coords, all_vdm_bbcoords):
      flattened_cg_coords = flatten_cg_coords(_cg)
      flattened_vdmbbs = flatten_vdg_bbs(_vdmbb)

      cg_and_vdmbb = flattened_cg_coords + flattened_vdmbbs
      cg_and_vdmbb = np.array(cg_and_vdmbb)
      all_cg_and_vdmbb_coords.append(cg_and_vdmbb)
   # Cluster these potential redundant vdGs by rmsd of CG + vdm bb coords. If their rmsds 
   # are larger than the rmsd thershold, then they are automatically considered different 
   # vdGs.
   rmsd_dist_matrix = create_rmsd_dist_matrix(all_cg_and_vdmbb_coords)
   condensed_dist_matrix = squareform(rmsd_dist_matrix)
   Z = linkage(condensed_dist_matrix, method='single')
   clusters = fcluster(Z, rmsd_thresh, criterion='distance')
   cluster_assignments = {} # key = clusnum, value = indices corresponding to that cluster
   for all_vdgs_index, clusnum in enumerate(clusters):
      if clusnum not in cluster_assignments.keys():
         cluster_assignments[clusnum] = []
      cluster_assignments[clusnum].append(all_vdgs_index)
   
def create_rmsd_dist_matrix(coords):
   n = len(coords)
   distance_matrix = np.zeros((n, n))
   for i in range(n):
      for j in range(i + 1, n):
         # Align and then calc rmsd 
         moved_coords_i, transf = pr.superpose(coords[i], coords[j])
         rmsd = pr.calcRMSD(moved_coords_i, coords[j])
         distance_matrix[i][j] = rmsd
         distance_matrix[j][i] = rmsd  # Symmetric matrix
   return distance_matrix

def flatten_cg_coords(_cg):
   coords = []
   for atom in _cg:
      coords.append(atom)
   return coords

def flatten_vdg_bbs(_vdmbb):
   # Flatten backbones to a single list where each element is a numpy array of the x, y, z 
   # coords of each vdm backbone atom
   bbs = []
   for res in _vdmbb:
      for atom in res:
         bbs.append(atom)
   return bbs

def add_vdgs_to_dict(vdm_combos, vdm_subset, re_ordered_aas, re_ordered_bbcoords, 
                     re_ordered_flankingseqs, re_ordered_CAs, re_ordered_scrr, cg_coords, pdbpath):
   # Construct the `vdm_combos` dict by categorizing subset of vdgs based on the number of vdms in
   # each subset, as well as the AA compositions of the vdms. Store each vdg subset as a list
   # containing [CG coords, backbone N-CA-C coords of the vdms, sequences flanking the vdms,
   # CA coordinates flanking the vdms, the PDB path, PDB seg/chain/resnum of the vdms]
   num_vdms_in_subset = len(vdm_subset)
   if num_vdms_in_subset not in vdm_combos.keys():
      vdm_combos[num_vdms_in_subset] = {}
   re_ord_aas_tup = tuple(re_ordered_aas)
   if re_ord_aas_tup not in vdm_combos[num_vdms_in_subset].keys():
      vdm_combos[num_vdms_in_subset][re_ord_aas_tup] = []
   # Add the vdg to `vdm_combos`
   vdm_combos[num_vdms_in_subset][re_ord_aas_tup].append([cg_coords, re_ordered_bbcoords,
      re_ordered_flankingseqs, re_ordered_CAs, pdbpath, re_ordered_scrr])

   return vdm_combos

def reorder_vdg_subset(vdm_subset, vdms_dict):
   #### Step 1) Recording 
   aas_of_vdms_in_order = []
   bb_coords_of_vdms_in_order = []
   flankingseqs_of_vdms_in_order = []
   flanking_CA_coords_of_vdms_in_order = []
   seg_ch_res_of_vdms_in_order = []
   for _vdmresind in vdm_subset:
      vdmAA, vdm_features = vdms_dict[_vdmresind]
      aas_of_vdms_in_order.append(vdmAA)
      vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict = vdm_features
      seg_ch_res_of_vdms_in_order.append(vdm_seg_chain_resnum_resname)
      bb_coords_of_vdms_in_order.append(bb_coords)
      flankingseqs = []
      flankingCAs = []
      sorted_flank_indices = sorted(list(flanking_seq_dict.keys()))
      # Decompress the flanking AA and CA info
      for flank_ind in sorted_flank_indices:
         flank_resname, flank_ca = flanking_seq_dict[flank_ind]
         flankingseqs.append(flank_resname)
         flankingCAs.append(flank_ca)
      flankingseqs_of_vdms_in_order.append(flankingseqs)
      flanking_CA_coords_of_vdms_in_order.append(flankingCAs)
   #### Step 2) Re-ordering based on alphabetical order
   super_list = [aas_of_vdms_in_order, bb_coords_of_vdms_in_order, 
      flankingseqs_of_vdms_in_order, flanking_CA_coords_of_vdms_in_order, 
      seg_ch_res_of_vdms_in_order]
   return sort_vdGs_by_AA(super_list)

def sort_vdGs_by_AA(super_list):
    # Check if all sublists have equal length
    assert all(len(sublist) == len(super_list[0]) for sublist in super_list)
    # Combine the sublists into a list of tuples, where each tuple corresponds to the 
    # elements at the same index
    combined = list(zip(*super_list))
    # Sort the combined list based on the first element (from the first sublist)
    sorted_combined = sorted(combined, key=lambda x: x[0])
    # Unzip the sorted combined list back into sublists
    sorted_sublists = list(zip(*sorted_combined))
    return [list(sublist) for sublist in sorted_sublists]

def get_vdm_res_features(prody_obj, pdbpath):
   # Identify the vdM residues (occ == 2). To be safe, select > 1.5 and < 2.5.
   vdm_residues = prody_obj.select('(occupancy) > 1.5 and (occupancy < 2.5)')
   vdm_resinds = set(vdm_residues.getResindices())
   # Record features of the vdm residues (bb coords, flanking residues, pdb paths, etc.)
   vdms_dict = {}
   for vdm_resind in vdm_resinds:
      vdm_obj = vdm_residues.select(f'resindex {vdm_resind}')
      # BB coords
      bb_coords = get_bb_coords(vdm_obj)

      # Sequence of (contiguous) flanking residues
      flanking_seq_dict = {} # key = relative flank num (-1, +1, etc.), 
                             # value = list(CA coords, AA identity)
      # Walk up and down the flanking residues and check that they are actually
      # neighboring the vdM, and not jumped through a chain break. If there's a
      # chain break, report the AA as "X".
      # > First, store the CA coords.
      for flank_num in range(1, num_flanking + 1):
         negative_flank_num = -1 * flank_num
         positive_flank_num = flank_num
         # Adding the flank_num (neg. or pos.) gives you the residue index
         for f in [negative_flank_num, positive_flank_num]:
            # Get the AA identity and CA coords of the "current" resindex
            current_resindex = vdm_resind + f
            AA, CA_coords = get_AA_and_CA_coords(prody_obj, current_resindex)
            flanking_seq_dict[f] = [AA, CA_coords]
      # Get the CA coords of the central vdM as well
      flanking_seq_dict[0] = ['vdm', bb_coords[1]] # label as 'vdm' for easy exclusion
                                                   # when calculating seq. similarity
      # > Then, go through the dict for chain breaks. Start in the fwd direction, and
      #   then progress backward. 
      fwd_rel_inds = range(1, num_flanking + 1)
      back_rel_inds = [-1 * i for i in fwd_rel_inds]
      central_vdm_CA = flanking_seq_dict[0][1]
      for list_indices in [fwd_rel_inds, back_rel_inds]:
         for ind in list_indices:
            # Get distance between current flanking res (relative to central vdm)
            # and res prior to it.
            if ind == 1 or ind == -1:
               prev_CA = central_vdm_CA
            curr_CA = flanking_seq_dict[ind][1]
            if curr_CA is None:
               flanking_seq_dict = found_chain_break(flanking_seq_dict, ind)
               break 
            dist = pr.calcDistance(np.array(prev_CA), np.array(curr_CA))
            # If the dist is > 4.5A, then it's a chain break.
            if dist > 4.5:
               flanking_seq_dict = found_chain_break(flanking_seq_dict, ind)
               break 
            # Otherwise, continue walking.
            prev_CA = curr_CA
      # Get this vdm resind's PDB identifier (seg, chain, resnum)
      vdm_seg_chain_resnum_resname = get_res_iden(vdm_obj)
      # Store all into dict
      vdm_descript = [vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict]
      vdm_AA = vdm_seg_chain_resnum_resname[-1]
      assert vdm_resind not in list(vdms_dict.keys())
      vdms_dict[vdm_resind] = [vdm_AA, vdm_descript]
   return vdms_dict

def get_res_iden(vdm_obj):
   seg =     list(set(vdm_obj.getSegnames()))
   chain =   list(set(vdm_obj.getChids()))
   resnum =  list(set(vdm_obj.getResnums()))
   resname = list(set(vdm_obj.getResnames()))
   assert len(seg) == 1
   assert len(chain) == 1
   assert len(resnum) == 1
   assert len(resname) == 1
   return seg[0], chain[0], resnum[0], resname[0]

def found_chain_break(flanking_seq_dict, chain_break_ind):
   # If a chain break is found, overwrite all the subsequent (or preceding) flanking
   # residues beyond the chain break as AA "X".
   flank_indices = list(flanking_seq_dict.keys())
   if chain_break_ind > 0 :
      inds_to_overwrite = [n for n in flank_indices if n >= chain_break_ind]
   elif chain_break_ind < 0:
      inds_to_overwrite = [n for n in flank_indices if n <= chain_break_ind]
   
   for overwrite_ind in inds_to_overwrite:
      flanking_seq_dict[overwrite_ind] = ['X', None]
   return flanking_seq_dict

def get_AA_and_CA_coords(prody_obj, current_resindex):
   if current_resindex < 0: 
      sel_str = f'resindex `{current_resindex}`'
   else:
      sel_str = f'resindex {current_resindex}'

   curr_resindex_obj = prody_obj.select(sel_str)
   if curr_resindex_obj is None:
      AA = 'X'
      CA_coords = None
   elif curr_resindex_obj.protein is None: # not a residue
      AA = 'X'
      CA_coords = None
   else:
      curr_res_AA = list(set(curr_resindex_obj.getResnames()))
      assert len(curr_res_AA) == 1
      AA = curr_res_AA[0]
      CA_obj = curr_resindex_obj.select(sel_str + ' and name CA')
      assert len(CA_obj) == 1
      CA_coords = CA_obj.getCoords()[0]
   return AA, CA_coords

def get_cg_coords(prody_obj):
   # The CG atoms have their occupancies set to >= 3.0, with unique values (e.g., 3.0, 
   # 3.1, 3.2, etc.) to allow a 1:1 correspondence of equivalent atoms between different
   # ligands.
   cg = prody_obj.select('occupancy >= 3.0')
   num_atoms = len(cg)
   cg_coords = []
   for ind in range(num_atoms):
      occ = f'3.{ind}'
      atom = cg.select(f'occupancy == {occ}')
      assert len(atom) == 1
      atom_coords = atom.getCoords()[0]
      cg_coords.append(atom_coords)

   return np.array(cg_coords)

def get_bb_coords(obj):
   bb_coords = []
   for atom in ['N', 'CA', 'C']:
      atom_obj = obj.select(f'name {atom}')
      assert len(atom_obj) == 1
      coord = atom_obj.getCoords()[0]
      bb_coords.append(coord)
   return bb_coords

def get_vdm_subsets(input_list):
    # Group by quadruples, triples, pairs, and singles.
    # Initialize an empty list to store all subsets
    all_subsets = []
    # Loop through subset sizes 1 to 4, generate combos, then add to list
    for r in range(1, 5):
        subsets = combinations(input_list, r)
        all_subsets.extend(subsets)
    return all_subsets

if __name__ == "__main__":
    main()
