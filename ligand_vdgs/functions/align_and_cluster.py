from itertools import permutations, combinations, product
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import prody as pr
import os

def get_vdg_AA_permutations(reordered_AAs, _vdgs):
   permuted_indices = permute_AA_duplicates(reordered_AAs)
   all_AA_cg_perm_cg_coords = []
   all_AA_cg_perm_vdm_bbcoords = []
   all_AA_cg_perm_flankingseqs = []
   all_AA_cg_perm_flankingCAs = []
   all_AA_cg_perm_pdbpaths = []
   all_AA_cg_perm_vdm_scrr = []

   # Iterate over all AA permutations of each vdg
   for _vdg in _vdgs:
      # CG coords and pdbpaths remain unchanged, but vdmbbs, flankingseqs, 
      # flankingCAs, and scrrs need to be permuted.
      for permutation in permuted_indices:
         all_AA_cg_perm_cg_coords.append(_vdg[0])
         all_AA_cg_perm_pdbpaths.append(_vdg[4])
         nonpermuted_vdmbb = _vdg[1]
         nonpermuted_flankingseqs = _vdg[2]
         nonpermuted_flankingCAs = _vdg[3]
         nonpermuted_vdm_scrr = _vdg[5]
         vdmbb_permutation = [nonpermuted_vdmbb[ix] for ix in permutation]
         flankingseqs_permutation = [nonpermuted_flankingseqs[ix] for ix in 
                                     permutation]
         flankingCAs_permutation = [nonpermuted_flankingCAs[ix] for ix in permutation]
         vdm_scrrs_permutation = [nonpermuted_vdm_scrr[ix] for ix in permutation]
         all_AA_cg_perm_vdm_bbcoords.append(vdmbb_permutation)
         all_AA_cg_perm_flankingseqs.append(flankingseqs_permutation)
         all_AA_cg_perm_flankingCAs.append(flankingCAs_permutation)
         all_AA_cg_perm_vdm_scrr.append(vdm_scrrs_permutation)

   return all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords, \
      all_AA_cg_perm_flankingseqs, all_AA_cg_perm_flankingCAs, \
         all_AA_cg_perm_pdbpaths, all_AA_cg_perm_vdm_scrr

def permute_AA_duplicates(seq):
    # Dictionary to store indices for each element in the sequence
    seen = {}
    for i, item in enumerate(seq):
        if item not in seen:
            seen[item] = []
        seen[item].append(i)

    permute_groups = [] # list of all index groups that have duplicates
    for indices in seen.values():
        if len(indices) > 1:  # only interested in elements with duplicates
            permute_groups.append(indices)

    if not permute_groups: # no duplicates, so return the original index list
        return [list(range(len(seq)))]
    
    # generate all possible permutations of indices within each group of duplicates
    permuted_idx_lists = []
    for perm_combination in product(*[permutations(group) for group in 
                                      permute_groups]):
        # start with the list of original indices
        permuted_idx = list(range(len(seq)))

        # flatten the product of permutations and assign them to the corresponding 
        # positions
        for group_idx, perm in zip(permute_groups, perm_combination):
            for orig_idx, new_idx in zip(group_idx, perm):
                permuted_idx[orig_idx] = new_idx

        permuted_idx_lists.append(permuted_idx)
    
    return permuted_idx_lists

def combine_cg_and_vdmbb_coords(all_AA_cg_perm_cg_coords, 
                                all_AA_cg_perm_vdm_bbcoords):
   all_cg_and_vdmbb_coords = []
   for _cg, _vdmbb_per_res in zip(all_AA_cg_perm_cg_coords, 
                                  all_AA_cg_perm_vdm_bbcoords):
      flattened_cg_coords = np.array(flatten_cg_coords(_cg))
      flattened_vdmbbs =    np.array(flatten_vdg_bbs(_vdmbb_per_res))
      cg_and_vdmbb = np.vstack([flattened_cg_coords, flattened_vdmbbs])
      all_cg_and_vdmbb_coords.append(cg_and_vdmbb)
   return all_cg_and_vdmbb_coords

def flatten_cg_coords(_cg):
   coords = []
   for atom in _cg:
      coords.append(atom)
   return coords

def flatten_vdg_bbs(_vdmbb):
   # Flatten backbones to a single list where each element is a numpy array of the 
   # x, y, z coords of each vdm backbone atom
   bbs = []
   for res in _vdmbb:
      for atom in res:
         bbs.append(atom)
   return bbs

def get_hierarchical_clusters(data, rmsd_cut,  calc_seq=False, calc_rmsd=False, 
                              seq_sim_thresh=None):
   # Returns dict where key = cluster number and value = indices of the list elements 
   # that belong in that cluster. This dict is not ordered by size.
   # `data` is either a list of coords, or a list of sequences.
   if sum([calc_seq, calc_rmsd]) != 1:
      raise ValueError('Either calc_seq or calc_rmsd, but not both, must be True.')
   
   dist_matrix = create_dist_matrix(data, calc_seq, calc_rmsd) 
   condensed_dist_matrix = squareform(dist_matrix)
   
   # If the distancee matrix is empty (because there's only one sample in `data`), then
   # linkage() cannot be calculated on it. Return a single cluster with the only sample
   # as its centroid.
   if len(condensed_dist_matrix) == 0:
       return {1: [0]}, {1: 0} 
   
   Z = linkage(condensed_dist_matrix, method='single')
   if seq_sim_thresh:
      seq_dissimilarity_thresh = 1 - seq_sim_thresh
      clusters = fcluster(Z, seq_dissimilarity_thresh, criterion='distance')
   else:
      clusters = fcluster(Z, rmsd_cut, criterion='distance')
   cluster_assignments = {} # key = clusnum, value = indices of `coords` 
                            # belonging to that cluster
   for all_vdgs_index, clusnum in enumerate(clusters):
      if clusnum not in cluster_assignments.keys():
         cluster_assignments[clusnum] = []
      cluster_assignments[clusnum].append(all_vdgs_index)
   
   # Identify and label centroids
   centroids = {} # key = clusnum, val = index of `coords` corresponding
              # to the centroid
   for clusnum, indices in cluster_assignments.items():
      cluster_data = [data[i] for i in indices]
      cluster_dist_matrix = dist_matrix[np.ix_(indices, indices)]
      centroid_index = np.argmin(cluster_dist_matrix.sum(axis=0))
      centroids[clusnum] = indices[centroid_index]
   
   return cluster_assignments, centroids

def create_dist_matrix(data, calc_seq=False, calc_rmsd=True):
   n = len(data)
   distance_matrix = np.zeros((n, n))
   for i in range(n):
      for j in range(i + 1, n):
         # Align and then calc rmsd or seq dissimilarity. If "None" is in the list 
         # of coords (which is what would happen if rmsd is being calculated on 
         # stretches of backbone and terminal residues are missing), then select 
         # only the overlapping residues where a coordinate could be extracted for 
         # both elements being measured.
         data_i = data[i]
         data_j = data[j]
         if calc_rmsd:  
            data_i, data_j = get_overlapping_res_coords(data_i, data_j)
            moved_coords_i, transf = pr.superpose(data_i, data_j)
            rmsd = pr.calcRMSD(moved_coords_i, data_j)
            value = rmsd
         elif calc_seq: # calc seq similarity
            seq_sim = calc_seq_similarity(data_i, data_j)
            value = (100 - seq_sim) / 100 # clustering is distance-based, so the 
                                          # metric is DISsimilarity
         distance_matrix[i][j] = value 
         distance_matrix[j][i] = value   # Symmetric matrix
   return distance_matrix

def get_overlapping_res_coords(list1, list2):
    list1_overlap = []
    list2_overlap = []
    
    # Loop over both lists and check if both elements at the same index are not None
    for i in range(len(list1)):
        if list1[i] is not None and list2[i] is not None:
            list1_overlap.append(list1[i])
            list2_overlap.append(list2[i])
     
    return np.array(list1_overlap), np.array(list2_overlap)

def calc_seq_similarity(list1, list2):
    # Count how many residues match
    list1 = [i for i in list1 if i != 'vdm']
    list2 = [i for i in list2 if i != 'vdm']
    assert len(list1) == len(list2)
    # Exclude residues if at least one of them is an "X"
    indices_to_exclude = []
    for ind in range(len(list1)):
         if list1[ind] == 'X' or list2[ind] == 'X':
            indices_to_exclude.append(ind)
    list1 = [item for idx, item in enumerate(list1) if idx not in indices_to_exclude]
    list2 = [item for idx, item in enumerate(list2) if idx not in indices_to_exclude]
    matches = sum(1 for a, b in zip(list1, list2) if a == b)

    # Calculate the percentage of matches
    if len(list1) == 0:
        match_percentage = 0
    else:
        match_percentage = (matches / len(list1)) * 100
    return match_percentage

def get_vdm_res_features(prody_obj, pdbpath, num_flanking):
   # Identify the vdM residues (occ == 2). To be safe, select > 1.5 and < 2.5.
   vdm_residues = prody_obj.select('(occupancy) > 1.5 and (occupancy < 2.5)')
   vdm_resinds = set(vdm_residues.getResindices())
   # Record features of the vdm residues (bb coords, flanking residues, pdb paths, 
   # etc.)
   vdms_dict = {}
   for vdm_resind in vdm_resinds:
      vdm_obj = vdm_residues.select(f'resindex {vdm_resind}')
      # BB coords
      bb_coords = get_bb_coords(vdm_obj)
      if bb_coords is None:
         continue

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
   # The CG atoms have their occupancies set to >= 3.0, with unique values (e.g., 
   # 3.0, 3.1, 3.2, etc.) to allow a 1:1 correspondence of equivalent atoms between 
   # different ligands.
   cg = prody_obj.select('occupancy >= 3.0')
   num_atoms = len(cg)
   cg_coords = []
   for ind in range(num_atoms):
      occ = f'3.{ind}'
      atom = cg.select(f'occupancy == {occ}')
      if len(atom) != 1:
         print(f'get_cg_coords: {len(atom)} atoms are occupancy {occ}.')
         return None
      assert len(atom) == 1
      atom_coords = atom.getCoords()[0]
      cg_coords.append(atom_coords)

   return np.array(cg_coords)

def get_bb_coords(obj):
   bb_coords = []
   for atom in ['N', 'CA', 'C']:
      atom_obj = obj.select(f'name {atom}')
      if atom_obj is None:
         print('get_bb_coords: atom_obj is None.')
         return None
      if len(atom_obj) > 1: # not sure why, but sometimes each atom is duplicated,
                            # despite having identical coords.
         ambiguous_coords = None
         for atom in atom_obj:
            if ambiguous_coords is None:
               ambiguous_coords = atom.getCoords()
            else:
               dist = pr.calcDistance(ambiguous_coords, atom.getCoords())
               if dist > 0.2:
                  print('get_bb_coords: >1 atom defined in atom_obj.')
                  return None
      coord = atom_obj.getCoords()[0]
      bb_coords.append(coord)
   return bb_coords

def get_vdg_subsets(input_list):
    # Group by quadruples, triples, pairs, and singles.
    # Initialize an empty list to store all subsets
    all_subsets = []
    # Loop through subset sizes 1 to 4, generate combos, then add to list
    for r in range(1, 5):
        subsets = combinations(input_list, r)
        all_subsets.extend(subsets)
    return all_subsets

def add_vdgs_to_dict(vdm_combos, vdg_subset, re_ordered_aas, re_ordered_bbcoords, 
                     re_ordered_flankingseqs, re_ordered_CAs, re_ordered_scrr, 
                     cg_coords, pdbpath):
   # Construct the `vdm_combos` dict by categorizing subset of vdgs based on the 
   # number of vdms in each subset, as well as the AA compositions of the vdms. 
   # Store each vdg subset as a list containing [CG coords, backbone N-CA-C coords 
   # of the vdms, sequences flanking the vdms, CA coordinates flanking the vdms, 
   # the PDB path, PDB seg/chain/resnum of the vdms]
   num_vdms_in_subset = len(vdg_subset)
   if num_vdms_in_subset not in vdm_combos.keys():
      vdm_combos[num_vdms_in_subset] = {}
   re_ord_aas_tup = tuple(re_ordered_aas)
   if re_ord_aas_tup not in vdm_combos[num_vdms_in_subset].keys():
      vdm_combos[num_vdms_in_subset][re_ord_aas_tup] = []
   # Add the vdg to `vdm_combos`
   vdm_combos[num_vdms_in_subset][re_ord_aas_tup].append([cg_coords, 
      re_ordered_bbcoords, re_ordered_flankingseqs, re_ordered_CAs, pdbpath, 
      re_ordered_scrr])

   return vdm_combos

def reorder_vdg_subset(vdg_subset, vdms_dict):
   '''Reorders vdms alphabetically.'''
   # First, record
   aas_of_vdms_in_order = []
   bb_coords_of_vdms_in_order = []
   flankingseqs_of_vdms_in_order = []
   flanking_CA_coords_of_vdms_in_order = []
   seg_ch_res_of_vdms_in_order = []
   for _vdmresind in vdg_subset:
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
   # Then, re-order based on alphabetical order
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

def permute_on_indices(symmetry_classes, coords_to_permute):
    # Given a list of indices to permute on (symmetry_classes), permute a list or 
    # array of coordinates.

   assert len(symmetry_classes) == len(coords_to_permute)

   # Group coordinates by their indices
   grouped_coords = {}
   for idx, coord in zip(symmetry_classes, coords_to_permute):
       if idx not in grouped_coords:
           grouped_coords[idx] = []
       grouped_coords[idx].append(coord)

   # Generate all permutations of coordinates within each group
   grouped_permutations = {}
   for key in grouped_coords:
       grouped_permutations[key] = list(permutations(grouped_coords[key]))

   # Generate all combinations of permutations (one for each group)
   result = []
   for perm_combination in product(*grouped_permutations.values()):
       # Rebuild the permuted list from the permuted groups
       permuted_list = []
       # For each original index, we need to pick the corresponding permuted item
       perm_idx = {key: 0 for key in grouped_permutations}  # Keep track of which 
                                            # element in each permutation we're at
       for idx in symmetry_classes:
           # Find the group corresponding to the current index (grouped_permutations)
           group_perm = perm_combination[list(
              grouped_permutations.keys()).index(idx)]
           permuted_list.append(group_perm[perm_idx[idx]])
           perm_idx[idx] += 1  # Move to the next item in this group's permutation
       
       result.append(np.array(permuted_list))

   return result

def get_vdg_AA_and_cg_perms(all_AA_perm_cg_coords,    all_AA_perm_vdm_bbcoords, 
                            all_AA_perm_flankingseqs, all_AA_perm_flankingCAs, 
                            all_AA_perm_pdbpaths,     all_AA_perm_vdm_scrr, 
                            symmetry_classes):
   if symmetry_classes is None:
      return [all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords,          
              all_AA_perm_flankingseqs, all_AA_perm_flankingCAs, 
              all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr]

   # Recompile the vdgs with the permuted symmetric cg indices
   all_AA_cg_perm_cg_coords = []
   all_AA_cg_perm_vdm_bbcoords = []
   all_AA_cg_perm_flankingseqs = []
   all_AA_cg_perm_flankingCAs = []
   all_AA_cg_perm_pdbpaths = []
   all_AA_cg_perm_vdm_scrr_cg_perm = []
   
   for cgcoords, vdmcoords, seq, CAs, pdbpath, scrr in zip(all_AA_perm_cg_coords,
            all_AA_perm_vdm_bbcoords, all_AA_perm_flankingseqs, 
            all_AA_perm_flankingCAs, all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr):
      
      perms = permute_on_indices(symmetry_classes, cgcoords)

      for perm_ind, perm in enumerate(perms):
         all_AA_cg_perm_cg_coords.append(perm)
         all_AA_cg_perm_vdm_bbcoords.append(vdmcoords)
         all_AA_cg_perm_flankingseqs.append(seq)
         all_AA_cg_perm_flankingCAs.append(CAs)
         all_AA_cg_perm_pdbpaths.append(pdbpath)
         # add permutation index to the scrr label
         vdm_scrr_perm = [scrr, f'cg_perm_{perm_ind+1}']
         all_AA_cg_perm_vdm_scrr_cg_perm.append(vdm_scrr_perm)
   
   return [all_AA_cg_perm_cg_coords,    all_AA_cg_perm_vdm_bbcoords, 
           all_AA_cg_perm_flankingseqs, all_AA_cg_perm_flankingCAs, 
           all_AA_cg_perm_pdbpaths,     all_AA_cg_perm_vdm_scrr_cg_perm]
   
def elements_in_clusters(indices_of_elements_in_cluster, cg_coords, vdm_bbcoords,
                         cgvdmbb_coords, flankingseqs, flankingCAs, pdbpaths, 
                         vdm_scrrs_cg_perms):
   assert len(cg_coords) == len(pdbpaths)
   clus_cg_coords = [cg_coords[index] for index in indices_of_elements_in_cluster]
   clus_vdmbb_coords = [vdm_bbcoords[index] for index in 
                        indices_of_elements_in_cluster]
   clus_cgvdmbb_coords = [cgvdmbb_coords[index] for index in 
                        indices_of_elements_in_cluster]
   clus_flankingseqs = [flankingseqs[index] for index in 
                        indices_of_elements_in_cluster]
   clus_flankingCAs = [flankingCAs[index] for index in indices_of_elements_in_cluster]
   clus_pdbpaths = [pdbpaths[index] for index in indices_of_elements_in_cluster]
   clus_vdm_scrr_cg_perms = [vdm_scrrs_cg_perms[index] for index in 
                    indices_of_elements_in_cluster]
   return clus_cg_coords, clus_vdmbb_coords, clus_cgvdmbb_coords, clus_flankingseqs, \
      clus_flankingCAs, clus_pdbpaths, clus_vdm_scrr_cg_perms

def write_out_clusters(clusdir, clus_assignments, centroid_assignments, all_cg_coords, 
                       all_pdbpaths, all_scrr_cg_perm, all_cg_and_vdmbb_coords, 
                       all_flankbb_coords, num_flanking, first_pdb_out, 
                       first_pdb_cg_vdmbb_coords, weights, atomgroup_dict): 

   assert len(all_pdbpaths) == len(all_cg_and_vdmbb_coords)
   ref = np.array([[0, 0, 0], [-1, 0, 1], [1, -1, 0]]) # it's just to align the 
        # first 3 atoms of CG. then, all atoms of all subsequent CGs will be 
        # aligned to that first CG.
   # Iterate over each cluster
   failed = []
   for clusnum, clus_mem_indices in \
                clus_assignments.items():
      clusnum_dir = os.path.join(clusdir, f'clus_{clusnum}')
      os.makedirs(clusnum_dir)
      centroid_ind = centroid_assignments[clusnum]
      # First, output the centroid before the other cluster members. Align centroids
      # on cg+vdmbb.
      data = get_clus_mem_data(centroid_ind, all_cg_coords, 
         all_cg_and_vdmbb_coords, all_flankbb_coords, all_pdbpaths, all_scrr_cg_perm, 
         num_flanking, atomgroup_dict)
      if data is None:
         failed.append(all_pdbpaths[centroid_ind])
         continue
      (cent_cg_coords, cent_cg_vdmbb_coords, cent_flankbb_coords, cent_pdbpath, 
         cent_scrr_cg_perm, cent_pr_obj) = data
      cent_pdb_outpath = get_clus_pdb_outpath(clusnum, clusnum_dir, cent_pdbpath, 
                                       cent_scrr_cg_perm, centroid_ind, is_centroid=True)
      if first_pdb_out is None:
         # If no pdb has been written out yet, align the first 3 CG atoms to the
         # reference 3 atoms and output the pdb.
         first_pdb_out, first_pdb_cg_vdmbb_coords, centroid_transf = print_out_first_pdb_of_clus(
            cent_pdb_outpath, cent_cg_coords, cent_cg_vdmbb_coords, cent_pr_obj, ref)
      else:
         # If a pdb has already been written out, align this vdg's cg+vdmbb onto that 
         # first pdb. "target_coords" to align to is cg+vdmbb of the first PDB (which 
         # itself is a centroid). Return the transformation b/c this moved centroid 
         # will be the "target" for all the other members of the cluster.
         target_coords = first_pdb_cg_vdmbb_coords
         centroid_transf = write_out_subsequent_clus_pdbs(cent_cg_vdmbb_coords, 
            cent_pr_obj, cent_pdb_outpath, target_coords, weights)
      target_coords = [] # transform the bb atoms individually bc some have None
      for b_ in cent_flankbb_coords:
         if b_ is None:
            target_coords.append(None)
         else:
            transformed_b = pr.applyTransformation(centroid_transf, b_)
            assert len(transformed_b) == 1
            # b_ is shape (3,) and transformed_b is shape (1, 3) but they need to match
            transformed_b = transformed_b[0] 
            target_coords.append(transformed_b)

      for ind in clus_mem_indices:
         if ind == centroid_ind:
            continue
         data = get_clus_mem_data(
            ind, all_cg_coords, all_cg_and_vdmbb_coords, all_flankbb_coords, all_pdbpaths, 
            all_scrr_cg_perm, num_flanking, atomgroup_dict)
         if data is None:
            failed.append(all_pdbpaths[ind])
            continue
         (clusmem_cg_coords, clusmem_cg_vdmbb_coords, clusmem_flankbb_coords, 
            clusmem_pdbpath, clusmem_scrr_cg_perm, clusmem_pr_obj) = data

         # Define atoms to align. 
         clusmem_coords = clusmem_flankbb_coords

         clusmem_pdb_outpath = get_clus_pdb_outpath(clusnum, clusnum_dir, 
                        clusmem_pdbpath, clusmem_scrr_cg_perm, ind, is_centroid=False)
         write_out_subsequent_clus_pdbs(clusmem_coords, clusmem_pr_obj, 
                           clusmem_pdb_outpath, target_coords, weights)
   
   return first_pdb_out, first_pdb_cg_vdmbb_coords, failed

def write_out_subsequent_clus_pdbs(clusmem_coords, clusmem_pr_obj, clusmem_pdb_outpath, 
                                   target_coords, weights=None):
   if weights is None:
      weights = np.array([1/len(clusmem_coords) for i in clusmem_coords])

   transf = get_transformation_for_flankbb(clusmem_coords, target_coords)
   pr_obj_copy = clusmem_pr_obj.copy() # b/c of mutability
   pr.applyTransformation(transf, pr_obj_copy)
   pr.writePDB(clusmem_pdb_outpath, pr_obj_copy)
   return transf

def get_transformation_for_flankbb(flankbb_mob, flankbb_tar):
   # Treat differently from cg+vdmbb because this is a list of coords that 
   # occassionally contains Nones.
   data_i, data_j = get_overlapping_res_coords(flankbb_mob, flankbb_tar)
   moved_coords_i, transf = pr.superpose(data_i, data_j)
   return transf

def get_clus_pdb_outpath(clusnum, clusnum_dir, pdbpath, scrr_cg_perm, clusmem_ind, 
                         is_centroid):
   pdbname = os.path.basename(pdbpath).removesuffix('.pdb')
   scrrs, perm = scrr_cg_perm
   perm = perm.removeprefix('cg_perm_')
   vdg_scrr_str = '_'.join(['_'.join([str(z) for z in v_s]) for v_s in scrrs])
   if is_centroid:
      pdb_str = f'clus{clusnum}_{pdbname}_{vdg_scrr_str}_CGperm{perm}_ix{clusmem_ind}_centroid.pdb.gz'
   else:
      pdb_str = f'clus{clusnum}_{pdbname}_{vdg_scrr_str}_CGperm{perm}_ix{clusmem_ind}.pdb.gz'

   return os.path.join(clusnum_dir, pdb_str)

def write_out_first_pdb(mobile, target, pr_obj, outpath, clusmem_cg_coords,
                        clusmem_cg_vdmbb_coords):
   # No weights because it's only aligning 3 atoms of CG
   mobile, target = np.array(mobile), np.array(target)
   moved_3atom_coords, transf = pr.superpose(mobile, target)
   pr_obj_copy = pr_obj.copy() # b/c of mutability
   pr.applyTransformation(transf, pr_obj_copy)
   moved_cg_vdmbb_coords = pr.applyTransformation(transf, clusmem_cg_vdmbb_coords)
   pr.writePDB(outpath, pr_obj_copy)
   # superposed on 3 reference atoms, but need to return coords of the entire CG.
   moved_cg_coords = pr.applyTransformation(transf, clusmem_cg_coords)
   return moved_cg_coords, moved_cg_vdmbb_coords, transf

def rewrite_temp_clusters(clusdir, clean_dir, size_subset, subset_AAs):
   # Clean up the cluster directory by merging degenerate vdGs (based on diff
   # AA perms and CG perms of the same PDB), deleting degenerate pdbs/clusters, 
   # reassigning cluster nums (and sort by cluster size, and removing the temp dir
   # when everything is done (outside this function). 
   # Must merge before deleting duplicates b/c need the duplicate names to determine 
   # which clusters are equivalent.
   
   temp_clus_assignments = {}
   ''' Structure:
   temp_clus_assignments[tuple([flankseq clusnum, flankbb clusnum, 
      cgvdmbb clusnum])] = ALL vdgs w/ sorted scrrs. Don't remove identical vdgs yet.
   '''
   clean_dir = os.path.join(clean_dir, str(size_subset), subset_AAs)
   clusdir = os.path.join(clusdir, str(size_subset), subset_AAs)
   for flankseq_clus in os.listdir(clusdir):
      flankseq_clus_dir = os.path.join(clusdir, flankseq_clus)
      for flankbb_clus in os.listdir(flankseq_clus_dir):
         flankbb_clus_dir = os.path.join(flankseq_clus_dir, flankbb_clus)
         for cgvdmbb_clus in os.listdir(flankbb_clus_dir):
            cgvdmbb_clus_dir = os.path.join(flankbb_clus_dir, cgvdmbb_clus)
            for pdb in os.listdir(cgvdmbb_clus_dir):
               pdbpath = os.path.join(cgvdmbb_clus_dir, pdb)
               temp_clus_assignments = get_temp_clus_assignments(flankseq_clus,
                  flankbb_clus, cgvdmbb_clus, pdb, pdbpath, temp_clus_assignments)

   temp_clus_assignments = merge_equivalent_clusters(temp_clus_assignments)
   pruned_clusters = delete_redun_vdgs(temp_clus_assignments)
   reassigned = reassign_temp_clusters(pruned_clusters)

   # Move the deduplicated and reassigned PDBs to the new clus dir.
   for new_seq_clusnum, new_bb_clusters in reassigned.items():
      for new_bb_clusnum, new_cgvdmbb_clusters in new_bb_clusters.items():
         for new_cgvdmbb_clusnum, vdgs in new_cgvdmbb_clusters.items():
            reassigned_clusdir = os.path.join(clean_dir, 
               f'flankseqclus_{new_seq_clusnum}', f'flankbbclus_{new_bb_clusnum}', 
               f'cgvdmbbclus_{new_cgvdmbb_clusnum}')
            os.makedirs(reassigned_clusdir)
            
            # Copy over the vdgs
            for vdg in vdgs:
               pdbbase, vdmscrrs, ix, vdgpath = vdg
               # rename the pdb for clarity
               scrr_str = '_'.join(['_'.join([str(z) for z in v_s]) for v_s in vdmscrrs])
               if 'centroid' in vdgpath:
                  new_pdbname = f'{pdbbase}_{scrr_str}_centroid.pdb.gz'
               else:
                  new_pdbname = f'{pdbbase}_{scrr_str}.pdb.gz'
               assert not os.path.exists(new_pdbname)
               os.rename(vdgpath, os.path.join(reassigned_clusdir, new_pdbname))
   
   return reassigned
         
def determine_redundant_temp_vdg(already_seen_temp_vdgs, 
                                 candidate_pdbbase, candidate_scrr):
   # Return True if the candidate vdg has already been seen.
   for a in already_seen_temp_vdgs:
      seen_pdbbase = a[0]
      seen_scrr = a[1]
      if candidate_pdbbase == seen_pdbbase and candidate_scrr == seen_scrr:
         return True
   return False

def determine_cluster_redundancy(clusnumA_vdgs, clusnumB_vdgs): 
   # Return yes if the vdg in clus A is in clus B
   for vdgA in clusnumA_vdgs:
      vdgA_pdbbase, vdgA_sorted_scrrs, vdgA_ix, vdgA_pdbpath = vdgA
      for vdgB in clusnumB_vdgs:
         vdgB_pdbbase, vdgB_sorted_scrrs, vdgB_ix, vdgB_pdbpath = vdgB
         if vdgA_pdbbase == vdgB_pdbbase and vdgA_sorted_scrrs == vdgB_sorted_scrrs:
            return True
   return False

def get_temp_clus_assignments(flankseq_clus, flankbb_clus, cgvdmbb_clus, pdb, 
                              pdbpath, temp_clus_assignments):
   pdb_base = pdb.split('.pdb.gz')[0]
   pdb_base = '_'.join(pdb_base.split('_')[1:]) # removes "clus{x}_" from name
   scrrs = pdb.split('.pdb.gz_')[1].removesuffix('.pdb.gz')
   scrrs = scrrs.split('_CGperm')[0]
   scrrs = scrrs.split('_')
   assert len(scrrs) % 4 == 0
   grouped_scrrs = [scrrs[i:i+4] for i in range(0, len(scrrs), 4)]
   sorted_scrrs = sorted(grouped_scrrs)

   clus_mem_index = pdb.split('_ix')[1].split('_')[0].removesuffix(
      '.pdb.gz') # index for this specific pdb of specific AA perm and cg perm for 
                 # mapping it back to the original "all_AA_cg_perm...etc..." lists. 
   vdg = tuple([pdb_base, sorted_scrrs, clus_mem_index, pdbpath]) # determine 
      # redundancy just on pdb_base and sorted_scrrs, but record pdb path so that 
      # one file can be referenced for copying over to a clean cluster dir.
   # Check for redundancy (does this vdG already exist in temp_clus_assignments?)
   clus_tuple = tuple([flankseq_clus, flankbb_clus, cgvdmbb_clus])
   if clus_tuple not in temp_clus_assignments.keys():
      temp_clus_assignments[clus_tuple] = []
   existing_in_clus = [i[:2] for i in temp_clus_assignments[clus_tuple]]
   assert vdg[:2] not in existing_in_clus
   temp_clus_assignments[clus_tuple].append(vdg)
   return temp_clus_assignments

def merge_equivalent_clusters(temp_clus_assignments):
   # After gathering all the pdb and sorted scrrs, find the degenerate vdGs and merge them.
   already_seen_temp_clustup_A = []
   clusnumbers = list(temp_clus_assignments.keys())
   for clustupA in clusnumbers:
      if clustupA not in temp_clus_assignments.keys(): # b/c dict dynamically changing
         continue
      clusnumA_vdgs = temp_clus_assignments[clustupA]
      already_seen_temp_clustup_A.append(clustupA)
      for clustupB in clusnumbers:
         if clustupA == clustupB:
            continue
         if clustupB not in temp_clus_assignments.keys(): # b/c dict dynamically changing
            continue
         if clustupB in already_seen_temp_clustup_A:
            continue # ensures that only the upper triangle is calculated 
         clusnumB_vdgs = temp_clus_assignments[clustupB]
         # Are any of clusnumA vdgs in clusnumB?
         if determine_cluster_redundancy(clusnumA_vdgs, clusnumB_vdgs): 
            # If so, merge the clusters
            temp_clus_assignments[clustupA] += temp_clus_assignments[clustupB]
            del temp_clus_assignments[clustupB]
   return temp_clus_assignments

def delete_redun_vdgs(temp_clus_assignments):
   # Next, delete the redundant vdgs within the merged clusters
   already_seen_temp_vdgs = []
   pruned_clusters = {}
   for clusnum, vdgs in temp_clus_assignments.items():
      for _ind, vdg in enumerate(vdgs):
         pdbbase, scrr, clusmem_ix, pdbpath = vdg
         if determine_redundant_temp_vdg(already_seen_temp_vdgs, pdbbase, scrr):
            continue
         else:
            already_seen_temp_vdgs.append(vdg)
            # add to pruned_clusters
            if clusnum not in pruned_clusters.keys():
               pruned_clusters[clusnum] = []
            pruned_clusters[clusnum].append(vdg)
   return pruned_clusters

def reassign_temp_clusters(pruned_clusters):
   '''
   Reassign cluster numbers, sorting cluster numbers based on size.
   Each subdict is re-ordered on the "clustering" level (meaning flank seq clusters
   are reordered first, followed by flank bb clusters, followed by cgvdmbb clusters).
   '''
   clus_tuples = list(pruned_clusters.keys())
   original_clus_labels = {}
   for tup in clus_tuples:
      seq_clus, bb_clus, cgvdmbb_clus = tup
      if seq_clus not in original_clus_labels.keys():
         original_clus_labels[seq_clus] = {}
      if bb_clus not in original_clus_labels[seq_clus].keys():
         original_clus_labels[seq_clus][bb_clus] = {}
      assert cgvdmbb_clus not in original_clus_labels[seq_clus][bb_clus].keys()
      original_clus_labels[seq_clus][bb_clus][cgvdmbb_clus] = pruned_clusters[tup]
   
   sorted_seq_clusters = sorted(original_clus_labels.items(), 
                                key=lambda x: len(x[1]), reverse=True)
   
   # Re-order/renumber seq clus nums based on # of flank bb clus nums
   reassigned_seq_clusters = {}
   for new_seq_clus_num, (seq_clusters, bb_clusters) in enumerate(
                        sorted_seq_clusters, start=1):
      reassigned_seq_clusters[new_seq_clus_num] = bb_clusters

   # Re-order/renumber flank bb clus nums based on # of cgvdmbb clus nums 
   reassigned_bb_clusters = {}
   for new_seq_clus_num, orig_bb_clusters in reassigned_seq_clusters.items():
      sorted_bb_clusters = sorted(orig_bb_clusters.items(), 
                                   key=lambda x: len(x[1]), reverse=True)
      reassigned_bb_clusters[new_seq_clus_num] = {}
      for new_bb_clus_num, (bb_clus_, cgvdmbb_clusters) in enumerate(
                           sorted_bb_clusters, start=1):
         reassigned_bb_clusters[new_seq_clus_num][new_bb_clus_num] = cgvdmbb_clusters
   
   # Re-order/renumber cgvdmbb clus nums based on # of vdgs
   reassigned_cgvdmbb_clusters = {}
   for new_seq_clus_num, orig_bb_clusters in reassigned_bb_clusters.items():
      reassigned_cgvdmbb_clusters[new_seq_clus_num] = {}
      for new_bb_clus_num, orig_cgvdmbb_clusters in orig_bb_clusters.items():
         sorted_cgvdmbb_clusters = sorted(orig_cgvdmbb_clusters.values(), 
                                          key=lambda x: len(x), reverse=True)
         reassigned_cgvdmbb_clusters[new_seq_clus_num][new_bb_clus_num] = {}
         for new_cgvdmbb_clus_num, vdgs in enumerate(sorted_cgvdmbb_clusters, start=1):
            reassigned_cgvdmbb_clusters[new_seq_clus_num][new_bb_clus_num][
                                                          new_cgvdmbb_clus_num] = vdgs

   return reassigned_cgvdmbb_clusters

def get_pr_obj_to_print(clusmem_pdbpath, vdg_scrr_cg_perm,
                                  num_flanking, atomgroup_dict):
   try:
      if clusmem_pdbpath in atomgroup_dict.keys():
         par = atomgroup_dict[clusmem_pdbpath]
      else:
         par = pr.parsePDB(clusmem_pdbpath)
   except:
      raise ValueError(f'Could not parse {clusmem_pdbpath}.')
      return None
   # Add CG to pr obj
   pr_obj = par.select('occupancy > 2.8') # occ >= 3 is CG
   vdg_scrrs, cg_perm = vdg_scrr_cg_perm
   # Iterate over vdms
   for vdm_scrr in vdg_scrrs:
      flank_obj = add_flank_obj_for_cluslevel_flank(par, vdm_scrr, num_flanking,
                                                    clusmem_pdbpath)
      pr_obj += flank_obj

   return pr_obj.copy()

def add_flank_obj_for_cluslevel_flank(par, vdm_scrr, num_flanking, pdbpath):
   # Check +/- flank num to see whether to include it in the pr obj or not.
   # If it doesn't exist or is a chain break, then don't include it.
   s, c, resnum, resname = vdm_scrr
   if resnum < 0:
      res_sel = f'resnum `{resnum}`'
   else:
      res_sel = f'resnum {resnum}'
   if s == '':
      vdm_ca = par.select(f'chain {c} and {res_sel} and name CA')
   else: 
      vdm_ca = par.select(f'segment {s} and chain {c} and {res_sel} and name CA')

   resnums_to_print = [resnum]
   # Start at vdm. Then, walk up num_flanking residues, breaking if chain break.
   # Then, start at vdm again and walk down, breaking if chain break.
   # walk up
   walk_up_resnums = [i for i in range(resnum+1, resnum+num_flanking+1)]
   walk_up_consec_res = determine_flank_resnums(par, walk_up_resnums, s, c, vdm_ca)
   resnums_to_print.extend(walk_up_consec_res) 
   # walk down
   walk_down_resnums = [i for i in range(resnum-1, resnum-num_flanking-1, -1)]
   walk_down_consec_res = determine_flank_resnums(par, walk_down_resnums, s, c, vdm_ca)
   resnums_to_print.extend(walk_down_consec_res)
   assert len(resnums_to_print) == len(set(resnums_to_print))
   resnums_to_print = [f'`{i}`' if i < 0 else str(i) for i in resnums_to_print]
   num_resnums_to_print = len(set(resnums_to_print))
   resnums_to_print = ' '.join(resnums_to_print)
   if s == '':
      selstr = f"chain {c} and resnum {resnums_to_print}"
   else:
      selstr = f"segment {s} and chain {c} and resnum {resnums_to_print}"

   flank_obj = par.select(selstr)
   num_flank_obj_resindices = len(set(flank_obj.getResindices()))
   if num_flank_obj_resindices != num_resnums_to_print:
      print(f'PDB {pdbpath} flank_obj {resname} has {num_flank_obj_resindices} '
            f'resindices, but resnums_to_print is {resnums_to_print}: '
            f'num_resnums_to_print = ({num_resnums_to_print}).')
   return flank_obj

def get_clus_mem_data(ind, all_cg_coords, all_cg_and_vdmbb_coords, all_flankbb_coords,
                      all_pdbpaths, all_scrr_cg_perm, num_flanking, atomgroup_dict):
   # Given the index of a cluster member, return the cg coords, cg+vdmbb coords, 
   # pdbpath, scrr_cg_perm, pr_obj, and output pdb name.
      
   clusmem_cg_coords = all_cg_coords[ind]
   clusmem_cg_vdmbb_coords = all_cg_and_vdmbb_coords[ind]
   clusmem_flankbb_coords = all_flankbb_coords[ind]
   clusmem_pdbpath = all_pdbpaths[ind]
   clusmem_scrr_cg_perm = all_scrr_cg_perm[ind]
   clusmem_pr_obj = get_pr_obj_to_print(clusmem_pdbpath, clusmem_scrr_cg_perm, 
                                        num_flanking, atomgroup_dict)
   if clusmem_pr_obj is None:
      return None
   
   return [clusmem_cg_coords, clusmem_cg_vdmbb_coords, clusmem_flankbb_coords, 
           clusmem_pdbpath, clusmem_scrr_cg_perm, clusmem_pr_obj]

def print_out_first_pdb_of_clus(pdb_outpath, cg_coords, cg_vdmbb_coords, pr_obj, ref):
   # Align this first pdb being output to the 3-atom reference.
   # Return the name of this cluster member, b/c it's the first one being output
   first_pdb_out = pdb_outpath
   mobile, target = cg_coords[:3], ref
   moved_cg_coords, moved_cg_vdmbb_coords, centroid_transf = write_out_first_pdb(
            mobile, target, pr_obj, pdb_outpath, cg_coords, cg_vdmbb_coords)
   _cg_coords = moved_cg_coords
   first_pdb_cg_vdmbb_coords = moved_cg_vdmbb_coords
   return first_pdb_out, first_pdb_cg_vdmbb_coords, centroid_transf

def get_weights(num_cg_atoms, num_bb_atoms, align_cg_weight):
   cg_weights_arr = np.array([align_cg_weight/num_cg_atoms for i in range(num_cg_atoms)])
   align_vdmbb_weight = 1-align_cg_weight
   num_all_atoms = num_cg_atoms + num_bb_atoms
   vdmbb_weights_arr = np.array([align_vdmbb_weight/num_bb_atoms 
                                 for i in range(num_bb_atoms)])
   weights = np.concatenate((cg_weights_arr, vdmbb_weights_arr))
   weights = weights.reshape((num_all_atoms, 1))
   if not np.abs(weights.sum() - 1.) < 1e-8:
      raise ValueError('Weights do not sum to 1.')
   return weights

def add_vdm_obj_for_cluslevel_cgvdmbb(par, vdm_scrr):
   r = vdm_scrr[2] 
   if r < 0:
      res_sel = f'resnum `{r}`'
   else:
      res_sel = f'resnum {r}'

   if vdm_scrr[0]: # has a pdb segment defined
      res_obj = par.select(
         f'segment {vdm_scrr[0]} and chain {vdm_scrr[1]} and {res_sel}'
         )
   else: # no segment
      res_obj = par.select(
         f'chain {vdm_scrr[1]} and {res_sel}')
   assert len(set(res_obj.getResindices())) == 1
   vdm_res_name = list(set(res_obj.getResnames()))
   assert len(vdm_res_name)  == 1
   assert vdm_res_name[0] == vdm_scrr[3]
   return res_obj

def determine_flank_resnums(par, list_resnums, s, c, vdm_ca):
   # Determine which flanking residues to add to pr obj based on chain breaks. 
   resnums_to_include = []
   assert len(vdm_ca) == 1
   vdm_ca = vdm_ca[0]
   prev_ca = vdm_ca
   # Iterate over resnums 
   for r in list_resnums:
      if r < 0:
         res_sel = f'resnum `{r}`'
      else:
         res_sel = f'resnum {r}'
      if s == '':
         sel = par.select(f'chain {c} and {res_sel} and name CA')
      else:
         sel = par.select(f'segment {s} and chain {c} and {res_sel} and name CA')
      # Return if this residue's CA doesn't exist
      if sel is None:
         return resnums_to_include
      # Check if there's a chain break b/n prev_ca and current CA
      assert len(sel) == 1
      curr_ca = sel[0]
      dist = pr.calcDistance(prev_ca, curr_ca)
      if dist > 4.5: # found chain break
         return resnums_to_include
      else: # otherwise, continue walking
         resnums_to_include.append(r)
         prev_ca = curr_ca
   return resnums_to_include

def flatten_flanking_CAs(cgvdmbb_clus_flankingCAs):
   cgvdmbb_clus_flat_flankCAs = []
   for vdg_flankingCAs in cgvdmbb_clus_flankingCAs:
      # `vdg_flankingCAs` is a list of vdm residues, so flatten it
      flat_flanking_CAs = []
      for res in vdg_flankingCAs:
         for CA_coord in res:
            flat_flanking_CAs.append(CA_coord)
      cgvdmbb_clus_flat_flankCAs.append(flat_flanking_CAs)
   return cgvdmbb_clus_flat_flankCAs

def flatten_flanking_seqs(flankingCAs_clus_flankingseqs):
   flattened_flankingseqs_for_vdgs_in_flankingCA_clus = []
   for vdg_flankingseq in flankingCAs_clus_flankingseqs:
      flat_vdg_flankingseq = []
      for vdm_res in vdg_flankingseq:
         flat_vdg_flankingseq += vdm_res
      flattened_flankingseqs_for_vdgs_in_flankingCA_clus.append(
         flat_vdg_flankingseq)
   return flattened_flankingseqs_for_vdgs_in_flankingCA_clus