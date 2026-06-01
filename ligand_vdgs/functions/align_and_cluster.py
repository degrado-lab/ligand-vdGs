# align_and_cluster.py

import os
from itertools import permutations, product
import numpy as np
import prody as pr
from functools import lru_cache
from clus_helpers import (calc_seq_similarity, get_res_iden, found_chain_break,
                          get_AA_and_CA_coords, get_bb_coords)
from utils import kabsch

_RMSD_CACHE_MAXSIZE = 100_000
_SEQ_CACHE_MAXSIZE  = 50_000

EPS = 1e-9   # strict-improvement epsilon for reassignment; tune to 1e-8 if needed

def get_leader_clusters(
   data_to_clus, threshold, AA_subset, size_subset, vdglib_dir,
   seq_weight=0.5,
   refresh_medoid_every=64,                  # periodic refresh for big clusters
   small_refresh_max=16,                     # robust for small clusters
   final_exact_medoid_pass=True,             # polish small clusters cheaply
   final_reassign_once=True                  # one refinement pass
):
   metrics, datasets = [], []
   for data, metric in data_to_clus:
      metrics.append(metric); datasets.append(data)
   metric_to_data = {m: d for m, d in zip(metrics, datasets)}
   n = len(datasets[0])
   for d in datasets: assert len(d) == n

   if n == 1:
      return {1: [0]}, {1: 0}

   # Implement persistent LRU memoization to reuse distances across leader pass, 
   # exact-medoid, and reassignment without per-call sizing or cache resets.

   # Initialize once per process; keep references on the function object
   if not hasattr(get_leader_clusters, "_seqsim_cached"):
      # small global registries mapping dataset ids -> actual arrays
      get_leader_clusters._SEQ_DATASETS  = {}
      get_leader_clusters._RMSD_DATASETS = {}

      # cache sequence similarity (0..100). key: (id(seq_dataset), i, j)
      @lru_cache(maxsize=_SEQ_CACHE_MAXSIZE)
      def _seqsim_cached(ds_id, i, j):
          seq = get_leader_clusters._SEQ_DATASETS[ds_id]
          return calc_seq_similarity(seq[i], seq[j]) # returns percent identity

      # cache RMSD by metric and dataset identity
      @lru_cache(maxsize=_RMSD_CACHE_MAXSIZE)
      def _rmsd_cached(metric_name, ds_id, i, j):
          a, b = _ord_pair(i, j)
          data = get_leader_clusters._RMSD_DATASETS[(metric_name, ds_id)]
          X = data[a]; Y = data[b]
          return _rmsd_pair(X, Y)

      # stash them on the function so they persist across calls
      get_leader_clusters._seqsim_cached = _seqsim_cached
      get_leader_clusters._rmsd_cached   = _rmsd_cached
 
   # aliases
   _seqsim_cached = get_leader_clusters._seqsim_cached
   _rmsd_cached   = get_leader_clusters._rmsd_cached
   _SEQ_DATASETS  = get_leader_clusters._SEQ_DATASETS
   _RMSD_DATASETS = get_leader_clusters._RMSD_DATASETS

   # register the current datasets by identity for caching 
   dsid_seq     = id(metric_to_data['flankseq']) if 'flankseq' in metric_to_data else None
   dsid_flankbb = id(metric_to_data['flankbb'])  if 'flankbb'  in metric_to_data else None
   dsid_cgvdmbb = id(metric_to_data['cgvdmbb'])  if 'cgvdmbb'  in metric_to_data else None

   if dsid_seq is not None and dsid_seq not in _SEQ_DATASETS:
      _SEQ_DATASETS[dsid_seq] = metric_to_data['flankseq']
   if dsid_flankbb is not None and ('flankbb', dsid_flankbb) not in _RMSD_DATASETS:
      _RMSD_DATASETS[('flankbb', dsid_flankbb)] = metric_to_data['flankbb']
   if dsid_cgvdmbb is not None and ('cgvdmbb', dsid_cgvdmbb) not in _RMSD_DATASETS:
      _RMSD_DATASETS[('cgvdmbb', dsid_cgvdmbb)] = metric_to_data['cgvdmbb']

   def _ord_pair(i, j):
      # Order indices so (i, j) and (j, i) share the same cache entry.
      return (i, j) if i <= j else (j, i)

   # ---- distance = seq_dissim * seq_weight + RMSD(flankbb?) + RMSD(cgvdmbb?)
   def _dist_idx_idx(i, j, early_stop=True, cap=None):
      # cap: if provided, we early-exit once total > min(cap, threshold)
      cutoff = threshold if cap is None else min(cap, threshold)
      total = 0.0
      a, b = _ord_pair(i, j)

      if dsid_seq is not None:
         # pull seq sim from persistent cache instead of recomputing
         sim = _seqsim_cached(dsid_seq, a, b)  # 0..100
         total += ((100.0 - sim) / 100.0) * seq_weight
         if early_stop and total > cutoff:
            return total

      if dsid_flankbb is not None:
         # cached RMSD for 'flankbb'
         total += _rmsd_cached('flankbb', dsid_flankbb, a, b)
         if early_stop and total > cutoff:
            return total

      if dsid_cgvdmbb is not None:
         # cached RMSD for 'cgvdmbb'
         total += _rmsd_cached('cgvdmbb', dsid_cgvdmbb, a, b)

      return total

   # ---- exact medoid (O(m^2), but m is small for small clusters) ----
   def _exact_medoid(members):
      if len(members) <= 2:
          return members[0]
      best, best_sum = members[0], float('inf')
      for c in members:
         s = 0.0
         for o in members:
            if o == c: continue
            # same distance, hits persistent caches
            s += _dist_idx_idx(c, o, early_stop=False)  # <- full distance
            if s >= best_sum:  # early break on the SUM, not on threshold
               break
         if s < best_sum:
             best_sum, best = s, c
      return best

   # ---- Leader pass ----
   reps = [0]            # representative indices
   members = [[0]]       # cluster memberships (lists of indices)

   for i in range(1, n):
      best_j, best_d = -1, float('inf')
      for j, r in enumerate(reps):
         # Use persistent caches for distance calls
         d = _dist_idx_idx(i, r, early_stop=True, cap=best_d)
         if d < best_d:
            best_d, best_j = d, j

      if best_d <= threshold and best_j >= 0:
         # assign to existing cluster
         members[best_j].append(i)

         # small-cluster refresh
         size_now = len(members[best_j])
         # robust for tiny clusters
         if size_now <= small_refresh_max:
            reps[best_j] = _exact_medoid(members[best_j])
         # periodic refresh for big clusters
         elif refresh_medoid_every and (size_now % refresh_medoid_every == 0):
            reps[best_j] = _exact_medoid(members[best_j])

      else:
         # start new cluster
         reps.append(i)
         members.append([i])

   # ---- Final polish: exact medoid for every cluster ----
   if final_exact_medoid_pass:
      for j in range(len(reps)):
         reps[j] = _exact_medoid(members[j])

      if final_reassign_once:
         # map each item -> its current cluster
         item2clus = {}
         for j, mem in enumerate(members):
            for idx in mem:
               item2clus[idx] = j

         # single reassignment step: move only if STRICTLY closer to another medoid and within threshold
         moved_any = False
         for i in range(n):
            cur_j = item2clus[i]
            # Distances to medoids hit persistent caches
            d_cur = _dist_idx_idx(i, reps[cur_j], early_stop=True, cap=threshold)
            best_j, best_d = cur_j, d_cur

            # If nothing can beat ~0, skip quickly
            if best_d <= EPS:
               continue

            for j, r in enumerate(reps):
               if j == cur_j:
                  continue
               strict_cap = min(best_d - EPS, threshold)  # if best_d==d_cur initially, this is d_cur-EPS
               if strict_cap <= 0.0:
                  continue  # can't possibly beat current
               d = _dist_idx_idx(i, r, early_stop=True, cap=strict_cap)
               if d < best_d:
                  best_d, best_j = d, j

            # Reassign only if strictly closer and within threshold
            if best_j != cur_j and best_d <= threshold and best_d + EPS < d_cur:
               # move i
               members[cur_j].remove(i)
               members[best_j].append(i)
               item2clus[i] = best_j
               moved_any = True

         # drop any empties and recompute medoids once more
         if moved_any:
            new_reps, new_members = [], []
            for j, mem in enumerate(members):
               if len(mem) == 0: continue
               new_members.append(mem)
               new_reps.append(_exact_medoid(mem))
            members, reps = new_members, new_reps

   clus_assignments = {cnum+1: mem for cnum, mem in enumerate(members)}
   centroids        = {cnum+1: reps[cnum] for cnum in range(len(reps))}
   return clus_assignments, centroids

def _rmsd_pair(X, Y,):
   '''RMSD between two coordinate arrays (N,3); NaN rows (missing flanking residues) are skipped.'''
   X = np.asarray(X, dtype=np.float64)
   Y = np.asarray(Y, dtype=np.float64)
   if X.shape != Y.shape:
      raise ValueError(f"RMSD pair got mismatched shapes: {X.shape} vs {Y.shape}")
   # Drop positions where either array has NaN (chain breaks / missing residues).
   # Mirrors calc_seq_similarity, which already drops 'X' flanking positions.
   valid = ~(np.isnan(X).any(axis=1) | np.isnan(Y).any(axis=1))
   if not valid.any():
      return float('inf')
   X = X[valid]
   Y = Y[valid]
   _, _, ssd = kabsch(X[None, ...], Y[None, ...], chunk_size=1)
   n_atoms = X.shape[0]
   return float(np.sqrt(ssd[0] / n_atoms))

def get_vdm_res_features(prody_obj, pdbpath, num_flanking):
   # Identify the vdM residues (occ == 2). To be safe, select > 1.5 and < 2.5.
   vdm_residues = prody_obj.select('(occupancy) > 1.5 and (occupancy < 2.5)')
   vdm_resinds = set(vdm_residues.getResindices())
   # Record features of the vdm residues (bb coords, flanking residues, pdb paths, etc.)
   vdms_dict = {}
   for vdm_resind in vdm_resinds:
      vdm_obj = vdm_residues.select(f'resindex {vdm_resind}')
      # BB coords
      bb_coords = get_bb_coords(vdm_obj)
      if bb_coords is None:  # missing N, CA, and/or C
         continue

      # Add heavy atoms for this vdM residue
      vdm_heavy_sel = vdm_obj.select('not element H D')
      if vdm_heavy_sel is None or len(vdm_heavy_sel) == 0:
         continue
      vdm_heavy_coords = vdm_heavy_sel.getCoords().astype(np.float32)
      vdm_heavy_names = vdm_heavy_sel.getNames()
      vdm_heavy_elems = vdm_heavy_sel.getElements()

      # Store heavy-atom metadata as a dict so it can be serialized cleanly
      vdm_heavy = {
         "coords": vdm_heavy_coords,
         "names": vdm_heavy_names,
         "elements": vdm_heavy_elems,}

      # Sequence of (contiguous) flanking residues
      flanking_seq_dict = {}  # key = relative flank num (-1, +1, etc.),
                              # value = list(AA identity, CA coords)
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
      flanking_seq_dict[0] = ['vdm', bb_coords[1]]  # label as 'vdm' for easy exclusion
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
            if np.any(np.isnan(curr_CA)):  # curr_CA is None:
               # If chain break, overwrite the residues preceding (if N-term) or
               # succeeding (if C-term) the chain break as AA "X".
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
      vdm_descript = [
          vdm_seg_chain_resnum_resname,
          bb_coords,
          flanking_seq_dict,
          vdm_heavy,]
      vdm_AA = vdm_seg_chain_resnum_resname[-1]
      assert vdm_resind not in list(vdms_dict.keys())
      vdms_dict[vdm_resind] = [vdm_AA, vdm_descript]
   return vdms_dict

def reorder_vdg_subset(vdg_subset, vdms_dict, cg_obj, prody_obj):
   ''' Return per-vdG features reordered by alphabetical AA name, after assigning each 
   vdM as 'bb' (backbone contact) or its sidechain AA. 
   '''
   aas_of_vdms_in_order = []
   bb_coords_of_vdms_in_order = []
   flankingseqs_of_vdms_in_order = []
   flanking_CA_coords_of_vdms_in_order = []
   seg_ch_res_of_vdms_in_order = []
   vdm_heavycoords_of_vdms_in_order = []

   for _vdmresind in vdg_subset:
      vdmAA, vdm_features = vdms_dict[_vdmresind]
      vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict, vdm_heavy_coords = vdm_features
      seg_ch_res_of_vdms_in_order.append(vdm_seg_chain_resnum_resname)
      bb_coords_of_vdms_in_order.append(bb_coords)
      vdm_heavycoords_of_vdms_in_order.append(vdm_heavy_coords)

      flankingseqs = []
      flankingCAs = []
      sorted_flank_indices = sorted(list(flanking_seq_dict.keys()))
      # Check if it's a backbone-only contact. If it is, instead of recording vdmAA, 
      # record "bb" as the vdm resname.
      _vdg_seg, _vdg_ch, _vdg_resnum, _vdg_resname = vdm_seg_chain_resnum_resname

      if _vdg_resnum < 0:
         res_sel = f'resnum `{_vdg_resnum}`'
      else:
         res_sel = f'resnum {_vdg_resnum}'

      if _vdg_seg: 
         vdm_sel = f'segname {_vdg_seg} and chain {_vdg_ch} and {res_sel}'
      else:
         vdm_sel = f'chain {_vdg_ch} and {res_sel}'
      vdm_res_obj = prody_obj.select(f'{vdm_sel} and not element H D')
      vdm_sc = vdm_res_obj.select(f'{vdm_sel} and sidechain') 
      
      if vdm_sc is None:  # then Gly
         aas_of_vdms_in_order.append('bb') # assign gly as 'bb'
      elif len(vdm_sc) == 0: # also Gly
         aas_of_vdms_in_order.append('bb')
      else: # not Gly, so check if it has a sc contact
         has_sidechain_contact = False
         # Determine whether the vdm is a bb or sc contact. H is excluded in cg_obj b/c 
         # a CG H within 4.5A of a sc atom is not necessarily a sc contact and will 
         # give poor results.
         for cg_atom in cg_obj:
            # sc condition
            dists_cg_atom_to_sc = pr.calcDistance(cg_atom, vdm_sc)
            if np.any(dists_cg_atom_to_sc <= 4.4):
               has_sidechain_contact = True 
               break 

         if has_sidechain_contact:
            aas_of_vdms_in_order.append(vdmAA)
         else: 
            # Even if it doesn't have a close sc contact, it may have a weaker sc 
            # contact, so don't rule it out completely. See whether the closest AA atom 
            # is bb or sc.
            vdm_bb = vdm_res_obj.select(f'{vdm_sel} and backbone') 
            min_dist_to_sc = None
            min_dist_to_bb = None
            for cg_atom in cg_obj:
               # get closest dist to AA sc
               dist_to_sc = pr.calcDistance(cg_atom, vdm_sc)
               if min_dist_to_sc is None:
                  min_dist_to_sc = min(dist_to_sc)
               else:
                  min_dist_to_sc = min(min_dist_to_sc, min(dist_to_sc))
               # get closest dist to AA bb
               dist_to_bb = pr.calcDistance(cg_atom, vdm_bb)
               if min_dist_to_bb is None:
                  min_dist_to_bb = min(dist_to_bb)
               else:
                  min_dist_to_bb = min(min_dist_to_bb, min(dist_to_bb))
            
            # Is bb closer to lig by at least 0.3A compared to sc? If yes, then bb.
            if ((min_dist_to_bb < min_dist_to_sc) and 
                (min_dist_to_sc - min_dist_to_bb > 0.3)):
               aas_of_vdms_in_order.append('bb')
            else: # then sc
               aas_of_vdms_in_order.append(vdmAA)

      # Decompress the flanking AA and CA info
      for flank_ind in sorted_flank_indices:
         flank_resname, flank_ca = flanking_seq_dict[flank_ind]
         flankingseqs.append(flank_resname)
         flankingCAs.append(flank_ca)
      flankingseqs_of_vdms_in_order.append(flankingseqs)
      flanking_CA_coords_of_vdms_in_order.append(flankingCAs)
   
   assert len(aas_of_vdms_in_order) == len(bb_coords_of_vdms_in_order)
   # Then, re-order based on alphabetical order
   super_list = [aas_of_vdms_in_order,
                 bb_coords_of_vdms_in_order, 
                 flankingseqs_of_vdms_in_order,
                 flanking_CA_coords_of_vdms_in_order, 
                 seg_ch_res_of_vdms_in_order,
                 vdm_heavycoords_of_vdms_in_order]
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

def get_vdg_AA_and_cg_perms(all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords,
    all_AA_perm_flankingseqs, all_AA_perm_flankingCAs, all_AA_perm_pdbpaths, 
    all_AA_perm_vdm_scrr, all_AA_perm_vdm_heavycoords, all_AA_perm_cg_names, 
    all_AA_perm_cg_elements, all_AA_perm_cg_seg, all_AA_perm_cg_chain,
    all_AA_perm_cg_resnum, all_AA_perm_cg_resname, symmetry_classes,):
    """
    Expand over CG symmetry permutations, carrying vdM features and CG atom metadata.
    Inputs are parallel lists, length = #vdGs-after-AA-perm.
    """

    if symmetry_classes is None:
        return (all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords,
            all_AA_perm_flankingseqs, all_AA_perm_flankingCAs,
            all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr,
            all_AA_perm_vdm_heavycoords,
            all_AA_perm_cg_names, all_AA_perm_cg_elements,
            all_AA_perm_cg_seg, all_AA_perm_cg_chain,
            all_AA_perm_cg_resnum, all_AA_perm_cg_resname,)

    # Recompile vdGs with permuted symmetric CG indices
    all_AA_cg_perm_cg_coords = []
    all_AA_cg_perm_vdm_bbcoords = []
    all_AA_cg_perm_flankingseqs = []
    all_AA_cg_perm_flankingCAs = []
    all_AA_cg_perm_pdbpaths = []
    all_AA_cg_perm_vdm_scrr_cg_perm = []
    all_AA_cg_perm_vdm_heavycoords = []

    all_AA_cg_perm_cg_names = []
    all_AA_cg_perm_cg_elements = []
    all_AA_cg_perm_cg_seg = []
    all_AA_cg_perm_cg_chain = []
    all_AA_cg_perm_cg_resnum = []
    all_AA_cg_perm_cg_resname = []

    for (cgcoords, vdmcoords, seq, CAs, pdbpath, scrr, vdm_heavycoords,
         cg_names, cg_elements, cg_seg, cg_chain, cg_resnum, cg_resname) in zip(
        all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords, all_AA_perm_flankingseqs, 
        all_AA_perm_flankingCAs, all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr,
        all_AA_perm_vdm_heavycoords, all_AA_perm_cg_names, all_AA_perm_cg_elements,
        all_AA_perm_cg_seg, all_AA_perm_cg_chain, all_AA_perm_cg_resnum, 
        all_AA_perm_cg_resname,):
        perms = permute_on_indices(symmetry_classes, cgcoords)

        for perm_ind, perm in enumerate(perms):
            all_AA_cg_perm_cg_coords.append(perm)
            all_AA_cg_perm_vdm_bbcoords.append(vdmcoords)
            all_AA_cg_perm_flankingseqs.append(seq)
            all_AA_cg_perm_flankingCAs.append(CAs)
            all_AA_cg_perm_pdbpaths.append(pdbpath)

            vdm_scrr_perm = [scrr, f"cg_perm_{perm_ind + 1}"]
            all_AA_cg_perm_vdm_scrr_cg_perm.append(vdm_scrr_perm)

            all_AA_cg_perm_vdm_heavycoords.append(vdm_heavycoords)

            # CG metadata duplicated across symmetry-equivalent perms
            all_AA_cg_perm_cg_names.append(cg_names)
            all_AA_cg_perm_cg_elements.append(cg_elements)
            all_AA_cg_perm_cg_seg.append(cg_seg)
            all_AA_cg_perm_cg_chain.append(cg_chain)
            all_AA_cg_perm_cg_resnum.append(cg_resnum)
            all_AA_cg_perm_cg_resname.append(cg_resname)

    return (
        all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords, 
        all_AA_cg_perm_flankingseqs, all_AA_cg_perm_flankingCAs, all_AA_cg_perm_pdbpaths, 
        all_AA_cg_perm_vdm_scrr_cg_perm, all_AA_cg_perm_vdm_heavycoords,
        all_AA_cg_perm_cg_names, all_AA_cg_perm_cg_elements, all_AA_cg_perm_cg_seg, 
        all_AA_cg_perm_cg_chain, all_AA_cg_perm_cg_resnum, all_AA_cg_perm_cg_resname,)

def reassign_cgvdmbb_clusters(cgvdmbb_clus_assignments,
    all_pdbpaths, all_scrr_cg_perm):
   """
   Merge equivalent clusters with degenerate vdGs (diff AA perms and CG perms of the same 
   PDB) and reassign cluster numbers by size. Must merge before deleting duplicates b/c 
   need the duplicate names to determine which clusters are equivalent. 

   Returns: 
       {new_cluster_number: [indices_into_all_*_arrays, ...]}
   where indices correspond to the entries in all_* arrays
   (e.g., all_AA_cg_perm_cg_coords, all_AA_cg_perm_pdbpaths, etc.)
   for mapping between clusters.
   """

   # Dict architecture containing _all_ vdGs with sorted scrrs:
   #   temp_clus_assignments[(clusnum,)] = [
   #       (pdb_base, sorted_scrrs, clus_mem_index, pdbpath), ...]
   temp_clus_assignments = {}

   for cgvdmbb_clusnum, member_indices in cgvdmbb_clus_assignments.items():
      clus_label = f'clus_{cgvdmbb_clusnum}'
      clus_tup = (clus_label,)

      if clus_tup not in temp_clus_assignments:
         temp_clus_assignments[clus_tup] = []

      for ix in member_indices:
         pdbpath = all_pdbpaths[ix]
         pdbname = os.path.basename(pdbpath)

         if pdbname.endswith('.pdb.gz'):
            pdb_base = pdbname[:-len('.pdb.gz')]
         elif pdbname.endswith('.pdb'):
            pdb_base = pdbname[:-len('.pdb')]
         else:
            pdb_base = pdbname

         scrrs, cg_perm = all_scrr_cg_perm[ix]
         # scrrs is a list of [seg, chain, resnum, resname]
         grouped_scrrs = [
            [str(seg), str(chain), str(resnum), str(resname)]
            for (seg, chain, resnum, resname) in scrrs]
         sorted_scrrs = sorted(grouped_scrrs)

         # index for this specific pdb of specific AA perm and cg perm for mapping it 
         # back to the original "all_AA_cg_perm"...etc... lists. 
         clus_mem_index = str(ix) 

         vdg = (pdb_base, sorted_scrrs, clus_mem_index, pdbpath)

         # Remove duplicates
         existing_pairs = [v[:2] for v in temp_clus_assignments[clus_tup]]
         if vdg[:2] not in existing_pairs:
            temp_clus_assignments[clus_tup].append(vdg)

   # Merge_equivalent_clusters -> delete_redun_vdgs -> reassign_temp_clusters
   temp_clus_assignments = merge_equivalent_clusters(temp_clus_assignments)
   pruned_clusters = delete_redun_vdgs(temp_clus_assignments)
   reassigned = reassign_temp_clusters(pruned_clusters)

   # Convert from vdg tuples back to index lists. Reassign clusters.
   renumbered_clusters = {}
   for new_cgvdmbb_clusnum, vdgs in reassigned.items():
      renumbered_clusters[new_cgvdmbb_clusnum] = [int(v[2]) for v in vdgs]

   return renumbered_clusters 

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
   '''
   clus_tuples = list(pruned_clusters.keys())
   original_clus_labels = {}
   for tup in clus_tuples:
      cgvdmbb_clus = tup
      assert cgvdmbb_clus not in original_clus_labels.keys()
      original_clus_labels[cgvdmbb_clus] = pruned_clusters[tup]
   
   # Re-order/renumber cgvdmbb clus nums based on # of vdgs
   reassigned_cgvdmbb_clusters = {}
   sorted_cgvdmbb_clusters = sorted(original_clus_labels.values(), 
                                    key=lambda x: len(x), reverse=True)
   for new_cgvdmbb_clus_num, vdgs in enumerate(sorted_cgvdmbb_clusters, start=1):
      reassigned_cgvdmbb_clusters[new_cgvdmbb_clus_num] = vdgs

   return reassigned_cgvdmbb_clusters