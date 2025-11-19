# align_and_cluster.py

import os
import gzip
from itertools import permutations, combinations, product
import numpy as np
import gc
import prody as pr
import time
from functools import lru_cache

_RMSD_CACHE_MAXSIZE = 100_000
_SEQ_CACHE_MAXSIZE  = 50_000
_ATOMGROUP_CACHE_MAXSIZE = 1000

EPS = 1e-9   # strict-improvement epsilon for reassignment; tune to 1e-8 if needed

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
          return calc_seq_similarity(seq[i], seq[j])

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

def kabsch(X, Y, chunk_size=30000):
   """
   Rotate and translate X into Y to minimize SSD (Kabsch, 1976).
   X, Y: arrays of shape [M, N, 3]
   Returns:
       R:   [M, 3, 3]
       t:   [M, 3]
       ssd: [M]
   """
   R_chunks, t_chunks, ssd_chunks = [], [], []
   M = len(X)

   for start in range(0, M, chunk_size):
      stop = min(start + chunk_size, M)
      Xc_in = X[start:stop]
      Yc_in = Y[start:stop]

      mask = np.logical_or(np.isnan(Xc_in), np.isnan(Yc_in))
      valid_atom = ~np.any(mask, axis=2)
      N_atoms = np.sum(valid_atom, axis=1, keepdims=True)
      zero_rows = (N_atoms == 0).squeeze(-1)
      safeN = np.where(N_atoms == 0, 1, N_atoms)

      X_nonan = np.where(mask, 0.0, Xc_in).astype(np.float32)
      Y_nonan = np.where(mask, 0.0, Yc_in).astype(np.float32)

      valid_coords = np.repeat(valid_atom[:, :, None], 3, axis=2)
      Xbar = (np.sum(X_nonan * valid_coords, axis=1, keepdims=True) /
              safeN[:, None, :].astype(np.float32))
      Ybar = (np.sum(Y_nonan * valid_coords, axis=1, keepdims=True) /
              safeN[:, None, :].astype(np.float32))

      Xc = X_nonan - Xbar
      Yc = Y_nonan - Ybar
      Xc[mask] = 0.0
      Yc[mask] = 0.0

      H = np.matmul(np.transpose(Xc, (0, 2, 1)), Yc)
      U, S, Vt = np.linalg.svd(H, full_matrices=False)
      d = np.sign(np.linalg.det(np.matmul(U, Vt))).astype(np.float32)

      D = np.zeros((H.shape[0], 3, 3), dtype=np.float32)
      D[:, 0, 0] = 1.0
      D[:, 1, 1] = 1.0
      D[:, 2, 2] = d

      R = np.matmul(U.astype(np.float32), np.matmul(D, Vt.astype(np.float32)))
      t = (Ybar - np.matmul(Xbar, R)).reshape(-1, 3)

      XRmY = np.matmul(Xc, R) - Yc
      ssd = np.sum(XRmY ** 2, axis=(1, 2)).astype(np.float64)

      if np.any(zero_rows):
         ssd[zero_rows] = np.inf

      R_chunks.append(R)
      t_chunks.append(t.astype(np.float32))
      ssd_chunks.append(ssd)

   R = np.concatenate(R_chunks) if R_chunks else np.empty((0, 3, 3), np.float32)
   t = np.concatenate(t_chunks) if t_chunks else np.empty((0, 3), np.float32)
   ssd = np.concatenate(ssd_chunks) if ssd_chunks else np.empty((0,), np.float64)
   return R, t, ssd

def _rmsd_pair(X, Y,):
   '''RMSD between two coordinate arrays (N,3)'''
   X = np.asarray(X, dtype=np.float64)
   Y = np.asarray(Y, dtype=np.float64)
   if X.shape != Y.shape:
      raise ValueError(f"RMSD pair got mismatched shapes: {X.shape} vs {Y.shape}")
   # Use vectorized Kabsch with a single pair (M=1)
   _, _, ssd = kabsch(X[None, ...], Y[None, ...], chunk_size=1)
   n_atoms = X.shape[0]
   return float(np.sqrt(ssd[0] / n_atoms))

def calc_seq_similarity(list1, list2):
   # Percent identity of flanking residue names after dropping positions labeled 'vdm' 
   # or 'X'. Returns identity=0 (max dissimilarity) if all positions drop.
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

   # Calculate the percentage of matches. If all residues were dropped out, 
   # it should be treated as max dissimilarity.
   if len(list1) == 0:
      match_percentage = 0
   else:
      match_percentage = (matches / len(list1)) * 100
   return match_percentage

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
      if bb_coords is None:
         continue

      # Sequence of (contiguous) flanking residues
      flanking_seq_dict = {} # key = relative flank num (-1, +1, etc.), 
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
            if np.any(np.isnan(curr_CA)): # curr_CA is None:
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
      flanking_seq_dict[overwrite_ind] = \
         ['X', np.array([np.nan, np.nan, np.nan])] # ['X', None]
   return flanking_seq_dict

def get_AA_and_CA_coords(prody_obj, current_resindex):
   '''
   Return (AA, CA_coords) for a residue index.
   If residue is missing/non-protein, returns AA='X' and CA_coords = [nan, nan, nan].
   '''
   # Build selection for this residue index
   if current_resindex < 0: 
      sel_str = f'resindex `{current_resindex}`'
   else:
      sel_str = f'resindex {current_resindex}'

   curr_resindex_obj = prody_obj.select(sel_str)
   if curr_resindex_obj is None or curr_resindex_obj.protein is None:
      return 'X', np.array([np.nan, np.nan, np.nan])

   # Resolve residue name
   curr_res_AA = list(set(curr_resindex_obj.getResnames()))
   if len(curr_res_AA) != 1:
      print(f'[WARNING] get_AA_and_CA_coords: \nResindex {current_resindex} in '
            f'{prody_obj.getTitle()} contains >1 AA: {curr_res_AA}. Defaulting to AA="X".')
      return 'X', np.array([np.nan, np.nan, np.nan])
   AA = curr_res_AA[0]

   # Select CA(s)
   CA_sel = curr_resindex_obj.select(sel_str + ' and name CA')
   if CA_sel is None or len(CA_sel) == 0:
      return 'X', np.array([np.nan, np.nan, np.nan])

   try:
      ca_atom = _pick_single_ca(CA_sel, prev_ca=None)
      CA_coords = ca_atom.getCoords()
      # Ensure shape is (3,)
      CA_coords = np.asarray(CA_coords, dtype=float).reshape(3,)
   except Exception:
      # print warning and be specific about error
      print(f'[WARNING] Could not resolve CA for resindex {current_resindex} in '
            f'{prody_obj.getTitle()}. Defaulting to AA="X".')
      return 'X', np.array([np.nan, np.nan, np.nan])

   return AA, CA_coords

def get_cg_coords(prody_obj, pdbpath):
   # The CG atoms have their occupancies set to >= 3.0, with unique values (e.g., 
   # 3.0, 3.1, 3.2, etc.) to allow a 1:1 correspondence of equivalent atoms between 
   # different ligands.
   cg = prody_obj.select('occupancy > 2.9')

   if cg is None or len(cg) == 0:
      print(f'[WARNING] get_cg_coords: no atoms with occupancy > 2.9 in {pdbpath}')
      return None

   num_atoms = len(cg)
   cg_coords = []
   for ind in range(num_atoms):
      occ = f'3.{ind}'
      atom = cg.select(f'occupancy == {occ}')
      if atom is None or len(atom) != 1:
         print(f'[WARNING] get_cg_coords: {0 if atom is None else len(atom)} atoms '
               f'are occupancy {occ} in {pdbpath}.')
         return None
      atom_coords = atom.getCoords()[0]
      cg_coords.append(atom_coords)

   cg_coords = np.asarray(cg_coords, dtype=float)

   if not np.isfinite(cg_coords).all(): # reject NaN or inf coords
      return None

   return cg_coords

def get_bb_coords(obj):
   bb_coords = []
   for atom in ['N', 'CA', 'C']:
      atom_obj = obj.select(f'name {atom}')
      if atom_obj is None:
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
                  return None
      coord = atom_obj.getCoords()[0]
      bb_coords.append(coord)
   return bb_coords

def get_vdg_subsets(input_list, target_size):
    # Generate combos of residues containing `target_size` elements.
   if len(input_list) < target_size:
        return []  
   return list(combinations(input_list, target_size)) 

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

def reorder_vdg_subset(vdg_subset, vdms_dict, cg_obj, prody_obj):
   ''' Return per-vdG features reordered by alphabetical AA name, after assigning each 
   vdM as 'bb' (backbone contact) or its sidechain AA. 
   '''
   aas_of_vdms_in_order = []
   bb_coords_of_vdms_in_order = []
   flankingseqs_of_vdms_in_order = []
   flanking_CA_coords_of_vdms_in_order = []
   seg_ch_res_of_vdms_in_order = []
   for _vdmresind in vdg_subset:
      vdmAA, vdm_features = vdms_dict[_vdmresind]
      vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict = vdm_features
      seg_ch_res_of_vdms_in_order.append(vdm_seg_chain_resnum_resname)
      bb_coords_of_vdms_in_order.append(bb_coords)
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
   all_pdbpaths, all_scrr_cg_perm, all_cg_and_vdmbb_coords, all_flankbb_coords, 
   num_flanking, first_pdb_out, first_pdb_cg_vdmbb_coords, weights, atomgroup_dict, 
   print_flankbb, symmetry_classes, reordered_AAs, write_cluster_members, 
   clusterlabel=None):
    
   assert len(all_pdbpaths) == len(all_cg_and_vdmbb_coords)
   ref = np.array([[0, 0, 0], [-1, 0, 1], [1, -1, 0]])  # reference for aligning CG

   failed = []
   for clusnum, clus_mem_indices in sorted(clus_assignments.items()):
      subdir = f'{clusterlabel}clus_{clusnum}' if clusterlabel else f'clus_{clusnum}'
      clusnum_dir = os.path.join(clusdir, subdir)
      os.makedirs(clusnum_dir, exist_ok=False)  # Intentionally fail if it's a re-run
      centroid_ind = centroid_assignments[clusnum]
      # Output the centroid first
      data = get_clus_mem_data(centroid_ind, all_cg_coords, all_cg_and_vdmbb_coords, 
         all_flankbb_coords, all_pdbpaths, all_scrr_cg_perm, num_flanking, 
         atomgroup_dict)
      if data is None:
         failed.append(all_pdbpaths[centroid_ind])
         continue
      (cent_cg_coords, cent_cg_vdmbb_coords, cent_flankbb_coords, cent_pdbpath, 
       cent_scrr_cg_perm, cent_pr_obj) = data
      cent_pdb_outpath = get_clus_pdb_outpath(clusnum, clusnum_dir, cent_pdbpath, 
                                  cent_scrr_cg_perm, centroid_ind, is_centroid=True)

      if first_pdb_out is None:
         # Write the very first PDB aligned to the 3-atom reference
         try:
            first_pdb_out, first_pdb_cg_vdmbb_coords = print_out_first_pdb_of_clus(
                cent_pdb_outpath, cent_cg_coords, cent_cg_vdmbb_coords, cent_pr_obj, 
                ref, cent_scrr_cg_perm, print_flankbb)
            # For the first cluster, the "moved centroid coords" are the first-PDB coords
            moved_cent_coords = first_pdb_cg_vdmbb_coords
         except Exception as e:
            print(f'[ERROR] {reordered_AAs} cluster {clusnum} centroid, the first '
                  f'PDB, failed ({cent_pdbpath}): \n{e}\n'
                  f'Skipping entire AA bucket; must troubleshoot.')
            first_pdb_out, first_pdb_cg_vdmbb_coords = None, None
            continue
      else:
         # Align this cluster's centroid to the very first PDB’s coords (target_coords). 
         # `get_transf_and_coords` handles symmetry permutations.
         try:
            moved_cent_transf, moved_cent_coords = get_transf_and_coords(
                     cent_cg_vdmbb_coords, first_pdb_cg_vdmbb_coords, weights, 
                     cent_pr_obj, cent_scrr_cg_perm, symmetry_classes)
         except Exception as e:
            print(f'[ERROR] {reordered_AAs} cluster {clusnum} centroid failed ({cent_pdbpath}): {e}')
            failed.append(cent_pdbpath)
            raise

         write_out_subsequent_clus_pdbs(cent_pr_obj, cent_pdb_outpath, 
                             cent_scrr_cg_perm, print_flankbb, moved_cent_transf)

      if not write_cluster_members: 
         continue
      # Now, process non-centroid members
      for ind in clus_mem_indices:
         if ind == centroid_ind:
            continue

         try:
            mem_data = get_clus_mem_data(
               ind, all_cg_coords, all_cg_and_vdmbb_coords, all_flankbb_coords,
               all_pdbpaths, all_scrr_cg_perm, num_flanking, atomgroup_dict)
            if mem_data is None:
               failed.append(all_pdbpaths[ind])
               continue

            (clusmem_cg_coords, clusmem_cg_vdmbb_coords, clusmem_flankbb_coords,
             clusmem_pdbpath, clusmem_scrr_cg_perm, clusmem_pr_obj) = mem_data

            clusmem_pdb_outpath = get_clus_pdb_outpath(
               clusnum, clusnum_dir, clusmem_pdbpath, clusmem_scrr_cg_perm, ind, 
               is_centroid=False)

            if clusnum == 1: # the centroid to align onto _is_ the first pdb instead of a 
               # cluster cent that's aligned onto the first pdb.
               target_coords = first_pdb_cg_vdmbb_coords
            else:
               target_coords = moved_cent_coords

            moved_transf, _ = get_transf_and_coords(
               clusmem_cg_vdmbb_coords, target_coords, weights, clusmem_pr_obj,
               clusmem_scrr_cg_perm, symmetry_classes)

            write_out_subsequent_clus_pdbs(
               clusmem_pr_obj, clusmem_pdb_outpath, clusmem_scrr_cg_perm,
               print_flankbb, moved_transf)

         except Exception as e:
            print(f"[WARNING] {reordered_AAs} cluster member {ind} of clusnum "
                  f"{clusnum} ({all_pdbpaths[ind]}): {e}")
            failed.append(all_pdbpaths[ind])

   return first_pdb_out, first_pdb_cg_vdmbb_coords, failed

def get_transf_and_coords(mobile_coords, target_coords, weights, obj, scrrs, 
                          symmetry_classes):
   if weights is None:
      weights = np.array([1/len(mobile_coords) for i in mobile_coords])
   # Superpose each permutation
   cgcoords = mobile_coords[:len(symmetry_classes)]
   perms = permute_on_indices(symmetry_classes, cgcoords)
   lowest_rmsd, best_transf, best_moved_coords = None, None, None 
   # Iterate through each permutation
   for perm in perms:
      # combining cg and vdmbb
      rearranged_coords = np.vstack([perm, mobile_coords[len(symmetry_classes):]])
      assert weights.shape[0] == rearranged_coords.shape[0], \
            "weights/coords length mismatch (rearranged)"
      moved_coords, transf = pr.superpose(rearranged_coords, target_coords, weights.ravel())
      rmsd = pr.calcRMSD(moved_coords, target_coords, weights.ravel())
      if lowest_rmsd is None or rmsd < lowest_rmsd:
         lowest_rmsd = rmsd
         best_transf = transf
         best_moved_coords = moved_coords
   return best_transf, best_moved_coords

def write_out_subsequent_clus_pdbs(pr_obj, pdb_outpath, scrrs, print_flankbb, transf):
   # Align to the first PDB if it's a centroid or to its cluster centroid if it's a cluster mem
   pr_obj_copy = pr_obj.copy() # b/c of mutability
   pr.applyTransformation(transf, pr_obj_copy)
   # Set occupancies so that only the vdMs being clustered are 2.0 (whereas originally, 
   # all probe contacts were 2.0)
   pr_obj_copy = reset_occs(pr_obj_copy, scrrs)
   # Write out PDB
   if not print_flankbb: 
      manually_write_pdb(pdb_outpath, pr_obj_copy.select('occupancy > 1.9')) # selects occ>=2.0
   else:
      manually_write_pdb(pdb_outpath, pr_obj_copy)
   
   del pr_obj_copy  # sometimes, deep copying leaves open handles 
   gc.collect()     

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

def manually_write_pdb(outpath, pr_obj):
   # Workaround because of sge's "Too many open files" error.
   with gzip.open(outpath, 'wt') as f:
      pr.writePDBStream(f, pr_obj)

def write_out_first_pdb(mobile, target, pr_obj, outpath, clusmem_cg_vdmbb_coords, scrrs, 
                        print_flankbb):
   # No weights because it's only aligning 3 atoms of CG
   mobile, target = np.array(mobile), np.array(target)
   moved_3atom_coords, transf = pr.superpose(mobile, target)
   pr_obj_copy = pr_obj.copy() # b/c of mutability
   pr.applyTransformation(transf, pr_obj_copy)
   moved_cg_vdmbb_coords = pr.applyTransformation(transf, clusmem_cg_vdmbb_coords)
   # Set occupancies so that only the vdMs being clustered are 2.0
   pr_obj_copy = reset_occs(pr_obj_copy, scrrs)
   # Write out PDB
   if not print_flankbb: 
      manually_write_pdb(outpath, pr_obj_copy.select('occupancy > 1.9')) # selects occ=2.0
   else:
      manually_write_pdb(outpath, pr_obj_copy)
   
   del pr_obj_copy  # sometimes, deep copying leaves open handles 
   gc.collect()     
   
   return moved_cg_vdmbb_coords

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

def get_pr_obj_to_print(clusmem_pdbpath, vdg_scrr_cg_perm,
                        num_flanking, atomgroup_dict):
    exceptions = []  # collect all exceptions we hit

    for attempt in range(6):  # 1 initial try + 5 retries
        try:
            # Use cache only if enabled AND the key is present
            if _ATOMGROUP_CACHE_MAXSIZE > 0 and clusmem_pdbpath in atomgroup_dict:
                par = atomgroup_dict[clusmem_pdbpath]
            else:
                par = pr.parsePDB(clusmem_pdbpath)
                if _ATOMGROUP_CACHE_MAXSIZE > 0:
                    # Evict oldest (FIFO) if at capacity
                    if len(atomgroup_dict) >= _ATOMGROUP_CACHE_MAXSIZE:
                        try:
                            atomgroup_dict.pop(next(iter(atomgroup_dict)))
                        except StopIteration:
                            pass
                    atomgroup_dict[clusmem_pdbpath] = par

            # if we got here, it worked
            return par

        except Exception as e:
            if attempt < 5:  # still have retries left
                time.sleep(30)
            else:
                print(f"[WARNING] get_pr_obj_to_print failed on {clusmem_pdbpath}")
                return None

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
      print(f'[WARNING] {pdbpath} flank_obj {resname} has {num_flank_obj_resindices} '
            f'resindices, but resnums_to_print is {resnums_to_print}: '
            f'num_resnums_to_print = ({num_resnums_to_print}).')
   return flank_obj

def get_clus_mem_data(ind, all_cg_coords, all_cg_and_vdmbb_coords, all_flankbb_coords,
                      all_pdbpaths, all_scrr_cg_perm, num_flanking, atomgroup_dict):
   # Given the index of a cluster member, return its CG coords, CG+vdM bb coords, 
   # flanking bb coords, pdb path, scrr/cg-perm label, and ProDy AtomGroup.
      
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

def print_out_first_pdb_of_clus(pdb_outpath, cg_coords, cg_vdmbb_coords, pr_obj, ref, 
                                scrrs, print_flankbb):
   # Align this first pdb being output to the 3-atom reference.
   # Return the name of this cluster member, b/c it's the first one being output
   first_pdb_out = pdb_outpath
   if len(cg_coords) >= 3:
      mobile, target = cg_coords[:3], ref
   moved_cg_vdmbb_coords = write_out_first_pdb(
      mobile, target, pr_obj, pdb_outpath, cg_vdmbb_coords, scrrs, print_flankbb)
   first_pdb_cg_vdmbb_coords = moved_cg_vdmbb_coords
   return first_pdb_out, first_pdb_cg_vdmbb_coords

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

def _pick_single_ca(ca_sel, prev_ca):
   # Return a single CA atom object from a ProDy selection, handling altlocs/dups
   if len(ca_sel) == 1:
      return ca_sel[0]

   # If we have multiple CAs (altlocs/dups), prefer the one closest to prev_ca.
   if prev_ca is None:
      return ca_sel[0]

   best_atom, best_dist = None, float('inf')
   for atom in ca_sel:
      d = pr.calcDistance(prev_ca, atom)
      if d < best_dist:
         best_atom, best_dist = atom, d
   return best_atom

def determine_flank_resnums(par, list_resnums, s, c, vdm_ca):
   # Determine which flanking residues to add to pr obj based on chain breaks. 
   resnums_to_include = []
   assert len(vdm_ca) == 1
   vdm_ca = vdm_ca[0]
   prev_ca = vdm_ca
   # Iterate over resnums 
   for r in list_resnums:
      res_sel = f'resnum `{r}`' if r < 0 else f'resnum {r}'
      if s == '':
         ca_sel = par.select(f'chain {c} and {res_sel} and name CA')
      else:
         ca_sel = par.select(f'segment {s} and chain {c} and {res_sel} and name CA')

      # If missing, stop (chain break).
      if ca_sel is None or len(ca_sel) == 0:
         return resnums_to_include

      # Collapse to one CA robustly
      try:
         curr_ca = _pick_single_ca(ca_sel, prev_ca)
      except Exception:
         return resnums_to_include

      # Check chain break distance to previous CA
      dist = pr.calcDistance(prev_ca, curr_ca)
      if dist > 4.5:  # chain break
          return resnums_to_include

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

def reset_occs(pr_obj_copy, scrrs):
   # Set occupancies so that only the vdMs being clustered are 2.0
   pr_obj_copy.protein.setOccupancies(1.0) # reset protein only; don't touch ligand
   scrr_list, perm = scrrs
   for scrr in scrr_list: 
      seg, chain, resnum, resname = scrr
      if resnum < 0:
         res_sel = f'resnum `{resnum}`'
      else:
         res_sel = f'resnum {resnum}'
      if seg == '':
         sel = f'chain {chain} and {res_sel} and resname {resname}'
      else:
         sel = f'segname {seg} and chain {chain} and {res_sel} and resname {resname}'
      
      target = pr_obj_copy.select(sel)
      if target is not None:
          target.setOccupancies(2.0)

   return pr_obj_copy