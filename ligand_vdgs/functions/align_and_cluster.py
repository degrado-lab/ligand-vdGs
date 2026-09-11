# align_and_cluster.py

import time
import numpy as np
from collections import OrderedDict, namedtuple
from scipy import sparse
from functools import lru_cache
from ligand_vdgs.functions.clus_helpers import calc_seq_similarity
from ligand_vdgs.functions import vdg_struct_utils as struct_utils
from ligand_vdgs.functions.vdg_struct_utils import (get_res_iden, found_chain_break,
    get_AA_and_CA_coords, get_bb_coords, get_bb_o_coords, build_flank_lookup_index)
from ligand_vdgs.functions.utils import kabsch_ssd
from ligand_vdgs.functions.vdg_fp_utils import (
    INTERNAL_DISTANCE_EPS, fp_tolerances,
    full_row_permutations, internal_distance_bound_matrix,
    precompute_bucket_fingerprints, precompute_internal_distance_descriptors)

_SEQ_CACHE_MAXSIZE  = 100_000   # sequence similarity (Stage 2 sub-clusters only)

# Stage-2 centre selection is O(m^2) in Python-level pair distances, so both the
# cluster size at which it stops being exact and how often it is refreshed have
# to be bounded; see _exact_medoid and get_leader_clusters' leader pass.
_EXACT_MEDOID_MAX = 1_000       # above this, the medoid is picked from a subsample
_MEDOID_SUBSAMPLE = 512

# Sized so a cluster clustered exactly (m <= _EXACT_MEDOID_MAX) has all m(m-1)/2
# of its pairs resident at once. Below that the repeated medoid refreshes of one
# growing cluster evict each other's pairs and every refresh recomputes a fresh
# Kabsch per pair, which turns the leader pass from quadratic-with-reuse into
# cubic. The bound still matters -- an unbounded dict is an OOM on a big bucket --
# so it is tied to the exact-path cap rather than raised freely.
_RMSD_CACHE_MAXSIZE = _EXACT_MEDOID_MAX * (_EXACT_MEDOID_MAX - 1) // 2

# pose_minimax_prototype: candidate cap, and the pair budget per Kabsch block.
_MINIMAX_CANDIDATE_CAP = 4_096
_MINIMAX_BLOCK_PAIRS   = 200_000

EPS = 1e-9   # strict-improvement epsilon for reassignment; tune to 1e-8 if needed

_CacheInfo = namedtuple('_CacheInfo', 'hits misses evictions maxsize currsize')



class _BoundedRmsdCache:
   """Bounded LRU for pair RMSDs, with the membership test ``lru_cache`` lacks.

   Batching only pays if it computes *just* the misses, which means asking what
   is already known -- ``functools.lru_cache`` exposes no such query. The bound
   is kept because it is load-bearing: a bucket of n vdGs has n(n-1)/2 pairs
   (367k for one guanidine bucket alone, far more for a large high-symmetry
   fragment), so an unbounded dict is an out-of-memory waiting to happen.

   Keys are ``(i, j, p)``: an ordered index pair plus the slot ordering applied
   to the lower-indexed record. The stored value is the flank RMSD only, so the
   metric is not part of the key. Indices are relative to the *current* call's
   datasets and collide across calls, so the ``cache_clear`` on every
   ``get_leader_clusters`` entry is what keeps this correct.
   """

   __slots__ = ("_d", "_maxsize", "hits", "misses", "evictions")

   def __init__(self, maxsize):
      self._d = OrderedDict()
      self._maxsize = maxsize
      self.hits = self.misses = self.evictions = 0

   def __contains__(self, key):
      return key in self._d

   def get(self, key):
      d = self._d
      if key in d:
         d.move_to_end(key)
         self.hits += 1
         return d[key]
      self.misses += 1
      return None

   def put(self, key, value):
      d = self._d
      if key in d:
         d.move_to_end(key)
      d[key] = value
      if len(d) > self._maxsize:
         d.popitem(last=False)
         # Nonzero means pairs are being recomputed; the bound needs revisiting.
         self.evictions += 1

   def cache_clear(self):
      self._d.clear()
      self.hits = self.misses = self.evictions = 0

   def cache_info(self):
      return _CacheInfo(self.hits, self.misses, self.evictions,
                        self._maxsize, len(self._d))


# Rows (pair x group element) per batched Kabsch call. Kabsch is ~99% of
# Stage-1 time, and calling it per candidate pair costs ~12x what batching does,
# so the graph builder defers survivors and flushes them in blocks. The cap is on
# rows rather than pairs because a 24-automorphism CG at subset size 2 expands
# each pair up to 48-fold, and the temporaries have to stay bounded regardless.
_GRAPH_BATCH_ROWS = 500_000


def _batched_min_rmsd(data, ia, ib, full_perms, n_total,
                      batch_rows=_GRAPH_BATCH_ROWS):
   """Symmetry-minimised RMSD for pairs (ia[k], ib[k]), in batched Kabsch calls.

   One `kabsch_ssd` call per block of pairs x group elements, rather than one
   per pair: identical values, and the difference between a Stage-1 bucket
   taking seconds and taking minutes.
   """
   n_perm = len(full_perms)
   out = np.empty(ia.size, dtype=np.float64)
   per_block = max(1, batch_rows // n_perm)
   for start in range(0, ia.size, per_block):
      stop = min(start + per_block, ia.size)
      m = stop - start
      A = data[ia[start:stop]]
      B = data[ib[start:stop]]
      X = np.empty((m * n_perm, n_total, 3), dtype=np.float32)
      for k, perm in enumerate(full_perms):
         X[k * m:(k + 1) * m] = A[:, perm]
      Y = np.empty_like(X)
      for k in range(n_perm):
         Y[k * m:(k + 1) * m] = B
      ssd = kabsch_ssd(X, Y).reshape(n_perm, m).min(axis=0)
      out[start:stop] = np.sqrt(ssd / n_total)
   return out


def _masked_min_rmsd(data, ia, ib, elem_mask, full_perms, n_total,
                     batch_rows=_GRAPH_BATCH_ROWS):
   """RMSD for pairs (ia[k], ib[k]) minimised over the group elements their
   `elem_mask[k]` row admits. Every row of `elem_mask` must admit at least one.

   Equal to `_batched_min_rmsd` wherever the pair is within the cutoff, provided
   the elements masked out were rejected by an admissible bound: the element
   realising the minimum then survives. Outside the cutoff the value may be
   larger, which no caller distinguishes.
   """
   pair_idx, perm_idx = np.nonzero(elem_mask)
   full_perms = np.asarray(full_perms, dtype=np.intp)
   ssd = np.empty(pair_idx.size, dtype=np.float64)
   for start in range(0, pair_idx.size, batch_rows):
      stop = min(start + batch_rows, pair_idx.size)
      p = pair_idx[start:stop]
      A = data[ia[p]]
      X = np.take_along_axis(A, full_perms[perm_idx[start:stop]][:, :, None], axis=1)
      ssd[start:stop] = kabsch_ssd(X, data[ib[p]])
   starts = np.zeros(ia.size, dtype=np.intp)
   np.cumsum(elem_mask.sum(axis=1)[:-1], out=starts[1:])
   best = np.minimum.reduceat(ssd, starts) if ssd.size else ssd
   return np.sqrt(best / n_total)


def stage1_edges(data, threshold, n_cg, perm_group, row_range=None,
                 counters=None):
   """Exact symmetry-aware within-cutoff pairs (i, j), i < j, for i in `row_range`.

   Returns `(qi, qj)` as int32 arrays sorted lexicographically, plus per-stage
   counters. Rows are independent, so splitting `range(n - 1)` across calls and
   concatenating the results in row order gives exactly the single-call graph.

   Per row the admissible cascade is: CA/CG-COM fingerprints; the internal-
   distance bound per group element; at subset size 2 the six-atom backbone fit
   per distinct backbone relabeling; then the exact fit, only under the group
   elements the bounds left standing. Every stage sees the same `perm_group` the
   exact RMSD minimises over; screening a subset of the group would make the
   bounds inadmissible, because the true distance is itself a minimum over all
   of it.
   """
   n, n_total = len(data), data.shape[1]
   n_res = (n_total - n_cg) // 3
   full_perms = np.stack(full_row_permutations(perm_group, n_cg))
   n_perm = len(full_perms)
   start, stop = (0, n - 1) if row_range is None else row_range
   stop = min(stop, n - 1)

   # Distinct backbone relabelings in the group, and which one each element
   # uses, for the stages that see backbone atoms but not CG atoms.
   bb_perms, bb_of_elem = [], np.empty(n_perm, dtype=np.intp)
   for e, (_cg_perm, bb_perm) in enumerate(perm_group):
      for k, seen in enumerate(bb_perms):
         if np.array_equal(bb_perm, seen):
            bb_of_elem[e] = k
            break
      else:
         bb_of_elem[e] = len(bb_perms)
         bb_perms.append(np.asarray(bb_perm, dtype=np.intp))
   elems_of_bb = [np.flatnonzero(bb_of_elem == k) for k in range(len(bb_perms))]
   res_orders = [tuple(bb[::3] // 3) for bb in bb_perms]

   fp = precompute_bucket_fingerprints(data, n_cg)
   fp_tol = fp_tolerances(threshold, n_total, n_cg, n_res) if fp else None
   desc = precompute_internal_distance_descriptors(data, n_cg)
   if fp:
      fp_ca = (np.stack([fp['fp0'], fp['fp1']], axis=1) if n_res == 2
               else fp['fp0'][:, None])
      fp_ca_ca = fp.get('fp2')
   bb_arr = data[:, n_cg:]
   cutoff = threshold + INTERNAL_DISTANCE_EPS

   stats = {'fp': 0, 'internal_lb': 0, 'bb_lb': 0, 'exact': 0, 'exact_rows': 0,
            'edges': 0, 'scan_s': 0.0, 'exact_s': 0.0}
   out_i, out_j = [], []
   pend_i, pend_j, pend_mask, pend_rows = [], [], [], 0

   def _flush():
      nonlocal pend_rows
      if not pend_i:
         return
      t0 = time.perf_counter()
      qi = np.concatenate(pend_i)
      qj = np.concatenate(pend_j)
      mask = np.concatenate(pend_mask)
      keep = _masked_min_rmsd(data, qi, qj, mask, full_perms, n_total) <= threshold
      out_i.append(qi[keep].astype(np.int32))
      out_j.append(qj[keep].astype(np.int32))
      stats['edges'] += int(keep.sum())
      stats['exact_s'] += time.perf_counter() - t0
      pend_i.clear(); pend_j.clear(); pend_mask.clear()
      pend_rows = 0

   t_scan = time.perf_counter()
   for i in range(start, stop):
      cand = np.arange(i + 1, n, dtype=np.intp)
      if fp:
         # A slot relabeling swaps which residue each CA fingerprint belongs to,
         # so a pair survives if ANY residue order is within tolerance.
         mask = np.zeros(cand.size, dtype=bool)
         for order in res_orders:
            sub = np.ones(cand.size, dtype=bool)
            for k, src in enumerate(order):
               sub &= np.abs(fp_ca[i, src] - fp_ca[cand, k]) <= fp_tol[k]
            mask |= sub
         if fp_ca_ca is not None:   # CA-CA distance is slot-order invariant
            mask &= np.abs(fp_ca_ca[i] - fp_ca_ca[cand]) <= fp_tol[2]
         cand = cand[mask]
      stats['fp'] += cand.size

      if cand.size and desc:
         elem_ok = internal_distance_bound_matrix(
            i, cand, desc, n_total, perm_group) <= cutoff
         alive = elem_ok.any(axis=1)
         cand, elem_ok = cand[alive], elem_ok[alive]
      else:
         elem_ok = np.ones((cand.size, n_perm), dtype=bool)
      stats['internal_lb'] += cand.size

      if cand.size and n_res == 2:
         # A single N/CA/C triplet is nearly rigid and does not prune; six atoms
         # across two residues do. Rejecting a backbone relabeling rejects every
         # group element that uses it.
         for bb_perm, elems in zip(bb_perms, elems_of_bb):
            need = elem_ok[:, elems].any(axis=1)
            rows = np.flatnonzero(need)
            if rows.size == 0:
               continue
            ssd = kabsch_ssd(bb_arr[i][bb_perm], bb_arr[cand[rows]])
            reject = rows[np.sqrt(ssd / n_total) > cutoff]
            elem_ok[np.ix_(reject, elems)] = False
         alive = elem_ok.any(axis=1)
         cand, elem_ok = cand[alive], elem_ok[alive]
      stats['bb_lb'] += cand.size

      if cand.size:
         stats['exact'] += cand.size
         rows = int(elem_ok.sum())
         stats['exact_rows'] += rows
         pend_i.append(np.full(cand.size, i, dtype=np.intp))
         pend_j.append(cand)
         pend_mask.append(elem_ok)
         pend_rows += rows
         if pend_rows >= _GRAPH_BATCH_ROWS:
            _flush()
   _flush()
   stats['scan_s'] = time.perf_counter() - t_scan - stats['exact_s']
   if counters is not None:
      counters.update(stats)
   if out_i:
      return np.concatenate(out_i), np.concatenate(out_j)
   return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)


def neighbor_csr(n, qi, qj):
   """Symmetric CSR neighbour lists from lexicographically sorted edges (qi < qj).

   Returns int32 ``(indptr, indices)``. Each row's neighbours come out
   ascending: the reverse edges are listed first, then the forward ones, and
   scipy's coo->csr is a counting sort by row that keeps input order within a
   row. That order is load-bearing -- Stage 2 walks a cluster's members in the
   order Butina emits them -- and is pinned in tests/test_butina_clustering.py.
   """
   qi = np.asarray(qi, dtype=np.int32)
   qj = np.asarray(qj, dtype=np.int32)
   rows = np.concatenate([qj, qi])
   cols = np.concatenate([qi, qj])
   graph = sparse.coo_matrix(
      (np.ones(rows.size, dtype=np.int8), (rows, cols)), shape=(n, n)).tocsr()
   return graph.indptr.astype(np.int32), graph.indices.astype(np.int32)


def butina_partition(indptr, indices):
   """Sphere-exclusion partition of a CSR neighbour graph.

   Butina (1999) sphere exclusion -- the same rule as the GROMOS/Daura
   conformational clustering algorithm: repeatedly take the unassigned vertex
   with the most neighbours as a representative and assign it together with
   its still-unassigned neighbours. Returns a list of int32 member arrays, the
   representative first, in the order clusters were formed.

   Every member lies within the cutoff of its representative, and the result
   does not depend on input order (ties on degree go to the lower index), so
   two fragments' cluster counts are comparable. The guarantee is on the
   *radius*, not the diameter: two members may be up to 2*cutoff apart.
   """
   n = indptr.size - 1
   degree = np.diff(indptr)
   order = np.argsort(-degree, kind='stable')
   assigned = np.zeros(n, dtype=bool)
   clusters = []
   for seed in order.tolist():
      if assigned[seed]:
         continue
      nbrs = indices[indptr[seed]:indptr[seed + 1]]
      nbrs = nbrs[~assigned[nbrs]]
      assigned[seed] = True
      assigned[nbrs] = True
      members = np.empty(nbrs.size + 1, dtype=np.int32)
      members[0] = seed
      members[1:] = nbrs
      clusters.append(members)
   return clusters


def get_butina_clusters(cgvdmbb_data, threshold, n_cg_atoms, perm_group,
                        counters=None):
   """Stage-1 pose clustering of one bucket in a single process.

   `stage1_edges` over every row, `neighbor_csr`, then `butina_partition`;
   the driver splits the first of those across processes for large buckets
   and this is the one-shot path for small ones and for tests. `perm_group`
   is the caller's, built once per bucket (`build_perm_group`): every stage
   must minimise over the same group, so no stage builds its own.

   Every Stage-1 atom is mandatory, so non-finite input raises here rather than
   being screened out: the fingerprint prefilter compares with ``<=``, which is
   False against NaN, so such a record would be pruned from every candidate list
   and emerge as a *singleton cluster* -- a plausible-looking result that is
   written to the library and only fails much later, inside `kabsch_ssd` at hit
   finding. Generation already filters these upstream; this makes the invariant
   local to the algorithm that depends on it.
   """
   data = np.asarray(cgvdmbb_data, dtype=np.float32)
   n = len(data)
   if data.size and not np.isfinite(data).all():
      bad = np.flatnonzero(~np.isfinite(data).all(axis=(1, 2)))
      raise ValueError(
         f"get_butina_clusters received non-finite Stage-1 coordinates in "
         f"{bad.size} of {n} records (first at index {int(bad[0])}); every "
         "Stage-1 atom is mandatory")
   if n == 0:
      return []
   if n == 1:
      return [np.zeros(1, dtype=np.int32)]
   qi, qj = stage1_edges(data, threshold, n_cg_atoms, perm_group, counters=counters)
   return butina_partition(*neighbor_csr(n, qi, qj))


def pose_minimax_prototype(cgvdmbb_data, member_indices, n_cg_atoms,
                        perm_group=None):
   """Minimax prototype of a group under pose RMSD, plus its exact pose radius.

   The published minimax prototype (Bien & Tibshirani 2011): the member
   minimising its greatest distance to the rest. Not a medoid -- a medoid
   minimises the *mean*, which does not bound the maximum, so storing one would
   discard the radius guarantee (measured 0.78 A against a 0.50 A cutoff).

   Returns the member minimising its greatest symmetry-aware RMSD to the rest,
   which is the choice that makes the stored radius as small as the group
   allows. Stage 2 partitions a pose cluster by flanking context, and a
   subgroup's members are only guaranteed to be within the cutoff of the
   *Stage-1* seed -- so once the stored row is a subgroup's own member, the
   radius has to be recomputed on pose geometry rather than assumed.
   """
   members = list(member_indices)
   if not members:
      raise ValueError('pose_minimax_prototype needs at least one member')
   if len(members) == 1:
      return members[0], 0.0
   data = np.asarray(cgvdmbb_data, dtype=np.float32)
   n_total = data.shape[1]
   if perm_group is None:
      perm_group = ((np.arange(n_cg_atoms, dtype=np.intp),
                     np.arange(n_total - n_cg_atoms, dtype=np.intp)),)
   full_perms = full_row_permutations(perm_group, n_cg_atoms)
   m = len(members)
   idx = np.asarray(members, dtype=np.intp)

   # Candidates are subsampled past the cap, but each candidate's greatest
   # distance is still measured against *every* member, so the returned radius is
   # exact for whichever row is returned. Only the optimality of the choice
   # weakens: the true minimax member may sit outside the sample, in which case
   # the stored radius is honest but larger than it had to be. Deterministic
   # stride rather than RNG, so a rebuilt library reproduces the same rows.
   if m > _MINIMAX_CANDIDATE_CAP:
      cand_pos = np.linspace(0, m - 1, _MINIMAX_CANDIDATE_CAP).astype(np.intp)
      cand_pos = np.unique(cand_pos)
   else:
      cand_pos = np.arange(m, dtype=np.intp)

   # Accumulated in row blocks. Materializing the m x m matrix is ~10 GB at
   # m=20k, past the per-job memory ceiling with several workers on a node, and
   # only the per-row maximum is ever read from it.
   worst = np.empty(cand_pos.size, dtype=np.float64)
   rows_per_block = max(1, _MINIMAX_BLOCK_PAIRS // m)
   for start in range(0, cand_pos.size, rows_per_block):
      stop = min(start + rows_per_block, cand_pos.size)
      rows = cand_pos[start:stop]
      ia = np.repeat(idx[rows], m)
      ib = np.tile(idx, rows.size)
      d = _batched_min_rmsd(data, ia, ib, full_perms, n_total).reshape(rows.size, m)
      # The self-pair is a zero on this row and cannot be the maximum, so it
      # needs no masking.
      worst[start:stop] = d.max(axis=1)
   best = int(np.argmin(worst))
   return members[int(cand_pos[best])], float(worst[best])


def get_leader_clusters(
   data_to_clus, threshold,
   seq_weight=0.5,
   missing_seq_similarity=0.0,             # fraction used as prior for unknown flanks
   refresh_medoid_every=64,                  # periodic refresh for big clusters
   small_refresh_max=16,                     # robust for small clusters
   final_exact_medoid_pass=True,             # polish small clusters cheaply
   final_reassign_once=True,                 # one refinement pass
   slot_orders=None,                         # interchangeable vdM slot orderings
):
   """Partition one pose cluster by flanking context. Returns the partition only.

   Distance is ``flanking-sequence dissimilarity * seq_weight + flanking-CA
   RMSD``, over whichever of the two datasets is supplied.

   The medoids this computes internally are cluster *centres* -- they decide
   which subgroup each vdG joins -- but they are not returned, because they are
   not what gets stored. The stored row for each subgroup is its
   ``pose_minimax_prototype``: a flank-metric medoid says nothing about pose, and
   the radius recorded next to the row is a pose radius.

   ``slot_orders`` are the vdM slot orderings that permute only same-label slots
   (``vdg_fp_utils.slot_orders`` on the bucket's label list). The distance is
   minimised over them, because Stage 1 already treats two vdGs as the same pose
   if they match under ANY of those orderings: comparing slot 1's flank to slot
   1's flank positionally would then pit the flanks of a swapped pair against
   each other crosswise, inflate the distance, and split one interaction mode
   into two subgroups -- each carrying about half the real cluster_size and
   cluster_num_parents. Over-splitting only, so no stored geometry is wrong, but
   the support counts downstream analysis reads would be. Mixed-label buckets
   (ARG_bb) admit only the identity and cost nothing extra.

   Stage-1 pose geometry does **not** come through here. It is sphere-exclusion
   clustering in :func:`get_butina_clusters`, which needs the whole within-cutoff
   neighbour graph rather than a single pass over current representatives, and
   which -- unlike this pass -- guarantees every member lies within the cutoff of
   the representative that is kept.

   **One call at a time per process.** The datasets and caches are attributes on
   this function object and are keyed by index into *this* call's datasets, so
   they are cleared on entry to every call -- not once per bucket, since a bucket
   makes one call per Stage-1 cluster. The pipeline parallelizes over processes,
   where that state is private.
   """
   metrics, datasets = [], []
   for data, metric in data_to_clus:
      metrics.append(metric); datasets.append(data)
   metric_to_data = {m: d for m, d in zip(metrics, datasets)}
   unknown = set(metrics) - {'flankseq', 'flankbb'}
   if unknown:
      raise ValueError(f'get_leader_clusters handles flankseq/flankbb only, got {sorted(unknown)}')

   if not datasets:
      return {}
   n = len(datasets[0])
   for d in datasets: assert len(d) == n

   if n == 0:
      return {}
   if n == 1:
      return {1: [0]}

   def _ord_pair(i, j):
      # Order indices so (i, j) and (j, i) share the same cache entry.
      return (i, j) if i <= j else (j, i)

   # Persistent memoization across leader pass, exact-medoid, reassignment
   if not hasattr(get_leader_clusters, "_seqsim_cached"):
      get_leader_clusters._SEQ_DATA      = None
      get_leader_clusters._FLANKBB_DATA  = None

      @lru_cache(maxsize=_SEQ_CACHE_MAXSIZE)
      def _seqsim_cached(i, j, missing_similarity, p):
          # p indexes the slot ordering applied to the LOWER-indexed record.
          # Applying the group to one side is enough: it is closed under
          # inverse, so min_g d(g.i, j) == min_g d(i, g.j).
          seq = get_leader_clusters._SEQ_DATA_PERMS[p]
          return calc_seq_similarity(
             seq[i], get_leader_clusters._SEQ_DATA[j],
             missing_similarity=missing_similarity * 100.0)

      _rmsd_store = _BoundedRmsdCache(_RMSD_CACHE_MAXSIZE)

      def _rmsd_cached(a, b, p):
          # Callers pass an already-ordered pair, so (i, j) and (j, i) hit the
          # same entry.
          key = (a, b, p)
          hit = _rmsd_store.get(key)
          if hit is not None:
             return hit
          val = _rmsd_pair(get_leader_clusters._FLANKBB_DATA_PERMS[p][a],
                           get_leader_clusters._FLANKBB_DATA[b])
          _rmsd_store.put(key, val)
          return val
      _rmsd_cached.cache_clear = _rmsd_store.cache_clear
      _rmsd_cached.cache_info = _rmsd_store.cache_info
      _rmsd_cached.store = _rmsd_store
      get_leader_clusters._seqsim_cached = _seqsim_cached
      get_leader_clusters._rmsd_cached   = _rmsd_cached

   # Indices are call-local, so clear on entry to every call.
   get_leader_clusters._seqsim_cached.cache_clear()
   get_leader_clusters._rmsd_cached.cache_clear()

   get_leader_clusters._SEQ_DATA = metric_to_data.get('flankseq')
   get_leader_clusters._FLANKBB_DATA = metric_to_data.get('flankbb')

   has_seq = get_leader_clusters._SEQ_DATA is not None
   has_flankbb = get_leader_clusters._FLANKBB_DATA is not None

   # Both datasets are flattened in slot order with one contiguous, equal-length
   # block per slot, so permuting slots is a block reorder. Materialized once per
   # bucket rather than reordered per pair: a pair distance is evaluated many
   # times across the leader, medoid and reassignment passes.
   orders = tuple(slot_orders) if slot_orders else ((0,),)
   if len(orders) > 1:
      n_slots = len(orders[0])
      def _permuted(dataset):
         out = [dataset]
         for order in orders[1:]:
            per_slot = len(dataset[0]) // n_slots if len(dataset[0]) else 0
            out.append([[rec[s * per_slot + k] for s in order
                         for k in range(per_slot)] for rec in dataset])
         return out
      get_leader_clusters._SEQ_DATA_PERMS = (
         _permuted(get_leader_clusters._SEQ_DATA) if has_seq else None)
      get_leader_clusters._FLANKBB_DATA_PERMS = (
         _permuted(get_leader_clusters._FLANKBB_DATA) if has_flankbb else None)
   else:
      get_leader_clusters._SEQ_DATA_PERMS = [get_leader_clusters._SEQ_DATA]
      get_leader_clusters._FLANKBB_DATA_PERMS = [get_leader_clusters._FLANKBB_DATA]
   n_orders = len(orders)

   _seqsim_cached = get_leader_clusters._seqsim_cached
   _rmsd_cached   = get_leader_clusters._rmsd_cached

   if has_seq and not 0.0 <= missing_seq_similarity <= 1.0:
      raise ValueError("missing_seq_similarity must be between 0 and 1")

   def _dist_one(a, b, p, early_stop, cutoff):
      total = 0.0
      if has_seq:
         sim = _seqsim_cached(a, b, float(missing_seq_similarity), p)  # 0..100
         total += ((100.0 - sim) / 100.0) * seq_weight
         if early_stop and total > cutoff:
            return total
      if has_flankbb:
         total += _rmsd_cached(a, b, p)
      return total

   def _dist_idx_idx(i, j, early_stop=True, cap=None):
      cutoff = threshold if cap is None else min(cap, threshold)
      a, b = _ord_pair(i, j)
      best = _dist_one(a, b, 0, early_stop, cutoff)
      # Identity first, then the rest against the running best: a same-label
      # bucket costs |orders| distance evaluations, a mixed one exactly the
      # original count.
      for p in range(1, n_orders):
         if best <= 0.0:
            break
         best = min(best, _dist_one(a, b, p, early_stop, min(cutoff, best)))
      return best

   def _exact_medoid(members):
      # Full pairwise: no early-exit, warms cache for the reassignment step
      m = len(members)
      if m <= 2:
         # m == 2: the two members have identical row sums, so there is no medoid
         # to find; the first is the deterministic tie-break.
         return members[0]
      if m > _EXACT_MEDOID_MAX:
         return _sampled_medoid(members)
      row_sums = np.zeros(m, dtype=np.float64)
      for ii in range(m):
         for jj in range(ii + 1, m):
            d = _dist_idx_idx(members[ii], members[jj], early_stop=False)
            row_sums[ii] += d
            row_sums[jj] += d
      return members[int(np.argmin(row_sums))]

   def _sampled_medoid(members):
      """Medoid of a deterministic subsample, for clusters too large to do exactly.

      This picks a cluster *centre*, which decides only which subgroup a vdG
      joins during the leader pass -- it is never stored. The row that gets
      written is ``pose_minimax_prototype``'s, chosen on pose geometry and
      re-measured there, so an approximate centre costs assignment quality at the
      margin and nothing in the recorded data.

      Exact is O(m^2) *Python-level* pair distances (measured: n=200 9 s, n=400
      36 s, n=800 158 s), which is cubic overall once the leader pass refreshes a
      growing cluster, and does not finish for the largest bb_bb/ARG buckets.
      """
      m = len(members)
      pos = np.unique(np.linspace(0, m - 1, _MEDOID_SUBSAMPLE).astype(np.intp))
      sample = [members[int(p)] for p in pos]
      k = len(sample)
      row_sums = np.zeros(k, dtype=np.float64)
      for ii in range(k):
         for jj in range(ii + 1, k):
            d = _dist_idx_idx(sample[ii], sample[jj], early_stop=False)
            row_sums[ii] += d
            row_sums[jj] += d
      return sample[int(np.argmin(row_sums))]

   def _should_refresh(size_now):
      """Whether a cluster that just grew to `size_now` should re-pick its centre.

      Every append while small, then a fixed stride, then powers of two. The
      stride alone is O(m / refresh_medoid_every) refreshes of a cluster growing
      to m, each itself quadratic -- cubic overall, and the reason a large
      bb_bb/ARG bucket did not finish. Geometric spacing past the exact-path cap
      makes it O(log m) refreshes, and a centre that is re-picked on every
      doubling is no staler in relative terms than one re-picked every 64.
      """
      if size_now <= small_refresh_max:
         return True
      if not refresh_medoid_every:   # caller disabled periodic refresh entirely
         return False
      if size_now <= _EXACT_MEDOID_MAX:
         return size_now % refresh_medoid_every == 0
      return size_now & (size_now - 1) == 0

   # ---- Leader pass ----
   reps = [0]
   members = [[0]]
   for i in range(1, n):
      best_j, best_d = -1, float('inf')
      for j, r in enumerate(reps):
         d = _dist_idx_idx(i, r, early_stop=True, cap=best_d)
         if d < best_d:
            best_d, best_j = d, j

      if best_d <= threshold and best_j >= 0:
         members[best_j].append(i)
         size_now = len(members[best_j])
         if _should_refresh(size_now):
            reps[best_j] = _exact_medoid(members[best_j])
      else:
         reps.append(i)
         members.append([i])

   # ---- Final polish: exact medoid for every cluster ----
   if final_exact_medoid_pass:
      for j in range(len(reps)):
         reps[j] = _exact_medoid(members[j])

      if final_reassign_once:
         item2clus = {idx: j for j, mem in enumerate(members) for idx in mem}
         moved_any = False
         for i in range(n):
            cur_j = item2clus[i]
            d_cur = _dist_idx_idx(i, reps[cur_j], early_stop=False)
            best_j, best_d = cur_j, d_cur

            if best_d <= EPS:
               continue

            for j, r in enumerate(reps):
               if j == cur_j:
                  continue
               strict_cap = min(best_d - EPS, threshold)
               if strict_cap <= 0.0:
                  continue
               d = _dist_idx_idx(i, r, early_stop=True, cap=strict_cap)
               if d < best_d:
                  best_d, best_j = d, j

            # Reassign only if strictly closer and within threshold
            if best_j != cur_j and best_d <= threshold and best_d + EPS < d_cur:
               item2clus[i] = best_j   # defer list mutation; update dict only
               moved_any = True

         # rebuild members from item2clus in one pass, then recompute medoids
         if moved_any:
            new_members_lists = [[] for _ in reps]
            for i, j in item2clus.items():
               new_members_lists[j].append(i)
            new_reps, new_members = [], []
            for j, mem in enumerate(new_members_lists):
               if not mem:
                  continue
               new_members.append(mem)
               new_reps.append(_exact_medoid(mem))
            members, reps = new_members, new_reps

   return {cnum + 1: mem for cnum, mem in enumerate(members)}

def _rmsd_pair(X, Y):
   """Stage-2 RMSD over coordinate rows that are finite in both structures.

   Stage 1 has no counterpart: its atoms are all mandatory, and records missing
   any of them are dropped upstream by
   `clus_and_deduplicate_vdgs._has_complete_stage1_coords`.
   """
   X = np.asarray(X, dtype=np.float32)
   Y = np.asarray(Y, dtype=np.float32)
   if X.shape != Y.shape:
      raise ValueError(f"RMSD pair got mismatched shapes: {X.shape} vs {Y.shape}")
   if X.ndim != 2 or X.shape[1] != 3:
      raise ValueError(f"RMSD pair expected shape (N, 3), got {X.shape}")
   valid = np.isfinite(X).all(axis=1) & np.isfinite(Y).all(axis=1)
   if not valid.any():
      return float('inf')
   X = X[valid]
   Y = Y[valid]
   ssd = kabsch_ssd(X[None, ...], Y[None, ...], chunk_size=1)
   return float(np.sqrt(ssd[0] / len(X)))

def _permuted_rmsd_pair(X, Y, cg_index_perms, n_cg):
   """Min RMSD over explicit, graph-validated CG automorphisms via Kabsch.

   Test-facing: production Stage-1 code takes the batched path
   (`_batched_min_rmsd`), and this is the readable single-pair reference the
   automorphism and slot-symmetry tests check that path's semantics against.
   """
   X = np.asarray(X, dtype=np.float32)
   Y = np.asarray(Y, dtype=np.float32)
   if X.shape != Y.shape:
      raise ValueError(f"Permuted RMSD got mismatched shapes: {X.shape} vs {Y.shape}")
   bb_X = X[n_cg:]
   n_atoms = len(X)
   P = len(cg_index_perms)
   X_perms = np.empty((P, n_atoms, 3), dtype=np.float32)
   for k, perm_idx in enumerate(cg_index_perms):
      X_perms[k, :n_cg] = X[np.asarray(perm_idx, dtype=np.intp)]
      X_perms[k, n_cg:] = bb_X
   # kabsch_ssd fast path: fixed Y [n_atoms, 3] vs P permuted X [P, n_atoms, 3]
   ssds = kabsch_ssd(Y, X_perms)
   return float(np.sqrt(ssds.min() / n_atoms))


def _mark_flanking_chain_breaks(flanking_seq_dict, num_flanking):
   """Mark discontinuous flanks on each side of the vdM independently.

   The walk ends for two different reasons and they are recorded separately.
   A flank with no usable CA keeps its own FLANK_MISSING label -- the residue is
   unreadable, which is not a statement about the chain -- and only the positions
   *past* it become FLANK_CHAIN_BREAK, since continuity can no longer be checked
   through it. A CA-CA step over 4.5 A is a real break, so that position and
   everything beyond it becomes FLANK_CHAIN_BREAK.

   Mutates and returns the caller's flanking_seq_dict.
   """
   central_CA = np.asarray(flanking_seq_dict[0][1], dtype=np.float32)
   for direction in (1, -1):
      prev_CA = central_CA
      for flank_num in range(1, num_flanking + 1):
         ind = direction * flank_num
         curr_CA = np.asarray(flanking_seq_dict[ind][1], dtype=np.float32)
         if not np.isfinite(curr_CA).all():
            found_chain_break(flanking_seq_dict, ind + direction)
            break
         # prev_CA is finite by construction: it seeds from the vdM's own CA,
         # which get_bb_coords guarantees, and every later value already passed
         # the finiteness check above.
         if np.linalg.norm(curr_CA - prev_CA) > 4.5:
            found_chain_break(flanking_seq_dict, ind)
            break
         prev_CA = curr_CA
   return flanking_seq_dict


def get_vdm_res_features(prody_obj, pdbpath, num_flanking):
   # Identify the vdM residues (occ == 2). To be safe, select > 1.5 and < 2.5.
   vdm_residues = prody_obj.select('(occupancy) > 1.5 and (occupancy < 2.5)')
   if vdm_residues is None or len(vdm_residues) == 0:
      return {}
   vdm_resinds = set(vdm_residues.getResindices())
   # Built once for the whole environment: the flank walk below asks for
   # 2*num_flanking residues per vdM, and each miss-free lookup here saves three
   # ProDy selections (~1 ms). Ambiguous residues are absent from it, so those
   # still take the selection path.
   flank_index = build_flank_lookup_index(prody_obj)
   vdms_dict = {}
   for vdm_resind in vdm_resinds:
      vdm_obj = vdm_residues.select(f'resindex {vdm_resind}')
      bb_coords = get_bb_coords(vdm_obj)
      if bb_coords is None:
         continue

      # Walk flanking residues fwd/bwd, checking for chain breaks; positions past
      # a break are marked FLANK_CHAIN_BREAK, unreadable ones FLANK_MISSING.
      flanking_seq_dict = {}
      for flank_num in range(1, num_flanking + 1):
         for f in [-flank_num, flank_num]:
            AA, CA_coords = get_AA_and_CA_coords(
               prody_obj, vdm_resind + f, flank_index=flank_index)
            flanking_seq_dict[f] = [AA, CA_coords]
      flanking_seq_dict[0] = ['vdm', bb_coords[1]]
      flanking_seq_dict = _mark_flanking_chain_breaks(
         flanking_seq_dict, num_flanking)
      vdm_seg_chain_resnum_resname = get_res_iden(vdm_obj)
      if vdm_seg_chain_resnum_resname is None:
          continue
      vdm_descript = [vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict,
                      get_bb_o_coords(vdm_obj)]
      vdm_AA = vdm_seg_chain_resnum_resname[-1]
      if vdm_resind in vdms_dict:
          raise ValueError(f"Duplicate vdM resindex {vdm_resind} encountered in vdG assembly.")
      vdms_dict[vdm_resind] = [vdm_AA, vdm_descript]
   return vdms_dict

def _coords_or_empty(sel):
   """Heavy-atom coords of a ProDy selection, or an empty (0, 3) array if None."""
   if sel is None or len(sel) == 0:
      return np.zeros((0, 3), dtype=np.float32)
   return np.asarray(sel.getCoords(), dtype=np.float32).reshape(-1, 3)

def _min_dist_to_cg(cg_coords, other_coords):
   """Min distance from any CG atom to any atom in `other_coords`. Both heavy-only."""
   diff = cg_coords[:, None, :] - other_coords[None, :, :]
   return float(np.sqrt(np.min(np.sum(diff * diff, axis=2))))

def reorder_vdg_subset(vdg_subset, vdms_dict, cg_coords, prody_obj):
   # Reorder by alphabetical AA name; label each vdM slot by the moiety that
   # contacts the CG -- its resname, the backbone label BB_LABEL, or 'X' when
   # non-canonical atoms are what touch the CG. Also emits a
   # per-slot flag (SLOT_* in vdg_struct_utils) carrying what the label does not.
   cg_coords = np.asarray(cg_coords, dtype=np.float32)
   aas_of_vdms_in_order = []
   bb_coords_of_vdms_in_order = []
   flankingseqs_of_vdms_in_order = []
   flanking_CA_coords_of_vdms_in_order = []
   seg_ch_res_of_vdms_in_order = []
   slot_flags_of_vdms_in_order = []
   bb_o_coords_of_vdms_in_order = []

   for _vdmresind in vdg_subset:
      vdmAA, vdm_features = vdms_dict[_vdmresind]
      vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict, bb_o = vdm_features
      seg_ch_res_of_vdms_in_order.append(vdm_seg_chain_resnum_resname)
      bb_coords_of_vdms_in_order.append(bb_coords)
      bb_o_coords_of_vdms_in_order.append(bb_o)

      sorted_flank_indices = sorted(list(flanking_seq_dict.keys()))
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
      if vdm_res_obj is None:
         raise ValueError(f'reorder_vdg_subset: no atoms found for "{vdm_sel}"')
      # (seg, chain, resnum) is not unique when insertion codes are in use: the
      # selection pulls every icode variant at once, and this would build one slot
      # out of two residues' atoms. Skip rather than merge -- same policy as
      # dock_utils.get_bsr_combinations on the query side. A lone icode is fine;
      # only a collision is not. The raise is the skip: the caller logs it and
      # moves to the next vdG subset.
      _icodes = vdm_res_obj.getIcodes()   # None when the source set no icodes at all
      if _icodes is not None and len(set(_icodes)) > 1:
         raise ValueError(
            f'reorder_vdg_subset: residue {_vdg_seg}:{_vdg_ch}:{_vdg_resnum} has '
            'multiple insertion codes (insertion codes are not part of the vdM key)')

      # Split by atom name, not by ProDy's flags: the motivating case is a residue
      # that wears a canonical resname while carrying atoms that resname does not
      # have (a GFP chromophore deposited as GLY, an oxidized CYS, an alkylated
      # LYS), and ProDy would hand those extras back as an ordinary sidechain.
      res_bb, res_sc, res_extra = struct_utils.split_residue_heavy_atoms(
         vdm_res_obj, _vdg_resname)
      if res_bb is None:
         # Resname outside the 20, so there is no reference atom set. Fall back to
         # ProDy's own split and treat nothing as non-canonical. Not reachable from
         # the current pipeline (vdG-miner keeps only the 20, and prep renames the
         # rest), kept so a future database change degrades instead of crashing.
         res_sc = _coords_or_empty(vdm_res_obj.select(f'{vdm_sel} and sidechain'))
         res_bb = _coords_or_empty(vdm_res_obj.select(f'{vdm_sel} and backbone'))
         res_extra = np.zeros((0, 3), dtype=np.float32)

      # Which moiety of the *canonical* residue is nearest the CG. Computed for
      # every slot including 'X', where it describes the real residue underneath
      # the modification -- something the label cannot say.
      if len(res_sc) == 0:
         reason = struct_utils.SLOT_NO_SC
      else:
         min_dist_to_sc = _min_dist_to_cg(cg_coords, res_sc)
         if min_dist_to_sc <= 4.5:
            reason = struct_utils.SLOT_SC
         elif len(res_bb) == 0:
            print(f'[WARNING] reorder_vdg_subset: no backbone atoms for "{vdm_sel}"; '
                  'treating as SC contact', flush=True)
            reason = struct_utils.SLOT_SC
         else:
            # Is bb closer to the lig by at least 0.3 A than sc? If yes, then bb.
            min_dist_to_bb = _min_dist_to_cg(cg_coords, res_bb)
            reason = (struct_utils.SLOT_BB_CLOSER
                      if min_dist_to_sc - min_dist_to_bb > 0.3
                      else struct_utils.SLOT_SC)

      is_modified = len(res_extra) > 0
      if is_modified and _min_dist_to_cg(cg_coords, res_extra) <= 4.5:
         # The atoms doing the contacting are not part of the residue this slot is
         # named for, so neither the resname nor a backbone label describes it.
         # Label it 'X' rather than dropping it: the geometry is real and worth
         # keeping, it just must not be counted as an observation of the residue
         # whose name it wears.
         label = struct_utils.NONCANONICAL_AA_LABEL
      elif reason == struct_utils.SLOT_SC:
         label = vdmAA
      else:
         label = struct_utils.bb_label_for(_vdg_resname)

      aas_of_vdms_in_order.append(label)
      slot_flags_of_vdms_in_order.append(
         reason | (struct_utils.SLOT_MODIFIED if is_modified else 0))

      # Decompress flanking info
      flankingseqs = [flanking_seq_dict[i][0] for i in sorted_flank_indices]
      flankingCAs = [flanking_seq_dict[i][1] for i in sorted_flank_indices]
      flankingseqs_of_vdms_in_order.append(flankingseqs)
      flanking_CA_coords_of_vdms_in_order.append(flankingCAs)

   assert len(aas_of_vdms_in_order) == len(bb_coords_of_vdms_in_order)
   assert len(slot_flags_of_vdms_in_order) == len(aas_of_vdms_in_order)
   # Re-order by alphabetical AA name. slot_flags goes last so the sort key
   # (entries 0 and 1) is unaffected.
   super_list = [aas_of_vdms_in_order, bb_coords_of_vdms_in_order,
                 flankingseqs_of_vdms_in_order, flanking_CA_coords_of_vdms_in_order,
                 seg_ch_res_of_vdms_in_order, slot_flags_of_vdms_in_order,
                 bb_o_coords_of_vdms_in_order]
   return sort_vdGs_by_AA(super_list)

def sort_vdGs_by_AA(super_list):
   assert all(len(sublist) == len(super_list[0]) for sublist in super_list)
   combined = list(zip(*super_list))
   def _sort_key(entry):
      # Primary key: AA name; secondary: CA coords for deterministic duplicate-AA ordering
      ca = np.asarray(entry[1][1], dtype=np.float32).reshape(-1)
      if ca.size != 3 or not np.isfinite(ca).all():
         # A vdM CA is mandatory (get_bb_coords rejects the residue otherwise), so
         # this means the caller built the entry itself and the tie-break is gone.
         raise ValueError(f'sort_vdGs_by_AA: vdM {entry[0]} has no usable CA coords '
                          'for the duplicate-AA tie-break')
      return (entry[0], float(ca[0]), float(ca[1]), float(ca[2]))
   sorted_combined = sorted(combined, key=_sort_key)
   return [list(sublist) for sublist in zip(*sorted_combined)]

def clear_caches():
   """Clear Stage-2's memoized distance caches and dataset registries.

   Call between buckets so one bucket's index-keyed cache cannot be read by the
   next. Stage 1 holds no module-level state -- get_butina_clusters is
   self-contained -- so there is nothing to clear on its side."""
   if not hasattr(get_leader_clusters, '_rmsd_cached'):
      return
   get_leader_clusters._rmsd_cached.cache_clear()
   get_leader_clusters._seqsim_cached.cache_clear()
   get_leader_clusters._SEQ_DATA = None
   get_leader_clusters._FLANKBB_DATA = None
