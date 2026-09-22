import time
import numpy as np
from collections import OrderedDict, namedtuple
from functools import lru_cache
from ligand_vdgs.functions.clus_helpers import calc_seq_similarity
from ligand_vdgs.functions import vdg_struct_utils as struct_utils
from ligand_vdgs.functions.vdg_struct_utils import (get_res_iden, found_chain_break,
    get_AA_and_CA_coords, get_bb_coords, get_bb_o_coords, build_flank_lookup_index)
from ligand_vdgs.functions.utils import _det3, kabsch_ssd
from ligand_vdgs.functions.vdg_fp_utils import (
    INTERNAL_DISTANCE_EPS, full_row_permutations, internal_distance_bound_matrix,
    precompute_internal_distance_descriptors)

_SEQ_CACHE_MAXSIZE  = 100_000

_EXACT_MEDOID_MAX = 1_000
_MEDOID_SUBSAMPLE = 512

_PAIR_BLOCK = 200_000

_MINIMAX_CANDIDATE_CAP = 4_096
_MINIMAX_BLOCK_PAIRS   = 200_000

EPS = 1e-9

_GRAPH_BATCH_ROWS = 500_000

_PIVOT_K, _PIVOT_K_SMALL, _PIVOT_SMALL_N = 32, 8, 8_000
_PIVOT_SAMPLE = 1_500
_PIVOT_CONTIG_COLS = 4

def _batched_min_rmsd(data, ia, ib, full_perms, n_total,
                      batch_rows=_GRAPH_BATCH_ROWS):
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

def assert_perm_group_closed(perm_group):
   """The pivot bound needs `perm_group` to be a group, not an arbitrary set.

   `identify_mol_automorphisms` raises past its cap rather than truncating, so this
   should never fire; if it ever did, the bound would be inadmissible and edges
   would vanish silently instead of failing.
   """
   if not {(tuple(cg1[cg2].tolist()), tuple(bb1[bb2].tolist()))
           for cg1, bb1 in perm_group for cg2, bb2 in perm_group} <= {
              (tuple(cg.tolist()), tuple(bb.tolist())) for cg, bb in perm_group}:
      raise ValueError('perm_group is not closed under composition, so the pivot '
                       'lower bound would be inadmissible')

def pivot_count(n):
   return _PIVOT_K_SMALL if n < _PIVOT_SMALL_N else _PIVOT_K

def select_pivots(data, n_cg, perm_group, k, cutoff, sample=_PIVOT_SAMPLE, seed=0):
   """Max-min pivot rows, ordered by measured prune rate on the sample (best first)."""
   n = len(data)
   cand = (np.arange(n, dtype=np.intp) if n <= sample else np.sort(
      np.random.default_rng(seed).choice(n, sample, replace=False)).astype(np.intp))
   full_perms = np.stack(full_row_permutations(perm_group, n_cg))
   to_cand = lambda p: _batched_min_rmsd(
      data, np.full(cand.size, int(p), dtype=np.intp), cand, full_perms, data.shape[1])
   chosen, cols = [int(cand[0])], [to_cand(cand[0])]
   best = cols[0].copy()
   while len(chosen) < min(k, cand.size) and float(best.max()) > 0.0:
      chosen.append(int(cand[int(np.argmax(best))]))
      cols.append(to_cand(chosen[-1]))
      np.minimum(best, cols[-1], out=best)
   ii, jj = np.triu_indices(cand.size, 1)
   return np.asarray(chosen, dtype=np.int64)[np.argsort(
      [np.count_nonzero(np.abs(c[ii] - c[jj]) <= cutoff) for c in cols], kind='stable')]

def pivot_distances(data, pivot_ids, n_cg, perm_group, row_range=None, out=None):
   """Exact group-minimised RMSD from every row to each pivot, shaped (k, rows).

   Transposed so `stage1_edges` scans each pivot column contiguously. Uses
   `_batched_min_rmsd`, never `_masked_min_rmsd`: the latter minimises over a pruned
   subset of group elements, which would make the derived bound inadmissible.
   """
   start, stop = (0, len(data)) if row_range is None else row_range
   rows = np.arange(start, stop, dtype=np.intp)
   out = np.empty((len(pivot_ids), rows.size), dtype=np.float32) if out is None else out
   if out.shape != (len(pivot_ids), rows.size):
      raise ValueError(f'pivot_distances: out has shape {out.shape}, expected '
                       f'{(len(pivot_ids), rows.size)}; a partial fill would leave '
                       'zeros that silently disable the prefilter for those rows')
   out[:] = np.stack([
      _batched_min_rmsd(data, np.full(rows.size, int(p), dtype=np.intp), rows,
                        np.stack(full_row_permutations(perm_group, n_cg)), data.shape[1])
      for p in pivot_ids])
   return out

def build_pivot_embedding(data, n_cg, perm_group, cutoff, k=None):
   """Admissible pivot lower-bound embedding for one bucket, shaped (k, n) float32.

   `d(x,y) = min_g sqrt(kabsch_ssd(x, g.y)/n_total)` is a quotient metric (SE(3) and
   `perm_group` both act by isometries), so `d(x,y) >= max_j abs(d(x,p_j) - d(y,p_j))`
   and a pair rejected on that bound cannot be an edge.
   """
   assert_perm_group_closed(perm_group)
   k = pivot_count(len(data)) if k is None else k
   return pivot_distances(data, select_pivots(data, n_cg, perm_group, k, cutoff),
                          n_cg, perm_group)

def stage1_edges(data, threshold, n_cg, perm_group, row_range=None,
                 counters=None, pivots=None):
   n, n_total = len(data), data.shape[1]
   n_res = (n_total - n_cg) // 3
   full_perms = np.stack(full_row_permutations(perm_group, n_cg))
   n_perm = len(full_perms)
   start, stop = (0, n - 1) if row_range is None else row_range
   stop = min(stop, n - 1)
   if stop <= start:
      return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)

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

   desc = precompute_internal_distance_descriptors(data, n_cg)
   bb_arr = data[:, n_cg:]
   cutoff = threshold + INTERNAL_DISTANCE_EPS
   piv = np.ascontiguousarray(
      build_pivot_embedding(data, n_cg, perm_group, cutoff) if pivots is None else pivots,
      dtype=np.float32)
   n_contig = min(_PIVOT_CONTIG_COLS, piv.shape[0])
   dbuf, mbuf, tbuf = (np.empty(n, dtype=np.float32), np.empty(n, dtype=bool),
                       np.empty(n, dtype=bool))

   stats = {'pivot': 0, 'internal_lb': 0, 'bb_lb': 0, 'exact': 0, 'exact_rows': 0,
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
      lo, m = i + 1, n - i - 1
      d, keep_m = dbuf[:m], mbuf[:m]
      np.subtract(piv[0, lo:], piv[0, i], out=d)
      np.less_equal(np.abs(d, out=d), cutoff, out=keep_m)
      for j in range(1, n_contig):
         np.subtract(piv[j, lo:], piv[j, i], out=d)
         np.less_equal(np.abs(d, out=d), cutoff, out=tbuf[:m])
         np.logical_and(keep_m, tbuf[:m], out=keep_m)
      cand = np.flatnonzero(keep_m) + lo
      for j in range(n_contig, piv.shape[0]):
         if cand.size == 0:
            break
         cand = cand[np.abs(piv[j, cand] - piv[j, i]) <= cutoff]
      stats['pivot'] += cand.size

      if cand.size and desc:
         elem_ok = internal_distance_bound_matrix(
            i, cand, desc, n_total, perm_group) <= cutoff
         alive = elem_ok.any(axis=1)
         cand, elem_ok = cand[alive], elem_ok[alive]
      else:
         elem_ok = np.ones((cand.size, n_perm), dtype=bool)
      stats['internal_lb'] += cand.size

      if cand.size and n_res == 2:
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

_CSR_CHUNK_EDGES = 1 << 20

def _scatter_by_row(indices, cursor, rows, values):
   order = np.argsort(rows, kind='stable')
   rows_sorted = rows[order]
   counts = np.bincount(rows_sorted, minlength=cursor.size)
   starts = np.zeros(counts.size + 1, dtype=np.int64)
   np.cumsum(counts, out=starts[1:])
   rank = np.arange(rows_sorted.size, dtype=np.int64) - starts[rows_sorted]
   indices[cursor[rows_sorted] + rank] = values[order]
   cursor += counts

def neighbor_csr_from_chunks(n, load_chunks):
   deg_rev = np.zeros(n, dtype=np.int64)
   deg_fwd = np.zeros(n, dtype=np.int64)
   total = 0
   for qi, qj in load_chunks():
      deg_rev += np.bincount(np.asarray(qj), minlength=n)
      deg_fwd += np.bincount(np.asarray(qi), minlength=n)
      total += len(qi)

   indptr = np.zeros(n + 1, dtype=np.int64)
   np.cumsum(deg_rev + deg_fwd, out=indptr[1:])
   indices = np.empty(2 * total, dtype=np.int32)
   rev_at = indptr[:-1].copy()
   fwd_at = indptr[:-1] + deg_rev
   for qi, qj in load_chunks():
      qi = np.asarray(qi, dtype=np.int32)
      qj = np.asarray(qj, dtype=np.int32)
      _scatter_by_row(indices, rev_at, qj, qi)
      _scatter_by_row(indices, fwd_at, qi, qj)
   return indptr.astype(np.int32, copy=False), indices

def neighbor_csr(n, qi, qj, chunk=_CSR_CHUNK_EDGES):
   qi = np.asarray(qi, dtype=np.int32)
   qj = np.asarray(qj, dtype=np.int32)

   def _chunks():
      return ((qi[a:a + chunk], qj[a:a + chunk])
              for a in range(0, qi.size, chunk))

   return neighbor_csr_from_chunks(n, _chunks)

def butina_partition(indptr, indices):
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

def batched_pair_row_sums(members, n_orders, flankbb_arr, seq_terms):
   k = len(members)
   idx = np.asarray(members, dtype=np.intp)
   ii, jj = np.triu_indices(k, 1)
   row_sums = np.zeros(k, dtype=np.float64)
   for start in range(0, ii.size, _PAIR_BLOCK):
      sl = slice(start, start + _PAIR_BLOCK)
      x, y = idx[ii[sl]], idx[jj[sl]]
      a, b = np.minimum(x, y), np.maximum(x, y)
      best = None
      for p in range(n_orders):
         d = np.zeros(a.size, dtype=np.float64) + seq_terms(a, b, p)
         if flankbb_arr is not None:
            d += _rmsd_rows(flankbb_arr[p][a], flankbb_arr[0][b])
         best = d if best is None else np.minimum(best, d)
      np.add.at(row_sums, ii[sl], best)
      np.add.at(row_sums, jj[sl], best)
   return row_sums

def pose_minimax_prototype(cgvdmbb_data, member_indices, n_cg_atoms,
                        perm_group=None):
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

   if m > _MINIMAX_CANDIDATE_CAP:
      cand_pos = np.linspace(0, m - 1, _MINIMAX_CANDIDATE_CAP).astype(np.intp)
      cand_pos = np.unique(cand_pos)
   else:
      cand_pos = np.arange(m, dtype=np.intp)

   worst = np.empty(cand_pos.size, dtype=np.float64)
   rows_per_block = max(1, _MINIMAX_BLOCK_PAIRS // m)
   for start in range(0, cand_pos.size, rows_per_block):
      stop = min(start + rows_per_block, cand_pos.size)
      rows = cand_pos[start:stop]
      ia = np.repeat(idx[rows], m)
      ib = np.tile(idx, rows.size)
      d = _batched_min_rmsd(data, ia, ib, full_perms, n_total).reshape(rows.size, m)
      worst[start:stop] = d.max(axis=1)
   best = int(np.argmin(worst))
   return members[int(cand_pos[best])], float(worst[best])

def get_leader_clusters(
   data_to_clus, threshold,
   seq_weight=0.5,
   missing_seq_similarity=0.0,
   refresh_medoid_every=64,
   small_refresh_max=16,
   final_exact_medoid_pass=True,
   final_reassign_once=True,
   slot_orders=None,
   _batched=True,
):
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
      return (i, j) if i <= j else (j, i)

   if not hasattr(get_leader_clusters, "_seqsim_cached"):
      get_leader_clusters._SEQ_DATA      = None
      get_leader_clusters._FLANKBB_DATA  = None

      @lru_cache(maxsize=_SEQ_CACHE_MAXSIZE)
      def _seqsim_cached(i, j, missing_similarity, p):
          seq = get_leader_clusters._SEQ_DATA_PERMS[p]
          return calc_seq_similarity(
             seq[i], get_leader_clusters._SEQ_DATA[j],
             missing_similarity=missing_similarity * 100.0)

      get_leader_clusters._seqsim_cached = _seqsim_cached

   get_leader_clusters._seqsim_cached.cache_clear()

   get_leader_clusters._SEQ_DATA = metric_to_data.get('flankseq')
   get_leader_clusters._FLANKBB_DATA = metric_to_data.get('flankbb')

   has_seq = get_leader_clusters._SEQ_DATA is not None
   has_flankbb = get_leader_clusters._FLANKBB_DATA is not None

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
   flankbb_arr = ([np.asarray(d, dtype=np.float32)
                   for d in get_leader_clusters._FLANKBB_DATA_PERMS]
                  if has_flankbb else None)

   _seqsim_cached = get_leader_clusters._seqsim_cached

   _primed = {}

   def _prime(i, targets):
      _primed.clear()
      if not (_batched and has_flankbb and targets):
         return
      t = np.asarray(targets, dtype=np.intp)
      a, b = np.minimum(t, i), np.maximum(t, i)
      for q in range(n_orders):
         for key, val in zip(zip(a.tolist(), b.tolist(), [q] * t.size),
                             _rmsd_rows(flankbb_arr[q][a], flankbb_arr[0][b]).tolist()):
            _primed[key] = val

   def _rmsd_lookup(a, b, p):
      val = _primed.get((a, b, p))
      if val is None:
         val = _rmsd_pair(get_leader_clusters._FLANKBB_DATA_PERMS[p][a],
                          get_leader_clusters._FLANKBB_DATA[b])
      return val

   if has_seq and not 0.0 <= missing_seq_similarity <= 1.0:
      raise ValueError("missing_seq_similarity must be between 0 and 1")

   def _dist_one(a, b, p, early_stop, cutoff):
      total = 0.0
      if has_seq:
         sim = _seqsim_cached(a, b, float(missing_seq_similarity), p)
         total += ((100.0 - sim) / 100.0) * seq_weight
         if early_stop and total > cutoff:
            return total
      if has_flankbb:
         total += _rmsd_lookup(a, b, p)
      return total

   def _dist_idx_idx(i, j, early_stop=True, cap=None):
      cutoff = threshold if cap is None else min(cap, threshold)
      a, b = _ord_pair(i, j)
      best = _dist_one(a, b, 0, early_stop, cutoff)
      for p in range(1, n_orders):
         if best <= 0.0:
            break
         best = min(best, _dist_one(a, b, p, early_stop, min(cutoff, best)))
      return best

   def _row_sums(members):
      if _batched:
         return batched_pair_row_sums(members, n_orders, flankbb_arr, _seq_terms)
      k = len(members)
      sums = np.zeros(k, dtype=np.float64)
      for ii in range(k):
         for jj in range(ii + 1, k):
            d = _dist_idx_idx(members[ii], members[jj], early_stop=False)
            sums[ii] += d
            sums[jj] += d
      return sums

   def _seq_terms(a, b, p):
      if not has_seq:
         return 0.0
      scale = seq_weight / 100.0
      missing = float(missing_seq_similarity)
      return np.fromiter(
         ((100.0 - _seqsim_cached(u, v, missing, p)) * scale
          for u, v in zip(a.tolist(), b.tolist())),
         dtype=np.float64, count=a.size)

   def _exact_medoid(members):
      m = len(members)
      if m <= 2:
         return members[0]
      if m > _EXACT_MEDOID_MAX:
         return _sampled_medoid(members)
      return members[int(np.argmin(_row_sums(members)))]

   def _sampled_medoid(members):
      m = len(members)
      pos = np.unique(np.linspace(0, m - 1, _MEDOID_SUBSAMPLE).astype(np.intp))
      sample = [members[int(p)] for p in pos]
      return sample[int(np.argmin(_row_sums(sample)))]

   def _should_refresh(size_now):
      if size_now <= small_refresh_max:
         return True
      if not refresh_medoid_every:
         return False
      if size_now <= _EXACT_MEDOID_MAX:
         return size_now % refresh_medoid_every == 0
      return size_now & (size_now - 1) == 0

   reps = [0]
   members = [[0]]
   for i in range(1, n):
      _prime(i, reps)
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

   if final_exact_medoid_pass:
      for j in range(len(reps)):
         reps[j] = _exact_medoid(members[j])

      if final_reassign_once:
         item2clus = {idx: j for j, mem in enumerate(members) for idx in mem}
         moved_any = False
         for i in range(n):
            _prime(i, reps)
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

            if best_j != cur_j and best_d <= threshold and best_d + EPS < d_cur:
               item2clus[i] = best_j
               moved_any = True

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

def masked_kabsch_ssd(X, Y, chunk_size=30000):
   X = np.asarray(X, dtype=np.float32)
   Y = np.asarray(Y, dtype=np.float32)
   if X.shape != Y.shape:
      raise ValueError(f"masked_kabsch_ssd got mismatched shapes: {X.shape} vs {Y.shape}")
   if X.ndim != 3 or X.shape[2] != 3:
      raise ValueError(f"masked_kabsch_ssd expected (M, n, 3), got {X.shape}")

   M = X.shape[0]
   if M == 0:
      return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.int64)

   valid = np.isfinite(X).all(axis=2) & np.isfinite(Y).all(axis=2)
   n_eff = valid.sum(axis=1).astype(np.int64)

   ssd = np.empty(M, dtype=np.float64)
   for start in range(0, M, chunk_size):
      stop = min(start + chunk_size, M)
      keep = valid[start:stop, :, None]
      inv = 1.0 / np.maximum(n_eff[start:stop], 1)[:, None, None]
      out = []
      for arr in (X[start:stop], Y[start:stop]):
         c = np.where(keep, arr, 0.0).astype(np.float64)
         c -= np.add.reduce(c, axis=1)[:, None, :] * inv
         out.append(np.where(keep, c, 0.0))
      Xc, Yc = out
      H = np.matmul(np.transpose(Xc, (0, 2, 1)), Yc)
      sv = np.linalg.svd(H, compute_uv=False)
      d = np.where(_det3(H) < 0.0, -1.0, 1.0)
      norms = (np.add.reduce(np.add.reduce(Xc * Xc, axis=2), axis=1)
               + np.add.reduce(np.add.reduce(Yc * Yc, axis=2), axis=1))
      trace = sv[:, 0] + sv[:, 1] + d * sv[:, 2]
      ssd[start:stop] = np.maximum(norms - 2.0 * trace, 0.0)
   return ssd, n_eff

def _rmsd_rows(X, Y):
   ssd, n_eff = masked_kabsch_ssd(X, Y)
   out = np.full(ssd.shape, np.inf, dtype=np.float64)
   ok = n_eff > 0
   out[ok] = np.sqrt(ssd[ok] / n_eff[ok])
   return out

def _rmsd_pair(X, Y):
   X = np.asarray(X, dtype=np.float32)
   Y = np.asarray(Y, dtype=np.float32)
   if X.shape != Y.shape:
      raise ValueError(f"RMSD pair got mismatched shapes: {X.shape} vs {Y.shape}")
   if X.ndim != 2 or X.shape[1] != 3:
      raise ValueError(f"RMSD pair expected shape (N, 3), got {X.shape}")
   return float(_rmsd_rows(X[None, ...], Y[None, ...])[0])

def _permuted_rmsd_pair(X, Y, cg_index_perms, n_cg):
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
   ssds = kabsch_ssd(Y, X_perms)
   return float(np.sqrt(ssds.min() / n_atoms))

def _mark_flanking_chain_breaks(flanking_seq_dict, num_flanking):
   central_CA = np.asarray(flanking_seq_dict[0][1], dtype=np.float32)
   for direction in (1, -1):
      prev_CA = central_CA
      for flank_num in range(1, num_flanking + 1):
         ind = direction * flank_num
         curr_CA = np.asarray(flanking_seq_dict[ind][1], dtype=np.float32)
         if not np.isfinite(curr_CA).all():
            found_chain_break(flanking_seq_dict, ind + direction)
            break
         if np.linalg.norm(curr_CA - prev_CA) > 4.5:
            found_chain_break(flanking_seq_dict, ind)
            break
         prev_CA = curr_CA
   return flanking_seq_dict

def get_vdm_res_features(prody_obj, pdbpath, num_flanking):
   vdm_residues = prody_obj.select('(occupancy) > 1.5 and (occupancy < 2.5)')
   if vdm_residues is None or len(vdm_residues) == 0:
      return {}
   vdm_resinds = set(vdm_residues.getResindices())
   flank_index = build_flank_lookup_index(prody_obj)
   vdms_dict = {}
   for vdm_resind in vdm_resinds:
      vdm_obj = vdm_residues.select(f'resindex {vdm_resind}')
      bb_coords = get_bb_coords(vdm_obj)
      if bb_coords is None:
         continue

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
   if sel is None or len(sel) == 0:
      return np.zeros((0, 3), dtype=np.float32)
   return np.asarray(sel.getCoords(), dtype=np.float32).reshape(-1, 3)

def _min_dist_to_cg(cg_coords, other_coords):
   diff = cg_coords[:, None, :] - other_coords[None, :, :]
   return float(np.sqrt(np.min(np.sum(diff * diff, axis=2))))

def reorder_vdg_subset(vdg_subset, vdms_dict, cg_coords, prody_obj):
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
      _icodes = vdm_res_obj.getIcodes()
      if _icodes is not None and len(set(_icodes)) > 1:
         raise ValueError(
            f'reorder_vdg_subset: residue {_vdg_seg}:{_vdg_ch}:{_vdg_resnum} has '
            'multiple insertion codes (insertion codes are not part of the vdM key)')

      res_bb, res_sc, res_extra = struct_utils.split_residue_heavy_atoms(
         vdm_res_obj, _vdg_resname)
      if res_bb is None:
         res_sc = _coords_or_empty(vdm_res_obj.select(f'{vdm_sel} and sidechain'))
         res_bb = _coords_or_empty(vdm_res_obj.select(f'{vdm_sel} and backbone'))
         res_extra = np.zeros((0, 3), dtype=np.float32)

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
            min_dist_to_bb = _min_dist_to_cg(cg_coords, res_bb)
            reason = (struct_utils.SLOT_BB_CLOSER
                      if min_dist_to_sc - min_dist_to_bb > 0.3
                      else struct_utils.SLOT_SC)

      is_modified = len(res_extra) > 0
      if is_modified and _min_dist_to_cg(cg_coords, res_extra) <= 4.5:
         label = struct_utils.NONCANONICAL_AA_LABEL
      elif reason == struct_utils.SLOT_SC:
         label = vdmAA
      else:
         label = struct_utils.bb_label_for(_vdg_resname)

      aas_of_vdms_in_order.append(label)
      slot_flags_of_vdms_in_order.append(
         reason | (struct_utils.SLOT_MODIFIED if is_modified else 0))

      flankingseqs = [flanking_seq_dict[i][0] for i in sorted_flank_indices]
      flankingCAs = [flanking_seq_dict[i][1] for i in sorted_flank_indices]
      flankingseqs_of_vdms_in_order.append(flankingseqs)
      flanking_CA_coords_of_vdms_in_order.append(flankingCAs)

   assert len(aas_of_vdms_in_order) == len(bb_coords_of_vdms_in_order)
   assert len(slot_flags_of_vdms_in_order) == len(aas_of_vdms_in_order)
   super_list = [aas_of_vdms_in_order, bb_coords_of_vdms_in_order,
                 flankingseqs_of_vdms_in_order, flanking_CA_coords_of_vdms_in_order,
                 seg_ch_res_of_vdms_in_order, slot_flags_of_vdms_in_order,
                 bb_o_coords_of_vdms_in_order]
   return sort_vdGs_by_AA(super_list)

def sort_vdGs_by_AA(super_list):
   assert all(len(sublist) == len(super_list[0]) for sublist in super_list)
   combined = list(zip(*super_list))
   def _sort_key(entry):
      ca = np.asarray(entry[1][1], dtype=np.float32).reshape(-1)
      if ca.size != 3 or not np.isfinite(ca).all():
         raise ValueError(f'sort_vdGs_by_AA: vdM {entry[0]} has no usable CA coords '
                          'for the duplicate-AA tie-break')
      return (entry[0], float(ca[0]), float(ca[1]), float(ca[2]))
   sorted_combined = sorted(combined, key=_sort_key)
   return [list(sublist) for sublist in zip(*sorted_combined)]

def clear_caches():
   if not hasattr(get_leader_clusters, '_seqsim_cached'):
      return
   get_leader_clusters._seqsim_cached.cache_clear()
   get_leader_clusters._SEQ_DATA = None
   get_leader_clusters._FLANKBB_DATA = None
