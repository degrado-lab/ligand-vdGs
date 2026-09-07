"""Admissible internal-distance prefilters for vdG pair comparisons."""
import math

import numpy as np

# The fingerprint scheme is defined for these subset sizes only.
_MIN_SUBSET_SIZE, _MAX_SUBSET_SIZE = 1, 2

# Multiplier applied to the proven bound in fp_tolerances. 1.0 is already
# provably safe -- the bound is tight, and the worst |delta fp| seen over 85k
# real true-match pairs was 0.78x it -- so this is pure insurance against the
# derivation being wrong, not a tuning parameter.
#
# It must stay a *multiplier* rather than a flat tolerance in Å: the bound grows
# with the RMSD threshold and the atom count, so a flat value that looks
# generous for the library's own buckets (bounds 1.4-3.4 Å) goes inadmissible
# for a hit-finding run with --rmsd-threshold above ~0.85, which is exactly the
# way the old FP_TOL = 5.0 was unsafe.
#
# Cost, measured as the fraction of pairs surviving the prefilter (lower is
# better): 1.0x -> 0.60, 1.25x -> 0.71, 1.5x -> 0.79, 2.0x -> 0.90.
FP_SAFETY_FACTOR = 1.5

# Float32 distance descriptors can round upward by a few ulps. The analytic
# lower bound is exact; this padding prevents descriptor quantization from
# turning a boundary match into a false rejection.
INTERNAL_DISTANCE_EPS = 1e-4


def fp_tolerances(rmsd_threshold, n_total, n_cg, n_res):
    """Admissible per-fingerprint tolerances for an `rmsd_threshold` Å cutoff.

    Returns one tolerance per fingerprint, in the order
    `precompute_bucket_fingerprints` emits them: `(fp0,)` for one residue,
    `(fp0, fp1, fp2)` for two.

    Every fingerprint is an *internal* distance, so it is invariant to the
    superposition; only the per-atom deviations matter. If two n-atom structures
    superpose at RMSD <= T then sum(d_i^2) <= n*T^2, and:

      fp0/fp1 = ||CA - CG_com||.  |delta| <= d_CA + ||delta CG_com||, and
        ||delta CG_com|| <= sqrt(S_cg / n_cg) by Cauchy-Schwarz, so maximizing
        under d_CA^2 + S_cg <= n*T^2 gives  T * sqrt(n * (1 + 1/n_cg)).
      fp2 = ||CA1 - CA2||.  |delta| <= d_CA1 + d_CA2 <= T * sqrt(2n).

    Both are **tight**: displacing the CA outward while every CG atom moves
    inward (fp0), or pushing the two CAs apart along their own axis (fp2),
    reaches the bound exactly. Anything the prefilter rejects at these
    tolerances is therefore provably outside the RMSD threshold -- rejecting a
    true match is not possible, which is what lets the caller skip the Kabsch
    fit outright.

    The returned values carry an extra `FP_SAFETY_FACTOR` on top of that bound.
    Loosening an admissible tolerance cannot lose a match, so the guarantee is
    unaffected; it only costs pruning.

    Each tolerance assumes the *whole* SSD budget lands on the atoms its own
    fingerprint touches, so the three are individually necessary conditions and
    are correct applied one at a time. They are **not** simultaneously
    achievable: do not try to tighten them into a joint constraint.

    Note the tolerance always exceeds `rmsd_threshold` itself (n >= 4 gives a
    factor of at least ~2.2), which is what makes `total + tol` a safe
    "reject" sentinel for callers accumulating a distance.
    """
    if n_cg <= 0 or n_total <= 0:
        raise ValueError(f'fp_tolerances: n_total={n_total}, n_cg={n_cg} must be positive')
    scale = rmsd_threshold * FP_SAFETY_FACTOR
    tol_com = scale * math.sqrt(n_total * (1.0 + 1.0 / n_cg))
    if n_res == 1:
        return (tol_com,)
    if n_res == 2:
        return (tol_com, tol_com, scale * math.sqrt(2.0 * n_total))
    raise ValueError(f'fp_tolerances: subset size {n_res} not supported '
                     f'(defined for {_MIN_SUBSET_SIZE}-{_MAX_SUBSET_SIZE})')


def precompute_bucket_fingerprints(cgvdmbb_data, n_cg):
    """
    Vectorized fingerprint precomputation for all vdGs in one bucket.

    cgvdmbb layout: [CG_atoms (n_cg) | N1, CA1, C1 | (N2, CA2, C2) ...]

    Returns a dict:
      size_subset==1: {'fp0': (N,) float32}
      size_subset==2: {'fp0': (N,), 'fp1': (N,), 'fp2': (N,)}
        fp0 = ||CA1 - CG_COM||,  fp1 = ||CA2 - CG_COM||,  fp2 = ||CA1 - CA2||
    Returns an empty dict when there is nothing to fingerprint (no data, n_cg
    == 0, or no backbone atoms); callers test truthiness before indexing.
    Raises ValueError on a backbone block that is not 1 or 2 whole residues.
    """
    if cgvdmbb_data is None or len(cgvdmbb_data) == 0 or n_cg == 0:
        return {}
    arr = cgvdmbb_data.astype(np.float32, copy=False)
    n_res, rem = divmod(arr.shape[1] - n_cg, 3)
    if n_res == 0 and rem == 0:   # no backbone atoms: nothing to fingerprint
        return {}
    if rem != 0 or n_res > _MAX_SUBSET_SIZE:
        raise ValueError(f'precompute_bucket_fingerprints: {arr.shape[1] - n_cg} backbone '
                         f'atoms is not a supported subset size (defined for '
                         f'{_MIN_SUBSET_SIZE}-{_MAX_SUBSET_SIZE} residues, 3 atoms each)')
    cg_com = arr[:, :n_cg, :].mean(axis=1)
    ca1    = arr[:, n_cg + 1, :]
    fp0    = np.linalg.norm(ca1 - cg_com, axis=1).astype(np.float32)
    if n_res == 2:
        ca2 = arr[:, n_cg + 4, :]
        return {
            'fp0': fp0,
            'fp1': np.linalg.norm(ca2 - cg_com, axis=1).astype(np.float32),
            'fp2': np.linalg.norm(ca1 - ca2,    axis=1).astype(np.float32),
        }
    return {'fp0': fp0}


def cross_residue_mask(n_res):
    """Upper-triangular mask over backbone atoms, selecting cross-residue pairs.

    Within-residue N-CA/CA-C/N-C distances are near-invariant between any two
    vdGs, so they cost descriptor width without pruning anything; the mask is
    what drops them. It is applied *after* a slot permutation, and stays correct
    under one because a slot permutation maps whole residue blocks onto whole
    residue blocks, so cross-residue pairs permute among themselves.
    """
    n_bb = 3 * n_res
    res_of = np.repeat(np.arange(n_res), 3)
    i, j = np.triu_indices(n_bb, k=1)
    keep = res_of[i] != res_of[j]
    mask = np.zeros((n_bb, n_bb), dtype=bool)
    mask[i[keep], j[keep]] = True
    return mask


def precompute_internal_distance_descriptors(cgvdmbb_data, n_cg):
    """Precompute the lean Stage-1 internal-distance descriptor.

    Coordinates are laid out as ``CG | N1,CA1,C1 | [N2,CA2,C2]``. The descriptor
    holds every CG-atom-to-backbone-atom distance (``cg_bb``) and, for a
    two-residue vdG, the whole backbone-to-backbone distance matrix (``bb_bb``)
    plus the ``cross_mask`` that selects its cross-residue entries. CG-CG
    distances are omitted: they added at most ~0.13 percentage points of
    rejection in the validation benchmark.

    ``bb_bb`` is kept as a full symmetric matrix rather than the nine
    cross-residue distances alone because a vdM slot permutation relabels
    backbone atoms, and the bound has to index it on *both* axes. A
    pre-flattened cross-only vector cannot be permuted.

    Values are stored as float32; lower-bound accumulation is float64 in
    :func:`internal_distance_lower_bounds`.
    """
    if cgvdmbb_data is None or len(cgvdmbb_data) == 0 or n_cg <= 0:
        return {}
    arr = np.asarray(cgvdmbb_data, dtype=np.float32)
    if arr.ndim != 3 or arr.shape[2] != 3:
        raise ValueError(
            f'precompute_internal_distance_descriptors expected (N, atoms, 3), '
            f'got {arr.shape}')
    if n_cg >= arr.shape[1]:
        return {}
    n_res, rem = divmod(arr.shape[1] - n_cg, 3)
    if rem != 0 or not _MIN_SUBSET_SIZE <= n_res <= _MAX_SUBSET_SIZE:
        raise ValueError(
            f'precompute_internal_distance_descriptors: {arr.shape[1] - n_cg} '
            f'backbone atoms is not a supported subset size (defined for '
            f'{_MIN_SUBSET_SIZE}-{_MAX_SUBSET_SIZE} residues, 3 atoms each)')

    cg = arr[:, :n_cg]
    bb = arr[:, n_cg:]
    n_records, n_bb = len(arr), bb.shape[1]
    cg_bb = np.empty((n_records, n_cg, n_bb), dtype=np.float32)
    # Loop over the three/six backbone atoms to avoid an N*n_cg*n_bb*3
    # temporary, which matters for the largest buckets on Wynton.
    for j in range(n_bb):
        delta = cg - bb[:, None, j]
        cg_bb[:, :, j] = np.sqrt(
            np.sum(delta * delta, axis=2, dtype=np.float32))

    out = {'cg_bb': cg_bb}
    if n_res == 2:
        bb_bb = np.zeros((n_records, n_bb, n_bb), dtype=np.float32)
        for i in range(n_bb):
            for j in range(i + 1, n_bb):
                delta = bb[:, i] - bb[:, j]
                d = np.sqrt(np.sum(delta * delta, axis=1, dtype=np.float32))
                bb_bb[:, i, j] = d
                bb_bb[:, j, i] = d
        out['bb_bb'] = bb_bb
        out['cross_mask'] = cross_residue_mask(n_res)
    return out


def internal_distance_lower_bounds(query_index, target_indices, descriptors,
                                   n_total, perm_group=None):
    """Return symmetry-minimised RMSD lower bounds for one query vs targets.

    For selected internal-distance changes ``delta`` over ``n_total`` atoms,
    both ``||delta||_2 / n_total`` and ``max(abs(delta)) / sqrt(2*n_total)``
    lower-bound the best-fit RMSD. The maximum of those is taken per group
    element, then minimised over the group. A pair may therefore be rejected
    when the returned value exceeds the RMSD cutoff, with no Kabsch fit.

    ``perm_group`` is the vdG's *full* symmetry group as ``(cg_perm, bb_perm)``
    index pairs -- CG automorphisms crossed with interchangeable vdM slot
    orders. Minimising over all of it is what keeps the bound admissible: the
    exact RMSD is itself a minimum over the same group, so a bound taken over
    any subset could exceed the true distance and reject a real match.
    """
    targets = np.asarray(target_indices, dtype=np.intp).reshape(-1)
    if targets.size == 0:
        return np.empty(0, dtype=np.float64)
    if not descriptors or 'cg_bb' not in descriptors:
        return np.zeros(targets.size, dtype=np.float64)
    if n_total <= 0:
        raise ValueError(f'n_total must be positive, got {n_total}')

    cg_bb = descriptors['cg_bb']
    if not 0 <= query_index < len(cg_bb):
        raise IndexError(f'query_index {query_index} outside 0..{len(cg_bb) - 1}')
    n_cg, n_bb = cg_bb.shape[1], cg_bb.shape[2]
    if perm_group is None:
        perm_group = ((np.arange(n_cg, dtype=np.intp),
                       np.arange(n_bb, dtype=np.intp)),)
    else:
        perm_group = validate_perm_group(perm_group, n_cg, n_bb)

    target_cg = cg_bb[targets].reshape(targets.size, -1)
    query_cg = cg_bb[query_index]
    bb_bb = descriptors.get('bb_bb')
    if bb_bb is not None:
        mask = descriptors['cross_mask']
        target_cross = bb_bb[targets][:, mask]
        query_bb = bb_bb[query_index]

    best = np.full(targets.size, np.inf, dtype=np.float64)
    max_denom = math.sqrt(2.0 * n_total)
    for cg_perm, bb_perm in perm_group:
        query = np.asarray(query_cg[cg_perm][:, bb_perm],
                           dtype=np.float64).reshape(-1)
        delta = query - target_cg
        sum_sq = np.sum(delta * delta, axis=1, dtype=np.float64)
        max_abs = np.max(np.abs(delta), axis=1)
        if bb_bb is not None:
            cross_delta = np.asarray(
                query_bb[bb_perm][:, bb_perm][mask], dtype=np.float64) - target_cross
            sum_sq += np.sum(cross_delta * cross_delta, axis=1, dtype=np.float64)
            max_abs = np.maximum(max_abs,
                                 np.max(np.abs(cross_delta), axis=1))
        bound = np.maximum(np.sqrt(sum_sq) / n_total,
                           max_abs / max_denom)
        np.minimum(best, bound, out=best)
    return best


def validate_descriptor_permutations(permutations, n_cg):
    """Validate CG permutations without importing the higher-level utils module."""
    normalized = []
    expected = tuple(range(n_cg))
    for permutation in permutations:
        perm = tuple(int(i) for i in permutation)
        if len(perm) != n_cg or tuple(sorted(perm)) != expected:
            raise ValueError(f'Invalid CG descriptor permutation {perm!r}')
        normalized.append(np.asarray(perm, dtype=np.intp))
    if not normalized:
        raise ValueError('At least one CG descriptor permutation is required')
    return tuple(normalized)


def _validate_index_permutation(perm, size, label):
    perm = np.asarray(perm, dtype=np.intp).reshape(-1)
    if perm.size != size or not np.array_equal(np.sort(perm), np.arange(size)):
        raise ValueError(f'Invalid {label} permutation {perm.tolist()!r} for size {size}')
    return perm


def validate_perm_group(perm_group, n_cg, n_bb):
    """Validate a symmetry group given as ``(cg_perm, bb_perm)`` index pairs."""
    normalized = []
    for entry in perm_group:
        try:
            cg_perm, bb_perm = entry
        except (TypeError, ValueError):
            raise ValueError(
                f'perm_group entries must be (cg_perm, bb_perm) pairs, got {entry!r}')
        normalized.append((_validate_index_permutation(cg_perm, n_cg, 'CG'),
                           _validate_index_permutation(bb_perm, n_bb, 'backbone')))
    if not normalized:
        raise ValueError('At least one symmetry group element is required')
    return tuple(normalized)


def slot_orders(aa_bucket_parts):
    """Residue orderings that permute only vdM slots carrying the same label.

    A bucket's label list is what says which slots are interchangeable -- two
    slots may be swapped exactly when their labels match, so ``ASP_ASP`` and
    ``bb_bb`` admit both orders while ``ARG_bb`` admits only the identity.

    Delegates to ``utils.group_preserving_permutations``, the same routine
    ``vdg_npz_utils.aa_perm_indices`` uses on the read path. Generation and hit
    finding have to agree on which slots are interchangeable, so they share one
    implementation rather than two that happen to match.
    """
    from ligand_vdgs.functions.utils import group_preserving_permutations
    return sorted(tuple(order)
                  for order in group_preserving_permutations(list(aa_bucket_parts)))


def build_perm_group(cg_symm_perms, n_cg, aa_bucket_parts):
    """The vdG's full symmetry group: CG automorphisms x interchangeable slots.

    Returns ``(cg_perm, bb_perm)`` index pairs. Handling slot symmetry here --
    inside the distance -- rather than by replicating each vdG once per slot
    ordering keeps one physical environment per record, so a vdG cannot end up
    in two clusters and no post-hoc union is needed to put it back together.
    """
    n_res = len(aa_bucket_parts)
    # `cg_symm_perms` arrives as a numpy array from the library's stored
    # automorphism table, so emptiness has to be tested by length, not truthiness.
    if cg_symm_perms is None or len(cg_symm_perms) == 0:
        cg_perms = (np.arange(n_cg, dtype=np.intp),)
    else:
        cg_perms = validate_descriptor_permutations(cg_symm_perms, n_cg)
    group = []
    for order in slot_orders(aa_bucket_parts):
        bb_perm = np.concatenate(
            [np.arange(3 * r, 3 * r + 3) for r in order]).astype(np.intp) \
            if n_res else np.empty(0, dtype=np.intp)
        for cg_perm in cg_perms:
            group.append((cg_perm, bb_perm))
    return tuple(group)


def full_row_permutations(perm_group, n_cg):
    """Whole-row atom permutations for the exact Kabsch fit.

    The cgvdmbb row is ``CG | N1,CA1,C1 | [N2,CA2,C2]``, so a group element
    relabels the CG block and the backbone block independently; concatenating
    them gives the single index array a batched fit can apply.
    """
    return tuple(np.concatenate([cg_perm, n_cg + bb_perm]).astype(np.intp)
                 for cg_perm, bb_perm in perm_group)


# ---- Vectorised helpers used by score_poses/hit_finder_core.py ------------
# Tolerances come from fp_tolerances(); each fingerprint gets its own.

def prefilter_query_indices_single(fp_q0, fp_v0, tol0):
    """Return flatnonzero of |fp_q0 - fp_v0| <= tol0, or None if empty."""
    mask = np.abs(fp_q0 - fp_v0) <= tol0
    return np.flatnonzero(mask) if np.any(mask) else None


def prefilter_query_indices_pair(fp_q0, fp_q1, fp_q2,
                                  fp_v0, fp_v1, fp_v2,
                                  tols):
    """Vectorised prefilter for size_subset==2 (hit-finder usage).

    `tols` is the 3-tuple from `fp_tolerances`. fp2 is tested first: it is the
    CA-CA distance, which discriminates best in practice.
    """
    tol0, tol1, tol2 = tols
    mask = np.abs(fp_q2 - fp_v2) <= tol2
    if not np.any(mask):
        return None
    mask &= np.abs(fp_q0 - fp_v0) <= tol0
    if not np.any(mask):
        return None
    mask &= np.abs(fp_q1 - fp_v1) <= tol1
    return np.flatnonzero(mask) if np.any(mask) else None
