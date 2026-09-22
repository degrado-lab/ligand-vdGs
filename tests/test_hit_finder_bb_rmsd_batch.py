"""B4 perf fix (2026-09-15): the backbone-only lower-bound RMSD used inside
score_one_model's scoring loop was one kabsch_ssd call per (nr_idx, aa_perm)
pair -- profiled as the dominant cost for 4ij8/SAM (60% of runtime, ~15M
calls). Batching it into one call per bucket must reproduce the exact
per-(nr_idx, aa_perm) value, not just an aggregate."""
import numpy as np
from ligand_vdgs.functions.utils import kabsch_ssd

def _naive_bb_rmsd(bucket_bb, resind_perms, bb_flat, n_atoms):
    """Exactly the per-item code being replaced (hit_finder_core.py:1019-1020)."""
    n_nr_vdgs = bucket_bb.shape[0]
    out = np.empty((n_nr_vdgs, len(resind_perms)), dtype=np.float64)
    for nr_idx in range(n_nr_vdgs):
        for aa_perm_idx, resind_perm in enumerate(resind_perms):
            out[nr_idx, aa_perm_idx] = np.sqrt(kabsch_ssd(
                bucket_bb[nr_idx][np.asarray(resind_perm, dtype=np.intp)].reshape(-1, 3),
                bb_flat[None])[0] / n_atoms)
    return out

def _batched_bb_rmsd(bucket_bb, resind_perms, bb_flat, n_atoms):
    """The replacement: one kabsch_ssd call for the whole bucket."""
    n_nr_vdgs, n_perms = bucket_bb.shape[0], len(resind_perms)
    return np.sqrt(kabsch_ssd(
        bb_flat,
        bucket_bb[:, np.asarray(resind_perms, dtype=np.intp), :, :].reshape(
            n_nr_vdgs * n_perms, -1, 3),
    ) / n_atoms).reshape(n_nr_vdgs, n_perms)

def test_batched_bb_rmsd_matches_naive_per_item():
    rng = np.random.default_rng(0)
    n_nr_vdgs, subset_size, atoms_per_res = 5, 2, 3
    # Deliberately asymmetric per residue-slot coords -- a bug that ignored
    # aa_perm, or transposed the (nr_idx, aa_perm) axes, disagrees per-cell
    # here, not just in aggregate.
    bucket_bb = rng.normal(size=(n_nr_vdgs, subset_size, atoms_per_res, 3)).astype(np.float32)
    resind_perms = [(0, 1), (1, 0)]
    bb_flat = rng.normal(size=(subset_size * atoms_per_res, 3)).astype(np.float32)
    n_atoms = subset_size * atoms_per_res + 4  # stand-in for the real bb+CG atom count

    naive = _naive_bb_rmsd(bucket_bb, resind_perms, bb_flat, n_atoms)
    batched = _batched_bb_rmsd(bucket_bb, resind_perms, bb_flat, n_atoms)

    assert naive.shape == batched.shape == (n_nr_vdgs, len(resind_perms))
    np.testing.assert_allclose(batched, naive, atol=1e-5)
    # Non-vacuity: the two aa_perms must give genuinely different values for
    # at least one nr_idx, or a permutation-ignoring bug (e.g. broadcasting
    # perm 0's slice for every column) would pass the allclose above too.
    assert not np.allclose(naive[:, 0], naive[:, 1]), (
        "fixture's two aa_perms give identical bb_rmsd for every nr_idx -- "
        "cannot discriminate a bug that ignores aa_perm; regenerate with a "
        "different seed/asymmetry.")
