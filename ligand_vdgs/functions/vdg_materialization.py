"""Alignment and selection helpers for materialized vdG structures."""
import numpy as np

from ligand_vdgs.functions import utils
from ligand_vdgs.functions import vdg_npz_utils as vdg_npz

FIRST_FRAME = np.array([[0., 0., 0.], [-1., 0., 1.], [1., -1., 0.]])
RES_FIELDS = ("resname", "seg", "chain", "resnum")

def weights(n_cg, n_bb, cg_weight):
    """Spread the requested weight evenly over CG and backbone atoms."""
    return np.r_[np.full(n_cg, cg_weight / max(n_cg, 1)),
                 np.full(n_bb, max(1. - cg_weight, 0.) / max(n_bb, 1))]

def fit(mobile, target, weights=None):
    """Weighted Kabsch fit using row-vector coordinates: mobile @ R + t."""
    w = np.ones(len(mobile)) if weights is None else np.asarray(weights, float)
    w = w / w.sum()
    mob_cent, tar_cent = w @ mobile, w @ target
    R = utils._proper_kabsch_rotations((mobile - mob_cent).T @ ((target - tar_cent) * w[:, None]))
    return R, tar_cent - mob_cent @ R

def best_fit(cg, bb, target, atom_weights, perms):
    """Fit CG plus optional backbone, minimizing over CG automorphisms."""
    cg, target = np.asarray(cg, float), np.asarray(target, float)
    w = np.ones(len(target)) if atom_weights is None else np.asarray(atom_weights, float)
    def one(perm):
        mobile = np.vstack([cg[list(perm)], np.reshape([] if bb is None else bb, (-1, 3))])
        R, t = fit(mobile, target, w)
        return w @ np.sum((mobile @ R + t - target) ** 2, axis=1), R, t
    return min((one(p) for p in perms or [range(len(cg))]), key=lambda result: result[0])[1:]

def move(atomgroup, coords, R, t):
    """Transform an atom group and its fitting coordinates."""
    vdg_npz.apply_rigid_transform(atomgroup, R, t)
    return np.asarray(coords, float) @ R + t

def place_on_global(atomgroup, cg, coords, global_ref, perms):
    """Place the first structure in a fixed frame and later structures on its CG."""
    if global_ref is not None:
        return move(atomgroup, coords, *best_fit(cg, None, global_ref, None, perms)), global_ref
    if len(cg) < 3:
        raise ValueError("[ERROR] Need at least 3 CG atoms for first-centroid alignment.")
    coords = move(atomgroup, coords, *fit(np.asarray(cg, float)[:3], FIRST_FRAME))
    return coords, atomgroup.getCoords()[:len(cg)].copy()

def rank_clusters(num_parents, cluster_size, min_cluster_size=1, top_n=None):
    """Rank by distinct parent support, with stored order breaking ties."""
    order = np.argsort(-np.asarray(num_parents), kind="stable")
    if min_cluster_size > 1:
        order = order[np.asarray(cluster_size)[order] >= min_cluster_size]
    return order[:top_n]

def cap_across_buckets(bucket_supports, max_files=None):
    """Count each bucket's leading picks retained by a global support cap."""
    if max_files is None:
        return [len(s) for s in bucket_supports]
    return np.bincount([b for _, b, _ in sorted(
        (-int(support), b, rank) for b, values in enumerate(bucket_supports)
        for rank, support in enumerate(values))[:max_files]],
        minlength=len(bucket_supports)).tolist()
