"""Geometry helpers shared by vdG hit scoring modes."""

import numpy as np

from ligand_vdgs.functions.utils import kabsch

PLACEMENT_CHUNK = 2048
BB_SLOT_SIDECHAIN_CLASH = 3.4
PRO_NH_DONOR_CUTOFF = 3.5
ACCEPTOR_ELEMENTS = ("N", "O", "S")

def cg_center(cg_coords):
    return np.asarray(cg_coords, dtype=np.float32).mean(axis=0, dtype=np.float32)

def _dist(a, b):
    d = a - b
    return float(np.sqrt(np.sum(d * d, dtype=np.float32), dtype=np.float32))

def fp_single_from_ca_and_cgcom(ca_coords, cg_com):
    return _dist(ca_coords[0], cg_com)

def fp_pair_from_ca_and_cgcom(ca_coords, cg_com):
    return (_dist(ca_coords[0], cg_com), _dist(ca_coords[1], cg_com),
            _dist(ca_coords[0], ca_coords[1]))

def _min_dist2(cgs, pts):
    """Minimum squared distance from each CG to any point in ``pts``."""
    cgs, pts = np.asarray(cgs, np.float32), np.asarray(pts, np.float32).reshape(-1, 3)
    return np.concatenate([np.sum((cgs[s:s + PLACEMENT_CHUNK, :, None] - pts) ** 2, axis=-1,
                                  dtype=np.float32).min(axis=(1, 2))
                           for s in range(0, len(cgs), PLACEMENT_CHUNK)] or [np.empty(0, np.float32)])

def backbone_slots_host_mask(blockers, cgs, cg_is_acceptor):
    """Return a mask for CGs that fit the backbone-only slot blockers."""
    cgs = np.asarray(cgs, np.float32)
    acc = None if cg_is_acceptor is None else np.broadcast_to(cg_is_acceptor, cgs.shape[:2])
    checks = [np.ones(len(cgs), bool)]
    checks.extend(_min_dist2(cgs, sc) >= BB_SLOT_SIDECHAIN_CLASH ** 2
                  for sc, _ in blockers or () if len(sc))
    checks.extend(np.where(acc, np.sum((cgs - n) ** 2, axis=-1), np.inf).min(axis=1)
                  >= PRO_NH_DONOR_CUTOFF ** 2 for _sc, n in blockers or ()
                  if n is not None and acc is not None and acc.any())
    return np.logical_and.reduce(checks)

def backbone_slots_can_host(blockers, cg_coords, cg_is_acceptor):
    return bool(backbone_slots_host_mask(blockers, np.asarray(cg_coords)[None], cg_is_acceptor)[0])

def has_any_bsr_atom_cg_contact(bsr_atom_coords, cg_coords, cutoff=3.8):
    if bsr_atom_coords is None or len(bsr_atom_coords) == 0:
        return False
    cg_coords = np.asarray(cg_coords, dtype=np.float32)
    if cg_coords.size == 0:
        return False
    diff = bsr_atom_coords[:, None, :] - cg_coords[None, :, :]
    return bool(np.any(np.sum(diff * diff, axis=2, dtype=np.float32) <= np.float32(cutoff * cutoff)))

def cg_bb_dist(q_cg, bb_flat):
    """Mean CG atom distance to the nearest fit backbone atom."""
    return float(np.linalg.norm(q_cg[:, None, :] - bb_flat[None, :, :], axis=-1).min(axis=1).mean())

def held_out_placement(vdg_bb, vdg_cg, bb_flat, q_cgs):
    """Fit vdG backbones, then compare their placed CGs with query CG labelings."""
    R, t, _ = kabsch(vdg_bb, np.broadcast_to(bb_flat, vdg_bb.shape))
    placed = (np.einsum("mij,mjk->mik", np.asarray(vdg_cg, np.float32), R) + t[:, None, :]).astype(np.float32)
    q_cgs = np.asarray(q_cgs, np.float32)
    ssd = np.concatenate([np.sum((placed[s:s + PLACEMENT_CHUNK, None] - q_cgs[None]) ** 2,
                                 axis=(2, 3), dtype=np.float64)
                          for s in range(0, len(placed), PLACEMENT_CHUNK)])
    k = ssd.argmin(axis=1)
    return R, t, np.sqrt(ssd[np.arange(len(k)), k] / q_cgs.shape[1]), k, placed

def proper_rotation_mask(R):
    return ((np.abs(np.linalg.det(R.astype(np.float64)) - 1.0) <= 1e-5)
            & (np.abs(np.einsum("mji,mjk->mik", R, R) - np.eye(3)).max(axis=(1, 2)) <= 1e-4))
