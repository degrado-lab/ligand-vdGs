"""
held_out_cg_rmsd must be a genuinely held-out quantity (DR-17).

Falsifier: this file goes red if `compute_held_out_cg_rmsd` reads the joint
bb+CG transform (R/t) instead of the backbone-only one (Rbb/tbb), or if it
scores the hit's recorded `q_cg_perm_idx` instead of minimizing over the CG
automorphism group.

Discriminating input: a carboxylate-like CG with a 2-fold automorphism, whose
recorded perm is deliberately the WRONG relabeling, and a fixture in which the
joint and backbone-only fits provably diverge.
"""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ligand_vdgs.functions.utils import kabsch
from ligand_vdgs.tools.benchmark_pose_recovery import compute_held_out_cg_rmsd

FRAG, SUBSET, SIGN, BUCKET, VDG_IDX, SITE = "CG_frag", 1, "neut", "ASP", 0, 0

# Carboxylate-like CG: the O1/O2 swap is the automorphism.
LIB_CG = np.array([[0., 0., 0.], [1.25, 0.6, 0.], [1.25, -0.6, 0.]], np.float32)

def _rec(R, t, R_bb, t_bb, q_cg_perm_idx):
    rec = dict(frag=FRAG, subset_size=SUBSET, charge_sign=SIGN, aa_bucket=BUCKET,
               vdg_index=VDG_IDX,
               q_site_idx=SITE, q_cg_perm_idx=q_cg_perm_idx)
    for i in range(3):
        for j in range(3):
            rec[f"R{i}{j}"] = f"{R[i, j]:.4f}"
            rec[f"Rbb{i}{j}"] = f"{R_bb[i, j]:.4f}"
    for k in range(3):
        rec[f"t{k}"] = f"{t[k]:.4f}"
        rec[f"tbb{k}"] = f"{t_bb[k]:.4f}"
    return rec

def _rot_z(deg):
    a = np.radians(deg); c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]], np.float32)

def _rot_axis(axis, deg):
    """Rodrigues rotation about an arbitrary (non-principal) axis."""
    k = axis / np.linalg.norm(axis)
    a = np.radians(deg)
    K = np.array([[0., -k[2], k[1]], [k[2], 0., -k[0]], [-k[1], k[0], 0.]], np.float32)
    return (np.eye(3, dtype=np.float32) + np.sin(a) * K
            + (1 - np.cos(a)) * (K @ K)).astype(np.float32)

def _rmsd(a, b): return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))

@pytest.fixture
def seeded_bucket():
    return {(FRAG, SUBSET, SIGN, BUCKET): {"cg": LIB_CG[None]}}

def test_uses_backbone_only_transform_not_the_joint_one(seeded_bucket):
    """R maps the library CG exactly onto the crystal CG (residual 0); Rbb does not.
    Reading R instead of Rbb returns 0.0 and is the bug this catches."""
    R_joint, t_joint = _rot_z(0.0), np.zeros(3, np.float32)
    crystal = (LIB_CG @ R_joint + t_joint).astype(np.float32)
    # Oblique axis + a y-offset: LIB_CG is mirror-symmetric under y -> -y, so a
    # pure z-rotation with an x-only shift gives R and R.T the same RMSD and the
    # test could not see a transposed rotation (CLAUDE.md: using R.T is silent).
    R_bb = _rot_axis(np.array([0.3, 0.5, 0.81], np.float32), 25.0)
    t_bb = np.array([0.4, 0.55, -0.2], np.float32)

    cache = {FRAG: {SITE: {0: crystal}}}
    got = compute_held_out_cg_rmsd(_rec(R_joint, t_joint, R_bb, t_bb, 0),
                                   cache, "/unused", seeded_bucket)
    expected = _rmsd(LIB_CG @ R_bb + t_bb, crystal)
    # Vacuity clause: the two transforms must give genuinely different answers,
    # otherwise this asserts nothing.
    assert expected > 0.1, f"fixture degenerate: bb-only residual {expected:.4f}"
    # Vacuity clause: the fixture must be able to SEE a transposed rotation.
    transposed = _rmsd(LIB_CG @ R_bb.T + t_bb, crystal)
    assert abs(transposed - expected) > 0.05, (
        f"fixture cannot distinguish R from R.T ({transposed:.4f} vs {expected:.4f})")
    assert got == pytest.approx(expected, abs=1e-4), (
        f"got {got:.4f}, expected {expected:.4f} (0.0 means it read the joint R/t)")

def test_minimizes_over_cg_automorphisms_not_the_recorded_perm(seeded_bucket):
    """The recorded q_cg_perm_idx is the BAD relabeling. Indexing by it instead of
    minimizing over the site's perms reports the inflated RMSD."""
    R_bb, t_bb = _rot_z(0.0), np.zeros(3, np.float32)
    good = LIB_CG.copy()                      # matches the placed CG
    bad = LIB_CG[[0, 2, 1]].copy()            # O1/O2 swapped -> inflated RMSD
    bad[1] += np.array([0., 0., 3.0], np.float32)   # force a large gap

    cache = {FRAG: {SITE: {0: bad, 1: good}}}
    got = compute_held_out_cg_rmsd(_rec(R_bb, t_bb, R_bb, t_bb, q_cg_perm_idx=0),
                                   cache, "/unused", seeded_bucket)
    inflated = _rmsd(LIB_CG @ R_bb + t_bb, bad)
    # Vacuity clause: the wrong perm must actually be wrong.
    assert inflated > 1.0, f"fixture degenerate: bad perm residual {inflated:.4f}"
    assert got == pytest.approx(0.0, abs=1e-4), (
        f"got {got:.4f}; {inflated:.4f} means it indexed by q_cg_perm_idx")

def test_joint_and_backbone_only_fits_provably_diverge():
    """The premise the metric rests on: when the CG is displaced relative to the
    backbone, the joint fit trades backbone error for CG error, so its CG residual
    is strictly smaller than the backbone-only fit's. If these two fits coincided,
    'held out' would be a distinction without a difference."""
    lib_bb = np.array([[0., 0., 0.], [1.5, 0., 0.], [2.0, 1.4, 0.],
                       [0.3, 2.1, 0.6]], np.float32)
    R_true, t_true = _rot_z(35.0), np.array([2., -1., 0.5], np.float32)
    q_bb = (lib_bb @ R_true + t_true).astype(np.float32)
    # CG placed by a DIFFERENT transform: no rigid motion fits both exactly.
    q_cg = (LIB_CG @ _rot_z(70.0) + t_true + np.array([0.9, 0.4, -0.7], np.float32)
            ).astype(np.float32)

    X_joint = np.concatenate((lib_bb, LIB_CG), axis=0).astype(np.float32)
    Y_joint = np.concatenate((q_bb, q_cg), axis=0).astype(np.float32)[None]
    R_j, t_j, _ = kabsch(X_joint, Y_joint)
    R_b, t_b, _ = kabsch(lib_bb, q_bb[None])

    cg_joint = _rmsd(LIB_CG @ R_j[0] + t_j[0], q_cg)
    cg_bb = _rmsd(LIB_CG @ R_b[0] + t_b[0], q_cg)

    # Vacuity clause: assert the two FITS differ, not merely the two numbers.
    assert np.linalg.norm(R_j[0] - R_b[0]) > 1e-3, (
        "joint and bb-only rotations coincide; fixture does not exercise the split")
    assert cg_joint < cg_bb - 1e-4, (
        f"in-sample CG residual {cg_joint:.4f} not strictly below held-out {cg_bb:.4f}")

def test_every_emitted_record_key_is_writable_by_the_fixed_column_writer():
    """vdg_hit_finder writes hits with a csv.DictWriter over a hardcoded
    RESULT_FIELDS list, whose default extrasaction='raise'. Adding a field to
    score_one_model's record without adding it there crashes every real run at
    the write step -- which is exactly what the Rbb/tbb addition did."""
    import csv, io as _io, re
    from ligand_vdgs.score_poses.vdg_hit_finder import RESULT_FIELDS

    src = _io.open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "ligand_vdgs/score_poses/hit_finder_core.py"),
        encoding="utf-8").read()
    block = src[src.index("rec = dict("):src.index("match_records.append(rec)")]
    keys = set(re.findall(r"^\s+(\w+)=", block, re.M))
    keys |= {f"R{i}{j}" for i in range(3) for j in range(3)}
    keys |= {f"t{k}" for k in range(3)}
    keys |= {f"Rbb{i}{j}" for i in range(3) for j in range(3)}
    keys |= {f"tbb{k}" for k in range(3)}

    # Vacuity clause: if the scrape found nothing, the assertions below are free.
    assert "Rbb00" in keys and "vdg_rmsd" in keys and len(keys) >= 40, (
        f"record-key scrape looks wrong ({len(keys)} keys); test would be vacuous")

    missing = sorted(keys - set(RESULT_FIELDS))
    assert not missing, f"emitted but absent from RESULT_FIELDS: {missing}"
    csv.DictWriter(_io.StringIO(), fieldnames=RESULT_FIELDS,
                   delimiter="\t").writerow({k: "0" for k in keys})

def test_rejects_a_corrupted_backbone_rotation_instead_of_extrapolating(seeded_bucket):
    """_is_proper_rotation only catches a malformed R (bad det / non-orthogonal),
    e.g. from TSV read/write corruption -- NOT a genuinely rank-deficient fit.
    utils.kabsch's reflection correction always returns a proper rotation by
    construction, so a real rank-deficient (collinear) backbone fit passes this
    guard silently; that hazard is untested and unclosed (private/project_planning/
    vacuity_audit_2026-09-12.md, S2). This test only covers the corrupted-matrix
    case: held_out applies Rbb to the CG, which lies outside the fit, so a
    corrupted Rbb must yield None, not a plausible-looking number."""
    from ligand_vdgs.tools.benchmark_pose_recovery import _is_proper_rotation
    crystal = LIB_CG.copy()
    cache = {FRAG: {SITE: {0: crystal}}}
    good = _rot_axis(np.array([0.3, 0.5, 0.81], np.float32), 25.0)
    assert _is_proper_rotation(good), "a valid rotation must not be rejected"

    for name, bad in [("scaled", good * 1.05), ("reflection", good * np.array(
            [[1., 1., 1.], [1., 1., 1.], [-1., -1., -1.]], np.float32))]:
        assert not _is_proper_rotation(bad), f"{name} matrix accepted as a rotation"
        got = compute_held_out_cg_rmsd(
            _rec(good, np.zeros(3, np.float32), bad, np.zeros(3, np.float32), 0),
            cache, "/unused", seeded_bucket)
        assert got is None, f"{name}: returned {got}, should refuse to extrapolate"

def test_storage_rounding_alone_does_not_trip_the_rotation_guard():
    """Vacuity guard on the tolerance: it must be loose enough for 4dp rounding but
    tight enough to still reject a 5% scale. If rounding tripped it, every real hit
    would silently return None and held_out_cg_rmsd would be all-NA."""
    from ligand_vdgs.tools.benchmark_pose_recovery import _is_proper_rotation
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(3000):
        A = rng.normal(size=(3, 3))
        U, _, Vt = np.linalg.svd(A)
        R = U @ Vt
        if np.linalg.det(R) < 0: R = U @ np.diag([1., 1., -1.]) @ Vt
        Rr = np.round(R, 4)
        worst = max(worst, abs(np.linalg.det(Rr) - 1.0))
        assert _is_proper_rotation(Rr), "4dp-rounded proper rotation was rejected"
    assert worst > 5e-5, f"rounding envelope not exercised (worst {worst:.2e})"
