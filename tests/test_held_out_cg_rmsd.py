"""
DR-36 placement mode and the per-hit held-out CG quantities (hit_finder_core).

Falsifiers:
- held_out_cg_rmsd applies the joint bb+CG transform, or R.T, instead of the backbone-only
  fit -> wrong per-vdG value (checked against an independent SVD Kabsch).
- held_out_cg_rmsd scores one stored labeling instead of the min over query CG labelings ->
  the exact vdG reports > 0 (the first labeling is deliberately the wrong one).
- placement matching still reads the query CG (fp prefilter, joint RMSD) -> the vdG with an
  exact backbone and a displaced CG goes missing from placement hits.
- joint hits are not a subset of placement hits at equal tau (DR-36 cutoff).
"""
import csv, io
import numpy as np
import pytest

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions.utils import kabsch, normalize_rmsd
from ligand_vdgs.score_poses import hit_finder_core as hfc
from ligand_vdgs.score_poses.vdg_hit_finder import RESULT_FIELDS
from tests.vacuity import assert_discriminates

LIB_BB = np.array([[0., 0., 0.], [1.46, 0., 0.], [2.0, 1.4, 0.1]], np.float32)
LIB_CG = np.array([[3.0, 3.0, 0.5], [4.2, 3.2, 0.6], [2.6, 4.1, 0.9], [2.2, 2.3, 1.4]], np.float32)
SWAP = [0, 2, 1, 3]  # the labeling automorphism
N_BB, N_CG = 3, 4
N_ATOMS = N_BB + N_CG
TAU = normalize_rmsd(N_ATOMS, "cgvdmbb")
CG_SHIFT = np.array([1.5, -1.2, 1.6], np.float32)
BSR = [("", "A", 10)]

def _rot(axis, deg):
    k = np.asarray(axis, float) / np.linalg.norm(axis); a = np.radians(deg)
    K = np.array([[0., -k[2], k[1]], [k[2], 0., -k[0]], [-k[1], k[0], 0.]])
    return (np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K).astype(np.float32)

R_Q, T_Q = _rot([0.3, 0.5, 0.81], 40.0), np.array([5.0, -3.0, 2.0], np.float32)
Q_BB, Q_CG = LIB_BB @ R_Q + T_Q, LIB_CG @ R_Q + T_Q

def _svd_fit(X, Y):
    """Independent Kabsch for Y ~ X @ R + t."""
    xc, yc = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((X - xc).T @ (Y - yc))
    D = np.diag([1., 1., np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt
    return R, yc - xc @ R

def _rmsd(a, b): return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))

def _fit_rmsd(X, Y):
    R, t = _svd_fit(X, Y)
    return _rmsd(X @ R + t, Y)

def _bb_distortion(target):
    """LIB_BB perturbed so its optimal-fit RMSD to LIB_BB lands near `target`."""
    d = np.array([[0.5, -0.3, 0.2], [-0.4, 0.1, -0.5], [0.2, 0.4, 0.3]], np.float32)
    s = min(np.linspace(0.05, 5, 400), key=lambda s: abs(_fit_rmsd(LIB_BB + s * d, LIB_BB) - target))
    return (LIB_BB + s * d).astype(np.float32)

# vdG rows: 0 exact; 1 exact bb + displaced CG; 2 bb far outside; 3 bb between tau and the
# placement cutoff tau*sqrt(n_atoms/N_bb).
VDG_BB = np.stack([LIB_BB, LIB_BB, _bb_distortion(1.3), _bb_distortion(0.62)])
VDG_CG = np.stack([LIB_CG, LIB_CG + CG_SHIFT, LIB_CG, LIB_CG])

def _bucket():
    n = len(VDG_BB)
    return dict(aa_bucket_parts=["ASP"], charge_signs=np.full(n, "neg", dtype="U11"),
                partition_indices=np.arange(n, dtype=np.int32),
                cluster_id=np.arange(n, dtype=np.int32),
                cluster_num_parents=np.array([10, 40, 30, 20], np.int32),
                parent_biounit=np.array([f"p{i}" for i in range(n)]),
                cg_elements=np.broadcast_to(np.array(["C", "O", "N", "C"]), (n, 4)),
                cg=VDG_CG.astype(np.float32), bb=VDG_BB[:, None].astype(np.float32))

def _run(match_mode):
    cache = hfc._BoundedBucketCache()
    cache.put(("/lib", "F", 1, vdg_npz.make_aa_bucket(["ASP"])), _bucket())
    site = [(Q_CG[p].astype(np.float32), Q_CG.mean(0), np.array([0, 1, 1, 0], bool), tuple(p),
             (5, 6, 7, 8)) for p in (SWAP, list(range(4)))]
    recs = hfc._score_one_bsr_combo(
        "q.pdb", None, "F", "F", [site], (["ASP"], BSR, ["ASP"], [Q_BB]), "/lib", cache,
        {}, None, None, "", match_mode)
    return {int(r["vdg_index"]): r for r in recs}

def _expected(i):
    R, t = _svd_fit(VDG_BB[i], Q_BB)
    placed = VDG_CG[i] @ R + t
    return dict(bb=_rmsd(VDG_BB[i] @ R + t, Q_BB),
                held=min(_rmsd(placed, Q_CG[p]) for p in (SWAP, list(range(4)))),
                transposed=_rmsd(VDG_CG[i] @ R.T + t, Q_CG),
                lever=float(np.linalg.norm(Q_CG[:, None] - Q_BB[None], axis=-1).min(1).mean()))

@pytest.fixture(scope="module")
def hits():
    return _run("joint"), _run("bb")

def test_fixture_spans_the_match_window():
    bb = [_expected(i)["bb"] for i in range(4)]
    cutoff = TAU * np.sqrt(N_ATOMS / N_BB)
    assert bb[0] < 1e-4 and bb[1] < 1e-4, bb
    assert bb[2] > cutoff + 0.2, (bb[2], cutoff)
    assert TAU < bb[3] < cutoff, (TAU, bb[3], cutoff)

def test_placement_matches_on_backbone_alone(hits):
    joint, bb = hits
    assert set(bb) == {0, 1, 3}, sorted(bb)
    assert 1 not in joint and 2 not in joint, sorted(joint)
    assert set(joint) <= set(bb), (sorted(joint), sorted(bb))
    assert {r["match_mode"] for r in bb.values()} == {"bb"}
    assert {r["match_mode"] for r in joint.values()} == {"joint"}

def test_per_vdg_values_match_an_independent_fit(hits):
    _, bb = hits
    for i, rec in bb.items():
        e = _expected(i)
        assert float(rec["bb_rmsd"]) == pytest.approx(e["bb"], abs=2e-4), i
        assert float(rec["held_out_cg_rmsd"]) == pytest.approx(e["held"], abs=2e-4), i
        assert float(rec["cg_bb_dist"]) == pytest.approx(e["lever"], abs=2e-4), i
    assert float(bb[1]["held_out_cg_rmsd"]) == pytest.approx(
        float(np.linalg.norm(CG_SHIFT)), abs=2e-4)

def test_held_out_rejects_joint_transform_transpose_and_stored_labeling(hits):
    _, bb = hits
    got1, got0 = float(bb[1]["held_out_cg_rmsd"]), float(bb[0]["held_out_cg_rmsd"])
    X = np.concatenate((VDG_BB[1], VDG_CG[1])); Y = np.concatenate((Q_BB, Q_CG))
    R_j, t_j = _svd_fit(X, Y)
    assert_discriminates(lambda v: abs(v - got1) < 1e-3, [_expected(1)["held"]],
                         [_rmsd(VDG_CG[1] @ R_j + t_j, Q_CG), _expected(1)["transposed"]],
                         "held-out uses the bb-only fit, not joint or R.T")
    assert_discriminates(lambda v: v < 1e-3, [got0], [_rmsd(LIB_CG @ R_Q + T_Q, Q_CG[SWAP])],
                         "held-out minimizes over query CG labelings")

def test_joint_vs_bb_inequalities_hold_on_every_row(hits):
    for rows in hits:
        for i, rec in rows.items():
            v, b, h = (float(rec[k]) for k in ("vdg_rmsd", "bb_rmsd", "held_out_cg_rmsd"))
            assert b <= v * np.sqrt(N_ATOMS / N_BB) + 1e-3, (rec["match_mode"], i, b, v)
            assert N_CG * h ** 2 + N_BB * b ** 2 >= N_ATOMS * v ** 2 - 1e-3, (rec["match_mode"], i)

def test_single_site_joint_rows_carry_the_same_held_out_values(hits):
    joint, bb = hits
    for i in joint:
        for k in ("bb_rmsd", "held_out_cg_rmsd", "cg_bb_dist"):
            assert joint[i][k] == bb[i][k], (i, k)

def test_joint_row_held_out_stays_on_the_joint_site():
    """Two sites: A = CG rotated 16 deg about CA (the joint fit absorbs it), B = CG with 0.6 Å
    internal noise (it cannot). Joint picks A; the bb-only placement lies nearer B. A joint row
    reporting B's 0.6 Å would pair held_out with a q_site_idx it was not measured at."""
    rng = np.random.default_rng(3)
    noise = rng.normal(size=(4, 3)); noise -= noise.mean(0)
    noise /= np.sqrt((noise ** 2).sum(1).mean())
    site_a = (Q_CG - Q_BB[1]) @ _rot([0.2, -0.7, 0.4], 16.0) + Q_BB[1]
    site_b = Q_CG + 0.6 * noise
    cache = hfc._BoundedBucketCache()
    cache.put(("/lib", "F", 1, "ASP"), {k: (v[:1] if isinstance(v, np.ndarray) and v.ndim else v)
                                        for k, v in _bucket().items()})
    sites = [[(s.astype(np.float32), s.mean(0), np.zeros(4, bool), tuple(range(4)), atoms)]
             for s, atoms in ((site_a, (0, 1, 2, 3)), (site_b, (10, 11, 12, 13)))]
    run = lambda adm: hfc._score_one_bsr_combo(
        "q.pdb", None, "F", "F", sites, (["ASP"], BSR, ["ASP"], [Q_BB]), "/lib", cache,
        {}, 1.0, None, "", adm)[0]
    joint, bb = run("joint"), run("bb")
    placed = LIB_CG @ R_Q + T_Q
    assert (joint["q_site_idx"], bb["q_site_idx"]) == (0, 1), "fixture: sites must diverge"
    assert_discriminates(lambda v: abs(v - float(joint["held_out_cg_rmsd"])) < 1e-3,
                         [_rmsd(placed, site_a)], [_rmsd(placed, site_b)],
                         "joint held_out measured at the joint site")
    assert float(bb["held_out_cg_rmsd"]) == pytest.approx(_rmsd(placed, site_b), abs=2e-4)

def test_records_match_result_fields_exactly(hits):
    for rows in hits:
        for rec in rows.values():
            assert set(rec) == set(RESULT_FIELDS), sorted(set(rec) ^ set(RESULT_FIELDS))
            csv.DictWriter(io.StringIO(), fieldnames=RESULT_FIELDS, delimiter="\t").writerow(rec)

def test_joint_and_backbone_only_fits_provably_diverge():
    """When the CG is displaced relative to the backbone the joint fit trades backbone error
    for CG error, so 'held out' is a real distinction."""
    X, Y = np.concatenate((VDG_BB[1], VDG_CG[1])), np.concatenate((Q_BB, Q_CG))
    R_j, t_j, _ = kabsch(X, Y[None])
    R_b, t_b, _ = kabsch(VDG_BB[1], Q_BB[None])
    assert np.linalg.norm(R_j[0] - R_b[0]) > 1e-3
    assert _rmsd(VDG_CG[1] @ R_j[0] + t_j[0], Q_CG) < _rmsd(VDG_CG[1] @ R_b[0] + t_b[0], Q_CG) - 1e-3

def _run_compact(match_mode):
    cache = hfc._BoundedBucketCache()
    cache.put(("/lib", "F", 1, vdg_npz.make_aa_bucket(["ASP"])), _bucket())
    site = [(Q_CG[p].astype(np.float32), Q_CG.mean(0), np.array([0, 1, 1, 0], bool), tuple(p),
             (5, 6, 7, 8)) for p in (SWAP, list(range(4)))]
    (rec,) = hfc._score_one_bsr_combo(
        "q.pdb", None, "F", "F", [site], (["ASP"], BSR, ["ASP"], [Q_BB]), "/lib", cache,
        {}, None, None, "", match_mode, True)
    return rec

def _joint_cg_residual(i):
    """In-sample CG residual of the best joint bb+CG fit over the query labelings."""
    fits = []
    for p in (SWAP, list(range(4))):
        R, t = _svd_fit(np.concatenate((VDG_BB[i], VDG_CG[i])), np.concatenate((Q_BB, Q_CG[p])))
        fits.append((_rmsd(np.concatenate((VDG_BB[i], VDG_CG[i])) @ R + t, np.concatenate((Q_BB, Q_CG[p]))),
                     _rmsd(VDG_CG[i] @ R + t, Q_CG[p])))
    return min(fits)[1]

@pytest.mark.parametrize("mode", ["joint", "bb"])
def test_compact_records_carry_the_dict_records_values(hits, mode):
    rows, rec = dict(zip(("joint", "bb"), hits))[mode], _run_compact(mode)
    assert rec["match_mode"] == mode and sorted(rec["vdg_index"].tolist()) == sorted(rows)
    for j, i in enumerate(rec["vdg_index"].tolist()):
        for k in ("bb_rmsd", "held_out_cg_rmsd", "cg_bb_dist", "vdg_rmsd"):
            assert float(rec[k][j]) == pytest.approx(float(rows[i][k]), abs=1e-4), (i, k)
        for k in ("vdg_cluster_num_parents", "vdg_cluster_id"):
            assert int(rec[k][j]) == int(rows[i][k]), (i, k)
        assert rec["nr_parent_biounit"][j] == f"p{i}" and rec["charge_sign"][j] == rows[i]["charge_sign"]
        assert rec["placed_cg_element"][j].tolist() == ["C", "O", "N", "C"]
        R, t = _svd_fit(VDG_BB[i], Q_BB)  # placed_cg is the bb-only placement, not the joint one
        np.testing.assert_allclose(rec["placed_cg"][j], VDG_CG[i] @ R + t, atol=1e-3)
        # q_lig_atom_idx names the ligand atoms the held-out error was measured against
        assert _rmsd(rec["placed_cg"][j], Q_CG[rec["q_lig_atom_idx"][j]]) == pytest.approx(
            float(rec["held_out_cg_rmsd"][j]), abs=2e-4), i
        if mode == "joint":
            assert float(rec["in_sample_cg_rmsd"][j]) == pytest.approx(_joint_cg_residual(i), abs=2e-4), i

def _compact(frag, s, b, h, n, mode="bb", n_cg=4):
    return dict(frag=frag, subset_size=s, match_mode=mode, aa_bucket="ASP", charge_sign=np.array(["neut"]),
                vdg_index=np.array([7]), bb_rmsd=np.array([b], np.float32), vdg_rmsd=np.array([0.5], np.float32),
                held_out_cg_rmsd=np.array([h], np.float32), cg_bb_dist=np.array([1.0], np.float32),
                vdg_cluster_num_parents=np.array([n]), nr_parent_biounit=np.array(["1abc_1"]), lig_instance="A_401",
                placed_cg=np.arange(n_cg * 3, dtype=np.float32).reshape(1, n_cg, 3),
                q_lig_atom_idx=np.arange(n_cg, dtype=np.int32)[None], placed_cg_element=np.full((1, n_cg), "C"))

def _tree(pts=((0., 0., 100.),)):
    from scipy.spatial import cKDTree
    return cKDTree(np.asarray(pts, np.float32))

def test_benchmark_selection_rules_pick_by_the_right_key():
    from ligand_vdgs.tools.benchmark_pose_recovery import _hit_arrays, _select
    rows = [(1, 0.10, 9.0, 5), (2, 0.10, 3.0, 50), (2, 0.40, 0.2, 80), (1, 0.05, 7.0, 80)]
    picks = _select(_hit_arrays([_compact("F", *r) for r in rows], _tree())["F"])
    # oracle: min held_out (2); top_bb: min bb (3); top_support: max support, tie -> lower bb (3)
    assert picks == dict(oracle=2, top_bb=3, top_support=3), picks
    rows[3] = (1, 0.5, 7.0, 80)  # now top_bb tie at 0.10 between 0 and 1 -> more support (1)
    assert _select(_hit_arrays([_compact("F", *r) for r in rows], _tree())["F"]) == dict(oracle=2, top_bb=1, top_support=2)

def test_hits_npz_matches_the_pinned_contract(tmp_path):
    from ligand_vdgs.tools.benchmark_pose_recovery import NPZ_FIELDS, _hit_arrays, write_hits_npz
    recs = [_compact("F4", 2, 0.1, 0.3, 9, "bb", 4), _compact("F5", 1, 0.2, 0.4, 3, "joint", 5)]
    recs[1]["placed_cg_element"] = np.array([["C", "O", "N", "C", "S"]])
    write_hits_npz(tmp_path / "hits.npz", _hit_arrays(recs, _tree([[0., 1., 2.5], [50., 50., 50.]])))
    z = np.load(tmp_path / "hits.npz")
    assert tuple(z.files) == NPZ_FIELDS
    assert (z["subset_size"].dtype, z["nr_idx"].dtype, z["cluster_num_parents"].dtype) == (np.int8, np.int32, np.int32)
    assert all(z[k].dtype == np.float32 for k in ("bb_rmsd", "vdg_rmsd", "held_out_cg_rmsd", "cg_bb_dist", "placed_cg",
                                                  "min_bb_dist"))
    assert all(z[k].dtype.kind == "U" for k in ("frag", "match_mode", "charge_sign", "aa_bucket", "nr_parent_biounit",
                                                "lig_instance", "placed_cg_element"))
    assert z["q_lig_atom_idx"].dtype == np.int32 and z["q_lig_atom_idx"].shape == (2, 5)
    i4, i5 = list(z["frag"]).index("F4"), list(z["frag"]).index("F5")
    assert np.isnan(z["vdg_rmsd"][i4]) and z["vdg_rmsd"][i5] == np.float32(0.5), "vdg_rmsd NaN only in bb mode"
    assert z["placed_cg"].shape == (2, 5, 3) and np.isnan(z["placed_cg"][i4, 4]).all()
    np.testing.assert_array_equal(z["placed_cg"][i4, :4], recs[0]["placed_cg"][0])
    np.testing.assert_array_equal(z["placed_cg"][i5], recs[1]["placed_cg"][0])
    assert z["q_lig_atom_idx"][i4].tolist() == [0, 1, 2, 3, -1] and z["q_lig_atom_idx"][i5].tolist() == [0, 1, 2, 3, 4]
    assert z["placed_cg_element"][i4].tolist() == ["C"] * 4 + [""] and z["placed_cg_element"][i5].tolist()[-1] == "S"
    # placed atoms are arange(3 n_cg) rows: (0,1,2) is 0.5 Å from the receptor point (0,1,2.5)
    assert z["min_bb_dist"][i4] == pytest.approx(0.5) and z["min_bb_dist"][i5] == pytest.approx(0.5)
    write_hits_npz(tmp_path / "empty.npz", {})
    z = np.load(tmp_path / "empty.npz")
    assert tuple(z.files) == NPZ_FIELDS and z["placed_cg"].shape == (0, 0, 3) and z["bb_rmsd"].dtype == np.float32
    assert z["q_lig_atom_idx"].shape == (0, 0) and z["q_lig_atom_idx"].dtype == np.int32
    assert z["placed_cg_element"].shape == (0, 0) and z["min_bb_dist"].dtype == np.float32
