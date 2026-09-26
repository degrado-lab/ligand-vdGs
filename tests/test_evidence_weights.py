"""M9 N_eff sidecar (functions/evidence_weights.py) and its falsifier geometry (tools/evidence_falsifiers.py)."""
import json

import numpy as np
import pytest

from ligand_vdgs.functions import evidence_weights as ew
from ligand_vdgs.functions.vdg_npz_utils import BUCKET_SCHEMA_VERSION
from ligand_vdgs.tools import evidence_falsifiers as ef
from tests.vacuity import assert_discriminates

def _bucket(path, s1, sizes, nr_bu, mem_cid, mem_bu):
    np.savez(path, schema=np.asarray(json.dumps({"schema_version": BUCKET_SCHEMA_VERSION})),
             first_stage_cluster_id=np.asarray(s1, np.int32), cluster_id=np.arange(1, len(s1) + 1, dtype=np.int32),
             cluster_size=np.asarray(sizes, np.int32), cluster_num_parents=np.asarray(sizes, np.int32),
             nr_parent_biounit=np.asarray(nr_bu), mem_cluster_id=np.asarray(mem_cid, np.int32),
             mem_parent_biounit=np.asarray(mem_bu, dtype="U32"))
    return path

def test_n_eff_counts_local_environments_per_stage1_cluster(tmp_path):
    # Stage-1 cluster 7: 300 redundant copies in one subgroup (row 0) -> 1, however many parents it has.
    # Stage-1 cluster 9: three distinct local environments (rows 1-3) -> 3 each, including the singleton row 3.
    mem_bu = [f"{i:04d}_1" for i in range(299)] + ["2abc_1", "3abc_1"]
    rows = ew.bucket_rows(_bucket(tmp_path / "b.npz", [7, 9, 9, 9], [300, 2, 2, 1],
                                  ["9zzz_1", "1abc_1", "1abc_2", "4abc_1"], [1] * 299 + [2, 3], mem_bu))
    assert rows["n_eff"].tolist() == [1, 3, 3, 3]
    assert rows["cluster_num_parents"].tolist()[0] == 300  # the redundancy M9 replaces
    assert rows["stage1_n_obs"].tolist() == [300, 5, 5, 5]
    # 1abc_1/1abc_2 are one entry: 1abc, 2abc, 3abc, 4abc
    assert rows["stage1_n_entries"].tolist() == [300, 4, 4, 4]

def test_schema_mismatch_is_refused(tmp_path):
    p = tmp_path / "old.npz"
    np.savez(p, first_stage_cluster_id=np.asarray([1], np.int32))
    with pytest.raises(ValueError, match="schema"):
        ew.bucket_rows(p)

def _sidecar(chunks):
    tables = {}
    coded = [ew.encode({k: np.asarray(v) for k, v in c.items()} | {f: np.zeros(len(c["nr_idx"]), np.int32)
                       for f in ew.VALUE_FIELDS if f not in c}, tables) for c in chunks]
    return ew.assemble(coded, tables)

def test_lookup_maps_each_hit_to_its_own_row():
    # two shards encoded against one growing table, as `merge` does; B arrives in the second shard
    sc = _sidecar([dict(frag=["A", "A"], subset_size=[1, 1], charge_sign=["neut", "neut"], aa_bucket=["ALA", "ALA"],
                        nr_idx=[0, 1], n_eff=[5, 2]),
                   dict(frag=["B"], subset_size=[1], charge_sign=["neg"], aa_bucket=["ALA"], nr_idx=[0], n_eff=[9])])
    hits = dict(frag=["B", "A", "A"], subset_size=[1, 1, 1], charge_sign=["neg", "neut", "neut"],
                aa_bucket=["ALA", "ALA", "ALA"], nr_idx=[0, 1, 0])
    assert ew.lookup(sc, hits).tolist() == [9, 2, 5]
    hits["charge_sign"][0] = "neut"  # sign is part of the key: B/neut/0 does not exist
    with pytest.raises(KeyError):
        ew.lookup(sc, hits)

def test_duplicate_and_out_of_range_keys_are_refused():
    row = dict(frag=["A"], subset_size=[1], charge_sign=["neut"], aa_bucket=["ALA"], nr_idx=[0], n_eff=[1])
    with pytest.raises(ValueError, match="duplicate"):
        _sidecar([row, row])
    with pytest.raises(ValueError, match="nr_idx outside"):
        _sidecar([row | dict(nr_idx=[1 << 27])])

def _obs(rng, n_obs, shift=None):
    """n_obs rigidly moved copies of one CG+vdM geometry (bb (v=1, 3, 3)); `shift` moves the CG of the last."""
    bb, cg = rng.normal(size=(1, 3, 3)) * 1.5, rng.normal(size=(4, 3)) * 2
    out = [(cg @ _rot(0.3 * k) + k, bb @ _rot(0.3 * k) + k) for k in range(n_obs)]
    if shift is not None:
        c, b = out[-1]
        out[-1] = (c + shift @ _rot(0.3 * (n_obs - 1)), b)
    return np.stack([c for c, _ in out]), np.stack([b for _, b in out])

WT, MUT = (("ALA", "GLY", "LEU", "SER", "THR"),), (("ALA", "GLY", "LEU", "SER", "VAL"),)

def test_subgroup_n_eff_kernel_invariants():
    rng, ident = np.random.default_rng(1), [(0, 1, 2, 3)]
    cg, bb = _obs(rng, 5)
    assert ew.subgroup_n_eff([WT] * 5, cg, bb, ident, [(0,)], 1.2) == pytest.approx(1.0)  # identical copies
    # same-pose mutant adds ~0: 4 WT + 1 MUT, all within t_CG -> K all ones -> 1
    assert ew.subgroup_n_eff([WT] * 4 + [MUT], cg, bb, ident, [(0,)], 1.2) == pytest.approx(1.0)
    cg_s, bb_s = _obs(rng, 5, shift=np.array([3.0, 0, 0]))
    # shifted mutant is its own environment: 4 WT (1) + 1 MUT (1) -> 2
    assert ew.subgroup_n_eff([WT] * 4 + [MUT], cg_s, bb_s, ident, [(0,)], 1.2) == pytest.approx(2.0)
    # the identical-flank rule merges the same shifted copy: flank identity overrides geometry
    assert ew.subgroup_n_eff([WT] * 5, cg_s, bb_s, ident, [(0,)], 1.2) == pytest.approx(1.0)
    assert_discriminates(lambda t: ew.subgroup_n_eff([WT] * 4 + [MUT], cg_s, bb_s, ident, [(0,)], t) == pytest.approx(2.0),
                         [1.2, 2.5], [3.5], "t_CG separates the shifted mutant only below its shift")

def test_subgroup_n_eff_slot_order_symmetry():
    # two-slot bucket (e.g. ALA_ALA): member b lists its slots swapped; flanks and backbone must be matched
    rng, ident = np.random.default_rng(2), [(0, 1, 2, 3)]
    bb, cg = rng.normal(size=(2, 3, 3)) * 1.5, rng.normal(size=(4, 3)) * 2
    f1, f2 = ("ALA",) * 5, ("GLY",) * 5
    flanks, cgs, bbs = [(f1, f2), (f2, f1)], np.stack([cg, cg]), np.stack([bb, bb[::-1]])
    assert ew.subgroup_n_eff(flanks, cgs, bbs, ident, [(0, 1), (1, 0)], 0.5) == pytest.approx(1.0)
    assert ew.subgroup_n_eff(flanks, cgs, bbs, ident, [(0, 1)], 0.5) == pytest.approx(2.0)

def _rot(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], float)

def test_falsifier_rmsds_are_frame_invariant_and_minimize_over_automorphisms():
    rng = np.random.default_rng(0)
    bb, cg = rng.normal(size=(3, 3)) * 1.5, rng.normal(size=(4, 3)) * 2
    move = lambda x: x @ _rot(0.7) + np.array([5.0, -3.0, 1.0])
    swapped = move(cg)[[1, 0, 2, 3]]  # same geometry, CG atoms 0/1 listed in the other order
    ident, sym = [(0, 1, 2, 3)], [(0, 1, 2, 3), (1, 0, 2, 3)]
    assert ef.joint_rmsd(cg, bb, move(cg), move(bb), ident) < 1e-3
    assert ef.in_frame_cg_rmsd(cg, bb, move(cg), move(bb), ident) < 1e-3
    assert ef.joint_rmsd(cg, bb, swapped, move(bb), sym) < 1e-3
    assert ef.joint_rmsd(cg, bb, swapped, move(bb), ident) > 0.5
    shifted = move(cg + np.array([1.0, 0, 0]))
    assert_discriminates(lambda y: ef.in_frame_cg_rmsd(cg, bb, y, move(bb), ident) < 0.1, [move(cg)], [shifted],
                         "in-frame CG RMSD sees a CG shift relative to the backbone")

def test_t_cg_calibration_is_per_n_p90_of_identical_pairs_with_pooled_fallback():
    four, three = "[C][C][C][O]", "[C][C][O]"
    rows = ([dict(kind="identical_flank", frag=four, cg_rmsd=str(v)) for v in np.linspace(0, 1, 201)]
            + [dict(kind="identical_flank", frag=three, cg_rmsd="5.0")] * 10
            + [dict(kind="mutant", frag=four, cg_rmsd="9.0")] * 500)  # mutants never calibrate
    t = ef.calibrate_t_cg(rows)
    assert t["t_cg"] == {"4": 0.9} and t["n_pairs"] == {"4": 201, "3": 10}
    assert ef.t_cg_threshold(t, three) == t["pooled"] > 0.9  # the small stratum falls back to all pairs

def test_sidecar_bucket_paths_decode_the_key(tmp_path):
    sc = _sidecar([dict(frag=["A", "B", "B"], subset_size=[1, 2, 2], charge_sign=["neut", "neg", "neg"],
                        aa_bucket=["ALA", "bb_bb", "bb_bb"], nr_idx=[3, 0, 7], stage1_n_entries=[1, 30, 30])])
    got = ef.sidecar_bucket_paths(sc, str(tmp_path), sc["stage1_n_entries"] >= 20)
    assert [k for k, _ in got] == [("B", 2, "neg", "bb_bb")] and got[0][1].endswith("bb_bb.npz")

def test_f1_counts_family_subgroups_per_stage1_cluster():
    # Stage-1 1: subgroups 10 (fam a,b + non-fam z), 11 (fam c only), 12 (non-fam only); Stage-1 2: fam d only
    d = dict(first_stage_cluster_id=np.array([1, 1, 1, 2]), cluster_id=np.array([10, 11, 12, 13]),
             nr_parent_biounit=np.array(["1aaa_1", "1ccc_1", "1zzz_1", "1ddd_1"]),
             mem_cluster_id=np.array([10, 10, 11, 12]), mem_parent_biounit=np.array(["1bbb_2", "1zzz_2", "1ccc_2", "1yyy_1"]))
    rows = ef.f1_bucket(d, {"F": {"1aaa", "1bbb", "1ccc", "1ddd"}}, min_entries=2)
    assert rows == [dict(family="F", s1=1, fam_entries=3, n_subgroups=3, fam_subgroups_any=2, fam_subgroups_only=1)]
    assert_discriminates(lambda m: len(ef.f1_bucket(d, {"F": {"1aaa", "1bbb", "1ccc", "1ddd"}}, m)) == 2, [1], [2, 4],
                         "min_entries gates which Stage-1 clusters are reported")

def test_parent_list_csr_roundtrip_and_duplicate_refusal():
    out = ew.assemble_parents(np.array([9, 7, 9, 7], np.int64), np.array(["2abc", "1abc", "1abc", "3abc"]))
    assert ew.s1_parents(out, 7).tolist() == ["1abc", "3abc"] and ew.s1_parents(out, 9).tolist() == ["1abc", "2abc"]
    with pytest.raises(KeyError):
        ew.s1_parents(out, 8)
    with pytest.raises(ValueError, match="duplicate"):
        ew.assemble_parents(np.array([7, 7], np.int64), np.array(["1abc", "1abc"]))

def test_s1_entry_pairs_are_mem_level():
    d = dict(first_stage_cluster_id=np.array([1, 2]), cluster_id=np.array([5, 6]), nr_parent_biounit=np.array(["1abc_1", "1abc_2"]),
             mem_cluster_id=np.array([5, 6, 6]), mem_parent_biounit=np.array(["1abc_3", "2abc_1", "3abc_1"]))
    s1, e = ew.s1_entry_pairs(d)
    assert list(zip(s1.tolist(), e.tolist())) == [(1, "1abc"), (2, "1abc"), (2, "2abc"), (2, "3abc")]

def test_sampled_n_eff_bounds_and_unbiasedness():
    rng, ident, m = np.random.default_rng(3), [(0, 1, 2, 3)], 40
    cg, bb = _obs(rng, m)
    fl = [((f"R{k % 4}",) * 5,) for k in range(m)]  # 4 flank groups; all CGs in the same frame -> K all ones
    assert ew.subgroup_n_eff(fl, cg, bb, ident, [(0,)], 1.2, sample=rng.choice(m, 8, replace=False)) == pytest.approx(1.0)
    cg_far = np.stack([cg[k] + np.array([5.0 * k, 0, 0]) @ _rot(0.3 * k) for k in range(m)])  # each CG 5k A off in its own frame
    distinct = [((f"R{k}",) * 5,) for k in range(m)]
    assert ew.subgroup_n_eff(distinct, cg_far, bb, ident, [(0,)], 1.2, sample=rng.choice(m, 8, replace=False)) == pytest.approx(m)
    # mixture: 10 copies of one environment + 30 distinct; sampled mean matches exact
    mix_fl = [WT] * 10 + distinct[10:]
    exact = ew.subgroup_n_eff(mix_fl, cg_far, bb, ident, [(0,)], 1.2)
    assert exact == pytest.approx(31.0)
    est = [ew.subgroup_n_eff(mix_fl, cg_far, bb, ident, [(0,)], 1.2, sample=rng.choice(m, 8, replace=False)) for _ in range(400)]
    assert abs(np.mean(est) - exact) < 0.03 * exact and min(est) >= 1 and max(est) <= m

def test_pair_kernel_matches_the_f2_instrument():
    rng = np.random.default_rng(4)
    sym = [(0, 1, 2, 3), (1, 0, 2, 3)]
    for _ in range(20):
        bb_a, bb_b = rng.normal(size=(1, 3, 3)), rng.normal(size=(1, 3, 3))
        cg_a, cg_b = rng.normal(size=(4, 3)) * 2, rng.normal(size=(4, 3)) * 2
        want = ef.in_frame_cg_rmsd(cg_a, bb_a.reshape(-1, 3), cg_b, bb_b.reshape(-1, 3), sym)
        got = ew.in_frame_cg_rmsd_pairs(cg_a[None], bb_a[None], cg_b[None], bb_b[None], sym)[0]
        assert got == pytest.approx(want, abs=1e-4)
    # two slots: swapping b's slots is undone only when both orders are allowed
    bb2, cg = rng.normal(size=(2, 3, 3)), rng.normal(size=(4, 3))
    move = lambda x: x @ _rot(1.1) + 2.0
    args = (cg[None], bb2[None], move(cg)[None], move(bb2)[::-1][None], [(0, 1, 2, 3)])
    assert ew.in_frame_cg_rmsd_pairs(*args, [(0, 1), (1, 0)])[0] < 1e-6 < ew.in_frame_cg_rmsd_pairs(*args, [(0, 1)])[0]
