"""Exact-kernel M9 reconstruction (tools/evidence_reconstruct.py) on synthetic parents."""
import pickle

import numpy as np
import prody as pr
import pytest

from ligand_vdgs.tools import evidence_falsifiers as ef, evidence_reconstruct as er
from tests.vacuity import assert_discriminates

def _parent(gap_after=None, altloc_b_better=False):
    """Chain A residues 1-7 (N/CA/C) plus ligand LIG 101 (C1, C2); residue 4's CA has altlocs A/B."""
    rows = []
    for k, rn in enumerate(("ALA", "GLY", "LEU", "SER", "THR", "VAL", "ILE")):
        x0 = 3.8 * k + (5.0 if gap_after is not None and k > gap_after else 0.0)
        rows += [(rn, k + 1, "N", x0, " ", 1.0), (rn, k + 1, "CA", x0 + 1.2, " ", 1.0), (rn, k + 1, "C", x0 + 2.4, " ", 1.0)]
    rows = [r for r in rows if not (r[1] == 4 and r[2] == "CA")]
    rows += [("SER", 4, "CA", 12.6, "A", 0.4 if altloc_b_better else 0.6), ("SER", 4, "CA", 99.0, "B", 0.6 if altloc_b_better else 0.4)]
    rows += [("LIG", 101, "C1", 0.0, " ", 1.0), ("LIG", 101, "C2", 1.5, " ", 1.0)]
    ag = pr.AtomGroup("t")
    ag.setCoords(np.array([[r[3], 0.0, 0.0] for r in rows]))
    ag.setNames([r[2] for r in rows]); ag.setResnames([r[0] for r in rows]); ag.setResnums([r[1] for r in rows])
    ag.setChids(["A"] * len(rows)); ag.setSegnames([""] * len(rows)); ag.setAltlocs([r[4] for r in rows])
    ag.setOccupancies([r[5] for r in rows])
    return ag

def test_parent_index_matches_instrument_and_survives_pickling():
    for gap, alt_b in ((None, False), (2, False), (None, True)):
        st = _parent(gap, alt_b)
        px = pickle.loads(pickle.dumps(er.ParentIndex(st).precompute()))
        for rn in range(1, 8):
            assert px.flank("", "A", rn) == ef.flank_seq(st, "", "A", rn)
        ca = px.coords("", "A", 4, ["CA"])[0, 0]
        assert ca == pytest.approx(99.0 if alt_b else 12.6)  # max occupancy wins, then altloc A
        assert px.coords("", "A", 101, ["C1", "C3"]) is None
    assert pickle.loads(pickle.dumps(er.ParentIndex(_parent(2)).precompute())).flank("", "A", 3) == ("ALA", "GLY", "LEU", "-", "-")

def test_slot_orders_follow_bucket_parts():
    assert er.slot_orders(["bb", "bb"]) == [(0, 1), (1, 0)]
    assert er.slot_orders(["ALA", "LEU"]) == [(0, 1)]
    assert er.slot_orders(["SER"]) == [(0,)]

def test_bucket_kernel_counts_environments_and_rescales_unresolved():
    rng = np.random.default_rng(0)
    cg, bb = rng.normal(size=(2, 3)) * 2, rng.normal(size=(3, 3))
    class Fake:  # one parent per biounit; "far" parents place the CG 3 A off in the vdM frame
        def __init__(self, far, flank): self.far, self.fl = far, flank
        def coords(self, seg, chain, resnum, names):
            return bb.copy() if tuple(names) == ("N", "CA", "C") else cg + (np.array([3.0, 0, 0]) if self.far else 0)
        def flank(self, seg, chain, resnum): return self.fl
    parents = {"1aaa_1": Fake(False, ("A",) * 5), "1bbb_1": Fake(False, ("B",) * 5), "1ccc_1": Fake(True, ("C",) * 5), "1ddd_1": None}
    # row 0: aaa (nr) + bbb (same pose, other flank) + ccc (shifted) + ddd (unresolved) -> 2 envs * 4/3; row 1: singleton
    d = dict(cluster_id=np.array([1, 2]), cluster_size=np.array([4, 1]), mem_cluster_id=np.array([1, 1, 1]), aa_bucket_parts=np.array(["SER"]),
             nr_parent_biounit=np.array(["1aaa_1", "1bbb_1"]), mem_parent_biounit=np.array(["1bbb_1", "1ccc_1", "1ddd_1"]))
    for p, n in (("nr_", 2), ("mem_", 3)):
        d |= {p + "cg_seg": np.array([""] * n), p + "cg_chain": np.array(["L"] * n), p + "cg_resnum": np.ones(n, int),
              p + "cg_names": np.array([["C1", "C2"]] * n), p + "scrr_seg": np.array([[""]] * n), p + "scrr_chain": np.array([["A"]] * n),
              p + "scrr_resnum": np.ones((n, 1), int)}
    loads = []
    recs = er.bucket_records([d], "", loader=lambda cache_dir, stem: loads.append(stem) or parents[stem])[0]
    assert sorted(loads) == sorted(parents)  # each parent loaded once
    run = lambda t: er.bucket_kernel(d, recs, [(0, 1)], t, cap=64, n_sample=64, rng=rng)
    n_eff, n_mem, n_unres = run(1.0)
    assert n_eff.tolist() == pytest.approx([2 * 4 / 3, 1.0]) and n_mem.tolist() == [4, 1] and n_unres.tolist() == [1, 0]
    assert_discriminates(lambda t: run(t)[0][0] == pytest.approx(8 / 3), [0.5, 2.0], [3.5], "t_CG merges the shifted member only above 3 A")
