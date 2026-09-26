"""Overlap graph (functions/fragment_overlap.py): topology and placed-geometry edges agree on real nesting."""
import numpy as np

from ligand_vdgs.functions import fragment_overlap as fo
from tests.vacuity import assert_discriminates

def test_topology_edges_containment_and_jaccard():
    e = fo.topology_edges([{0, 1, 2}, {0, 1, 2, 3, 4}, {4, 5}, {7, 8}])
    got = {(i, j): (s, round(jac, 3), a, b) for i, j, s, jac, a, b in
           zip(e["i"].tolist(), e["j"].tolist(), e["n_shared"].tolist(), e["jaccard"].tolist(), e["i_in_j"].tolist(), e["j_in_i"].tolist())}
    assert got == {(0, 1): (3, 0.6, True, False), (1, 2): (1, 0.167, False, False)}
    assert fo.maximal_mask(4, e).tolist() == [False, True, True, True]

LIG = np.array([[0, 0, 0], [1.5, 0, 0], [2.2, 1.3, 0], [3.7, 1.3, 0], [4.4, 2.6, 0], [5.9, 2.6, 0]], float)
ELEM = np.array(["C", "C", "N", "C", "O", "C"])

def _site(idx, jitter=0.0, rng=np.random.default_rng(0), n=5):
    """Placed site over ligand atoms `idx`, padded to n with '' / NaN."""
    x = np.full((n, 3), np.nan)
    x[:len(idx)] = LIG[idx] + rng.normal(scale=jitter, size=(len(idx), 3))
    return x, np.array(list(ELEM[idx]) + [""] * (n - len(idx)))

def test_placed_edges_recover_topology_from_geometry_alone():
    idx = [[0, 1, 2], [0, 1, 2, 3, 4], [4, 5], [1, 2, 3]]
    sites = [_site(i, 0.1) for i in idx]
    pe = fo.placed_edges(np.stack([x for x, _ in sites]), np.stack([e for _, e in sites]))
    te = fo.topology_edges(idx)
    key = lambda e: {(i, j): (s, a, b) for i, j, s, a, b in zip(*(e[k].tolist() for k in ("i", "j", "n_shared", "i_in_j", "j_in_i")))}
    assert key(pe) == key(te)
    assert fo.maximal_mask(4, pe).tolist() == [False, True, True, False]

def test_placed_shared_is_element_matched_one_to_one_within_tol():
    a, ea = _site([0, 1, 2])
    b, eb = _site([0, 1, 2])
    eb = np.where(eb == "N", "O", eb)  # same position, other element: not shared
    assert fo.placed_shared(a[None], ea[None], b[None], eb[None]).tolist() == [2]
    # one b atom 0.75 A from two a atoms pairs with exactly one of them, in both directions
    mid, emid = np.array([[[0.75, 0, 0]] + [[np.nan] * 3] * 4]), np.array([["C", "", "", "", ""]])
    assert fo.placed_shared(a[None], ea[None], mid, emid).tolist() == fo.placed_shared(mid, emid, a[None], ea[None]).tolist() == [1]
    shifted = a + np.array([0.6, 0, 0])
    assert_discriminates(lambda t: fo.placed_shared(a[None], ea[None], shifted[None], ea[None], t)[0] == 3,
                         [0.7, 1.0], [0.5], "tol decides whether a 0.6 A shift is still the same atom")

def test_cluster_observations_join_nr_and_mem_rows():
    d = dict(cluster_id=np.array([4, 9]), mem_cluster_id=np.array([9, 4, 9]))
    for f, nr, mem in (("parent_biounit", ["1abc_1", "2abc_1"], ["3abc_1", "4abc_1", "5abc_1"]), ("cg_seg", ["", ""], ["", "", ""]),
                       ("cg_chain", ["A", "B"], ["A", "A", "C"]), ("cg_resnum", [301, 302], [1, 2, 3])):
        d["nr_" + f], d["mem_" + f] = np.array(nr), np.array(mem)
    for f, v in (("scrr_seg", ""), ("scrr_chain", "A"), ("scrr_resnum", 10)):
        d["nr_" + f], d["mem_" + f] = np.full((2, 1), v), np.full((3, 1), v)
    d["nr_cg_names"], d["mem_cg_names"] = np.array([["C1", "C2"]] * 2), np.array([["C1", "C2"]] * 3)
    obs = fo.cluster_observations(d, 1)
    assert [k[0] for k, _ in obs] == ["2abc_1", "3abc_1", "5abc_1"]
    assert obs[0] == (("2abc_1", "", "B", 302, (("", "A", 10),)), frozenset({"C1", "C2"}))
