"""
Placement (bb) mode reads neither the query CG nor side chains (coordinator Q1 ruling, 2026-09-25).

Falsifiers:
- the slot gate still reads the query (crystal) CG -> a clashing query CG empties the hit pool, and
  a clean query CG lets a clashing placement through.
- the gate reads real side-chain atoms instead of the virtual CB -> a far side-chain atom blocks.
- the virtual CB has the wrong handedness (sign of the cross term) -> > 2 Å from the real L-CB.
- the Pro N acceptor check reads query-CG elements instead of the library's cg_elements.
- joint mode reads side chains: its slot gate or contact filter uses real atoms, not N/CA/C + vCB.
- min_bb_dist's receptor includes side-chain atoms or a Gly virtual CB.
"""
import inspect

import numpy as np
import prody as pr
import pytest

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions import vdg_struct_utils as su
from ligand_vdgs.functions.vdg_struct_utils import BB_LABEL
from ligand_vdgs.score_poses import hit_finder_core as hfc
from tests.vacuity import assert_discriminates

# Engh-Huber ideal Ala (L): N, CA, C and the real CB.
BB = np.array([[-0.525, 1.363, 0.0], [0.0, 0.0, 0.0], [1.526, 0.0, 0.0]], np.float32)
CB = np.array([-0.529, -0.774, -1.205], np.float32)
CLEAN = np.array([[3, 4, 3], [4, 4, 3], [3.5, 5, 3], [3, 3.5, 4]], np.float32)
CLASH = np.concatenate([CB[None] + [1.0, 0.5, 0.5], CLEAN[1:]]).astype(np.float32)
NEAR_PRO_N = np.concatenate([BB[:1] + [0.0, 3.0, 0.0], CLEAN[1:]]).astype(np.float32)  # 3.0 Å from N, far from CB
BSR = [("", "A", 10)]

def _run(resname, q_cg, vdg_cgs, elements=("C",) * 4, q_acceptor=np.zeros(4, bool), contact_cutoff=None, mode="bb"):
    """Hits of an exact-backbone bucket (placed CG == stored CG); returns hit vdg indices."""
    n = len(vdg_cgs)
    cache = hfc._BoundedBucketCache()
    cache.put(("/lib", "F", 1, vdg_npz.make_aa_bucket([BB_LABEL])), dict(aa_bucket_parts=[BB_LABEL], charge_signs=np.full(n, "neut", "U11"),
                  partition_indices=np.arange(n, dtype=np.int32), cluster_id=np.arange(n, dtype=np.int32),
                  cluster_num_parents=np.ones(n, np.int32), parent_biounit=np.array(["p"] * n),
                  cg_elements=np.broadcast_to(np.array(elements, "U2"), (n, 4)),
                  cg=np.asarray(vdg_cgs, np.float32), bb=np.broadcast_to(BB, (n, 1, 3, 3)).copy()))
    return sorted(int(r["vdg_index"]) for r in hfc._score_one_bsr_combo(
        "q.pdb", None, "F", "F", [[(np.asarray(q_cg, np.float32), q_cg.mean(0), q_acceptor, (0, 1, 2, 3), (0, 1, 2, 3))]],
        ([BB_LABEL], BSR, [resname], [BB]), "/lib", cache, {}, None, contact_cutoff, "", mode))

def test_virtual_cb_is_the_l_cb():
    vcb = su.virtual_cb(BB)
    assert_discriminates(lambda v: np.linalg.norm(v - CB) < 0.1, [vcb],
                         [BB[1] + 0.58273431 * np.cross(BB[1] - BB[0], BB[2] - BB[1])  # mirrored (D) CB
                          + 0.56802827 * (BB[1] - BB[0]) - 0.54067466 * (BB[2] - BB[1])], "vCB is the L-enantiomer CB")
    assert np.linalg.norm(vcb - BB[1]) == pytest.approx(1.52, abs=0.03)

def test_fixture_geometry():
    vcb = su.virtual_cb(BB)
    d = lambda cg, p: np.linalg.norm(cg - p, axis=1).min()
    assert d(CLASH, vcb) < hfc.BB_SLOT_SIDECHAIN_CLASH <= d(CLEAN, vcb)
    assert d(NEAR_PRO_N[:1], BB[0]) < hfc.PRO_NH_DONOR_CUTOFF <= d(NEAR_PRO_N[1:], BB[0])
    assert d(NEAR_PRO_N, vcb) >= hfc.BB_SLOT_SIDECHAIN_CLASH and d(CLEAN, BB[0]) >= hfc.PRO_NH_DONOR_CUTOFF

def test_slot_gate_acts_on_the_placed_cg_not_the_query_cg():
    # vdG 0 places a clean CG, vdG 1 a CG on the virtual CB; the query CG decides nothing.
    for q_cg in (CLEAN, CLASH):
        assert _run("ALA", q_cg, [CLEAN, CLASH]) == [0], "query CG must not gate placement hits"
    assert _run("GLY", CLASH, [CLEAN, CLASH]) == [0, 1], "Gly has no (virtual) CB"

def test_pro_n_check_reads_library_elements():
    assert _run("PRO", CLEAN, [NEAR_PRO_N], elements=("O", "C", "C", "C")) == []
    assert _run("PRO", CLEAN, [NEAR_PRO_N], elements=("C",) * 4, q_acceptor=np.ones(4, bool)) == [0]
    assert _run("ALA", CLEAN, [NEAR_PRO_N], elements=("O", "C", "C", "C")) == [0], "Pro N check is Pro-only"

def test_contact_cutoff_is_refused_in_placement_mode():
    with pytest.raises(ValueError, match="joint-mode only"):
        _run("ALA", CLEAN, [CLEAN], contact_cutoff=3.8)

def test_receptor_is_backbone_plus_virtual_cb():
    atoms = [("ALA", "N", BB[0]), ("ALA", "CA", BB[1]), ("ALA", "C", BB[2]), ("ALA", "O", [2.1, 1.0, 0.0]),
             ("ALA", "CB", [9.0, 9.0, 9.0])]  # displaced real CB: side chains must not enter
    atoms += [("GLY", n, np.asarray(c) + 10.0) for _r, n, c in atoms[:4]]
    ag = pr.AtomGroup("r")
    ag.setCoords(np.array([a[2] for a in atoms], float)); ag.setNames([a[1] for a in atoms])
    ag.setResnames([a[0] for a in atoms]); ag.setResnums([1] * 5 + [2] * 4); ag.setChids(["A"] * 9)
    ag.setElements([a[1][0] for a in atoms])
    np.testing.assert_allclose(su.receptor_bb_vcb_coords(ag), np.concatenate(
        [np.array([a[2] for a in atoms if a[1] != "CB"], np.float32), su.virtual_cb(BB)[None]]), atol=1e-5)

def test_joint_mode_gates_the_query_cg_on_the_virtual_cb():
    assert _run("ALA", CLASH, [CLASH], mode="joint") == [] and _run("GLY", CLASH, [CLASH], mode="joint") == [0]
    assert _run("PRO", NEAR_PRO_N, [NEAR_PRO_N], q_acceptor=np.array([1, 0, 0, 0], bool), mode="joint") == []

def test_joint_contact_filter_uses_backbone_and_virtual_cb_only():
    np.testing.assert_allclose(su.bsr_contact_atoms(BB[None], ["ALA"]), np.concatenate([BB, su.virtual_cb(BB)[None]]))
    np.testing.assert_allclose(su.bsr_contact_atoms(BB[None], ["GLY"]), BB)
    d = np.linalg.norm(CLEAN[:, None] - su.bsr_contact_atoms(BB[None], ["ALA"]), axis=-1).min()
    assert _run("ALA", CLEAN, [CLEAN], contact_cutoff=d - 0.1, mode="joint") == []
    assert _run("ALA", CLEAN, [CLEAN], contact_cutoff=d + 0.1, mode="joint") == [0]

def test_bb_output_is_invariant_to_the_truth_labeling():
    """bb q_atom_indices are the query-CG (truth) argmin. Permuting which labeling carries which truth
    coordinates changes them but not the placements, so nothing may key on them (the deleted
    deduplicate_hits kept 1 placement under one labeling and 2 under the other)."""
    shift, far = CLEAN + [0.0, 0.0, 2.0], CLEAN + [20.0, 0.0, 0.0]
    def run(l2, l3):
        cache = hfc._BoundedBucketCache()
        cache.put(("/lib", "F", 1, vdg_npz.make_aa_bucket([BB_LABEL])), dict(aa_bucket_parts=[BB_LABEL], charge_signs=np.full(2, "neut", "U11"),
                  partition_indices=np.arange(2, dtype=np.int32), cluster_id=np.arange(2, dtype=np.int32),
                  cluster_num_parents=np.ones(2, np.int32), parent_biounit=np.array(["p"] * 2),
                  cg_elements=np.full((2, 4), "C", "U2"), cg=np.stack([CLEAN, shift]).astype(np.float32),
                  bb=np.broadcast_to(BB, (2, 1, 3, 3)).copy()))
        sites = [[(np.asarray(x, np.float32), x.mean(0), np.zeros(4, bool), (0, 1, 2, 3), atoms)]
                 for x, atoms in ((CLEAN, (0, 1, 2, 3)), (l2, (1, 2, 3, 4)), (l3, (5, 6, 7, 8)))]
        return hfc._score_one_bsr_combo("q.pdb", None, "F", "F", sites, ([BB_LABEL], BSR, ["ALA"], [BB]), "/lib", cache,
                                        {}, None, None, "", "bb")
    t1, t2 = run(shift, far), run(far, shift)  # truth coordinates of labelings 2 and 3 swapped
    placements = lambda recs: sorted(tuple(r[k] for k in ("vdg_index", "bb_rmsd", *(f"Rbb{i}{j}" for i in range(3) for j in range(3)),
                                                          "tbb0", "tbb1", "tbb2")) for r in recs)
    assert [r["q_atom_indices"] for r in t1] != [r["q_atom_indices"] for r in t2], "truth must reach q_atom_indices"
    assert placements(t1) == placements(t2) and len(t1) == 2
    assert not hasattr(hfc, "deduplicate_hits") and "deduplicate" not in inspect.signature(hfc.score_one_model).parameters
