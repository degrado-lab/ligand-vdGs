"""The query-side guard that replaced the bbGLY/bbPRO label split.

A vdG stores only N/CA/C per vdM, so a backbone-labeled slot says nothing about the side chain the
query residue carries. Hit finding asks whether the query residue can host the geometry, using only
its backbone and residue identity (no side-chain coordinates in any mode): a virtual CB (none for
Gly) must be clear of the CG, and a Pro N must not face a CG acceptor (Pro has no NH).
"""
import unittest

import numpy as np

from ligand_vdgs.score_poses.hit_finder_core import (
    BB_SLOT_SIDECHAIN_CLASH, backbone_slots_can_host)
from ligand_vdgs.functions.vdg_struct_utils import BB_LABEL, virtual_cb, vcb_slot_blockers

# N, CA, C shared by every fixture; only the residue identity varies.
_BB = np.array([[0.0, 0.0, 0.0], [1.46, 0.0, 0.0], [2.0, 1.42, 0.0]], np.float32)
_VCB = virtual_cb(_BB)
_FAR = _VCB + [0.0, 0.0, -(BB_SLOT_SIDECHAIN_CLASH + 0.6)]
# 3.0 A from N (donor test fires) but > 3.4 A from the virtual CB (clash test does not), so the
# Pro tests cannot pass for the wrong reason.
_NEAR_N = [(-1.8, 1.8, 1.6)]

def _can_host(resname, cg_coords, cg_acceptor_mask=None):
    cg = np.asarray(cg_coords, np.float32)
    return backbone_slots_can_host(vcb_slot_blockers(_BB[None], [BB_LABEL], [resname]), cg,
                                   np.zeros(len(cg), bool) if cg_acceptor_mask is None else cg_acceptor_mask)

class VirtualCbClashTests(unittest.TestCase):
    def test_fixture_geometry(self):
        self.assertAlmostEqual(np.linalg.norm(np.subtract(_NEAR_N[0], _BB[0])), 3.0, places=1)
        self.assertGreater(np.linalg.norm(_NEAR_N[0] - _VCB), BB_SLOT_SIDECHAIN_CLASH)

    def test_glycine_hosts_a_cg_where_a_cb_would_be(self):
        self.assertTrue(_can_host("GLY", [_VCB]))

    def test_residue_with_cb_rejects_a_cg_on_its_virtual_cb(self):
        self.assertFalse(_can_host("ALA", [_VCB]))

    def test_residue_with_cb_accepts_a_cg_clear_of_it(self):
        self.assertTrue(_can_host("ALA", [_FAR]))

    def test_any_cg_atom_clashing_is_enough(self):
        self.assertFalse(_can_host("ALA", [_FAR, _VCB]))

    def test_sidechain_labeled_slot_is_not_screened(self):
        # Under a side-chain label the contact IS the claim; the slot must not block itself.
        self.assertIsNone(vcb_slot_blockers(_BB[None], ["ALA"], ["ALA"]))

class ProlineDonorTests(unittest.TestCase):
    def test_proline_rejects_a_geometry_needing_the_backbone_nh(self):
        self.assertFalse(_can_host("PRO", _NEAR_N, np.array([True])))

    def test_proline_accepts_the_same_geometry_when_the_atom_is_not_an_acceptor(self):
        self.assertTrue(_can_host("PRO", _NEAR_N, np.array([False])))

    def test_non_proline_may_donate_its_backbone_nh(self):
        self.assertTrue(_can_host("ALA", _NEAR_N, np.array([True])))

    def test_proline_accepts_a_geometry_clear_of_both_cb_and_nh(self):
        self.assertTrue(_can_host("PRO", [(6.0, 4.0, -3.0)], np.array([True])))

if __name__ == "__main__":
    unittest.main()
