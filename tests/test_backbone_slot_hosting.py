"""The query-side guard that replaced the bbGLY/bbPRO label split.

A vdG stores only N/CA/C per vdM, so a backbone-labeled slot says nothing about
the sidechain the query residue carries. Rather than partitioning the library by
donor residue, hit finding asks the query residue directly whether it can host the
geometry -- against its own real atoms, which is strictly more informative than a
label describing the *library* residue.
"""
import unittest

import numpy as np
import prody as pr

from ligand_vdgs.score_poses.hit_finder_core import (
    BB_SLOT_SIDECHAIN_CLASH, PRO_NH_DONOR_CUTOFF,
    backbone_slot_blockers, backbone_slots_can_host)
from ligand_vdgs.functions.vdg_struct_utils import BB_LABEL


def _residue(resname, atoms):
    """One residue as an AtomGroup. `atoms` is [(name, element, (x,y,z)), ...]."""
    ag = pr.AtomGroup("q")
    ag.setCoords(np.array([a[2] for a in atoms], dtype=float))
    ag.setNames([a[0] for a in atoms])
    ag.setElements([a[1] for a in atoms])
    ag.setResnames([resname] * len(atoms))
    ag.setResnums(np.array([5] * len(atoms), dtype=int))
    ag.setChids(["A"] * len(atoms))
    ag.setSegnames([""] * len(atoms))
    ag.setOccupancies(np.ones(len(atoms)))
    ag.setAltlocs([" "] * len(atoms))
    return ag


# Backbone shared by every fixture; the sidechain is what varies.
_BB = [("N", "N", (0.0, 0.0, 0.0)),
       ("CA", "C", (1.46, 0.0, 0.0)),
       ("C", "C", (2.0, 1.42, 0.0)),
       ("O", "O", (1.3, 2.4, 0.0))]

_COMBO = [("", "A", 5)]
_LABELS = [BB_LABEL]


def _can_host(resname, sidechain, cg_coords, cg_acceptor_mask=None):
    struct = _residue(resname, _BB + sidechain)
    blockers = backbone_slot_blockers(struct, _COMBO, _LABELS, [resname])
    cg = np.asarray(cg_coords, dtype=np.float32)
    if cg_acceptor_mask is None:
        cg_acceptor_mask = np.zeros(len(cg), dtype=bool)
    return backbone_slots_can_host(blockers, cg, cg_acceptor_mask)


class SidechainClashTests(unittest.TestCase):
    # CB sits near (2.0, -0.8, 1.2) for this backbone; the CG is placed relative
    # to that rather than to the backbone, since CB is what does the blocking.
    _CB = ("CB", "C", (2.0, -0.8, 1.2))

    def test_glycine_hosts_a_cg_where_a_cb_would_be(self):
        # No CB, so nothing blocks -- this is the geometry bbGLY made unreachable
        # for every other residue, and which the clash test now admits case by case.
        self.assertTrue(_can_host("GLY", [], [(2.0, -0.8, 1.2)]))

    def test_residue_with_cb_rejects_a_cg_on_top_of_it(self):
        self.assertFalse(_can_host("ALA", [self._CB], [(2.0, -0.8, 1.2)]))

    def test_residue_with_cb_accepts_a_cg_clear_of_it(self):
        # Same residue, CG moved beyond the clash cutoff: admitted. This is the
        # 54% of glycine-derived geometries the label split discarded wholesale.
        far = (2.0, -0.8, 1.2 + BB_SLOT_SIDECHAIN_CLASH + 0.6)
        self.assertTrue(_can_host("ALA", [self._CB], [far]))

    def test_any_cg_atom_clashing_is_enough(self):
        # The CG is rigid, so one overlapping atom disqualifies the whole placement.
        far = (2.0, -0.8, 1.2 + BB_SLOT_SIDECHAIN_CLASH + 0.6)
        self.assertFalse(_can_host("ALA", [self._CB], [far, (2.0, -0.8, 1.2)]))

    def test_sidechain_labeled_slot_is_not_screened(self):
        # Under a sidechain label the contact IS the claim, so the sidechain must
        # not be treated as an obstacle -- otherwise every sidechain vdG self-rejects.
        struct = _residue("ALA", _BB + [self._CB])
        self.assertIsNone(backbone_slot_blockers(struct, _COMBO, ["ALA"], ["ALA"]))


class ProlineDonorTests(unittest.TestCase):
    _RING = [("CB", "C", (2.2, -1.0, 1.4)),
             ("CG", "C", (1.4, -2.2, 1.6)),
             ("CD", "C", (0.2, -1.4, 1.0))]

    # 3.0 A from N, so the donor test fires, but >3.4 A from every sidechain atom
    # of both fixtures, so the clash test does not. Isolating the two matters:
    # placed naively the CG trips the clash test first and these tests pass for
    # the wrong reason.
    _NEAR_N = [(-1.8, 1.8, 1.6)]

    def test_proline_rejects_a_geometry_needing_the_backbone_nh(self):
        # An acceptor in H-bond range of N: proline has no amide hydrogen to
        # donate, so it cannot reproduce this interaction.
        self.assertFalse(_can_host("PRO", self._RING, self._NEAR_N,
                                   cg_acceptor_mask=np.array([True])))

    def test_proline_accepts_the_same_geometry_when_the_atom_is_not_an_acceptor(self):
        # Identical coordinates, carbon instead of N/O/S: no donation implied.
        self.assertTrue(_can_host("PRO", self._RING, self._NEAR_N,
                                  cg_acceptor_mask=np.array([False])))

    def test_non_proline_may_donate_its_backbone_nh(self):
        self.assertTrue(_can_host("ALA", [("CB", "C", (2.0, -0.8, 1.2))],
                                  self._NEAR_N,
                                  cg_acceptor_mask=np.array([True])))

    def test_proline_accepts_a_geometry_clear_of_both_ring_and_nh(self):
        # ~81% of the backbone pool, against the 2.75% bbPRO exposed.
        self.assertTrue(_can_host("PRO", self._RING, [(6.0, 4.0, -3.0)],
                                  cg_acceptor_mask=np.array([True])))


if __name__ == "__main__":
    unittest.main()
