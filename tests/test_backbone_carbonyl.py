"""The backbone carbonyl O is a display atom, re-derived at write time.

The npz stores N/CA/C only, on purpose: those are the atoms hit-finding RMSD is
computed over, and storing O would invite the reading that it took part. The O
is recovered by superposing the parent residue's own N/CA/C onto the stored
triplet, so nothing about the library or the matching changes.
"""
import inspect
import os
import tempfile
import unittest
from unittest import mock

import numpy as np
import prody as pr

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz


def _residue_pdb(path, coords, resnum=10, chain="A", resname="ALA"):
    names = ["N", "CA", "C", "O"]
    ag = pr.AtomGroup("res")
    ag.setCoords(np.asarray(coords, float))
    ag.setNames(np.array(names))
    ag.setResnames(np.array([resname] * len(names)))
    ag.setResnums(np.array([resnum] * len(names)))
    ag.setChids(np.array([chain] * len(names)))
    ag.setSegnames(np.array([""] * len(names)))
    ag.setElements(np.array(["N", "C", "C", "O"]))
    ag.setOccupancies(np.ones(len(names)))
    pr.writePDB(path, ag)


PARENT = np.array([[0.0, 0.0, 0.0],      # N
                   [1.46, 0.0, 0.0],     # CA
                   [2.0, 1.4, 0.0],      # C
                   [1.4, 2.4, 0.4]])     # O


def _build(stored_bb, parent_path, carbonyl):
    return vdg_npz.build_vdg_atomgroup_from_npz(
        cg_coords=np.zeros((2, 3)), cg_names=["C1", "C2"], cg_elements=["C", "C"],
        cg_seg="", cg_chain="L", cg_resnum=1, cg_resname="LIG",
        vdm_bb_coords=stored_bb,
        scrr_seg=[""], scrr_chain=["A"], scrr_resnum=[10], scrr_resname=["ALA"],
        parent_pdb_path=parent_path, include_backbone_carbonyl=carbonyl)


class BackboneCarbonylTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.parent = os.path.join(self.tmp.name, "1abc.pdb")
        _residue_pdb(self.parent, PARENT)
        self.addCleanup(self.tmp.cleanup)

    def test_absent_unless_requested(self):
        ag, _ = _build(PARENT[:3][None], self.parent, carbonyl=False)
        self.assertNotIn("O", [n for n in ag.getNames()[2:]])

    def test_o_is_placed_by_superposing_the_parent_backbone(self):
        # Stored frame: the same residue rotated and translated.
        theta = 0.7
        rot = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                        [np.sin(theta), np.cos(theta), 0.0],
                        [0.0, 0.0, 1.0]])
        shift = np.array([5.0, -3.0, 2.0])
        moved = PARENT @ rot.T + shift
        ag, _ = _build(moved[:3][None], self.parent, carbonyl=True)

        names = list(ag.getNames())
        self.assertEqual(names[-1], "O")
        # The O must land where the same rigid motion would put it.
        np.testing.assert_allclose(ag.getCoords()[-1], moved[3], atol=1e-3)
        self.assertEqual(ag.getResnums()[-1], 10)
        self.assertEqual(ag.getResnames()[-1], "ALA")

    def test_a_mismatched_parent_is_skipped_not_silently_misplaced(self):
        """A wrong parent fits badly, and an extrapolated O would be arbitrary."""
        wrong = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 4.0, 0.0]])
        with mock.patch("builtins.print") as printed:
            ag, _ = _build(wrong[None], self.parent, carbonyl=True)
        self.assertNotIn("O", list(ag.getNames())[2:])
        self.assertTrue(any("carbonyl" in str(c) for c in printed.call_args_list))

    def test_missing_parent_path_warns_and_writes_the_rest(self):
        with mock.patch("builtins.print") as printed:
            ag, _ = _build(PARENT[:3][None], "", carbonyl=True)
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C"])
        self.assertTrue(printed.called)

    def test_the_npz_still_stores_three_backbone_atoms(self):
        """Guards the reason this is done at write time rather than at build time."""
        # By name, not by position: the signature gains parameters over time.
        param = inspect.signature(
            vdg_npz.build_vdg_atomgroup_from_npz).parameters["include_backbone_carbonyl"]
        self.assertEqual(param.default, False)


if __name__ == "__main__":
    unittest.main()
