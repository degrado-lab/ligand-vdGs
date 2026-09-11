"""The vdM sidechain is a display atom set, re-derived at write time.

Like the carbonyl O, sidechain atoms are recovered by superposing the parent
residue's own N/CA/C onto the stored triplet -- so nothing about the library or
about any alignment changes: the fit block build_vdg_atomgroup_from_npz returns
is CG + N/CA/C whether or not the sidechain was asked for.
"""
import inspect
import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import prody as pr

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.generate_vdgs import materialize_vdg_pdbs as mat
from ligand_vdgs.functions.vdg_struct_utils import VDM_OCC


# A SER-like residue: backbone triplet, carbonyl, sidechain heavy atoms, an
# attached H (must be dropped) and an OXT (backbone, must be dropped).
NAMES = ["N", "CA", "C", "O", "CB", "OG", "HB2", "OXT"]
ELEMS = ["N", "C", "C", "O", "C", "O", "H", "O"]
COORDS = np.array([[0.0, 0.0, 0.0],     # N
                   [1.46, 0.0, 0.0],    # CA
                   [2.0, 1.4, 0.0],     # C
                   [1.4, 2.4, 0.4],     # O
                   [2.0, -0.8, 1.2],    # CB
                   [3.4, -0.6, 1.4],    # OG
                   [1.8, -1.9, 1.0],    # HB2
                   [3.2, 1.6, -0.4]])   # OXT


def _write_residue(path, names, elems, coords, occupancies=None, resname="SER",
                   resnum=10, chain="A"):
    n = len(names)
    ag = pr.AtomGroup("res")
    ag.setCoords(np.asarray(coords, float))
    ag.setNames(np.array(names))
    ag.setResnames(np.array([resname] * n))
    ag.setResnums(np.array([resnum] * n))
    ag.setChids(np.array([chain] * n))
    ag.setSegnames(np.array([""] * n))
    ag.setElements(np.array(elems))
    ag.setOccupancies(np.ones(n) if occupancies is None else np.asarray(occupancies, float))
    pr.writePDB(path, ag)


def _build(stored_bb, parent_path, sidechain=False, carbonyl=False, resname="SER"):
    return vdg_npz.build_vdg_atomgroup_from_npz(
        cg_coords=np.zeros((2, 3)), cg_names=["C1", "C2"], cg_elements=["C", "C"],
        cg_seg="", cg_chain="L", cg_resnum=1, cg_resname="LIG",
        vdm_bb_coords=stored_bb,
        scrr_seg=[""], scrr_chain=["A"], scrr_resnum=[10], scrr_resname=[resname],
        parent_pdb_path=parent_path, include_backbone_carbonyl=carbonyl,
        include_sidechain=sidechain)


def _moved(coords, theta=0.7, shift=(5.0, -3.0, 2.0)):
    rot = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                    [np.sin(theta), np.cos(theta), 0.0],
                    [0.0, 0.0, 1.0]])
    return coords @ rot.T + np.asarray(shift, float)


class SidechainTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.parent = os.path.join(self.tmp.name, "1abc.pdb")
        _write_residue(self.parent, NAMES, ELEMS, COORDS)

    def test_absent_unless_requested(self):
        ag, _ = _build(COORDS[:3][None], self.parent)
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C"])

    def test_placed_by_superposing_the_parent_backbone(self):
        moved = _moved(COORDS)
        ag, _ = _build(moved[:3][None], self.parent, sidechain=True)
        names = list(ag.getNames())
        self.assertEqual(names[2:], ["N", "CA", "C", "CB", "OG"])  # no O/OXT/H
        coords = ag.getCoords()
        np.testing.assert_allclose(coords[-2], moved[4], atol=1e-3)  # CB
        np.testing.assert_allclose(coords[-1], moved[5], atol=1e-3)  # OG
        self.assertEqual(list(ag.getElements())[-2:], ["C", "O"])
        self.assertEqual(list(ag.getResnames())[-2:], ["SER", "SER"])
        self.assertEqual(list(ag.getOccupancies())[-2:], [VDM_OCC, VDM_OCC])

    def test_alignment_block_is_unchanged_by_the_extras(self):
        """The claim the flag rests on: display atoms never enter a fit."""
        plain = _build(COORDS[:3][None], self.parent)[1]
        for kw in ({"sidechain": True}, {"carbonyl": True},
                   {"sidechain": True, "carbonyl": True}):
            with self.subTest(**kw):
                np.testing.assert_array_equal(
                    _build(COORDS[:3][None], self.parent, **kw)[1], plain)

    def test_carbonyl_and_sidechain_together(self):
        ag, _ = _build(COORDS[:3][None], self.parent, sidechain=True, carbonyl=True)
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C", "O", "CB", "OG"])

    def test_glycine_yields_nothing_and_says_nothing(self):
        gly = os.path.join(self.tmp.name, "2abc.pdb")
        _write_residue(gly, ["N", "CA", "C", "O"], ["N", "C", "C", "O"],
                       COORDS[:4], resname="GLY")
        with mock.patch("builtins.print") as printed:
            ag, _ = _build(COORDS[:3][None], gly, sidechain=True, resname="GLY")
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C"])
        self.assertFalse(printed.called)

    def test_a_stripped_non_glycine_residue_warns(self):
        bare = os.path.join(self.tmp.name, "3abc.pdb")
        _write_residue(bare, ["N", "CA", "C", "O"], ["N", "C", "C", "O"], COORDS[:4])
        with mock.patch("builtins.print") as printed:
            ag, _ = _build(COORDS[:3][None], bare, sidechain=True)
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C"])
        self.assertTrue(any("sidechain" in str(c) for c in printed.call_args_list))

    def test_a_mismatched_parent_is_skipped_not_silently_misplaced(self):
        wrong = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 4.0, 0.0]])
        with mock.patch("builtins.print") as printed:
            ag, _ = _build(wrong[None], self.parent, sidechain=True)
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C"])
        self.assertTrue(any("sidechain" in str(c) for c in printed.call_args_list))

    def test_repeated_atom_names_take_the_higher_occupancy_copy(self):
        """Real case (2y1x's SAH): two atoms share a name at a blank altloc."""
        dup = os.path.join(self.tmp.name, "4abc.pdb")
        names = NAMES + ["OG"]
        elems = ELEMS + ["O"]
        decoy = COORDS[5] + np.array([0.9, 0.0, 0.0])
        coords = np.vstack([COORDS, decoy])
        _write_residue(dup, names, elems, coords, occupancies=[1.0] * 8 + [0.3])
        ag, _ = _build(COORDS[:3][None], dup, sidechain=True)
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C", "CB", "OG"])
        np.testing.assert_allclose(ag.getCoords()[-1], COORDS[5], atol=1e-3)

    def test_blank_element_column_falls_back_to_the_atom_name(self):
        blank = os.path.join(self.tmp.name, "5abc.pdb")
        _write_residue(blank, NAMES, [""] * len(NAMES), COORDS)
        ag, _ = _build(COORDS[:3][None], blank, sidechain=True)
        self.assertEqual(list(ag.getNames())[2:], ["N", "CA", "C", "CB", "OG"])
        self.assertEqual(list(ag.getElements())[-2:], ["C", "O"])

    def test_non_canonical_residue_keeps_its_own_atoms(self):
        """`X` slots hold real residues; the rule is name-based, not a table."""
        mse = os.path.join(self.tmp.name, "6abc.pdb")
        _write_residue(mse, ["N", "CA", "C", "O", "CB", "CG", "SE", "CE"],
                       ["N", "C", "C", "O", "C", "C", "SE", "C"],
                       np.vstack([COORDS[:6], COORDS[6] + 2.0, COORDS[7] + 2.0]),
                       resname="MSE")
        ag, _ = _build(COORDS[:3][None], mse, sidechain=True, resname="MSE")
        self.assertEqual(list(ag.getNames())[2:],
                         ["N", "CA", "C", "CB", "CG", "SE", "CE"])
        self.assertEqual(list(ag.getElements())[-2:], ["SE", "C"])


class SidechainDefaultTests(unittest.TestCase):
    """The CLI writes sidechains unless told not to; the npz builder does not."""

    def _args(self, *extra):
        argv = ["materialize_vdg_pdbs.py", "-c", "x/nr_vdgs", "-o", "out", *extra]
        with mock.patch.object(sys, "argv", argv):
            return mat.parse_args()

    def test_cli_defaults_on_and_can_be_turned_off(self):
        self.assertTrue(self._args().sidechain)
        self.assertFalse(self._args("--no-sidechain").sidechain)

    def test_the_npz_builder_still_defaults_off(self):
        """Both CLIs pass the flag explicitly; the raw builder stays off so an
        importer does not silently acquire a parent database dependency."""
        param = inspect.signature(
            vdg_npz.build_vdg_atomgroup_from_npz).parameters["include_sidechain"]
        self.assertEqual(param.default, False)

    def test_a_missing_mirror_drops_the_extras_instead_of_failing(self):
        """The average user has no parent database; backbone-only vdGs still
        come out of the npz, and the warning says what was lost."""
        vdg_npz._warned_missing_parent_db.clear()  # warn-once is process-wide
        with mock.patch.object(vdg_npz, "require_parent_pdb_dir",
                               side_effect=vdg_npz.ParentPdbDirError("no mirror")), \
             mock.patch("builtins.print") as printed:
            self.assertFalse(vdg_npz.parent_extras_available(["--sidechain"]))
        said = " ".join(str(c) for c in printed.call_args_list)
        self.assertIn("--sidechain", said)
        self.assertIn("no mirror", said)

    def test_the_message_is_printed_once_per_run(self):
        """Both writers check twice (before the walk, then on the first bucket);
        the user should not see the block twice."""
        vdg_npz._warned_missing_parent_db.clear()
        with mock.patch.object(vdg_npz, "require_parent_pdb_dir",
                               side_effect=vdg_npz.ParentPdbDirError("no mirror")), \
             mock.patch("builtins.print") as printed:
            vdg_npz.parent_extras_available(["--sidechain"])
            vdg_npz.parent_extras_available(["--sidechain"])
        self.assertEqual(printed.call_count, 1)

    def test_a_reachable_mirror_keeps_them(self):
        with mock.patch.object(vdg_npz, "require_parent_pdb_dir",
                               return_value="/mirror"), \
             mock.patch("builtins.print") as printed:
            self.assertTrue(vdg_npz.parent_extras_available(["--sidechain", "--carbonyl"]))
        self.assertFalse(printed.called)


if __name__ == "__main__":
    unittest.main()
