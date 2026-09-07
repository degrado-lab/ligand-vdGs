"""The parent biounit is parsed once per vdG, not once per optional extra.

Both extras that build_vdg_atomgroup_from_npz can add -- the non-CG ligand atoms
and the vdM backbone carbonyl O -- are re-derived from the *same* parent PDB.
Parsing dominates on a network filesystem, so the parse happens in the caller and
the struct is shared. These tests pin both halves: the count, and the fact that
both extras still land when the parse is shared.
"""
import os
import tempfile
import unittest
from unittest import mock

import numpy as np
import prody as pr

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions.vdg_struct_utils import NONCG_LIGAND_OCC

# Residue 10 chain A: a backbone with its carbonyl O, the source of the extra O.
RES = np.array([[0.0, 0.0, 0.0],     # N
                [1.46, 0.0, 0.0],    # CA
                [2.0, 1.4, 0.0],     # C
                [1.4, 2.4, 0.4]])    # O
# Ligand 1 chain L: four CG atoms (non-collinear, so the fit is well conditioned)
# plus one atom outside the CG, which is what include_full_ligand must recover.
LIG = np.array([[10.0, 0.0, 0.0],    # C1  (CG)
                [11.4, 0.0, 0.0],    # C2  (CG)
                [11.9, 1.3, 0.0],    # C3  (CG)
                [11.0, 2.3, 0.6],    # C4  (CG)
                [13.3, 1.6, 0.2]])   # O9  (NOT in the CG)
CG_NAMES = ["C1", "C2", "C3", "C4"]


def _parent_pdb(path):
    coords = np.vstack([RES, LIG])
    names = ["N", "CA", "C", "O", "C1", "C2", "C3", "C4", "O9"]
    ag = pr.AtomGroup("parent")
    ag.setCoords(coords)
    ag.setNames(np.array(names))
    ag.setResnames(np.array(["ALA"] * 4 + ["LIG"] * 5))
    ag.setResnums(np.array([10] * 4 + [1] * 5))
    ag.setChids(np.array(["A"] * 4 + ["L"] * 5))
    ag.setSegnames(np.array([""] * 9))
    ag.setElements(np.array(["N", "C", "C", "O", "C", "C", "C", "C", "O"]))
    ag.setOccupancies(np.ones(9))
    pr.writePDB(path, ag)


def _build(parent_path, full_ligand, carbonyl, parent_struct=None):
    # Stored frame == parent frame, so every fit below is the identity and any
    # displaced atom is a real error rather than fit noise.
    return vdg_npz.build_vdg_atomgroup_from_npz(
        cg_coords=LIG[:4], cg_names=CG_NAMES, cg_elements=["C"] * 4,
        cg_seg="", cg_chain="L", cg_resnum=1, cg_resname="LIG",
        vdm_bb_coords=RES[:3][None],
        scrr_seg=[""], scrr_chain=["A"], scrr_resnum=[10], scrr_resname=["ALA"],
        parent_pdb_path=parent_path, include_full_ligand=full_ligand,
        include_backbone_carbonyl=carbonyl, parent_struct=parent_struct)


class SharedParentParseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.parent = os.path.join(self.tmp.name, "1abc.pdb")
        _parent_pdb(self.parent)
        self.addCleanup(self.tmp.cleanup)

    def _parse_count(self, **kwargs):
        real = vdg_npz.pr.parsePDB
        with mock.patch.object(vdg_npz.pr, "parsePDB",
                               side_effect=real) as parsed:
            ag, _ = _build(self.parent, **kwargs)
        return ag, parsed.call_count

    def test_both_extras_share_one_parse(self):
        """The case the shared parse exists for; two parses would be the old bug."""
        ag, n = self._parse_count(full_ligand=True, carbonyl=True)
        self.assertEqual(n, 1)

        names = list(ag.getNames())
        occs = list(ag.getOccupancies())
        # The non-CG ligand atom, at its own occupancy band (below the CG band).
        self.assertIn("O9", names)
        self.assertAlmostEqual(occs[names.index("O9")], NONCG_LIGAND_OCC)
        np.testing.assert_allclose(ag.getCoords()[names.index("O9")], LIG[4],
                                   atol=1e-3)
        # The carbonyl O is last, and is the residue's O, not the ligand's.
        self.assertEqual(names[-1], "O")
        self.assertEqual(ag.getResnames()[-1], "ALA")
        np.testing.assert_allclose(ag.getCoords()[-1], RES[3], atol=1e-3)

    def test_each_extra_alone_still_parses_once(self):
        for full_ligand, carbonyl in ((True, False), (False, True)):
            with self.subTest(full_ligand=full_ligand, carbonyl=carbonyl):
                _, n = self._parse_count(full_ligand=full_ligand, carbonyl=carbonyl)
                self.assertEqual(n, 1)

    def test_no_extras_means_no_parse(self):
        _, n = self._parse_count(full_ligand=False, carbonyl=False)
        self.assertEqual(n, 0)

    def test_a_caller_supplied_struct_is_not_reparsed(self):
        """The member loop in materialize_vdg_pdbs hands in its own parse."""
        struct = pr.parsePDB(self.parent)
        with mock.patch.object(vdg_npz.pr, "parsePDB") as parsed:
            ag, _ = _build(self.parent, full_ligand=True, carbonyl=True,
                           parent_struct=struct)
        self.assertEqual(parsed.call_count, 0)
        self.assertIn("O9", list(ag.getNames()))
        self.assertEqual(list(ag.getNames())[-1], "O")

    def test_an_unparseable_parent_drops_both_extras_without_raising(self):
        missing = os.path.join(self.tmp.name, "nope.pdb")
        with mock.patch("builtins.print"):
            ag, _ = _build(missing, full_ligand=True, carbonyl=True)
        names = list(ag.getNames())
        self.assertNotIn("O9", names)
        self.assertEqual(names, CG_NAMES + ["N", "CA", "C"])


if __name__ == "__main__":
    unittest.main()
