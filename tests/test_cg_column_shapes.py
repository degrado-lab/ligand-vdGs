"""The CG's residue is one value per record and its elements one row per bucket.

Both were once stored per CG atom per record. The shapes are the invariant, so
they are asserted here rather than left to the writer's comments.
"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ligand_vdgs.functions import clus_helpers
from ligand_vdgs.functions.vdg_npz_utils import build_vdg_atomgroup_from_npz
from ligand_vdgs.functions import vdg_struct_utils as su


class _Atom:
    def __init__(self, name, element, seg, chain, resnum, resname, occ):
        self._v = (name, element, seg, chain, resnum, resname, occ)

    def getCoords(self): return np.zeros(3, dtype=np.float32)
    def getName(self): return self._v[0]
    def getElement(self): return self._v[1]
    def getSegname(self): return self._v[2]
    def getChid(self): return self._v[3]
    def getResnum(self): return self._v[4]
    def getResname(self): return self._v[5]
    def getOccupancy(self): return self._v[6]


class GetCgAtomsCollapsesResidue(unittest.TestCase):
    def _run(self, atoms):
        # Bypass ProDy selection; exercise only the collapse/validation logic.
        real_select, real_sort = su.select_cg_atoms, su.sort_cg_atoms_by_slot
        su.select_cg_atoms = lambda obj: object()
        su.sort_cg_atoms_by_slot = lambda cg: atoms
        try:
            return su.get_cg_atoms(None, "test.pdb")
        finally:
            su.select_cg_atoms, su.sort_cg_atoms_by_slot = real_select, real_sort

    def test_one_residue_returns_scalars(self):
        atoms = [_Atom(f"O{i}", "O", "", "A", 4, "SO4", 1.0) for i in range(3)]
        result = self._run(atoms)
        self.assertIsNotNone(result)
        coords, names, elements, seg, chain, resnum, resname = result
        self.assertEqual(coords.shape, (3, 3))
        self.assertEqual(len(names), 3)  # per atom
        self.assertEqual(len(elements), 3)  # per atom
        self.assertEqual((seg, chain, resnum, resname), ("", "A", 4, "SO4"))

    def test_straddling_two_residues_is_refused(self):
        atoms = [_Atom("O1", "O", "", "A", 4, "SO4", 1.0),
                 _Atom("O2", "O", "", "A", 5, "SO4", 1.0)]
        self.assertIsNone(self._run(atoms))


class SmartSlotOrderSurvivesRecordUnpacking(unittest.TestCase):
    def test_cg_fields_are_copied_without_relabeling(self):
        coords = np.arange(9, dtype=np.float32).reshape(3, 3)
        bb = [np.zeros((3, 3), dtype=np.float32),
              np.ones((3, 3), dtype=np.float32)]
        record = [
            coords, bb, [["A"], ["B"]], [np.zeros(3), np.ones(3)], "1abc.pdb",
            [["", "A", 10, "GLY"], ["", "A", 11, "ALA"]],
            ["O1", "P", "O2"], ["O", "P", "O"], "", "L", 1, "LIG", [0, 1],
            18.4, 1.0, 22.7, 1.0,   # cg_max_b, cg_min_occ, vdm_max_b, vdm_min_occ
        ]

        out = clus_helpers.unpack_vdg_records([record, record])
        out_coords, out_names = out[0], out[6]

        # One record in, one record out: interchangeable slots are handled inside
        # the Stage-1 distance, not by replicating the vdG here.
        self.assertEqual(len(out_coords), 2)
        for cg_coords, cg_names in zip(out_coords, out_names):
            np.testing.assert_array_equal(cg_coords, coords)
            self.assertEqual(cg_names, ["O1", "P", "O2"])
        self.assertIsNot(out_coords[0], out_coords[1])
        self.assertIsNot(out_names[0], out_names[1])

    def test_a_wrong_width_record_is_rejected(self):
        with self.assertRaises(ValueError):
            clus_helpers.unpack_vdg_records([[1, 2, 3]])


class BuildAtomGroupTakesScalarResidue(unittest.TestCase):
    def test_every_cg_atom_gets_the_one_residue(self):
        n_cg = 4
        ag, _ = build_vdg_atomgroup_from_npz(
            cg_coords=np.arange(n_cg * 3, dtype=float).reshape(n_cg, 3),
            cg_names=[f"C{i}" for i in range(n_cg)], cg_elements=["C"] * n_cg,
            cg_seg="", cg_chain="A", cg_resnum=7, cg_resname="LIG",
            vdm_bb_coords=np.zeros((1, 3, 3)),
            scrr_seg=[""], scrr_chain=["B"], scrr_resnum=[10], scrr_resname=["ASP"])
        self.assertEqual(list(ag.getResnums()[:n_cg]), [7] * n_cg)
        self.assertEqual(list(ag.getResnames()[:n_cg]), ["LIG"] * n_cg)
        self.assertEqual(list(ag.getChids()[:n_cg]), ["A"] * n_cg)


class CentroidIsStoredExactlyOnce(unittest.TestCase):
    """The medoid lives in nr_*, never also in mem_*."""

    def test_extract_members_drops_the_medoid(self):
        from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as clus

        n = 5
        args = ([f"pdb{i}" for i in range(n)],                     # pdbpaths
                [[("", "A", 10 + i, "ASP")] for i in range(n)],    # scrr
                [["O1", "O2"] for _ in range(n)],                  # cg_names
                [""] * n, ["A"] * n, list(range(n)), ["SO4"] * n,  # cg residue
                [(10.0, 1.0, 20.0, 1.0)] * n)                      # quality
        members = clus._extract_members(*args, member_idxs=[0, 1, 2, 3, 4],
                                        nr_idx=2)

        self.assertEqual(len(members), n - 1)
        self.assertNotIn("pdb2", [m["pdbpath"] for m in members])
        # cluster_size counts the medoid, so the rows are one short of it.
        self.assertEqual(len(members) + 1, n)
        self.assertNotIn("is_centroid", members[0])
        self.assertEqual(members[0]["quality"], (10.0, 1.0, 20.0, 1.0))


if __name__ == "__main__":
    unittest.main()
