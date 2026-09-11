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
    @staticmethod
    def _record(names=("O1", "P", "O2")):
        return {
            "cg_coords": np.arange(9, dtype=np.float32).reshape(3, 3),
            "bbcoords": [np.zeros((3, 3), dtype=np.float32),
                         np.ones((3, 3), dtype=np.float32)],
            "flankseqs": [["A", "-", "vdm", "!", "GLY"], ["B", "C", "vdm", "D", "E"]],
            "flankCAs": [np.zeros((5, 3)), np.ones((5, 3))],
            "biounit": "1abc",
            "scrr": [["", "A", 10, "GLY"], ["", "A", 11, "ALA"]],
            "cg_names": list(names), "cg_elements": ["O", "P", "O"],
            "cg_seg": "", "cg_chain": "L", "cg_resnum": 1, "cg_resname": "LIG",
            "slot_flags": [0, 1],
            "quality": (18.4, 1.0, 22.7, 1.0),
            "bbo": [np.zeros(3, dtype=np.float32), np.full(3, np.nan, dtype=np.float32)],
            # Synthetic until session 3's gate lands. Required fields: the
            # writer refuses a record without them rather than filling in.
            "cg_heavy_degree": [1, 4, 1], "cg_num_h": [0, 0, 1],
            "cg_formal_charge": [-1, 0, 0], "cg_nbr_elems": ["", "OC", ""],
            "perception": 1,
            "buried_area": [12.5, 4.0], "shared_area": [1.5, 0.0],
            "n_atom_pairs": [6, 2],
            "min_heavy_dist": [3.1, 4.8],
        }

    def test_cg_fields_are_copied_without_relabeling(self):
        cols = clus_helpers.records_to_columns([self._record(), self._record()])
        # One record in, one row out: interchangeable slots are handled inside
        # the Stage-1 distance, not by replicating the vdG here.
        self.assertEqual(cols["cgvdmbb"].shape, (2, 3 + 6, 3))
        np.testing.assert_array_equal(cols["cgvdmbb"][1, :3],
                                      np.arange(9).reshape(3, 3))
        np.testing.assert_array_equal(cols["cgvdmbb"][0, 6:], np.ones((3, 3)))
        self.assertEqual(cols["cg_names"][0].tolist(), ["O1", "P", "O2"])
        self.assertEqual(cols["flank_seq"][0].tolist(),
                         ["A", "-", "vdm", "!", "GLY", "B", "C", "vdm", "D", "E"])
        self.assertEqual(cols["flank_ca"].shape, (2, 10, 3))
        self.assertEqual(cols["scrr_resnum"].tolist(), [[10, 11], [10, 11]])
        self.assertEqual(cols["quality"].shape, (2, 4))

    def test_a_name_filling_its_width_is_kept_and_one_over_is_refused(self):
        cols = clus_helpers.records_to_columns([self._record(names=("HO3'", "P", "O2"))])
        self.assertEqual(cols["cg_names"][0, 0], "HO3'")
        with self.assertRaisesRegex(ValueError, "cg_names.*5-char"):
            clus_helpers.records_to_columns([self._record(names=("HO3''", "P", "O2"))])
        long_stem = dict(self._record(), biounit="x" * 33)
        with self.assertRaisesRegex(ValueError, "biounit.*33-char"):
            clus_helpers.records_to_columns([long_stem])

    def test_a_ragged_record_is_rejected(self):
        short = self._record(names=("O1", "P"))
        with self.assertRaises(ValueError):
            clus_helpers.records_to_columns([self._record(), short])


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
    """The nr row lives in nr_*, never also in mem_*."""

    def test_writer_drops_the_nr_row_from_the_members(self):
        import tempfile
        from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as clus

        n = 5
        recs = [SmartSlotOrderSurvivesRecordUnpacking._record() for _ in range(n)]
        for i, rec in enumerate(recs):
            rec["biounit"] = f"1a{i:02d}" if i < 4 else "1a00_2"   # 4 entries
        cols = clus_helpers.records_to_columns(recs)
        members = np.arange(n, dtype=np.int32)
        subgroups = [clus.Subgroup(1, 1, 2, members, 0.25)]
        with tempfile.TemporaryDirectory() as tmp:
            clus._write_bucket_npz(tmp, 2, ("GLY", "ALA"), cols, subgroups, "/db")
            with np.load(os.path.join(tmp, "nr_vdgs", "2", "GLY_ALA.npz")) as z:
                self.assertEqual(z["cluster_size"].tolist(), [n])
                self.assertEqual(z["cluster_num_parents"].tolist(), [4])
                self.assertEqual(z["nr_parent_biounit"].tolist(), ["1a02"])
                mem = z["mem_parent_biounit"].tolist()
                self.assertEqual(len(mem), n - 1)
                self.assertNotIn("1a02", mem)
                self.assertEqual(z["mem_cluster_id"].tolist(), [1] * (n - 1))
                self.assertEqual(z["nr_vdm_bb_coords"].shape, (1, 2, 3, 3))
                self.assertEqual(z["cg_elements"].tolist(), ["O", "P", "O"])
                self.assertEqual(str(z["parent_pdb_dir"]), "/db")
                self.assertAlmostEqual(float(z["nr_cg_max_b"][0]), 18.4, places=5)


if __name__ == "__main__":
    unittest.main()
