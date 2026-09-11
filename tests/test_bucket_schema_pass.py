"""The columns and provenance added in the rebuild's single schema pass.

Three properties, each with a failure mode that is silent rather than loud:
per-observation annotations must reach *both* row sets (counting observations
needs mem_ rows), contact strength must never be invented, and a bucket from an
older writer must be refused rather than read as if its numbers meant the same
thing.
"""
import json
import os
import tempfile
import unittest

import numpy as np

from ligand_vdgs.functions import clus_helpers, vdg_npz_utils
from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as clus


def _record(**overrides):
    rec = {
        "cg_coords": np.arange(9, dtype=np.float32).reshape(3, 3),
        "bbcoords": [np.zeros((3, 3), dtype=np.float32),
                     np.ones((3, 3), dtype=np.float32)],
        "flankseqs": [["A", "-", "vdm", "!", "GLY"], ["B", "C", "vdm", "D", "E"]],
        "flankCAs": [np.zeros((5, 3)), np.ones((5, 3))],
        "biounit": "1abc",
        "scrr": [["", "A", 10, "GLY"], ["", "A", 11, "ALA"]],
        "cg_names": ["O1", "P", "O2"], "cg_elements": ["O", "P", "O"],
        "cg_seg": "", "cg_chain": "L", "cg_resnum": 1, "cg_resname": "LIG",
        "slot_flags": [0, 1],
        "quality": (18.4, 1.0, 22.7, 1.0),
        "bbo": [np.zeros(3, dtype=np.float32), np.full(3, np.nan, dtype=np.float32)],
        # Synthetic until session 3's SASA gate lands.
        "cg_heavy_degree": [1, 4, 1], "cg_num_h": [0, 0, 1],
        "cg_formal_charge": [-1, 0, 0], "cg_nbr_elems": ["", "OC", ""],
        "perception": vdg_npz_utils_perception(),
        "buried_area": [12.5, 4.0], "shared_area": [1.5, 0.0],
        "n_atom_pairs": [6, 2],
        "min_heavy_dist": [3.1, 5.8],
    }
    rec.update(overrides)
    return rec


def vdg_npz_utils_perception():
    from ligand_vdgs.functions import ligand_perception
    return ligand_perception.PERCEPTION_OPENBABEL


def _write(tmp, recs, nr_idx=0):
    cols = clus_helpers.records_to_columns(recs)
    members = np.arange(len(recs), dtype=np.int32)
    subgroups = [clus.Subgroup(1, 1, nr_idx, members, 0.25)]
    clus._write_bucket_npz(tmp, 2, ("GLY", "ALA"), cols, subgroups, "/db")
    return os.path.join(tmp, "nr_vdgs", "2", "GLY_ALA.npz")


class AnnotationColumnsReachBothRowSets(unittest.TestCase):
    def test_every_new_column_is_written_for_nr_and_mem(self):
        recs = [_record(), _record(cg_num_h=[2, 2, 2], perception=0)]
        with tempfile.TemporaryDirectory() as tmp:
            with np.load(_write(tmp, recs)) as z:
                for key in ("cg_heavy_degree", "cg_num_h", "cg_formal_charge",
                            "cg_nbr_elems", "perception", "vdm_buried_area",
                            "vdm_shared_area", "vdm_n_atom_pairs",
                            "vdm_min_heavy_dist"):
                    self.assertIn(f"nr_{key}", z.files, key)
                    self.assertIn(f"mem_{key}", z.files, key)
                # Counting observations means nr + mem; a mem_ row that lost its
                # annotation would make any read-time pooling wrong by exactly
                # the members, which is most of the data.
                self.assertEqual(z["nr_cg_num_h"].tolist(), [[0, 0, 1]])
                self.assertEqual(z["mem_cg_num_h"].tolist(), [[2, 2, 2]])
                self.assertEqual(z["nr_perception"].tolist(), [1])
                self.assertEqual(z["mem_perception"].tolist(), [0])
                self.assertEqual(z["nr_cg_heavy_degree"].dtype, np.int8)
                self.assertEqual(z["nr_cg_nbr_elems"].dtype, np.uint32)
                # Multiplicity and identity survive the packing; the P atom's
                # out-of-match neighbours were one O and one C.
                self.assertEqual(
                    vdg_npz_utils.decode_nbr_elems(z["nr_cg_nbr_elems"][0, 1]),
                    ("C", "O"))
                self.assertEqual(z["nr_vdm_n_atom_pairs"].dtype, np.int16)
                self.assertEqual(z["nr_vdm_buried_area"].dtype, np.float32)

    def test_contact_strength_is_per_slot_not_per_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            with np.load(_write(tmp, [_record(), _record()])) as z:
                self.assertEqual(z["nr_vdm_buried_area"].shape, (1, 2))
                np.testing.assert_allclose(z["nr_vdm_buried_area"][0],
                                           [12.5, 4.0], rtol=1e-6)
                # 5.8 A: past the deleted 4.5 A writer guard and inside the
                # SASA gate's reach. A surviving guard would drop this record.
                np.testing.assert_allclose(z["nr_vdm_min_heavy_dist"][0],
                                           [3.1, 5.8], rtol=1e-6)

    def test_a_missing_field_is_an_error_not_a_fill(self):
        for field in ("buried_area", "shared_area", "n_atom_pairs",
                      "min_heavy_dist", "cg_heavy_degree", "perception"):
            rec = _record()
            del rec[field]
            with self.assertRaises(KeyError, msg=field):
                clus_helpers.records_to_columns([rec])

    def test_heavy_degree_zero_is_refused_at_write(self):
        # No atom of a connected >=4-atom fragment has zero heavy neighbours, so
        # a 0 is a value the perception never computed. -1 is the legal
        # "unreadable"; 0 would be indistinguishable from a measurement.
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                _write(tmp, [_record(cg_heavy_degree=[1, 0, 1])])
        with tempfile.TemporaryDirectory() as tmp:
            self.assertTrue(os.path.isfile(
                _write(tmp, [_record(cg_heavy_degree=[-1, -1, -1])])))


class ContactStrengthIsRequired(unittest.TestCase):
    ENV = [["1abc", "", "A", 900, 1], ["1abc", "", "A", 10], ["1abc", "", "A", 11]]

    def _rec(self, **over):
        rec = {"buried_area": [1.0, 2.0], "shared_area": [0.5, 0.0],
               "n_atom_pairs": [3, 4], "min_heavy_dist": [3.0, 4.0]}
        rec.update(over)
        return rec

    def test_values_are_keyed_by_residue_not_by_position(self):
        out = clus._contact_strength_by_residue(self._rec(), self.ENV)
        self.assertEqual(out[("", "A", 11)], (2.0, 0.0, 4, 4.0))

    def test_a_missing_field_raises(self):
        for field in ("buried_area", "shared_area", "n_atom_pairs",
                      "min_heavy_dist"):
            rec = self._rec()
            del rec[field]
            with self.assertRaises(ValueError, msg=field):
                clus._contact_strength_by_residue(rec, self.ENV)

    def test_a_length_mismatch_raises(self):
        # Silently truncating would attach slot 0's strength to slot 1.
        with self.assertRaises(ValueError):
            clus._contact_strength_by_residue(self._rec(buried_area=[1.0]), self.ENV)

    def test_a_non_finite_distance_raises(self):
        # inf is the gate's "never measured" marker and cannot describe a
        # residue it admitted as a member.
        with self.assertRaises(ValueError):
            clus._contact_strength_by_residue(
                self._rec(min_heavy_dist=[3.0, float("inf")]), self.ENV)


class SchemaVersionIsEnforced(unittest.TestCase):
    def test_a_fresh_bucket_carries_the_current_version_and_loads(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, [_record(), _record()])
            with np.load(path) as z:
                block = json.loads(str(z["schema"]))
                self.assertEqual(block["schema_version"],
                                 vdg_npz_utils.BUCKET_SCHEMA_VERSION)
                self.assertEqual(block["parent_pdb_dir"], "/db")
                self.assertIn("build_date", block)
            self.assertIsNotNone(vdg_npz_utils.load_vdg_bucket(
                tmp, "", 2, "GLY_ALA"))

    def test_an_old_bucket_is_refused_not_skipped(self):
        # Refusing matters more than warning: the columns changed meaning, so a
        # stale bucket does not fail to load, it answers with wrong numbers.
        # And the refusal must not be swallowed by CORRUPT_NPZ_ERRORS, which
        # downgrades a bad file to a warning and None.
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, [_record(), _record()])
            with np.load(path) as z:
                arrays = {k: z[k] for k in z.files}
            del arrays["schema"]          # what a pre-rebuild bucket looks like
            with open(path, "wb") as handle:
                np.savez_compressed(handle, **arrays)
            # Both readers, not just the hit-finding one: load_bucket_npz is
            # what load_cluster_members, the bioisostere analysis and
            # h_class_diagnostic go through, so a check in one reader only
            # means a stale library is refused by hit finding and read
            # silently by everything else.
            with self.assertRaises(vdg_npz_utils.BucketSchemaMismatch):
                vdg_npz_utils.load_vdg_bucket(tmp, "", 2, "GLY_ALA")
            with self.assertRaises(vdg_npz_utils.BucketSchemaMismatch):
                vdg_npz_utils.load_bucket_npz(path)
        self.assertNotIsInstance(vdg_npz_utils.BucketSchemaMismatch("x"),
                                 vdg_npz_utils.CORRUPT_NPZ_ERRORS)


if __name__ == "__main__":
    unittest.main()


class NeighbourElementPacking(unittest.TestCase):
    """Multiplicity and out-of-vocabulary elements must both survive."""

    def test_multiplicity_is_kept(self):
        # A carbon with two carbon neighbours outside the match is not the same
        # environment as one with a single carbon; a presence bitmask loses that.
        one = vdg_npz_utils.encode_nbr_elems(["C"])
        two = vdg_npz_utils.encode_nbr_elems(["C", "C"])
        self.assertNotEqual(one, two)
        self.assertEqual(vdg_npz_utils.decode_nbr_elems(two), ("C", "C"))

    def test_order_does_not_matter_but_composition_does(self):
        self.assertEqual(vdg_npz_utils.encode_nbr_elems(["O", "C"]),
                         vdg_npz_utils.encode_nbr_elems(["C", "O"]))
        self.assertNotEqual(vdg_npz_utils.encode_nbr_elems(["C", "O"]),
                            vdg_npz_utils.encode_nbr_elems(["C", "N"]))

    def test_an_unlisted_element_is_recorded_not_dropped(self):
        # Se/metals: "something heavy is attached here" is the load-bearing part.
        code = vdg_npz_utils.encode_nbr_elems(["Se"])
        self.assertEqual(vdg_npz_utils.decode_nbr_elems(code), ("X",))
        self.assertNotEqual(code, vdg_npz_utils.encode_nbr_elems([]))

    def test_empty_is_zero(self):
        self.assertEqual(vdg_npz_utils.encode_nbr_elems([]), 0)
        self.assertEqual(vdg_npz_utils.decode_nbr_elems(0), ())

    def test_halogens_stay_distinct_from_carbon(self):
        # 'Cl' must not be read as a C followed by an l.
        self.assertEqual(vdg_npz_utils.decode_nbr_elems(
            vdg_npz_utils.encode_nbr_elems(["Cl"])), ("Cl",))
        self.assertNotEqual(vdg_npz_utils.encode_nbr_elems(["Cl"]),
                            vdg_npz_utils.encode_nbr_elems(["C"]))

    def test_the_packed_column_fits_a_uint32(self):
        # Four heavy neighbours is the physical maximum; the count field
        # saturates at 7, so no real atom can overflow a slot.
        code = vdg_npz_utils.encode_nbr_elems(["C"] * 4 + ["N"] * 4)
        self.assertLess(code, 2 ** 32)
        self.assertEqual(vdg_npz_utils.decode_nbr_elems(code),
                         ("C",) * 4 + ("N",) * 4)
