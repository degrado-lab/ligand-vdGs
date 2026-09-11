"""One module owns how a parent structure's stem, entry and paths relate."""
import os
import tempfile
import unittest

from ligand_vdgs.functions import parent_db


class StemTests(unittest.TestCase):
    def test_stem_strips_every_structure_extension(self):
        for path, stem in (("/db/ab/1abc.pdb", "1abc"), ("1abc_2.pdb.gz", "1abc_2"),
                           ("x/1ABC.cif", "1ABC"), ("1abc.cif.gz", "1abc"),
                           ("1abc", "1abc")):
            self.assertEqual(parent_db.stem_of(path), stem)

    def test_entry_is_shared_by_assemblies_and_plinder_systems(self):
        self.assertEqual(parent_db.entry_of("1abc"), "1abc")
        self.assertEqual(parent_db.entry_of("1abc_2"), "1abc")
        # A PLINDER system ID: double underscores and dotted chain IDs.
        self.assertEqual(parent_db.entry_of("1abc__1__1.A__1.B"), "1abc")
        self.assertEqual(parent_db.stem_of("/db/ab/1abc__1__1.A__1.B.pdb"),
                         "1abc__1__1.A__1.B")

    def test_shard_is_the_lowercased_middle_two(self):
        self.assertEqual(parent_db.shard("1ABC"), "ab")
        self.assertEqual(parent_db.shard("1abc_2"), "ab")

    def test_paths_round_trip_through_the_stem(self):
        path = parent_db.structure_path("/db", "1abc_2")
        self.assertEqual(path, os.path.join("/db", "ab", "1abc_2.pdb"))
        self.assertEqual(parent_db.stem_of(path), "1abc_2")
        self.assertEqual(parent_db.probe_path("/probe", "1abc_2"),
                         os.path.join("/probe", "ab", "1abc_2.probe.gz"))


class MirrorTests(unittest.TestCase):
    def test_iter_and_is_mirror(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(parent_db.is_mirror(tmp))
            for stem in ("2xyz", "1abc", "1abc_2"):
                path = parent_db.structure_path(tmp, stem)
                os.makedirs(os.path.dirname(path), exist_ok=True)
                open(path, "w").write("x")
            # A file in the wrong shard is listed but does not make a mirror.
            os.makedirs(os.path.join(tmp, "zz"))
            open(os.path.join(tmp, "zz", "3def.pdb"), "w").write("x")
            stems = [stem for stem, _ in parent_db.iter_structures(tmp)]
            self.assertEqual(stems, ["1abc", "1abc_2", "2xyz", "3def"])
            self.assertTrue(parent_db.is_mirror(tmp))


if __name__ == "__main__":
    unittest.main()
