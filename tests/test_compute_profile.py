"""The compute profile is opt-in and must cost nothing when it is off."""
import json
import os
import tempfile
import unittest

from ligand_vdgs.functions.compute_profile import ComputeProfile


class DisabledProfileTests(unittest.TestCase):
    def test_every_method_is_a_no_op(self):
        p = ComputeProfile(enabled=False)
        with p.phase("x"):
            p.add("n", 3)
            p.set("k", 1)
            p.merge({"a": 1}, prefix="s.")
            p.section("sub", {"phases": {}})
        self.assertEqual(p.to_dict()["phases"], {})
        self.assertEqual(p.to_dict()["counters"], {})
        self.assertEqual(p.to_dict()["sections"], {})

    def test_write_is_skipped_and_reports_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "p.json")
            self.assertIsNone(ComputeProfile(enabled=False).write(path))
            self.assertFalse(os.path.exists(path))


class EnabledProfileTests(unittest.TestCase):
    def test_phases_accumulate_across_repeated_names(self):
        p = ComputeProfile(enabled=True)
        for _ in range(3):
            with p.phase("cluster"):
                pass
        self.assertEqual(p.phases["cluster"]["calls"], 3)
        self.assertGreaterEqual(p.phases["cluster"]["wall_s"], 0.0)
        self.assertIn("cpu_s", p.phases["cluster"])

    def test_a_raising_block_is_still_timed(self):
        p = ComputeProfile(enabled=True)
        with self.assertRaises(RuntimeError):
            with p.phase("boom"):
                raise RuntimeError
        self.assertEqual(p.phases["boom"]["calls"], 1)

    def test_counters_add_and_merge_with_a_prefix(self):
        p = ComputeProfile(enabled=True)
        p.add("n", 2)
        p.add("n", 3)
        p.merge({"fp": 10, "edges": 1}, prefix="stage1.")
        p.merge({"fp": 5}, prefix="stage1.")
        self.assertEqual(p.counters["n"], 5)
        self.assertEqual(p.counters["stage1.fp"], 15)
        self.assertEqual(p.counters["stage1.edges"], 1)

    def test_round_trips_through_a_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "nested", "p.json")
            p = ComputeProfile(enabled=True, cg="CC(=O)O")
            p.add("buckets", 4)
            self.assertEqual(p.write(path), path)
            with open(path) as handle:
                on_disk = json.load(handle)
            self.assertEqual(on_disk["meta"]["cg"], "CC(=O)O")
            self.assertEqual(on_disk["counters"]["buckets"], 4)
            self.assertEqual(ComputeProfile.load(path), on_disk)

    def test_loading_a_missing_or_corrupt_file_is_not_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(ComputeProfile.load(os.path.join(tmp, "absent.json")), {})
            bad = os.path.join(tmp, "bad.json")
            with open(bad, "w") as handle:
                handle.write("{not json")
            self.assertEqual(ComputeProfile.load(bad), {})
            self.assertEqual(ComputeProfile.load(None), {})


if __name__ == "__main__":
    unittest.main()
