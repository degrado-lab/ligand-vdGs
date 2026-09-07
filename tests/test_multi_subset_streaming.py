import json
import os
import pickle
import tempfile
import unittest
from unittest import mock

import numpy as np

from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as pipeline


class MultiSubsetStreamingTests(unittest.TestCase):
    def test_one_environment_builds_both_subset_sizes_once(self):
        environment = [["1abc", "", "A", 10, 1],
                       ["1abc", "", "A", 20],
                       ["1abc", "", "A", 30]]
        cg_result = (
            np.zeros((3, 3), dtype=np.float32),
            ["C1", "O1", "O2"], ["C", "O", "O"],
            "", "A", 10, "LIG")

        with tempfile.TemporaryDirectory() as temp_dir:
            environments_dir = os.path.join(temp_dir, "environments", "ab")
            os.makedirs(environments_dir)
            with open(os.path.join(environments_dir, "1abc.jsonl"), "w") as handle:
                json.dump({"env": environment, "cg_max_b": 15.0,
                           "cg_min_occ": 1.0, "vdm_max_b": 25.0,
                           "vdm_min_occ": 1.0}, handle)
                handle.write("\n")

            worker_root = os.path.join(temp_dir, "worker")
            args = (
                [("ab", "1abc.jsonl")],
                os.path.join(temp_dir, "environments"),
                os.path.join(temp_dir, "pdb"),
                "test_cg", None, [1, 0, 2],
                os.path.join(temp_dir, "log"), 3, ["C", "O", "O"], 2,
                worker_root, (0, 0), (1, 2), pipeline._FLUSH_RECORDS_THRESHOLD)

            def reorder(subset, _vdms, _cg, _atomgroup):
                labels = ["ALA"] if len(subset) == 1 else ["ALA", "SER"]
                n = len(labels)
                return (labels, [np.zeros((3, 3))] * n,
                        ["AAAAA"] * n, [np.zeros((5, 3))] * n,
                        [("", "A", 20 + i, label)
                         for i, label in enumerate(labels)], [0] * n)

            with mock.patch.object(
                    pipeline, "_get_atomgroup_for_env", return_value=object()), \
                    mock.patch.object(
                        pipeline, "_resolve_duplicate_ligand_occupancies",
                        side_effect=lambda atomgroup, _label: atomgroup), \
                    mock.patch.object(pipeline, "get_cg_atoms", return_value=cg_result), \
                    mock.patch.object(
                        pipeline.clust, "get_vdm_res_features",
                        return_value={20: object(), 30: object()}) as get_features, \
                    mock.patch.object(
                        pipeline, "get_vdg_subsets_target_size",
                        side_effect=lambda _indices, size:
                            [(20,)] if size == 1 else [(20, 30)]), \
                    mock.patch.object(
                        pipeline.clust, "reorder_vdg_subset",
                        side_effect=reorder):
                pipeline._stream_one_chunk(args)

            self.assertEqual(get_features.call_count, 1)
            for size, aa_key in ((1, "ALA"), (2, "ALA_SER")):
                bucket_dir = os.path.join(worker_root, str(size))
                paths = [os.path.join(bucket_dir, name)
                         for name in os.listdir(bucket_dir)]
                self.assertEqual(len(paths), 1)
                self.assertEqual(
                    pipeline._aa_key_from_bucket_fname(os.path.basename(paths[0])),
                    aa_key)
                with open(paths[0], "rb") as handle:
                    self.assertEqual(len(pickle.load(handle)), 1)

    def test_cg_element_mismatch_is_skipped_and_counted(self):
        """A CG whose elements disagree with the SMARTS must be dropped at stream
        time, not written out and raised on at npz-write time."""
        environment = [["1abc", "", "A", 10, 1],
                       ["1abc", "", "A", 20],
                       ["1abc", "", "A", 30]]
        # SMARTS multiset is C,O,O; OpenBabel picked an N.
        cg_result = (
            np.zeros((3, 3), dtype=np.float32),
            ["C1", "N1", "O2"], ["C", "N", "O"],
            "", "A", 10, "LIG")

        with tempfile.TemporaryDirectory() as temp_dir:
            environments_dir = os.path.join(temp_dir, "environments", "ab")
            os.makedirs(environments_dir)
            with open(os.path.join(environments_dir, "1abc.jsonl"), "w") as handle:
                json.dump({"env": environment, "cg_max_b": 15.0,
                           "cg_min_occ": 1.0, "vdm_max_b": 25.0,
                           "vdm_min_occ": 1.0}, handle)
                handle.write("\n")

            worker_root = os.path.join(temp_dir, "worker")
            args = (
                [("ab", "1abc.jsonl")],
                os.path.join(temp_dir, "environments"),
                os.path.join(temp_dir, "pdb"),
                "test_cg", None, [1, 0, 2],
                os.path.join(temp_dir, "log"), 3, ["C", "O", "O"], 2,
                worker_root, (0, 0), (1, 2),
                pipeline._FLUSH_RECORDS_THRESHOLD)

            with mock.patch.object(
                    pipeline, "_get_atomgroup_for_env", return_value=object()), \
                    mock.patch.object(
                        pipeline, "_resolve_duplicate_ligand_occupancies",
                        side_effect=lambda atomgroup, _label: atomgroup), \
                    mock.patch.object(pipeline, "get_cg_atoms", return_value=cg_result), \
                    mock.patch.object(
                        pipeline.clust, "get_vdm_res_features",
                        return_value={20: object(), 30: object()}):
                _, skips, _warns = pipeline._stream_one_chunk(args)

            self.assertEqual(skips["cg_elements_mismatch"], 1)
            for size in (1, 2):
                bucket_dir = os.path.join(worker_root, str(size))
                self.assertEqual(os.listdir(bucket_dir), [])

    def test_generic_smarts_slot_acceptance_is_chunk_independent(self):
        """A slot the SMARTS does not pin must not be validated against whatever
        record a worker happened to see first: acceptance would then depend on
        how --num-procs partitions the shards."""
        # '[#7,#8]=[C;!R][O;!R]' -- slot 0 is an OR query, so it is unpinned.
        expected_element_seq = (None, "C", "O")
        # Two shards that differ only at the unpinned slot; both are legal matches.
        cg_by_label = {
            "1abc__A_10_1": (np.zeros((3, 3), dtype=np.float32),
                             ["N1", "C1", "O1"], ["N", "C", "O"],
                             "", "A", 10, "LIG"),
            "2xyz__A_10_1": (np.zeros((3, 3), dtype=np.float32),
                             ["O1", "C1", "O2"], ["O", "C", "O"],
                             "", "A", 10, "LIG"),
        }

        def run(chunks, temp_dir, tag):
            """Stream `chunks` (one call per chunk, as separate workers would) and
            return (records written, skip counts)."""
            total_records, total_skips = 0, 0
            for i, chunk in enumerate(chunks):
                worker_root = os.path.join(temp_dir, f"worker_{tag}_{i}")
                args = (
                    chunk, os.path.join(temp_dir, "environments"),
                    os.path.join(temp_dir, "pdb"), "test_cg", None, [1, 0, 2],
                    os.path.join(temp_dir, "log"), 3, expected_element_seq, 2,
                    worker_root, (0, 0), (1,), pipeline._FLUSH_RECORDS_THRESHOLD)

                def reorder(subset, _vdms, _cg, _atomgroup):
                    return (["ALA"], [np.zeros((3, 3))], ["AAAAA"],
                            [np.zeros((5, 3))], [("", "A", 20, "ALA")], [0])

                with mock.patch.object(
                        pipeline, "_get_atomgroup_for_env", return_value=object()), \
                        mock.patch.object(
                            pipeline, "_resolve_duplicate_ligand_occupancies",
                            side_effect=lambda atomgroup, _label: atomgroup), \
                        mock.patch.object(
                            pipeline, "get_cg_atoms",
                            side_effect=lambda _ag, label: cg_by_label[label]), \
                        mock.patch.object(
                            pipeline.clust, "get_vdm_res_features",
                            return_value={20: object()}), \
                        mock.patch.object(
                            pipeline, "get_vdg_subsets_target_size",
                            side_effect=lambda _indices, size: [(20,)]), \
                        mock.patch.object(
                            pipeline.clust, "reorder_vdg_subset", side_effect=reorder):
                    _, skips, _warns = pipeline._stream_one_chunk(args)

                total_skips += skips["cg_elements_mismatch"]
                bucket_dir = os.path.join(worker_root, "1")
                for name in os.listdir(bucket_dir):
                    with open(os.path.join(bucket_dir, name), "rb") as handle:
                        total_records += len(pickle.load(handle))
            return total_records, total_skips

        with tempfile.TemporaryDirectory() as temp_dir:
            shards = []
            for biounit in ("1abc", "2xyz"):
                subdir = biounit[1:3]
                environments_dir = os.path.join(temp_dir, "environments", subdir)
                os.makedirs(environments_dir, exist_ok=True)
                with open(os.path.join(environments_dir, f"{biounit}.jsonl"),
                          "w") as handle:
                    json.dump({"env": [[biounit, "", "A", 10, 1],
                                       [biounit, "", "A", 20]],
                               "cg_max_b": 15.0, "cg_min_occ": 1.0,
                               "vdm_max_b": 25.0, "vdm_min_occ": 1.0}, handle)
                    handle.write("\n")
                shards.append((subdir, f"{biounit}.jsonl"))

            one_chunk = run([shards], temp_dir, "single")
            two_chunks = run([[shards[0]], [shards[1]]], temp_dir, "split")

        # Both records are kept either way; nothing is rejected on the unpinned slot.
        self.assertEqual(one_chunk, (2, 0))
        self.assertEqual(one_chunk, two_chunks)

    def test_warning_counts_returned_are_per_chunk_deltas(self):
        """_WARN_COUNTS is module state and the executor reuses a worker across
        chunks, so returning the running totals would re-report every earlier
        chunk's warnings and inflate the end-of-streaming line. Also covers the
        keys nothing else counts: vdm_not_in_contact drops one residue and lets
        the environment through, so it never reaches a skips[...] bucket."""
        pipeline._WARN_COUNTS.clear()
        self.addCleanup(pipeline._WARN_COUNTS.clear)

        with tempfile.TemporaryDirectory() as temp_dir:
            shards = []
            for biounit, subdir in (("1abc", "ab"), ("2xyz", "xy")):
                environments_dir = os.path.join(temp_dir, "environments", subdir)
                os.makedirs(environments_dir, exist_ok=True)
                with open(os.path.join(environments_dir, f"{biounit}.jsonl"),
                          "w") as handle:
                    json.dump({"env": [[biounit, "", "A", 10, 1],
                                       [biounit, "", "A", 20]],
                               "cg_max_b": 15.0, "cg_min_occ": 1.0,
                               "vdm_max_b": 25.0, "vdm_min_occ": 1.0}, handle)
                    handle.write("\n")
                shards.append((subdir, f"{biounit}.jsonl"))

            logfile = os.path.join(temp_dir, "log")

            def run(chunk, tag):
                args = (
                    chunk, os.path.join(temp_dir, "environments"),
                    os.path.join(temp_dir, "pdb"), "test_cg", None, [1, 0, 2],
                    logfile, 3, ["C", "O", "O"], 2,
                    os.path.join(temp_dir, f"worker_{tag}"), (0, 0), (1,),
                    pipeline._FLUSH_RECORDS_THRESHOLD)

                # One residue-level warning per environment, then the environment
                # proceeds -- exactly the vdm_not_in_contact shape.
                def fake_atomgroup(*_a, **_kw):
                    pipeline._log_warn_capped(
                        logfile, "vdm_not_in_contact", "[WARNING] dropped.\n")
                    return None

                with mock.patch.object(
                        pipeline, "_get_atomgroup_for_env",
                        side_effect=fake_atomgroup):
                    _, skips, warns = pipeline._stream_one_chunk(args)
                return skips, warns

            skips_a, warns_a = run([shards[0]], "a")
            skips_b, warns_b = run([shards[1]], "b")

        # Each chunk reports only its own warning, not the running total.
        self.assertEqual(warns_a, {"vdm_not_in_contact": 1})
        self.assertEqual(warns_b, {"vdm_not_in_contact": 1})
        # ... and the running total really did advance, so the deltas are not
        # just a counter that never incremented.
        self.assertEqual(pipeline._WARN_COUNTS["vdm_not_in_contact"], 2)
        # The warning key is absent from the skip counters on purpose: summing
        # the two sets would double count the environment-level reasons.
        self.assertNotIn("vdm_not_in_contact", skips_a)
        self.assertNotIn("vdm_not_in_contact", skips_b)


if __name__ == "__main__":
    unittest.main()
