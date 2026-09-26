"""Charge-sign partitioning in hit_finder_core.py: load_vdg_bucket is now per charge-sign
partition (no all-signs mode), and _combo_worker's task-level exception
handling must still abort score_one_model's nprocs>1 path on a real failure
(e.g. BucketSchemaMismatch) instead of silently degrading to a clean-looking
0-hit result.
"""
import os
import tempfile
import unittest

import numpy as np

from ligand_vdgs.functions import clus_helpers, vdg_npz_utils
from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as clus
from ligand_vdgs.functions.vdg_npz_utils import load_vdg_bucket_all_signs
from ligand_vdgs.score_poses.hit_finder_core import _raise_if_combo_task_errors

from tests.test_bucket_schema_pass import _record

def _write(tmp, recs, sign, aa_tuple=("GLY", "ALA")):
    # One single-member cluster per record, so nr row count == len(recs) --
    # a Subgroup spanning multiple records would collapse them into one
    # cluster representative, undercounting nr rows.
    clus._write_bucket_npz(
        tmp, 2, sign, aa_tuple, clus_helpers.records_to_columns(recs),
        [clus.Subgroup(1, 1, i, np.array([i], dtype=np.int32), 0.25)
         for i in range(len(recs))], "/db")

class LoadAllSignsMerges(unittest.TestCase):
    def test_rows_from_every_present_sign_are_concatenated(self):
        # 'neg' gets 1 nr row, 'pos' gets 2 -- a merge that dropped or
        # duplicated a sign's rows would disagree on this count, not just the
        # aggregate presence/absence.
        with tempfile.TemporaryDirectory() as tmp:
            _write(tmp, [_record()], sign="neg")
            _write(tmp, [_record(), _record()], sign="pos")
            bucket = load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA")
            self.assertEqual(bucket["cg"].shape[0], 3)
            self.assertEqual(bucket["cluster_id"].shape[0], 3)
            self.assertEqual(bucket["aa_bucket_parts"], ["GLY", "ALA"])
            self.assertEqual(
                list(zip(bucket["charge_signs"].tolist(), bucket["partition_indices"].tolist())),
                [("pos", 0), ("pos", 1), ("neg", 0)])

    def test_missing_sign_is_not_padded_with_anything(self):
        # Only 'neg' written; 'pos'/'neut'/'unreadable' contribute zero rows,
        # not zero-filled placeholders that would look like real vdGs.
        with tempfile.TemporaryDirectory() as tmp:
            _write(tmp, [_record(), _record()], sign="neg")
            self.assertEqual(
                load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA")["cg"].shape[0], 2)

    def test_no_sign_present_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA"))

    def test_fields_subset_and_per_row_cg_elements(self):
        with tempfile.TemporaryDirectory() as tmp:
            _write(tmp, [_record()], sign="neg")
            _write(tmp, [_record(), _record()], sign="pos")
            full = load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA")
            self.assertEqual(full["cg_elements"].shape, full["cg"].shape[:2])
            self.assertEqual(set(load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA", fields=("cg",))),
                             {"cg", "aa_bucket_parts", "charge_signs", "partition_indices"})

class IterBucketFiles(unittest.TestCase):
    def test_yields_every_sign_file_only_for_a_completed_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            _write(tmp, [_record()], sign="neg")
            _write(tmp, [_record()], sign="pos")
            self.assertEqual(list(vdg_npz_utils.iter_bucket_files(tmp, "")), [], "incomplete job must be skipped")
            with open(os.path.join(tmp, "_log"), "w") as f:
                f.write("Job completed.\n")
            got = [(k, sign, b, os.path.exists(p)) for k, sign, b, p in vdg_npz_utils.iter_bucket_files(tmp, "")]
            self.assertEqual(got, [(2, "pos", "GLY_ALA", True), (2, "neg", "GLY_ALA", True)])
            self.assertEqual(list(vdg_npz_utils.iter_bucket_files(tmp, "", subset=1)), [])

    def test_legacy_v2_bucket_raises_loud_not_swallowed_to_none(self):
        # A flat schema-v2 bucket (no <sign>/ level) must still refuse loud
        # through this merge helper -- silently returning None here would
        # read as "zero matches" for every fragment in a v2 library.
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "nr_vdgs", "2"))
            open(os.path.join(tmp, "nr_vdgs", "2", "GLY_ALA.npz"), "w").close()
            with self.assertRaises(vdg_npz_utils.BucketSchemaMismatch):
                load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA")

class ComboTaskErrorsAbort(unittest.TestCase):
    def test_no_errors_is_a_no_op(self):
        _raise_if_combo_task_errors([], n_tasks=5)  # must not raise

    def test_any_error_raises_with_count_and_traceback(self):
        with self.assertRaises(RuntimeError) as ctx:
            _raise_if_combo_task_errors(["tb-of-first-failure"], n_tasks=5)
        self.assertIn("1/5", str(ctx.exception))
        self.assertIn("tb-of-first-failure", str(ctx.exception))

if __name__ == "__main__":
    unittest.main()
