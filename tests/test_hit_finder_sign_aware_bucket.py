"""DR-61 fallout in hit_finder_core.py: load_vdg_bucket is now per charge-sign
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
from ligand_vdgs.score_poses.hit_finder_core import (
    _load_vdg_bucket_all_signs, _raise_if_combo_task_errors)

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
            bucket = _load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA")
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
                _load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA")["cg"].shape[0], 2)

    def test_no_sign_present_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(_load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA"))

    def test_legacy_v2_bucket_raises_loud_not_swallowed_to_none(self):
        # A pre-DR-61 flat bucket (no <sign>/ level) must still refuse loud
        # through this merge helper -- silently returning None here would
        # read as "zero matches" for every fragment in a v2 library.
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "nr_vdgs", "2"))
            open(os.path.join(tmp, "nr_vdgs", "2", "GLY_ALA.npz"), "w").close()
            with self.assertRaises(vdg_npz_utils.BucketSchemaMismatch):
                _load_vdg_bucket_all_signs(tmp, "", 2, "GLY_ALA")

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
