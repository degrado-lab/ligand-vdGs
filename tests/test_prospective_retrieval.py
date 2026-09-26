"""Falsifiers for crystal exclusion and reciprocal retrieval ranks."""
import os
import tempfile
import unittest

import numpy as np

from ligand_vdgs.identify_bioisosteres.bioisostere_benchmark import LONG_CARB, MAIN_CARB, MAIN_TET
from ligand_vdgs.identify_bioisosteres.prospective_retrieval import (
    _bucket_counts, _case_targets, _entry, _fold, _rank, _write)
from tests.vacuity import assert_discriminates

class TestProspectiveRetrieval(unittest.TestCase):
    def test_fold_counts_exclude_crystal_entries_and_deduplicate_members(self):
        stems = {fold: [] for fold in (0, 1)}
        for i in range(100):
            stem = f'{i:04x}'
            if len(stems[_fold(stem)]) < 2: stems[_fold(stem)].append(stem)
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, 'ARG.npz')
            np.savez(path, cluster_id=np.array([1, 2]),
                     nr_parent_biounit=np.array([stems[0][0], stems[1][0]]),
                     mem_cluster_id=np.array([1, 1, 2, 2]),
                     mem_parent_biounit=np.array([stems[0][0], stems[0][1], stems[1][0], stems[1][1]]))
            counts, entries = _bucket_counts(path, {stems[0][0]})
            self.assertEqual(counts, [1, 2])
            self.assertEqual(entries, [{stems[0][1]}, set(stems[1])])
            self.assertFalse(entries[0] & entries[1])
            assert_discriminates(lambda x: x == [1, 2], accepts=[counts], rejects=[[3, 3]],
                                 label='exclude crystal entry and deduplicate repeated parent')
            self.assertEqual(_entry(f'{stems[0][1]}_assembly'), stems[0][1])

    def test_tie_inclusive_rank_and_missing_target(self):
        pool = [('target', .8), ('tie', .8), ('higher', .9), ('lower', .2)]
        self.assertEqual(_rank(pool, 'target'), (3, 4, .75, .8))
        self.assertIsNone(_rank(pool, 'absent'))
        assert_discriminates(lambda x: x == 3, accepts=[_rank(pool, 'target')[0]], rejects=[2],
                             label='tie-inclusive reciprocal rank')

    def test_crystal_targets_require_both_full_heads(self):
        good = {'library_occurrence': {'entries': {
            '2qbp': {'head_fragment_names': [MAIN_CARB]},
            '2nt7': {'head_fragment_names': [MAIN_TET]},
            '7g01': {'head_fragment_names': [MAIN_CARB]},
            '5hz5': {'head_fragment_names': [MAIN_TET]}}}}
        targets = _case_targets(good, {MAIN_TET, MAIN_CARB, LONG_CARB})
        self.assertEqual(targets[(MAIN_TET, MAIN_CARB)], {'PTP1B', 'FABP5'})
        bad = {'library_occurrence': {'entries': dict(good['library_occurrence']['entries'])}}
        bad['library_occurrence']['entries']['5hz5'] = {'head_fragment_names': []}
        self.assertRaises(ValueError, _case_targets, bad, {MAIN_TET, MAIN_CARB})
        assert_discriminates(lambda x: x == {'PTP1B', 'FABP5'},
                             accepts=[targets[(MAIN_TET, MAIN_CARB)]], rejects=[{'PTP1B'}],
                             label='both crystal cases represented')

    def test_empty_success_artifact_is_refused(self):
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, 'rows.tsv')
            self.assertRaises(ValueError, _write, path, [])
            self.assertFalse(os.path.exists(path))

if __name__ == '__main__':
    unittest.main()
