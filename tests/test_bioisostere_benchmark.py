"""Falsifiers for pairwise-complete rankings, sign partitions, and parent folds."""
import os
import tempfile
import unittest

import numpy as np

from ligand_vdgs.identify_bioisosteres.bioisostere_benchmark import _fold, _fold_counts, _rank, _sign_counts
from ligand_vdgs.identify_bioisosteres.compare_aa_profiles import build_sim_matrix, fast_single_aa_pearson
from tests.vacuity import assert_discriminates

def _profile(name, labels, values):
    return name, np.array(values, dtype=float), labels, np.ones(len(labels), dtype=int), len(labels)

class TestBioisostereBenchmark(unittest.TestCase):
    def test_fast_pearson_preserves_pairwise_missingness_and_coverage(self):
        profiles = [_profile('A', ['a', 'b', 'c', 'd', 'e'], [2, 1, 3, 5, 4]),
                    _profile('B', ['a', 'b', 'd', 'e'], [5, 2, 3, 1]),
                    _profile('C', ['a', 'b'], [3, 5])]
        slow, _ = build_sim_matrix(profiles, 'pearson', min_shared=3, mode='single', min_coverage=.5)
        fast = fast_single_aa_pearson(profiles, min_shared=3, min_coverage=.5)
        np.testing.assert_allclose(fast, slow, equal_nan=True, atol=1e-12)
        self.assertTrue(np.isnan(fast[0, 2]))
        self.assertTrue(np.isfinite(fast[0, 1]))
        assert_discriminates(lambda x: np.isclose(x, slow[0, 1]), accepts=[fast[0, 1]],
                             rejects=[np.corrcoef(np.array([2, 1, 3, 5, 4]),
                                                  np.array([5, 2, 0, 3, 1]))[0, 1]],
                             label='pairwise missing-category Pearson')

    def test_rank_counts_ties_and_refuses_missing_partner(self):
        pool = [('partner', .5), ('tie', .5), ('better', .8), ('worse', -.1)]
        self.assertEqual(_rank(pool, 'partner'), (3, 4))
        self.assertRaises(AssertionError, _rank, pool, 'absent')
        assert_discriminates(lambda rank: rank == 3, accepts=[_rank(pool, 'partner')[0]],
                             rejects=[2], label='tie-inclusive partner rank')

    def test_sign_counts_do_not_pool_charge_or_x(self):
        with tempfile.TemporaryDirectory() as root:
            for sign, value in (('neut', 7), ('neg', 11)):
                base = os.path.join(root, 'frag', 'nr_vdgs', '1', sign)
                os.makedirs(base)
                np.savez(os.path.join(base, 'ARG.npz'), cluster_num_parents=np.array([value]))
                np.savez(os.path.join(base, 'X.npz'), cluster_num_parents=np.array([100]))
            neutral = _sign_counts(root, 'frag', 'neut')
            self.assertEqual(neutral, {'ARG': 7})
            self.assertEqual(_sign_counts(root, 'frag', 'neg'), {'ARG': 11})
            assert_discriminates(lambda x: x == {'ARG': 7}, accepts=[neutral],
                                 rejects=[{'ARG': 18}], label='neutral-sign isolation')

    def test_nested_fold_uses_members_and_disjoint_parent_entries(self):
        stems = {fold: [] for fold in (0, 1)}
        for i in range(100):
            stem = f'{i:04x}'
            if len(stems[_fold(stem)]) < 2: stems[_fold(stem)].append(stem)
        with tempfile.TemporaryDirectory() as root:
            for key in ('short', 'long'):
                base = os.path.join(root, key, 'nr_vdgs', '1', 'neut')
                os.makedirs(base)
                np.savez(os.path.join(base, 'ARG.npz'), cluster_id=np.array([1, 2]),
                         nr_parent_biounit=np.array([stems[0][0], stems[1][0]]),
                         mem_cluster_id=np.array([1, 2]),
                         mem_parent_biounit=np.array([stems[0][1], stems[1][1]]))
            counts_a, parents_a = _fold_counts(root, 'short', 0)
            counts_b, parents_b = _fold_counts(root, 'long', 1)
            self.assertEqual((counts_a['ARG'], counts_b['ARG']), (2, 2))
            self.assertFalse(parents_a & parents_b)
            assert_discriminates(lambda x: not x, accepts=[parents_a & parents_b],
                                 rejects=[set(stems[0] + stems[1])], label='parent-entry fold separation')

if __name__ == '__main__':
    unittest.main()
