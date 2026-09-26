"""fig2_clustermap: NaN-rate gate on row distance, and the enrichment-matrix column mapping.

Falsifiers: a NaN rate under --nan-rate-max is silently left NaN (breaks scipy linkage)
instead of filled to distance 1; a NaN rate over the max is filled instead of raising; a
category absent from the fixed column order leaks into the matrix instead of being dropped.
"""
import unittest

import numpy as np

from ligand_vdgs.identify_bioisosteres.fig2_clustermap import build_enrichment_matrix, row_distance

class TestRowDistance(unittest.TestCase):
    def test_low_nan_rate_is_filled_not_left_nan(self):
        dist = row_distance(np.array([[1.0, 0.5, np.nan], [0.5, 1.0, 0.2], [np.nan, 0.2, 1.0]]),
                            nan_rate_max=0.5)
        self.assertTrue(np.isfinite(dist).all())
        self.assertAlmostEqual(dist[0, 2], 1.0)

    def test_high_nan_rate_raises_instead_of_filling(self):
        self.assertRaises(ValueError, row_distance,
                          np.array([[1.0, np.nan, np.nan], [np.nan, 1.0, np.nan],
                                   [np.nan, np.nan, 1.0]]), 0.1)

class TestEnrichmentMatrix(unittest.TestCase):
    def test_unlisted_category_is_dropped_not_leaked(self):
        mat = build_enrichment_matrix(
            [('fragA', np.array([1.0, 2.0]), ['ASP', 'NOT_A_CATEGORY'], np.array([9, 9]), 9)],
            ['ASP', 'GLY'])
        self.assertEqual(mat.shape, (1, 2))
        self.assertAlmostEqual(mat[0, 0], 1.0)
        self.assertTrue(np.isnan(mat[0, 1]))

if __name__ == '__main__':
    unittest.main()
