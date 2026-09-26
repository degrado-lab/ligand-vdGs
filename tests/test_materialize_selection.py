"""Centroid selection in materialize_vdg_pdbs: rank by support, cap across buckets.

Falsifier (TODO §2 gripe): `-n 10 --max-files 15` over a neg and a neut GLN bucket
wrote neg's 10 then 5 neut, because the cap was spent in charge-sign iteration order
and clusters were ranked by cluster_size rather than cluster_num_parents.
"""
import unittest

from ligand_vdgs.generate_vdgs.materialize_vdg_pdbs import cap_across_buckets, rank_clusters
from tests.vacuity import assert_discriminates

# Supports from the phosphate-ester fragment's nr_vdgs/1/{neg,neut}/GLN.npz top 10.
NEG = [59, 57, 44, 33, 24, 21, 18, 18, 18, 18]
NEUT = [195, 192, 139, 113, 112, 92, 88, 70, 68, 67]

def no_dropped_pick_outranks_a_kept_one(case):
    supports, kept = case
    kept_min = min((s for b, n in zip(supports, kept) for s in b[:n]), default=None)
    dropped_max = max((s for b, n in zip(supports, kept) for s in b[n:]), default=None)
    return kept_min is None or dropped_max is None or dropped_max <= kept_min

class TestRankClusters(unittest.TestCase):
    def test_ranks_by_support_not_size(self):
        # Index 0 is the largest cluster but the second-worst supported.
        # By size it would be [0, 3, 2, 1]; the 9-9 tie keeps stored order, not size.
        self.assertEqual(rank_clusters([5, 9, 9, 2], [100, 10, 30, 50]).tolist(), [1, 2, 0, 3])

    def test_min_size_filters_after_support_ranking(self):
        self.assertEqual(rank_clusters([5, 9, 9, 2], [100, 10, 30, 50], min_cluster_size=30,
                                       top_n=2).tolist(), [2, 0])

class TestCapAcrossBuckets(unittest.TestCase):
    def test_gripe_split_goes_to_best_supported(self):
        kept = cap_across_buckets([NEG, NEUT], max_files=15)
        self.assertEqual(kept, [5, 10])
        # The old iteration-order cap produced [10, 5].
        assert_discriminates(no_dropped_pick_outranks_a_kept_one,
                             accepts=[([NEG, NEUT], kept)], rejects=[([NEG, NEUT], [10, 5])],
                             label='support-ordered cap')

    def test_ties_go_to_earlier_bucket_and_keep_a_prefix(self):
        self.assertEqual(cap_across_buckets([[7, 3], [7, 7], [3]], max_files=3), [1, 2, 0])
        self.assertEqual(cap_across_buckets([[3], [3], [3]], max_files=2), [1, 1, 0])

    def test_no_cap_or_loose_cap_keeps_everything(self):
        self.assertEqual(cap_across_buckets([NEG, NEUT]), [10, 10])
        self.assertEqual(cap_across_buckets([NEG, NEUT], max_files=100), [10, 10])
        self.assertEqual(cap_across_buckets([], max_files=5), [])

if __name__ == '__main__':
    unittest.main()
