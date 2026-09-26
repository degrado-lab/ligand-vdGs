"""fig3_geometry_overlay: pre-registered bucket pick from the mmd2 TSV.

Falsifiers: the pooled 'bb' row wins despite being excluded by spec; a bucket is picked by
lowest mmd2 (tuning on the result) instead of highest min(N_A, N_B) support.
"""
import csv
import os
import tempfile
import unittest

from ligand_vdgs.identify_bioisosteres.fig3_geometry_overlay import AMIDE, ARYL_CL, ESTER, pick_bucket, stats_for_bucket

ROWS = [
    dict(frag_A=AMIDE, frag_B=ESTER, subset_size='1', aa_bucket='bb', N_A='9999', N_B='9999', mmd2='9.0'),
    dict(frag_A=AMIDE, frag_B=ESTER, subset_size='1', aa_bucket='LEU', N_A='1170', N_B='1098', mmd2='0.05'),
    dict(frag_A=AMIDE, frag_B=ESTER, subset_size='1', aa_bucket='PHE', N_A='200', N_B='200', mmd2='-9.0'),
    dict(frag_A=AMIDE, frag_B=ARYL_CL, subset_size='1', aa_bucket='LEU', N_A='1170', N_B='500', mmd2='0.10')]

class TestPickBucket(unittest.TestCase):
    def _tsv(self, rows):
        f = tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False)
        w = csv.DictWriter(f, fieldnames=['frag_A', 'frag_B', 'subset_size', 'aa_bucket', 'N_A', 'N_B', 'mmd2'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(rows)
        f.close()
        return f.name

    def test_excludes_bb_and_picks_highest_support_not_lowest_mmd2(self):
        path = self._tsv(ROWS)
        try:
            self.assertEqual(pick_bucket(path), 'LEU')
        finally:
            os.unlink(path)

    def test_stats_for_bucket_labels_all_three_pairs(self):
        path = self._tsv(ROWS)
        try:
            stats = stats_for_bucket(path, 'LEU')
        finally:
            os.unlink(path)
        self.assertAlmostEqual(stats['mmd2_amide_ester'], 0.05)
        self.assertAlmostEqual(stats['mmd2_amide_arylCl'], 0.10)
        self.assertEqual(set(stats), {'mmd2_amide_ester', 'mmd2_amide_arylCl'})

if __name__ == '__main__':
    unittest.main()
