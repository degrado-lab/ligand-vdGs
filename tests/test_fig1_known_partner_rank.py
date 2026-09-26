"""fig1_known_partner_rank: nested exclusion and percentile-rank math.

Falsifiers: letting AMIDE_SUB (nested in AMIDE, common.fragment_keys_nested) into K_A moves
ESTER's pct from 0.5 to 0.667; a non-finite s(A,B) is treated as ranked instead of raising;
a tie is dropped by a strict '>' instead of counted by '>='.
"""
import unittest
import os
import tempfile

import numpy as np

from ligand_vdgs.identify_bioisosteres.fig1_known_partner_rank import known_partner_percentile, load_profiles

AMIDE = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[N;!R;D2][C;!R;!H0]'
AMIDE_SUB = '[C;!R;!H0][C;!R;H0]([N;!R;D2])=[O;!R;D1]'
ESTER = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D2][C;!R;!H0]'
ARYL_CL = '[c;!H0][c;!H0][c;H0][Cl;!R;D1]'
LABELS = [AMIDE, ESTER, AMIDE_SUB, ARYL_CL]

def _sim(values):
    sim = np.eye(len(LABELS))
    for (i, j), v in values.items():
        sim[i, j] = sim[j, i] = v
    return sim

class TestKnownPartnerPercentile(unittest.TestCase):
    def test_excludes_nested_fragment_from_K_A(self):
        pct, n_ka, s = known_partner_percentile(
            _sim({(0, 1): 0.5, (0, 2): 0.9, (0, 3): 0.1}), LABELS, AMIDE, ESTER)
        self.assertEqual(n_ka, 2)
        self.assertAlmostEqual(pct, 0.5)
        self.assertAlmostEqual(s, 0.5)

    def test_raises_when_partner_similarity_is_nonfinite(self):
        self.assertRaises(AssertionError, known_partner_percentile,
                          _sim({(0, 1): float('nan'), (0, 2): 0.9, (0, 3): 0.1}), LABELS, AMIDE, ESTER)

    def test_tie_counted_by_ge_not_strict_gt(self):
        pct, n_ka, _ = known_partner_percentile(
            _sim({(0, 1): 0.5, (0, 3): 0.5, (0, 2): 0.9}), LABELS, AMIDE, ESTER)
        self.assertEqual(n_ka, 2)
        self.assertAlmostEqual(pct, 1.0)

    def test_figure_profiles_require_pooled_backbone_mode(self):
        with tempfile.TemporaryDirectory() as root:
            paths = [os.path.join(root, f'{key}_single_aa_freq.npz') for key in ('A', 'B')]
            for path in paths:
                np.savez(path, enrichments=np.arange(5), aa_labels=np.array(list('ABCDE')),
                         counts=np.ones(5), total_count=5, bb_mode='pooled')
            self.assertEqual(len(load_profiles(root)), 2)
            np.savez(paths[1], enrichments=np.arange(5), aa_labels=np.array(list('ABCDE')),
                     counts=np.ones(5), total_count=5, bb_mode='per-residue')
            self.assertRaisesRegex(ValueError, 'expected only pooled', load_profiles, root)

if __name__ == '__main__':
    unittest.main()
