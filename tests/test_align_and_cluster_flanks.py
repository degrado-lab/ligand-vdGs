import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster
from ligand_vdgs.functions.vdg_struct_utils import (FLANK_CHAIN_BREAK,
    FLANK_MISSING)


class FlankingChainBreakTests(unittest.TestCase):
    def test_chain_breaks_are_scanned_independently_on_both_sides(self):
        flanks = {
            -2: ["VAL", np.array([-9.0, 0.0, 0.0], dtype=np.float32)],
            -1: ["GLY", np.array([-3.8, 0.0, 0.0], dtype=np.float32)],
             0: ["vdm", np.array([0.0, 0.0, 0.0], dtype=np.float32)],
             1: ["ALA", np.array([5.0, 0.0, 0.0], dtype=np.float32)],
             2: ["SER", np.array([6.0, 0.0, 0.0], dtype=np.float32)],
        }

        result = align_and_cluster._mark_flanking_chain_breaks(flanks, 2)

        self.assertEqual(result[-1][0], "GLY")
        self.assertEqual(result[-2][0], FLANK_CHAIN_BREAK)
        self.assertEqual(result[1][0], FLANK_CHAIN_BREAK)
        self.assertEqual(result[2][0], FLANK_CHAIN_BREAK)
        self.assertTrue(np.isnan(result[-2][1]).all())
        self.assertTrue(np.isnan(result[1][1]).all())
        self.assertTrue(np.isnan(result[2][1]).all())

    def test_unreadable_flank_is_not_relabeled_as_a_chain_break(self):
        # -1 has no usable CA, so it stays FLANK_MISSING: the residue is
        # unreadable, which says nothing about chain continuity. Only -2, whose
        # continuity can no longer be checked through it, becomes a break.
        flanks = {
            -2: ["VAL", np.array([-7.6, 0.0, 0.0], dtype=np.float32)],
            -1: [FLANK_MISSING, np.array([np.nan] * 3, dtype=np.float32)],
             0: ["vdm", np.array([0.0, 0.0, 0.0], dtype=np.float32)],
             1: ["ALA", np.array([3.8, 0.0, 0.0], dtype=np.float32)],
             2: ["SER", np.array([7.6, 0.0, 0.0], dtype=np.float32)],
        }

        result = align_and_cluster._mark_flanking_chain_breaks(flanks, 2)

        self.assertEqual(result[-1][0], FLANK_MISSING)
        self.assertEqual(result[-2][0], FLANK_CHAIN_BREAK)
        # The other side is continuous and untouched.
        self.assertEqual(result[1][0], "ALA")
        self.assertEqual(result[2][0], "SER")

    def test_break_at_the_outermost_flank_overwrites_nothing_further(self):
        # Missing residue at the last flank: there is nothing past it, so the
        # walk must not raise or relabel the readable positions.
        flanks = {
            -1: ["GLY", np.array([-3.8, 0.0, 0.0], dtype=np.float32)],
             0: ["vdm", np.array([0.0, 0.0, 0.0], dtype=np.float32)],
             1: [FLANK_MISSING, np.array([np.nan] * 3, dtype=np.float32)],
        }

        result = align_and_cluster._mark_flanking_chain_breaks(flanks, 1)

        self.assertEqual(result[-1][0], "GLY")
        self.assertEqual(result[1][0], FLANK_MISSING)


if __name__ == "__main__":
    unittest.main()
