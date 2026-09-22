"""DR-62 `cg_placed_h`: per-CG-atom placed-H status, read geometrically from
the already-open parent atomgroup rather than the CCD template."""
import unittest

import numpy as np
import prody as pr

from ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs import (
    _cg_placed_h_in_slot_order)

def _atomgroup(rows):
    """`rows`: list of (name, element, seg, chain, resnum, coord)."""
    ag = pr.AtomGroup("synthetic")
    ag.setNames([r[0] for r in rows])
    ag.setElements([r[1] for r in rows])
    ag.setSegnames([r[2] for r in rows])
    ag.setChids([r[3] for r in rows])
    ag.setResnums(np.array([r[4] for r in rows], dtype=int))
    ag.setCoords(np.array([r[5] for r in rows], dtype=float))
    return ag

class CgPlacedHTests(unittest.TestCase):
    def test_mixed_placed_h_within_and_outside_cutoff(self):
        # O1 has a placed H 1.0 A away (within the 1.3 A cutoff); P has none;
        # O2 has an H but 2.0 A away (outside the cutoff).
        out = _cg_placed_h_in_slot_order(
            _atomgroup([("O1", "O", "", "A", 10, [0.0, 0.0, 0.0]),
                        ("P", "P", "", "A", 10, [3.0, 0.0, 0.0]),
                        ("O2", "O", "", "A", 10, [6.0, 0.0, 0.0]),
                        ("H1", "H", "", "A", 10, [1.0, 0.0, 0.0]),
                        ("H2", "H", "", "A", 10, [8.0, 0.0, 0.0])]),
            np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0]]),
            "", "A", 10)
        self.assertEqual(out.tolist(), [1, 0, 0])
        self.assertEqual(out.dtype, np.int8)

    def test_h_on_a_different_residue_does_not_count(self):
        # A geometrically close H belongs to a different resnum -- e.g. a
        # crystallographic water or a neighboring copy in the trimmed
        # atomgroup -- and must not be attributed to this CG atom.
        self.assertEqual(_cg_placed_h_in_slot_order(
            _atomgroup([("O1", "O", "", "A", 10, [0.0, 0.0, 0.0]),
                        ("H1", "H", "", "A", 11, [0.2, 0.0, 0.0])]),
            np.array([[0.0, 0.0, 0.0]]), "", "A", 10).tolist(), [0])

    def test_no_hydrogens_anywhere_is_a_real_zero_not_unresolvable(self):
        # An unprotonated source structure carries no H atoms at all; that is
        # a genuine "no placed H" answer (0), not the -1 sentinel.
        self.assertEqual(_cg_placed_h_in_slot_order(
            _atomgroup([("O1", "O", "", "A", 10, [0.0, 0.0, 0.0])]),
            np.array([[0.0, 0.0, 0.0]]), "", "A", 10).tolist(), [0])

if __name__ == "__main__":
    unittest.main()
