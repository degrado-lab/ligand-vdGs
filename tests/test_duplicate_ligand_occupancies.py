"""Guards for the vdG-miner submodule's CG-occupancy reader."""
import os
import sys
import tempfile
import unittest

import numpy as np
import prody as pr

PROGRAMS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "external", "vdG-miner",
    "vdg_miner", "programs")
if PROGRAMS_DIR not in sys.path:
    sys.path.insert(0, PROGRAMS_DIR)

from fingerprint_helpers import _resolve_duplicate_ligand_occupancies


def _atomgroup(occupancies):
    ag = pr.AtomGroup("synthetic")
    n_atoms = len(occupancies)
    ag.setCoords(np.arange(n_atoms * 3, dtype=float).reshape(n_atoms, 3))
    ag.setNames([f"C{i}" for i in range(n_atoms)])
    ag.setResnames(["LIG"] * n_atoms)
    ag.setResnums(np.full(n_atoms, 900, dtype=int))
    ag.setChids(["A"] * n_atoms)
    ag.setSegnames([""] * n_atoms)
    ag.setElements(["C"] * n_atoms)
    ag.setOccupancies(np.asarray(occupancies, dtype=float))
    ag.setAltlocs([" "] * n_atoms)
    return ag


class DuplicateLigandOccupancyTests(unittest.TestCase):
    def test_duplicate_slot_is_skipped_and_logged(self):
        ag = _atomgroup([3.00, 3.00, 2.00])
        with tempfile.NamedTemporaryFile("r+", delete=False) as fh:
            logfile = fh.name
        try:
            self.assertIsNone(
                _resolve_duplicate_ligand_occupancies(ag, "1abc_A_900", logfile))
            with open(logfile) as fh:
                warning = fh.read()
        finally:
            os.unlink(logfile)

        self.assertIn("1abc_A_900", warning)
        self.assertIn("duplicate CG slot occupancy 3.00 (2 atoms)", warning)
        self.assertIn("skipping environment", warning)

    def test_unique_slots_are_preserved_without_pruning(self):
        ag = _atomgroup([3.00, 3.01, 2.00])
        coords_before = ag.getCoords().copy()
        result = _resolve_duplicate_ligand_occupancies(ag, "1abc_A_900")

        self.assertIs(result, ag)
        self.assertEqual(result.numAtoms(), 3)
        np.testing.assert_array_equal(result.getCoords(), coords_before)


if __name__ == "__main__":
    unittest.main()
