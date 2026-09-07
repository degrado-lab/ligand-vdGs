"""build_flank_lookup_index must never change what get_AA_and_CA_coords answers.

The index is a pure speedup on the flank walk, and the flank CAs it feeds go
into Stage-2 subclustering -- so a single residue it answers differently from
the selection path silently changes which vdGs get deduplicated, with nothing
in the output to signal it. The test that matters is therefore differential:
every resindex, indexed vs. not, on an atomgroup built to contain the cases the
fast path must refuse to decide.
"""

import unittest

import numpy as np
import prody as pr

from ligand_vdgs.functions.vdg_struct_utils import (
    FLANK_MISSING, build_flank_lookup_index, get_AA_and_CA_coords)


# (resname, [(atom name, occupancy, altloc), ...]). Order matters only in that
# each entry becomes one resindex, in this order.
RESIDUES = [
    ("ALA", [("N", 1.0, ""), ("CA", 1.0, ""), ("C", 1.0, ""), ("O", 1.0, "")]),
    ("GLY", [("N", 1.0, ""), ("CA", 1.0, ""), ("C", 1.0, ""), ("O", 1.0, "")]),
    # Non-protein but carrying an atom named CA: the one case where a naive
    # name-only index would hand back a real coordinate for a ligand.
    ("LIG", [("C1", 1.0, ""), ("CA", 1.0, ""), ("C2", 1.0, ""), ("O1", 1.0, "")]),
    ("HOH", [("O", 1.0, "")]),
    # Altloc CAs: resolved by occupancy, then by altloc 'A' among ties. The fast
    # path must not pick either itself.
    ("SER", [("N", 1.0, ""), ("CA", 0.4, "A"), ("CA", 0.6, "B"), ("C", 1.0, "")]),
    ("THR", [("N", 1.0, ""), ("CA", 0.5, "A"), ("CA", 0.5, "B"), ("C", 1.0, "")]),
    ("VAL", [("N", 1.0, ""), ("C", 1.0, ""), ("O", 1.0, "")]),          # no CA
    ("MSE", [("N", 1.0, ""), ("CA", 1.0, ""), ("C", 1.0, ""), ("SE", 1.0, "")]),
    ("UNK", [("N", 1.0, ""), ("CA", 1.0, ""), ("C", 1.0, "")]),
    ("PRO", [("N", 1.0, ""), ("CA", 1.0, ""), ("C", 1.0, ""), ("O", 1.0, "")]),
]


def _build_atomgroup():
    names, resnames, resnums, occs, altlocs = [], [], [], [], []
    for i, (resname, atoms) in enumerate(RESIDUES):
        for name, occ, altloc in atoms:
            names.append(name)
            resnames.append(resname)
            resnums.append(i + 1)
            occs.append(occ)
            altlocs.append(altloc)
    n = len(names)
    ag = pr.AtomGroup("flank_index_test")
    ag.setCoords((np.random.default_rng(0).random((n, 3)) * 30).astype(np.float32))
    ag.setNames(np.array(names))
    ag.setResnames(np.array(resnames))
    ag.setResnums(np.array(resnums))
    ag.setChids(np.array(["A"] * n))
    ag.setOccupancies(np.array(occs, dtype=float))
    ag.setAltlocs(np.array(altlocs))
    return ag


class FlankLookupIndexTests(unittest.TestCase):
    def setUp(self):
        self.ag = _build_atomgroup()
        self.index = build_flank_lookup_index(self.ag)

    def test_indexed_and_selection_paths_agree_everywhere(self):
        # Range runs past both ends: the flank walk asks for vdm_resind +/- f,
        # which goes negative at the start of an environment atomgroup. A
        # positional index would wrap -1 to the last residue and return a real
        # neighbour there.
        for resindex in range(-3, len(RESIDUES) + 3):
            with self.subTest(resindex=resindex):
                slow_aa, slow_ca = get_AA_and_CA_coords(self.ag, resindex)
                fast_aa, fast_ca = get_AA_and_CA_coords(
                    self.ag, resindex, flank_index=self.index)
                self.assertEqual(slow_aa, fast_aa)
                self.assertEqual(slow_ca.dtype, np.float32)
                self.assertEqual(fast_ca.dtype, np.float32)
                if np.isnan(slow_ca).all():
                    self.assertTrue(np.isnan(fast_ca).all())
                else:
                    np.testing.assert_array_equal(slow_ca, fast_ca)

    def test_only_unambiguous_residues_are_indexed(self):
        # ALA, GLY, MSE, PRO -- not the ligand, water, altloc pairs, CA-less
        # residue or UNK, all of which must fall through.
        self.assertEqual(sorted(self.index), [0, 1, 7, 9])

    def test_out_of_range_resindices_fall_through(self):
        for resindex in (-1, -2, len(RESIDUES), len(RESIDUES) + 1):
            aa, ca = get_AA_and_CA_coords(self.ag, resindex, flank_index=self.index)
            self.assertEqual(aa, FLANK_MISSING)
            self.assertTrue(np.isnan(ca).all())

    def test_returned_coords_do_not_alias_the_index(self):
        # The selection path hands back a fresh array; a caller that writes to
        # one must not corrupt every later lookup in the environment.
        _, first = get_AA_and_CA_coords(self.ag, 0, flank_index=self.index)
        first[:] = 999.0
        _, second = get_AA_and_CA_coords(self.ag, 0, flank_index=self.index)
        self.assertFalse(np.allclose(second, 999.0))

    def test_unindexable_atomgroup_falls_through(self):
        # An empty AtomGroup has no arrays to read; the caller must still work.
        empty = pr.AtomGroup("empty")
        self.assertIsNone(build_flank_lookup_index(empty))


if __name__ == "__main__":
    unittest.main()
