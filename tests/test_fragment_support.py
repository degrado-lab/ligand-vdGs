"""The ratified support definition (DR-5), tested on its two falsifiers.

    support(f) = |{ biounit b : some ligand instance in b has a match of f with
                    every atom of m observed }|

Four units have been in circulation for this quantity and three of them are wrong.
Two implementations that count a wrong unit still pass the obvious test, so the
inputs here are chosen to separate them:

  * an NCS-heavy biounit -- a fragment in 1000 copies must score 1, not 1000.
    An implementation that counts ligand copies or CG sites fails only here.
  * one biounit holding two DIFFERENT CCD codes that share a fragment (ATP and
    ADP sharing adenine) -- it must score 1, not 2. An implementation that SUMS
    over codes instead of unioning passes the NCS test and fails only here, which
    is why the sum is computed alongside and asserted to differ.
  * an observed-atom mask that omits one atom of the fragment -- the fragment must
    vanish, because a match reaching an unobserved (phantom) atom can never be
    mined. 41.4% of ligands have at least one unobserved heavy atom.

Synthetic rosters throughout: the join is pure bookkeeping and must be testable
without the CCD store or the parent database.
"""
import unittest
from collections import Counter

from ligand_vdgs.generate_vdgs.fragment_database_ligs import support_from_roster

ADENINE = '[n;D2][c;!H0][n;D2][c;H0][c;H0]'
RIBOSE = '[C;!H0][O;r5;D2][C;!H0][C;!H0][O;!R;D1]'

# ATP and ADP observe different atom sets; both carry the adenine fragment.
ATP_NAMES = ('C2', 'C4', 'C5', 'N1', 'N3', 'N7', 'N9', "C1'", "O4'", "C2'", "O2'")
ADP_NAMES = ('C2', 'C4', 'C5', 'N1', 'N3', 'N7', 'N9')

KEYS_BY_RESNAME = {
    'ATP': {ADENINE: [('C2', 'C4', 'N1', 'N3', 'N7')],
            RIBOSE: [("C1'", "C2'", "O2'", "O4'", 'C4')]},
    'ADP': {ADENINE: [('C2', 'C4', 'N1', 'N3', 'N7')]},
}


def _roster(types, biounits):
    return {'types': types, 'biounits': biounits, 'db_identity': {'sha256': 'test'}}


class SupportUnitTests(unittest.TestCase):

    def test_ncs_copies_of_one_ligand_count_once(self):
        """1000 copies in one biounit -> support 1."""
        types = [('ATP', ATP_NAMES)]
        # The same type id repeated, which is what an NCS-heavy structure produces
        # before dedup. Passing duplicates in is the point: it tests np.unique
        # rather than the roster builder's own set().
        roster = _roster(types, {'1abc': [0] * 1000})
        support, _ = support_from_roster(roster, KEYS_BY_RESNAME)
        self.assertEqual(support[ADENINE], 1)
        self.assertEqual(support[RIBOSE], 1)

    def test_two_ccd_codes_in_one_biounit_union_rather_than_sum(self):
        """ATP and ADP together contribute 1, not 2, to a shared fragment."""
        types = [('ATP', ATP_NAMES), ('ADP', ADP_NAMES)]
        roster = _roster(types, {'1abc': [0, 1]})
        support, _ = support_from_roster(roster, KEYS_BY_RESNAME)
        self.assertEqual(support[ADENINE], 1)

        # Vacuity clause: a summing implementation must give a DIFFERENT answer on
        # this input, or the test would pass against the bug it exists to catch.
        summed = Counter()
        for type_ids in roster['biounits'].values():
            for type_idx in type_ids:
                resname, names = types[type_idx]
                for key, name_sets in KEYS_BY_RESNAME[resname].items():
                    if any(set(names).issuperset(n) for n in name_sets):
                        summed[key] += 1
        self.assertEqual(summed[ADENINE], 2)
        self.assertNotEqual(summed[ADENINE], support[ADENINE])

    def test_distinct_biounits_do_add_up(self):
        """The union is per biounit, not global -- two biounits give 2."""
        types = [('ATP', ATP_NAMES), ('ADP', ADP_NAMES)]
        roster = _roster(types, {'1abc': [0], '2xyz': [1]})
        support, _ = support_from_roster(roster, KEYS_BY_RESNAME)
        self.assertEqual(support[ADENINE], 2)

    def test_unobserved_atom_removes_the_fragment(self):
        """A match reaching an atom the structure never observed cannot be mined."""
        partial = tuple(n for n in ATP_NAMES if n != 'N3')
        roster = _roster([('ATP', partial)], {'1abc': [0]})
        support, _ = support_from_roster(roster, KEYS_BY_RESNAME)
        self.assertEqual(support.get(ADENINE, 0), 0)
        # Discriminating: the OTHER fragment of the same ligand, whose atoms are
        # all still observed, must be untouched. Without this the test would also
        # pass for an implementation that drops the ligand entirely.
        self.assertEqual(support[RIBOSE], 1)

        restored = _roster([('ATP', ATP_NAMES)], {'1abc': [0]})
        self.assertEqual(support_from_roster(restored, KEYS_BY_RESNAME)[0][ADENINE], 1)

    def test_a_fragment_can_be_reached_through_either_of_two_atom_sets(self):
        """Two matches of one key in one ligand: either being observed suffices."""
        keys = {'LIG': {ADENINE: [('A1', 'A2', 'A3', 'A4', 'A5'),
                                  ('B1', 'B2', 'B3', 'B4', 'B5')]}}
        only_second = _roster([('LIG', ('B1', 'B2', 'B3', 'B4', 'B5'))], {'1abc': [0]})
        self.assertEqual(support_from_roster(only_second, keys)[0][ADENINE], 1)
        neither = _roster([('LIG', ('A1', 'A2', 'A3', 'B4', 'B5'))], {'1abc': [0]})
        self.assertEqual(support_from_roster(neither, keys)[0].get(ADENINE, 0), 0)

    def test_biounit_with_no_recorded_type_contributes_nothing(self):
        roster = _roster([('ATP', ATP_NAMES)], {'1abc': [], '2xyz': [0]})
        support, stats = support_from_roster(roster, KEYS_BY_RESNAME)
        self.assertEqual(support[ADENINE], 1)
        self.assertEqual(stats['num_biounits'], 2)


if __name__ == '__main__':
    unittest.main()
