"""`collapse_keys` must fold degenerate spellings onto one representative.

Two spellings of one fragment still occur with canonical SMILES: the annotation
(`r<n>`, `D<n>`, `H0`) is overlaid on the written string afterwards and is
invisible to RDKit's canonical ranking, so when a graph automorphism is broken
only by the annotation, either placement can be written. Measured on a
3,000-structure roster: 655 of 24,944 raw keys were degenerate.

The collapse must happen at KEY level and before the roster join. Collapsing after
it would leave two biounit counts to merge and neither merge is right -- summing
double-counts a biounit carrying both spellings, taking the max undercounts one
carrying only one of each.

Discriminating inputs, in order:
  * a reversal pair, which is the common degeneracy;
  * a key in the same element bucket that is a genuinely DIFFERENT fragment and
    must stay its own representative (a test with only the pair would pass for an
    implementation that collapses everything);
  * the `_key_token_signature` pre-grouping, which claims to be a necessary
    condition for equivalence -- asserted by running the collapse with the
    grouping disabled and requiring the identical alias map, plus a clause that
    fails if nothing collapsed at all.
"""
import unittest

from ligand_vdgs.generate_vdgs import fragment_database_ligs as F
from ligand_vdgs.functions.utils import (fragment_key_query_mol,
                                         fragment_query_mols_equivalent)

KEY_A = '[C;!R][C;!R](=[O;!R])[N;!R]'
KEY_B = '[N;!R][C;!R](=[O;!R])[C;!R]'   # same fragment, written backwards
KEY_C = '[C;!R][C;!R](=[O;!R])[O;!R]'   # a different fragment, different elements
KEY_D = '[C;!R][N;!R]([C;!R])[C;!R]'    # CCCN, same elements as A/B, different graph

ELEMENTS = {KEY_A: 'CCNO', KEY_B: 'CCNO', KEY_C: 'CCOO', KEY_D: 'CCCN'}


class CollapseKeyTests(unittest.TestCase):

    def setUp(self):
        # Guard the premise: if these stopped comparing equivalent, every test
        # below would pass vacuously without ever entering the collapse branch.
        self.assertTrue(fragment_query_mols_equivalent(
            fragment_key_query_mol(KEY_A), fragment_key_query_mol(KEY_B)))
        self.assertFalse(fragment_query_mols_equivalent(
            fragment_key_query_mol(KEY_A), fragment_key_query_mol(KEY_D)))

    def test_reversal_pair_collapses_onto_one_representative(self):
        alias, n_collapsed = F.collapse_keys(ELEMENTS, {}, log=lambda *a: None)
        self.assertEqual(n_collapsed, 1)
        self.assertEqual(alias[KEY_A], alias[KEY_B])
        # Deterministic representative: sorted iteration, so the same input always
        # yields the same directory name downstream.
        self.assertEqual(alias[KEY_A], min(KEY_A, KEY_B))

    def test_a_different_fragment_stays_its_own_representative(self):
        alias, _ = F.collapse_keys(ELEMENTS, {}, log=lambda *a: None)
        self.assertEqual(alias[KEY_C], KEY_C)
        self.assertEqual(alias[KEY_D], KEY_D)
        self.assertEqual(len(set(alias.values())), 3)

    def test_token_grouping_is_a_necessary_condition(self):
        """Grouping must not separate a pair the full check would have merged."""
        fast, n_fast = F.collapse_keys(ELEMENTS, {}, log=lambda *a: None)
        original = F._key_token_signature
        try:
            F._key_token_signature = lambda key: 0     # group by element string only
            slow, n_slow = F.collapse_keys(ELEMENTS, {}, log=lambda *a: None)
        finally:
            F._key_token_signature = original
        self.assertEqual(fast, slow)
        # Vacuity clause: comparing two alias maps in which nothing collapsed
        # would prove nothing.
        self.assertGreater(n_slow, 0)
        self.assertEqual(n_fast, n_slow)

    def test_signature_ignores_bond_and_ring_symbols_but_not_atoms(self):
        """Two spellings of one fragment must land in the same group."""
        self.assertEqual(F._key_token_signature(KEY_A), F._key_token_signature(KEY_B))
        self.assertNotEqual(F._key_token_signature(KEY_A),
                            F._key_token_signature(KEY_C))
        # A cyclic and an acyclic spelling of the same atoms share a signature --
        # they must, or a ring-closing key could never be compared to its
        # path-written twin.
        cyclic = '[C;r5;!H0]1[C;r5;!H0][C;r5;!H0][N;r5;D2][C;r5;!H0]1'
        acyclic = '[C;r5;!H0][C;r5;!H0][C;r5;!H0][N;r5;D2][C;r5;!H0]'
        self.assertEqual(F._key_token_signature(cyclic),
                         F._key_token_signature(acyclic))

    def test_an_unparseable_key_is_its_own_representative(self):
        elements = dict(ELEMENTS)
        elements['[C;!R][C;!R'] = 'CCNO'      # unbalanced bracket
        alias, _ = F.collapse_keys(elements, {}, log=lambda *a: None)
        self.assertEqual(alias['[C;!R][C;!R'], '[C;!R][C;!R')


if __name__ == '__main__':
    unittest.main()
