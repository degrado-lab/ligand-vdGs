"""record_frag's alias cache must be a pure speedup.

The dedup scan is O(bucket^2) and, before the cache, re-ran in full for every
ligand carrying a degenerate key -- such a key never becomes a bucket key of its
own, so the fast `smiles in frag_dict[elements]` path never catches it. The risk
of caching is that the resolved key drifts, so the discriminating input here is a
degenerate key seen repeatedly, interleaved with its canonical partner and with a
repeated ligand name: the resulting dict must be byte-identical to the uncached
one, and the degenerate key must still not appear as a key.
"""
import unittest

from rdkit import Chem

from ligand_vdgs.generate_vdgs.fragment_database_ligs import (fragment_elements,
                                                              record_frag)
from ligand_vdgs.functions.utils import (fragment_key_query_mol,
                                         fragment_query_mols_equivalent)

KEY_A = '[C;!R][C;!R](=[O;!R])[N;!R]'
KEY_B = '[N;!R][C;!R](=[O;!R])[C;!R]'   # same fragment, written backwards


def _submol(key):
    mol = Chem.MolFromSmiles(key, sanitize=False)
    return mol if mol is not None else Chem.MolFromSmarts(key)


class AliasCacheTests(unittest.TestCase):
    def setUp(self):
        # Guard the premise: if these two stopped comparing equivalent the test
        # below would pass vacuously, never entering the alias branch at all.
        self.assertTrue(fragment_query_mols_equivalent(
            fragment_key_query_mol(KEY_A), fragment_key_query_mol(KEY_B)))
        self.seq = [(KEY_A, 'LIG1'), (KEY_B, 'LIG2'), (KEY_B, 'LIG3'),
                    (KEY_B, 'LIG2'), (KEY_A, 'LIG4'), (KEY_B, 'LIG5')]

    def _run(self, use_cache):
        frag_dict, smarts_cache, errors, alias = {}, {}, {}, {}
        for key, lig in self.seq:
            # record_frag takes the element string, not the Mol: the worker that
            # produces it in the parallel path has already discarded the Mol.
            record_frag(fragment_elements(_submol(key)), frag_dict, key, lig,
                        smarts_cache, errors, alias if use_cache else {})
        return frag_dict

    def test_cache_does_not_change_the_dict(self):
        self.assertEqual(self._run(True), self._run(False))

    def test_degenerate_key_stays_folded_into_the_first_one(self):
        frag_dict = self._run(True)
        keys = {k for sub in frag_dict.values() for k in sub}
        self.assertEqual(keys, {KEY_A})
        ligs = next(iter(frag_dict.values()))[KEY_A]
        # Every ligand once, in first-seen order; the repeat of LIG2 is not
        # appended twice.
        self.assertEqual(ligs, ['LIG1', 'LIG2', 'LIG3', 'LIG4', 'LIG5'])


if __name__ == '__main__':
    unittest.main()
