"""`Frags.KEY_SCHEMA` must be bumped whenever the key vocabulary changes.

The failure this closes: `resources/database_frags_dict.pkl` and
`resources/frag_cost_estimate.tsv` were both written in the pre-annotation
vocabulary (`[O;!R]=[C;!R]([O;!R])...` where the code now emits
`[O;!R;D1]=[C;!R;H0](...)`). `utils.fragment_keys_equivalent(old, new)` is False,
so every lookup missed -- and nothing caught it, because the two artifacts were
sha-consistent with EACH OTHER and no guard compared either to the code.

The version string alone does not close it: someone changes a primitive and
forgets to bump. So the golden keys below pin the vocabulary itself. Changing
`ring_query_for_atom`, the `D<n>` rule or the `H0`/`!H0` rule changes a golden
key and fails here, and the fix is to update both the key AND `KEY_SCHEMA` --
at which point every reader of an old artifact refuses it.

Golden set chosen for coverage of the vocabulary, not for being easy: aromatic
(bare, no ring primitive), saturated 5-ring and 6-ring (`r5`/`r6`, which must stay
distinct -- THF and THP are different fragments while pyridine and pyrrole pool),
acyclic (`!R`), a terminal vs. bridging heteroatom pair (the `D<n>` distinction
that restricts the automorphism group), and a carbon `H0`/`!H0` pair.
"""
import unittest

from rdkit import Chem

from ligand_vdgs.functions import Frags

# (name, SMILES, atom-count window) -> expected keys.
GOLDEN = {
    # r5 vs r6 on the same heteroatom: THF and THP must not pool.
    'THF (whole 5-ring)': (
        'C1CCOC1', 5, 5,
        {'[C;r5;!H0]1[C;r5;!H0][C;r5;!H0][O;r5;D2][C;r5;!H0]1'}),
    'pyrrolidine (whole 5-ring)': (
        'C1CCNC1', 5, 5,
        {'[C;r5;!H0]1[C;r5;!H0][C;r5;!H0][N;r5;D2][C;r5;!H0]1'}),
    # Aromatic atoms stay bare; a 6-ring is only ever seen as its arcs at size 5,
    # and the all-carbon arc is dropped.
    'pyridine (arcs of a 6-ring)': (
        'c1ccncc1', 5, 5,
        {'[c;!H0][c;!H0][c;!H0][c;!H0][n;D2]',
         '[c;!H0][c;!H0][c;!H0][n;D2][c;!H0]',
         '[c;!H0][c;!H0][n;D2][c;!H0][c;!H0]'}),
    # Terminal (D1) vs bridging (D2) oxygen in one fragment: the distinction that
    # keeps a bridging O from mapping onto a terminal one.
    'methyl acetate': (
        'COC(C)=O', 5, 5,
        {'[C;!R;!H0][O;!R;D2][C;!R;H0]([C;!R;!H0])=[O;!R;D1]'}),
}


def _keys(smiles, lo, hi):
    mol = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(mol)
    return set(Frags.enumerate_induced_fragments(mol, lo, hi))


class KeySchemaTests(unittest.TestCase):

    def test_key_schema_is_a_nonempty_version_string(self):
        self.assertIsInstance(Frags.KEY_SCHEMA, str)
        self.assertTrue(Frags.KEY_SCHEMA)

    def test_golden_keys_pin_the_vocabulary(self):
        for name, (smiles, lo, hi, expected) in GOLDEN.items():
            with self.subTest(fragment=name, key_schema=Frags.KEY_SCHEMA):
                self.assertEqual(
                    _keys(smiles, lo, hi), expected,
                    f'{name} no longer keys as it did under KEY_SCHEMA '
                    f'{Frags.KEY_SCHEMA!r}. If the change is intended, update the '
                    f'golden key AND bump Frags.KEY_SCHEMA, so every artifact '
                    f'written under the old vocabulary is refused.')

    def test_ring_sizes_are_not_interchangeable(self):
        """Vacuity clause: the golden set must actually discriminate r5 from r6.

        If `ring_query_for_atom` collapsed to a bare `R`, THF and THP would key
        alike and the golden set above would still pass on its own entries.
        """
        thf = _keys('C1CCOC1', 5, 5)
        thp_ring = _keys('C1CCOCC1', 6, 6)
        self.assertTrue(thf)
        self.assertTrue(thp_ring)
        self.assertEqual(thf & thp_ring, set())
        self.assertTrue(any('r5' in k for k in thf))
        self.assertTrue(any('r6' in k for k in thp_ring))


if __name__ == '__main__':
    unittest.main()
