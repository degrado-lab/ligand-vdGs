"""Ring-context annotation of fragment keys.

The discriminating cases first: fragments whose atom sequence is identical but
whose ring context is not (piperidine vs an acyclic amine; THF vs THP; ribose vs
glucose), aromatic rings that must *not* split by size, bridged systems where
RDKit and OpenBabel must agree on smallest ring or a query silently never finds
its library, and the two places the annotation is easy to break silently -- the
submol self-match (the ring is already cut there) and charge normalization (the
charge is no longer at the end of the bracket).
"""
import unittest

from rdkit import Chem, RDLogger

from ligand_vdgs.functions.Frags import get_fragments, ring_query_for_atom
from ligand_vdgs.functions.utils import (identify_mol_automorphisms,
                                         mol_from_fragment,
                                         fragment_keys_equivalent)
from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
    charge_normalized_fragment)

RDLogger.DisableLog('rdApp.*')


def frag_keys(smiles):
    return set(get_fragments(2, Chem.MolFromSmiles(smiles), 4, 5))


class RingQueryTests(unittest.TestCase):
    def test_fused_aromatic_is_left_bare(self):
        # Indole's fusion carbons are in both a 5- and a 6-ring. Aromatic atoms
        # get no annotation at all, so nothing here can shatter by ring size.
        mol = Chem.MolFromSmiles('c1ccc2[nH]ccc2c1')
        self.assertEqual({ring_query_for_atom(a) for a in mol.GetAtoms()}, {None})

    def test_saturated_ring_atoms_carry_smallest_ring_size(self):
        # Fused aliphatic: the shared bond's atoms are in both rings; smallest wins.
        mol = Chem.MolFromSmiles('C1CC2CC2C1')   # bicyclo[3.1.0]hexane
        self.assertEqual(sorted(ring_query_for_atom(a) for a in mol.GetAtoms()),
                         ['r3', 'r3', 'r3', 'r5', 'r5', 'r5'])

    def test_aromatic_atoms_are_left_bare(self):
        # Lower case already implies a ring; annotating adds a constant.
        mol = Chem.MolFromSmiles('c1ccncc1')
        self.assertEqual({ring_query_for_atom(a) for a in mol.GetAtoms()}, {None})

    def test_macrocycle_gets_its_exact_ring_size_never_plain_R(self):
        # SMARTS `R` matches any ring atom, so a plain-R macrocycle key would mine
        # every THF and sugar in the mirror (estimated 27,399 structures for 42
        # ligands' worth of key). `r12` mines macrocycles only.
        mol = Chem.MolFromSmiles('C1COCCOCCOCCO1')
        self.assertEqual({ring_query_for_atom(a) for a in mol.GetAtoms()}, {'r12'})

    def test_ring_keys_are_disjoint_as_queries_in_both_toolkits(self):
        # The property that makes fragment libraries partition rather than nest.
        from openbabel import openbabel as ob
        conv = ob.OBConversion(); conv.SetInFormat('smi')
        # Every crown key, not an arbitrary one: only some of them catch THF
        # under a plain-`R` regression, so `next(iter(...))` passed at coin-flip
        # rate depending on set iteration order.
        pairs = [(k, 'C1COCCOCCOCCO1', 'C1CCOC1') for k in frag_keys('C1COCCOCCOCCO1')]
        pairs += [(k, 'C1CCOC1', 'C1COCCOCCOCCO1')
                  for k in frag_keys('C1CCOC1') if k.count('[') == 5]
        self.assertTrue(pairs)
        for key, should_match, should_not in pairs:
            q = mol_from_fragment(key)
            sp = ob.OBSmartsPattern(); self.assertTrue(sp.Init(key), key)
            for smi, expect in [(should_match, True), (should_not, False)]:
                self.assertEqual(Chem.MolFromSmiles(smi).HasSubstructMatch(q), expect, (key, smi))
                m = ob.OBMol(); conv.ReadString(m, smi)
                self.assertEqual(bool(sp.Match(m)), expect, (key, smi))

    def test_acyclic_atoms_are_marked_not_in_ring(self):
        mol = Chem.MolFromSmiles('CCOCC')
        self.assertEqual({ring_query_for_atom(a) for a in mol.GetAtoms()}, {'!R'})


class FragmentKeyTests(unittest.TestCase):
    def test_ring_and_chain_of_the_same_atoms_get_different_keys(self):
        # The case RMSD cannot separate: 100% of chain-context CCNC vdGs have a
        # ring-context one inside the hit threshold.
        self.assertEqual(frag_keys('C1CCNCC1') & frag_keys('CCCCNC'), set())

    def test_five_and_six_membered_saturated_rings_get_different_keys(self):
        # Furanose vs pyranose is the largest protein-side divergence in the
        # library (26x null); THF/THP and pyrrolidine/piperidine ride along.
        self.assertEqual(frag_keys('C1CCOC1') & frag_keys('C1CCOCC1'), set())
        self.assertEqual(frag_keys('C1CCNC1') & frag_keys('C1CCNCC1'), set())
        ribose, glucose = frag_keys('OCC1OC(O)C(O)C1O'), frag_keys('OCC1OC(O)C(O)C(O)C1O')
        self.assertTrue(ribose and glucose)
        self.assertEqual(ribose & glucose, set())

    def test_five_and_six_membered_aromatics_share_keys(self):
        # Deliberate: pyridine and pyrrole pool. They stay separate pose clusters
        # inside the shared bucket (~0.63 A apart against a 0.5 A cutoff), so the
        # geometry survives and only the ring-size label is given up -- which is
        # the trade that keeps rare fragments from shattering.
        self.assertTrue(frag_keys('c1ccncc1') & frag_keys('c1cc[nH]c1'))

    def test_ring_derived_fragment_survives_the_submol_self_match(self):
        # PathToSubmol cuts the ring, so an annotated key matched back against the
        # submol finds nothing. If the two roles are ever merged again, this
        # returns empty rather than raising. Saturated ring, since aromatics are
        # deliberately left bare.
        keys = frag_keys('C1CCNCC1')
        self.assertTrue(keys)
        self.assertTrue(all(';r6]' in k for k in keys), keys)

    def test_rdkit_and_openbabel_agree_on_smallest_ring_in_bridged_systems(self):
        # Query keys come from RDKit ring perception; mining matches with
        # OpenBabel [r<n>]. A disagreement is a query that silently never finds
        # its library. SSSR is where toolkits are most likely to differ.
        from openbabel import openbabel as ob
        conv = ob.OBConversion(); conv.SetInFormat('smi')
        pats = {}
        for n in range(3, 9):
            pats[n] = ob.OBSmartsPattern(); pats[n].Init(f'[r{n}]')
        for name, smi in [('norbornane', 'C1CC2CCC1C2'), ('quinuclidine', 'C1CN2CCC1CC2'),
                          ('adamantane', 'C1C2CC3CC1CC(C2)C3'), ('tropane', 'CN1C2CCC1CC(O)C2'),
                          ('cubane', 'C12C3C4C1C5C2C3C45'), ('camphor', 'CC1(C)C2CCC1(C)C(=O)C2'),
                          ('morphinan', 'C1CCC2C3CCCCC3CC2C1'), ('spiro[4.5]decane', 'C1CCC2(C1)CCCCC2')]:
            mol = Chem.MolFromSmiles(smi)
            rd = [next((n for n in range(3, 9) if a.IsInRingSize(n)), 0) if a.IsInRing() else 0
                  for a in mol.GetAtoms()]
            om = ob.OBMol(); conv.ReadString(om, smi)
            obsz = [0] * om.NumAtoms()
            for n in range(8, 2, -1):
                if pats[n].Match(om):
                    for (i,) in pats[n].GetUMapList():
                        obsz[i - 1] = n
            self.assertEqual(rd, obsz, name)

    def test_partially_cyclic_fragment_annotates_per_atom(self):
        # Benzoate: ring carbon plus an exocyclic carboxylate.
        keys = frag_keys('c1ccccc1C(=O)[O-]')
        self.assertIn('c[C;!R](=[O;!R])[O-;!R]', keys)

    def test_keys_parse_as_smarts_in_both_toolkits(self):
        from openbabel import openbabel as ob
        for key in sorted(frag_keys('C1CCNCC1') | frag_keys('c1ccccc1C(=O)[O-]')):
            self.assertIsNotNone(mol_from_fragment(key), key)
            pattern = ob.OBSmartsPattern()
            self.assertTrue(pattern.Init(key), key)

    def test_automorphisms_still_found_on_annotated_keys(self):
        # The carboxylate oxygens stay interchangeable through the annotation.
        perms = identify_mol_automorphisms(
            mol_from_fragment('c[C;!R](=[O;!R])[O-;!R]'))
        self.assertEqual(len(perms), 2)


class KeyEquivalenceTests(unittest.TestCase):
    """Degenerate spellings must merge; different ring context must not.

    Regression: dedup compared a key against the *submol* it was cut from,
    which has no RingInfo, so every annotated key raised and the duplicate was
    recorded as a distinct fragment -- 190 of 1029 selected fragments.
    """

    def test_degenerate_spellings_of_one_fragment_are_equivalent(self):
        for a, b in [('[C;!R][C;!R][O;!R][C;!R]', '[C;!R][O;!R][C;!R][C;!R]'),
                     ('[O;!R][C;r5]([O;r5])[O;r5]', '[O;r5][C;r5]([O;r5])[O;!R]'),
                     ('cc([C;!R])o', '[C;!R]c(c)o')]:
            self.assertTrue(fragment_keys_equivalent(a, b), (a, b))

    def test_equivalence_is_reflexive_for_ring_keys(self):
        # Direct matching of two annotated keys is not: `r<n>` is unsatisfiable
        # on the acyclic graph a key parses to.
        for key in ['[C;r6][C;r6][O;r6][C;r6][C;r6]', '[C;r5][C;r5]([C;r5])[O;!R]']:
            self.assertTrue(fragment_keys_equivalent(key, key), key)

    def test_different_ring_context_is_not_equivalent(self):
        for a, b in [('[C;r5][C;r5][O;r5][C;r5][C;r5]', '[C;r6][C;r6][O;r6][C;r6][C;r6]'),
                     ('[C;!R][C;!R][O;!R][C;!R][C;!R]', '[C;r5][C;r5][O;r5][C;r5][C;r5]'),
                     ('[C;r5][C;r5]([C;r5])[O;!R]', '[C;r5][C;r5]([C;!R])[O;r5]'),
                     ('[C;!R][C;!R](=[O;!R])[O;!R]', '[C;!R][C;!R](=[O;!R])[O-;!R]')]:
            self.assertFalse(fragment_keys_equivalent(a, b), (a, b))


class ChargeNormalizationTests(unittest.TestCase):
    def test_charge_is_stripped_from_an_annotated_bracket(self):
        # The charge is mid-bracket now; a $-anchored pattern silently misses it,
        # which stops every protonation variant from collapsing.
        self.assertEqual(
            charge_normalized_fragment('c[C;!R](=[O;!R])[O-;!R]'),
            'c[C;!R](=[O;!R])[O;!R]')

    def test_annotation_is_preserved_while_charge_is_dropped(self):
        self.assertEqual(charge_normalized_fragment('[O-;r5][P;!R]'), '[O;r5][P;!R]')

    def test_uncharged_annotated_key_is_unchanged(self):
        self.assertIsNone(charge_normalized_fragment('[C;r6][N;r6]'))

    def test_unannotated_keys_still_normalize(self):
        self.assertEqual(charge_normalized_fragment('CC(=O)[O-]'), 'CC(=O)O')


if __name__ == '__main__':
    unittest.main()
