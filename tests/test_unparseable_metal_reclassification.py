"""The undesired-element filter must parse atoms, never match strings.

REMOVED BEHAVIOUR, recorded here because the removal is deliberate and is FORCED BY
THE SUPPORT UNIT, not a convenience. Support is counted in distinct parent biounits
(DR-5), so a ligand that exists only as a SMILES string has no biounit to count and
can never clear any threshold; and an instance whose resname has no CCD template
falls back to OpenBabel perception, which the roster excludes, so it contributes no
support either. Neither can reach the vocabulary, so neither needs a SMILES path.
The enumerator therefore no longer reads a CCD SMILES file at all -- ligand chemistry
comes from the CCD template store, which is also what the miner and the instance
side use. The whole
`'|'`-dative-bond reclassification (`undesired_elements_in_dative_smiles`) existed
only because the CCD writes dative bonds as `|`, which RDKit rejects, so metal
complexes died at `MolFromSmiles` before `undesired_elements` ever saw them. From
a template there is no SMILES to fail to parse: the metal is an atom in the graph
and `undesired_elements` sees it directly. The end-to-end half of this file went
with it; the classification it protected is now this file's subject.

What remains is the part that was never about SMILES syntax: `undesired_elements`
decides on PARSED atom symbols, so a raw-string matcher's false positives must not
reappear. `In1cccc1` is iodine bonded to an aromatic N, not indium; `CSnc1ccccc1`
is C-S-n, not tin; `[Ala]` is not aluminium.
"""
import unittest

from rdkit import Chem

from ligand_vdgs.generate_vdgs.fragment_database_ligs import (UNDESIRED_ELEMENTS,
                                                              undesired_elements)


def _has_undesired(smiles):
    """Parse unsanitized, as the enumerator's element filter does.

    Unsanitized on purpose: a metal complex that fails RDKit's valence check must
    be classified as non-druglike, not as a parse failure.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    return undesired_elements(mol)


class UndesiredElementTests(unittest.TestCase):
    """The element check must run on atoms, not on characters."""

    def test_unbracketed_metal_spellings_do_not_match(self):
        # The exact false positives the raw-string approach was abandoned for.
        # `CCa1ccccc1` does not parse at all, so it can only be asserted not to
        # read as a metal; the other two must parse and read False, which is what
        # keeps this from passing merely because nothing parsed.
        for smi in ('In1cccc1', 'CSnc1ccccc1'):
            self.assertIs(_has_undesired(smi), False, smi)
        self.assertIsNot(_has_undesired('CCa1ccccc1'), True)

    def test_bracket_syntax_is_not_mistaken_for_an_element(self):
        for smi in ('[13C]', '[2H+]', '[C@@H](N)C', '[nH]1cccc1', '[N+](=O)[O-]'):
            self.assertIs(_has_undesired(smi), False, smi)

    def test_bracketed_undesired_elements_match(self):
        # Isotope/charge decoration and the lowercase aromatic metalloids.
        for smi in ('[Fe]', '[Sn+2]', '[Y]', '[K+]', '[Yb+3]',
                    '[se]1cccc1', '[te]1cccc1', '[as]1ccccc1'):
            self.assertIs(_has_undesired(smi), True, smi)

    def test_a_real_metal_complex_is_caught_even_unsanitized(self):
        """Zinc acetate fails a valence check; it must still read as non-druglike."""
        self.assertIs(_has_undesired('CC(=O)O[Zn]OC(C)=O'), True)

    def test_druglike_molecules_pass(self):
        for smi in ('CC(=O)Nc1ccc(O)cc1', 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C'):
            self.assertIs(_has_undesired(smi), False, smi)

    def test_boron_is_not_undesired(self):
        # Boronic-acid warheads are deliberately mined; guards the element table.
        self.assertNotIn('B', UNDESIRED_ELEMENTS)
        self.assertIs(_has_undesired('OB(O)c1ccccc1'), False)

    def test_the_filter_reaches_templated_metal_ligands(self):
        """The path that replaced the SMILES one: a template carrying a metal.

        Skipped rather than failed when the CCD store is unavailable, but when it
        is available this is the real check -- an Fe-S cluster must be excluded by
        the element filter alone, with no SMILES involved.
        """
        from ligand_vdgs.functions import ligand_perception
        try:
            ligand_perception.require_template_store()
        except Exception as exc:                      # pragma: no cover
            self.skipTest(f'CCD template store unavailable: {exc}')
        perceived = ligand_perception.perceive_ligand_graph('FES')
        self.assertIsNotNone(perceived.mol)
        self.assertTrue(undesired_elements(perceived.mol))
        # Discriminating: an ordinary templated ligand must pass, or this would
        # also hold for a filter that rejected every template.
        atp = ligand_perception.perceive_ligand_graph('ATP')
        self.assertIsNotNone(atp.mol)
        self.assertFalse(undesired_elements(atp.mol))


if __name__ == '__main__':
    unittest.main()
