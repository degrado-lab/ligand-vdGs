import os
import unittest

from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (is_solvent_artifact,
    select_fragments)
from ligand_vdgs.functions.utils import mol_from_fragment


def _is_artifact(smarts):
    return is_solvent_artifact(mol_from_fragment(smarts))


class SolventArtifactFilterTests(unittest.TestCase):
    def test_halogen_oxyanions_are_rejected(self):
        for smarts in ("O=[Cl](=O)(=O)[O-]", "O=[Cl](=O)[O-]",
                       "O=I(=O)(=O)[O-]", "O=[Br](=O)O"):
            self.assertTrue(_is_artifact(smarts), smarts)

        # Detection is structural, so an alternate drawing is caught too.
        self.assertTrue(_is_artifact("[O-][Cl]([O-])([O-])[O-]"))

    def test_ligand_chemistry_is_kept(self):
        # Phosphates and sulfates are real binding moieties despite also being
        # common buffer components.
        for smarts in ("O=P(O)(O)O", "O=S(=O)(O)O", "[O-]S([O-])([O-])O",
                       "NS(=O)(=O)O", "CS(=O)(=O)O"):
            self.assertFalse(_is_artifact(smarts), smarts)

        # A halogen with fewer than two oxygen neighbors is ordinary chemistry.
        for smarts in ("CC(S)[Cl]", "cOP(=O)[O-]", "SC([Cl])([Cl])[Br]"):
            self.assertFalse(_is_artifact(smarts), smarts)

        # A carbon-bound hypervalent halogen oxide is real organohalogen
        # chemistry, not a crystallization salt.
        for smarts in ("C[I](=O)=O", "c1ccccc1[I](=O)=O"):
            self.assertFalse(_is_artifact(smarts), smarts)

    def test_select_fragments_applies_the_filter(self):
        frags_dict = {
            "ClOOOO": {"O=[Cl](=O)(=O)[O-]": ["LCP", "AAA", "BBB"]},
            "COO": {"CC(=O)O": ["AAA", "BBB", "CCC"]},
        }
        self.assertEqual(
            select_fragments(frags_dict, counts_threshold=2, max_size=5),
            ["CC(=O)O"],
        )


if __name__ == "__main__":
    unittest.main()


class ProtonationVariantCollapseTests(unittest.TestCase):
    """Charge-strict selection vs charge-loose SMARTS matching (docs/pitfalls.md)."""

    def test_charge_is_stripped_without_touching_aromaticity(self):
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            charge_normalized_fragment)
        self.assertEqual(charge_normalized_fragment('CC(=O)[O-]'), 'CC(=O)O')
        self.assertEqual(charge_normalized_fragment('CC(C)[N+]'), 'CC(C)N')
        self.assertEqual(charge_normalized_fragment('c[n+](c)C'), 'cn(c)C')
        self.assertEqual(charge_normalized_fragment('O=P([O-])([O-])O'), 'O=P(O)(O)O')
        # Uncharged fragments have no key of their own.
        self.assertIsNone(charge_normalized_fragment('CC(=O)O'))

    def test_aromatic_atoms_are_never_rewritten_as_aliphatic(self):
        """Fragment keys match exactly on aromaticity; nothing may bridge that.

        Round-tripping through RDKit sanitization would re-perceive aromaticity
        and turn `cC(c)N` into `CC(C)N` on any fragment too small to close a
        ring, merging chemically distinct fragments.
        """
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            charge_normalized_fragment)
        for aromatic in ('cC(c)N', 'cCCCN', 'nCCO', 'cc(c)n', 'c1cc[nH]c1'):
            self.assertIsNone(charge_normalized_fragment(aromatic), aromatic)
        # A charged aromatic keeps its aromatic case when neutralized.
        self.assertEqual(charge_normalized_fragment('cc(c)[n+]'), 'cc(c)n')

    def test_a_group_collapses_only_when_its_neutral_form_is_present(self):
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            group_protonation_variants)
        both = {'CC(=O)O': {'LIG1', 'LIG2'}, 'CC(=O)[O-]': {'LIG2', 'LIG3'}}
        reps, aliases = group_protonation_variants(both)
        self.assertEqual(aliases, {'CC(=O)[O-]': 'CC(=O)O'})
        # Ligand names pool, so the threshold sees prevalence rather than
        # which way the CCD happened to draw the acid.
        self.assertEqual(reps['CC(=O)O'], {'LIG1', 'LIG2', 'LIG3'})

        # Neutral form absent: promoting the charged one would change which
        # structures get mined, so nothing is merged.
        anion_only = {'CC(=O)[O-]': {'LIG2'}}
        reps, aliases = group_protonation_variants(anion_only)
        self.assertEqual(aliases, {})
        self.assertEqual(set(reps), {'CC(=O)[O-]'})

    def test_aliases_round_trip_through_the_library(self):
        import tempfile
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            write_fragment_aliases)
        from ligand_vdgs.functions.vdg_npz_utils import resolve_fragment_alias
        with tempfile.TemporaryDirectory() as lib:
            write_fragment_aliases(os.path.join(lib, 'fragment_aliases.tsv'),
                                   {'CC(=O)[O-]': 'CC(=O)O'})
            self.assertEqual(resolve_fragment_alias(lib, 'CC(=O)[O-]'), 'CC(=O)O')
            self.assertEqual(resolve_fragment_alias(lib, 'CC(=O)O'), 'CC(=O)O')
            self.assertEqual(resolve_fragment_alias(lib, 'CCO'), 'CCO')

    def test_missing_alias_file_is_not_an_error(self):
        import tempfile
        from ligand_vdgs.functions.vdg_npz_utils import resolve_fragment_alias
        with tempfile.TemporaryDirectory() as lib:
            self.assertEqual(resolve_fragment_alias(lib, 'CC(=O)[O-]'), 'CC(=O)[O-]')
