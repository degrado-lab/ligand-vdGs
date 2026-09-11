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

        # Neutral form absent and only one charged spelling: nothing to pool, so
        # the key that exists in the dict stays the representative (a promoted
        # key needs >= 2 spellings; see PromotedKeyTests).
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
                                   {'CC(=O)[O-]': 'CC(=O)O'}, {'CC(=O)O', 'CC(=O)[O-]'})
            self.assertEqual(resolve_fragment_alias(lib, 'CC(=O)[O-]'), 'CC(=O)O')
            self.assertEqual(resolve_fragment_alias(lib, 'CC(=O)O'), 'CC(=O)O')
            self.assertEqual(resolve_fragment_alias(lib, 'CCO'), 'CCO')

    def test_missing_alias_file_is_not_an_error(self):
        import tempfile
        from ligand_vdgs.functions.vdg_npz_utils import resolve_fragment_alias
        with tempfile.TemporaryDirectory() as lib:
            self.assertEqual(resolve_fragment_alias(lib, 'CC(=O)[O-]'), 'CC(=O)[O-]')


class PromotedKeyTests(unittest.TestCase):
    """Two charged spellings with no neutral twin pool under the charge-stripped key.

    Aromatic nitro is the real case: RDKit writes it as `[N+](=O)[O-]` or
    `[N+](=O)O` and never neutral, so the 2026-09-06 build mined both, one nested
    in the other.
    """
    NITRO_A = 'c[N+;!R](=[O;!R])[O-;!R]'
    NITRO_B = 'c[N+;!R](=[O;!R])[O;!R]'
    NITRO_PROMOTED = 'c[N;!R](=[O;!R])[O;!R]'

    def test_no_twin_pair_is_promoted_and_pooled(self):
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            group_protonation_variants)
        reps, aliases = group_protonation_variants(
            {self.NITRO_A: {'NIT', 'DNP'}, self.NITRO_B: {'DNP', 'PNP'}})
        self.assertEqual(set(reps), {self.NITRO_PROMOTED})
        self.assertEqual(reps[self.NITRO_PROMOTED], {'NIT', 'DNP', 'PNP'})
        self.assertEqual(aliases, {self.NITRO_A: self.NITRO_PROMOTED,
                                   self.NITRO_B: self.NITRO_PROMOTED})

    def test_differently_ordered_spellings_still_meet(self):
        """Grouping compares queries, not strings: the second spelling is written
        with the ring atom last, so its charge-stripped text differs from the first's."""
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            group_protonation_variants)
        from ligand_vdgs.functions.utils import fragment_keys_equivalent
        a, b = '[C;r5][N;!R]=[N+;!R]=[N-;!R]', '[N-;!R]=[N+;!R]=[N;!R][C;r5]'
        reps, aliases = group_protonation_variants({a: {'AZ1'}, b: {'AZ2'}})
        (rep,) = reps
        self.assertEqual(reps[rep], {'AZ1', 'AZ2'})
        self.assertEqual(set(aliases), {a, b})
        # The promoted key keeps the ring annotation and is the members' loose form.
        self.assertIn(';r5]', rep)
        self.assertTrue(fragment_keys_equivalent(rep, '[C;r5][N;!R]=[N;!R]=[N;!R]'))

    def test_twin_beats_promotion_and_lone_spelling_is_untouched(self):
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            group_protonation_variants)
        # Twin present (in a different atom order): it is the representative, not
        # a promoted key, and the twin's own ligands pool in.
        qualifying = {'[C;!R][C;!R](=[O;!R])[O-;!R]': {'L1'},
                      '[O-;!R]C(=[O;!R])[C;!R]'.replace('C(', '[C;!R]('): {'L2'},
                      '[O;!R]=[C;!R]([O;!R])[C;!R]': {'L3'},
                      'C[N+;!R](C)(C)C': {'Q1'}}
        reps, aliases = group_protonation_variants(qualifying)
        self.assertEqual(set(reps), {'[O;!R]=[C;!R]([O;!R])[C;!R]', 'C[N+;!R](C)(C)C'})
        self.assertEqual(reps['[O;!R]=[C;!R]([O;!R])[C;!R]'], {'L1', 'L2', 'L3'})
        self.assertEqual(set(aliases.values()), {'[O;!R]=[C;!R]([O;!R])[C;!R]'})
        self.assertEqual(reps['C[N+;!R](C)(C)C'], {'Q1'})

    def test_promoted_key_parses_and_covers_each_member(self):
        from ligand_vdgs.functions.utils import (fragment_key_query_mol,
                                                 fragment_query_mols_equivalent)
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            charge_normalized_fragment)
        promoted = fragment_key_query_mol(self.NITRO_PROMOTED)
        self.assertIsNotNone(promoted)
        for member in (self.NITRO_A, self.NITRO_B):
            self.assertTrue(fragment_query_mols_equivalent(
                promoted, fragment_key_query_mol(charge_normalized_fragment(member))))

    def test_resolve_fragment_key_finds_variants_and_reorderings(self):
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            resolve_fragment_key)
        candidates = [self.NITRO_A, self.NITRO_B, '[C;!R][O;!R][C;!R][C;!R]', 'cC(=O)O']
        self.assertEqual(resolve_fragment_key(self.NITRO_PROMOTED, candidates),
                         [(self.NITRO_A, 'charge_variant'),
                          (self.NITRO_B, 'charge_variant')])
        self.assertEqual(resolve_fragment_key('[C;!R][C;!R][O;!R][C;!R]', candidates),
                         [('[C;!R][O;!R][C;!R][C;!R]', 'equivalent')])
        self.assertEqual(resolve_fragment_key(self.NITRO_A, candidates),
                         [(self.NITRO_A, 'exact'), (self.NITRO_B, 'charge_variant')])
        self.assertEqual(resolve_fragment_key('CCCC', candidates), [])

    def test_alias_file_records_kind_and_still_loads(self):
        import tempfile
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            write_fragment_aliases)
        from ligand_vdgs.functions.vdg_npz_utils import load_fragment_aliases
        aliases = {self.NITRO_A: self.NITRO_PROMOTED, self.NITRO_B: self.NITRO_PROMOTED,
                   'CC(=O)[O-]': 'CC(=O)O'}
        dict_keys = {self.NITRO_A, self.NITRO_B, 'CC(=O)[O-]', 'CC(=O)O'}
        with tempfile.TemporaryDirectory() as lib:
            path = os.path.join(lib, 'fragment_aliases.tsv')
            write_fragment_aliases(path, aliases, dict_keys)
            rows = [l.rstrip('\n').split('\t') for l in open(path) if not l.startswith('#')]
            kinds = {r[0]: r[2] for r in rows}
            self.assertEqual(kinds, {self.NITRO_A: 'promoted', self.NITRO_B: 'promoted',
                                     'CC(=O)[O-]': 'neutral_twin'})
            # Third column must not leak into the representative.
            self.assertEqual(load_fragment_aliases(lib), aliases)


if __name__ == "__main__":
    unittest.main()
