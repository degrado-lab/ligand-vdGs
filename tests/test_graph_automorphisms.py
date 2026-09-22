import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from ligand_vdgs.functions import align_and_cluster
from ligand_vdgs.functions import utils

def _automorphisms(smarts):
    return utils.identify_mol_automorphisms(utils.mol_from_fragment(smarts))

class GraphAutomorphismTests(unittest.TestCase):
    def test_exact_maps_do_not_expand_to_independent_orbit_permutations(self):
        automorphisms = _automorphisms("CCCC")

        self.assertEqual(automorphisms, (
            (0, 1, 2, 3),
            (3, 2, 1, 0),
        ))
        self.assertNotIn((3, 1, 2, 0), automorphisms)

    def test_generic_graph_symmetry_is_not_limited_to_resonance_groups(self):
        self.assertEqual(len(_automorphisms("CC(C)C")), 6)
        self.assertEqual(len(_automorphisms("CS(=O)(=O)C")), 4)

    def test_every_returned_mapping_preserves_the_labeled_graph(self):
        mol = Chem.MolFromSmarts("c1ccccc1")
        atom_labels, adjacency = utils._automorphism_graph(mol)
        automorphisms = utils.identify_mol_automorphisms(mol)

        self.assertEqual(len(automorphisms), 12)
        for permutation in automorphisms:
            for i in range(mol.GetNumAtoms()):
                self.assertEqual(atom_labels[i], atom_labels[permutation[i]])
                for j in range(mol.GetNumAtoms()):
                    self.assertEqual(adjacency[i][j], adjacency[permutation[i]][permutation[j]])

    def test_resonance_normalization_keeps_intended_terminal_oxygen_symmetry(self):
        self.assertEqual(_automorphisms("CC(=O)[O-]"), (
            (0, 1, 2, 3),
            (0, 1, 3, 2),
        ))
        self.assertEqual(_automorphisms("CC(=O)O"), (
            (0, 1, 2, 3),
            (0, 1, 3, 2),
        ))

        phosphate = _automorphisms("COP(=O)([O-])[O-]")
        self.assertEqual(len(phosphate), 6)
        self.assertTrue(all(p[1:3] == (1, 2) for p in phosphate))
        self.assertEqual({p[3:] for p in phosphate}, {
            (3, 4, 5), (3, 5, 4), (4, 3, 5),
            (4, 5, 3), (5, 3, 4), (5, 4, 3),
        })
        protonated_phosphate = _automorphisms("COP(=O)([O-])O")
        self.assertEqual(len(protonated_phosphate), 6)
        self.assertEqual({p[3:] for p in protonated_phosphate}, {
            (3, 4, 5), (3, 5, 4), (4, 3, 5),
            (4, 5, 3), (5, 3, 4), (5, 4, 3),
        })
        self.assertEqual(_automorphisms("COP(=O)O"), (
            (0, 1, 2, 3, 4),
            (0, 1, 2, 4, 3),
        ))

        phosphorimidate = _automorphisms("N=P([O-])(O)O")
        self.assertEqual(len(phosphorimidate), 6)
        self.assertTrue(all(p[:2] == (0, 1) for p in phosphorimidate))
        self.assertEqual({p[2:] for p in phosphorimidate}, {
            (2, 3, 4), (2, 4, 3), (3, 2, 4),
            (3, 4, 2), (4, 2, 3), (4, 3, 2),
        })

        self.assertEqual(len(_automorphisms("O=[Cl](=O)(=O)[O-]")), 24)
        self.assertEqual(len(_automorphisms("O=I(=O)(=O)[O-]")), 24)
        self.assertEqual(len(_automorphisms("O=[Cl](=O)[O-]")), 6)
        self.assertEqual(len(_automorphisms("cB(=O)O")), 2)
        self.assertEqual(len(_automorphisms("CB(O)(O)[O+]")), 6)

    def test_phosphate_drawings_all_merge_all_four_oxygens(self):
        for fragment in ("O=P(O)(O)O", "O=P([O-])(O)O",
                         "O=P([O-])([O-])O", "O=P([O-])([O-])[O-]"):
            self.assertEqual(len(_automorphisms(fragment)), 24, fragment)

    def test_saturated_center_with_mixed_charges_keeps_its_drawing(self):
        bisulfate = _automorphisms("[O-]S([O-])([O-])O")
        self.assertEqual(len(bisulfate), 6)
        self.assertTrue(all(p[1] == 1 and p[4] == 4 for p in bisulfate))
        self.assertEqual({(p[0], p[2], p[3]) for p in bisulfate}, {
            (0, 2, 3), (0, 3, 2), (2, 0, 3),
            (2, 3, 0), (3, 0, 2), (3, 2, 0),
        })

        self.assertEqual(len(_automorphisms("cS([O-])(O)O")), 6)
        self.assertEqual(len(_automorphisms("O=S(=O)([O-])O")), 24)
        self.assertEqual(len(_automorphisms("O=P([O-])([O-])[O-]")), 24)
        self.assertEqual(len(_automorphisms("OS(O)(O)O")), 24)
        self.assertEqual(len(_automorphisms("O[B-](O)(O)O")), 24)

    def test_terminal_sulfur_exchanges_only_with_terminal_sulfur(self):
        self.assertEqual(_automorphisms("CC(=S)S"), (
            (0, 1, 2, 3),
            (0, 1, 3, 2),
        ))
        self.assertEqual(len(_automorphisms("cC(=S)S")), 2)
        self.assertEqual(len(_automorphisms("S=C(S)S")), 6)

        thiophosphate = _automorphisms("OP(O)(=S)[S-]")
        self.assertEqual(len(thiophosphate), 4)
        self.assertEqual({(p[0], p[2]) for p in thiophosphate}, {(0, 2), (2, 0)})
        self.assertEqual({(p[3], p[4]) for p in thiophosphate}, {(3, 4), (4, 3)})

        self.assertEqual(_automorphisms("CC(=O)S"), ((0, 1, 2, 3),))

    def test_terminal_group_ignores_drawn_bond_order_and_charge(self):
        self.assertEqual(_automorphisms("CCS(=O)O"), (
            (0, 1, 2, 3, 4),
            (0, 1, 2, 4, 3),
        ))
        self.assertEqual(len(_automorphisms("cS(=O)O")), 2)

        self.assertEqual(len(_automorphisms("cS([O-])(O)O")), 6)

        self.assertEqual(_automorphisms("[N+]C(N)=O"), (
            (0, 1, 2, 3),
            (2, 1, 0, 3),
        ))

    def test_aromatic_terminal_atoms_are_not_resonance_forms(self):
        self.assertEqual(_automorphisms("cc(o)O"), ((0, 1, 2, 3),))
        self.assertEqual(_automorphisms("cc(S)s"), ((0, 1, 2, 3),))
        self.assertEqual(len(_automorphisms("NS(n)(=O)=O")), 2)

    def test_permissive_phosphorus_fragment_terminal_policy(self):
        self.assertEqual(_automorphisms("CCP(=O)O"), (
            (0, 1, 2, 3, 4),
            (0, 1, 2, 4, 3),
        ))
        self.assertEqual(_automorphisms("CP(=O)(O)F"), (
            (0, 1, 2, 3, 4),
            (0, 1, 3, 2, 4),
        ))

        self.assertEqual(set(_automorphisms("CP(C)(=O)O")), {
            (0, 1, 2, 3, 4),
            (0, 1, 2, 4, 3),
            (2, 1, 0, 3, 4),
            (2, 1, 0, 4, 3),
        })

        thiophosphate = _automorphisms("[O-]P(O)(O)=S")
        self.assertEqual(len(thiophosphate), 6)
        self.assertTrue(all(p[1] == 1 and p[4] == 4 for p in thiophosphate))
        self.assertEqual({(p[0], p[2], p[3]) for p in thiophosphate}, {
            (0, 2, 3), (0, 3, 2), (2, 0, 3),
            (2, 3, 0), (3, 0, 2), (3, 2, 0),
        })

        self.assertEqual(_automorphisms("CP(=O)(OC)F"), (
            (0, 1, 2, 3, 4, 5),
        ))

    def test_amidine_and_guanidine_resonance_normalization(self):
        expected_amidine = (
            (0, 1, 2, 3),
            (0, 1, 3, 2))
        for smarts in ("CC(=N)N", "CC(=[NH2+])N"):
            with self.subTest(smarts=smarts):
                self.assertEqual(_automorphisms(smarts), expected_amidine)

        expected_three_n_permutations = {
            (0, 2, 3), (0, 3, 2), (2, 0, 3),
            (2, 3, 0), (3, 0, 2), (3, 2, 0)}
        for smarts in ("N=C(N)N", "NC(=[NH2+])N"):
            with self.subTest(smarts=smarts):
                automorphisms = _automorphisms(smarts)
                self.assertEqual(len(automorphisms), 6)
                self.assertTrue(all(p[1] == 1 for p in automorphisms))
                self.assertEqual(
                    {(p[0], p[2], p[3]) for p in automorphisms},
                    expected_three_n_permutations)

        self.assertEqual(_automorphisms("CNC(=[NH2+])N"), (
            (0, 1, 2, 3, 4),
            (0, 1, 2, 4, 3),
        ))
        self.assertEqual(_automorphisms("CC(=N)NC"), (
            (0, 1, 2, 3, 4),
        ))

    def test_aromatic_n_resonance_is_component_scoped(self):
        expected_counts = {
            "Cc([n+])n": 2,
            "Cc([n-])n": 2,
            "[n+]c([n+])n": 6,
            "[n+]c(n)N": 2,
            "[n+]c(n)O": 2,
            "[n+]c(n)n": 6,
            "[n+]cncn": 2,
            "c[n-]c[n+]c": 2,
            "cc([n+])n": 2,
            "cn([n+])n": 2,
            "cnc[n+]c": 2,
            "nc([n+])F": 2,
            "nc([n+])s": 2}
        for smarts, expected_count in expected_counts.items():
            with self.subTest(smarts=smarts):
                self.assertEqual(len(_automorphisms(smarts)), expected_count)

        self.assertEqual(utils._find_resonance_aromatic_N_atoms(
            utils.mol_from_fragment("Cc([n+])n")), {2, 3})
        self.assertEqual(_automorphisms("c([nH])n"), (
            (0, 1, 2),
            (0, 2, 1),
        ))

        self.assertEqual(_automorphisms("[n+]CCCn"), (
            (0, 1, 2, 3, 4),
        ))
        self.assertEqual(_automorphisms("NCC[N+]"), (
            (0, 1, 2, 3),
        ))

    def test_n_o_resonance_normalization(self):
        expected = (
            (0, 1, 2, 3),
            (0, 1, 3, 2),
        )
        for smarts in ("C=[N+]([O-])O", "cN(=O)O"):
            with self.subTest(smarts=smarts):
                self.assertEqual(_automorphisms(smarts), expected)

        self.assertEqual(_automorphisms("C[N+](=O)OC"), (
            (0, 1, 2, 3, 4),
        ))

    def test_s_n_resonance_normalization(self):
        self.assertEqual(_automorphisms("CS(=N)(N)=O"), (
            (0, 1, 2, 3, 4),
            (0, 1, 3, 2, 4),
        ))
        self.assertEqual(_automorphisms("N=S(N)(=O)F"), (
            (0, 1, 2, 3, 4),
            (2, 1, 0, 3, 4),
        ))
        self.assertEqual(_automorphisms("cS(=N)(N)=O"), (
            (0, 1, 2, 3, 4),
            (0, 1, 3, 2, 4),
        ))

    def test_conjugated_n_chain_resonance_normalization(self):
        self.assertEqual(_automorphisms("cN=NNc"), (
            (0, 1, 2, 3, 4),
            (4, 3, 2, 1, 0),
        ))
        self.assertEqual(_automorphisms("cN=NNC"), (
            (0, 1, 2, 3, 4),
        ))

    def test_cn_rule_keeps_a_different_third_element_fixed(self):
        self.assertEqual(_automorphisms("N=C(N)n"), (
            (0, 1, 2, 3),
            (2, 1, 0, 3),
        ))
        self.assertEqual(_automorphisms("NC(=[N+])S"), (
            (0, 1, 2, 3),
            (2, 1, 0, 3),
        ))

    def test_cn_normalization_boundary_keeps_ordinary_graph_symmetry(self):
        self.assertEqual(_automorphisms("NC(=O)N"), (
            (0, 1, 2, 3),
            (3, 1, 2, 0),
        ))
        self.assertEqual(_automorphisms("NCN"), (
            (0, 1, 2),
            (2, 1, 0),
        ))

        charged_aminal = utils.mol_from_fragment("[NH3+]C[NH2]")
        self.assertEqual(utils._find_resonance_CN_groups(charged_aminal), {})
        self.assertEqual(utils.identify_mol_automorphisms(charged_aminal), (
            (0, 1, 2),
            (2, 1, 0),
        ))

        self.assertEqual(
            utils._find_resonance_CN_groups(
                utils.mol_from_fragment("N#CN")),
            {},
        )
        self.assertEqual(
            len(utils._find_resonance_CN_groups(
                utils.mol_from_fragment("[C+](N)N"))),
            2,
        )

    def test_cn_normalization_recovers_resonance_swap_rmsd(self):
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.2, 0.0],
            [0.1, 1.4, 0.3],
            [0.3, 0.4, 2.0],
        ], dtype=np.float32)
        swapped = coords[[0, 1, 3, 2]]

        identity_only_rmsd = align_and_cluster._permuted_rmsd_pair(
            coords, swapped, ((0, 1, 2, 3),), n_cg=4
        )
        normalized_rmsd = align_and_cluster._permuted_rmsd_pair(
            coords, swapped, _automorphisms("CC(=N)N"), n_cg=4
        )

        self.assertGreater(identity_only_rmsd, 0.8)
        self.assertAlmostEqual(normalized_rmsd, 0.0, places=6)

    def test_explicit_automorphisms_prevent_false_low_rmsd(self):
        x = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.2, 0.0],
            [0.1, 1.4, 0.3],
            [0.3, 0.4, 2.0],
        ], dtype=np.float32)
        y = x[[3, 1, 2, 0]]

        rmsd = align_and_cluster._permuted_rmsd_pair(
            x, y, _automorphisms("CCCC"), n_cg=4
        )

        self.assertAlmostEqual(rmsd, 0.8134896, places=6)
        self.assertGreater(rmsd, 0.5)

    def test_fragment_parser_is_unsanitized_smarts_only(self):
        mol = utils.mol_from_fragment("CC(=O)O")

        self.assertIsNotNone(mol)
        self.assertTrue(all(atom.HasQuery() for atom in mol.GetAtoms()))

    def test_inplace_ligand_rmsd_matches_rdkit_calcrms(self):
        ref = Chem.MolFromSmiles("CCCC")
        query = Chem.Mol(ref)
        ref_conf = Chem.Conformer(4)
        query_conf = Chem.Conformer(4)
        for atom_idx, xyz in enumerate([
            (0.0, 0.0, 0.0),
            (1.0, 0.2, 0.0),
            (0.1, 1.4, 0.3),
            (0.3, 0.4, 2.0),
        ]):
            ref_conf.SetAtomPosition(atom_idx, xyz)
            query_conf.SetAtomPosition(atom_idx, xyz)
        ref.AddConformer(ref_conf)
        query.AddConformer(query_conf)

        expected = rdMolAlign.CalcRMS(query, ref, maxMatches=0)
        observed = utils.best_inplace_symmetry_rmsd(ref, query)

        self.assertAlmostEqual(observed, expected)

    def test_inplace_ligand_rmsd_search_is_not_truncated(self):
        ref = Chem.MolFromSmiles("CCCC")
        query = Chem.Mol(ref)
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.2, 0.0],
            [0.1, 1.4, 0.3],
            [0.3, 0.4, 2.0],
        ])
        for mol, atom_coords in ((ref, coords), (query, coords[::-1])):
            conf = Chem.Conformer(4)
            for atom_idx, xyz in enumerate(atom_coords):
                conf.SetAtomPosition(atom_idx, xyz)
            mol.AddConformer(conf)

        self.assertGreater(
            rdMolAlign.CalcRMS(query, ref, maxMatches=1), 1.0)
        self.assertAlmostEqual(
            utils.best_inplace_symmetry_rmsd(ref, query), 0.0)

    def test_inplace_ligand_rmsd_allows_carboxyl_oxygen_exchange(self):
        for smiles in ("CC(=O)O", "CC(=O)[O-]"):
            with self.subTest(smiles=smiles):
                ref = Chem.MolFromSmiles(smiles)
                query = Chem.Mol(ref)
                coords = np.array([
                    [0.0, 0.0, 0.0],
                    [1.1, 0.3, 0.2],
                    [-0.4, 1.5, 0.7],
                    [2.3, -0.8, 1.1],
                ])
                for mol, atom_coords in (
                        (ref, coords), (query, coords[[0, 1, 3, 2]])):
                    conf = Chem.Conformer(4)
                    for atom_idx, xyz in enumerate(atom_coords):
                        conf.SetAtomPosition(atom_idx, xyz)
                    mol.AddConformer(conf)

                self.assertAlmostEqual(
                    utils.best_inplace_symmetry_rmsd(ref, query), 0.0)

    def test_inplace_ligand_rmsd_allows_terminal_phosphate_oh_exchange(self):
        ref = Chem.MolFromSmiles("COP(=O)([O-])O")
        query = Chem.Mol(ref)
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.2, 0.1, 0.3],
            [-0.4, 1.7, 0.2],
            [0.3, -0.6, 2.1],
            [2.4, 1.1, -0.8],
            [-1.3, 0.7, 1.5],
        ])
        for mol, atom_coords in ((ref, coords), (query, coords[[0, 1, 2, 5, 4, 3]])):
            conf = Chem.Conformer(mol.GetNumAtoms())
            for atom_idx, xyz in enumerate(atom_coords):
                conf.SetAtomPosition(atom_idx, xyz)
            mol.AddConformer(conf)

        self.assertAlmostEqual(utils.best_inplace_symmetry_rmsd(ref, query), 0.0)

        translated = Chem.Mol(ref)
        translated_conf = translated.GetConformer()
        for atom_idx in range(translated.GetNumAtoms()):
            point = translated_conf.GetAtomPosition(atom_idx)
            translated_conf.SetAtomPosition(
                atom_idx, (point.x + 5.0, point.y, point.z)
            )
        self.assertAlmostEqual(
            utils.best_inplace_symmetry_rmsd(ref, translated), 5.0
        )

    def test_explicit_hydrogen_atom_nodes_are_rejected(self):
        methane_with_hydrogens = Chem.AddHs(Chem.MolFromSmiles("C"))

        with self.assertRaisesRegex(ValueError, "hydrogen-free molecular graph"):
            utils.identify_mol_automorphisms(methane_with_hydrogens)
        with self.assertRaisesRegex(ValueError, "hydrogen-free molecular graph"):
            utils.best_inplace_symmetry_rmsd(
                methane_with_hydrogens, Chem.Mol(methane_with_hydrogens)
            )

    def test_automorphism_limit_fails_instead_of_returning_an_incomplete_set(self):
        with self.assertRaisesRegex(ValueError, "more than max_automorphisms"):
            utils.identify_mol_automorphisms(
                Chem.MolFromSmarts("c1ccccc1"), max_automorphisms=5
            )

class DeclaredAnnotationTests(unittest.TestCase):
    def test_declared_degree_keeps_bridging_atoms_out_of_terminal_groups(self):
        self.assertEqual(
            len(_automorphisms("[O;!R][P;!R]([O;!R])([O;!R])[O;!R]")), 24)
        self.assertEqual(
            len(_automorphisms(
                "[O;D1;!R][P;!R]([O;D1;!R])([O;D1;!R])[O;D2;!R]")), 6)
        self.assertEqual(
            len(_automorphisms(
                "[O;D1;!R][P;!R]([O;D1;!R])([O;D2;!R])[O;D2;!R]")), 4)
        self.assertEqual(
            len(_automorphisms(
                "[O;D1;!R][P;!R]([O;D2;!R])([O;D2;!R])[O;D2;!R]")), 6)
        self.assertEqual(
            len(_automorphisms("[O;D1;!R][P;!R]([O;D1;!R])[O;D2;!R]")), 2)

    def test_declared_degree_survives_resonance_label_nulling(self):
        self.assertEqual(
            len(_automorphisms("[O;D1;!R][C;!R](=[O;D1;!R])[O;D1;!R]")), 6)
        self.assertEqual(
            len(_automorphisms("[O;D1;!R][C;!R](=[O;D1;!R])[O;D2;!R]")), 2)
        self.assertEqual(
            len(_automorphisms("[N;D1;!R][C;!R](=[N;D1;!R])[N;D1;!R]")), 6)
        self.assertEqual(
            len(_automorphisms("[N;D1;!R][C;!R](=[N;D1;!R])[N;D2;!R]")), 2)

    def test_every_returned_mapping_preserves_declared_degree_and_h(self):
        for smarts in (
            "[O;D1;!R][P;!R]([O;D1;!R])([O;D1;!R])[O;D2;!R]",
            "[O;D1;!R][P;!R]([O;D1;!R])([O;D2;!R])[O;D2;!R]",
            "[N;D1;!R][C;!R](=[N;D1;!R])[N;D2;!R]",
            "[c;!H0]1[c;!H0][c;!H0][c;H0][c;!H0][c;!H0]1",
        ):
            mol = utils.mol_from_fragment(smarts)
            declared = [
                (utils._query_primitive(atom, utils._QUERY_DEGREE_LINE),
                 utils._query_primitive(atom, utils._QUERY_HCOUNT_LINE))
                for atom in mol.GetAtoms()
            ]
            for permutation in utils.identify_mol_automorphisms(mol):
                for i, image in enumerate(permutation):
                    self.assertEqual(declared[i], declared[image], smarts)

    def test_carbon_h_flag_splits_otherwise_identical_positions(self):
        self.assertEqual(
            len(_automorphisms("[c;!H0]1[c;!H0][c;!H0][c;!H0][c;!H0][c;!H0]1")),
            12)
        self.assertEqual(
            len(_automorphisms("[c;!H0]1[c;!H0][c;!H0][c;H0][c;!H0][c;!H0]1")), 2)
        self.assertEqual(
            len(_automorphisms("[C;!R;!H0][O;D2;!R][C;!R;!H0]")), 2)
        self.assertEqual(
            len(_automorphisms("[C;!R;!H0][O;D2;!R][C;!R;H0]")), 1)

    def test_declared_h_on_heteroatoms_stays_normalized(self):
        self.assertEqual(
            len(_automorphisms("[O;D1;!R;H1][C;!R](=[O;D1;!R;H0])[C;!R]")), 2)
        self.assertEqual(
            len(_automorphisms("[N;D2;!R;H1][C;!R](=[N;D1;!R;H0])[N;D1;!R;H0]")),
            2)

class ExtractElementsTests(unittest.TestCase):
    def _assert_matches_rdkit(self, smiles):
        mol = utils.mol_from_fragment(smiles)
        self.assertIsNotNone(mol, smiles)
        self.assertEqual([e.capitalize() for e in utils.extract_elements(smiles)],
                         [a.GetSymbol() for a in mol.GetAtoms()], smiles)

    def test_bracketed_and_aromatic_atoms(self):
        for smiles in ("C[Se]C", "C[As](C)C", "cc[nH]c", "bc(c)c", "cb(o)O",
                       "Br[13CH2]Cl", "C[C@@H](N)C(=O)[O-]", "O=P([O-])(O)O",
                       "c1ccc2c(c1)[nH]cn2"):
            self._assert_matches_rdkit(smiles)

if __name__ == "__main__":
    unittest.main()
