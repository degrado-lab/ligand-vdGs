"""Tests for the symmetry-aware Stage-1 internal-distance lower bound.

The bound is minimised over the vdG's *whole* symmetry group -- CG
automorphisms crossed with orderings of interchangeable vdM slots. Screening
over a subset of that group would make it inadmissible, so the admissibility
tests below cover both kinds of symmetry together.
"""

import itertools
import math
import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster as ac
from ligand_vdgs.functions.utils import kabsch_ssd
from ligand_vdgs.functions.vdg_fp_utils import (
    INTERNAL_DISTANCE_EPS,
    build_perm_group,
    full_row_permutations,
    internal_distance_lower_bounds,
    precompute_internal_distance_descriptors,
)


def _phosphate_permutations():
    """All four terminal oxygens exchange; phosphorus (slot 1) stays fixed."""
    out = []
    for oxygens in itertools.permutations((0, 2, 3, 4)):
        permutation = [0, 1, 2, 3, 4]
        for slot, source in zip((0, 2, 3, 4), oxygens):
            permutation[slot] = source
        out.append(tuple(permutation))
    return tuple(out)


def _group(cg_permutations, n_cg, n_res, same_label=False):
    """Symmetry group for a bucket of `n_res` slots, optionally interchangeable."""
    parts = ['bb'] * n_res if same_label else [f'AA{i}' for i in range(n_res)]
    return build_perm_group(cg_permutations, n_cg, parts)


def _exact_symmetry_rmsds(data, query_index, target_indices, n_cg, perm_group):
    """Exact Kabsch RMSD minimised over the same group the bound uses."""
    data = np.asarray(data, dtype=np.float32)
    targets = np.asarray(target_indices, dtype=np.intp)
    n_total = data.shape[1]
    best = np.full(len(targets), np.inf, dtype=np.float64)
    for full in full_row_permutations(perm_group, n_cg):
        ssd = kabsch_ssd(data[query_index, full], data[targets])
        np.minimum(best, np.sqrt(ssd / n_total), out=best)
    return best


class DescriptorTests(unittest.TestCase):
    def setUp(self):
        self.cg = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
        ], dtype=np.float32)
        self.bb1 = np.array([
            [4.0, 0.0, 0.0],
            [0.0, 0.0, 12.0],
            [0.0, 4.0, 3.0],
        ], dtype=np.float32)
        self.bb2 = np.array([
            [0.0, 0.0, 5.0],
            [0.0, -12.0, 0.0],
            [8.0, 0.0, 0.0],
        ], dtype=np.float32)

    def _translated_pair(self, record):
        shift = np.array([21.0, -7.0, 3.0], dtype=np.float32)
        return np.stack((record, record + shift))

    def test_subset_one_shape_dtype_and_manual_distances(self):
        data = self._translated_pair(np.concatenate((self.cg, self.bb1)))
        descriptors = precompute_internal_distance_descriptors(data, n_cg=2)

        self.assertEqual(set(descriptors), {'cg_bb'})
        self.assertEqual(descriptors['cg_bb'].shape, (2, 2, 3))
        self.assertEqual(descriptors['cg_bb'].dtype, np.float32)
        expected = np.array([
            [4.0, 12.0, 5.0],
            [5.0, math.sqrt(153.0), math.sqrt(10.0)],
        ], dtype=np.float32)
        np.testing.assert_allclose(descriptors['cg_bb'][0], expected, atol=1e-6)
        np.testing.assert_allclose(descriptors['cg_bb'][1], expected, atol=1e-6)

    def test_subset_two_stores_the_whole_backbone_matrix(self):
        record = np.concatenate((self.cg, self.bb1, self.bb2))
        data = self._translated_pair(record)
        descriptors = precompute_internal_distance_descriptors(data, n_cg=2)

        self.assertEqual(set(descriptors), {'cg_bb', 'bb_bb', 'cross_mask'})
        self.assertEqual(descriptors['cg_bb'].shape, (2, 2, 6))
        # Whole 6x6 matrix, not the nine cross distances alone: a slot
        # permutation indexes it on both axes.
        self.assertEqual(descriptors['bb_bb'].shape, (2, 6, 6))
        self.assertEqual(descriptors['bb_bb'].dtype, np.float32)
        self.assertEqual(int(descriptors['cross_mask'].sum()), 9)

        bb = record[2:]
        expected = np.linalg.norm(bb[:, None, :] - bb[None, :, :], axis=2)
        for record_index in (0, 1):
            np.testing.assert_allclose(
                descriptors['bb_bb'][record_index], expected, atol=1e-4)
            # Translation-invariant: both records give the same descriptor.
            np.testing.assert_allclose(descriptors['bb_bb'][record_index],
                                       descriptors['bb_bb'][0], atol=1e-6)

    def test_cross_mask_selects_only_inter_residue_pairs(self):
        from ligand_vdgs.functions.vdg_fp_utils import cross_residue_mask
        mask = cross_residue_mask(2)
        rows, cols = np.nonzero(mask)
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(r < 3 <= c for r, c in zip(rows, cols)))

    def test_invalid_coordinate_layouts_raise(self):
        invalid = (
            np.zeros((7, 3), dtype=np.float32),
            np.zeros((2, 7, 2), dtype=np.float32),
            np.zeros((2, 6, 3), dtype=np.float32),   # 2 CG + 4 backbone
            np.zeros((2, 11, 3), dtype=np.float32),  # 2 CG + 9 backbone
        )
        for data in invalid:
            with self.subTest(shape=data.shape):
                with self.assertRaises(ValueError):
                    precompute_internal_distance_descriptors(data, n_cg=2)


class AdmissibilityTests(unittest.TestCase):
    def test_randomized_bounds_do_not_exceed_exact_symmetry_rmsd(self):
        rng = np.random.default_rng(20260903)
        carboxylate = ((0, 1, 2, 3), (0, 1, 3, 2))
        cases = (
            (4, 1, ((0, 1, 2, 3),), False),
            (4, 2, carboxylate, False),
            (4, 2, carboxylate, True),                    # + slot swap
            (5, 1, _phosphate_permutations(), False),
            (5, 2, _phosphate_permutations(), False),
            (5, 2, _phosphate_permutations(), True),      # 24 x 2 group
        )
        for n_cg, n_res, permutations, same_label in cases:
            perm_group = _group(permutations, n_cg, n_res, same_label)
            with self.subTest(n_cg=n_cg, n_res=n_res,
                              n_group=len(perm_group), slots=same_label):
                n_total = n_cg + 3 * n_res
                query = rng.normal(scale=3.0, size=(1, n_total, 3))
                close = query + rng.normal(scale=0.3, size=(32, n_total, 3))
                far = rng.normal(scale=3.0, size=(32, n_total, 3))
                data = np.concatenate((query, close, far)).astype(np.float32)
                targets = np.arange(1, len(data), dtype=np.intp)
                descriptors = precompute_internal_distance_descriptors(data, n_cg)

                lower = internal_distance_lower_bounds(
                    0, targets, descriptors, n_total, perm_group)
                exact = _exact_symmetry_rmsds(data, 0, targets, n_cg, perm_group)

                self.assertTrue(np.all(lower <= exact + INTERNAL_DISTANCE_EPS),
                                (lower - exact).max())
                for cutoff in (0.5, 0.65, 1.0):
                    false_reject = ((lower > cutoff + INTERNAL_DISTANCE_EPS)
                                    & (exact <= cutoff))
                    self.assertFalse(np.any(false_reject))

    def test_ignoring_slot_symmetry_would_be_inadmissible(self):
        """The reason the group has to include slot orders, not just CG ones.

        A pair related by swapping two interchangeable vdM slots has exact RMSD
        0. Bounding it over CG automorphisms alone reports a large distance and
        would reject the pair outright.
        """
        rng = np.random.default_rng(11)
        n_cg, n_res = 4, 2
        n_total = n_cg + 3 * n_res
        base = rng.normal(scale=3.0, size=(n_total, 3)).astype(np.float32)
        swapped = base.copy()
        swapped[n_cg:n_cg + 3] = base[n_cg + 3:]
        swapped[n_cg + 3:] = base[n_cg:n_cg + 3]
        data = np.stack((base, swapped))
        descriptors = precompute_internal_distance_descriptors(data, n_cg)

        cg_only = _group(None, n_cg, n_res, same_label=False)
        with_slots = _group(None, n_cg, n_res, same_label=True)
        exact = _exact_symmetry_rmsds(data, 0, [1], n_cg, with_slots)[0]
        self.assertLess(exact, 1e-5)

        blind = internal_distance_lower_bounds(0, [1], descriptors, n_total, cg_only)[0]
        aware = internal_distance_lower_bounds(0, [1], descriptors, n_total, with_slots)[0]
        self.assertGreater(blind, 0.5)       # would reject a zero-RMSD pair
        self.assertAlmostEqual(aware, 0.0, places=6)


class SymmetryTests(unittest.TestCase):
    @staticmethod
    def _base_record(n_cg):
        cg = np.array([
            [2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ], dtype=np.float32)[:n_cg]
        bb = np.array([
            [7.0, 1.0, 0.0],
            [6.0, -2.0, 1.0],
            [8.0, 0.0, -1.0],
        ], dtype=np.float32)
        return np.concatenate((cg, bb))

    def _assert_permutation_is_recognized(self, base, n_cg, permutation,
                                          allowed_permutations):
        target = base.copy()
        target[:n_cg] = base[np.asarray(permutation, dtype=np.intp)]
        data = np.stack((base, target))
        descriptors = precompute_internal_distance_descriptors(data, n_cg)
        identity = _group((tuple(range(n_cg)),), n_cg, 1)
        allowed = _group(allowed_permutations, n_cg, 1)

        identity_bound = internal_distance_lower_bounds(
            0, [1], descriptors, len(base), identity)[0]
        symmetry_bound = internal_distance_lower_bounds(
            0, [1], descriptors, len(base), allowed)[0]
        exact = _exact_symmetry_rmsds(data, 0, [1], n_cg, allowed)[0]

        self.assertGreater(identity_bound, 0.05)
        self.assertAlmostEqual(symmetry_bound, 0.0, places=7)
        self.assertLessEqual(exact, 1e-6)

    def test_phosphate_four_oxygen_permutation(self):
        permutation = (4, 1, 0, 2, 3)
        self._assert_permutation_is_recognized(
            self._base_record(5), 5, permutation, _phosphate_permutations())

    def test_carboxylate_terminal_oxygen_swap(self):
        base = self._base_record(4)
        permutation = (0, 1, 3, 2)
        self._assert_permutation_is_recognized(
            base, 4, permutation, ((0, 1, 2, 3), permutation))


class Float32PaddingTests(unittest.TestCase):
    def test_one_ulp_boundary_bound_is_kept_for_exact_decision(self):
        # This collinear stretch makes the max-component lower bound tight.
        # Rounding the intended 0.65-A displacement through float32 puts both
        # the bound and exact RMSD just above 0.65 A; padding must retain it as
        # a candidate for the exact calculation rather than rejecting on noise.
        cutoff = 0.65
        displacement = np.float32(cutoff * math.sqrt(8.0))
        query = np.array([
            [-10.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ], dtype=np.float32)
        target = query.copy()
        target[0, 0] -= displacement / 2.0
        target[1, 0] += displacement / 2.0
        data = np.stack((query, target))
        descriptors = precompute_internal_distance_descriptors(data, n_cg=1)
        lower = internal_distance_lower_bounds(
            0, [1], descriptors, n_total=4)[0]

        self.assertGreater(lower, cutoff)
        self.assertLessEqual(lower, cutoff + INTERNAL_DISTANCE_EPS)

    def test_graph_keeps_a_pair_whose_bound_sits_inside_the_padding(self):
        cutoff = 0.5
        data = np.zeros((2, 4, 3), dtype=np.float32)
        perm_group = build_perm_group(None, 1, ['ALA'])
        qi, qj = ac.stage1_edges(data, cutoff, 1, perm_group)
        self.assertEqual(list(zip(qi.tolist(), qj.tolist())), [(0, 1)])


if __name__ == '__main__':
    unittest.main()
