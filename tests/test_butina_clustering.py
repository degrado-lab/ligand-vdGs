"""Tests for Stage-1 sphere-exclusion (Butina) clustering.

The properties asserted here are the ones the previous moving-medoid leader
pass did not have: every member lies within the cutoff of the representative
that is actually stored, every record lands in exactly one cluster, and the
result does not depend on input order.
"""

import itertools
import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster as ac
from ligand_vdgs.functions.utils import kabsch_ssd
from ligand_vdgs.functions.vdg_fp_utils import (
    build_perm_group, full_row_permutations)


def _phosphate_permutations():
    out = []
    for oxygens in itertools.permutations((0, 2, 3, 4)):
        permutation = [0, 1, 2, 3, 4]
        for slot, source in zip((0, 2, 3, 4), oxygens):
            permutation[slot] = source
        out.append(tuple(permutation))
    return tuple(out)


def _brute_force_distances(data, n_cg, perm_group):
    """Dense symmetry-minimised RMSD matrix, the reference the graph must match."""
    data = np.asarray(data, dtype=np.float32)
    n, n_total = len(data), data.shape[1]
    full_perms = full_row_permutations(perm_group, n_cg)
    out = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        stacked = np.concatenate([data[:, p] for p in full_perms], axis=0)
        repeated = np.empty_like(stacked)
        for k in range(len(full_perms)):
            repeated[k * n:(k + 1) * n] = data[i]
        ssd = kabsch_ssd(repeated, stacked).reshape(len(full_perms), n).min(axis=0)
        out[i] = np.sqrt(ssd / n_total)
    return np.minimum(out, out.T)


def _clustered_blob(rng, n_cg, n_res, n_centers=5, per_center=9, spread=0.22):
    n_total = n_cg + 3 * n_res
    centers = rng.normal(scale=4.0, size=(n_centers, n_total, 3))
    return np.concatenate([
        center + rng.normal(scale=spread, size=(per_center, n_total, 3))
        for center in centers]).astype(np.float32)


def _edge_set(qi, qj):
    return set(zip(qi.tolist(), qj.tolist()))


def _partition(clusters, remap=None):
    return {frozenset(int(remap[m]) if remap is not None else int(m) for m in members)
            for members in clusters}


CASES = (
    ('asymmetric CG, 1 slot',   4, ['ALA'],        None),
    ('carboxylate, 1 slot',     4, ['ALA'],        ((0, 1, 2, 3), (0, 1, 3, 2))),
    ('carboxylate, 2 slots',    4, ['ARG', 'bb'],  ((0, 1, 2, 3), (0, 1, 3, 2))),
    ('carboxylate, like slots', 4, ['bb', 'bb'],   ((0, 1, 2, 3), (0, 1, 3, 2))),
    ('phosphate, like slots',   5, ['bb', 'bb'],   _phosphate_permutations()),
)


class ButinaGuaranteeTests(unittest.TestCase):
    def test_every_member_is_within_the_cutoff_of_its_representative(self):
        rng = np.random.default_rng(4)
        cutoff = 0.6
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                group = build_perm_group(perms, n_cg, parts)
                exact = _brute_force_distances(data, n_cg, group)
                clusters = ac.get_butina_clusters(data, cutoff, n_cg, build_perm_group(perms, n_cg, parts))
                for members in clusters:
                    # The representative is emitted first.
                    worst = exact[members, members[0]].max()
                    self.assertLessEqual(worst, cutoff + 1e-6)

    def test_every_record_belongs_to_exactly_one_cluster(self):
        rng = np.random.default_rng(5)
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                clusters = ac.get_butina_clusters(data, 0.6, n_cg, build_perm_group(perms, n_cg, parts))
                assigned = sorted(i for m in clusters for i in m.tolist())
                self.assertEqual(assigned, list(range(len(data))))

    def test_neighbour_graph_is_exactly_the_within_cutoff_graph(self):
        """The prefilter cascade may only skip work, never drop a true pair."""
        rng = np.random.default_rng(6)
        cutoff = 0.6
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                group = build_perm_group(perms, n_cg, parts)
                exact = _brute_force_distances(data, n_cg, group)
                got = _edge_set(*ac.stage1_edges(data, cutoff, n_cg, group))
                want = {(i, j) for i in range(len(data))
                        for j in range(i + 1, len(data)) if exact[i, j] <= cutoff}
                self.assertEqual(got, want)

    def test_boundary_pairs_are_kept_under_every_group_element(self):
        """The discriminating case for per-element pruning: a pair within the
        cutoff only under a non-identity element, whose bound under the
        identity exceeds the cutoff. Failure looks like a dropped edge."""
        rng = np.random.default_rng(12)
        n_cg, parts = 5, ['bb', 'bb']
        base = rng.normal(scale=3.0, size=(n_cg + 6, 3)).astype(np.float32)
        # Relabel two oxygens AND swap the two slots; the identity element sees
        # a large RMSD, the combined element sees ~0.
        other = base[[0, 1, 3, 2, 4, 8, 9, 10, 5, 6, 7]] + \
            rng.normal(scale=0.02, size=base.shape).astype(np.float32)
        data = np.stack([base, other])
        group = build_perm_group(_phosphate_permutations(), n_cg, parts)
        exact = _brute_force_distances(data, n_cg, group)
        self.assertLess(exact[0, 1], 0.1)
        self.assertEqual(_edge_set(*ac.stage1_edges(data, 0.5, n_cg, group)), {(0, 1)})
        identity_only = build_perm_group(None, n_cg, ['ARG', 'bb'])
        self.assertEqual(_edge_set(*ac.stage1_edges(data, 0.5, n_cg, identity_only)), set())


class ButinaDeterminismTests(unittest.TestCase):
    def test_result_is_invariant_to_input_order(self):
        """Cluster counts are compared across fragments, so order must not matter."""
        rng = np.random.default_rng(7)
        n_cg, parts = 4, ['ARG']
        data = _clustered_blob(rng, n_cg, 1)
        order = rng.permutation(len(data))
        base = ac.get_butina_clusters(data, 0.6, n_cg, build_perm_group(None, n_cg, parts))
        shuffled = ac.get_butina_clusters(data[order], 0.6, n_cg, build_perm_group(None, n_cg, parts))
        self.assertEqual(_partition(base), _partition(shuffled, remap=order))

    def test_singleton_and_empty_inputs(self):
        self.assertEqual(ac.get_butina_clusters(
            np.zeros((0, 7, 3), np.float32), 0.5, 4, build_perm_group(None, 4, ['ALA'])), [])
        clusters = ac.get_butina_clusters(
            np.zeros((1, 7, 3), np.float32), 0.5, 4, build_perm_group(None, 4, ['ALA']))
        self.assertEqual([m.tolist() for m in clusters], [[0]])

    def test_empty_and_complete_graphs(self):
        n = 6
        indptr, indices = ac.neighbor_csr(n, np.empty(0, np.int32), np.empty(0, np.int32))
        self.assertEqual(len(ac.butina_partition(indptr, indices)), n)
        ii, jj = np.triu_indices(n, k=1)
        indptr, indices = ac.neighbor_csr(n, ii.astype(np.int32), jj.astype(np.int32))
        self.assertEqual(np.diff(indptr).tolist(), [n - 1] * n)
        clusters = ac.butina_partition(indptr, indices)
        self.assertEqual([m.tolist() for m in clusters], [list(range(n))])


class MandatoryCoordinateTests(unittest.TestCase):
    """Non-finite Stage-1 coordinates must raise, not become a bogus singleton.

    Every Stage-1 atom is mandatory. Without the entry guard a NaN record is not
    an error but a *plausible* result: the fingerprint prefilter compares with
    ``<=``, which is False against NaN, so the record is pruned from every
    candidate list and comes back as its own cluster.
    """

    @staticmethod
    def _two_tight_pairs_plus(bad_value):
        """Four records that cluster into two pairs, with a fifth carrying `bad_value`.

        Deliberately not a lone bad record: the pre-guard failure was silent, so
        the discriminating input is one where the rest of the bucket clusters
        normally and the bad record merely adds a cluster nobody notices.
        """
        # n_cg=4 plus two vdM slots x 3 backbone atoms = 10 rows per record.
        rng = np.random.default_rng(11)
        a = rng.normal(scale=2.0, size=(10, 3)).astype(np.float32)
        b = rng.normal(scale=2.0, size=(10, 3)).astype(np.float32)
        data = np.stack([a, a + 0.01, b, b + 0.01, a.copy()])
        data[4, 2, 1] = bad_value
        return data

    def test_non_finite_record_raises_instead_of_forming_a_cluster(self):
        for label, bad in (("nan", np.nan), ("+inf", np.inf), ("-inf", -np.inf)):
            with self.subTest(label):
                data = self._two_tight_pairs_plus(bad)
                with self.assertRaisesRegex(ValueError, "non-finite Stage-1"):
                    ac.get_butina_clusters(data, 0.5, 4, build_perm_group(None, 4, ["bb", "bb"]))

    def test_the_bad_record_is_what_raises_not_the_shape(self):
        """Same input with the bad value repaired must cluster, or the test is vacuous."""
        data = self._two_tight_pairs_plus(np.nan)
        data[4, 2, 1] = float(data[0, 2, 1])
        clusters = ac.get_butina_clusters(data, 0.5, 4, build_perm_group(None, 4, ["bb", "bb"]))
        assigned = sorted(i for m in clusters for i in m.tolist())
        self.assertEqual(assigned, list(range(5)))
        # Records 0, 1 and 4 are the same pose; 2 and 3 are the other one.
        self.assertEqual(_partition(clusters), {frozenset({0, 1, 4}), frozenset({2, 3})})

    def test_a_single_non_finite_record_raises_before_the_n_eq_1_shortcut(self):
        # n == 1 returns before any distance is computed, so this row would go
        # straight to disk as an nr vdG if the guard sat any later.
        data = np.zeros((1, 7, 3), dtype=np.float32)
        data[0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "non-finite Stage-1"):
            ac.get_butina_clusters(data, 0.5, 4, build_perm_group(None, 4, ["ALA"]))


class SlotSymmetryTests(unittest.TestCase):
    def test_slot_swapped_duplicate_joins_one_cluster_without_expansion(self):
        """Slot symmetry lives in the distance, so one environment is one record.

        Previously each vdG was replicated once per slot ordering and the copies
        were stitched back together by a transitive whole-cluster union. Here the
        swapped pair simply has RMSD 0 and clusters together.
        """
        rng = np.random.default_rng(8)
        n_cg = 4
        base = rng.normal(scale=3.0, size=(n_cg + 6, 3)).astype(np.float32)
        swapped = base.copy()
        swapped[n_cg:n_cg + 3] = base[n_cg + 3:]
        swapped[n_cg + 3:] = base[n_cg:n_cg + 3]
        data = np.stack((base, swapped))

        clusters = ac.get_butina_clusters(data, 0.5, n_cg, build_perm_group(None, n_cg, ['bb', 'bb']))
        self.assertEqual([m.tolist() for m in clusters], [[0, 1]])

        # With distinguishable slots the same pair is genuinely different.
        apart = ac.get_butina_clusters(data, 0.5, n_cg, build_perm_group(None, n_cg, ['ARG', 'bb']))
        self.assertEqual(len(apart), 2)


class PoseMinimaxPrototypeTests(unittest.TestCase):
    def test_minimax_member_and_exact_radius(self):
        rng = np.random.default_rng(9)
        n_cg, parts = 4, ['ALA']
        data = _clustered_blob(rng, n_cg, 1, n_centers=1, per_center=12, spread=0.4)
        group = build_perm_group(None, n_cg, parts)
        exact = _brute_force_distances(data, n_cg, group)
        members = list(range(len(data)))

        rep, radius = ac.pose_minimax_prototype(data, members, n_cg, group)
        worst = exact[members].max(axis=1)
        self.assertEqual(rep, members[int(np.argmin(worst))])
        self.assertAlmostEqual(radius, worst.min(), places=6)
        # Minimax is what the stored radius is minimised against.
        self.assertLessEqual(radius, exact[members, members[0]].max() + 1e-9)

    def test_single_member_group(self):
        data = np.zeros((3, 7, 3), dtype=np.float32)
        self.assertEqual(ac.pose_minimax_prototype(data, [2], 4), (2, 0.0))


class BatchedRmsdTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(10)
        self.n_cg, parts = 5, ['bb', 'bb']
        self.data = _clustered_blob(rng, self.n_cg, 2, n_centers=3, per_center=8)
        self.group = build_perm_group(_phosphate_permutations(), self.n_cg, parts)
        self.full_perms = full_row_permutations(self.group, self.n_cg)
        self.ii, self.jj = np.triu_indices(len(self.data), k=1)

    def test_batched_matches_pairwise_across_chunk_boundaries(self):
        n_total = self.data.shape[1]
        one_shot = ac._batched_min_rmsd(self.data, self.ii, self.jj, self.full_perms, n_total)
        chunked = ac._batched_min_rmsd(self.data, self.ii, self.jj, self.full_perms,
                                       n_total, batch_rows=len(self.full_perms) * 3)
        np.testing.assert_allclose(one_shot, chunked, atol=1e-6)
        np.testing.assert_allclose(
            one_shot, _brute_force_distances(self.data, self.n_cg, self.group)[self.ii, self.jj],
            atol=1e-5)

    def test_masked_fit_equals_the_full_fit_when_the_argmin_element_survives(self):
        n_total = self.data.shape[1]
        full = ac._batched_min_rmsd(self.data, self.ii, self.jj, self.full_perms, n_total)
        n_perm = len(self.full_perms)
        all_on = np.ones((self.ii.size, n_perm), dtype=bool)
        np.testing.assert_allclose(
            ac._masked_min_rmsd(self.data, self.ii, self.jj, all_on, self.full_perms,
                                n_total, batch_rows=7), full, atol=1e-6)
        # Keep only the argmin element of every pair: still the same value.
        per_elem = np.stack([
            ac._batched_min_rmsd(self.data, self.ii, self.jj, (perm,), n_total)
            for perm in self.full_perms], axis=1)
        only_best = np.zeros_like(all_on)
        only_best[np.arange(self.ii.size), per_elem.argmin(axis=1)] = True
        np.testing.assert_allclose(
            ac._masked_min_rmsd(self.data, self.ii, self.jj, only_best, self.full_perms,
                                n_total), full, atol=1e-6)


class NeighbourOrderTests(unittest.TestCase):
    """Stage-1 neighbour order is load-bearing, so it is pinned here.

    Stage 2 (`get_leader_clusters`) is a leader algorithm: it walks a Stage-1
    cluster's members in the order Butina emitted them, so that order decides the
    Stage-2 partition and, through it, the stored `cluster_size` and
    `cluster_num_parents`. The other tests in this file compare the graph as a
    *set* and partitions as *frozensets*, which is deliberate for the properties
    they assert but leaves order entirely unpinned.

    The order must not depend on where `_flush` happens to land (driven by
    `_GRAPH_BATCH_ROWS`, a performance knob) nor on how the rows were split
    across processes -- otherwise tuning either would silently re-partition the
    library. Those sweeps are the discriminating cases, not the ascending check.
    """

    def _edges(self, data, group, cutoff, batch_rows=None, blocks=None):
        original = ac._GRAPH_BATCH_ROWS
        if batch_rows is not None:
            ac._GRAPH_BATCH_ROWS = batch_rows
        try:
            if blocks is None:
                return ac.stage1_edges(data, cutoff, group[0][0].size, group)
            parts = [ac.stage1_edges(data, cutoff, group[0][0].size, group, row_range=b)
                     for b in blocks]
            return (np.concatenate([p[0] for p in parts]),
                    np.concatenate([p[1] for p in parts]))
        finally:
            ac._GRAPH_BATCH_ROWS = original

    def test_neighbour_lists_are_ascending_and_duplicate_free(self):
        rng = np.random.default_rng(21)
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                group = build_perm_group(perms, n_cg, parts)
                indptr, indices = ac.neighbor_csr(len(data), *self._edges(data, group, 0.6))
                for i in range(len(data)):
                    nbrs = indices[indptr[i]:indptr[i + 1]].tolist()
                    self.assertEqual(nbrs, sorted(set(nbrs)), f'row {i}')
                    self.assertNotIn(i, nbrs)

    def test_edges_do_not_depend_on_the_flush_boundary_or_the_row_split(self):
        """Tiny batches force a flush per row and row blocks split the scan;
        the edge arrays must be byte-identical."""
        rng = np.random.default_rng(22)
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                group = build_perm_group(perms, n_cg, parts)
                base = self._edges(data, group, 0.6)
                self.assertGreater(base[0].size, 0, 'no edges: test would be vacuous')
                n = len(data)
                for batch_rows, blocks in ((1, None), (7, None), (64, None),
                                           (None, [(0, 5), (5, 20), (20, n - 1)]),
                                           (3, [(0, n // 2), (n // 2, n - 1)])):
                    got = self._edges(data, group, 0.6, batch_rows, blocks)
                    np.testing.assert_array_equal(got[0], base[0])
                    np.testing.assert_array_equal(got[1], base[1])

    def test_cluster_member_order_is_stable_across_flush_boundaries(self):
        """What Stage 2 actually consumes: member *lists*, not member sets."""
        rng = np.random.default_rng(23)
        original = ac._GRAPH_BATCH_ROWS
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                runs = []
                for batch_rows in (original, 1, 7):
                    ac._GRAPH_BATCH_ROWS = batch_rows
                    try:
                        clusters = ac.get_butina_clusters(data, 0.6, n_cg, build_perm_group(perms, n_cg, parts))
                    finally:
                        ac._GRAPH_BATCH_ROWS = original
                    runs.append([m.tolist() for m in clusters])
                self.assertTrue(any(len(m) > 1 for m in runs[0]),
                                'all singletons: order test would be vacuous')
                for other in runs[1:]:
                    self.assertEqual(other, runs[0])


if __name__ == '__main__':
    unittest.main()
