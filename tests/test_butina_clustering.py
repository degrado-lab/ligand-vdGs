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
                clusters, reps, radii = ac.get_butina_clusters(
                    data, cutoff, n_cg, perms, parts)
                for cnum, members in clusters.items():
                    worst = exact[members, reps[cnum]].max()
                    self.assertLessEqual(worst, cutoff + 1e-6)
                    self.assertAlmostEqual(radii[cnum], worst, places=6)

    def test_every_record_belongs_to_exactly_one_cluster(self):
        rng = np.random.default_rng(5)
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                clusters, _reps, _radii = ac.get_butina_clusters(
                    data, 0.6, n_cg, perms, parts)
                assigned = sorted(i for m in clusters.values() for i in m)
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
                adj, _dists = ac._stage1_neighbor_graph(
                    data, cutoff, n_cg, group)
                got = {(min(i, j), max(i, j))
                       for i in range(len(data)) for j in adj[i]}
                want = {(i, j) for i in range(len(data))
                        for j in range(i + 1, len(data)) if exact[i, j] <= cutoff}
                self.assertEqual(got, want)


class ButinaDeterminismTests(unittest.TestCase):
    def test_result_is_invariant_to_input_order(self):
        """Cluster counts are compared across fragments, so order must not matter."""
        rng = np.random.default_rng(7)
        n_cg, parts = 4, ['ARG']
        data = _clustered_blob(rng, n_cg, 1)
        order = rng.permutation(len(data))
        base, _r, _rad = ac.get_butina_clusters(data, 0.6, n_cg, None, parts)
        shuffled, _r2, _rad2 = ac.get_butina_clusters(
            data[order], 0.6, n_cg, None, parts)

        def as_partition(clusters, remap=None):
            out = set()
            for members in clusters.values():
                out.add(frozenset(int(remap[m]) if remap is not None else int(m)
                                  for m in members))
            return out

        self.assertEqual(as_partition(base), as_partition(shuffled, remap=order))

    def test_singleton_and_empty_inputs(self):
        self.assertEqual(ac.get_butina_clusters(
            np.zeros((0, 7, 3), np.float32), 0.5, 4, None, ['ALA']), ({}, {}, {}))
        clusters, reps, radii = ac.get_butina_clusters(
            np.zeros((1, 7, 3), np.float32), 0.5, 4, None, ['ALA'])
        self.assertEqual((clusters, reps, radii), ({1: [0]}, {1: 0}, {1: 0.0}))


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
                    ac.get_butina_clusters(data, 0.5, 4, None, ["bb", "bb"])

    def test_the_bad_record_is_what_raises_not_the_shape(self):
        """Same input with the bad value repaired must cluster, or the test is vacuous."""
        data = self._two_tight_pairs_plus(np.nan)
        data[4, 2, 1] = float(data[0, 2, 1])
        clusters, _reps, _radii = ac.get_butina_clusters(
            data, 0.5, 4, None, ["bb", "bb"])
        assigned = sorted(i for m in clusters.values() for i in m)
        self.assertEqual(assigned, list(range(5)))
        # Records 0, 1 and 4 are the same pose; 2 and 3 are the other one.
        partition = {frozenset(m) for m in clusters.values()}
        self.assertEqual(partition, {frozenset({0, 1, 4}), frozenset({2, 3})})

    def test_a_single_non_finite_record_raises_before_the_n_eq_1_shortcut(self):
        # n == 1 returns before any distance is computed, so this row would go
        # straight to disk as an nr vdG if the guard sat any later.
        data = np.zeros((1, 7, 3), dtype=np.float32)
        data[0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "non-finite Stage-1"):
            ac.get_butina_clusters(data, 0.5, 4, None, ["ALA"])

    def test_empty_input_is_still_accepted(self):
        self.assertEqual(ac.get_butina_clusters(
            np.zeros((0, 7, 3), np.float32), 0.5, 4, None, ["ALA"]), ({}, {}, {}))


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

        clusters, _reps, radii = ac.get_butina_clusters(
            data, 0.5, n_cg, None, ['bb', 'bb'])
        self.assertEqual(len(clusters), 1)
        self.assertEqual(sorted(clusters[1]), [0, 1])
        self.assertAlmostEqual(radii[1], 0.0, places=6)

        # With distinguishable slots the same pair is genuinely different.
        apart, _r, _rad = ac.get_butina_clusters(
            data, 0.5, n_cg, None, ['ARG', 'bb'])
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
    def test_batched_matches_pairwise_across_chunk_boundaries(self):
        rng = np.random.default_rng(10)
        n_cg, parts = 5, ['bb', 'bb']
        data = _clustered_blob(rng, n_cg, 2, n_centers=3, per_center=8)
        group = build_perm_group(_phosphate_permutations(), n_cg, parts)
        full_perms = full_row_permutations(group, n_cg)
        n_total = data.shape[1]
        ii, jj = np.triu_indices(len(data), k=1)

        one_shot = ac._batched_min_rmsd(data, ii, jj, full_perms, n_total)
        chunked = ac._batched_min_rmsd(data, ii, jj, full_perms, n_total,
                                       batch_rows=len(full_perms) * 3)
        np.testing.assert_allclose(one_shot, chunked, atol=1e-6)
        np.testing.assert_allclose(
            one_shot, _brute_force_distances(data, n_cg, group)[ii, jj], atol=1e-5)


class NeighbourOrderTests(unittest.TestCase):
    """Stage-1 neighbour order is load-bearing, so it is pinned here.

    Stage 2 (`get_leader_clusters`) is a leader algorithm: it walks a Stage-1
    cluster's members in the order Butina emitted them, so that order decides the
    Stage-2 partition and, through it, the stored `cluster_size` and
    `cluster_num_parents`. The other tests in this file compare the graph as a
    *set* and partitions as *frozensets*, which is deliberate for the properties
    they assert but leaves order entirely unpinned.

    That matters for two reasons. `adj[x]` being ascending is what lets the
    Python adjacency lists be replaced by a CSR edge array without moving a
    single vdG between subgroups. And the order currently falls out of where
    `_flush` happens to land, which is driven by `_GRAPH_BATCH_ROWS` -- a
    performance knob. If the order were flush-dependent, tuning that constant
    would silently re-partition the library, so the batch-size sweep below is the
    discriminating case, not the ascending check.
    """

    def _graph(self, data, n_cg, parts, perms, cutoff, batch_rows=None):
        group = build_perm_group(perms, n_cg, parts)
        original = ac._GRAPH_BATCH_ROWS
        if batch_rows is not None:
            ac._GRAPH_BATCH_ROWS = batch_rows
        try:
            return ac._stage1_neighbor_graph(data, cutoff, n_cg, group)
        finally:
            ac._GRAPH_BATCH_ROWS = original

    def test_neighbour_lists_are_ascending_and_duplicate_free(self):
        rng = np.random.default_rng(21)
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                adj, _d = self._graph(data, n_cg, parts, perms, 0.6)
                for i, nbrs in enumerate(adj):
                    self.assertEqual(nbrs, sorted(nbrs), f'adj[{i}] not ascending')
                    self.assertEqual(len(nbrs), len(set(nbrs)),
                                     f'adj[{i}] has duplicates')

    def test_neighbour_order_does_not_depend_on_the_flush_boundary(self):
        """Tiny batches force a flush per row; the graph must be byte-identical.

        Failure looks like the same edge *set* with a different order inside some
        `adj[x]` -- which every other test in this file would pass.
        """
        rng = np.random.default_rng(22)
        for label, n_cg, parts, perms in CASES:
            with self.subTest(label):
                data = _clustered_blob(rng, n_cg, len(parts))
                base_adj, base_d = self._graph(data, n_cg, parts, perms, 0.6)
                self.assertTrue(any(base_adj), 'no edges: test would be vacuous')
                for batch_rows in (1, 7, 64):
                    adj, dists = self._graph(
                        data, n_cg, parts, perms, 0.6, batch_rows=batch_rows)
                    self.assertEqual(adj, base_adj, f'batch_rows={batch_rows}')
                    for a, b in zip(dists, base_d):
                        np.testing.assert_allclose(a, b)

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
                        clusters, _r, _rad = ac.get_butina_clusters(
                            data, 0.6, n_cg, perms, parts)
                    finally:
                        ac._GRAPH_BATCH_ROWS = original
                    runs.append({c: list(m) for c, m in clusters.items()})
                self.assertTrue(any(len(m) > 1 for m in runs[0].values()),
                                'all singletons: order test would be vacuous')
                for other in runs[1:]:
                    self.assertEqual(other, runs[0])


if __name__ == '__main__':
    unittest.main()
