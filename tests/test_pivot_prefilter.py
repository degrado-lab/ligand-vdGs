"""Stage-1's pivot (LAESA) prefilter must prune hard and never drop a true edge.

The bound is `d(x,y) >= max_j abs(d(x,p_j) - d(y,p_j))`, admissible because the
symmetry-minimised optimal-superposition RMSD is a quotient metric.
"""
import itertools
import os
import tempfile
import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster as ac
from ligand_vdgs.functions.vdg_fp_utils import INTERNAL_DISTANCE_EPS, build_perm_group
from vacuity import assert_discriminates

THRESHOLD = 0.65
CUTOFF = THRESHOLD + INTERNAL_DISTANCE_EPS

def _phosphate_permutations():
    """4! relabelings of the four oxygens around P -- the 24-permutation CG case."""
    return tuple(tuple(1 if i == 1 else oxy[(0, 2, 3, 4).index(i)] for i in range(5))
                 for oxy in itertools.permutations((0, 2, 3, 4)))

def _two_terminal_oxygens():
    """The real `[O;!R;D2][P;!R;D4](=[O;!R;D1])[O;!R;D1]` group: swap the two D1 oxygens."""
    return ((0, 1, 2, 3), (0, 1, 3, 2))

# (label, n_cg, aa_parts, cg_perms) -- mirrors the symmetry cases in test_butina_clustering.
CASES = (('no symmetry', 4, ['ARG', 'bb'], None),
         ('like slots', 4, ['bb', 'bb'], None),
         ('phosphate D1 pair', 4, ['bb', 'bb'], _two_terminal_oxygens()),
         ('phosphate 24-perm', 5, ['bb', 'bb'], _phosphate_permutations()))

def _clustered_data(n_cg, n_res, seed, groups=6, per_group=25, spread=0.25):
    """Loose clusters, so the prefilter sees both accepts and rejects."""
    rng = np.random.default_rng(seed)
    n = n_cg + 3 * n_res
    return np.concatenate([c + rng.normal(scale=spread, size=(per_group, n, 3))
                           for c in rng.normal(scale=6.0, size=(groups, n, 3))]).astype(np.float32)

def _brute_force(data, n_cg, perm_group):
    """Dense symmetry-minimised RMSD, computed straight from the metric."""
    ii, jj = np.triu_indices(len(data), 1)
    return ii, jj, ac._batched_min_rmsd(
        np.asarray(data, dtype=np.float32), ii, jj,
        np.stack(ac.full_row_permutations(perm_group, n_cg)), data.shape[1])

def _edges(data, threshold, n_cg, perm_group, pivots=None):
    qi, qj = ac.stage1_edges(np.asarray(data, dtype=np.float32), threshold, n_cg,
                             perm_group, pivots=pivots)
    return set(zip(qi.tolist(), qj.tolist()))

class PivotAdmissibilityTests(unittest.TestCase):
    def test_pivot_prefilter_never_drops_a_true_edge(self):
        """Recall against the metric itself, not against another implementation."""
        for label, n_cg, parts, perms in CASES:
            for threshold in (0.65, 1.0):
                with self.subTest(case=label, threshold=threshold):
                    data = _clustered_data(n_cg, len(parts), seed=11)
                    group = build_perm_group(perms, n_cg, parts)
                    ii, jj, dist = _brute_force(data, n_cg, group)
                    truth = {(int(a), int(b)) for a, b, d in zip(ii, jj, dist) if d <= threshold}
                    self.assertTrue(truth, f'{label}: no true edges, comparison would be vacuous')
                    self.assertEqual(_edges(data, threshold, n_cg, group), truth)

    def test_pivot_prefilter_matches_a_neutralised_prefilter(self):
        """A zero embedding prunes nothing, so it is the no-prefilter reference."""
        for label, n_cg, parts, perms in CASES:
            with self.subTest(case=label):
                data = _clustered_data(n_cg, len(parts), seed=5)
                group = build_perm_group(perms, n_cg, parts)
                self.assertEqual(
                    _edges(data, 0.65, n_cg, group),
                    _edges(data, 0.65, n_cg, group,
                           pivots=np.zeros((4, len(data)), dtype=np.float32)))

    def test_the_prefilter_is_actually_doing_work(self):
        """Non-vacuity: the real embedding rejects pairs, the neutralised one cannot."""
        n_cg, parts = 4, ['bb', 'bb']
        data = _clustered_data(n_cg, len(parts), seed=6, groups=8, per_group=15)
        ii, jj = np.triu_indices(len(data), 1)
        def prunes_a_tenth(embedding):
            return bool((np.abs(embedding[:, ii] - embedding[:, jj]).max(axis=0)
                         > CUTOFF).mean() > 0.1)
        assert_discriminates(
            prunes_a_tenth,
            [ac.build_pivot_embedding(data, n_cg, build_perm_group(None, n_cg, parts), CUTOFF)],
            [np.zeros((4, len(data)), dtype=np.float32)],
            'pivot prefilter prunes')

class PivotFalsifierTests(unittest.TestCase):
    def test_an_inflated_embedding_drops_edges(self):
        """Falsifier for the edge-identity check: over-prune and the check must refuse."""
        n_cg, parts = 4, ['bb', 'bb']
        data = _clustered_data(n_cg, len(parts), seed=11)
        group = build_perm_group(None, n_cg, parts)
        truth = _edges(data, THRESHOLD, n_cg, group)
        self.assertTrue(truth, 'no true edges, the falsifier would be vacuous')
        self.assertTrue(_edges(data, THRESHOLD, n_cg, group,
                               pivots=ac.build_pivot_embedding(data, n_cg, group, CUTOFF) * 5.0)
                        < truth,
                        'inflating the bound did not drop any edge, so an edge-set '
                        'comparison could not detect an inadmissible prefilter')

    def test_closure_check_refuses_a_truncated_group(self):
        """A truncated automorphism set is not a group and breaks the triangle inequality."""
        full = build_perm_group(_phosphate_permutations(), 5, ['bb', 'bb'])
        ac.assert_perm_group_closed(full)
        with self.assertRaises(ValueError):
            ac.assert_perm_group_closed(tuple(full[:3]))

class PivotMemmapFillTests(unittest.TestCase):
    """Disjoint column spans are filled by separate workers; a partial fill would
    leave zeros, which keep the bound admissible but silently disable pruning."""

    def _ids_and_reference(self, data, n_cg, group):
        ids = ac.select_pivots(data, n_cg, group, ac.pivot_count(len(data)), CUTOFF)
        return ids, ac.pivot_distances(data, ids, n_cg, group)

    def test_spans_assemble_to_the_single_process_embedding(self):
        n_cg, parts = 4, ['bb', 'bb']
        data = _clustered_data(n_cg, len(parts), seed=3, groups=5, per_group=12)
        group = build_perm_group(None, n_cg, parts)
        ids, reference = self._ids_and_reference(data, n_cg, group)
        spans = [(0, 17), (17, 41), (41, len(data))]
        self.assertGreater(len(spans), 1, 'a single span would not test concurrent fill')
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'pivots.npy')
            np.lib.format.open_memmap(path, mode='w+', dtype=np.float32,
                                      shape=(len(ids), len(data)))
            for start, stop in spans:
                out = np.lib.format.open_memmap(path, mode='r+')
                ac.pivot_distances(data, ids, n_cg, group, row_range=(start, stop),
                                   out=out[:, start:stop])
                out.flush()
            np.testing.assert_array_equal(np.load(path), reference)

    def test_a_skipped_span_is_detectable(self):
        """Falsifier: the zeros a missed span leaves must be visible to the check above."""
        n_cg, parts = 4, ['bb', 'bb']
        data = _clustered_data(n_cg, len(parts), seed=3, groups=5, per_group=12)
        group = build_perm_group(None, n_cg, parts)
        ids, reference = self._ids_and_reference(data, n_cg, group)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'pivots.npy')
            np.lib.format.open_memmap(path, mode='w+', dtype=np.float32,
                                      shape=(len(ids), len(data)))
            out = np.lib.format.open_memmap(path, mode='r+')
            ac.pivot_distances(data, ids, n_cg, group, row_range=(0, 17), out=out[:, 0:17])
            out.flush()
            self.assertFalse(np.array_equal(np.load(path), reference),
                             'skipping a span left an embedding indistinguishable from '
                             'a complete one, so the assembly check cannot fire')

    def test_a_mismatched_out_shape_raises(self):
        """The guard that turns a partial fill into a failure instead of slow zeros."""
        n_cg, parts = 4, ['bb', 'bb']
        data = _clustered_data(n_cg, len(parts), seed=3, groups=5, per_group=12)
        group = build_perm_group(None, n_cg, parts)
        ids = ac.select_pivots(data, n_cg, group, ac.pivot_count(len(data)), CUTOFF)
        with self.assertRaises(ValueError):
            ac.pivot_distances(data, ids, n_cg, group, row_range=(0, 17),
                               out=np.zeros((len(ids), 16), dtype=np.float32))

if __name__ == '__main__':
    unittest.main()
