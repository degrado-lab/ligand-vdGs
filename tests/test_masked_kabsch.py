"""The batched Stage-2 flank fit must return what the scalar one did.

`kabsch_ssd` rejects non-finite input and centres over a fixed point count, so
Stage 2 used to drop rows per pair and call it one pair at a time. The batch
version carries the mask as weights instead; if it disagrees anywhere, the
Stage-2 partition silently changes and cannot be recovered (`mem_` rows hold no
coordinates).
"""
import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster as ac
from ligand_vdgs.functions.utils import kabsch_ssd


def _reference_rmsd(X, Y):
    """The pre-batch implementation: drop rows non-finite in either, then fit."""
    X = np.asarray(X, dtype=np.float32)
    Y = np.asarray(Y, dtype=np.float32)
    valid = np.isfinite(X).all(axis=1) & np.isfinite(Y).all(axis=1)
    if not valid.any():
        return float("inf")
    X, Y = X[valid], Y[valid]
    return float(np.sqrt(kabsch_ssd(X[None, ...], Y[None, ...], chunk_size=1)[0]
                         / len(X)))


class MaskedKabschTests(unittest.TestCase):
    TOL = 1e-9   # A, against a 0.5 A clustering threshold

    def setUp(self):
        self.rng = np.random.default_rng(7)
        self.X = self.rng.normal(size=(10, 3)).astype(np.float32)
        self.Y = self.rng.normal(size=(10, 3)).astype(np.float32)

    def _assert_matches(self, label, X, Y):
        want, got = _reference_rmsd(X, Y), ac._rmsd_pair(X, Y)
        if np.isinf(want):
            self.assertTrue(np.isinf(got), label)
        else:
            self.assertLessEqual(abs(want - got), self.TOL, f"{label}: {want} vs {got}")

    def test_missing_rows_are_dropped_exactly_as_before(self):
        nan_y = self.Y.copy()
        nan_y[3] = np.nan
        nan_x = self.X.copy()
        nan_x[7] = np.nan
        same_row = self.X.copy()
        same_row[3] = np.nan
        two_valid = np.full((10, 3), np.nan, dtype=np.float32)
        two_valid[:2] = self.X[:2]
        one_valid = np.full((10, 3), np.nan, dtype=np.float32)
        one_valid[:1] = self.X[:1]
        for label, X, Y in (
            ("all finite", self.X, self.Y),
            ("one missing row in Y", self.X, nan_y),
            ("disjoint missing rows", nan_x, nan_y),
            ("same row missing in both", same_row, nan_y),
            ("two rows left -- rank-deficient fit", two_valid, self.Y),
            ("one row left -- fit is exact, RMSD 0", one_valid, self.Y),
        ):
            with self.subTest(label):
                self._assert_matches(label, X, Y)

    def test_no_shared_rows_is_infinite_not_zero(self):
        """Zero would read as a perfect match and merge two unrelated vdGs."""
        blank = np.full((10, 3), np.nan, dtype=np.float32)
        self.assertTrue(np.isinf(ac._rmsd_pair(blank, self.Y)))
        self.assertTrue(np.isinf(ac._rmsd_pair(blank, blank)))

    def test_degenerate_fits(self):
        """A rank-deficient cross-covariance must not change the sign handling."""
        t = np.linspace(0.0, 1.0, 10, dtype=np.float32)
        collinear = np.stack([t, 2.0 * t, 3.0 * t], axis=1)
        planar = self.X.copy()
        planar[:, 2] = 0.0
        mirrored = self.X * np.array([-1.0, 1.0, 1.0], dtype=np.float32)
        for label, X, Y in (
            ("collinear vs general", collinear, self.Y),
            ("collinear vs scaled collinear", collinear,
             (collinear * 1.7 + 0.3).astype(np.float32)),
            ("planar vs planar", planar, planar[::-1].copy()),
            # Reflection is not a rotation: the fit must not report 0 here.
            ("mirror image", self.X, mirrored),
            ("identical", self.X, self.X.copy()),
        ):
            with self.subTest(label):
                self._assert_matches(label, X, Y)
        self.assertGreater(ac._rmsd_pair(self.X, mirrored), 0.1)
        self.assertEqual(ac._rmsd_pair(self.X, self.X.copy()), 0.0)

    def test_batch_matches_pair_by_pair_over_random_masks(self):
        rng = np.random.default_rng(11)
        n_pairs, n_rows = 2000, 10
        A = rng.normal(size=(n_pairs, n_rows, 3)).astype(np.float32)
        B = rng.normal(size=(n_pairs, n_rows, 3)).astype(np.float32)
        A[rng.random((n_pairs, n_rows)) < 0.2] = np.nan
        B[rng.random((n_pairs, n_rows)) < 0.2] = np.nan
        batched = ac._rmsd_rows(A, B)
        scalar = np.array([_reference_rmsd(A[i], B[i]) for i in range(n_pairs)])
        np.testing.assert_array_equal(np.isinf(batched), np.isinf(scalar))
        finite = np.isfinite(scalar)
        self.assertGreater(finite.sum(), n_pairs // 2)   # the masks did not kill everything
        self.assertLessEqual(np.abs(batched[finite] - scalar[finite]).max(), self.TOL)

    def test_chunking_does_not_change_the_answer(self):
        rng = np.random.default_rng(12)
        A = rng.normal(size=(97, 10, 3)).astype(np.float32)
        B = rng.normal(size=(97, 10, 3)).astype(np.float32)
        A[rng.random((97, 10)) < 0.2] = np.nan
        whole, _ = ac.masked_kabsch_ssd(A, B)
        split, _ = ac.masked_kabsch_ssd(A, B, chunk_size=7)
        np.testing.assert_array_equal(whole, split)


class BatchedRowSumTests(unittest.TestCase):
    """`batched_pair_row_sums` is what picks every Stage-2 cluster centre."""

    @staticmethod
    def _scalar_row_sums(members, n_orders, flankbb_arr, seq_terms):
        """The double Python loop the batch replaced, written out."""
        k = len(members)
        sums = np.zeros(k, dtype=np.float64)
        for i in range(k):
            for j in range(i + 1, k):
                a, b = sorted((members[i], members[j]))
                best = min(
                    float(seq_terms(np.array([a]), np.array([b]), p)[0])
                    + (ac._rmsd_pair(flankbb_arr[p][a], flankbb_arr[0][b])
                       if flankbb_arr is not None else 0.0)
                    for p in range(n_orders))
                sums[i] += best
                sums[j] += best
        return sums

    def test_matches_the_scalar_double_loop(self):
        rng = np.random.default_rng(13)
        n, rows = 40, 10
        data = rng.normal(size=(n, rows, 3)).astype(np.float32)
        data[rng.random((n, rows)) < 0.15] = np.nan
        flankbb_arr = [data, data[:, ::-1].copy()]           # two slot orderings
        seq = rng.random((n, n)) * 0.4
        seq = np.minimum(seq, seq.T)

        def seq_terms(a, b, p):
            return seq[a, b] + 0.01 * p

        members = list(range(3, 33))
        want = self._scalar_row_sums(members, 2, flankbb_arr, seq_terms)
        got = ac.batched_pair_row_sums(members, 2, flankbb_arr, seq_terms)
        np.testing.assert_allclose(got, want, atol=1e-9)
        # The medoid is an argmin over these, so agreeing on the winner is the
        # part that actually reaches the library.
        self.assertEqual(int(np.argmin(got)), int(np.argmin(want)))

    def test_slot_ordering_minimises_the_summed_distance_not_either_term(self):
        """Ordering 0 wins on sequence, ordering 1 on geometry; the sum decides.

        Minimising the two terms separately would give the best sequence term
        (0.1) plus the best geometric one (0.0); the true distance is ordering
        1's 0.4 + 0.0, and `_dist_idx_idx` minimises the sum.
        """
        rng = np.random.default_rng(15)
        rows = 6
        # Two genuinely non-superimposable flanks, so the identity ordering has
        # a large geometric term; ordering 1 maps record 0's flank onto B, so
        # its geometric term is 0 but its sequence term is worse.
        A = rng.normal(size=(rows, 3)).astype(np.float32)
        B = rng.normal(size=(rows, 3)).astype(np.float32)
        flankbb_arr = [np.stack([A, B]), np.stack([B, B])]
        seq_by_p = {0: 0.1, 1: 0.4}

        def seq_terms(a, b, p):
            return np.full(a.size, seq_by_p[p], dtype=np.float64)

        got = ac.batched_pair_row_sums([0, 1], 2, flankbb_arr, seq_terms)
        d0 = seq_by_p[0] + ac._rmsd_pair(flankbb_arr[0][0], flankbb_arr[0][1])
        d1 = seq_by_p[1] + ac._rmsd_pair(flankbb_arr[1][0], flankbb_arr[0][1])
        self.assertLess(d1, d0)                      # the case is the one intended
        np.testing.assert_allclose(got, [min(d0, d1)] * 2, atol=1e-9)
        # Minimising each term on its own would give 0.1 + 0.0; the real answer
        # is 0.4, so that shortcut is excluded by this assertion.
        self.assertGreater(got[0], seq_by_p[0] + 1e-6)

    def test_missing_flanks_do_not_silently_become_zero_distance(self):
        rows = 8
        data = np.full((4, rows, 3), np.nan, dtype=np.float32)
        rng = np.random.default_rng(14)
        data[0] = rng.normal(size=(rows, 3))
        data[1] = rng.normal(size=(rows, 3))
        # Records 2 and 3 share no readable flank row with anything.
        sums = ac.batched_pair_row_sums([0, 1, 2, 3], 1, [data],
                                        lambda a, b, p: np.zeros(a.size))
        self.assertTrue(np.isfinite(sums[0]) or np.isinf(sums[0]))
        self.assertTrue(np.isinf(sums[2]) and np.isinf(sums[3]))


class ScalarBatchedEquivalenceTests(unittest.TestCase):
    """The batched Stage-2 path must partition exactly as the scalar one did.

    `get_leader_clusters(..., _batched=False)` is the scalar reference: same
    control flow, same early exits, distances computed one pair at a time. The
    batched path only moves where the arithmetic happens, so every cluster id,
    member list and stored prototype has to come out identical. A disagreement
    here is unrecoverable in a real build -- `mem_` rows carry no coordinates.
    """

    @staticmethod
    def _bucket(rng, n, n_res, n_flank, missing, n_cg=4):
        """A synthetic bucket with the features that break naive batching:
        missing flanks, duplicate records, and tight clusters."""
        n_centres = -(-n // 8)   # ceil, so the blobs cover every record
        centres = rng.normal(scale=3.0, size=(n_centres, n_flank, 3))
        flank_ca = np.concatenate(
            [c + rng.normal(scale=0.25, size=(8, n_flank, 3)) for c in centres]
        )[:n].astype(np.float32)
        # Exact duplicates: ties, which is where a changed tie-break shows up.
        flank_ca[1] = flank_ca[0]
        flank_ca[rng.random((n, n_flank)) < missing] = np.nan
        tokens = np.array(["ALA", "GLY", "SER", "-", "!"])
        seq = tokens[rng.integers(0, 5, size=(n, n_flank))].astype("U4")
        seq[:, n_flank // 2 :: n_flank] = "vdm"
        n_total = n_cg + 3 * n_res
        pose = np.concatenate(
            [c + rng.normal(scale=0.2, size=(8, n_total, 3))
             for c in rng.normal(scale=3.0, size=(n_centres, n_total, 3))]
        )[:n].astype(np.float32)
        return seq.tolist(), list(flank_ca), pose

    def _run(self, seq, flank_ca, threshold, orders, batched):
        ac.clear_caches()
        out = ac.get_leader_clusters(
            zip([seq, flank_ca], ["flankseq", "flankbb"]),
            threshold, missing_seq_similarity=0.4,
            final_exact_medoid_pass=True, final_reassign_once=True,
            slot_orders=orders, _batched=batched)
        ac.clear_caches()
        return out

    def test_seeded_buckets_partition_identically(self):
        for seed in range(6):
            rng = np.random.default_rng(seed)
            n_res = 1 + seed % 2
            n_flank = n_res * 5
            orders = ((0, 1), (1, 0)) if n_res == 2 else ((0,),)
            n = 60 + 17 * seed
            seq, flank_ca, pose = self._bucket(
                rng, n, n_res, n_flank, missing=0.05 * (seed % 4))
            threshold = 0.5 + 0.1 * (seed % 3)
            with self.subTest(seed=seed, n=n, n_res=n_res):
                want = self._run(seq, flank_ca, threshold, orders, False)
                got = self._run(seq, flank_ca, threshold, orders, True)
                # Cluster ids and member lists, not just the partition: the ids
                # become first/second_stage_cluster_id in the library.
                self.assertEqual(sorted(want), sorted(got))
                for cid in want:
                    self.assertEqual(list(want[cid]), list(got[cid]), f"cluster {cid}")
                # The subgroups actually cover the bucket.
                self.assertEqual(sorted(i for mem in got.values() for i in mem),
                                 list(range(n)))
                # ... and the stored row for each, which is what the library keeps.
                perm_group = ((np.arange(4, dtype=np.intp),
                               np.arange(pose.shape[1] - 4, dtype=np.intp)),)
                for cid in want:
                    a = ac.pose_minimax_prototype(pose, want[cid], 4, perm_group)
                    b = ac.pose_minimax_prototype(pose, got[cid], 4, perm_group)
                    self.assertEqual(a[0], b[0], f"prototype row, cluster {cid}")
                    self.assertAlmostEqual(a[1], b[1], places=9)

    def test_the_reference_path_is_actually_exercised(self):
        """Guard against the switch being ignored and both arms running batched."""
        rng = np.random.default_rng(99)
        seq, flank_ca, _ = self._bucket(rng, 40, 1, 5, missing=0.1)
        # `_rmsd_pair` is itself a one-pair `_rmsd_rows`, so the signal is the
        # *batch size*, not the call count.
        batched_calls = []
        real = ac._rmsd_rows

        def counting(X, Y):
            if len(X) > 1:
                batched_calls.append(len(X))
            return real(X, Y)

        try:
            ac._rmsd_rows = counting
            self._run(seq, flank_ca, 0.6, ((0,),), False)
            self.assertEqual(batched_calls, [], "scalar reference used a batched fit")
            self._run(seq, flank_ca, 0.6, ((0,),), True)
            self.assertGreater(len(batched_calls), 0, "batched path never batched")
        finally:
            ac._rmsd_rows = real


if __name__ == "__main__":
    unittest.main()
