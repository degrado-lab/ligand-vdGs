"""The fingerprint prefilter's tolerances must be admissible.

`fp_tolerances` claims that anything the prefilter rejects is provably outside
the RMSD threshold. That is a correctness claim, not a tuning choice: if it is
wrong, clustering and hit finding silently drop true matches without ever
computing the fit that would have found them. These tests pin the closed form,
show it is not exceeded, show it is tight (so it is not loose enough to be
useless), and check that tightening it does not move any cluster.
"""

import math
import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster as ac
from ligand_vdgs.functions.utils import kabsch_ssd
from ligand_vdgs.functions.vdg_fp_utils import (FP_SAFETY_FACTOR, build_perm_group, fp_tolerances,
                                                precompute_bucket_fingerprints)


def fingerprints(X, n_cg):
    """The same quantities precompute_bucket_fingerprints emits, for one vdG."""
    com = X[:n_cg].mean(axis=0)
    ca1 = X[n_cg + 1]
    out = [np.linalg.norm(ca1 - com)]
    if len(X) - n_cg >= 6:
        ca2 = X[n_cg + 4]
        out += [np.linalg.norm(ca2 - com), np.linalg.norm(ca1 - ca2)]
    return np.array(out)


def rmsd(X, Y):
    return float(np.sqrt(kabsch_ssd(X[None, ...], Y[None, ...])[0] / len(X)))


class ToleranceFormulaTests(unittest.TestCase):
    def test_closed_form(self):
        # fp0/fp1: T*sqrt(n*(1+1/n_cg));  fp2: T*sqrt(2n), each times the factor
        k = 0.5 * FP_SAFETY_FACTOR
        t = fp_tolerances(0.5, 10, 4, 2)
        self.assertAlmostEqual(t[0], k * math.sqrt(10 * 1.25), places=12)
        self.assertAlmostEqual(t[1], t[0], places=12)
        self.assertAlmostEqual(t[2], k * math.sqrt(20), places=12)
        self.assertEqual(len(fp_tolerances(0.5, 7, 4, 1)), 1)

    def test_safety_factor_may_only_loosen(self):
        # Below 1.0 the returned tolerance drops under the proven bound and the
        # prefilter starts discarding true matches. That is the one edit to this
        # constant that breaks correctness rather than just costing speed.
        self.assertGreaterEqual(FP_SAFETY_FACTOR, 1.0)

    def test_tolerance_always_exceeds_the_threshold(self):
        # _dist_idx_idx returns `total + tol` as a reject sentinel, which only
        # rejects if tol > threshold. Holds for every shape the library builds.
        for n_cg in range(1, 13):
            for n_res in (1, 2):
                n = n_cg + 3 * n_res
                for T in (0.5, 0.75, 1.0, 1.5):
                    self.assertGreater(min(fp_tolerances(T, n, n_cg, n_res)), T)

    def test_unsupported_subset_size_raises(self):
        with self.assertRaises(ValueError):
            fp_tolerances(0.5, 13, 4, 3)


class AdmissibilityTests(unittest.TestCase):
    """No pair within the RMSD threshold may exceed the tolerance."""

    def test_random_perturbations_never_exceed_the_bound(self):
        rng = np.random.default_rng(11)
        for n_cg, n_res, T in [(4, 1, 0.5), (4, 2, 0.65), (5, 2, 1.0), (12, 2, 1.0)]:
            n = n_cg + 3 * n_res
            tol = np.array(fp_tolerances(T, n, n_cg, n_res))
            for _ in range(400):
                X = rng.normal(scale=3.0, size=(n, 3))
                D = rng.normal(size=(n, 3))
                D *= T / np.sqrt((D * D).sum() / n)   # place the pair at RMSD = T
                Y = X + D
                if rmsd(X, Y) > T + 1e-6:
                    continue
                d = np.abs(fingerprints(X, n_cg) - fingerprints(Y, n_cg))
                # Checked against the unpadded bound: the guarantee is a property
                # of the derivation, not of FP_SAFETY_FACTOR.
                bound = tol / FP_SAFETY_FACTOR
                self.assertTrue((d <= bound + 1e-5).all(),
                                f'{d} exceeded {bound} at n_cg={n_cg} n_res={n_res}')

    def test_the_bound_is_tight_not_merely_valid(self):
        # A loose bound is admissible but prunes nothing. The worst case for fp0
        # displaces the CA outward while every CG atom moves inward; for fp2 it
        # pushes the two CAs apart along their own axis. Both reach the bound.
        n_cg, n_res, T = 5, 2, 1.0
        n = n_cg + 3 * n_res
        tol = fp_tolerances(T, n, n_cg, n_res)
        u = np.array([1.0, 0.0, 0.0])

        X = np.zeros((n, 3))
        X[:n_cg] = -8.0 * u
        X[n_cg + 1] = 8.0 * u
        X[n_cg + 4] = 4.0 * u
        eps = T * math.sqrt(n / (n_cg * (n_cg + 1)))
        D = np.zeros((n, 3))
        D[:n_cg] = -eps * u
        D[n_cg + 1] = n_cg * eps * u
        got = abs(fingerprints(X, n_cg)[0] - fingerprints(X + D, n_cg)[0])
        self.assertLessEqual(rmsd(X, X + D), T + 1e-6)
        self.assertGreater(got / (tol[0] / FP_SAFETY_FACTOR), 0.97)

        X2 = np.zeros((n, 3))
        X2[:n_cg] = np.linspace(-1, 1, n_cg)[:, None] * u
        X2[n_cg + 1] = -6.0 * u
        X2[n_cg + 4] = 6.0 * u
        d = T * math.sqrt(n / 2)
        D2 = np.zeros((n, 3))
        D2[n_cg + 1] = -d * u
        D2[n_cg + 4] = d * u
        got2 = abs(fingerprints(X2, n_cg)[2] - fingerprints(X2 + D2, n_cg)[2])
        self.assertLessEqual(rmsd(X2, X2 + D2), T + 1e-6)
        self.assertGreater(got2 / (tol[2] / FP_SAFETY_FACTOR), 0.97)


class ClusteringUnchangedTests(unittest.TestCase):
    """Tightening an admissible prefilter must not move a single cluster."""

    def _cluster(self, data, n_cg, threshold, disabled):
        # Patch the name in align_and_cluster, not in vdg_fp_utils: the module
        # imports it directly, so patching the source has no effect. If that
        # import style ever changes, this patch stops firing and the comparison
        # below goes vacuous -- test_the_prefilter_is_actually_doing_work is the
        # guard for that.
        original = ac.fp_tolerances
        if disabled:
            ac.fp_tolerances = lambda t, n, ncg, nres: (1e9,) * (1 if nres == 1 else 3)
        try:
            clusters = ac.get_butina_clusters(
                data, threshold, n_cg, build_perm_group(None, n_cg, ['ARG', 'bb']))
            return [members.tolist() for members in clusters]
        finally:
            ac.fp_tolerances = original

    def test_derived_tolerance_matches_no_prefilter(self):
        rng = np.random.default_rng(5)
        n_cg, n_res = 4, 2
        n = n_cg + 3 * n_res
        # Several loose groups so the prefilter has both accepts and rejects.
        centers = rng.normal(scale=6.0, size=(6, n, 3))
        data = np.concatenate(
            [c + rng.normal(scale=0.25, size=(25, n, 3)) for c in centers]
        ).astype(np.float32)
        for threshold in (0.65, 1.0):
            with self.subTest(threshold=threshold):
                self.assertEqual(self._cluster(data, n_cg, threshold, False),
                                 self._cluster(data, n_cg, threshold, True))

    def test_the_prefilter_is_actually_doing_work(self):
        # Guards the test above from passing because nothing is ever rejected.
        rng = np.random.default_rng(6)
        n_cg, n_res, T = 4, 2, 0.65
        n = n_cg + 3 * n_res
        data = rng.normal(scale=6.0, size=(120, n, 3)).astype(np.float32)
        fps = precompute_bucket_fingerprints(data, n_cg)
        tol = fp_tolerances(T, n, n_cg, n_res)
        rejected = (np.abs(fps['fp2'][:, None] - fps['fp2'][None, :]) > tol[2])
        self.assertGreater(rejected.mean(), 0.1)


if __name__ == '__main__':
    unittest.main()
