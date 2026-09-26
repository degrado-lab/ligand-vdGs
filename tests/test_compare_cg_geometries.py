"""compare_cg_geometries: mmd2 (unbiased weighted MMD) and the ideal-frame CG alignment.

Falsifiers: a diagonal-included mmd2 gives 17.0 not 14.0 on a hand-computed 2-vs-2 case;
an unweighted mmd2 disagrees with a reference double-loop reimplementation; a frame fit
that leaks the raw world frame fails to recover a known CG offset under a random
per-observation transform, or changes under a shared one; a self-split using randomized
str hash() (not a fixed hash) differs across interpreter processes; a fit that lets vdM
slot 1 leak in changes when slot 1 is perturbed with slot 0 and the CG held fixed.
"""
import os
import sys
import subprocess
import unittest

import numpy as np

from ligand_vdgs.identify_bioisosteres.compare_cg_geometries import (
    IDEAL_NCAC, _ideal_frame_points, compare_bucket_geometry, mmd2)

class TestMMD2Formula(unittest.TestCase):
    def test_excludes_diagonal_self_pairs(self):
        w = np.ones(2)
        # Hand-derived: offdiag_mean(A)=offdiag_mean(B)=3, cross_mean=10 -> mmd2=2*10-3-3=14.
        # A diagonal-included estimator would give offdiag_mean=1.5 -> mmd2=17 instead.
        self.assertAlmostEqual(
            mmd2(np.array([[0., 0., 0.], [3., 0., 0.]]), w,
                np.array([[10., 0., 0.], [13., 0., 0.]]), w), 14.0, places=10)

    def test_respects_nonuniform_weights_vs_reference_loop(self):
        rng = np.random.default_rng(0)
        pts_A, w_A = rng.normal(size=(6, 3)), rng.uniform(0.5, 5, 6)
        pts_B, w_B = rng.normal(size=(7, 3)) + 4, rng.uniform(0.5, 5, 7)

        def ref_mmd2(pA, wA, pB, wB):
            def offdiag(p, w):
                num = den = 0.0
                for i in range(len(w)):
                    for j in range(len(w)):
                        if i != j:
                            num += w[i] * w[j] * np.linalg.norm(p[i] - p[j])
                            den += w[i] * w[j]
                return num / den
            def cross(pA, wA, pB, wB):
                num = den = 0.0
                for i in range(len(wA)):
                    for j in range(len(wB)):
                        num += wA[i] * wB[j] * np.linalg.norm(pA[i] - pB[j])
                        den += wA[i] * wB[j]
                return num / den
            return 2 * cross(pA, wA, pB, wB) - offdiag(pA, wA) - offdiag(pB, wB)

        self.assertAlmostEqual(mmd2(pts_A, w_A, pts_B, w_B),
                               ref_mmd2(pts_A, w_A, pts_B, w_B), places=8)

    def test_nan_below_two_points(self):
        pts, w = np.zeros((1, 3)), np.ones(1)
        self.assertTrue(np.isnan(mmd2(pts, w, pts, w)))

    def test_same_population_split_gives_small_two_sided_value(self):
        pts, w = np.random.default_rng(4).normal(scale=5, size=(120, 3)), np.ones(120)
        half = np.arange(120) % 2 == 0
        # Unbiased estimator of 0 for a same-population split: small relative to the ~5-unit
        # point scale, and NOT required to be >= 0 (Gretton et al. 2012).
        self.assertLess(abs(mmd2(pts[half], w[half], pts[~half], w[~half])), 1.0)

class TestIdealFrame(unittest.TestCase):
    def _random_rigid(self, rng):
        R, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(R) < 0:
            R[:, 0] *= -1
        return R, rng.normal(size=3) * 20

    def test_recovers_known_cg_offset_under_random_per_observation_transform(self):
        rng = np.random.default_rng(1)
        n = 25
        cg_offset_ideal = np.array([2.0, -1.0, 3.0])
        bb_slot0, cg_raw = np.empty((n, 3, 3)), np.empty((n, 1, 3))
        for i in range(n):
            R, t = self._random_rigid(rng)
            # utils.kabsch(X, Y) returns R, t for X@R+t=Y, so the inverse maps ideal-frame
            # points back to this observation's raw (pre-fit) space.
            bb_slot0[i] = (IDEAL_NCAC - t) @ R.T
            cg_raw[i, 0] = (cg_offset_ideal - t) @ R.T
        np.testing.assert_allclose(_ideal_frame_points(bb_slot0, cg_raw),
                                   np.broadcast_to(cg_offset_ideal, (n, 3)), atol=1e-4)

    def test_rigid_invariant_to_a_shared_world_transform(self):
        rng = np.random.default_rng(2)
        n = 10
        bb_slot0 = IDEAL_NCAC + rng.normal(size=(n, 3, 3)) * 0.3
        cg_raw = rng.normal(size=(n, 1, 3)) * 2
        R, t = self._random_rigid(rng)
        np.testing.assert_allclose(_ideal_frame_points(bb_slot0, cg_raw),
                                   _ideal_frame_points(bb_slot0 @ R + t, cg_raw @ R + t),
                                   atol=1e-4)

class TestSubsetSize2Anchoring(unittest.TestCase):
    def test_slot1_does_not_affect_mmd2(self):
        def make_lib(n, seed):
            return (np.repeat(np.broadcast_to(IDEAL_NCAC, (n, 3, 3))[:, None], 2, axis=1).copy(),
                    np.random.default_rng(seed).normal(size=(n, 1, 3)), np.ones(n),
                    np.array([f'p{i}' for i in range(n)]))

        bb_A, cg_A, w_A, parent_A = make_lib(8, 10)
        bb_B, cg_B, w_B, parent_B = make_lib(8, 11)
        bb_A_perturbed = bb_A.copy()
        bb_A_perturbed[:, 1] = np.random.default_rng(3).normal(size=(8, 3, 3)) * 50

        self.assertAlmostEqual(
            compare_bucket_geometry(cg_A, bb_A, w_A, parent_A, cg_B, bb_B, w_B, parent_B)['mmd2'],
            compare_bucket_geometry(cg_A, bb_A_perturbed, w_A, parent_A,
                                    cg_B, bb_B, w_B, parent_B)['mmd2'], places=8)

class TestSelfSplitDeterminism(unittest.TestCase):
    def test_split_is_stable_across_hash_seeds(self):
        script = (
            "import numpy as np\n"
            "from ligand_vdgs.identify_bioisosteres.compare_cg_geometries import _self_split_mmd2\n"
            "rng = np.random.default_rng(5)\n"
            "n = 40\n"
            "pts = rng.normal(size=(n, 3))\n"
            "w = rng.uniform(0.5, 3, n)\n"
            "parent = np.array([f'parent_{i % 13}' for i in range(n)])\n"
            "print(repr(_self_split_mmd2(pts, w, parent)))\n")
        results = set()
        for seed in ('0', '1', '2'):
            out = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True,
                                 env={**os.environ, 'PYTHONHASHSEED': seed})
            self.assertEqual(out.returncode, 0, out.stderr)
            results.add(out.stdout.strip())
        self.assertEqual(len(results), 1, f"split differed across hash seeds: {results}")

if __name__ == '__main__':
    unittest.main()
