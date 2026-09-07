import unittest

import numpy as np

from ligand_vdgs.functions.utils import kabsch


class KabschDegeneracyTests(unittest.TestCase):
    rotation = np.array(
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=np.float32,
    )
    translation = np.array([2.5, -1.25, 4.0], dtype=np.float32)

    def assert_proper_rotation(self, rotation):
        np.testing.assert_allclose(
            rotation.T @ rotation,
            np.eye(3, dtype=np.float32),
            atol=5e-6,
            rtol=0.0,
        )
        self.assertAlmostEqual(float(np.linalg.det(rotation)), 1.0, places=5)

    def test_rank_zero_one_and_two_return_proper_rotations(self):
        point_sets = {
            "rank_zero": np.array(
                [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]],
                dtype=np.float32,
            ),
            "rank_one": np.array(
                [[-2.0, 0.0, 0.0], [-0.5, 0.0, 0.0], [1.0, 0.0, 0.0]],
                dtype=np.float32,
            ),
            "rank_two": np.array(
                [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 1.5, 0.0]],
                dtype=np.float32,
            ),
        }

        for rank_name, X in point_sets.items():
            Y = (X @ self.rotation + self.translation)[None, ...]
            for batched_x in (False, True):
                with self.subTest(rank=rank_name, batched_x=batched_x):
                    X_input = X[None, ...] if batched_x else X
                    rotations, translations, ssds = kabsch(X_input, Y, chunk_size=1)

                    self.assert_proper_rotation(rotations[0])
                    np.testing.assert_allclose(
                        X @ rotations[0] + translations[0],
                        Y[0],
                        atol=1e-5,
                        rtol=0.0,
                    )
                    self.assertLess(float(ssds[0]), 1e-8)

    def test_rank_two_reflection_uses_unconstrained_axis_for_proper_rotation(self):
        X = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 1.5, 0.0]],
            dtype=np.float32,
        )
        reflection = np.diag(np.array([-1.0, 1.0, 1.0], dtype=np.float32))
        Y = (X @ reflection + self.translation)[None, ...]

        rotations, translations, ssds = kabsch(X, Y)

        self.assert_proper_rotation(rotations[0])
        np.testing.assert_allclose(
            X @ rotations[0] + translations[0], Y[0], atol=1e-5, rtol=0.0
        )
        self.assertLess(float(ssds[0]), 1e-8)

    def test_non_finite_inputs_are_rejected_before_svd(self):
        X = np.eye(3, dtype=np.float32)
        Y = X[None, ...].copy()

        for value in (np.nan, np.inf, -np.inf):
            for invalid_input in ("X", "Y"):
                with self.subTest(value=value, invalid_input=invalid_input):
                    bad_X = X.copy()
                    bad_Y = Y.copy()
                    if invalid_input == "X":
                        bad_X[0, 0] = value
                    else:
                        bad_Y[0, 0, 0] = value
                    with self.assertRaisesRegex(ValueError, "Non-finite values"):
                        kabsch(bad_X, bad_Y)


if __name__ == "__main__":
    unittest.main()


class KabschSsdParityTests(unittest.TestCase):
    """kabsch_ssd must agree with kabsch's ssd wherever clustering relies on it.

    It reaches the same number by a different route -- singular values only,
    no rotation and no explicit residual -- so the agreement is algebraic, not
    bit-for-bit. The tolerance below is what stands between the two and a
    different clustering result.
    """

    TOL = 1e-5  # in A of RMSD, against a 0.5 A clustering threshold

    def rmsds(self, X, Y, n):
        from ligand_vdgs.functions.utils import kabsch_ssd
        _, _, ref = kabsch(X, Y)
        return np.sqrt(ref / n), np.sqrt(kabsch_ssd(X, Y) / n)

    def test_agrees_for_a_fixed_X_against_many_Y(self):
        rng = np.random.default_rng(0)
        for n in (4, 7, 8, 11, 14):
            X = rng.normal(size=(n, 3)).astype(np.float32)
            Y = (X[None] + rng.normal(scale=0.3, size=(32, n, 3))).astype(np.float32)
            a, b = self.rmsds(X, Y, n)
            np.testing.assert_allclose(a, b, atol=self.TOL, rtol=0.0)

    def test_agrees_for_batched_X_against_batched_Y(self):
        rng = np.random.default_rng(1)
        n, M = 8, 24
        X = rng.normal(size=(M, n, 3)).astype(np.float32)
        Y = rng.normal(size=(M, n, 3)).astype(np.float32)
        a, b = self.rmsds(X, Y, n)
        np.testing.assert_allclose(a, b, atol=self.TOL, rtol=0.0)

    def test_agrees_when_the_proper_rotation_correction_fires(self):
        # A reflected target: the smallest singular value must enter negatively.
        # Getting this sign wrong is silent and gives a plausible-looking RMSD.
        rng = np.random.default_rng(2)
        n = 6
        X = rng.normal(size=(n, 3)).astype(np.float32)
        Y = (X[None] * np.array([1, 1, -1], np.float32)
             + rng.normal(scale=0.05, size=(8, n, 3))).astype(np.float32)
        a, b = self.rmsds(X, Y, n)
        np.testing.assert_allclose(a, b, atol=self.TOL, rtol=0.0)

    def test_identical_structures_give_exactly_zero_not_a_negative(self):
        # The closed form can land a hair below zero where a residual cannot.
        from ligand_vdgs.functions.utils import kabsch_ssd
        rng = np.random.default_rng(3)
        X = rng.normal(size=(7, 3)).astype(np.float32)
        ssd = kabsch_ssd(X, X[None])
        self.assertGreaterEqual(float(ssd[0]), 0.0)
        self.assertLess(float(ssd[0]), 1e-8)

    def test_rank_deficient_inputs_match(self):
        # Three-point backbone fits are rank 2; collinear ones are rank 1.
        from ligand_vdgs.functions.utils import kabsch_ssd
        for pts in ([[0, 0, 0], [1, 0, 0], [0, 1, 0]],      # coplanar
                    [[0, 0, 0], [1, 0, 0], [2, 0, 0]]):     # collinear
            X = np.asarray(pts, np.float32)
            Y = (X[None] + np.float32(0.1)).astype(np.float32)
            _, _, ref = kabsch(X, Y)
            np.testing.assert_allclose(kabsch_ssd(X, Y), ref, atol=1e-5, rtol=0.0)

    def test_non_finite_inputs_are_rejected(self):
        from ligand_vdgs.functions.utils import kabsch_ssd
        X = np.zeros((4, 3), np.float32)
        Y = np.zeros((1, 4, 3), np.float32)
        Y[0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            kabsch_ssd(X, Y)

    def test_multi_chunk_path_matches_single_chunk(self):
        # chunk_size splits the batch; each chunk must stand alone.
        from ligand_vdgs.functions.utils import kabsch_ssd
        rng = np.random.default_rng(4)
        X = rng.normal(size=(6, 3)).astype(np.float32)
        Y = rng.normal(size=(200, 6, 3)).astype(np.float32)
        np.testing.assert_allclose(
            kabsch_ssd(X, Y, chunk_size=17), kabsch_ssd(X, Y),
            atol=1e-9, rtol=0.0)

    def test_hand_rolled_determinant_matches_numpy(self):
        # Only its sign is used, but a wrong expansion would flip that sign.
        from ligand_vdgs.functions.utils import _det3
        rng = np.random.default_rng(5)
        H = rng.normal(size=(64, 3, 3))
        np.testing.assert_allclose(_det3(H), np.linalg.det(H), atol=1e-9, rtol=1e-9)
        singular = H.copy()
        singular[:, 2] = singular[:, 0] + singular[:, 1]   # rank 2
        np.testing.assert_allclose(_det3(singular), 0.0, atol=1e-9)

    def test_empty_batch_returns_empty(self):
        from ligand_vdgs.functions.utils import kabsch_ssd
        out = kabsch_ssd(np.zeros((4, 3), np.float32), np.zeros((0, 4, 3), np.float32))
        self.assertEqual(out.shape, (0,))
        self.assertEqual(out.dtype, np.float64)


if __name__ == "__main__":
    unittest.main()
