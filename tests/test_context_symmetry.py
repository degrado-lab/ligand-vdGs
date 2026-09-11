"""Permuted-RMSD and the CG symmetry wiring in leader clustering.

Clustering compares every pair over the CG's full automorphism group.
These tests pin the full group reaching the RMSD.
"""

import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster
from ligand_vdgs.functions import utils
from ligand_vdgs.functions.vdg_fp_utils import build_perm_group

PHOSPHATE = "O=P(O)(O)O"  # slots: [O(=), P, O, O, O]


def _phosphate_perms():
    return utils.identify_mol_automorphisms(utils.mol_from_fragment(PHOSPHATE))


class PermutedRmsdTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.n_cg = 5
        self.X = rng.normal(size=(11, 3)).astype(np.float32)
        self.Y = rng.normal(size=(11, 3)).astype(np.float32)
        self.perms = _phosphate_perms()

    def test_phosphate_full_group_is_24(self):
        self.assertEqual(len(self.perms), 24)

    def test_rmsd_is_symmetric(self):
        forward = align_and_cluster._permuted_rmsd_pair(
            self.X, self.Y, self.perms, self.n_cg)
        reverse = align_and_cluster._permuted_rmsd_pair(
            self.Y, self.X, self.perms, self.n_cg)
        self.assertAlmostEqual(forward, reverse, places=5)

    def test_the_minimum_is_taken_over_every_permutation(self):
        # Any subset can only score worse, which is why nothing may be dropped.
        full = align_and_cluster._permuted_rmsd_pair(
            self.X, self.Y, self.perms, self.n_cg)
        for k in (1, 6, 12):
            subset = align_and_cluster._permuted_rmsd_pair(
                self.X, self.Y, self.perms[:k], self.n_cg)
            self.assertGreaterEqual(subset, full - 1e-6)

    def test_rmsd_covers_the_backbone_as_well_as_the_cg(self):
        # A permutation only relabels CG slots, so if the backbone were excluded
        # a near-symmetric CG would score ~0 under every permutation and the
        # choice would carry no information.
        moved = self.X.copy()
        moved[self.n_cg:] += 5.0
        self.assertGreater(
            align_and_cluster._permuted_rmsd_pair(moved, self.Y, self.perms, self.n_cg),
            align_and_cluster._permuted_rmsd_pair(self.X, self.Y, self.perms, self.n_cg))


class ButinaWiringTests(unittest.TestCase):
    """The automorphism group must actually reach the RMSD, not just be stored."""

    def setUp(self):
        rng = np.random.default_rng(7)
        self.n_cg = 5
        cg = rng.normal(size=(self.n_cg, 3)).astype(np.float32)
        bb = rng.normal(size=(3, 3)).astype(np.float32)
        # Second vdG: same CG with slots 2 and 3 relabelled. That swap is a valid
        # fragment automorphism, so the pair must merge.
        swapped = cg[[0, 1, 3, 2, 4]]
        self.data = np.stack([
            np.vstack([cg, bb]),
            np.vstack([swapped, bb]),
        ]).astype(np.float32)
        self.perms = _phosphate_perms()

    def _cluster(self, perms="default"):
        clusters = align_and_cluster.get_butina_clusters(
            self.data, 0.5, self.n_cg, build_perm_group(
                self.perms if perms == "default" else perms, self.n_cg, ["ALA"]))
        return clusters

    def test_the_relabelled_pair_merges(self):
        self.assertEqual(len(self._cluster()), 1)

    def test_without_the_group_the_same_pair_stays_apart(self):
        # Guards against the perms being accepted but silently unused.
        self.assertEqual(len(self._cluster(perms=None)), 2)


if __name__ == "__main__":
    unittest.main()
