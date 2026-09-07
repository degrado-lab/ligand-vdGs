import unittest

import numpy as np

from ligand_vdgs.functions import align_and_cluster
from ligand_vdgs.functions import clus_helpers
from ligand_vdgs.functions.vdg_struct_utils import (FLANK_CHAIN_BREAK,
    FLANK_MISSING, is_valid_backbone_coords)
from ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs import (
    _has_complete_stage1_coords,
)


class MissingDataClusteringTests(unittest.TestCase):
    def test_stage2_rmsd_uses_shared_finite_atom_count(self):
        x = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [np.nan, np.nan, np.nan],
        ], dtype=np.float32)
        y = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [np.nan, np.nan, np.nan],
        ], dtype=np.float32)

        # The two shared points have SSD=0.5, so RMSD=sqrt(0.5/2)=0.5.
        self.assertAlmostEqual(align_and_cluster._rmsd_pair(x, y), 0.5)

    def test_stage2_rmsd_masks_all_non_finite_rows(self):
        x = np.array([[0.0, 0.0, 0.0], [np.inf, 0.0, 0.0]], dtype=np.float32)
        y = np.array([[1.0, 2.0, 3.0], [0.0, 0.0, 0.0]], dtype=np.float32)
        self.assertAlmostEqual(align_and_cluster._rmsd_pair(x, y), 0.0)

    def test_missing_sequence_uses_configured_similarity_as_prior(self):
        reference = ["ALA", "GLY", "vdm", "SER", "THR"]
        M, B = FLANK_MISSING, FLANK_CHAIN_BREAK
        cases = [
            ([M, M, "vdm", M, M], 40.0),
            (["ALA", M, "vdm", M, M], 55.0),
            (["VAL", M, "vdm", M, M], 30.0),
            (["ALA", "GLY", "vdm", M, M], 70.0),
            # A chain break carries the same prior as an unreadable residue:
            # neither supplies a residue to match.
            ([M, B, "vdm", B, M], 40.0),
            (["ALA", B, "vdm", B, B], 55.0),
        ]
        for sequence, expected in cases:
            with self.subTest(sequence=sequence):
                similarity = clus_helpers.calc_seq_similarity(
                    reference, sequence, missing_similarity=40.0)
                self.assertAlmostEqual(similarity, expected)

    def test_noncanonical_label_is_scored_as_a_residue_not_as_absent_data(self):
        # 'X' is a vdM slot label (NONCANONICAL_AA_LABEL), never a flank marker.
        # It must not be silently absorbed into the missing-data prior.
        reference = ["ALA", "GLY", "vdm", "SER", "THR"]
        mismatch = ["X", "X", "vdm", "X", "X"]
        self.assertAlmostEqual(
            clus_helpers.calc_seq_similarity(
                reference, mismatch, missing_similarity=40.0),
            0.0)
        self.assertAlmostEqual(
            clus_helpers.calc_seq_similarity(
                mismatch, list(mismatch), missing_similarity=40.0),
            100.0)

    def test_missing_flanks_remain_clustered_with_expected_size_cutoff(self):
        reference_seq = ["ALA", "GLY", "vdm", "SER", "THR"]
        missing_seq = [FLANK_MISSING, FLANK_MISSING, "vdm",
                       FLANK_MISSING, FLANK_MISSING]
        x = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
        ], dtype=np.float32)
        y = np.array([
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
        ], dtype=np.float32)

        # RMSD=0.4 over the two shared atoms. With the neutral 40% sequence
        # prior the combined distance is 0.4 + 0.3 = 0.7, below the expected-
        # size threshold 0.5 + 0.3 = 0.8.
        assignments = align_and_cluster.get_leader_clusters(
            [([reference_seq, missing_seq], "flankseq"),
             ([x, y], "flankbb")],
            threshold=0.8,
            missing_seq_similarity=0.4,
            final_exact_medoid_pass=False,
            final_reassign_once=False,
        )
        self.assertEqual(assignments, {1: [0, 1]})
        align_and_cluster.clear_caches()

    def test_stage1_coordinate_validation_requires_full_finite_rank2_backbone(self):
        cg = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
        backbone = np.array([
            [-1.459, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.551, 1.422, 0.0],
        ], dtype=np.float32)

        self.assertTrue(is_valid_backbone_coords(backbone))
        self.assertTrue(_has_complete_stage1_coords(
            [cg, [backbone]], expected_n_cg=2, expected_num_vdms=1))

        collinear = backbone.copy()
        collinear[:, 1] = 0.0
        self.assertFalse(is_valid_backbone_coords(collinear))
        self.assertFalse(_has_complete_stage1_coords(
            [cg, [collinear]], expected_n_cg=2, expected_num_vdms=1))

        non_finite = backbone.copy()
        non_finite[0, 0] = np.nan
        self.assertFalse(_has_complete_stage1_coords(
            [cg, [non_finite]], expected_n_cg=2, expected_num_vdms=1))

        self.assertFalse(_has_complete_stage1_coords(
            [cg[:1], [backbone]], expected_n_cg=2, expected_num_vdms=1))


if __name__ == "__main__":
    unittest.main()
