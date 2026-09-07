"""Stage 2 must minimise over the same-label slot orderings Stage 1 quotients out.

Stage 1 treats two vdGs as the same pose if they match under ANY same-label slot
permutation (vdg_fp_utils.build_perm_group). If Stage 2 then compares slot 1's
flank to slot 1's flank positionally, a pair Stage 1 matched under a swap has its
flanks compared crosswise, the distance comes out too large, and one interaction
mode is split into two subgroups -- each carrying about half the real
cluster_size and cluster_num_parents.

The discriminating input is a pair that is IDENTICAL up to a slot swap: it must
land in one subgroup for a same-label bucket (ASP_ASP) and, because a mixed
bucket has no such symmetry to exploit, stay split for ARG_bb.
"""
import unittest

import numpy as np

from ligand_vdgs.functions.align_and_cluster import get_leader_clusters
from ligand_vdgs.functions.vdg_fp_utils import slot_orders

FLANK = 2
PER_SLOT = 2 * FLANK + 1          # [-2, -1, vdm, +1, +2] as flatten_* emit them

SLOT_A_SEQ = ["ALA", "GLY", "vdm", "SER", "THR"]
SLOT_B_SEQ = ["TRP", "TYR", "vdm", "PHE", "HIS"]


# Deliberately different SHAPES, not one translated copy: a translated pair can
# be brought back into register by the Kabsch fit inside _rmsd_pair, which would
# make the crosswise comparison cheap and the test non-discriminating.
SLOT_A_CAS = np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0],
                       [11.4, 0.0, 0.0], [15.2, 0.0, 0.0]], dtype=np.float32)
SLOT_B_CAS = np.array([[0.0, 20.0, 0.0], [2.0, 23.2, 0.0], [5.5, 24.0, 1.5],
                       [8.0, 21.5, 3.0], [6.5, 18.0, 4.5]], dtype=np.float32)

# Below the 0.5 a fully-mismatched flank sequence contributes on its own
# (dissimilarity 1.0 x seq_weight 0.5), so a crosswise comparison cannot pass it
# however the CA fit lands, while an aligned comparison scores 0.
THRESHOLD = 0.3


def _record(order):
    """One flattened record whose slots appear in `order` (0 = A, 1 = B)."""
    seqs = {0: SLOT_A_SEQ, 1: SLOT_B_SEQ}
    cas = {0: list(SLOT_A_CAS), 1: list(SLOT_B_CAS)}
    flat_seq, flat_ca = [], []
    for s in order:
        flat_seq += list(seqs[s])
        flat_ca += list(cas[s])
    return flat_seq, flat_ca


def _cluster(aa_parts, use_slot_orders):
    """Partition the swapped pair, with and without the slot-permutation fix."""
    a_seq, a_ca = _record((0, 1))
    b_seq, b_ca = _record((1, 0))      # identical vdG, slots written the other way
    seqs, cas = [a_seq, b_seq], [a_ca, b_ca]
    kwargs = {}
    if use_slot_orders:
        kwargs["slot_orders"] = slot_orders(aa_parts)
    return get_leader_clusters(
        zip([seqs, cas], ["flankseq", "flankbb"]), threshold=THRESHOLD,
        missing_seq_similarity=0.0, **kwargs)


class Stage2SlotPermutationTests(unittest.TestCase):
    def test_same_label_bucket_merges_a_slot_swapped_pair(self):
        parts = ("ASP", "ASP")
        self.assertEqual(sorted(slot_orders(parts)), [(0, 1), (1, 0)])
        merged = _cluster(parts, use_slot_orders=True)
        self.assertEqual(len(merged), 1,
                         f"slot-swapped pair should be one subgroup, got {merged}")
        self.assertEqual(sorted(next(iter(merged.values()))), [0, 1])

    def test_without_the_fix_the_same_pair_is_split(self):
        # Guards against the test passing for a reason other than the fix: the
        # identical input must still split when the orderings are withheld.
        split = _cluster(("ASP", "ASP"), use_slot_orders=False)
        self.assertEqual(len(split), 2,
                         "positional comparison should split the pair; if this "
                         "fails the input is not actually discriminating")

    def test_mixed_label_bucket_has_no_permutation_to_exploit(self):
        # ARG_bb admits only the identity, so the pair stays split even with the
        # fix enabled -- the fix must not merge slots that are not interchangeable.
        parts = ("ARG", "bb")
        self.assertEqual(slot_orders(parts), [(0, 1)])
        self.assertEqual(len(_cluster(parts, use_slot_orders=True)), 2)

    def test_identical_records_are_unaffected(self):
        # A pair needing no permutation must behave the same either way.
        seq, ca = _record((0, 1))
        for use in (True, False):
            kwargs = {"slot_orders": slot_orders(("ASP", "ASP"))} if use else {}
            parts = get_leader_clusters(
                zip([[seq, list(seq)], [ca, list(ca)]], ["flankseq", "flankbb"]),
                threshold=THRESHOLD, missing_seq_similarity=0.0, **kwargs)
            self.assertEqual(len(parts), 1)


if __name__ == "__main__":
    unittest.main()
