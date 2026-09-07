import os
import shutil
import tempfile
import unittest

import numpy as np

from ligand_vdgs.functions import vdg_struct_utils
from ligand_vdgs.functions.utils import (identify_mol_automorphisms,
                                         mol_from_fragment)
from ligand_vdgs.functions.vdg_npz_utils import (
    make_aa_bucket,
    build_vdg_atomgroup_from_npz,
    _MIN_ROTATION_ANCHOR_RATIO,
    _rotation_anchor_ratio,
    aa_perm_indices,
    cg_symmetry_path,
    load_cg_symmetry,
    write_cg_symmetry,
)


class RotationAnchorTests(unittest.TestCase):
    def test_collinear_points_do_not_determine_rotation(self):
        coords = np.array(
            [[-2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
        )
        self.assertEqual(_rotation_anchor_ratio(coords), 0.0)

    def test_non_collinear_points_determine_rotation(self):
        coords = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 1.5, 0.0]]
        )
        self.assertGreater(_rotation_anchor_ratio(coords), _MIN_ROTATION_ANCHOR_RATIO)

    def test_nearly_collinear_points_are_treated_as_unstable(self):
        coords = np.array(
            [[-2.0, 0.0, 0.0], [0.0, 1e-5, 0.0], [2.0, 0.0, 0.0]]
        )
        self.assertLessEqual(_rotation_anchor_ratio(coords), _MIN_ROTATION_ANCHOR_RATIO)


class CgSymmetrySidecarTests(unittest.TestCase):
    """The recorded group must survive a round trip and never be guessed."""

    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root, True)

    def test_round_trip_preserves_smarts_and_group(self):
        frag = "CC(=O)O"
        perms = identify_mol_automorphisms(mol_from_fragment(frag))
        write_cg_symmetry(os.path.join(self.root, frag), frag, perms)
        self.assertEqual(load_cg_symmetry(self.root, frag), (frag, perms))

    def test_label_may_differ_from_the_smarts_the_group_came_from(self):
        # The library dir is the -c label; deriving the group from it would give
        # the wrong answer here, which is the whole reason for the sidecar.
        frag, label = "O=P(O)(O)O", "phosphate"
        perms = identify_mol_automorphisms(mol_from_fragment(frag))
        write_cg_symmetry(os.path.join(self.root, label), frag, perms)
        smarts, loaded = load_cg_symmetry(self.root, label)
        self.assertEqual(smarts, frag)
        self.assertEqual(len(loaded), 24)
        self.assertIsNone(mol_from_fragment(label))

    def test_missing_sidecar_raises_rather_than_falling_back(self):
        # Every generation run writes the sidecar, so an absent one is a broken
        # library. Deriving the group from the directory name instead would
        # silently cluster under a different group.
        with self.assertRaises(FileNotFoundError):
            load_cg_symmetry(self.root, "no-such-frag")

    def test_rewrite_is_atomic_and_leaves_no_temp_files(self):
        # Subset sizes 1 and 2 share a fragment dir and race to write this.
        frag = "CC(=O)O"
        perms = identify_mol_automorphisms(mol_from_fragment(frag))
        fragdir = os.path.join(self.root, frag)
        write_cg_symmetry(fragdir, frag, perms)
        write_cg_symmetry(fragdir, frag, perms)
        outdir = os.path.dirname(cg_symmetry_path(self.root, frag))
        self.assertEqual(sorted(os.listdir(outdir)), ["cg_symmetry.npz"])


class AaSlotPermutationTests(unittest.TestCase):
    """Slots swap on bucket labels, never on an nr vdG's stored resnames.

    Regression tests for the hit finder having used ``nr_scrr_resname``:
    a nonstandard residue in a ``bb`` slot yielded no permutations at all (the
    nr vdG was dropped silently), and a duplicate resname invented a swap
    between a sidechain-contact slot and a backbone-contact slot.
    """

    def test_bb_slot_is_independent_of_the_residue_filling_it(self):
        # Whatever sits in the bb slot -- LEU, a second ASP, MSE, or a
        # GLY-named chromophore -- the labels differ, so nothing swaps.
        self.assertEqual(aa_perm_indices(["ASP", "bb"]), [[0, 1]])

    def test_like_labels_swap(self):
        self.assertEqual(sorted(aa_perm_indices(["bb", "bb"])), [[0, 1], [1, 0]])
        self.assertEqual(sorted(aa_perm_indices(["ASP", "ASP"])), [[0, 1], [1, 0]])

    def test_unlike_labels_do_not_swap(self):
        self.assertEqual(aa_perm_indices(["ASP", "HIS"]), [[0, 1]])

    def test_backbone_slots_swap_with_each_other_but_not_with_sidechains(self):
        # 'bb' is one role, so two backbone slots are interchangeable whatever
        # residues they came from -- which is why the per-nr-vdG resnames must
        # not be consulted when building permutations.
        su = vdg_struct_utils
        self.assertEqual(sorted(aa_perm_indices([su.BB_LABEL] * 2)),
                         [[0, 1], [1, 0]])
        self.assertEqual(aa_perm_indices([su.BB_LABEL, "SER"]), [[0, 1]])
        self.assertEqual(aa_perm_indices([su.BB_LABEL, su.NONCANONICAL_AA_LABEL]),
                         [[0, 1]])


class BucketNameTests(unittest.TestCase):
    """Bucket names join labels with '_', so no label may contain one."""

    def test_backbone_labels_round_trip_through_a_bucket_name(self):
        su = vdg_struct_utils
        for labels in ([su.BB_LABEL, "SER"],
                       [su.BB_LABEL, su.BB_LABEL],
                       [su.NONCANONICAL_AA_LABEL, "MET"]):
            bucket = make_aa_bucket(labels)
            self.assertEqual(bucket.split("_"), sorted(labels))

    def test_stored_parts_dtype_holds_the_longest_label(self):
        # aa_bucket_parts is written as U4, which fits every resname plus 'bb'
        # and 'X'. A longer label would silently truncate and then never match a
        # query label, so widen the dtype alongside adding one.
        su = vdg_struct_utils
        longest = max(su.BB_LABELS | {su.NONCANONICAL_AA_LABEL, "TRP"}, key=len)
        stored = np.asarray([longest], dtype="U4")
        self.assertEqual(str(stored[0]), longest)


class OccupancyProtocolTests(unittest.TestCase):
    """Every role marker stays in its own band and CG slot order survives a write."""

    def _atomgroup(self, n_cg, include_full_ligand=False, parent_pdb_path=None):
        rng = np.random.RandomState(0)
        ag, _ = build_vdg_atomgroup_from_npz(
            cg_coords=rng.rand(n_cg, 3) * 5,
            cg_names=[f"C{i}" for i in range(n_cg)], cg_elements=["C"] * n_cg,
            cg_seg="", cg_chain="A", cg_resnum=1,
            cg_resname="LIG", vdm_bb_coords=rng.rand(1, 3, 3) * 5,
            scrr_seg=[""], scrr_chain=["B"], scrr_resnum=[10], scrr_resname=["ASP"],
            include_full_ligand=include_full_ligand, parent_pdb_path=parent_pdb_path)
        return ag

    def test_other_role_markers_are_outside_the_cg_band(self):
        self.assertLess(vdg_struct_utils.NONCG_LIGAND_OCC,
                        vdg_struct_utils.CG_OCC_BASE)
        self.assertLess(vdg_struct_utils.VDM_OCC, vdg_struct_utils.CG_OCC_BASE)

    def test_cg_band_is_closed_below_four(self):
        # Everything from 4.0 up is reserved, so the top slot must stay under it.
        top = vdg_struct_utils.cg_slot_occupancy(
            vdg_struct_utils.CG_OCC_CAPACITY - 1)
        self.assertLess(top, 4.0)
        with self.assertRaises(ValueError):
            vdg_struct_utils.cg_slot_occupancy(vdg_struct_utils.CG_OCC_CAPACITY)

    def test_all_cg_slots_are_selected_past_slot_ten(self):
        ag = self._atomgroup(12)
        cg = vdg_struct_utils.select_cg_atoms(ag)
        self.assertEqual(len(cg), 12)
        atoms = vdg_struct_utils.sort_cg_atoms_by_slot(cg)
        self.assertEqual([a.getName() for a in atoms],
                         [f"C{i}" for i in range(12)])

    def test_duplicate_slot_occupancies_are_rejected(self):
        ag = self._atomgroup(4)
        occ = ag.getOccupancies()
        occ[1] = occ[0]
        ag.setOccupancies(occ)
        cg = vdg_struct_utils.select_cg_atoms(ag)
        self.assertIsNone(vdg_struct_utils.sort_cg_atoms_by_slot(cg))


if __name__ == "__main__":
    unittest.main()
