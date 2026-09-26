"""Net-charge-sign partition (clus_and_deduplicate_vdgs._charge_sign)
and its downstream path plumbing (vdg_npz_utils.vdg_npz_path/CHARGE_SIGNS,
_bucket_npz_path/_preexisting_bucket_outputs's sign directory level).

Discriminating case: ANNOT_UNREADABLE (-1) collides with a real -1 formal
charge (carboxylate O-). _charge_sign must key off heavy_degree, not
formal_charge, or a genuine anion gets mis-tagged 'unreadable' and a CG with
one genuinely-unreadable atom gets mis-tagged by whatever its other atoms sum
to.
"""
import os
import tempfile
import unittest

from ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs import (
    _charge_sign, _bucket_npz_path, _preexisting_bucket_outputs, ANNOT_UNREADABLE)
from ligand_vdgs.functions.vdg_npz_utils import vdg_npz_path, CHARGE_SIGNS
from tests.vacuity import assert_discriminates

def _annot(heavy_degree, formal_charge):
    return {"heavy_degree": heavy_degree, "formal_charge": formal_charge}

class TestChargeSign(unittest.TestCase):
    def test_unambiguous_signs(self):
        # Per-fixture mapping, not a count: each case must land in the ONE sign
        # named, including net charge from mixed-sign atoms within one CG.
        self.assertEqual(_charge_sign(_annot([2, 2, 2, 2], [0, 0, 0, 0])), 'neut')
        self.assertEqual(_charge_sign(_annot([2, 2, 2, 2], [1, 0, 0, 0])), 'pos')
        self.assertEqual(_charge_sign(_annot([2, 2, 2, 2], [0, -1, 0, 0])), 'neg')
        self.assertEqual(_charge_sign(_annot([2, 2, 2, 2], [1, -2, 0, 0])), 'neg')
        self.assertEqual(_charge_sign(_annot([2, 2, 2, 2], [2, -1, 0, 0])), 'pos')

    def test_unreadable_sentinel_does_not_collide_with_real_negative_charge(self):
        # The exact bug caught in this session: a clause keyed on
        # formal_charge == -1 would call BOTH cases below 'unreadable'.
        # assert_discriminates proves ours (keyed on heavy_degree) does not.
        assert_discriminates(
            lambda a: _charge_sign(a) == 'unreadable',
            accepts=[_annot([2, ANNOT_UNREADABLE, 1], [0, ANNOT_UNREADABLE, 0])],
            rejects=[_annot([2, 2, 1], [0, 0, -1])],
            label='_charge_sign unreadable flag')
        self.assertEqual(_charge_sign(_annot([2, 2, 1], [0, 0, -1])), 'neg')

    def test_partial_unreadable_cg_is_unreadable_despite_plausible_other_atoms(self):
        # One unreadable atom makes the whole CG undecidable, even though the
        # other atoms would net to a plausible-looking charge (-1) alone.
        self.assertEqual(
            _charge_sign(_annot([2, ANNOT_UNREADABLE, 2], [1, ANNOT_UNREADABLE, -2])),
            'unreadable')

class TestSignDirectoryLayout(unittest.TestCase):
    def test_writer_and_reader_paths_agree_on_every_sign(self):
        # The writer's path (_bucket_npz_path) and the reader's (vdg_npz_path)
        # must resolve to the identical file for every sign, or a bucket the
        # writer places is one the reader can never find.
        with tempfile.TemporaryDirectory() as lib_root:
            for sign in CHARGE_SIGNS:
                self.assertEqual(
                    _bucket_npz_path(os.path.join(lib_root, 'test_frag'), 1, sign,
                                     ['ALA', 'bb']),
                    vdg_npz_path(lib_root, 'test_frag', 1, sign, 'ALA_bb'))

    def test_preexisting_bucket_outputs_sees_sign_subdirectories(self):
        # Non-vacuity: an empty library reports nothing, and a real leftover
        # bucket under a sign subdir IS seen -- the guard that stops a re-run
        # from silently leaving a stale build mixed into a fresh one.
        with tempfile.TemporaryDirectory() as vdglib_dir:
            self.assertEqual(_preexisting_bucket_outputs(vdglib_dir, [1]), [])
            os.makedirs(os.path.join(vdglib_dir, 'nr_vdgs', '1', 'pos'))
            open(os.path.join(vdglib_dir, 'nr_vdgs', '1', 'pos', 'ALA_bb.npz'), 'w').close()
            self.assertEqual(_preexisting_bucket_outputs(vdglib_dir, [1]),
                             [(1, 'pos/ALA_bb.npz')])

if __name__ == '__main__':
    unittest.main()
