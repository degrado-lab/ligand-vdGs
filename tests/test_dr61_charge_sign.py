"""DR-61: net-charge-sign partition (clus_and_deduplicate_vdgs._charge_sign)
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
import numpy as np

from ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs import (
    _charge_sign, _bucket_npz_path, _preexisting_bucket_outputs, ANNOT_UNREADABLE)
from ligand_vdgs.functions import clus_helpers, vdg_npz_utils
from ligand_vdgs.functions.vdg_npz_utils import vdg_npz_path, CHARGE_SIGNS
from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as clus
from ligand_vdgs.tools.h_class_diagnostic import read_bucket
from tests.test_bucket_schema_pass import _record
from tests.vacuity import assert_discriminates

def _write_sign_bucket(frag_dir, sign, records):
    clus._write_bucket_npz(
        frag_dir, 2, sign, ('GLY', 'ALA'), clus_helpers.records_to_columns(records),
        [clus.Subgroup(1, 1, i, np.array([i], dtype=np.int32), 0.25)
         for i in range(len(records))], '/db')

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

    def test_h_class_reader_merges_current_signs_and_refuses_flat_bucket(self):
        # Falsifier: a reader that inspects only the first sign reports one row
        # instead of the three rows written across the two current partitions.
        with tempfile.TemporaryDirectory() as lib:
            frag_dir = os.path.join(lib, 'CG')
            _write_sign_bucket(frag_dir, 'neg', [_record()])
            _write_sign_bucket(frag_dir, 'pos', [_record(), _record()])
            with open(os.path.join(frag_dir, 'CG_log'), 'w') as handle:
                handle.write('Job completed.\n')
            self.assertEqual(len(read_bucket(lib, 'CG', 2, 'GLY_ALA')['cluster_id']), 3)

            legacy = os.path.join(lib, 'legacy', 'nr_vdgs', '2')
            os.makedirs(legacy)
            open(os.path.join(legacy, 'GLY_ALA.npz'), 'wb').close()
            with self.assertRaises(vdg_npz_utils.BucketSchemaMismatch):
                read_bucket(lib, 'legacy', 2, 'GLY_ALA')

if __name__ == '__main__':
    unittest.main()
