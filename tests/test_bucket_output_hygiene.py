"""Per-bucket outputs must describe the current run, and silent drops must log.

Both failures this covers are invisible ones: a bucket whose records are all
unusable used to vanish with nothing in the log, and a bucket rerun after a
failure used to leave the previous run's npz beside its new .FAILED (counted by
_subset_output_counts) or a previous .FAILED beside a npz that is now fine.
"""
import os
import pickle
import tempfile
import unittest

import numpy as np

from ligand_vdgs.functions.clus_helpers import VDG_FIELDS
import ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs as C

GOOD_BB = [[[0., 0., 0.], [1.5, 0., 0.], [2.5, 1.0, 0.]]]   # N, CA, C


def _record(cg_coords, pdbpath):
    rec = {k: [] for k in VDG_FIELDS}
    rec['cg_coords'] = cg_coords
    rec['bbcoords'] = GOOD_BB
    rec['pdbpath'] = pdbpath
    return rec


class BucketOutputHygieneTests(unittest.TestCase):
    def test_all_dropped_bucket_logs_a_count_and_examples(self):
        with tempfile.TemporaryDirectory() as tmp:
            bucket = os.path.join(tmp, 'bucket.pkl')
            with open(bucket, 'wb') as fh:
                pickle.dump([_record([[np.nan, 0., 0.], [1.5, 0., 0.]],
                                     f'/db/{i}.pdb') for i in range(7)], fh)
            logfile = os.path.join(tmp, 'log')
            label, total, stats = C._run_one_bucket_strict(
                (bucket, 0.4, ['ALA'], ((0, 1),), tmp, logfile, 1, None))
            self.assertEqual((label, total), ('ALA', 0))
            log = open(logfile).read()
        self.assertIn('dropped 7/7', log)
        self.assertIn('/db/0.pdb', log)
        # Sampled, not dumped in full: a systematic failure is millions of rows.
        self.assertEqual(log.count('/db/'), 5)

    def test_stale_npz_and_failed_marker_are_cleared_for_that_bucket_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            size_dir = os.path.join(tmp, 'nr_vdgs', '2')
            os.makedirs(size_dir)
            for name in ('ASP_GLN.npz', 'ASP_GLN.FAILED', 'ASP_GLU.npz',
                         'ASP_GLU.FAILED'):
                open(os.path.join(size_dir, name), 'w').write('stale')
            C._clear_bucket_outputs(tmp, 2, 'ASP_GLN')
            self.assertEqual(sorted(os.listdir(size_dir)),
                             ['ASP_GLU.FAILED', 'ASP_GLU.npz'])
            # Idempotent, and silent on a subset directory that does not exist:
            # it runs before the bucket has written anything.
            C._clear_bucket_outputs(tmp, 2, 'ASP_GLN')
            C._clear_bucket_outputs(tmp, 99, 'NOPE')


class CrashDurabilityTests(unittest.TestCase):
    """A killed writer must not leave a file the readers accept, and a file that
    cannot be read must not pass for an absent one."""

    def test_interrupted_pickle_leaves_no_readable_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, C._bucket_fname('ASP'))
            real_dump = C.pickle.dump

            def die(obj, fh, protocol=None):
                fh.write(b'\x80\x05partial')   # a real truncated pickle
                raise KeyboardInterrupt          # SIGTERM from the queue

            C.pickle.dump = die
            try:
                with self.assertRaises(KeyboardInterrupt):
                    C._dump_pickle_atomic([1, 2, 3], path)
            finally:
                C.pickle.dump = real_dump
            # Nothing on the final name, and no temp file the walkers could pick
            # up (they filter on '.pkl', which the temp name must not end in).
            self.assertEqual(os.listdir(tmp), [])
            C._dump_pickle_atomic([1, 2, 3], path)
            with open(path, 'rb') as fh:
                self.assertEqual(pickle.load(fh), [1, 2, 3])

    def test_merge_raises_naming_the_corrupt_shard_instead_of_dropping_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            wdir = os.path.join(tmp, 'w0')
            os.makedirs(wdir)
            good = os.path.join(wdir, C._bucket_fname('ASP', flush_idx=0))
            with open(good, 'wb') as fh:
                pickle.dump(['a', 'b'], fh)
            bad = os.path.join(wdir, C._bucket_fname('ASP', flush_idx=1))
            with open(bad, 'wb') as fh:
                fh.write(b'\x80\x05partial')
            final = os.path.join(tmp, 'final')
            with self.assertRaises(RuntimeError) as ctx:
                C._merge_worker_dirs([wdir], final)
            self.assertIn(bad, str(ctx.exception))
            # Both shards belong to the same aa_key: silently keeping the good
            # half would be the failure mode this replaces.
            self.assertNotIn(C._bucket_fname('ASP'), os.listdir(final))

    def test_merge_keeps_every_record_across_workers_and_flushes(self):
        with tempfile.TemporaryDirectory() as tmp:
            wdirs = []
            for w in range(2):
                wdir = os.path.join(tmp, f'w{w}')
                os.makedirs(wdir)
                wdirs.append(wdir)
                for flush in (0, None):
                    with open(os.path.join(wdir, C._bucket_fname('ASP_GLN', flush)),
                              'wb') as fh:
                        pickle.dump([(w, flush)], fh)
                # A leftover temp file from a killed writer must be ignored.
                open(os.path.join(wdir, C._bucket_fname('ASP_GLN', 9) + '.99.tmp'),
                     'wb').write(b'junk')
            final = os.path.join(tmp, 'final')
            C._merge_worker_dirs(wdirs, final)
            self.assertEqual(os.listdir(final), [C._bucket_fname('ASP_GLN')])
            with open(os.path.join(final, C._bucket_fname('ASP_GLN')), 'rb') as fh:
                merged = pickle.load(fh)
            self.assertEqual(sorted(merged, key=str),
                             sorted([(0, 0), (0, None), (1, 0), (1, None)], key=str))

    def test_corrupt_npz_is_reported_not_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            size_dir = os.path.join(tmp, 'nr_vdgs', '1')
            os.makedirs(size_dir)
            np.savez_compressed(os.path.join(size_dir, 'ALA.npz'),
                                cluster_id=np.arange(3),
                                cluster_size=np.array([4, 2, 1]))
            # Truncated npz: np.load raises only when the members are read.
            raw = open(os.path.join(size_dir, 'ALA.npz'), 'rb').read()
            with open(os.path.join(size_dir, 'GLN.npz'), 'wb') as fh:
                fh.write(raw[:len(raw) // 2])
            logfile = os.path.join(tmp, 'log')
            n_nr, n_in, unreadable = C._subset_output_counts(tmp, 1, logfile)
            self.assertEqual((n_nr, n_in), (3, 7))
            self.assertEqual(unreadable, ['GLN'])
            self.assertIn('[ERROR] Unreadable npz', open(logfile).read())

    def test_failed_markers_are_read_back_off_disk(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(C._scan_failed_markers(tmp, 1), [])
            size_dir = os.path.join(tmp, 'nr_vdgs', '1')
            os.makedirs(size_dir)
            open(os.path.join(size_dir, 'ALA.npz'), 'w').write('x')
            C._write_failed_marker(tmp, 1, 'bb_GLN', 'boom')
            C._write_failed_marker(tmp, 1, 'ALA', 'boom')
            self.assertEqual(C._scan_failed_markers(tmp, 1), ['ALA', 'bb_GLN'])

    def test_preexisting_outputs_are_detected_per_requested_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(C._preexisting_bucket_outputs(tmp, (1, 2)), [])
            os.makedirs(os.path.join(tmp, 'nr_vdgs', '2'))
            open(os.path.join(tmp, 'nr_vdgs', '2', 'ASP_GLN.FAILED'), 'w').write('x')
            # A size this run is not building is not the run's business.
            self.assertEqual(C._preexisting_bucket_outputs(tmp, (1,)), [])
            self.assertEqual(C._preexisting_bucket_outputs(tmp, (1, 2)),
                             [(2, 'ASP_GLN.FAILED')])


if __name__ == '__main__':
    unittest.main()
