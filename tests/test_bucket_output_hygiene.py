"""Files a killed writer leaves behind, and files a reader cannot read, must
never pass for good output."""
import os
import tempfile
import unittest

import numpy as np

from ligand_vdgs.functions import clus_helpers as ch
import ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs as C


def _cols(biounits):
    n = len(biounits)
    return {"biounit": np.asarray(biounits, dtype="U32"),
            "cgvdmbb": np.zeros((n, 7, 3), dtype=np.float32)}


class CrashDurabilityTests(unittest.TestCase):
    def test_interrupted_shard_leaves_no_readable_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, C._bucket_fname('ASP'))
            real_savez = ch.np.savez

            def die(handle, **arrays):
                handle.write(b'PK\x03\x04partial')   # a real truncated zip
                raise KeyboardInterrupt              # SIGTERM from the queue

            ch.np.savez = die
            try:
                with self.assertRaises(KeyboardInterrupt):
                    ch.save_shard(path, _cols(['1abc']))
            finally:
                ch.np.savez = real_savez
            # Nothing on the final name, and no temp file the merge could pick
            # up (it filters on '.npz', which the temp name must not end in).
            self.assertEqual(os.listdir(tmp), [])
            ch.save_shard(path, _cols(['1abc', '2xyz']))
            self.assertEqual(list(ch.load_shard(path)['biounit']), ['1abc', '2xyz'])

    def test_merge_raises_naming_the_corrupt_shard_instead_of_dropping_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            wdir = os.path.join(tmp, 'w0')
            os.makedirs(wdir)
            ch.save_shard(os.path.join(wdir, C._bucket_fname('ASP', flush_idx=0)),
                          _cols(['1abc']))
            bad = os.path.join(wdir, C._bucket_fname('ASP', flush_idx=1))
            with open(bad, 'wb') as fh:
                fh.write(b'PK\x03\x04partial')
            final = os.path.join(tmp, 'final')
            with self.assertRaises(RuntimeError) as ctx:
                C._merge_worker_dirs([wdir], final, None, os.path.join(tmp, 'log'))
            self.assertIn(bad, str(ctx.exception))
            # Both shards belong to the same aa_key: silently keeping the good
            # half would be the failure mode this replaces.
            self.assertNotIn('ASP', os.listdir(final))

    def test_merge_keeps_every_record_across_workers_and_flushes(self):
        with tempfile.TemporaryDirectory() as tmp:
            wdirs = []
            for w in range(2):
                wdir = os.path.join(tmp, f'w{w}')
                os.makedirs(wdir)
                wdirs.append(wdir)
                for flush in (0, None):
                    ch.save_shard(os.path.join(wdir, C._bucket_fname('ASP_GLN', flush)),
                                  _cols([f'{w}abc_{flush}']))
                # A leftover temp file from a killed writer must be ignored.
                open(os.path.join(wdir, C._bucket_fname('ASP_GLN', 9) + '.99.tmp'),
                     'wb').write(b'junk')
            final = os.path.join(tmp, 'final')
            counts = C._merge_worker_dirs(wdirs, final, None, os.path.join(tmp, 'log'))
            self.assertEqual(counts, {'ASP_GLN': 4})
            self.assertEqual(os.listdir(final), ['ASP_GLN'])
            merged = ch.load_columns(os.path.join(final, 'ASP_GLN'))
            self.assertEqual(sorted(merged['biounit']),
                             ['0abc_0', '0abc_None', '1abc_0', '1abc_None'])
            self.assertEqual(merged['cgvdmbb'].shape, (4, 7, 3))

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
