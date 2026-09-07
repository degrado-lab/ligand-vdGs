"""The 'streamed nothing' path must be loud without failing the job.

A fragment whose streaming produced no records used to exit 0 and get a
'Job completed.' line, making a wrong --pdb-dir indistinguishable from a
fragment that genuinely has no vdGs. It now exits EXIT_NO_VDGS, which the
wrapper turns into a warning and, crucially, withholds the completion marker.
"""
import os
import re
import tempfile
import unittest

from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as pipeline
from ligand_vdgs.generate_vdgs import vdg_generation_wrapper as wrapper
from ligand_vdgs.functions.Frags import check_vdg_job_status


class NoVdgsExitContractTests(unittest.TestCase):

    def test_exit_code_is_distinct_from_success_and_from_failure(self):
        # Must not collide with 0 (success) or 1 (the hard-fail path), or the
        # wrapper cannot tell the three apart.
        self.assertNotIn(pipeline.EXIT_NO_VDGS, (0, 1))

    def test_wrapper_imports_the_same_constant(self):
        self.assertEqual(wrapper.EXIT_NO_VDGS, pipeline.EXIT_NO_VDGS)

    def test_guard_fires_on_any_cause_not_only_element_mismatch(self):
        # The old guard was `if stream_skips.get("cg_elements_mismatch") and not
        # _streamed_any`, so a wrong --pdb-dir (which yields no skip reasons at
        # all) fell straight through to an empty library and exit 0.
        src = open(pipeline.__file__).read()
        self.assertIn('if not _streamed_any:', src)
        self.assertNotIn(
            'if stream_skips.get("cg_elements_mismatch") and not _streamed_any:',
            src)

    def test_a_log_without_the_marker_reads_as_unfinished(self):
        # The whole point of withholding the marker: the post-build sweep and
        # check_vdg_job_status must both still call the fragment incomplete.
        with tempfile.TemporaryDirectory() as d:
            cg = 'cn(c)[O;!R]'
            frag_dir = os.path.join(d, cg)
            os.makedirs(frag_dir)
            log = os.path.join(frag_dir, f'{cg}_log')
            with open(log, 'w') as fh:
                fh.write('[WARNING] Streamed 0 vdG records for CG cn(c)[O;!R]: '
                         'no environment produced a single record.\n')
            self.assertFalse(check_vdg_job_status(cg, d))
            with open(log, 'a') as fh:
                fh.write('=' * 79 + '\nJob completed.\n')
            self.assertTrue(check_vdg_job_status(cg, d))


    def test_warning_text_never_contains_the_completion_marker(self):
        """The message must not defeat the mechanism it describes.

        Frags.check_vdg_job_status substring-searches the whole log for
        'Job completed.', so a warning that quotes the phrase makes an empty
        fragment read as finished.
        """
        src = open(pipeline.__file__).read()
        emitted = re.findall(r'err_text = \((.*?)\)\n', src, flags=re.S)
        self.assertTrue(emitted, 'no err_text blocks found; test is stale')
        for block in emitted:
            body = '\n'.join(l for l in block.split('\n')
                              if not l.lstrip().startswith('#'))
            self.assertNotIn('Job completed.', body)

    def test_an_empty_fragment_log_does_not_read_as_complete(self):
        # End-to-end on the real reader, using the exact text the code emits.
        with tempfile.TemporaryDirectory() as d:
            cg = 'cn(c)[O;!R]'
            os.makedirs(os.path.join(d, cg))
            with open(os.path.join(d, cg, f'{cg}_log'), 'w') as fh:
                fh.write("[WARNING] Streamed 0 vdG records for CG cn(c)[O;!R]: no "
                         "environment produced a single record. nr_vdgs/ is EMPTY and "
                         "this fragment is NOT complete, so this log deliberately "
                         "carries no completion marker.\n")
            self.assertFalse(check_vdg_job_status(cg, d))


class CountOneArityTests(unittest.TestCase):
    """_count_one must return (counts, read_failures, unreadable) on EVERY path.

    The early return for a structure with no readable ligands was missed when the
    read-failure tally was added, and a 40-structure smoke sample happened to
    contain no such structure -- so the happy path passed while the real 3000-
    structure run would have died on the first ligand-less file.

    The third field also has to distinguish an unreadable file (None) from a
    structure that genuinely has no ligands ({}). Conflating them leaves dead
    files in the denominator of the extrapolation, which biases every count
    downward and under-requests resources -- jobs killed at h_rt.
    """

    def _call(self, ligands):
        from ligand_vdgs.generate_vdgs import estimate_frag_cost as efc
        efc.add_vdg_miner_paths()
        import cg as cg_mod
        real = cg_mod.read_ligand_blocks
        cg_mod.read_ligand_blocks = lambda *a, **k: ligands
        try:
            efc._init_worker(['[C;!R][O;!R]'])
            return efc._count_one('/nonexistent.pdb')
        finally:
            cg_mod.read_ligand_blocks = real

    def test_no_ligands_returns_a_three_tuple(self):
        for ligands in (None, {}, []):
            with self.subTest(ligands=ligands):
                result = self._call(ligands)
                self.assertIsInstance(result, tuple)
                self.assertEqual(len(result), 3)
                counts, failures, unreadable = result
                self.assertEqual((counts, failures), ({}, 0))
                # None means the file could not be read after retries; an empty
                # mapping means it was read and had no ligands.
                self.assertEqual(unreadable, ligands is None)

    def test_aggregator_survives_ligandless_structures(self):
        # The real failure: _tally unpacks the tuple, so an early return of the
        # wrong arity raised "not enough values to unpack" and killed the whole
        # estimate pass.
        from ligand_vdgs.generate_vdgs import estimate_frag_cost as efc
        struct_hits, occurrences, failures = [0], [0], [0]

        def tally(result):
            matched, n_failed, _unreadable = result
            failures[0] += n_failed
            for i, n in matched.items():
                struct_hits[i] += 1
                occurrences[i] += n

        tally(self._call(None))
        tally(({0: 3}, 1, False))
        self.assertEqual((struct_hits[0], occurrences[0], failures[0]), (1, 3, 1))


class PoolDeathContractTests(unittest.TestCase):

    def test_both_pools_use_an_executor_that_reports_a_killed_worker(self):
        # mp.Pool respawns an OOM-killed worker and never reports the lost task,
        # so imap_unordered blocks forever and the job burns its full h_rt.
        src = open(pipeline.__file__).read()
        # Calls, not mentions: the surrounding comments explain the swap and
        # legitimately name the old API.
        code = '\n'.join(l for l in src.split('\n')
                         if not l.lstrip().startswith('#'))
        self.assertNotIn('imap_unordered', code)
        self.assertNotIn('.Pool(', code)
        self.assertEqual(src.count('ProcessPoolExecutor('), 2)
        self.assertEqual(src.count('except BrokenProcessPool as e:'), 2)

    def test_broken_process_pool_is_imported_not_attribute_accessed(self):
        # concurrent.futures exposes only public class names via __getattr__ on
        # 3.10, so `concurrent.futures.process.X` raises AttributeError exactly
        # when a worker has died.
        src = open(pipeline.__file__).read()
        self.assertNotIn('concurrent.futures.process.', src)
        self.assertIn('from concurrent.futures.process import BrokenProcessPool',
                      src)


if __name__ == '__main__':
    unittest.main()
