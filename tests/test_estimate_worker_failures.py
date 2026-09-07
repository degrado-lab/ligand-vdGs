"""One bad structure must not abort the counting pass, nor silently bias it.

The pass costs minutes and writes nothing on the way, so an uncaught worker
exception threw away the whole run. But a dropped structure that stays in the
denominator of the extrapolation biases every fragment's count downward, which
under-requests resources -- so dropped structures must leave the denominator too.

Both branches route every result through the same _consume_result, so the handler
the cluster actually runs (num_procs > 1) is the one covered here -- a spawned
worker cannot see a monkeypatched _count_one, which is why the handler is shared
rather than duplicated per branch.
"""
import os
import tempfile
import unittest

from concurrent.futures.process import BrokenProcessPool

from ligand_vdgs.generate_vdgs import estimate_frag_cost as efc

FRAG = '[C;!R][O;!R]'


def _mirror(n):
    """A mirror of `n` empty-but-present PDB files."""
    d = tempfile.mkdtemp()
    for i in range(n):
        sub = os.path.join(d, f'{i:02d}')
        os.makedirs(sub, exist_ok=True)
        open(os.path.join(sub, f'{i}.pdb'), 'w').close()
    return d


class WorkerFailureTests(unittest.TestCase):

    def setUp(self):
        self.real = efc._count_one
        self.dir = _mirror(20)

    def tearDown(self):
        efc._count_one = self.real

    def _run(self, fake, **kw):
        efc._count_one = fake
        return efc.estimate_fragment_counts(
            [FRAG], self.dir, sample_size=20, num_procs=1, **kw)

    def test_a_raising_structure_does_not_abort_the_pass(self):
        calls = {'n': 0}

        def fake(path):
            calls['n'] += 1
            if calls['n'] == 1:
                raise RuntimeError('boom')
            return {0: 1}, 0, False

        structures, occurrences = self._run(fake)
        # 19 of 20 counted, and the scale divides by 19, not 20 -- so a fragment
        # in every counted structure still extrapolates to the whole mirror.
        self.assertEqual(structures[FRAG], 20)
        self.assertEqual(occurrences[FRAG], 20)

    def test_errored_structures_leave_the_denominator(self):
        # 2 of 20 errors -- under MAX_LOST_SAMPLE_FRACTION, so the pass proceeds.
        # The 18 survivors all contain the fragment; if the errors stayed in the
        # denominator this would come out at 18/20 of the mirror, not all of it.
        calls = {'n': 0}

        def fake(path):
            calls['n'] += 1
            if calls['n'] <= 2:
                raise RuntimeError('boom')
            return {0: 1}, 0, False

        structures, _ = self._run(fake)
        self.assertEqual(structures[FRAG], 20)

    def test_losing_too_much_of_the_sample_raises_rather_than_extrapolating(self):
        # Half the sample errors: scaling the survivors up by 2x would be a
        # confident, wrong estimate, so this is a failure rather than a warning.
        calls = {'n': 0}

        def fake(path):
            calls['n'] += 1
            if calls['n'] % 2:
                raise RuntimeError('boom')
            return {0: 1}, 0, False

        with self.assertRaisesRegex(RuntimeError, 'refusing to emit an estimate'):
            self._run(fake)

    def test_an_entirely_uncountable_sample_raises(self):
        def fake(path):
            raise RuntimeError('boom')

        with self.assertRaisesRegex(RuntimeError, 'refusing to emit an all-zero'):
            self._run(fake)

    def test_unreadable_files_also_leave_the_denominator(self):
        calls = {'n': 0}

        def fake(path):
            calls['n'] += 1
            if calls['n'] <= 2:
                return {}, 0, True  # read_ligand_blocks returned None
            return {0: 1}, 0, False

        structures, _ = self._run(fake)
        self.assertEqual(structures[FRAG], 20)

    def test_a_ligandless_structure_stays_in_the_denominator(self):
        # It was read successfully and genuinely has no ligands: that is real
        # information about the mirror, not a lost sample.
        calls = {'n': 0}

        def fake(path):
            calls['n'] += 1
            if calls['n'] <= 10:
                return {}, 0, False
            return {0: 1}, 0, False

        structures, _ = self._run(fake)
        self.assertEqual(structures[FRAG], 10)


class ConsumeResultTests(unittest.TestCase):
    """_consume_result is what both branches use, so cover it directly."""

    def setUp(self):
        self.tallied = []
        self.errored = []

    def test_a_raising_result_is_recorded_not_propagated(self):
        def boom():
            raise RuntimeError('boom')

        efc._consume_result('/a.pdb', boom, self.tallied.append, self.errored)
        self.assertEqual(self.tallied, [])
        self.assertEqual(len(self.errored), 1)
        self.assertEqual(self.errored[0][0], '/a.pdb')
        self.assertIn('boom', self.errored[0][1])

    def test_a_good_result_is_tallied(self):
        efc._consume_result('/a.pdb', lambda: ({0: 2}, 0, False),
                            self.tallied.append, self.errored)
        self.assertEqual(self.tallied, [({0: 2}, 0, False)])
        self.assertEqual(self.errored, [])

    def test_a_broken_pool_propagates_rather_than_being_counted(self):
        # A dead pool is not a bad structure: continuing would extrapolate from a
        # sample that silently stopped being processed.
        def dead():
            raise BrokenProcessPool('worker died')

        with self.assertRaises(BrokenProcessPool):
            efc._consume_result('/a.pdb', dead, self.tallied.append, self.errored)
        self.assertEqual(self.errored, [])

    def test_both_branches_route_through_it(self):
        # Guards the reason it is module level: a future edit that inlines the
        # handler back into one branch leaves the production path untested.
        src = open(efc.__file__).read()
        code = '\n'.join(l for l in src.split('\n')
                          if not l.lstrip().startswith('#'))
        self.assertEqual(code.count('_consume_result('), 3)  # def + 2 call sites


class PoolDeathTests(unittest.TestCase):

    def test_the_counting_pool_reports_a_killed_worker(self):
        # mp.Pool respawns an OOM-killed worker and never reports the lost task,
        # so imap_unordered blocks forever and the job burns its full h_rt.
        src = open(efc.__file__).read()
        code = '\n'.join(l for l in src.split('\n')
                         if not l.lstrip().startswith('#'))
        self.assertNotIn('imap_unordered', code)
        self.assertNotIn('.Pool(', code)
        self.assertIn('ProcessPoolExecutor(', code)
        self.assertIn('from concurrent.futures.process import BrokenProcessPool', src)
        self.assertNotIn('concurrent.futures.process.', code)


if __name__ == '__main__':
    unittest.main()
