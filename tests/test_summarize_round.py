"""scripts/summarize_round.py -- keeping A/B experiment deaths out of "failed".

The case that decides whether this script is worth having: a job killed at 0:29
in the short arm and a job that crashed must NOT land in the same bucket. Round
one exists to find out whether the pipeline works; folding the h_rt experiment's
expected deaths into the failure count destroys that signal.
"""
import importlib.util
import os
import sys

import pytest

from ligand_vdgs.functions import utils

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    'summarize_round', os.path.join(REPO, 'scripts', 'summarize_round.py'))
summarize_round = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(summarize_round)

DONE_FRAG = '[C;!R;!H0][O;!R;D2]'
SHORT_KILLED = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D2]'
LONG_KILLED = '[c;!H0][c;!H0][c;!H0][c;!H0][n;D2]'
CRASHED = '[N;!R;D3][C;!R;H0]=[O;!R;D1]'


def _build(tmp_path, jobs):
    """jobs: [(fragment, arm, h_rt, completed, sge_log_text_or_None)]."""
    lib, logs, sge = (tmp_path / n for n in ('lib', 'logs', 'sge'))
    for d in (lib, logs, sge):
        d.mkdir(parents=True, exist_ok=True)
    lines = ['# short_arm_h_rt\t0:29:00', '# long_arm_h_rt\t6:00:00',
             'order\tscript\tfragment\tslots\th_rt\tab_arm']
    for order, (frag, arm, h_rt, completed, log_text) in enumerate(jobs):
        label = utils.smiles_to_filename(frag)
        frag_dir = lib / label
        frag_dir.mkdir(exist_ok=True)
        (frag_dir / f'{label}_log').write_text(
            'work\nJob completed.\n' if completed else 'work\nstreaming...\n')
        if log_text is not None:
            log_path = logs / (utils.smiles_to_job_name(frag) + '.o123')
            log_path.write_text(log_text)
            # Strictly increasing mtimes: sge_log_for takes the newest match, so a
            # loose prefix match would deterministically pick a LATER job's log
            # rather than winning or losing on same-second ties.
            os.utime(log_path, (1_700_000_000 + order, 1_700_000_000 + order))
        lines.append(f'{order}\t{sge}/{label}.sh\t{frag}\t16\t{h_rt}\t{arm}')
    order_path = tmp_path / 'submission_order.tsv'
    order_path.write_text('\n'.join(lines) + '\n')
    return str(order_path), str(lib), str(logs)


KILLED_LOG = 'date\nhostname\nstreaming 12000 records\n'      # no DONE, no FAILED
CRASH_LOG = 'date\nTraceback...\nFAILED (exit 1)\n'
DONE_LOG = 'date\nrunning\nDONE\n'


def _run(monkeypatch, capsys, order, lib, logs, extra=()):
    monkeypatch.setattr(sys, 'argv', [
        'summarize_round.py', '--submission-order', order,
        '--vdg-lib-dir', lib, '--log-dir', logs, *extra])
    code = summarize_round.main()
    return code, capsys.readouterr().out


def test_a_crash_and_a_kill_are_not_the_same_outcome(monkeypatch, capsys, tmp_path):
    # Both leave an incomplete fragment directory. Only the SGE log distinguishes
    # them: the template prints "FAILED (exit N)" when the wrapper returns
    # non-zero, and prints nothing at all when the shell is killed by a signal.
    order, lib, logs = _build(tmp_path, [
        (CRASHED, 'long', '6:00:00', False, CRASH_LOG),
        (LONG_KILLED, 'long', '6:00:00', False, KILLED_LOG)])
    _, out = _run(monkeypatch, capsys, order, lib, logs, extra=('--list',))
    assert 'crashed            1' in out.replace('  crashed', 'crashed')
    assert 'exit 1' in out
    assert 'killed             1' in out.replace('  killed', 'killed')


def test_exit_zero_without_a_completion_line_is_not_reported_as_completed(
        monkeypatch, capsys, tmp_path):
    # The wrapper returning 0 while never writing 'Job completed.' is an
    # inconsistency, not a success and not a kill.
    order, lib, logs = _build(tmp_path, [(CRASHED, 'long', '6:00:00', False, DONE_LOG)])
    _, out = _run(monkeypatch, capsys, order, lib, logs)
    assert '0 completed' in out
    assert 'exit 0 but no completion line' in out


def test_a_populated_directory_is_not_taken_for_completion(
        monkeypatch, capsys, tmp_path):
    # Every fragment here has a directory and a log; only one has the completion
    # line. If existence were the test, all three would read as completed.
    order, lib, logs = _build(tmp_path, [
        (DONE_FRAG, 'short', '0:29:00', True, DONE_LOG),
        (SHORT_KILLED, 'short', '0:29:00', False, KILLED_LOG),
        (LONG_KILLED, 'long', '6:00:00', False, KILLED_LOG)])
    _, out = _run(monkeypatch, capsys, order, lib, logs)
    assert '1 completed' in out


def test_a_missing_sge_log_is_not_called_a_failure(monkeypatch, capsys, tmp_path):
    order, lib, logs = _build(tmp_path, [(LONG_KILLED, 'long', '6:00:00', False, None)])
    _, out = _run(monkeypatch, capsys, order, lib, logs)
    assert '0 completed, 0 real failure' in out
    assert 'no_sge_log' in out


def test_an_empty_submission_order_is_refused(monkeypatch, capsys, tmp_path):
    order, lib, logs = _build(tmp_path, [])
    with pytest.raises(SystemExit):
        _run(monkeypatch, capsys, order, lib, logs)


def test_a_job_name_that_prefixes_another_does_not_steal_its_log(
        monkeypatch, capsys, tmp_path):
    # 185 of the 5,642 keys in the shipped estimate have a job name that is a
    # strict prefix of another's, so log lookup by prefix is a real hazard here,
    # not a hypothetical one. SHORT_PREFIX crashed and LONGER was killed; if the
    # lookup matched loosely, one would inherit the other's verdict.
    short_prefix = '[C;!R;!H0][O;!R;D2]'
    longer = '[C;!R;!H0][O;!R;D2][C;!R;!H0]'
    assert utils.smiles_to_job_name(longer).startswith(
        utils.smiles_to_job_name(short_prefix)), 'fixture no longer exercises the hazard'
    order, lib, logs = _build(tmp_path, [
        (short_prefix, 'long', '6:00:00', False, CRASH_LOG),
        (longer, 'short', '0:29:00', False, KILLED_LOG)])
    _, out = _run(monkeypatch, capsys, order, lib, logs, extra=('--list',))
    # Assert WHICH fragment got which verdict, not just the totals: with a loose
    # match the short-prefix fragment inherits the longer one's (newer) killed log
    # and the counts alone can still come out right.
    def section(name):
        marker = f'  {name}:'
        assert marker in out, f'no {name} section at all:\n{out}'
        return out.split(marker, 1)[1].split('\n\n', 1)[0]

    crashed_section = section('crashed')
    assert short_prefix in crashed_section, out
    assert longer not in crashed_section, out
    assert longer in section('killed'), out
    # Both count now: the kill is no longer excused by an arm (2026-09-12).
    assert '0 completed, 2 real failure' in out


def test_every_kill_counts_as_a_real_failure_now(monkeypatch, capsys, tmp_path):
    """Until 2026-09-12 a kill in the A/B's 'short' arm was excluded from the
    failure count, because the 0:29:00 tier expected some. That tier is retired and
    both surviving bands ask 48:00:00, so nothing should hit its wall.

    ASSERT WHICH, NOT HOW MANY: the two fragments carry different h_rt values and
    both must appear in the killed section, so an exclusion that dropped either one
    cannot hide in the total. Falsifier: restore the `!= 'short'` filter in the
    verdict and this reports 1 real failure instead of 2.
    """
    killed_short = '[C;!R;!H0][O;!R;D2]'
    killed_long = '[N;!R;D2][C;!R;!H0]'
    order, lib, logs = _build(tmp_path, [
        (killed_short, 'short', '0:29:00', False, KILLED_LOG),
        (killed_long, 'long', '6:00:00', False, KILLED_LOG)])
    _, out = _run(monkeypatch, capsys, order, lib, logs, extra=('--list',))
    marker = '  killed:'
    assert marker in out, f'no killed section at all:\n{out}'
    killed_section = out.split(marker, 1)[1].split('\n\n', 1)[0]
    assert killed_short in killed_section, out
    assert killed_long in killed_section, out
    assert '0 completed, 2 real failure' in out, out
    assert 'EXPECTED BY DESIGN' not in out, out

