'''summarize_round.py must not report a RUNNING job as a failure.

WHY THIS EXISTS. The `killed` bucket means "the SGE log shows neither DONE nor
FAILED (exit N)". A log that is still being written has exactly that signature,
so before the in_flight split every running job was counted as a real failure.
On 2026-09-11 the script printed "297 completed, 60 real failure(s)" for a round
in which all 60 of those jobs were running normally and 0 had failed. The job ID
is already in the log filename as <job_name>.o<job_id>, so distinguishing the two
needs no manifest and no extra flag -- only a live job-ID set.

NAMED FALSIFIER: a fragment whose log lacks DONE/FAILED and whose job ID is still
in qstat must classify as 'in_flight'; one whose job ID is gone must classify as
'killed'. Reverting the split makes both come back 'killed' and reddens
test_which_fragment_gets_which_outcome.

TWO VACUITY HAZARDS, both deliberately armed here:

1. ASSERT WHICH, NOT HOW MANY. "one in_flight and one killed" is also what the
   swapped assignment produces. So the assertions name the fragment and its
   outcome individually.

2. THE FIXTURE MUST FORCE THE FAILURE, NOT MERELY ALLOW IT. The two logs are
   byte-identical in content, so the classification cannot come from the text --
   it can only come from the job-ID lookup. And the job names are chosen so that
   RUNNING is a strict prefix of RUNNING_LONGER, with the longer fragment's log
   given a strictly LATER mtime. `sge_log_for` matches on `prefix + '.o'`, which
   is prefix-safe; if anyone loosens it to `startswith(prefix)`, the `max(...,
   key=getmtime)` tie-break will then deterministically pick the WRONG log for
   the shorter fragment rather than coincidentally the right one. That is the
   exact trap that let an earlier version of this check pass while both
   individual verdicts were wrong (the two fixture logs shared an mtime).
'''
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ligand_vdgs.functions import utils
from scripts.summarize_round import classify, sge_log_for

# Strict prefix pair: job names are the fragment strings, so RUNNING's job name
# is a prefix of RUNNING_LONGER's. This arms vacuity hazard 2.
RUNNING = '[C;!R;!H0][O;!R;D1]'
RUNNING_LONGER = '[C;!R;!H0][O;!R;D1][C;!R;H0]'

RUNNING_JID = '5000001'
DEAD_JID = '5000002'

# No DONE, no 'FAILED (exit N)' -- the signature shared by a killed job and a
# job that is still writing. Identical for both fragments on purpose.
PARTIAL_LOG = 'Fri Sep 11 12:00:00 PDT 2026\nqb3-id100\nsome output\n'


def _build(tmp_path):
    '''Library with two INCOMPLETE fragments and one log each.

    Neither fragment gets a 'Job completed.' line, so both fall through to the
    log inspection -- otherwise the completed short-circuit would hide the
    branch under test.
    '''
    lib = tmp_path / 'lib'
    logs = tmp_path / 'logs'
    lib.mkdir()
    logs.mkdir()
    for frag in (RUNNING, RUNNING_LONGER):
        frag_dir = lib / utils.smiles_to_filename(frag)
        frag_dir.mkdir()
        label = utils.smiles_to_filename(frag)
        (frag_dir / f'{label}_log').write_text('started, never finished\n')

    short_log = logs / f'{utils.smiles_to_job_name(RUNNING)}.o{RUNNING_JID}'
    long_log = logs / f'{utils.smiles_to_job_name(RUNNING_LONGER)}.o{DEAD_JID}'
    short_log.write_text(PARTIAL_LOG)
    long_log.write_text(PARTIAL_LOG)
    # Strictly increasing mtimes: the longer (prefix-colliding) log is NEWER, so
    # a loosened prefix match would grab it for the shorter fragment.
    os.utime(short_log, (1_000_000, 1_000_000))
    os.utime(long_log, (2_000_000, 2_000_000))
    return str(lib), str(logs)


def test_which_fragment_gets_which_outcome(tmp_path):
    '''The live one is in_flight, the dead one is killed -- named individually,
    because the counts alone are identical under a swap.'''
    lib, logs = _build(tmp_path)
    live = {RUNNING_JID}  # DEAD_JID deliberately absent

    running_outcome, _ = classify(RUNNING, lib, logs, live_ids=live)
    dead_outcome, _ = classify(RUNNING_LONGER, lib, logs, live_ids=live)

    assert running_outcome == 'in_flight', running_outcome
    assert dead_outcome == 'killed', dead_outcome


def test_logs_are_identical_so_only_the_job_id_can_decide(tmp_path):
    '''Non-vacuity clause for hazard 1. If the two logs differed, the outcomes
    above could be explained by their contents rather than by the lookup.'''
    lib, logs = _build(tmp_path)
    a = open(sge_log_for(logs, RUNNING)).read()
    b = open(sge_log_for(logs, RUNNING_LONGER)).read()
    assert a == b, 'fixture logs differ; the in_flight verdict is not isolated'
    assert 'DONE' not in a and 'FAILED' not in a


def test_prefix_collision_resolves_to_the_right_log(tmp_path):
    '''Non-vacuity clause for hazard 2. The shorter job name is a strict prefix
    of the longer one and the longer log is newer, so a loose match would return
    the wrong file here.'''
    lib, logs = _build(tmp_path)
    assert utils.smiles_to_job_name(RUNNING_LONGER).startswith(
        utils.smiles_to_job_name(RUNNING)), 'fixture no longer arms the hazard'
    assert os.path.getmtime(sge_log_for(logs, RUNNING_LONGER)) > os.path.getmtime(
        sge_log_for(logs, RUNNING)), 'mtimes no longer force the wrong pick'

    assert sge_log_for(logs, RUNNING).endswith(f'.o{RUNNING_JID}')
    assert sge_log_for(logs, RUNNING_LONGER).endswith(f'.o{DEAD_JID}')


def test_without_a_live_set_nothing_is_called_in_flight(tmp_path):
    '''live_ids=None means "we could not ask qstat". The script must then fall
    back to the old, conservative labelling rather than guessing in_flight --
    and it warns, which is checked in the CLI, not here.'''
    lib, logs = _build(tmp_path)
    assert classify(RUNNING, lib, logs, live_ids=None)[0] == 'killed'
    assert classify(RUNNING_LONGER, lib, logs, live_ids=None)[0] == 'killed'


def test_empty_live_set_is_not_the_same_as_no_live_set(tmp_path):
    '''An empty set is a real answer ("nothing of yours is queued") and must
    classify as killed; None is the absence of an answer. Both land on killed
    here, so the distinction is asserted on the in_flight case instead: with the
    running job present, only the set form reports it.'''
    lib, logs = _build(tmp_path)
    assert classify(RUNNING, lib, logs, live_ids=set())[0] == 'killed'
    assert classify(RUNNING, lib, logs, live_ids={RUNNING_JID})[0] == 'in_flight'
