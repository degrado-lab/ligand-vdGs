"""--resume in make_sge_scripts_for_frags: which fragments a resumed build re-runs.

The case that matters is the one a naive check gets wrong. An h_rt kill leaves a
populated fragment directory behind, so "the directory exists" and "the fragment
finished" are different facts, and only the 'Job completed.' line distinguishes
them. Every case here is built around that difference.
"""
import os

import pytest

from ligand_vdgs.functions import utils
from ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags import (
    clear_partial_snippet, output_script, partition_by_completion)

DONE = 'c1ccccc1O'          # finished
KILLED = '[C;!R][O;!R]'     # ran, was killed: directory present, no completion line
FRESH = '[N;!R][C;!R]'      # never started: no directory at all


def _library(tmp_path):
    lib = tmp_path / 'lib'
    lib.mkdir()
    for smiles, log in [(DONE, 'stuff\nJob completed.\nTotal job time: 1 h\n'),
                        (KILLED, 'stuff\nProcessing 65598 PDBs...\n')]:
        label = utils.smiles_to_filename(smiles)
        frag = lib / label
        frag.mkdir()
        (frag / f'{label}_log').write_text(log)
        (frag / 'nr_vdgs').mkdir()          # partial output, as a kill leaves it
    return str(lib)


def test_finished_skipped_killed_rerun_fresh_rerun(tmp_path):
    lib = _library(tmp_path)
    unfinished, finished, partial = partition_by_completion([DONE, KILLED, FRESH], lib)
    assert finished == [DONE]
    assert unfinished == [KILLED, FRESH]
    # Only the killed one has a leftover directory; the fresh one has nothing to clear.
    assert partial == [KILLED]


def test_a_directory_is_not_evidence_of_completion(tmp_path):
    # The discriminating pair: DONE and KILLED are indistinguishable by directory
    # existence and differ only in the log line. If the verdict for both were the
    # same, this whole feature would be doing nothing.
    lib = _library(tmp_path)
    for smiles in (DONE, KILLED):
        label = utils.smiles_to_filename(smiles)
        assert os.path.isdir(os.path.join(lib, label))
    unfinished, finished, _ = partition_by_completion([DONE, KILLED], lib)
    assert finished != unfinished
    assert DONE in finished and KILLED in unfinished


def test_completion_line_alone_flips_the_verdict(tmp_path):
    # Same fragment, same directory contents, one line appended. Nothing else may
    # decide this.
    lib = _library(tmp_path)
    label = utils.smiles_to_filename(KILLED)
    log = os.path.join(lib, label, f'{label}_log')
    before = partition_by_completion([KILLED], lib)
    with open(log, 'a') as handle:
        handle.write('Job completed.\n')
    after = partition_by_completion([KILLED], lib)
    assert before[0] == [KILLED] and before[1] == []
    assert after[0] == [] and after[1] == [KILLED]


def test_an_empty_library_reruns_everything(tmp_path):
    lib = tmp_path / 'empty'
    lib.mkdir()
    unfinished, finished, partial = partition_by_completion([DONE, FRESH], str(lib))
    assert unfinished == [DONE, FRESH]
    assert finished == [] and partial == []


def test_clear_snippet_is_guarded_and_names_one_directory(tmp_path):
    lib = _library(tmp_path)
    label = utils.smiles_to_filename(KILLED)
    snippet = clear_partial_snippet(lib, label)
    target = os.path.join(os.path.abspath(lib), label)
    # It must delete exactly one path, and only after establishing that the path
    # is non-empty, is a directory, and holds this fragment's own log.
    assert snippet.count('rm -rf') == 1
    assert f'{label}_log' in snippet
    assert '[ -n "$RESUME_DIR" ]' in snippet and '[ -d "$RESUME_DIR" ]' in snippet
    assert target in snippet
    # A bare unguarded removal of the library root would be the catastrophic bug.
    assert f'rm -rf -- "{os.path.abspath(lib)}"\n' not in snippet


def test_clear_snippet_quotes_a_hostile_label(tmp_path):
    # smiles_to_filename encodes '/' but the guard must not depend on that.
    snippet = clear_partial_snippet('/tmp/lib dir', 'frag name')
    assert "'/tmp/lib dir/frag name'" in snippet


def test_pre_run_defaults_to_empty_so_a_normal_build_is_unchanged(tmp_path):
    template = ['#!/bin/bash\n', '$PRE_RUN\n', 'python run.py -c $CG\n']
    out = tmp_path / 'sge'
    out.mkdir()
    output_script(template, FRESH, str(out), {})
    written = (out / (utils.smiles_to_filename(FRESH) + '.sh')).read_text()
    assert '$PRE_RUN' not in written
    assert 'rm -rf' not in written


def test_pre_run_is_substituted_when_supplied(tmp_path):
    template = ['#!/bin/bash\n', '$PRE_RUN\n', 'python run.py -c $CG\n']
    out = tmp_path / 'sge'
    out.mkdir()
    snippet = clear_partial_snippet('/tmp/lib', 'frag')
    output_script(template, FRESH, str(out), {'$PRE_RUN': snippet})
    written = (out / (utils.smiles_to_filename(FRESH) + '.sh')).read_text()
    assert 'rm -rf' in written
    assert '$PRE_RUN' not in written
    # The snippet must land before the wrapper runs, or it clears the output the
    # wrapper just produced.
    assert written.index('rm -rf') < written.index('python run.py')


def test_shipped_template_carries_the_hook_before_the_wrapper():
    # Anchored on the $WRAPPER placeholder, not on the wrapper's filename: the
    # template stopped naming the file directly when the invocation was made
    # absolute (round one, 2026-09-11), and a filename anchor silently stopped
    # locating the invocation at all.
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    text = open(os.path.join(repo, 'resources', 'frag_sge_template.sh')).read()
    assert '$PRE_RUN' in text
    invocation = [line for line in text.splitlines()
                  if line.strip().startswith('python ')]
    assert len(invocation) == 1, invocation
    assert '$WRAPPER' in invocation[0], invocation
    assert text.index('$PRE_RUN') < text.index(invocation[0])
