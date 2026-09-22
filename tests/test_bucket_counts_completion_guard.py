'''load_bucket_counts must refuse a fragment whose generation job never finished.

WHY THIS EXISTS. load_bucket_counts walks the library with a raw np.load rather
than through load_vdg_bucket, and load_vdg_bucket is the only reader that calls
Frags.check_vdg_job_status on our behalf. A fragment directory can be populated
but partial while its SGE job is still running, and this function AGGREGATES
across buckets, so a partial fragment does not fail -- it returns counts that are
low by an unknown amount. On 2026-09-11, 57 of 357 fragments in frag_lib_annot1
were partial and in flight while this code was importable and callable.

NAMED FALSIFIER: a fragment whose <cg>_log exists but lacks 'Job completed.'
must raise. Before the guard, it returned counts instead. Reverting the guard
turns test_incomplete_fragment_raises red.

THE VACUITY HAZARD, and what is done about it. The cheap version of this test
builds an incomplete fixture with no buckets in it. That passes against ANY
implementation -- including one that never reads the log -- because there is
nothing to count either way. So here the two fixtures are byte-identical in
their bucket contents and differ ONLY in the log line, and
test_fixtures_differ_only_in_the_log asserts exactly that: the incomplete
fixture, forced through with allow_incomplete=True, must return the SAME
non-empty counts as the complete one. If that assertion fails the other two
tests prove nothing, because a difference in content could explain the verdict.

ASSERT WHICH, NOT HOW MANY: the complete case checks label -> count per bucket
(ASP -> 3, GLN -> 5), not the total 8. A total of 8 is also what a swapped pair
or a double-counted bucket produces.
'''
import os

import numpy as np
import pytest

from ligand_vdgs.identify_bioisosteres.common import load_bucket_counts

# Distinct per-bucket counts so a swap or a miscount cannot hide in the total.
BUCKETS = {'ASP': 3, 'GLN': 5}


def _make_fragment(lib_dir, cg_label, completed):
    '''Build <lib>/<cg>/nr_vdgs/1/pos/{ASP,GLN}.npz plus <lib>/<cg>/<cg>_log.
    DR-61: buckets live one level deeper, under a charge-sign subdirectory.

    The two variants differ ONLY in whether the log carries 'Job completed.';
    every bucket is written identically, which is what makes the guard the only
    possible source of a difference in behaviour.
    '''
    frag_dir = os.path.join(lib_dir, cg_label)
    npz_dir = os.path.join(frag_dir, 'nr_vdgs', '1')
    sign_dir = os.path.join(npz_dir, 'pos')
    os.makedirs(sign_dir)
    for label, n in BUCKETS.items():
        np.savez(os.path.join(sign_dir, f'{label}.npz'),
                 cluster_id=np.arange(n, dtype=np.int32),
                 cluster_num_parents=np.ones(n, dtype=np.int32))
    log = f'{"=" * 79}\nsome generation output\n'
    if completed:
        log += 'Job completed.\nTotal job time: 0 h, 1 mins, and 2 secs.\n'
    with open(os.path.join(frag_dir, f'{cg_label}_log'), 'w') as fh:
        fh.write(log)
    return npz_dir


def test_completed_fragment_counts_each_bucket(tmp_path):
    npz_dir = _make_fragment(str(tmp_path), 'CG_done', completed=True)
    counts, extras = load_bucket_counts(npz_dir, bb_mode='off')
    # WHICH label got WHICH count, not just the total.
    assert counts == BUCKETS, counts
    assert extras == {'X': 0, 'noncanonical_bb': 0}, extras


def test_incomplete_fragment_raises(tmp_path):
    npz_dir = _make_fragment(str(tmp_path), 'CG_partial', completed=False)
    with pytest.raises(ValueError, match='incomplete'):
        load_bucket_counts(npz_dir, bb_mode='off')


def test_fixtures_differ_only_in_the_log(tmp_path):
    '''The non-vacuity clause. Without this, the two tests above could both pass
    because the incomplete fixture happened to be empty or unreadable, which
    would make the guard untested rather than working.'''
    done = _make_fragment(str(tmp_path / 'a'), 'CG_x', completed=True)
    partial = _make_fragment(str(tmp_path / 'b'), 'CG_x', completed=False)

    forced, _ = load_bucket_counts(partial, bb_mode='off', allow_incomplete=True)
    reference, _ = load_bucket_counts(done, bb_mode='off')

    # Non-empty, so there was genuinely something to count in the partial case.
    assert forced, 'incomplete fixture counted to nothing; the guard is untested'
    # Identical, so the raise above can only have come from the log line.
    assert forced == reference, (forced, reference)


def test_missing_log_also_raises(tmp_path):
    '''check_vdg_job_status returns False for an absent log as well as a
    log without the marker; a fragment whose job died before writing one must
    not be counted either.'''
    npz_dir = _make_fragment(str(tmp_path), 'CG_nolog', completed=False)
    os.remove(os.path.join(str(tmp_path), 'CG_nolog', 'CG_nolog_log'))
    with pytest.raises(ValueError, match='incomplete'):
        load_bucket_counts(npz_dir, bb_mode='off')


def test_support_sums_cluster_num_parents_not_cluster_count(tmp_path):
    '''Support is sum(cluster_num_parents), not len(cluster_id): a fragment
    with 2 clusters whose parent counts are [1, 4] must report 5, not 2.
    FALSIFIER: counting cluster rows instead of summing parents returns 2.'''
    frag_dir = os.path.join(str(tmp_path), 'CG_parents')
    npz_dir = os.path.join(frag_dir, 'nr_vdgs', '1')
    os.makedirs(os.path.join(npz_dir, 'pos'))
    np.savez(os.path.join(npz_dir, 'pos', 'ASP.npz'),
             cluster_id=np.arange(2, dtype=np.int32),
             cluster_num_parents=np.asarray([1, 4], dtype=np.int32))
    with open(os.path.join(frag_dir, 'CG_parents_log'), 'w') as fh:
        fh.write('Job completed.\n')
    counts, _ = load_bucket_counts(npz_dir, bb_mode='off')
    assert counts == {'ASP': 5}, counts
    assert counts['ASP'] != 2, 'counted clusters instead of summing parents'

def test_unexpected_directory_shape_is_rejected_not_guessed(tmp_path):
    '''The guard derives <cg_label> and <vdg_lib_dir> from npz_dir by path
    arithmetic. session_C flagged that this silently assumes a shape: a caller
    passing a differently-nested directory would get a cg_label naming the wrong
    thing, and would then be told its job is "incomplete" for a reason unrelated
    to any job. So the shape is asserted rather than trusted.

    Falsifier: without the check, this input reports an incomplete job instead
    of an unusable path, which sends the reader looking for a dead SGE job that
    does not exist.
    '''
    odd = tmp_path / 'lib' / 'CG_x' / 'not_nr_vdgs' / '1'
    odd.mkdir(parents=True)
    np.savez(str(odd / 'ASP.npz'), cluster_id=np.arange(3, dtype=np.int32))
    with pytest.raises(ValueError, match='Expected npz_dir of the form'):
        load_bucket_counts(str(odd), bb_mode='off')
