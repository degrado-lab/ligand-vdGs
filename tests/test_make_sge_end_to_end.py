"""make_sge_scripts_for_frags.main() end to end, on a synthetic library.

Why this exists: file-level ownership across concurrent sessions cannot stop one
session changing a function another session calls. On 2026-09-11 this script sat
broken against a landed `select_fragments` -- `instance_counts=` had become
`support=` and the threshold unit had changed from CG occurrences to distinct
parent biounits -- and nothing failed until fleet generation would have. A test
that only compared parameter NAMES would not have caught a keyword-only change
either, so these cases run the real call.
"""
import json
import os
import pathlib
import pickle
import sys

import pytest

from ligand_vdgs.functions import Frags, utils
from ligand_vdgs.functions.db_identity import identity_of
from ligand_vdgs.generate_vdgs import make_sge_scripts_for_frags as mk

# Carbons carry !H0/H0, heteroatoms carry D<n>: the two share one isotope field
# and split_bracket_annotations rejects an atom carrying both.
ACID = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D2]'
PYRIDINE = '[c;!H0][c;!H0][c;!H0][c;!H0][n;D2]'
ETHER = '[C;!R;!H0][O;!R;D2]'
FRAGMENTS = [ACID, PYRIDINE, ETHER]
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Cost is now estimated inline by main() rather than read from a file, so every
# case fakes estimate_fragment_counts. main() discards its structures return
# value (only occurrences drive tiering), so STRUCTURES is just a filler.
STRUCTURES = 900

def _fake_estimate(occurrences, sample_scale=1.0):
    """A stand-in for estimate_fragment_counts, keyed by fragment SMARTS."""
    def fake(fragments, pdb_dir, sample_size=None, num_procs=1, seed=0):
        mk._LAST_SAMPLE_SCALE[0] = sample_scale
        return ({f: STRUCTURES for f in fragments},
                {f: occurrences.get(f, STRUCTURES) for f in fragments})
    return fake

def _setup(tmp_path, support=500, shipped_template=False):
    """A complete set of inputs for one main() run. Returns an argv list."""
    for name in ('lib', 'sge', 'logs', 'pdb'):
        (tmp_path / name).mkdir(parents=True, exist_ok=True)

    dict_path = tmp_path / 'frags.pkl'
    with open(dict_path, 'wb') as handle:
        pickle.dump({'key_schema': Frags.KEY_SCHEMA,
                     'frags': {'CCOO': {ACID: ['LIG']},
                               'CCCCN': {PYRIDINE: ['PYR']},
                               'CO': {ETHER: ['MET']}},
                     'support_pooled': {f: support for f in FRAGMENTS},
                     'db_identity': identity_of(str(tmp_path / 'pdb'))}, handle)

    if shipped_template:
        # The real template, so a regression in the shipped file is caught here
        # rather than only in a fixture that happens to mirror it.
        template = pathlib.Path(REPO) / 'resources' / 'frag_sge_template.sh'
    else:
        template = tmp_path / 'template.sh'
        template.write_text('# header\n#!/bin/bash\n#$ -pe smp $NUM_PROCS\n'
                            '#$ -l h_rt=$RUN_TIME\n'
                            'python $WRAPPER -s $SMILES -c $CG '
                            '-o $OUTPUT_DIR\n')

    return ['make_sge_scripts_for_frags.py',
            '--frags-dict', str(dict_path),
            '--template', str(template),
            '--sge-out-dir', str(tmp_path / 'sge'),
            '--vdg-lib-dir', str(tmp_path / 'lib'),
            '--log-dir', str(tmp_path / 'logs'),
            '--pdb-dir', str(tmp_path / 'pdb'),
            '--max-size', '5',
            '--min-support', '100']

def _run(monkeypatch, argv, occurrences=None, sample_scale=1.0):
    occurrences = occurrences or {f: STRUCTURES for f in FRAGMENTS}
    monkeypatch.setattr(mk, 'estimate_fragment_counts', _fake_estimate(occurrences, sample_scale))
    monkeypatch.setattr(sys, 'argv', argv)
    mk.main()

def _scripts(tmp_path):
    return sorted(f for f in os.listdir(tmp_path / 'sge') if f.endswith('.sh'))

def _manifest(tmp_path):
    """submission_order.tsv as a list of dicts, in submission order."""
    lines = [line for line in
             (tmp_path / 'sge' / 'submission_order.tsv').read_text().splitlines()
             if line and not line.startswith('#')]
    return [dict(zip(lines[0].split('\t'), line.split('\t'))) for line in lines[1:]]

def _finish(tmp_path, smiles):
    """Write the 'Job completed.' line a finished fragment leaves behind."""
    label = utils.smiles_to_filename(smiles)
    frag = tmp_path / 'lib' / label
    frag.mkdir(parents=True, exist_ok=True)
    (frag / f'{label}_log').write_text('work\nJob completed.\n')

def _kill(tmp_path, smiles):
    """The directory an h_rt kill leaves: populated, no completion line."""
    label = utils.smiles_to_filename(smiles)
    frag = tmp_path / 'lib' / label
    frag.mkdir(parents=True, exist_ok=True)
    (frag / f'{label}_log').write_text('work\nProcessing 65598 PDBs...\n')
    (frag / 'nr_vdgs').mkdir(exist_ok=True)

def test_a_full_run_writes_one_script_per_selected_fragment(monkeypatch, tmp_path):
    # The signature-compatibility case: this calls the real select_fragments with
    # the real keywords. A renamed, removed or keyword-only parameter fails here.
    _run(monkeypatch, _setup(tmp_path))
    assert len(_scripts(tmp_path)) == len(FRAGMENTS)
    assert os.path.isfile(tmp_path / 'lib' / 'library_provenance.json')

def test_support_below_the_threshold_selects_nothing(monkeypatch, tmp_path):
    # Vacuity guard on the case above: if `support` were ignored -- which is what
    # passing the wrong quantity would look like -- every fragment would still be
    # selected and the previous test would pass for the wrong reason.
    _run(monkeypatch, _setup(tmp_path, support=3))
    assert _scripts(tmp_path) == []

def test_a_stale_dict_is_refused_rather_than_silently_built(monkeypatch, tmp_path):
    argv = _setup(tmp_path)
    dict_path = argv[argv.index('--frags-dict') + 1]
    with open(dict_path, 'rb') as handle:
        payload = pickle.load(handle)
    payload.pop('key_schema')
    with open(dict_path, 'wb') as handle:
        pickle.dump(payload, handle)
    with pytest.raises(ValueError, match='key_schema'):
        _run(monkeypatch, argv)
    assert _scripts(tmp_path) == []

def test_rerunning_into_a_used_directory_needs_resume(monkeypatch, tmp_path):
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    with pytest.raises(FileExistsError):
        _run(monkeypatch, argv)

def test_resume_skips_finished_and_rewrites_the_rest(monkeypatch, tmp_path):
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    for name in _scripts(tmp_path):
        os.remove(tmp_path / 'sge' / name)
    (tmp_path / 'sge' / 'leftover.txt').write_text('makes the directory non-empty\n')

    _finish(tmp_path, ACID)
    _kill(tmp_path, PYRIDINE)          # ETHER never started
    _run(monkeypatch, argv + ['--resume'])

    written = _scripts(tmp_path)
    assert utils.smiles_to_filename(ACID) + '.sh' not in written
    assert sorted(written) == sorted(
        utils.smiles_to_filename(f) + '.sh' for f in (PYRIDINE, ETHER))
    rows = _manifest(tmp_path)
    assert len(rows) == 2
    assert all(os.path.isfile(row['script']) for row in rows)

def test_resume_leaves_a_partial_directory_alone(monkeypatch, tmp_path):
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    for name in _scripts(tmp_path):
        os.remove(tmp_path / 'sge' / name)
    _kill(tmp_path, PYRIDINE)
    _run(monkeypatch, argv + ['--resume'])
    assert 'rm -rf' not in (tmp_path / 'sge' /
        (utils.smiles_to_filename(PYRIDINE) + '.sh')).read_text()
    # and the leftover output is still there for the user to inspect
    assert os.path.isdir(tmp_path / 'lib' / utils.smiles_to_filename(PYRIDINE))

def test_na_sample_scale_does_not_crash(monkeypatch, capsys, tmp_path):
    # estimate_fragment_counts leaves _LAST_SAMPLE_SCALE at None whenever it
    # cannot compute a scale; main() must warn, not crash, and still proceed.
    _run(monkeypatch, _setup(tmp_path), sample_scale=None)
    assert 'no sample_scale' in capsys.readouterr().out
    assert len(_scripts(tmp_path)) == len(FRAGMENTS)

def test_resume_preserves_an_alias_a_top_up_merged_in(monkeypatch, tmp_path):
    # A3: --resume's `aliases` is only this run's full-selection recompute, which
    # cannot reproduce an alias an intervening --include-only top-up added for an
    # out-of-threshold charged variant. A plain overwrite here silently drops it
    # even though the top-up's vdGs are still in the library.
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    alias_path = tmp_path / 'lib' / 'fragment_aliases.tsv'
    with open(alias_path, 'a') as handle:
        handle.write('[C;!R][Cl;!R;D1]\tTOPPED_UP_REP\tcharge\n')
    for name in _scripts(tmp_path):
        os.remove(tmp_path / 'sge' / name)
    _finish(tmp_path, ACID)
    _run(monkeypatch, argv + ['--resume'])
    assert 'TOPPED_UP_REP' in alias_path.read_text()

def test_resume_does_not_rewrite_the_original_provenance(monkeypatch, tmp_path):
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    prov = tmp_path / 'lib' / 'library_provenance.json'
    original = json.loads(prov.read_text())
    for name in _scripts(tmp_path):
        os.remove(tmp_path / 'sge' / name)
    _finish(tmp_path, ACID)
    _run(monkeypatch, argv + ['--resume'])
    assert json.loads(prov.read_text()) == original

def test_tiers_key_on_occurrences_not_structures(monkeypatch, tmp_path):
    # DR-7. Every fragment shares one structure count (STRUCTURES, faked constant
    # and discarded by main()), so the tiers can only come apart if the OCCURRENCE
    # return value is what keys them.
    cheap, dear = mk.TIER_1[0] // 4, mk.RESOURCE_TIERS[-2][0] * 10
    _run(monkeypatch, _setup(tmp_path), occurrences={ACID: cheap, PYRIDINE: dear, ETHER: cheap})
    acid = (tmp_path / 'sge' / (utils.smiles_to_filename(ACID) + '.sh')).read_text()
    pyr = (tmp_path / 'sge' / (utils.smiles_to_filename(PYRIDINE) + '.sh')).read_text()
    # Positive on both sides: 'not in' alone would also pass on an empty file.
    assert f'h_rt={mk.TIER_1[2]}' in acid
    assert f'h_rt={mk.RESOURCE_TIERS[-1][2]}' in pyr
    # h_rt alone no longer discriminates: both bands ask '48:00:00' (2026-09-12), so
    # a resources_for() collapsed to a constant would still pass the two asserts
    # above. -pe smp is the axis that still varies between tiers (10 vs 20) -- key
    # on it too so a constant-resources regression actually fails this test.
    assert mk.TIER_1[1] != mk.RESOURCE_TIERS[-1][1], 'tiers no longer differ in slots either'
    assert f'-pe smp {mk.TIER_1[1]}' in acid
    assert f'-pe smp {mk.RESOURCE_TIERS[-1][1]}' in pyr

def test_manifest_orders_widest_jobs_first(monkeypatch, tmp_path):
    # The durable monitor submits this manifest directly. The old fragment-order
    # output [10, 20, 20] would strand wide reservations behind narrow jobs.
    _run(monkeypatch, _setup(tmp_path),
         occurrences={ACID: 900, PYRIDINE: 30_000, ETHER: 25_000})
    rows = _manifest(tmp_path)
    assert [int(row['slots']) for row in rows] == [20, 20, 10]
    # Within the tied 20-slot tier, the pricier fragment (higher occurrence count)
    # bubbles up first instead of falling back to alphabetical script path.
    assert [row['script'] for row in rows[:2]] == [
        str(tmp_path / 'sge' / (utils.smiles_to_filename(f) + '.sh'))
        for f in (PYRIDINE, ETHER)]
    assert [int(row['order']) for row in rows] == list(range(len(rows)))

def test_resume_refuses_a_changed_max_size(monkeypatch, tmp_path):
    # The silent-mix case DR-6 makes fatal: a resume at a different --max-size
    # would emit jobs mining a second vocabulary into the existing library, and
    # because a resume deliberately does not rewrite the provenance record,
    # nothing afterwards would show it.
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    for name in _scripts(tmp_path):
        os.remove(tmp_path / 'sge' / name)
    _finish(tmp_path, ACID)
    changed = list(argv)
    changed[changed.index('--max-size') + 1] = '4'
    # Matched on the reason: without this the case would also pass if --max-size 4
    # merely filtered every fragment out and main() exited for an unrelated cause.
    with pytest.raises(SystemExit, match='--max-size'):
        _run(monkeypatch, changed + ['--resume'])
    assert _scripts(tmp_path) == []

def test_min_support_is_required_for_a_build_but_not_for_a_top_up(
        monkeypatch, tmp_path):
    argv = _setup(tmp_path)
    without = [a for a in argv]
    index = without.index('--min-support')
    del without[index:index + 2]
    with pytest.raises(SystemExit, match='min-support'):
        _run(monkeypatch, without)
    # --include-only applies no threshold, so the flag is genuinely unused there and
    # must not be demanded. Seed current provenance first: old/unidentified libraries
    # are deliberately unsupported.
    _run(monkeypatch, argv)
    for path in (tmp_path / 'sge').iterdir():
        path.unlink()
    _run(monkeypatch, without + ['--include-only', '--include', ACID])
    assert _scripts(tmp_path) == [utils.smiles_to_filename(ACID) + '.sh']

def test_slot_hour_ceiling_is_arithmetically_right(monkeypatch, capsys, tmp_path):
    # The coordinator sizes the overnight round from this number, so it has to be
    # right rather than plausible. Occurrences are chosen to put one fragment in
    # each of three known bands; the expected total is computed from the tier table
    # rather than hard-coded, but the arithmetic is spelled out independently of
    # main()'s own expression.
    # Derived for ANY number of bands: one fragment inside the first band, one
    # just past its boundary, one far past the last finite boundary. The table
    # collapsed from four bands to two on 2026-09-12, and a `t0, t1, t2 = ...`
    # unpack here failed on the shape rather than on the arithmetic it tests.
    finite = [u for u in (t[0] for t in mk.RESOURCE_TIERS) if u != float('inf')]
    _run(monkeypatch, _setup(tmp_path),
         occurrences={ACID: finite[0] // 2, PYRIDINE: finite[0], ETHER: finite[-1] * 10})
    out = capsys.readouterr().out

    rows = _manifest(tmp_path)
    expected = sum(int(row['slots']) * mk._h_rt_to_hours(row['h_rt']) for row in rows)
    assert f'TOTAL {expected:,.0f} slot-hours' in out, (expected, out)
    peak = sum(int(row['slots']) for row in rows)
    assert f'peak {peak:,} slots' in out, (peak, out)

    # Non-vacuity: the fragments must really span more than one band, or the sum
    # would not exercise the per-band arithmetic at all. With two bands the most
    # distinct (slots, h_rt) pairs available is two, so assert against the table
    # rather than a literal 3.
    assert len({(row['slots'], row['h_rt']) for row in rows}) == min(3, len(mk.RESOURCE_TIERS)), rows
    # and the total must not coincidentally equal the single largest band
    largest = max(int(row['slots']) * mk._h_rt_to_hours(row['h_rt']) for row in rows)
    assert expected > largest, (expected, largest)

def test_a_library_in_another_vocabulary_is_refused(monkeypatch, tmp_path):
    # The 2026-09-11 near-miss: ~/docking/frag_lib held 673 pre-annotation fragment
    # directories plus a library_provenance.json with NO key_schema field. A full
    # build there would have overwritten both provenance files and pointed
    # annot-1 jobs at them, with no guard anywhere in that path.
    argv = _setup(tmp_path)
    lib = tmp_path / 'lib'
    (lib / 'library_provenance.json').write_text(json.dumps(
        {'frags_dict': 'old.pkl', 'frags_dict_sha256': 'd04d5f', 'max_size': 5,
         'min_instances': 250}))            # exactly the shipped old record's shape
    (lib / 'ccccn').mkdir()                 # an old-vocabulary fragment directory
    with pytest.raises(SystemExit, match='vocabulary mismatch'):
        _run(monkeypatch, argv)
    # Nothing may be written, including the provenance the run would have clobbered.
    assert _scripts(tmp_path) == []
    assert json.loads((lib / 'library_provenance.json').read_text())['frags_dict'] \
        == 'old.pkl'

def test_a_library_in_the_same_vocabulary_is_accepted(monkeypatch, tmp_path):
    # Discriminating pair for the case above: identical situation except the
    # recorded schema matches. If the guard refused any existing provenance, this
    # would fail and --resume would be unusable.
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    prov = json.loads((tmp_path / 'lib' / 'library_provenance.json').read_text())
    assert prov['key_schema'] == Frags.KEY_SCHEMA, prov
    for name in _scripts(tmp_path):
        os.remove(tmp_path / 'sge' / name)
    _finish(tmp_path, ACID)
    _run(monkeypatch, argv + ['--resume'])           # must not raise
    assert _scripts(tmp_path)

def test_the_vocabulary_guard_also_covers_a_top_up(monkeypatch, tmp_path):
    # --include-only writes into an existing library by design, so it is the path
    # where a vocabulary mix is easiest to cause and where check_provenance's
    # frags_dict comparison is deliberately only a warning.
    argv = _setup(tmp_path)
    lib = tmp_path / 'lib'
    _run(monkeypatch, argv)
    for path in (tmp_path / 'sge').iterdir():
        path.unlink()
    provenance = json.loads((lib / 'library_provenance.json').read_text())
    provenance['key_schema'] = 'older-0'
    (lib / 'library_provenance.json').write_text(json.dumps(provenance))
    with pytest.raises(SystemExit, match='vocabulary mismatch'):
        _run(monkeypatch, argv + ['--include-only', '--include', ACID])
    assert _scripts(tmp_path) == []

def _wrapper_invocation(script_text):
    """The `python ...` line a generated script runs, as a token list."""
    for line in script_text.splitlines():
        stripped = line.strip()
        if stripped.startswith('python ') and not stripped.startswith('#'):
            return stripped.split()
    raise AssertionError(f'no python invocation in script:\n{script_text}')

def test_generated_scripts_invoke_the_wrapper_by_absolute_path(
        monkeypatch, tmp_path):
    # Round one, 2026-09-11: the template invoked the wrapper by a path relative to
    # the repo root, and `#$ -cwd` runs the job in the SUBMIT directory, so all 153
    # dispatched jobs died with exit 2 and no traceback when submitted from the
    # script directory. A relative path here is not a style issue.
    _run(monkeypatch, _setup(tmp_path, shipped_template=True))
    for name in _scripts(tmp_path):
        tokens = _wrapper_invocation((tmp_path / 'sge' / name).read_text())
        target = tokens[1]
        assert os.path.isabs(target), (name, target)
        assert os.path.isfile(target), (name, target)
        assert target.endswith('vdg_generation_wrapper.py'), (name, target)

def test_the_generated_wrapper_call_works_from_an_unrelated_cwd(
        monkeypatch, tmp_path):
    # The assertion above is about the string; this one is about the behaviour that
    # actually broke. Run the generated invocation's target from a directory that
    # is not the repo root and require it to start. --help is enough: the failure
    # being guarded against is "python: can't open file", which happens before any
    # argument is parsed.
    import subprocess
    _run(monkeypatch, _setup(tmp_path, shipped_template=True))
    script = (tmp_path / 'sge' / _scripts(tmp_path)[0]).read_text()
    target = _wrapper_invocation(script)[1]
    elsewhere = tmp_path / 'unrelated'
    elsewhere.mkdir()
    result = subprocess.run([sys.executable, target, '--help'], cwd=str(elsewhere),
                            capture_output=True, text=True)
    assert result.returncode == 0, (result.returncode, result.stdout[-400:],
                                    result.stderr[-400:])
    assert '--subset-sizes' in result.stdout

    # Non-vacuity: the same call by the OLD relative path must fail from here, or
    # this test would pass whatever the template contained.
    relative = os.path.join('ligand_vdgs', 'generate_vdgs',
                            'vdg_generation_wrapper.py')
    broken = subprocess.run([sys.executable, relative, '--help'],
                            cwd=str(elsewhere), capture_output=True, text=True)
    assert broken.returncode != 0, broken.stdout[-200:]

def test_the_shipped_template_has_no_relative_python_invocation():
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    text = open(os.path.join(repo, 'resources', 'frag_sge_template.sh')).read()
    tokens = _wrapper_invocation(text)
    assert tokens[1] == '$WRAPPER', tokens
    assert 'python ligand_vdgs/' not in text

def test_manifest_has_no_ab_arm_column_and_the_flag_is_gone(monkeypatch, tmp_path):
    """The A/B experiment and the 0:29:00 tier it tested were removed 2026-09-12.

    Two falsifiers: restoring the trailing `ab_arm` column reddens the column
    assertion; re-adding the --short-queue-ab flag reddens the argparse one
    (argparse exits 2 on an unrecognised argument, so ACCEPTANCE means it is back).
    """
    argv = _setup(tmp_path)
    _run(monkeypatch, argv)
    rows = _manifest(tmp_path)
    assert rows, 'no manifest rows; the column assertion below would be vacuous'
    assert set(rows[0]) == {'order', 'script', 'fragment', 'slots', 'h_rt', 'tier_count'}, rows[0]
    with pytest.raises(SystemExit):
        _run(monkeypatch, argv + ['--short-queue-ab', '20260911'])
