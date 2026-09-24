"""The cost manifest, not script filenames, controls SGE submission."""
import sys

import pytest

from ligand_vdgs.generate_vdgs import submit_frag_jobs_by_size as submit
from tests.vacuity import assert_discriminates

def _fleet(tmp_path):
    scripts = [tmp_path / name for name in ('z_expensive.sh', 'a_cheaper.sh', 'b_small.sh')]
    for path, slots in zip(scripts, (20, 20, 10)):
        path.write_text(f'#!/bin/bash\n#$ -pe smp {slots}\n')
    manifest = tmp_path / 'submission_order.tsv'
    manifest.write_text('# submit in this order\norder\tscript\tfragment\tslots\th_rt\ttier_count\n' +
                        ''.join(f'{i}\t{path}\tfrag{i}\t{slots}\t48:00:00\t{cost}\n'
                                for i, (path, slots, cost) in enumerate(zip(scripts, (20, 20, 10),
                                                                             (1000, 500, 10)))))
    return scripts, manifest

def test_submitter_uses_manifest_for_preview_and_qsub(tmp_path, monkeypatch, capsys):
    scripts, _ = _fleet(tmp_path)
    expected = [str(path) for path in scripts]
    assert_discriminates(lambda paths: paths == expected, [expected],
                         [sorted(expected[:2]) + expected[2:]], 'cost before filename')
    monkeypatch.setattr(sys, 'argv', ['submit', str(tmp_path), '--print-order'])
    submit.main()
    assert capsys.readouterr().out.splitlines() == [f'{slots}\t{path}' for slots, path in
                                                    zip((20, 20, 10), expected)]
    calls = []
    monkeypatch.setattr(submit.subprocess, 'run', lambda args, check: calls.append(args))
    monkeypatch.setattr(sys, 'argv', ['submit', str(tmp_path)])
    submit.main()
    assert calls == [['qsub', path] for path in expected]

def test_submitter_refuses_script_slot_mismatch_before_any_qsub(tmp_path, monkeypatch):
    scripts, _ = _fleet(tmp_path)
    scripts[1].write_text('#!/bin/bash\n#$ -pe smp 10\n')
    calls = []
    monkeypatch.setattr(submit.subprocess, 'run', lambda args, check: calls.append(args))
    monkeypatch.setattr(sys, 'argv', ['submit', str(tmp_path)])
    with pytest.raises(SystemExit, match='slot request disagrees'):
        submit.main()
    assert calls == []
