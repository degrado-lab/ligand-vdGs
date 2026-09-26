"""B6 -- estimate_frag_cost.py must pin frags_dict/pdb_dir identity BEFORE the
counting pass, not after. A pin taken after a long pass attests to whatever the
dict/DB looked like when the pass finished, not what the counts actually came
from -- the exact "grouping changed but every key still present" case the pin
exists to catch, and this project runs concurrent sessions that can change either
mid-run.
"""
import pickle
import sys

from ligand_vdgs.functions import ligand_structure
from ligand_vdgs.functions.db_identity import identity_of
from ligand_vdgs.functions.utils import file_sha256
from ligand_vdgs.generate_vdgs import estimate_frag_cost as efc
from tests.vacuity import assert_discriminates

FRAGMENT = '[C;!R;!H0][O;!R;D2]'

def _setup(tmp_path):
    for name in ('pdb',):
        (tmp_path / name).mkdir()
    dict_path = tmp_path / 'frags.pkl'
    with open(dict_path, 'wb') as handle:
        pickle.dump({'key_schema': ligand_structure.KEY_SCHEMA,
                     'frags': {'CO': {FRAGMENT: ['MET']}},
                     'support_pooled': {FRAGMENT: 10}}, handle)
    return dict_path, ['estimate_frag_cost.py',
                       '--frags-dict', str(dict_path),
                       '--pdb-dir', str(tmp_path / 'pdb'),
                       '--max-size', '5',
                       '--output', str(tmp_path / 'out.tsv')]

def test_a_dict_mutated_mid_pass_records_its_pre_pass_hash(
        monkeypatch, tmp_path):
    dict_path, argv = _setup(tmp_path)
    original_sha = file_sha256(str(dict_path))

    def _fake_counts(fragments, pdb_dir, sample_size=None, num_procs=1, seed=0):
        # Simulates a concurrent session rewriting the dict while this pass runs.
        with open(dict_path, 'ab') as handle:
            handle.write(b'\x00mutated-mid-pass')
        efc._LAST_SAMPLE_SCALE[0] = 1.0
        return ({f: 1 for f in fragments}, {f: 1 for f in fragments})

    monkeypatch.setattr(efc, 'estimate_fragment_counts', _fake_counts)
    monkeypatch.setattr(sys, 'argv', argv)
    efc.main()

    assert file_sha256(str(dict_path)) != original_sha, \
        'test setup broken: the mutation did not land'
    with open(argv[argv.index('--output') + 1]) as handle:
        headers = dict(line[2:].rstrip('\n').split('\t', 1) for line in handle
                       if line.startswith('# '))
    assert headers['frags_dict_sha256'] == original_sha

def test_a_db_grown_mid_pass_records_its_pre_pass_identity(monkeypatch, tmp_path):
    """DB side of B6. No DB_IDENTITY exists, so a pin taken after the pass would walk
    the grown DB and record the post-pass identity."""
    _dict_path, argv = _setup(tmp_path)
    pdb_dir = str(tmp_path / 'pdb')
    pre_pass = identity_of(pdb_dir)['sha256']

    def _fake_counts(fragments, pdb_dir, sample_size=None, num_procs=1, seed=0):
        # Simulates a mirror sync adding a structure while this pass runs.
        (tmp_path / 'pdb' / 'ab').mkdir()
        (tmp_path / 'pdb' / 'ab' / '1abc.pdb').write_text('ATOM      1  N   ALA A   1\n')
        efc._LAST_SAMPLE_SCALE[0] = 1.0
        return ({f: 1 for f in fragments}, {f: 1 for f in fragments})

    monkeypatch.setattr(efc, 'estimate_fragment_counts', _fake_counts)
    monkeypatch.setattr(sys, 'argv', argv)
    efc.main()

    post_pass = identity_of(pdb_dir)['sha256']
    with open(argv[argv.index('--output') + 1]) as handle:
        headers = dict(line[2:].rstrip('\n').split('\t', 1) for line in handle
                       if line.startswith('# '))
    assert_discriminates(lambda sha: sha == pre_pass, [headers['pdb_db_identity']],
                         [post_pass], 'header pins the pre-pass DB identity')
