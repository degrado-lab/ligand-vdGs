"""Severe scheduler/roster gates from review_findings_20260912.md (A1, B1).

Each of these silently did nothing before the fix: A1 let a default full build mine the
wrong parent database, B1 let the roster exit 0 after losing a large fraction of the census.
Discriminating pairs throughout, per DR-52/gates.md: a broken artifact must be refused, a
clean one must not be.
"""
import types
import pytest

from ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags import check_frags_dict_identity
from ligand_vdgs.generate_vdgs.build_ligand_roster import check_loss_threshold

# ---------------------------------------------------------------------------
# A1 -- default full build must refuse a --pdb-dir that isn't the dict's database.
# ---------------------------------------------------------------------------

def _identity(sha256, n_structures=10):
    return {'sha256': sha256, 'n_structures': n_structures, 'version': 1}

def test_mismatched_pdb_dir_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(
        'ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags.identity_of',
        lambda pdb_dir: _identity('b' * 64))
    with pytest.raises(SystemExit, match='different --pdb-dir'):
        check_frags_dict_identity(
            {'db_identity': _identity('a' * 64)},
            types.SimpleNamespace(pdb_dir=str(tmp_path), frags_dict='dict.pkl'))

def test_matching_pdb_dir_is_not_refused(tmp_path, monkeypatch):
    # Discriminating pair with the case above: identical inputs except the identity
    # agrees. If this raised too, the check would be a no-op that always fires.
    monkeypatch.setattr(
        'ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags.identity_of',
        lambda pdb_dir: _identity('a' * 64))
    check_frags_dict_identity(
        {'db_identity': _identity('a' * 64)},
        types.SimpleNamespace(pdb_dir=str(tmp_path), frags_dict='dict.pkl'))

def test_dict_with_no_recorded_identity_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(
        'ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags.identity_of',
        lambda pdb_dir: _identity('a' * 64))
    with pytest.raises(SystemExit, match='records no parent-database identity'):
        check_frags_dict_identity(
            {}, types.SimpleNamespace(pdb_dir=str(tmp_path), frags_dict='dict.pkl'))

# ---------------------------------------------------------------------------
# B1 -- the roster must refuse when it silently lost a large fraction of the census.
# ---------------------------------------------------------------------------

def test_high_loss_roster_is_refused():
    with pytest.raises(SystemExit, match='20 of 100'):
        check_loss_threshold(
            {'num_structures': 100, 'num_unreadable_files': 15,
             'num_errored_structures': 5}, '/fake/pdb_dir', '/fake/log')

def test_zero_loss_roster_is_not_refused():
    check_loss_threshold(
        {'num_structures': 100, 'num_unreadable_files': 0,
         'num_errored_structures': 0}, '/fake/pdb_dir', '/fake/log')

def test_loss_exactly_at_the_limit_is_not_refused():
    # Boundary case, strict '>' -- matches estimate_frag_cost.py's MAX_LOST_SAMPLE_FRACTION
    # policy exactly, since B1's whole point is that the two guards must not drift apart.
    check_loss_threshold(
        {'num_structures': 100, 'num_unreadable_files': 10,
         'num_errored_structures': 0}, '/fake/pdb_dir', '/fake/log')

def test_loss_just_over_the_limit_is_refused():
    with pytest.raises(SystemExit, match='11 of 100'):
        check_loss_threshold(
            {'num_structures': 100, 'num_unreadable_files': 11,
             'num_errored_structures': 0}, '/fake/pdb_dir', '/fake/log')

def test_dr51_build_stats_still_pass_the_new_gate():
    # DR-51's real roster log: num_unreadable_files never appears (defaultdict key never
    # incremented, i.e. exactly 0) and num_errored_structures is 0. The new gate must not
    # retroactively flag a build that was already clean.
    check_loss_threshold(
        {'num_structures': 65598, 'num_errored_structures': 0},
        '/fake/pdb_dir', '/fake/log')
