"""Severe scheduler/roster gates from review_findings_20260912.md (A1, B1).

Each of these silently did nothing before the fix: A1 let a default full build mine the
wrong parent database, B1 let the roster exit 0 after losing a large fraction of the census.
Discriminating pairs throughout, per docs/directives/testing.md: a broken artifact
must be refused, a clean one must not be.
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

def _args(tmp_path, sha256):
    """main() pins the run's DB identity on args once, before any check."""
    return types.SimpleNamespace(pdb_dir=str(tmp_path), frags_dict='dict.pkl', db_identity=sha256)

def test_mismatched_pdb_dir_is_refused(tmp_path):
    with pytest.raises(SystemExit, match='different --pdb-dir'):
        check_frags_dict_identity({'db_identity': _identity('a' * 64)}, _args(tmp_path, 'b' * 64))

def test_matching_pdb_dir_is_not_refused(tmp_path):
    # Discriminating pair with the case above: identical inputs except the identity
    # agrees. If this raised too, the check would be a no-op that always fires.
    check_frags_dict_identity({'db_identity': _identity('a' * 64)}, _args(tmp_path, 'a' * 64))

def test_dict_with_no_recorded_identity_is_refused(tmp_path):
    with pytest.raises(SystemExit, match='records no parent-database identity'):
        check_frags_dict_identity({}, _args(tmp_path, 'a' * 64))

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
    # A real production roster log: num_unreadable_files never appears (defaultdict key never
    # incremented, i.e. exactly 0) and num_errored_structures is 0. The new gate must not
    # retroactively flag a build that was already clean.
    check_loss_threshold(
        {'num_structures': 65598, 'num_errored_structures': 0},
        '/fake/pdb_dir', '/fake/log')

# Instance-level loss: every structure reads, but perception drops ligands inside them.
def _instance_stats(kept, fallback, unreadable=0):
    return {'num_structures': 100, 'num_errored_structures': 0, 'num_instances': kept,
            'num_ob_fallback_instances': fallback, 'num_unreadable_instances': unreadable}

def test_mass_ob_fallback_is_refused_even_when_all_files_read():
    # Falsifier: a perception regression pushing 40% of instances to OpenBabel passes the
    # structure-level gate (0 files lost); only the instance gate can catch it.
    with pytest.raises(SystemExit, match='40 of 100 ligand instances'):
        check_loss_threshold(_instance_stats(60, 40), '/fake/pdb_dir', '/fake/log')

def test_unreadable_instances_count_toward_the_instance_gate():
    with pytest.raises(SystemExit, match='11 of 100 ligand instances'):
        check_loss_threshold(_instance_stats(89, 5, 6), '/fake/pdb_dir', '/fake/log')

def test_instance_loss_exactly_at_the_limit_is_not_refused():
    check_loss_threshold(_instance_stats(90, 5, 5), '/fake/pdb_dir', '/fake/log')

def test_dr51_instance_stats_pass_the_instance_gate():
    # A real production roster log: 1031 fallback of 143944 instances (0.7%), 0 unreadable.
    stats = _instance_stats(142913, 1031)
    stats['num_structures'] = 65598
    check_loss_threshold(stats, '/fake/pdb_dir', '/fake/log')
