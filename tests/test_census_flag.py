"""run_frag_cost_estimate.qsub must not hardcode --sample-size to today's exact structure
count, so the "census, not a sample" claim silently degraded into a sample the moment
the mirror grew past that literal. --census fixes this by passing sample_size=None
through to sample_pdb_paths, which takes every path regardless of the count.
"""
import os
import sys

from ligand_vdgs.generate_vdgs.estimate_frag_cost import sample_pdb_paths, parse_args

def _mirror(root, n):
    for i in range(n):
        d = os.path.join(root, str(i % 3))
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, f's{i}.pdb'), 'w').write('END\n')
    return root

def test_none_sample_size_takes_every_structure_regardless_of_count(tmp_path):
    # The mechanism --census relies on: unlike a hardcoded literal, None stays a
    # census no matter how many structures the mirror grows to hold.
    sample, total = sample_pdb_paths(str(_mirror(tmp_path / 'grown', 37)), None)
    assert total == 37
    assert len(sample) == 37

def test_a_sample_size_smaller_than_the_mirror_is_a_strict_subset(tmp_path):
    # Discriminating pair: same mirror, a real cap. If this returned everything
    # too, the case above would prove nothing about None being special.
    sample, total = sample_pdb_paths(str(_mirror(tmp_path / 'grown', 37)), 10)
    assert total == 37
    assert len(sample) == 10

def test_census_flag_defaults_to_false(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['estimate_frag_cost.py',
                                      '--pdb-dir', '/fake', '--sample-size', '5000'])
    assert parse_args().census is False

def test_census_flag_is_wired_through_argparse(monkeypatch):
    monkeypatch.setattr(sys, 'argv',
                        ['estimate_frag_cost.py', '--pdb-dir', '/fake', '--census'])
    assert parse_args().census is True
