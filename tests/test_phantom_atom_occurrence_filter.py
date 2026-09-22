"""B2 -- estimate_frag_cost.py's occurrence count must not count a match that
reaches a template-added phantom (unobserved heavy atom). The miner drops exactly
these matches (cg.py:355), so counting them here inflates the occurrence column
that keys the resource tiers (DR-7).

Fixture: acetate (ACT) with its methyl carbon unmodelled -- the same partial-density
fixture test_ccd_templates.py uses to establish that the phantom is a real, bonded
graph atom with no PDB name.
"""
from ligand_vdgs.generate_vdgs import estimate_frag_cost as efc

ACETATE = """\
HETATM    1  C   ACT A 900       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O   ACT A 900       1.250   0.000   0.000  1.00  0.00           O
HETATM    3  OXT ACT A 900      -0.700   1.210   0.000  1.00  0.00           O
HETATM    4  CH3 ACT A 900      -0.750  -1.290   0.000  1.00  0.00           C
END
"""
# Methyl carbon unmodelled -- routine partial density, same as test_ccd_templates.py.
ACETATE_PARTIAL = "\n".join(
    line for line in ACETATE.splitlines() if ' CH3 ' not in line) + "\n"

# Reaches the methyl carbon: real occurrence when CH3 is observed, phantom-only
# when it isn't.
THROUGH_METHYL = '[C][C](=[O])[O]'
# Never touches the methyl carbon at all.
CARBOXYLATE_ONLY = '[C](=[O])[O]'
FRAGMENTS = [THROUGH_METHYL, CARBOXYLATE_ONLY]

def _write(tmp_path, text):
    d = tmp_path / 'ac'
    d.mkdir()
    path = d / '1act.pdb'
    path.write_text(text)
    return str(path)

def test_a_match_reaching_the_phantom_is_not_counted(tmp_path):
    efc._init_worker(FRAGMENTS)
    counts, read_failures, unreadable, all_failed = efc._count_one(_write(tmp_path, ACETATE_PARTIAL))
    assert not unreadable and read_failures == 0 and not all_failed
    assert FRAGMENTS.index(THROUGH_METHYL) not in counts, counts

def test_a_match_away_from_the_phantom_is_still_counted(tmp_path):
    # Discriminating pair: same structure, a fragment that never reaches CH3. If
    # this were also absent, the case above would prove nothing about phantoms
    # specifically -- it would mean nothing counts at all.
    efc._init_worker(FRAGMENTS)
    counts, read_failures, unreadable, all_failed = efc._count_one(_write(tmp_path, ACETATE_PARTIAL))
    assert not unreadable and read_failures == 0 and not all_failed
    assert counts.get(FRAGMENTS.index(CARBOXYLATE_ONLY)) == 1, counts

def test_the_same_match_is_counted_once_fully_observed(tmp_path):
    # Proves THROUGH_METHYL is a real, matchable fragment -- its absence above is
    # the phantom filter working, not a SMARTS that can never match this mol.
    efc._init_worker(FRAGMENTS)
    counts, read_failures, unreadable, all_failed = efc._count_one(_write(tmp_path, ACETATE))
    assert not unreadable and read_failures == 0 and not all_failed
    assert counts.get(FRAGMENTS.index(THROUGH_METHYL)) == 1, counts
