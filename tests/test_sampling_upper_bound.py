"""sampling_upper_bound's Poisson widening must not fire when there is no sampling.

B3: scale == 1 means a --census pass with no read failures -- est_count is then an
exact count, not one draw from a random sample, so the "recover the raw sampled
count, inflate by z of it" rationale in the docstring does not apply. The shipped
2026-09 estimate had sample_scale 1.000000 and still widened a raw=4 count by +100%.
"""
from ligand_vdgs.generate_vdgs.estimate_frag_cost import sampling_upper_bound

def test_a_census_count_is_not_widened():
    assert sampling_upper_bound(4, 1.0) == 4

def test_a_real_sample_is_still_widened():
    # Discriminating pair: same est_count, scale < 1 (a genuine sample, so raw < 4
    # and the Poisson term is non-trivial). If this also returned 4 unchanged, the
    # case above would prove nothing about scale == 1 being special.
    assert sampling_upper_bound(4, 0.5) > 4

def test_zero_estimate_is_unaffected_by_the_census_case():
    assert sampling_upper_bound(0, 1.0) == 0
