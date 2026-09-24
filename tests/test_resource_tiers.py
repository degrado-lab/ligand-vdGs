"""Tier assignment in make_sge_scripts_for_frags, after the 2026-09-12 collapse to
two bands. NAMED FALSIFIERS, one per invariant:
  * a band above 48 h -> test_no_band_crosses_the_latency_cliff. Past 48 h the
    median queue wait goes 1.6 min -> 21.1 min (p90 8.6 h).
  * tiers keyed on the wrong estimate column -> test_boundaries.
  * a body that slices the table (a retired path dropped TIER_1 via
    `RESOURCE_TIERS[1:]`) -> test_boundaries. Injected after collection so
    expectations still came from the full table: [0-10] and [19999-10] failed
    `assert 30 == 10` (DR-50 addendum, DR-57).
"""
# Why the assertions key on SLOTS and not h_rt: the old table had four bands with
# four h_rt values, so boundary cases could assert on h_rt. Both surviving bands
# request 48:00:00, so an h_rt-keyed boundary test would now pass against ANY
# implementation, including one ignoring the estimate entirely (DR-20, DR-24,
# DR-27). test_boundary_cases_are_not_vacuous asserts the two bands really differ.
import pytest

from ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags import (
    RESOURCE_TIERS, TIER_1, _h_rt_to_hours, resources_for)

LATENCY_CLIFF_H = 48.0       # measured: requests past this wait 13x longer

def _boundary_cases():
    """(est_occurrences, expected slots) either side of every tier boundary.

    Derived from RESOURCE_TIERS, not written out: hard-coded counts went stale
    silently once DR-7 rekeyed the tiers from structures onto occurrences.
    """
    cases = [(0, RESOURCE_TIERS[0][1])]
    for i, (upper, slots, _h_rt) in enumerate(RESOURCE_TIERS[:-1]):
        cases.append((upper - 1, slots))                     # last count inside
        cases.append((upper, RESOURCE_TIERS[i + 1][1]))      # boundary is exclusive
    return cases

@pytest.mark.parametrize('est,expected_slots', _boundary_cases())
def test_boundaries(est, expected_slots):
    assert resources_for(est)[0] == expected_slots

def test_boundary_cases_are_not_vacuous():
    """The clause that fails if the parametrised test proves nothing.

    If every band carried the same slot count, test_boundaries would pass without
    distinguishing any input from any other. Require the bands to straddle a real
    change. This is the assertion that h_rt can no longer support.
    """
    cases = _boundary_cases()
    assert len({s for _, s in cases}) >= 2, cases
    for (_, below), (_, above) in zip(cases[1::2], cases[2::2]):
        assert below != above, (below, above)
    # And confirm h_rt is genuinely NOT discriminating, so nobody reinstates it
    # as the boundary axis without noticing it has gone flat.
    assert len({t[2] for t in RESOURCE_TIERS}) == 1, RESOURCE_TIERS

def test_no_band_crosses_the_latency_cliff():
    for upper, slots, h_rt in RESOURCE_TIERS:
        assert _h_rt_to_hours(h_rt) <= LATENCY_CLIFF_H, (upper, slots, h_rt)

def test_fixed_overrides_win():
    slots, h_rt = resources_for(0, fixed_num_procs=2, fixed_h_rt='1:00:00')
    assert (slots, h_rt) == (2, '1:00:00')

def test_tiers_are_monotone_and_ordered():
    uppers = [t[0] for t in RESOURCE_TIERS]
    assert uppers == sorted(uppers)
    assert RESOURCE_TIERS[0] is TIER_1
    assert RESOURCE_TIERS[-1][0] == float('inf')
