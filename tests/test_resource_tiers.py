"""Tier assignment in make_sge_scripts_for_frags, including the short-queue tier."""
import pytest

from ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags import (
    RESOURCE_TIERS, SHORT_QUEUE_TIER, _h_rt_to_hours, resources_for)

CEILING = '336:00:00'


def test_short_tier_fits_the_30_min_window():
    # The whole point of the tier: SGE admits a job to the short queue only if
    # its request is <= 30 min, so an off-by-a-minute h_rt silently loses it.
    assert _h_rt_to_hours(SHORT_QUEUE_TIER[2]) <= 0.5


@pytest.mark.parametrize('est,expected_h_rt', [
    (0, SHORT_QUEUE_TIER[2]),
    (SHORT_QUEUE_TIER[0] - 1, SHORT_QUEUE_TIER[2]),
    (SHORT_QUEUE_TIER[0], '6:00:00'),        # boundary is exclusive
    (999, '6:00:00'),
    (5000, CEILING),
])
def test_boundaries(est, expected_h_rt):
    assert resources_for(est, CEILING)[1] == expected_h_rt


def test_no_short_queue_falls_back_to_the_old_tier():
    assert resources_for(0, CEILING, short_queue=False) == (4, '6:00:00')


def test_max_h_rt_clamps_but_does_not_lengthen():
    # A ceiling below the tier clamps; a ceiling above it must not promote a
    # short job to a long h_rt and drop it out of the short queue.
    assert resources_for(10_000, '4:00:00')[1] == '4:00:00'
    assert resources_for(0, '40:00:00')[1] == SHORT_QUEUE_TIER[2]


def test_fixed_overrides_win():
    slots, h_rt = resources_for(0, CEILING, fixed_num_procs=2, fixed_h_rt='1:00:00')
    assert (slots, h_rt) == (2, '1:00:00')


def test_tiers_are_monotone_and_ordered():
    uppers = [t[0] for t in RESOURCE_TIERS]
    assert uppers == sorted(uppers)
    assert RESOURCE_TIERS[0] is SHORT_QUEUE_TIER
