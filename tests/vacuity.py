"""DR-52's rule -- a check is not trusted until its refusal has been observed --
applied to a clause INSIDE a test. A clause labelled non-vacuity discharges the
reviewer's obligation to look, so one that cannot fire is worse than none
(DR-49 addendum 1). Refusal demonstrated in test_vacuity_helper.py.
"""

def assert_discriminates(clause, accepts, rejects, label=''):
    """Fail unless `clause` accepts every `accepts` and rejects every `rejects`.
    `rejects` are the counterexamples the clause claims to exclude; naming them is
    what stops the clause from being a string heuristic nobody ever exercised.
    Reports WHICH input misbehaved, never a count (see the testing directive).
    """
    accepts, rejects = list(accepts), list(rejects)
    if not accepts or not rejects:
        raise AssertionError(
            f'{label}: assert_discriminates needs at least one accepted AND one '
            f'rejected input; a one-sided call is the vacuity it exists to catch')
    for case in accepts:
        assert clause(case), f'{label}: clause rejected {case!r}, which it must accept'
    for case in rejects:
        assert not clause(case), (
            f'{label}: clause accepted {case!r}, the state it claims to reject -- '
            f'it cannot fire, so the clause is itself vacuous')
