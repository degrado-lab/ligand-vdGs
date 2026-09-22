"""assert_discriminates must REFUSE the two clauses confirmed vacuous -- the bb_
clause (DR-49 addendum 1) and the aromatic clause (vacuity audit V7) -- and accept
their replacements. If this file goes green with the broken clauses in place, the
helper has stopped discriminating and every call to it is worthless.
"""
import pytest
from rdkit import Chem
from vacuity import assert_discriminates

# As shipped, and the part-wise replacement now at test_bb_slot_attribution.py:72.
def _bb_broken(key): return 'bb_' not in key.replace('-bb', '')
def _bb_fixed(key): return 'bb' not in key.split('_')

# V7 'is any key aromatic': the shipped form fires on the 'l' of Cl and on the
# lowercase ring primitive r, so it is True for entirely aliphatic key sets.
def _arom_broken(keys): return any(any(c.islower() for c in k if c.isalpha()) for k in keys)
def _arom_fixed(keys):
    return any(a.GetIsAromatic() for k in keys for a in Chem.MolFromSmarts(k).GetAtoms())

DECONVOLUTED = 'ALA-bb_TRP-bb'
UNDECONVOLUTED = 'ALA-bb_bb'        # what a first-slot-only reader emits
AROMATIC_KEYS = ('[c;D3]-[n;D2]',)
ALIPHATIC_KEYS = ('[C;D4]-[Cl;D1]', '[C;r6;D3]')

def test_helper_refuses_the_bb_clause_that_shipped():
    with pytest.raises(AssertionError, match='cannot fire'):
        assert_discriminates(_bb_broken, [DECONVOLUTED], [UNDECONVOLUTED], 'bb')

def test_helper_accepts_the_part_wise_replacement():
    assert_discriminates(_bb_fixed, [DECONVOLUTED], [UNDECONVOLUTED], 'bb')

def test_helper_refuses_the_aromatic_clause_that_shipped():
    with pytest.raises(AssertionError, match='cannot fire'):
        assert_discriminates(_arom_broken, [AROMATIC_KEYS], [ALIPHATIC_KEYS], 'aromatic')

def test_helper_accepts_the_parsed_query_replacement():
    assert_discriminates(_arom_fixed, [AROMATIC_KEYS], [ALIPHATIC_KEYS], 'aromatic')

def test_one_sided_use_is_itself_refused():
    """The likeliest misuse: handing the helper only inputs the clause accepts."""
    with pytest.raises(AssertionError, match='one-sided'):
        assert_discriminates(_bb_fixed, [DECONVOLUTED], [], 'bb')

def test_the_broken_and_fixed_clauses_genuinely_differ():
    """Non-vacuity for this file: if the two clauses agreed on the counterexamples,
    the four tests above would pass without saying anything about the helper.
    """
    assert _bb_broken(UNDECONVOLUTED) and not _bb_fixed(UNDECONVOLUTED)
    assert _arom_broken(ALIPHATIC_KEYS) and not _arom_fixed(ALIPHATIC_KEYS)
    assert _bb_fixed(DECONVOLUTED) and _arom_fixed(AROMATIC_KEYS)
