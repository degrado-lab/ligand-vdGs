"""fragment_keys_nested: which fragment pairs share observations (Fig 1 exclusion set).

Falsifier: a nesting test that ignores annotations would call the P1 ester nested in
the P1 amide's 4-atom sub-key or vice versa, and one that ignores topology would miss
the 4-atom amide key inside the 5-atom one.
"""
import unittest

from ligand_vdgs.identify_bioisosteres.common import fragment_keys_nested
from tests.vacuity import assert_discriminates

AMIDE = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[N;!R;D2][C;!R;!H0]'
AMIDE_SUB = '[C;!R;!H0][C;!R;H0]([N;!R;D2])=[O;!R;D1]'
ESTER = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D2][C;!R;!H0]'
CARBOXYLATE = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D1]'
# Same atoms as CARBOXYLATE but the ester oxygen is D2: shares no observations with it.
ESTER_SUB = '[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D2]'

class TestFragmentKeysNested(unittest.TestCase):
    def test_nesting_needs_topology_and_identical_annotations(self):
        assert_discriminates(lambda pair: fragment_keys_nested(*pair),
                             accepts=[(AMIDE, AMIDE_SUB), (AMIDE_SUB, AMIDE), (ESTER, ESTER_SUB),
                                      (AMIDE, AMIDE)],
                             rejects=[(AMIDE, ESTER), (AMIDE_SUB, ESTER), (CARBOXYLATE, ESTER_SUB),
                                      (CARBOXYLATE, ESTER)],
                             label='fragment_keys_nested')

if __name__ == '__main__':
    unittest.main()
