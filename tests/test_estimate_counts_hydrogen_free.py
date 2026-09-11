"""_count_one must perceive ligands hydrogen-free, or every degree-annotated key counts 0.

Both toolkits count explicit H atoms in the SMARTS `D` primitive, so on a protonated
parent database -- which is what the pipeline mines -- a hydroxyl O reads `D2` and
`[O;D1;!R]` matches nothing. The fragment then falls below --min-instances and is never
built. The estimate has to route through functions/ligand_perception, not build its own
OBMol, so it cannot drift from what the miner actually matches.

The discriminating input is an *explicitly protonated* ligand; a file without hydrogens
passes whether or not the deletion happens.
"""
import os
import tempfile
import unittest

from ligand_vdgs.generate_vdgs import estimate_frag_cost as efc

# Methanol: C, O, and the hydroxyl H, with the H written as a real atom record.
# O is D1 only once that H is gone; with it, O is D2 and C is D2 rather than D1.
METHANOL = """\
HETATM    1  C1  MOH A 301       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O1  MOH A 301       1.430   0.000   0.000  1.00  0.00           O
HETATM    3  HO1 MOH A 301       1.760   0.900   0.000  1.00  0.00           H
END
"""


class HydrogenFreeCountingTests(unittest.TestCase):

    FRAGMENTS = ['[C;D1;!R][O;D1;!R]', '[C;D1;!R][O;D2;!R]']

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        # parent_db mirror layout, so read_ligand_blocks/stem_of find it.
        d = os.path.join(self.tmp, 'me')
        os.makedirs(d)
        self.path = os.path.join(d, '1meo.pdb')
        with open(self.path, 'w') as f:
            f.write(METHANOL)
        efc._init_worker(self.FRAGMENTS)

    def test_the_hydroxyl_counts_as_D1_not_D2(self):
        counts, read_failures, unreadable = efc._count_one(self.path)
        self.assertFalse(unreadable)
        self.assertEqual(read_failures, 0)
        # index 0 is the D1 key, index 1 the D2 key it would be mistaken for.
        self.assertEqual(counts.get(0), 1,
                         f'[C;D1][O;D1] should match methanol once; got {counts}')
        self.assertNotIn(1, counts,
                         f'[C;D1][O;D2] must not match a hydroxyl; got {counts}')

    def test_an_unreadable_block_is_a_read_failure_not_a_zero_count(self):
        """Counting an unparsable ligand as zero occurrences undercounts the fragment
        and can drop it below --min-instances; it is tallied separately instead."""
        d = os.path.join(self.tmp, 'ba')
        os.makedirs(d)
        bad = os.path.join(d, '1bad.pdb')
        with open(bad, 'w') as f:
            f.write('HETATM    1        XXX A 301    not numbers at all\n')
        counts, read_failures, unreadable = efc._count_one(bad)
        self.assertFalse(unreadable)
        self.assertEqual(counts, {})
        self.assertEqual(read_failures, 1)


if __name__ == '__main__':
    unittest.main()
