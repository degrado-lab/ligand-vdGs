"""The roster records atom names the fragment side can actually join against.

Three ways this goes wrong, all silent -- a fragment simply draws less support,
with no error anywhere:

1. CCD *alt* names. `ccd_templates.index_by_name` accepts them (older depositions
   write `O1P` for `OP1`, `*` for `'` in nucleotides), so an instance is templated
   under a spelling the fragment side never emits -- enumeration walks
   `template.atoms` and writes `atom.name`. Recorded raw, the two vocabularies are
   disjoint for that instance.
2. Placed hydrogens. Their names come from the protonation tool, not the CCD
   (`HCHA` on a HEM), so looking them up in the template before filtering rejects
   the whole instance. Measured: 5,453 of 5,906 instances on a 3,000-structure
   roster.
3. OpenBabel atom IDs. `apply_template_to_obmol` regenerates the residue's atom
   IDs from element symbols as a side effect of its `AddBond` calls, so
   `pdb_atom_names` on a TEMPLATED instance can return `FE, FE, S, S` for a block
   that deposited `FE1, FE2, S1, S2`. That is why the names are read from the PDB
   block rather than from the perceived OBMol. Covered by
   `test_ob_atom_ids_are_destroyed_by_templating`, which fails if OpenBabel ever
   stops doing this -- at which point the block read is merely redundant, not
   wrong.

Unobserved atoms need no special handling and get a test anyway: they are simply
absent from the block, which is exactly the "every atom observed" clause.
"""
import unittest

from ligand_vdgs.functions import ccd_templates, ligand_perception
from ligand_vdgs.generate_vdgs.build_ligand_roster import _canonical_names


def _hetatm(serial, name, resname, element, xyz=(0.0, 0.0, 0.0)):
    """One HETATM line with `name` in columns 13-16 and `element` in 77-78."""
    return (f'HETATM{serial:>5d} {name:<4s} {resname:>3s} A 501    '
            f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          '
            f'{element:>2s}')


class CanonicalNameTests(unittest.TestCase):

    def test_alt_names_are_mapped_to_the_templates_own_spelling(self):
        """`O5*` must be recorded as `O5'`, the name enumeration emits."""
        block = '\n'.join([_hetatm(1, "O5*", 'ATP', 'O'),
                           _hetatm(2, "C5*", 'ATP', 'C'),
                           _hetatm(3, "N1", 'ATP', 'N')])
        names = _canonical_names('ATP', block)
        self.assertEqual(set(names), {"O5'", "C5'", 'N1'})
        # Discriminating: the raw spellings must NOT survive, or the join is still
        # against a vocabulary the fragment side never writes.
        self.assertNotIn('O5*', names)
        self.assertNotIn('C5*', names)

    def test_canonical_names_are_left_alone(self):
        block = '\n'.join([_hetatm(1, "O5'", 'ATP', 'O'), _hetatm(2, 'N1', 'ATP', 'N')])
        self.assertEqual(set(_canonical_names('ATP', block)), {"O5'", 'N1'})

    def test_placed_hydrogens_are_skipped_not_rejected(self):
        """A tool-named H must not take the whole instance down with it."""
        block = '\n'.join([_hetatm(1, 'N1', 'ATP', 'N'),
                           _hetatm(2, 'HZZ9', 'ATP', 'H'),
                           _hetatm(3, 'C2', 'ATP', 'C')])
        names = _canonical_names('ATP', block)
        self.assertIsNotNone(names)
        self.assertEqual(set(names), {'N1', 'C2'})

    def test_a_template_hydrogen_is_dropped_even_with_a_blank_element_column(self):
        """The template's element decides, so a blank element column still works."""
        template = ccd_templates.get_template('ATP')
        h_name = next(a.name for a in template.atoms if a.element in ('H', 'D'))
        line = _hetatm(2, h_name, 'ATP', 'H')
        blanked = line[:76] + '  '
        block = '\n'.join([_hetatm(1, 'N1', 'ATP', 'N'), blanked])
        self.assertEqual(set(_canonical_names('ATP', block)), {'N1'})

    def test_unobserved_heavy_atoms_are_simply_absent(self):
        block = _hetatm(1, 'N1', 'ATP', 'N')
        names = _canonical_names('ATP', block)
        self.assertEqual(set(names), {'N1'})
        self.assertNotIn('N3', names)

    def test_a_heavy_name_unknown_to_the_template_fails_loudly(self):
        block = '\n'.join([_hetatm(1, 'N1', 'ATP', 'N'),
                           _hetatm(2, 'ZZ9', 'ATP', 'C')])
        self.assertIsNone(_canonical_names('ATP', block))

    def test_non_hetatm_lines_are_ignored(self):
        """read_ligand_blocks appends CONECT records to the block."""
        block = '\n'.join([_hetatm(1, 'N1', 'ATP', 'N'),
                           'CONECT    1    2',
                           'ANISOU    1  N1  ATP A 501'])
        self.assertEqual(set(_canonical_names('ATP', block)), {'N1'})


class OpenBabelAtomIdTests(unittest.TestCase):
    """Why the names come from the block and not from `pdb_atom_names`."""

    FES_BLOCK = (
        'HETATM 9600  FE1 FES E 501       8.958  81.453  73.804  1.00 40.05          Fe\n'
        'HETATM 9601  FE2 FES E 501      10.333  81.857  76.137  1.00 40.33          Fe\n'
        'HETATM 9602  S1  FES E 501       8.252  81.008  75.864  1.00 39.83           S\n'
        'HETATM 9603  S2  FES E 501      11.097  82.079  74.050  1.00 38.58           S\n')

    def test_ob_atom_ids_are_destroyed_by_templating(self):
        """Templating rewrites `FE1/FE2/S1/S2` to `FE/FE/S/S`.

        If this ever starts passing the wrong way -- i.e. the IDs survive -- the
        block read in `_canonical_names` becomes redundant rather than wrong, and
        this test says so instead of silently protecting nothing.
        """
        perceived = ligand_perception.perceive_ligand_instance(self.FES_BLOCK, 'FES')
        self.assertIsNotNone(perceived)
        self.assertEqual(perceived.provenance,
                         ligand_perception.PERCEPTION_CCD_TEMPLATE)
        from_obmol = set(ligand_perception.pdb_atom_names(perceived.obmol).values())
        from_block = set(_canonical_names('FES', self.FES_BLOCK))
        self.assertEqual(from_block, {'FE1', 'FE2', 'S1', 'S2'})
        self.assertNotEqual(from_obmol, from_block)
        self.assertTrue(from_obmol.issubset({'FE', 'S'}),
                        f'unexpected OB atom IDs: {from_obmol}')


if __name__ == '__main__':
    unittest.main()
