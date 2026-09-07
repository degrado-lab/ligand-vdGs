"""Preprocessing filters that keep prepwizard from mangling parent PDBs.

Prepwizard cannot build residues it does not recognize as amino acids. It reclassifies
them, and re-emits the result under some *other* residue's identity, so the damage is
invisible downstream -- the output is still a well-formed PDB. Three behaviours are
covered here, all observed in a real 65,600-PDB database.
"""
import os
import shutil
import tempfile
import unittest

import prody as pr

from ligand_vdgs.preprocessing._prep_filters import (
    drop_prepwizard_hazard_residues, incomplete_backbone_resindices,
    restore_renamed_ligands, single_carbon_het_resindices,
    snapshot_restorable_ligands)


def _atom(serial, name, resname, chain, resnum, xyz, element, het=False):
    return (f"{'HETATM' if het else 'ATOM  '}{serial:5d} {name:<4s}{resname:>4s} {chain}"
            f"{resnum:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
            f"{element:>12s}\n")


def _residue(serial, resname, chain, resnum, origin, names_elements, het=False):
    """A residue whose atoms are laid out 1.5 A apart along x from *origin*."""
    lines, x0, y0, z0 = [], *origin
    for i, (name, element) in enumerate(names_elements):
        lines.append(_atom(serial + i, name, resname, chain, resnum,
                           (x0 + 1.5 * i, y0, z0), element, het))
    return lines


# A minimal but realistic structure: two complete alanines, one truncated to a lone
# backbone N, a free ligand well away from the protein, and a one-carbon HET group.
GLY_BB = [('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O')]
LIGAND = [('N1', 'N'), ('C2', 'C'), ('O2', 'O'), ('N3', 'N')]


def _write_structure(path):
    lines = []
    lines += _residue(1, 'ALA', 'A', 1, (0.0, 0.0, 0.0), GLY_BB)
    lines += _residue(10, 'ALA', 'A', 2, (0.0, 4.0, 0.0), GLY_BB)
    lines += _residue(20, 'ALA', 'A', 3, (0.0, 8.0, 0.0), [('N', 'N')])  # truncated
    lines += _residue(30, 'CYT', 'B', 201, (40.0, 0.0, 0.0), LIGAND, het=True)
    lines += _residue(40, 'CF0', 'B', 202, (60.0, 0.0, 0.0), [('C1', 'C')], het=True)
    lines += _residue(50, 'ZN', 'B', 203, (80.0, 0.0, 0.0), [('ZN', 'ZN')], het=True)
    with open(path, 'w') as fh:
        fh.writelines(lines)
        fh.write('END\n')


class PrepFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.pdb = os.path.join(self.tmp, 'test.pdb')
        _write_structure(self.pdb)
        self.ag = pr.parsePDB(self.pdb)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _resnames(self, resindices):
        hv = self.ag.getHierView()
        return sorted(f'{r.getResname()} {r.getChid()}{r.getResnum()}'
                      for r in hv.iterResidues() if r.getResindices()[0] in resindices)

    def test_lone_backbone_n_is_dropped(self):
        # The case behind 2y1x SAH A:1001 -- prepwizard turns this N into an ammonium
        # ion carrying a ligand's resname, 54 A from the ligand it is now labelled as.
        self.assertEqual(self._resnames(incomplete_backbone_resindices(self.ag)),
                         ['ALA A3'])

    def test_complete_residues_and_ions_survive(self):
        kept = self._resnames(drop_prepwizard_hazard_residues(self.ag))
        self.assertIn('ALA A1', kept)
        self.assertIn('ALA A2', kept)
        self.assertIn('CYT B201', kept)
        self.assertIn('ZN B203', kept)

    def test_single_carbon_het_is_dropped_but_monatomic_ion_is_not(self):
        # CF0/0QE are covalent-inhibitor warheads; ZN is a real ion and must survive.
        # The rule is restricted to carbon precisely to keep that distinction.
        dropped = self._resnames(single_carbon_het_resindices(self.ag))
        self.assertEqual(dropped, ['CF0 B202'])

    def test_renamed_free_ligand_is_restored(self):
        # Prepwizard rewrites CYT (free cytosine) into an ATOM record named CYS, which
        # hides it from find_cg_matches -- that reader only looks at HETATM lines.
        snapshot = snapshot_restorable_ligands(self.pdb)
        prepped = os.path.join(self.tmp, 'prepped.pdb')
        with open(self.pdb) as fh, open(prepped, 'w') as out:
            for line in fh:
                if line[17:20].strip() == 'CYT':
                    line = 'ATOM  ' + line[6:17] + 'CYS' + line[20:]
                out.write(line)

        self.assertEqual(restore_renamed_ligands(prepped, snapshot), 1)
        restored = [l for l in open(prepped) if l[22:26].strip() == '201']
        self.assertTrue(all(l.startswith('HETATM') for l in restored))
        self.assertTrue(all(l[17:20] == 'CYT' for l in restored))

    def test_chain_bonded_residue_is_not_restorable(self):
        # A chromophore or modified residue bonded into the chain (CR8 in 3tmr, MDO in
        # 2o7d) is protein. Turning it back into HETATM would make it mineable as a
        # ligand -- a change in library composition, not a bug fix.
        bonded = os.path.join(self.tmp, 'bonded.pdb')
        lines = open(self.pdb).readlines()
        # Place the HET residue in bonding contact with ALA A1's C.
        lines = [l for l in lines if l[22:26].strip() != '201']
        lines += _residue(30, 'CR8', 'B', 201, (4.4, 0.0, 0.0), LIGAND, het=True)
        with open(bonded, 'w') as fh:
            fh.writelines(lines)
        self.assertNotIn('CR8', set(snapshot_restorable_ligands(bonded).values()))

    def test_restore_is_a_no_op_when_nothing_was_renamed(self):
        snapshot = snapshot_restorable_ligands(self.pdb)
        untouched = os.path.join(self.tmp, 'untouched.pdb')
        shutil.copy(self.pdb, untouched)
        before = open(untouched).read()
        self.assertEqual(restore_renamed_ligands(untouched, snapshot), 0)
        self.assertEqual(open(untouched).read(), before)


if __name__ == '__main__':
    unittest.main()
