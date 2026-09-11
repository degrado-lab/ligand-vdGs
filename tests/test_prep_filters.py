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


class ModifiedResidueTests(unittest.TestCase):
    """Rebuild without prepwizard: modified residues must become protein records, not
    ligands (see _prep_filters.modified_residues_to_protein)."""

    CCD = {'KCX': 'L-PEPTIDE LINKING', 'DCY': 'D-PEPTIDE LINKING', 'PLP': 'NON-POLYMER',
           'SAM': 'NON-POLYMER', '04C': 'peptide-like', 'HOH': 'NON-POLYMER'}

    @staticmethod
    def _res(serial, resname, chain, resnum, x0, names, het, extra=()):
        """N, CA, C along x at 1.3 A so residue i's C is 1.3 A from residue i+1's N when
        consecutive residues start 3.9 A apart; *extra* = [(name, element, dx, dy)]."""
        lines = [_atom(serial + i, n, resname, chain, resnum, (x0 + 1.3 * i, 0.0, 0.0),
                       n[0], het) for i, n in enumerate(names)]
        for j, (n, el, dx, dy) in enumerate(extra):
            lines.append(_atom(serial + len(names) + j, n, resname, chain, resnum,
                               (x0 + dx, dy, 0.0), el, het))
        return lines

    def _structure(self):
        from ligand_vdgs.preprocessing._prep_filters import modified_residues_to_protein
        L = []
        bb = ['N', 'CA', 'C']
        # chain A: ALA - MSE(het) - KCX(het, carbamate on NZ) - ALA - SEC(het) - DCY(het)
        L += self._res(1, 'ALA', 'A', 1, 0.0, bb, False)
        L += self._res(10, 'MSE', 'A', 2, 3.9, bb, True, [('SE', 'SE', 1.3, 2.0)])
        L += self._res(20, 'KCX', 'A', 3, 7.8, bb, True, [('NZ', 'N', 1.3, 4.0),
                                                         ('CX', 'C', 1.3, 5.3)])
        L += self._res(30, 'ALA', 'A', 4, 11.7, bb, False)
        L += self._res(40, 'SEC', 'A', 5, 15.6, bb, True, [('SE', 'SE', 1.3, 2.0)])
        L += self._res(50, 'DCY', 'A', 6, 19.5, bb, True, [('SG', 'S', 1.3, 2.0)])
        # PLP covalently on KCX-adjacent lysine-like N: 1.5 A from KCX NZ, NON-POLYMER
        L += [_atom(60, 'C4A', 'PLP', 'A', 301, (9.1, 5.5, 0.0), 'C', True),
              _atom(61, 'N1', 'PLP', 'A', 301, (9.1, 6.8, 0.0), 'N', True)]
        # SAM ligand with N/CA/C names, 1.6 A from ALA 4's C: cofactor, stays HETATM
        L += self._res(70, 'SAM', 'A', 302, 14.3 + 1.6, ['N', 'CA', 'C'], True)
        # free MSE ligand with OXT far away: selenomethionine ligand, stays HETATM
        L += self._res(80, 'MSE', 'B', 1, 100.0, bb + ['OXT'], True,
                       [('SE', 'SE', 1.3, 2.0)])
        # free KCX ligand (not bonded to anything standard): stays HETATM
        L += self._res(90, 'KCX', 'C', 1, 200.0, bb + ['OXT'], True)
        # unknown resname bonded to ALA 1's N: left alone, reported
        L += [_atom(99, 'C1', 'ZZZ', 'A', 401, (-1.5, 0.0, 0.0), 'C', True)]
        # peptide-like ligand bonded to protein (covalent inhibitor): stays HETATM
        L += [_atom(100, 'C1', '04C', 'A', 402, (1.3, -1.5, 0.0), 'C', True)]
        L += [_atom(101, 'O', 'HOH', 'W', 1, (50.0, 50.0, 0.0), 'O', True)]
        return L, modified_residues_to_protein(L, self.CCD)

    def test_records_and_resnames(self):
        _, (new, stats) = self._structure()
        rec = {}
        for l in new:
            if l[:6] in ('ATOM  ', 'HETATM'):
                rec[(l[21], int(l[22:26]))] = (l[:6].strip(), l[17:20], l[12:16].strip())
        names = {k: set() for k in rec}
        for l in new:
            if l[:6] in ('ATOM  ', 'HETATM'):
                names[(l[21], int(l[22:26]))].add(l[12:16].strip())
        self.assertEqual(rec[('A', 2)][:2], ('ATOM', 'MET'))      # MSE bonded
        self.assertIn('SE', names[('A', 2)])                       # SE kept, not SD
        self.assertEqual(rec[('A', 3)][:2], ('ATOM', 'KCX'))      # own resname kept
        self.assertEqual(rec[('A', 5)][:2], ('ATOM', 'CYS'))      # SEC bonded
        self.assertEqual(rec[('A', 6)][:2], ('ATOM', 'DCY'))      # D-peptide linking
        self.assertEqual(rec[('A', 301)][:2], ('HETATM', 'PLP'))  # cofactor
        self.assertEqual(rec[('A', 302)][:2], ('HETATM', 'SAM'))  # N/CA/C names, non-polymer
        self.assertEqual(rec[('B', 1)][:2], ('HETATM', 'MET'))    # free Se-Met ligand
        self.assertEqual(rec[('C', 1)][:2], ('HETATM', 'KCX'))    # free ncAA ligand
        self.assertEqual(rec[('A', 401)][:2], ('HETATM', 'ZZZ'))  # unknown
        self.assertEqual(rec[('A', 402)][:2], ('HETATM', '04C'))  # peptide-like
        self.assertEqual(rec[('W', 1)][:2], ('HETATM', 'HOH'))
        self.assertEqual(stats, {'renamed': 3, 'amino_acids_to_atom': 2,
                                 'modres_to_atom': 2, 'unknown_resnames': ['ZZZ']})

    def test_only_record_and_resname_columns_change(self):
        old, (new, _) = self._structure()
        self.assertEqual(len(old), len(new))
        for a, b in zip(old, new):
            self.assertEqual(a[20:], b[20:])
            self.assertEqual(a[6:17], b[6:17])

    def test_missing_table_means_no_conversion_but_renames_still_happen(self):
        from ligand_vdgs.preprocessing._prep_filters import (
            load_ccd_polymer_types, modified_residues_to_protein)
        self.assertEqual(load_ccd_polymer_types('/nonexistent/ccd.tsv'), {})
        old, _ = self._structure()
        new, stats = modified_residues_to_protein(old, {})
        self.assertEqual(stats['renamed'], 3)
        self.assertEqual(stats['modres_to_atom'], 0)
        self.assertIn('KCX', stats['unknown_resnames'])
        kcx = [l for l in new if l[17:20] == 'KCX' and l[21] == 'A']
        self.assertTrue(all(l.startswith('HETATM') for l in kcx))

    def test_nothing_to_do_returns_equal_lines(self):
        from ligand_vdgs.preprocessing._prep_filters import modified_residues_to_protein
        lines = self._res(1, 'ALA', 'A', 1, 0.0, ['N', 'CA', 'C'], False)
        new, stats = modified_residues_to_protein(lines, self.CCD)
        self.assertEqual(new, lines)
        self.assertEqual(stats, {'renamed': 0, 'amino_acids_to_atom': 0,
                                 'modres_to_atom': 0, 'unknown_resnames': []})
