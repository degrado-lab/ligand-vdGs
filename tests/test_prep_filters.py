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
    lines, x0, y0, z0 = [], *origin
    for i, (name, element) in enumerate(names_elements):
        lines.append(_atom(serial + i, name, resname, chain, resnum,
                           (x0 + 1.5 * i, y0, z0), element, het))
    return lines

GLY_BB = [('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O')]
LIGAND = [('N1', 'N'), ('C2', 'C'), ('O2', 'O'), ('N3', 'N')]

def _write_structure(path):
    lines = []
    lines += _residue(1, 'ALA', 'A', 1, (0.0, 0.0, 0.0), GLY_BB)
    lines += _residue(10, 'ALA', 'A', 2, (0.0, 4.0, 0.0), GLY_BB)
    lines += _residue(20, 'ALA', 'A', 3, (0.0, 8.0, 0.0), [('N', 'N')])
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
        self.assertEqual(self._resnames(incomplete_backbone_resindices(self.ag)),
                         ['ALA A3'])

    def test_complete_residues_and_ions_survive(self):
        kept = self._resnames(drop_prepwizard_hazard_residues(self.ag))
        self.assertIn('ALA A1', kept)
        self.assertIn('ALA A2', kept)
        self.assertIn('CYT B201', kept)
        self.assertIn('ZN B203', kept)

    def test_single_carbon_het_is_dropped_but_monatomic_ion_is_not(self):
        dropped = self._resnames(single_carbon_het_resindices(self.ag))
        self.assertEqual(dropped, ['CF0 B202'])

    def test_renamed_free_ligand_is_restored(self):
        snapshot = snapshot_restorable_ligands(self.pdb)
        prepped = os.path.join(self.tmp, 'prepped.pdb')
        with open(self.pdb) as fh, open(prepped, 'w') as out:
            for line in fh:
                if line[17:20].strip() == 'CYT':
                    line = 'ATOM  ' + line[6:17] + 'CYS' + line[20:]
                out.write(line)

        self.assertEqual(restore_renamed_ligands(prepped, snapshot), 1)
        with open(prepped) as fh:
            restored = [l for l in fh if l[22:26].strip() == '201']
        self.assertTrue(all(l.startswith('HETATM') for l in restored))
        self.assertTrue(all(l[17:20] == 'CYT' for l in restored))

    def test_chain_bonded_residue_is_not_restorable(self):
        bonded = os.path.join(self.tmp, 'bonded.pdb')
        with open(self.pdb) as fh:
            lines = fh.readlines()
        lines = [l for l in lines if l[22:26].strip() != '201']
        lines += _residue(30, 'CR8', 'B', 201, (4.4, 0.0, 0.0), LIGAND, het=True)
        with open(bonded, 'w') as fh:
            fh.writelines(lines)
        self.assertNotIn('CR8', set(snapshot_restorable_ligands(bonded).values()))

    def test_restore_is_a_no_op_when_nothing_was_renamed(self):
        snapshot = snapshot_restorable_ligands(self.pdb)
        untouched = os.path.join(self.tmp, 'untouched.pdb')
        shutil.copy(self.pdb, untouched)
        with open(untouched) as fh:
            before = fh.read()
        self.assertEqual(restore_renamed_ligands(untouched, snapshot), 0)
        with open(untouched) as fh:
            self.assertEqual(fh.read(), before)

class ModifiedResidueTests(unittest.TestCase):
    CCD = {'KCX': 'L-PEPTIDE LINKING', 'DCY': 'D-PEPTIDE LINKING', 'PLP': 'NON-POLYMER',
           'SAM': 'NON-POLYMER', '04C': 'peptide-like', 'HOH': 'NON-POLYMER'}

    @staticmethod
    def _res(serial, resname, chain, resnum, x0, names, het, extra=()):
        lines = [_atom(serial + i, n, resname, chain, resnum, (x0 + 1.3 * i, 0.0, 0.0),
                       n[0], het) for i, n in enumerate(names)]
        return lines + [_atom(serial + len(names) + j, n, resname, chain, resnum,
                              (x0 + dx, dy, 0.0), el, het)
                        for j, (n, el, dx, dy) in enumerate(extra)]

    def _structure(self):
        from ligand_vdgs.preprocessing._prep_filters import modified_residues_to_protein
        L = []
        bb = ['N', 'CA', 'C']
        L += self._res(1, 'ALA', 'A', 1, 0.0, bb, False)
        L += self._res(10, 'MSE', 'A', 2, 3.9, bb, True, [('SE', 'SE', 1.3, 2.0)])
        L += self._res(20, 'KCX', 'A', 3, 7.8, bb, True, [('NZ', 'N', 1.3, 4.0),
                                                         ('CX', 'C', 1.3, 5.3)])
        L += self._res(30, 'ALA', 'A', 4, 11.7, bb, False)
        L += self._res(40, 'SEC', 'A', 5, 15.6, bb, True, [('SE', 'SE', 1.3, 2.0)])
        L += self._res(50, 'DCY', 'A', 6, 19.5, bb, True, [('SG', 'S', 1.3, 2.0)])
        L += [_atom(60, 'C4A', 'PLP', 'A', 301, (9.1, 5.5, 0.0), 'C', True),
              _atom(61, 'N1', 'PLP', 'A', 301, (9.1, 6.8, 0.0), 'N', True)]
        L += self._res(70, 'SAM', 'A', 302, 14.3 + 1.6, ['N', 'CA', 'C'], True)
        L += self._res(80, 'MSE', 'B', 1, 100.0, bb + ['OXT'], True,
                       [('SE', 'SE', 1.3, 2.0)])
        L += self._res(90, 'KCX', 'C', 1, 200.0, bb + ['OXT'], True)
        L += [_atom(99, 'C1', 'ZZZ', 'A', 401, (-1.5, 0.0, 0.0), 'C', True)]
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
        self.assertEqual(rec[('A', 2)][:2], ('ATOM', 'MET'))
        self.assertIn('SE', names[('A', 2)])
        self.assertEqual(rec[('A', 3)][:2], ('ATOM', 'KCX'))
        self.assertEqual(rec[('A', 5)][:2], ('ATOM', 'CYS'))
        self.assertEqual(rec[('A', 6)][:2], ('ATOM', 'DCY'))
        self.assertEqual(rec[('A', 301)][:2], ('HETATM', 'PLP'))
        self.assertEqual(rec[('A', 302)][:2], ('HETATM', 'SAM'))
        self.assertEqual(rec[('B', 1)][:2], ('HETATM', 'MET'))
        self.assertEqual(rec[('C', 1)][:2], ('HETATM', 'KCX'))
        self.assertEqual(rec[('A', 401)][:2], ('HETATM', 'ZZZ'))
        self.assertEqual(rec[('A', 402)][:2], ('HETATM', '04C'))
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

def _atom2(serial, name, resname, chain, resnum, xyz, element, het=False):
    return (f"{'HETATM' if het else 'ATOM  '}{serial:5d} {name:<4s}{resname:>4s}"
            f"{chain:>2s}{resnum:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00{element:>12s}\n")

class TextRepairOrderTests(unittest.TestCase):
    CCD = {'KCX': 'L-PEPTIDE LINKING', 'ATP': 'NON-POLYMER'}

    def _structure(self):
        bb = ['N', 'CA', 'C']
        lines = []
        for i, (resname, het, extra) in enumerate(
                [('ALA', False, ()), ('MSE', True, (('SE', 'SE'),)),
                 ('KCX', True, (('NZ', 'N'),))]):
            x0 = 3.9 * i
            for j, name in enumerate(bb):
                lines.append(_atom2(1 + 10 * i + j, name, resname, 'AA', i + 1,
                                    (x0 + 1.3 * j, 0.0, 0.0), name[0], het))
            for k, (name, el) in enumerate(extra):
                lines.append(_atom2(5 + 10 * i + k, name, resname, 'AA', i + 1,
                                    (x0 + 1.3, 2.0, 0.0), el, het))
        for i in range(3):
            for j, name in enumerate(bb):
                lines.append(_atom2(100 + 10 * i + j, name, 'ALA', 'A', i + 1,
                                    (50.0 + 3.9 * i + 1.3 * j, 0.0, 0.0), name[0]))
        return lines

    def _repaired(self):
        from ligand_vdgs.preprocessing._prep_filters import text_repairs
        return text_repairs(self._structure(), self.CCD)

    def test_all_three_rewrites_and_the_remap_land_together(self):
        new, mapping, stats = self._repaired()
        self.assertEqual(mapping, {'AA': 'B'})
        rec = {(l[20:22], l[22:26].strip()): (l[:6].strip(), l[17:20]) for l in new}
        self.assertEqual(rec[(' B', '2')], ('ATOM', 'MET'))
        self.assertEqual(rec[(' B', '3')], ('ATOM', 'KCX'))
        self.assertEqual(rec[(' A', '1')], ('ATOM', 'ALA'))
        self.assertIn('SE', {l[12:16].strip() for l in new if l[20:22] == ' B'
                             and l[22:26].strip() == '2'})
        self.assertEqual(stats, {'renamed': 1, 'amino_acids_to_atom': 1,
                                 'modres_to_atom': 1, 'unknown_resnames': []})

    def test_column_21_is_blank_everywhere_afterwards(self):
        new, _, _ = self._repaired()
        self.assertEqual({l[20] for l in new}, {' '})

    def test_remap_runs_first_so_chains_do_not_merge(self):
        new, _, _ = self._repaired()
        keys = {(l[20:22], l[22:26]) for l in new}
        self.assertIn((' B', '   1'), keys)
        self.assertIn((' A', '   1'), keys)

    def test_clean_structure_is_returned_unchanged(self):
        from ligand_vdgs.preprocessing._prep_filters import text_repairs
        lines = [_atom2(1 + i, n, 'ALA', 'A', 1, (1.3 * i, 0.0, 0.0), n[0])
                 for i, n in enumerate(['N', 'CA', 'C'])]
        new, mapping, stats = text_repairs(lines, self.CCD)
        self.assertEqual(new, lines)
        self.assertEqual(mapping, {})
        self.assertEqual(stats['renamed'] + stats['modres_to_atom'], 0)

def _repair_module():
    import importlib.util
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(repo, 'scripts', 'remap_chain_ids.py')
    spec = importlib.util.spec_from_file_location('_remap_chain_ids', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

class RepairPassTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.src = os.path.join(self.tmp, 'src')
        self.out = os.path.join(self.tmp, 'out')
        bb = ['N', 'CA', 'C']
        clean = [_atom2(1 + i, n, 'ALA', 'A', 1, (1.3 * i, 0.0, 0.0), n[0])
                 for i, n in enumerate(bb)]
        dirty = TextRepairOrderTests()._structure()
        for stem, lines in (('1abc', clean), ('2abd', clean), ('3xya', dirty)):
            d = os.path.join(self.src, stem[1:3].lower())
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, stem + '.pdb'), 'w') as f:
                f.writelines(lines)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, *extra):
        import subprocess
        import sys
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = [sys.executable, os.path.join(repo, 'scripts', 'remap_chain_ids.py'),
               '--pdb-dir', self.src, '--out-dir', self.out, *extra]
        return subprocess.run(cmd, capture_output=True, text=True, cwd=repo)

    def _manifest(self):
        path = os.path.join(self.tmp, 'm.tsv')
        with open(path) as f:
            rows = [l.rstrip('\n').split('\t') for l in f]
        return {r[0]: r[1:] for r in rows[1:]}

    def test_all_mode_writes_every_structure_and_copies_the_clean_ones_verbatim(self):
        res = self._run('--all', '--manifest', os.path.join(self.tmp, 'm.tsv'))
        self.assertEqual(res.returncode, 0, res.stderr)
        written = {f for _, _, fs in os.walk(self.out) for f in fs}
        self.assertEqual(written, {'1abc.pdb', '2abd.pdb', '3xya.pdb'})
        for stem in ('1abc', '2abd'):
            rel = os.path.join(stem[1:3].lower(), stem + '.pdb')
            with open(os.path.join(self.src, rel)) as a, open(os.path.join(self.out, rel)) as b:
                self.assertEqual(a.read(), b.read())
        man = self._manifest()
        self.assertEqual(man['1abc'][0], '0')
        self.assertEqual(man['3xya'][0], '1')
        self.assertEqual(man['3xya'][1], 'AA->B')

    def test_repaired_file_carries_the_remap_remark(self):
        self._run('--all')
        with open(os.path.join(self.out, 'xy', '3xya.pdb')) as f:
            head = f.readline()
        self.assertTrue(head.startswith('REMARK 900 CHAIN ID REMAPPED AA -> B'), head)

    def test_without_all_only_changed_structures_are_written(self):
        self._run()
        written = {f for _, _, fs in os.walk(self.out) for f in fs}
        self.assertEqual(written, {'3xya.pdb'})

    def test_shards_partition_the_database_exactly_once(self):
        for shard in range(3):
            self._run('--all', '--shard', str(shard), '--num-shards', '3',
                      '--manifest', os.path.join(self.tmp, f'm{shard}.tsv'))
        seen = []
        for shard in range(3):
            with open(os.path.join(self.tmp, f'm{shard}.tsv')) as f:
                seen += [l.split('\t')[0] for l in f.read().splitlines()[1:]]
        self.assertEqual(sorted(seen), ['1abc', '2abd', '3xya'])
        written = {f for _, _, fs in os.walk(self.out) for f in fs}
        self.assertEqual(written, {'1abc.pdb', '2abd.pdb', '3xya.pdb'})

    def test_an_unrepairable_structure_is_skipped_but_still_accounted_for(self):
        import itertools
        import string
        from ligand_vdgs.preprocessing._chain_ids import CHAIN_POOL
        ids = [a + b for a, b in itertools.islice(
            itertools.product(string.ascii_uppercase, repeat=2), len(CHAIN_POOL) + 1)]
        self.assertEqual(len(set(ids)), len(CHAIN_POOL) + 1)
        lines = [_atom2(i + 1, 'CA', 'ALA', cid, 1, (3.0 * i, 0.0, 0.0), 'C')
                 for i, cid in enumerate(ids)]
        d = os.path.join(self.src, 'ov')
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, '4ovf.pdb'), 'w') as f:
            f.writelines(lines)

        res = self._run('--all', '--manifest', os.path.join(self.tmp, 'm.tsv'))
        self.assertEqual(res.returncode, 0, res.stderr)
        written = {f for _, _, fs in os.walk(self.out) for f in fs}
        self.assertNotIn('4ovf.pdb', written)
        self.assertEqual(written, {'1abc.pdb', '2abd.pdb', '3xya.pdb'})
        man = self._manifest()
        self.assertEqual(set(man), {'1abc', '2abd', '3xya', '4ovf'})
        self.assertEqual(man['4ovf'][0], 'skipped_overflow')
        self.assertIn('single characters are available', man['4ovf'][5])
        self.assertIn('4ovf', res.stderr)

    def test_a_failed_write_costs_one_structure_not_the_whole_shard(self):
        rci = _repair_module()
        real = rci._copy_structure
        def flaky(src, out):
            if out.endswith('2abd.pdb'):
                raise FileExistsError(17, 'File exists')
            return real(src, out)
        rci._copy_structure = flaky
        try:
            rci.main(['--pdb-dir', self.src, '--out-dir', self.out, '--all',
                      '--manifest', os.path.join(self.tmp, 'm.tsv')])
        finally:
            rci._copy_structure = real

        man = self._manifest()
        self.assertEqual(set(man), {'1abc', '2abd', '3xya'})
        self.assertEqual(man['2abd'][0], 'skipped_error')
        self.assertIn('FileExistsError', man['2abd'][5])
        written = {f for _, _, fs in os.walk(self.out) for f in fs}
        self.assertEqual(written, {'1abc.pdb', '3xya.pdb'})

    def test_a_manifest_row_is_never_written_for_a_file_that_was_not(self):
        rci = _repair_module()
        real = rci._write_repaired
        def flaky(out, remarks, lines):
            if out.endswith('3xya.pdb'):
                raise OSError(28, 'No space left on device')
            return real(out, remarks, lines)
        rci._write_repaired = flaky
        try:
            rci.main(['--pdb-dir', self.src, '--out-dir', self.out, '--all',
                      '--manifest', os.path.join(self.tmp, 'm.tsv')])
        finally:
            rci._write_repaired = real

        man = self._manifest()
        self.assertEqual(man['3xya'][0], 'skipped_error')
        for stem, fields in man.items():
            path = os.path.join(self.out, stem[1:3].lower(), stem + '.pdb')
            if fields[0].startswith('skipped'):
                self.assertFalse(os.path.isfile(path), f'{stem} skipped but written')
            else:
                self.assertTrue(os.path.isfile(path), f'{stem} has a row but no file')

    def test_a_rerun_is_idempotent_and_leaves_no_scratch_files(self):
        self._run('--all', '--manifest', os.path.join(self.tmp, 'm1.tsv'))
        first = {}
        for root, _, fs in os.walk(self.out):
            for f in fs:
                with open(os.path.join(root, f)) as fh:
                    first[f] = fh.read()
        res = self._run('--all', '--manifest', os.path.join(self.tmp, 'm.tsv'))
        self.assertEqual(res.returncode, 0, res.stderr)
        second = {}
        for root, _, fs in os.walk(self.out):
            for f in fs:
                with open(os.path.join(root, f)) as fh:
                    second[f] = fh.read()
        self.assertEqual(first, second)
        self.assertFalse([f for f in second if f.endswith('.tmp')],
                         'atomic writes must not leave .tmp files behind')

    def test_a_racing_sibling_directory_creation_is_not_fatal(self):
        rci = _repair_module()
        real = os.makedirs
        state = {'raised': False}
        def racing(path, *a, **kw):
            real(path, *a, **kw)
            if not state['raised']:
                state['raised'] = True
                raise FileExistsError(17, 'File exists', path)
        os.makedirs = racing
        try:
            rci._ensure_dir(os.path.join(self.out, 'zz'))
        finally:
            os.makedirs = real
        self.assertTrue(state['raised'])
        self.assertTrue(os.path.isdir(os.path.join(self.out, 'zz')))

    def test_refuses_to_repair_in_place(self):
        import subprocess
        import sys
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        res = subprocess.run(
            [sys.executable, os.path.join(repo, 'scripts', 'remap_chain_ids.py'),
             '--pdb-dir', self.src, '--out-dir', self.src, '--all'],
            capture_output=True, text=True, cwd=repo)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('never repair in place', res.stderr)

if __name__ == '__main__':
    unittest.main()
