"""Multi-character chain IDs (columns 21-22) collide under every single-column reader.

Real case: 8any in the 2026-09 database, chain ``A4`` residue 66 (ASP, HETATM) and
chain ``4`` residue 66 (PHE, ATOM) -- ProDy 2.6 gives both one resindex, so two
residues merge into one and the contact gate measures the merged pair. The miner no
longer drops such a structure outright (that was the Probe-line hash mismatch), which
makes the residue-level collision the remaining defect, not a secondary one.
"""
import io
import unittest

import prody as pr

from ligand_vdgs.preprocessing._chain_ids import (
    CHAIN_POOL, ChainIdOverflow, assert_single_char_chains, chain_id_mapping,
    multichar_chain_ids, parse_remap_remarks, remap_chain_ids, remap_remarks)


def _atom(serial, name, resname, chain2, resnum, x=0.0, het=False):
    """*chain2* is written verbatim into columns 21-22 (right-justified)."""
    rec = 'HETATM' if het else 'ATOM  '
    return (f"{rec}{serial:5d} {name:<4s}{resname:>4s}{chain2:>2s}{resnum:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           {name[0]}\n")


class RemapTest(unittest.TestCase):

    def setUp(self):
        # PHE 4 66 (single-char chain) and ASP A4 66 (two-char chain), plus a
        # chain 'A' that must keep its letter, and a ligand on chain 'AA'.
        self.lines = [
            _atom(1, 'N', 'PHE', '4', 66, 0.0),
            _atom(2, 'CA', 'PHE', '4', 66, 1.5),
            _atom(3, 'N', 'ASP', 'A4', 66, 10.0, het=True),
            _atom(4, 'CA', 'ASP', 'A4', 66, 11.5, het=True),
            _atom(5, 'N', 'GLY', 'A', 1, 20.0),
            _atom(6, 'C1', 'NAD', 'AA', 500, 30.0, het=True),
            'TER       7      NADAA 500\n',
        ]

    def test_prody_merges_the_collision_without_remap(self):
        ag = pr.parsePDBStream(io.StringIO(''.join(self.lines)))
        sel = ag.select('chain 4 and resnum 66')
        self.assertEqual(set(sel.getResnames()), {'PHE', 'ASP'})
        self.assertEqual(len(set(sel.getResindices())), 1)

    def test_remap_gives_every_chain_one_column_and_no_collision(self):
        new, mapping = remap_chain_ids(self.lines)
        self.assertEqual(multichar_chain_ids(new), [])
        # 'A' and '4' are taken by real chains, so A4 -> B and AA -> C.
        self.assertEqual(mapping, {'A4': 'B', 'AA': 'C'})
        ag = pr.parsePDBStream(io.StringIO(''.join(new)))
        self.assertEqual(set(ag.select('chain 4 and resnum 66').getResnames()), {'PHE'})
        self.assertEqual(set(ag.select('chain B and resnum 66').getResnames()), {'ASP'})
        self.assertEqual(set(ag.select('chain C').getResnames()), {'NAD'})
        # Everything outside columns 21-22 is untouched, TER included.
        for old, upd in zip(self.lines, new):
            self.assertEqual(old[:20] + old[22:], upd[:20] + upd[22:])
        self.assertEqual(new[-1][20:22], ' C')

    def test_untouched_file_is_returned_as_is(self):
        clean = [l for l in self.lines if l[20] == ' ']
        new, mapping = remap_chain_ids(clean)
        self.assertIs(new, clean)
        self.assertEqual(mapping, {})
        assert_single_char_chains(clean)
        with self.assertRaises(ValueError):
            assert_single_char_chains(self.lines, 'x')

    def test_remarks_round_trip(self):
        _, mapping = remap_chain_ids(self.lines)
        self.assertEqual(parse_remap_remarks(remap_remarks(mapping)), mapping)

    def test_overflow_raises_instead_of_colliding(self):
        # Every single character in use, then one more two-character chain.
        lines = [_atom(i, 'CA', 'GLY', c, 1) for i, c in enumerate(CHAIN_POOL)]
        lines.append(_atom(99, 'CA', 'GLY', 'ZZ', 1))
        with self.assertRaises(ChainIdOverflow):
            chain_id_mapping(lines)
        # One short of the pool still fits.
        lines = [_atom(i, 'CA', 'GLY', c, 1) for i, c in enumerate(CHAIN_POOL[1:])]
        lines.append(_atom(99, 'CA', 'GLY', 'ZZ', 1))
        self.assertEqual(chain_id_mapping(lines), {'ZZ': CHAIN_POOL[0]})



from ligand_vdgs.preprocessing._prep_filters import embedded_amino_acid_hetatm_to_atom


def _residue(serial, resname, chain, resnum, x0, het, oxt=False):
    """N, CA, C (and OXT) along x, 1.3 A apart, so residue i's C is 1.3 A from residue
    i+1's N when consecutive residues are placed 3.9 A apart."""
    names = ['N', 'CA', 'C'] + (['OXT'] if oxt else [])
    return [_atom(serial + i, n, resname, chain, resnum, x0 + 1.3 * i, het=het)
            for i, n in enumerate(names)]


class HetatmProteinTest(unittest.TestCase):

    def test_bonded_and_fragment_residues_become_atom_but_free_ligand_stays(self):
        lines = []
        # A HETATM dipeptide (both residues bonded to each other): rewritten.
        lines += _residue(1, 'ALA', 'B', 1, 0.0, het=True)
        lines += _residue(5, 'TRP', 'B', 2, 3.9, het=True)
        # A HETATM residue bonded only to an ATOM residue: rewritten.
        lines += _residue(10, 'GLY', 'C', 5, 50.0, het=False)
        lines += _residue(15, 'ASP', 'C', 6, 53.9, het=True)
        # An isolated HETATM residue with no OXT: a trimmed chain fragment, rewritten.
        lines += _residue(20, 'ARG', 'E', 120, 100.0, het=True)
        # An isolated HETATM residue with OXT: a free amino-acid ligand, kept.
        lines += _residue(25, 'GLU', 'A', 104, 200.0, het=True, oxt=True)
        # A C-terminal residue (has OXT) bonded into a chain: rewritten.
        lines += _residue(30, 'LYS', 'C', 7, 57.8, het=True, oxt=True)
        # An isolated ATOM residue and a non-amino-acid HETATM: untouched.
        lines += _residue(35, 'SER', 'F', 1, 300.0, het=False)
        lines.append(_atom(40, 'C1', 'NAD', 'L', 500, 400.0, het=True))

        new, n = embedded_amino_acid_hetatm_to_atom(lines)
        self.assertEqual(n, 5)
        rec = {(l[21], int(l[22:26])): l[:6].strip() for l in new}
        self.assertEqual(rec[('B', 1)], 'ATOM')
        self.assertEqual(rec[('B', 2)], 'ATOM')
        self.assertEqual(rec[('C', 6)], 'ATOM')
        self.assertEqual(rec[('E', 120)], 'ATOM')
        self.assertEqual(rec[('C', 7)], 'ATOM')
        self.assertEqual(rec[('A', 104)], 'HETATM')
        self.assertEqual(rec[('F', 1)], 'ATOM')
        self.assertEqual(rec[('L', 500)], 'HETATM')
        for old, upd in zip(lines, new):
            self.assertEqual(old[6:], upd[6:])

    def test_nothing_to_do_returns_input(self):
        lines = _residue(1, 'ALA', 'A', 1, 0.0, het=False)
        new, n = embedded_amino_acid_hetatm_to_atom(lines)
        self.assertIs(new, lines)
        self.assertEqual(n, 0)


if __name__ == '__main__':
    unittest.main()
