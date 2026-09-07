"""Covers the '|'-gated reclassification of unparseable metal complexes.

The CCD writes dative bonds as '|', which RDKit rejects, so metal complexes die at
`MolFromSmiles` before `undesired_elements` ever sees them and land in num_ligs_failed
rather than not_druglike. `undesired_elements_in_dative_smiles` recovers the
classification by dropping the '|' marks and re-parsing, so RDKit still does the element
parsing -- removing bond marks cannot change the atom set, which is all the check needs.
Anything still unparseable afterwards stays a reported failure.

Exists because the parent ligand set is expanding past the CCD, so this has to hold by
construction rather than by the zero-false-positive rate measured on the CCD.
Discriminating inputs first, all of which a string-level element matcher gets wrong:
strings that spell a metal only through unbracketed atoms ("In1cccc1" is iodine + aromatic
N, not indium), and bracket contents that merely start with an element symbol ("[Ala]" is
not aluminium). All must stay in the failure log. A genuine complex is the floor case.
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, '..', 'ligand_vdgs', 'generate_vdgs',
                      'fragment_database_ligs.py')
sys.path.insert(0, os.path.join(HERE, '..'))
from ligand_vdgs.generate_vdgs.fragment_database_ligs import (  # noqa: E402
    UNDESIRED_ELEMENTS, undesired_elements_in_dative_smiles as undesired_in_smiles)

# Rows RDKit cannot parse that must be reclassified as not-druglike.
METAL_ROWS = [('[Pt+]|1|2(|NCCN|1)|I[Pt+]|3(|NCCN|3)I|2', '0JC'),
              ('O|[Cu++]', '1CU')]
# Rows RDKit cannot parse that must REMAIN parse failures: the first three spell a metal
# only via unbracketed atoms, the fourth has bracket contents that merely begin with an
# element symbol, and the last two carry no undesired element at all.
DECOY_ROWS = [('In1cccc1|', 'FP1'),        # iodine + aromatic N, not indium
              ('CSnc1ccccc1|', 'FP2'),     # C-S-n, not tin
              ('CCa1ccccc1|', 'FP3'),      # not calcium
              ('CC|[N,Na]C', 'FP4'),       # bracket content is not an atom
              ('CC(|)C(=O)NCC', 'FP5'),    # '|' but no metal
              ('c1ccccc1(', 'FP6')]        # no '|' at all; gate must not fire
# Parseable rows: must be filtered by the normal path, not the '|' fallback.
PARSEABLE_UNDESIRED = [('CC(=O)O[Zn]OC(C)=O', 'ZN1'),
                       ('N[CH](Cc1c[nH]c2[se]ccc12)C(O)=O', '23S')]  # aromatic Se
DRUGLIKE_ROWS = [('CC(=O)Nc1ccc(O)cc1', 'TYL'),
                 ('CN1C=NC2=C1C(=O)N(C)C(=O)N2C', 'CFF')]


def run_script(rows, tmp):
    ccd = os.path.join(tmp, 'ccd.smi')
    with open(ccd, 'w') as fh:
        fh.writelines(f'{smi}\t{code}\n' for smi, code in rows)
    logfile = os.path.join(tmp, 'log')
    proc = subprocess.run(
        [sys.executable, SCRIPT, '--ccd', ccd, '--outdir', tmp, '--logfile', logfile],
        capture_output=True, text=True)
    with open(logfile) as fh:
        return proc, fh.read()


class DativeSmilesElementMatchTests(unittest.TestCase):
    """Unit-level: RDKit must do the element parsing, not a string heuristic."""

    def test_unbracketed_metal_spellings_do_not_match(self):
        # The exact false positive the raw-string approach was abandoned for.
        for smi in ('In1cccc1', 'CSnc1ccccc1', 'CCa1ccccc1', 'CKc1ccc1', 'CUc1ccc1'):
            self.assertFalse(undesired_in_smiles(smi), smi)

    def test_bracket_syntax_is_not_mistaken_for_an_element(self):
        for smi in ('[13C]H4', '[2H+]', '[C@@H](N)C', '[nH]1cccc1', '[N+](=O)[O-]',
                    '[C:1]CC', '[#6][#7]', '[*]C', '[B](F)(F)n1cccc1'):
            self.assertFalse(undesired_in_smiles(smi), smi)

    def test_bracketed_undesired_elements_match(self):
        # Includes isotope/charge decoration and the lowercase aromatic metalloids.
        for smi in ('[Fe]', '[Sn+2]', '[Y]', '[K+]', '[Yb+3]',
                    '[se]1cccc1', '[te]1cccc1', '[as]1ccccc1'):
            self.assertTrue(undesired_in_smiles(smi), smi)

    def test_unparseable_garbage_is_never_reclassified(self):
        # Strings that stay unparseable after dropping '|' must fall through to the
        # failure log. '[Ala]' is the discriminating case: a string-level matcher
        # reading bracket contents sees "Al" and wrongly calls it aluminium.
        for smi in ('[N,Na]', '[Ala]|C', 'CC(|)C(=O)NCC', 'c1ccccc1('):
            self.assertFalse(undesired_in_smiles(smi), smi)

    def test_boron_is_not_undesired(self):
        # Boronic-acid warheads are deliberately mined; guards the element table.
        self.assertNotIn('B', UNDESIRED_ELEMENTS)


class ReclassificationEndToEndTests(unittest.TestCase):
    def test_metals_reclassified_and_decoys_stay_failures(self):
        rows = METAL_ROWS + DECOY_ROWS + PARSEABLE_UNDESIRED + DRUGLIKE_ROWS
        with tempfile.TemporaryDirectory() as tmp:
            proc, log = run_script(rows, tmp)
        self.assertEqual(proc.returncode, 0, proc.stderr[-2000:])

        failed_block, _, metal_block = log.partition(
            'Unparseable rows carrying an undesired element')
        # Every decoy must be visible as a parse failure...
        for _, code in DECOY_ROWS:
            self.assertIn(code, failed_block, f'{code} should be a logged failure')
            self.assertNotIn(code, metal_block, f'{code} must not be reclassified')
        # ...and every genuine complex must be reclassified, not counted as a failure.
        for _, code in METAL_ROWS:
            self.assertIn(code, metal_block, f'{code} should be reclassified')
            self.assertNotIn(code, failed_block, f'{code} must not be a parse failure')

        self.assertIn(f'Number of mols unsuccessfully parsed: {len(DECOY_ROWS)}', log)
        self.assertIn('...of which were unparseable (bad SMILES syntax): '
                      f'{len(METAL_ROWS)}', log)
        # Reclassified rows plus the two parseable undesired ligands.
        self.assertIn('Number of mols skipped b/c not druglike: '
                      f'{len(METAL_ROWS) + len(PARSEABLE_UNDESIRED)}', log)
        self.assertIn(f'Number of mols parsed: {len(DRUGLIKE_ROWS)}', log)

    def test_reclassified_rows_alone_do_not_write_an_empty_pickle(self):
        # The gate must not let a file of all-metal rows slip past the empty-output
        # guard by draining num_ligs_failed.
        with tempfile.TemporaryDirectory() as tmp:
            proc, log = run_script(METAL_ROWS, tmp)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(os.path.exists(os.path.join(tmp, 'database_frags_dict.pkl')))
        for _, code in METAL_ROWS:
            self.assertIn(code, log)

    def test_inert_without_the_dative_bond_notation(self):
        # Contract for a ligand source that does not use '|': an unparseable metal
        # complex stays a parse failure. Note this is guaranteed by the implementation
        # (dropping '|' from a string without one re-parses to the same None), not by
        # the "'|' in smiles" guard, which is only a short-circuit -- so this test will
        # not catch that guard being removed.
        rows = [('[Fe](C', 'BAD')] + DRUGLIKE_ROWS
        with tempfile.TemporaryDirectory() as tmp:
            proc, log = run_script(rows, tmp)
        self.assertEqual(proc.returncode, 0, proc.stderr[-2000:])
        self.assertIn('Number of mols unsuccessfully parsed: 1', log)
        self.assertNotIn('Unparseable rows carrying an undesired element', log)


if __name__ == '__main__':
    unittest.main()
