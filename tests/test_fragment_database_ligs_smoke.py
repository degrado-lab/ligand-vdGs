"""End-to-end smoke test for fragment_database_ligs on a four-ligand CCD file.

Exists because the dedup branch of `record_frag` (second ligand sharing a key
with the first) had no coverage: a missing import there crashed a real build
while the unit suite stayed green. Discriminating inputs first -- a protonation
pair whose keys must stay distinct here and collapse only in
`group_protonation_variants`, and a ring/chain pair that must NOT share a key --
then one ordinary ligand as a floor.
"""
import os
import pickle
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, '..', 'ligand_vdgs', 'generate_vdgs', 'fragment_database_ligs.py')


class FragmentDatabaseLigsSmokeTests(unittest.TestCase):
    def test_builds_dict_with_ring_annotated_keys(self):
        # Runs against an empty CCD store: with a template available the SMILES
        # column is ignored (rebuild-notes 6), and these codes are real CCD
        # entries whose templates are different molecules (PIP is the PIPES
        # buffer, not piperidine). Template-driven keys are covered separately
        # in test_ccd_templates.py.
        rows = [('CC(=O)O', 'ACY', 'acetic acid'),
                ('CC(=O)[O-]', 'ACT', 'acetate'),
                ('C1CCNCC1', 'PIP', 'piperidine'),
                ('CCCCNC', 'BAM', 'N-methylbutylamine')]
        with tempfile.TemporaryDirectory() as tmp:
            ccd = os.path.join(tmp, 'ccd.smi')
            with open(ccd, 'w') as fh:
                fh.writelines('\t'.join(r) + '\n' for r in rows)
            proc = subprocess.run(
                [sys.executable, SCRIPT, '--ccd', ccd, '--outdir', tmp,
                 '--logfile', os.path.join(tmp, 'log')],
                capture_output=True, text=True,
                env={**os.environ,
                     'LIGAND_VDGS_CCD_DIR': os.path.join(tmp, 'no_ccd')})
            self.assertEqual(proc.returncode, 0, proc.stderr[-2000:])
            with open(os.path.join(tmp, 'database_frags_dict.pkl'), 'rb') as fh:
                frags = pickle.load(fh)
        keys = {k for sd in frags.values() for k in sd}
        self.assertIn('[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D1]', keys)
        self.assertIn('[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O-;!R;D1]', keys)   # distinct here, pooled later
        ring = {k for k in keys if ';r6;' in k}
        chain = {k for k in keys if k.count(';!R;') >= 4 and 'N' in k}
        self.assertTrue(ring and chain)
        self.assertEqual(ring & chain, set())


if __name__ == '__main__':
    unittest.main()
