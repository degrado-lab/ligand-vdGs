"""One policy decides whether a recorded input mismatch is fatal, for both records.

Two things carry recorded inputs: the library's own provenance JSON (checked on a
--include-only top-up) and the --frag-cost-estimate TSV header (checked on every run).
They used to disagree about the same field, so the same mismatch raised on one path and
warned on the other. check_recorded_inputs is the single policy; these tests pin the
matrix so the two cannot drift apart again.

    field                  full build      --include-only
    db_identity            raise           raise
    pdb_dir                warn            warn
    frags_dict_sha256      raise           warn
    max_size               per caller: fatal for the library, a warning for the estimate
"""
import argparse
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

from ligand_vdgs.generate_vdgs.estimate_frag_cost import read_estimate_header
from ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags import (
    PROVENANCE_FILENAME, check_estimate_header, check_provenance, write_provenance)

from ligand_vdgs.functions.db_identity import identity_of


def _mirror(root, structures):
    for stem, text in structures.items():
        d = os.path.join(root, stem[1:3].lower())
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, stem + '.pdb'), 'w') as f:
            f.write(text)
    return root


class ProvenancePolicyTests(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.frags = os.path.join(self.tmp, 'frags.pkl')
        with open(self.frags, 'wb') as f:
            f.write(b'not really a pickle, only hashed')
        self.lib = os.path.join(self.tmp, 'lib')
        # Two real mirrors: the same database at two paths, and a different one.
        contents = {'1abc': 'ATOM      1  N   ALA A   1\n',
                    '2abd': 'ATOM      1  CA  GLY B   7\n'}
        self.old = _mirror(os.path.join(self.tmp, 'old'), contents)
        self.moved = _mirror(os.path.join(self.tmp, 'moved', 'old'), contents)
        self.other = _mirror(os.path.join(self.tmp, 'other'),
                             dict(contents, **{'3xyz': 'HETATM    1  C1  LIG C 301\n'}))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _args(self, pdb_dir=None, max_size=10, include_only=False, frags_dict=None):
        pdb_dir = self.old if pdb_dir is None else pdb_dir
        return argparse.Namespace(
            frag_cost_estimate=os.path.join(self.tmp, 'est.tsv'),
            frags_dict=frags_dict or self.frags, pdb_dir=pdb_dir, max_size=max_size,
            include_only=include_only, min_instances=100)

    def _header(self, **over):
        """An estimate TSV header recorded against self.old at --max-size 10."""
        path = os.path.join(self.tmp, 'est.tsv')
        fields = {'pdb_dir': self.old, 'max_size': '10',
                  'pdb_db_identity': identity_of(self.old)['sha256']}
        fields.update(over)
        with open(path, 'w') as f:
            for k, v in fields.items():
                f.write(f'# {k}\t{v}\n')
            f.write('[C;!R][C;!R][C;!R][C;!R]\t1\t1\n')
        return read_estimate_header(path)

    def _other_dict(self, content=b'a different fragment dict'):
        path = os.path.join(self.tmp, 'other.pkl')
        with open(path, 'wb') as f:
            f.write(content)
        return path

    @staticmethod
    def _warnings(fn):
        buf = io.StringIO()
        with redirect_stdout(buf):
            fn()
        return buf.getvalue()

    # --- db identity: fatal on both paths; the path itself only warns ------------

    def test_estimate_from_another_database_raises_on_a_full_build(self):
        with self.assertRaises(SystemExit) as cm:
            check_estimate_header(self._header(), self._args(pdb_dir=self.other))
        self.assertIn('parent database contents', str(cm.exception))

    def test_estimate_from_another_database_raises_under_include_only_too(self):
        """A top-up grows the fragment dict, never the parent database."""
        with self.assertRaises(SystemExit):
            check_estimate_header(self._header(),
                                  self._args(pdb_dir=self.other, include_only=True))

    def test_library_provenance_rejects_a_top_up_from_another_database(self):
        write_provenance(self.lib, self._args(pdb_dir=self.old))
        with self.assertRaises(SystemExit) as cm:
            check_provenance(self.lib,
                             self._args(pdb_dir=self.other, include_only=True))
        self.assertIn('parent database contents', str(cm.exception))

    def test_the_same_database_at_another_path_passes_with_a_warning(self):
        """The workflow this exists for: the library is copied to another machine and
        topped up against a local mirror of the same database at a different path."""
        write_provenance(self.lib, self._args(pdb_dir=self.old))
        out = self._warnings(lambda: check_provenance(
            self.lib, self._args(pdb_dir=self.moved, include_only=True)))
        self.assertIn('WARNING', out)
        self.assertIn('Informational only', out)

        out = self._warnings(
            lambda: check_estimate_header(self._header(),
                                          self._args(pdb_dir=self.moved)))
        self.assertIn('Informational only', out)

    def test_the_same_database_at_the_same_path_is_silent(self):
        """Path spelling must not produce noise: --pdb-dir's own default carries a
        trailing slash."""
        from ligand_vdgs.functions.utils import file_sha256
        complete = self._header(frags_dict_sha256=file_sha256(self.frags))
        for spelling in (self.old, self.old + '/', os.path.join(self.old, 'ab', '..')):
            out = self._warnings(lambda: check_estimate_header(
                complete, self._args(pdb_dir=spelling)))
            self.assertEqual(out, '', spelling)

    # --- frags_dict_sha256: fatal on a full build, a warning on a top-up ----------

    def test_a_changed_frags_dict_raises_on_a_full_build(self):
        header = self._header(frags_dict_sha256='0' * 64)
        with self.assertRaises(SystemExit) as cm:
            check_estimate_header(header, self._args())
        self.assertIn('sha256', str(cm.exception))

    def test_a_changed_frags_dict_only_warns_under_include_only(self):
        """Top-up mode grows the dict by design; select_fragments is the real guard,
        raising if the estimate does not cover an included fragment."""
        header = self._header(frags_dict_sha256='0' * 64)
        out = self._warnings(
            lambda: check_estimate_header(header, self._args(include_only=True)))
        self.assertIn('WARNING', out)

        write_provenance(self.lib, self._args())
        out = self._warnings(lambda: check_provenance(
            self.lib, self._args(include_only=True, frags_dict=self._other_dict())))
        self.assertIn('WARNING', out)

    # --- max_size: fatal for the library, a warning for the estimate --------------

    def test_max_size_is_fatal_for_the_library_and_a_warning_for_the_estimate(self):
        """--max-size changes prepare_fragments' representative/alias mapping, so the
        same SMARTS can resolve to a key the rest of the library does not use. The
        estimate is only a set of counts, and main()'s membership check already catches
        the case that actually breaks it."""
        out = self._warnings(
            lambda: check_estimate_header(self._header(max_size='99'), self._args()))
        self.assertIn('WARNING', out)

        write_provenance(self.lib, self._args(max_size=10))
        with self.assertRaises(SystemExit) as cm:
            check_provenance(self.lib, self._args(max_size=99, include_only=True))
        self.assertIn('max-size', str(cm.exception))

    # --- unknown, not mismatched --------------------------------------------------

    def test_a_record_written_before_identities_existed_warns_instead_of_raising(self):
        header = self._header()
        header.pop('pdb_db_identity')
        out = self._warnings(lambda: check_estimate_header(
            header, self._args(pdb_dir=self.other)))
        self.assertIn('no parent-database identity', out)

    def test_a_missing_provenance_file_does_not_block_an_older_library(self):
        out = self._warnings(
            lambda: check_provenance(self.lib, self._args(include_only=True)))
        self.assertIn('WARNING', out)

    def test_provenance_round_trips_and_a_matching_run_is_silent(self):
        write_provenance(self.lib, self._args())
        with open(os.path.join(self.lib, PROVENANCE_FILENAME)) as f:
            record = json.load(f)
        self.assertEqual(record['parent_pdb_dir'], os.path.abspath(self.old))
        self.assertEqual(record['parent_db_identity'], identity_of(self.old)['sha256'])
        out = self._warnings(
            lambda: check_provenance(self.lib, self._args(include_only=True)))
        self.assertEqual(out, '')


if __name__ == '__main__':
    unittest.main()
