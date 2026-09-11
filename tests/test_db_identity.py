"""A parent database is identified by its contents, not by where it sits.

The library is meant to be copyable between machines, so the recorded path legitimately
goes stale while the database stays the same; and a different database can be written to
the path the old one had. Comparing paths gets both of those backwards. These tests pin
the two directions: same database at two paths passes, different database at the same
path fails.
"""
import json
import os
import shutil
import tempfile
import unittest

from ligand_vdgs.functions.db_identity import (
    IDENTITY_FILENAME, compute_identity, identity_of, read_identity, write_identity)


def _mirror(root, structures):
    """{stem: text} written in parent_db mirror layout."""
    for stem, text in structures.items():
        d = os.path.join(root, stem[1:3].lower())
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, stem + '.pdb'), 'w') as f:
            f.write(text)
    return root


class DbIdentityTests(unittest.TestCase):

    STRUCTURES = {'1abc': 'ATOM      1  N   ALA A   1\n',
                  '2abd': 'ATOM      1  CA  GLY B   7\n',
                  '3xyz': 'HETATM    1  C1  LIG C 301\n'}

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.a = _mirror(os.path.join(self.tmp, 'a'), self.STRUCTURES)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_the_same_database_at_another_path_has_the_same_identity(self):
        b = os.path.join(self.tmp, 'somewhere', 'else', 'b')
        shutil.copytree(self.a, b)
        self.assertEqual(compute_identity(self.a)['sha256'],
                         compute_identity(b)['sha256'])

    def test_a_different_database_at_the_same_path_has_a_different_identity(self):
        before = compute_identity(self.a)['sha256']
        with open(os.path.join(self.a, 'ab', '2abd.pdb'), 'a') as f:
            f.write('ATOM      2  C   GLY B   7\n')      # same path, more content
        self.assertNotEqual(before, compute_identity(self.a)['sha256'])

    def test_adding_or_removing_a_structure_changes_the_identity(self):
        before = compute_identity(self.a)['sha256']
        _mirror(self.a, {'4wxy': 'ATOM      1  N   SER D   2\n'})
        after_add = compute_identity(self.a)['sha256']
        self.assertNotEqual(before, after_add)
        os.remove(os.path.join(self.a, 'wx', '4wxy.pdb'))
        self.assertEqual(before, compute_identity(self.a)['sha256'])

    def test_renaming_a_structure_changes_the_identity(self):
        """Same bytes, same count, different roster -- size alone must not decide it."""
        before = compute_identity(self.a)['sha256']
        d = os.path.join(self.a, 'ab')
        os.rename(os.path.join(d, '2abd.pdb'), os.path.join(d, '2abe.pdb'))
        self.assertNotEqual(before, compute_identity(self.a)['sha256'])

    def test_non_structure_files_alongside_the_mirror_do_not_count(self):
        """The identity file itself, manifests and logs must not change the identity, or
        stamping the database would invalidate the stamp."""
        before = compute_identity(self.a)['sha256']
        write_identity(self.a, source_dir='/somewhere/original')
        with open(os.path.join(self.a, 'manifest.tsv'), 'w') as f:
            f.write('stem\tchanged\n')
        self.assertEqual(before, compute_identity(self.a)['sha256'])

    def test_the_recorded_source_dir_is_provenance_and_not_part_of_the_hash(self):
        """A copy of a database is the same database, so where it was produced cannot
        feed the identity."""
        a = write_identity(self.a, source_dir='/one/place')
        b = write_identity(self.a, source_dir='/a/different/place')
        self.assertEqual(a['sha256'], b['sha256'])
        self.assertEqual(b['repaired_from'], os.path.abspath('/a/different/place'))

    def test_identity_of_caches_and_the_cache_agrees_with_a_fresh_walk(self):
        self.assertIsNone(read_identity(self.a))
        first = identity_of(self.a)
        self.assertTrue(os.path.isfile(os.path.join(self.a, IDENTITY_FILENAME)))
        self.assertEqual(first['sha256'], compute_identity(self.a)['sha256'])
        self.assertEqual(identity_of(self.a)['sha256'], first['sha256'])

    def test_a_corrupt_or_versionless_cache_is_recomputed_not_trusted(self):
        path = os.path.join(self.a, IDENTITY_FILENAME)
        for bad in ('{ truncated', json.dumps({'sha256': 'x' * 64}),
                    json.dumps({'version': 999, 'sha256': 'x' * 64})):
            with open(path, 'w') as f:
                f.write(bad)
            self.assertIsNone(read_identity(self.a), bad)
            self.assertEqual(identity_of(self.a)['sha256'],
                             compute_identity(self.a)['sha256'])

    def test_a_read_only_mirror_still_yields_an_identity(self):
        """A shared database another group owns cannot be stamped, and must still be
        usable rather than raising."""
        mode = os.stat(self.a).st_mode
        os.chmod(self.a, 0o555)
        try:
            record = identity_of(self.a)
            self.assertEqual(record['sha256'], compute_identity(self.a)['sha256'])
            self.assertFalse(os.path.isfile(os.path.join(self.a, IDENTITY_FILENAME)))
        finally:
            os.chmod(self.a, mode)

    def test_structure_count_is_reported(self):
        self.assertEqual(compute_identity(self.a)['n_structures'], 3)


if __name__ == '__main__':
    unittest.main()
