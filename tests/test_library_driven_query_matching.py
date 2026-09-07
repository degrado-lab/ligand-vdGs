"""Hit finding enumerates the *library's* fragment keys against the query ligand.

The discriminating cases are the ones the old path (re-fragmenting the query with
hardcoded radius 2 / size 4-5) could not reach: a key larger than that window, a
key whose ring primitive is only satisfiable on the intact ligand, a charged
query matched by the neutral key it was collapsed onto, and a symmetric CG whose
extra labelings come from the library's stored automorphism group rather than
from RDKit's matching.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions.utils import (identify_mol_automorphisms,
                                         mol_from_fragment, smiles_to_filename)
from ligand_vdgs.score_poses import hit_finder_core as hf

RDLogger.DisableLog('rdApp.*')


def ligand_with_coords(smiles):
    """A 3D-embedded, H-stripped ligand mol with ring perception, as the hit
    finder's own query-ligand builder returns."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=0xf00d) == 0
    mol = Chem.RemoveHs(mol)
    Chem.GetSymmSSSR(mol)
    return mol


class LibraryDrivenMatchingTests(unittest.TestCase):
    def setUp(self):
        self.lib = tempfile.mkdtemp(prefix='vdglib_')
        self.addCleanup(shutil.rmtree, self.lib)

    def add_frag(self, key, complete=True, automorphisms=None, cg_smarts=None):
        """Materialize the parts of a fragment directory that hit finding reads."""
        db_name = smiles_to_filename(key)
        frag_dir = os.path.join(self.lib, db_name)
        os.makedirs(frag_dir, exist_ok=True)
        with open(os.path.join(frag_dir, f'{db_name}_log'), 'w') as fh:
            fh.write('Job completed.\n' if complete else 'Started.\n')
        if automorphisms is None:
            automorphisms = identify_mol_automorphisms(mol_from_fragment(key))
        vdg_npz.write_cg_symmetry(frag_dir, cg_smarts or key, automorphisms)
        return db_name

    def match(self, lig_smiles, warn=None):
        entries = set(os.listdir(self.lib))
        # Rebuild the per-process pattern cache: these tests swap libraries
        # inside one process, which a real worker never does.
        hf.init_worker(entries)
        return hf.match_library_frags_to_query(
            ligand_with_coords(lig_smiles), self.lib, entries, warn=warn)

    def match_collecting_warnings(self, lig_smiles):
        """(result, warnings). The warn-once state is per process, so it has to be
        cleared or a second test in the same process silently sees no warning."""
        hf._warned_charge_only_frags.clear()
        hf._warned_incomplete_frags.clear()
        warnings = []
        return self.match(lig_smiles, warn=warnings.append), warnings

    def test_key_larger_than_the_old_hardcoded_window_is_found(self):
        """A 6-atom key from a library built with --max-frag-size 6. Query-side
        enumeration capped at 5 could never produce this key, so the fragment
        was unreachable no matter what the ligand contained."""
        db_name = self.add_frag('[C;!R][C;!R](=[O;!R])[N;!R][C;!R][C;!R]')
        frags, _, in_lib = self.match('CCC(=O)NCCO')
        self.assertTrue(in_lib[db_name])
        self.assertEqual(len(frags[db_name]), 1)

    def test_ring_key_matches_only_the_ring_and_needs_the_intact_ligand(self):
        """`r6` is satisfiable only on the whole ligand: the same atoms cut out
        as a fragment form a chain. A ligand carrying both a THP and an acyclic
        ether must contribute exactly the ring site."""
        ring_key = self.add_frag('[C;r6][C;r6][O;r6][C;r6][C;r6]')
        chain_key = self.add_frag('[C;!R][C;!R][O;!R][C;!R][C;!R]')
        frags, _, _ = self.match('CCOCCC1CCOCC1')
        self.assertEqual(len(frags[ring_key]), 1)
        self.assertEqual(len(frags[chain_key]), 1)
        # The ring site's atoms are all in a ring of the query; the chain
        # site's are not.
        lig = ligand_with_coords('CCOCCC1CCOCC1')
        ring_atoms = frags[ring_key][0][0][2]
        chain_atoms = frags[chain_key][0][0][2]
        self.assertTrue(all(lig.GetAtomWithIdx(i).IsInRingSize(6)
                            for i in ring_atoms))
        self.assertFalse(any(lig.GetAtomWithIdx(i).IsInRing()
                             for i in chain_atoms))

    def test_neutral_key_matches_a_charged_query(self):
        """Protonation-variant collapse means the library holds only the neutral
        acid; a ligand drawn as a carboxylate must still reach it."""
        db_name = self.add_frag('[C;!R][C;!R](=[O;!R])[O;!R]')
        frags, _, _ = self.match('CCC(=O)[O-]')
        self.assertEqual(len(frags[db_name]), 1)

    def test_symmetric_cg_gets_every_stored_labeling(self):
        """RDKit's matching is not resonance-normalized: it finds one labeling
        of a carboxylate, and the library clustered under two."""
        key = '[C;!R][C;!R](=[O;!R])[O;!R]'
        automorphisms = [(0, 1, 2, 3), (0, 1, 3, 2)]
        db_name = self.add_frag(key, automorphisms=automorphisms)
        frags, _, _ = self.match('CCC(=O)O')
        site = frags[db_name][0]
        self.assertEqual(len(site), 2)
        # The two labelings differ by swapping the two oxygens, and the coords
        # follow the labeling rather than staying put.
        coords = [np.asarray(sub.GetConformer().GetPositions()) for sub, _, _ in site]
        self.assertFalse(np.allclose(coords[0], coords[1]))
        self.assertTrue(np.allclose(coords[0][[0, 1, 3, 2]], coords[1]))

    def test_repeated_and_overlapping_occurrences(self):
        """Two chemically identical sites stay two sites (they are scored
        separately), while overlapping matches of one site collapse into one."""
        db_name = self.add_frag('[C;!R][C;!R](=[O;!R])[O;!R]')
        frags, _, _ = self.match('OC(=O)CCCCCCC(=O)O')
        self.assertEqual(len(frags[db_name]), 2)
        atom_sets = {site[0][2] for site in frags[db_name]}
        self.assertEqual(len(atom_sets), 2)

    def test_incomplete_job_is_excluded_but_reported(self):
        db_name = self.add_frag('[C;!R][C;!R](=[O;!R])[O;!R]', complete=False)
        frags, _, in_lib = self.match('CCC(=O)O')
        self.assertNotIn(db_name, frags)
        self.assertFalse(in_lib[db_name])

    def test_absent_fragment_is_not_reported_as_a_library_miss(self):
        """frags_in_lib now describes library fragments found in the ligand, so a
        key that does not occur must not appear at all."""
        db_name = self.add_frag('[N;!R][S;!R](=[O;!R])=[O;!R]')
        frags, _, in_lib = self.match('CCC(=O)O')
        self.assertEqual(frags, {})
        self.assertNotIn(db_name, in_lib)

    def test_ligand_without_ring_perception_still_matches(self):
        """A ring primitive against a mol with no RingInfo is a RuntimeError
        ("RingInfo not initialized"), not a no-match, so the matcher must not
        rely on its caller having perceived rings."""
        db_name = self.add_frag('[C;r6][C;r6][O;r6][C;r6][C;r6]')
        entries = set(os.listdir(self.lib))
        hf.init_worker(entries)
        unperceived = Chem.MolFromSmiles('C1CCOCC1', sanitize=False)
        AllChem.Compute2DCoords(unperceived)
        frags, _, _ = hf.match_library_frags_to_query(
            unperceived, self.lib, entries)
        self.assertEqual(len(frags[db_name]), 1)

    def test_unsearchable_fragment_directory_is_reported(self):
        """A directory name RDKit cannot parse is unreachable -- the searched key
        set is exactly what parses -- so it must not fail silently. Files at the
        library root are not fragments and must not be reported."""
        os.makedirs(os.path.join(self.lib, 'not-a-smarts]['))
        self.add_frag('[C;!R][C;!R](=[O;!R])[O;!R]')
        with open(os.path.join(self.lib, 'fragment_aliases.tsv'), 'w') as fh:
            fh.write('# alias\trepresentative\n')
        warnings = []
        entries = set(os.listdir(self.lib))
        hf.init_worker(entries)
        hf._warned_unsearchable_libs.discard(self.lib)
        hf.match_library_frags_to_query(
            ligand_with_coords('CCC(=O)O'), self.lib, entries,
            warn=warnings.append)
        self.assertEqual(len(warnings), 1)
        self.assertIn('not-a-smarts][', warnings[0])
        self.assertNotIn('fragment_aliases.tsv', warnings[0])

    def test_query_atom_order_follows_recorded_cg_smarts(self):
        """The stored cg_smarts, not the directory name, is the order the CG
        coords are in -- so the query fragment must be built in that order."""
        db_name = self.add_frag('[C;!R][C;!R](=[O;!R])[N;!R]',
                                cg_smarts='[N;!R][C;!R](=[O;!R])[C;!R]',
                                automorphisms=[(0, 1, 2, 3)])
        frags, query_frag_map, _ = self.match('CCC(=O)N')
        sub = frags[db_name][0][0][0]
        self.assertEqual([a.GetSymbol() for a in sub.GetAtoms()],
                         ['N', 'C', 'O', 'C'])
        self.assertEqual(query_frag_map[db_name], '[N;!R][C;!R](=[O;!R])[C;!R]')

    def test_charged_only_key_warns_when_the_ligand_is_drawn_neutral(self):
        """The one gap protonation collapse leaves: the library holds a charged
        key with no neutral twin, the ligand carries that moiety neutral, and
        nothing matches. Silent without this warning -- a fragment that matches
        nothing is indistinguishable from one the ligand does not contain."""
        db_name = self.add_frag('[C;!R][N+;!R]([C;!R])[C;!R]')
        (frags, _, _), warnings = self.match_collecting_warnings('CN(C)C')
        self.assertNotIn(db_name, frags)
        self.assertEqual(len(warnings), 1)
        self.assertIn(db_name, warnings[0])
        self.assertIn('[C;!R][N;!R]([C;!R])[C;!R]', warnings[0])
        self.assertIn('No usable neutral twin', warnings[0])

    def test_charge_only_warning_distinguishes_a_present_neutral_twin(self):
        """With the neutral twin in the library the moiety is still searched, so
        the warning must not claim the ligand's occurrences are unreachable --
        only this key's own vdGs are."""
        charged = self.add_frag('[C;!R][N+;!R]([C;!R])[C;!R]')
        neutral = self.add_frag('[C;!R][N;!R]([C;!R])[C;!R]')
        (frags, _, _), warnings = self.match_collecting_warnings('CN(C)C')
        self.assertIn(neutral, frags)
        self.assertNotIn(charged, frags)
        self.assertEqual(len(warnings), 1)
        self.assertIn('is searched instead', warnings[0])

    def test_charge_only_warning_finds_a_twin_written_another_way(self):
        """Fragment keys inherit the parent's perception, so one fragment has
        several spellings. Looking the twin up by exact name would call it
        missing and overstate the gap."""
        self.add_frag('[C;!R][C;!R][O;!R][N+;!R]')
        # Same fragment as the neutralized key, written from the other end.
        twin = self.add_frag('[N;!R][O;!R][C;!R][C;!R]')
        (frags, _, _), warnings = self.match_collecting_warnings('CCONC')
        self.assertIn(twin, frags)
        self.assertEqual(len(warnings), 1)
        self.assertIn('is searched instead', warnings[0])

    def test_charge_only_warning_treats_an_incomplete_twin_as_absent(self):
        """A crashed twin job leaves a directory behind; it is not a fallback."""
        self.add_frag('[C;!R][N+;!R]([C;!R])[C;!R]')
        self.add_frag('[C;!R][N;!R]([C;!R])[C;!R]', complete=False)
        (_, _, _), warnings = self.match_collecting_warnings('CN(C)C')
        self.assertEqual(len(warnings), 2)
        self.assertTrue(any('No usable neutral twin' in w for w in warnings))
        self.assertTrue(any('has not completed' in w for w in warnings))

    def test_charged_key_does_not_warn_when_the_ligand_lacks_the_moiety(self):
        """The check must fire on a lost hit, not on every charged key in the
        library -- otherwise it is noise on every model and gets ignored."""
        self.add_frag('[C;!R][N+;!R]([C;!R])[C;!R]')
        (_, _, _), warnings = self.match_collecting_warnings('CCCCO')
        self.assertEqual(warnings, [])

    def test_neutral_key_matching_a_charged_query_does_not_warn(self):
        """The direction that already works must stay quiet."""
        db_name = self.add_frag('[C;!R][C;!R](=[O;!R])[O;!R]')
        (frags, _, _), warnings = self.match_collecting_warnings('CCC(=O)[O-]')
        self.assertEqual(len(frags[db_name]), 1)
        self.assertEqual(warnings, [])

    def test_charge_only_warning_is_emitted_once_per_fragment(self):
        """Per process, not per model: a run scores thousands of models against
        one library."""
        self.add_frag('[C;!R][N+;!R]([C;!R])[C;!R]')
        (_, _, _), warnings = self.match_collecting_warnings('CN(C)C')
        self.assertEqual(len(warnings), 1)
        more = []
        self.match('CN(C)C', warn=more.append)
        self.assertEqual(more, [])


if __name__ == '__main__':
    unittest.main()
