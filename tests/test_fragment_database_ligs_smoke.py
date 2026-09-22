"""End-to-end smoke test for fragment_database_ligs on a synthetic roster.

Exists because the parts that only run in `main` -- the roster load, the two
joins, the collapse, the atomic write and the wrapper's provenance fields -- have
no coverage from the unit tests, and a missing import or a wrong key in the
payload crashes a multi-hour build while the unit suite stays green.

Discriminating inputs first, then a floor case:
  * BGC's 6-ring, cut to 5 atoms, must key as an open arc, never a closed ring;
    ADP's true 5-ring (ribofuranose) must key closed -- the rotation split this
    enumerator exists to fix;
  * two different CCD codes in ONE biounit sharing fragments, so a summing join
    would show up as doubled support;
  * the same code repeated in one biounit (NCS), which must not multiply support;
  * a partially observed instance, whose fragments touching the missing atom must
    not be counted.

The roster is synthetic (real CCD codes, hand-written observed-atom sets) so the
test needs no parent database, only the CCD template store.
"""
import os
import pickle
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, '..', 'ligand_vdgs', 'generate_vdgs',
                      'fragment_database_ligs.py')
sys.path.insert(0, os.path.join(HERE, '..'))

from ligand_vdgs.functions import Frags, ccd_templates, ligand_perception  # noqa: E402


def _heavy_names(resname):
    template = ccd_templates.get_template(resname)
    return tuple(sorted(a.name for a in template.atoms if a.element not in ('H', 'D')))


def _write_roster(path, types, biounits):
    with open(path, 'wb') as fh:
        pickle.dump({'db_identity': {'sha256': 'smoke', 'version': 1},
                     'ccd_identity': ccd_templates.store_identity(),
                     'pdb_dir': '/nonexistent',
                     'partial': False,
                     'types': types,
                     'biounits': biounits,
                     'stats': {}}, fh)


def _run(tmp, roster_path, extra=()):
    logfile = os.path.join(tmp, 'log')
    proc = subprocess.run(
        [sys.executable, SCRIPT, '--roster', roster_path, '--outdir', tmp,
         '--logfile', logfile, *extra],
        capture_output=True, text=True)
    payload = None
    out_path = os.path.join(tmp, 'database_frags_dict.pkl')
    if os.path.exists(out_path):
        with open(out_path, 'rb') as fh:
            payload = pickle.load(fh)
    log = open(logfile).read() if os.path.exists(logfile) else ''
    return proc, payload, log


class FragmentDatabaseLigsSmokeTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            ligand_perception.require_template_store()
        except Exception as exc:                      # pragma: no cover
            raise unittest.SkipTest(f'CCD template store unavailable: {exc}')

    def test_end_to_end_roster_build(self):
        atp, adp, bgc = _heavy_names('ATP'), _heavy_names('ADP'), _heavy_names('BGC')
        types = [('ATP', atp), ('ADP', adp), ('BGC', bgc)]
        biounits = {
            # ATP and ADP together: a summing join would double their shared
            # fragments here.
            '1aaa': [0, 1],
            # The same type three times (NCS): must not multiply.
            '2bbb': [0, 0, 0],
            '3ccc': [2],
        }
        with tempfile.TemporaryDirectory() as tmp:
            roster = os.path.join(tmp, 'roster.pkl')
            _write_roster(roster, types, biounits)
            proc, payload, log = _run(tmp, roster)

        self.assertEqual(proc.returncode, 0, proc.stderr[-3000:])
        self.assertIsNotNone(payload)
        self.assertEqual(payload['key_schema'], Frags.KEY_SCHEMA)
        self.assertEqual(payload['db_identity']['sha256'], 'smoke')
        self.assertEqual(payload['ccd_identity'], ccd_templates.store_identity())
        self.assertEqual(payload['params']['support_unit'], 'distinct parent biounits')
        self.assertIn('support', payload)
        self.assertIn('support_pooled', payload)
        self.assertIn('Support frontier', log)

        support = payload['support']
        self.assertTrue(support)
        # No fragment can exceed the number of biounits, which is the invariant a
        # count in any of the three wrong units would break: ATP appears in two
        # biounits and three instances, so a copy-counting join would report 4.
        self.assertLessEqual(max(support.values()), len(biounits))
        # And the bound must actually be approached, or the assertion above holds
        # trivially for an implementation that counts nothing.
        self.assertEqual(max(support.values()), 2)

    def test_5_atom_arc_never_closes_but_a_real_5_ring_does(self):
        """BGC's 6-ring cut to 5 atoms is an arc; ADP's ribofuranose is a true 5-ring."""
        from rdkit import Chem
        import re
        # A ring-closure digit sits OUTSIDE a bracket; `D1`/`r6` digits sit inside,
        # so a plain substring test would fire on every annotated key.
        closure = re.compile(r'\](%?\d)')
        bgc = ligand_perception.perceive_ligand_graph('BGC').mol
        Chem.SanitizeMol(bgc)
        bgc_ring_keys = [k for k in Frags.enumerate_induced_fragments(bgc, 5, 5) if 'r6' in k]
        self.assertTrue(bgc_ring_keys)
        for key in bgc_ring_keys:
            self.assertIsNone(closure.search(key),
                              f'{key} closes a ring, but a 5-atom cut of a 6-ring '
                              f'cannot -- max-frag-size 5 is settled (DR-5 addendum)')
        # Vacuity clause: a genuine 5-ring (not a hand-written string) must close.
        adp = ligand_perception.perceive_ligand_graph('ADP').mol
        Chem.SanitizeMol(adp)
        adp_ring_keys = [k for k in Frags.enumerate_induced_fragments(adp, 5, 5) if 'r5' in k]
        self.assertTrue(adp_ring_keys)
        self.assertTrue(any(closure.search(k) for k in adp_ring_keys),
                         'no r5 key closes a ring for a genuine 5-ring at max frag size 5')

    def test_unobserved_atom_costs_support(self):
        """Dropping one observed atom must lower support for fragments touching it."""
        atp = _heavy_names('ATP')
        partial = tuple(n for n in atp if n != 'N1')
        with tempfile.TemporaryDirectory() as tmp:
            full_roster = os.path.join(tmp, 'full.pkl')
            _write_roster(full_roster, [('ATP', atp)], {'1aaa': [0]})
            _, full, _ = _run(tmp, full_roster)
            os.remove(os.path.join(tmp, 'database_frags_dict.pkl'))
            part_roster = os.path.join(tmp, 'part.pkl')
            _write_roster(part_roster, [('ATP', partial)], {'1aaa': [0]})
            _, part, _ = _run(tmp, part_roster)

        self.assertIsNotNone(full)
        self.assertIsNotNone(part)
        lost = set(full['support']) - set(part['support'])
        self.assertTrue(lost, 'removing an observed atom cost no fragment its support')
        # Not everything is lost: fragments away from N1 must survive, or this
        # would also pass for a build that simply dropped the ligand.
        self.assertTrue(set(part['support']))

    def test_pooling_agrees_between_write_time_and_read_time(self):
        """Every representative the reader picks must have a `support_pooled` entry.

        The enumerator pools protonation variants over the keys it WRITES; the
        reader pools over the keys it FINDS. If the record floor is applied after
        pooling, a group whose neutral twin was dropped gets a different (promoted)
        representative at read time and `select_fragments` raises "support is
        missing", aborting the build. The floor is set above 1 here on purpose --
        at the default it drops only the support-0 keys and the mismatch is rarer.
        """
        from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
            prepare_fragments, select_fragments)
        atp, adp, bgc = _heavy_names('ATP'), _heavy_names('ADP'), _heavy_names('BGC')
        types = [('ATP', atp), ('ADP', adp), ('BGC', bgc)]
        biounits = {'1aaa': [0, 1], '2bbb': [0], '3ccc': [2], '4ddd': [1, 2]}
        with tempfile.TemporaryDirectory() as tmp:
            roster = os.path.join(tmp, 'roster.pkl')
            _write_roster(roster, types, biounits)
            proc, payload, _ = _run(tmp, roster, extra=('--record-min-support', '3'))
        self.assertEqual(proc.returncode, 0, proc.stderr[-3000:])
        reps, _aliases = prepare_fragments(payload['frags'], 5)
        self.assertTrue(reps)
        missing = [r for r in reps if r not in payload['support_pooled']]
        self.assertEqual(missing, [], f'{len(missing)} representative(s) have no '
                                      f'pooled support entry, e.g. {missing[:2]}')
        # And the selector itself must run, which is what actually aborts a build.
        selected = select_fragments(payload['frags'], 3, 5,
                                    support=payload['support_pooled'])
        self.assertTrue(selected)
        # Vacuity clause: the floor must have dropped something, or write-time and
        # read-time saw the same set trivially.
        self.assertGreater(payload['stats']['num_frags_below_record_floor'], 0)

    def test_existing_output_is_refused_without_force(self):
        atp = _heavy_names('ATP')
        with tempfile.TemporaryDirectory() as tmp:
            roster = os.path.join(tmp, 'roster.pkl')
            _write_roster(roster, [('ATP', atp)], {'1aaa': [0]})
            first, payload, _ = _run(tmp, roster)
            self.assertEqual(first.returncode, 0, first.stderr[-2000:])
            second, _, _ = _run(tmp, roster)
            self.assertNotEqual(second.returncode, 0)
            self.assertIn('--force', second.stderr + second.stdout)
            third, _, _ = _run(tmp, roster, extra=('--force',))
            self.assertEqual(third.returncode, 0, third.stderr[-2000:])

    def test_a_partial_roster_is_refused(self):
        atp = _heavy_names('ATP')
        with tempfile.TemporaryDirectory() as tmp:
            roster = os.path.join(tmp, 'roster.pkl')
            _write_roster(roster, [('ATP', atp)], {'1aaa': [0]})
            with open(roster, 'rb') as fh:
                payload = pickle.load(fh)
            payload['partial'] = True
            with open(roster, 'wb') as fh:
                pickle.dump(payload, fh)
            proc, _, _ = _run(tmp, roster)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn('--limit', proc.stderr + proc.stdout)

    def test_a_roster_from_a_different_ccd_store_is_refused(self):
        atp = _heavy_names('ATP')
        with tempfile.TemporaryDirectory() as tmp:
            roster = os.path.join(tmp, 'roster.pkl')
            _write_roster(roster, [('ATP', atp)], {'1aaa': [0]})
            with open(roster, 'rb') as fh:
                payload = pickle.load(fh)
            payload['ccd_identity']['template_db_mtime_ns'] -= 1
            with open(roster, 'wb') as fh:
                pickle.dump(payload, fh)
            proc, output, _ = _run(tmp, roster)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIsNone(output)
        self.assertIn('different CCD template store', proc.stderr + proc.stdout)

if __name__ == '__main__':
    unittest.main()
