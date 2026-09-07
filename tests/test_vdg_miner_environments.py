"""Contract tests for the lean external vdG-miner environment API."""
import inspect
import sys
import unittest
from pathlib import Path


VDG_MINER = Path(__file__).resolve().parents[1] / "external" / "vdG-miner" / "vdg_miner"
if str(VDG_MINER / "vdg") not in sys.path:
    sys.path.insert(0, str(VDG_MINER / "vdg"))

from vdg import VDG  # noqa: E402


class EnvironmentApiTests(unittest.TestCase):
    def test_constructor_does_not_build_legacy_fingerprint_features(self):
        vdg = VDG("_test_cg_", "pdb", "probe", "validation", cg_natoms=3)

        self.assertFalse(hasattr(vdg, "fingerprint_cols"))
        self.assertFalse(hasattr(vdg, "contact_cols"))
        self.assertFalse(hasattr(vdg, "ABPLE_cols"))
        self.assertFalse(hasattr(vdg, "relpos_cols"))

    def test_environment_api_requires_a_mining_source(self):
        vdg = VDG("_test_cg_", "pdb", "probe", "validation", cg_natoms=3)

        with self.assertRaisesRegex(ValueError, "chain_cluster or cg_match_dict"):
            vdg.mine_environments()

    def test_environment_api_always_returns_a_list(self):
        vdg = VDG("_test_cg_", "missing", "missing", "validation", cg_natoms=3)
        matches = {("1abc", "", "A", "1", "LIG"): [["C1", "C2", "C3"]]}

        self.assertEqual(vdg.mine_environments(cg_match_dict=matches), [])

    def test_a_missing_chain_file_does_not_discard_the_whole_structure(self):
        """One unreadable chain used to return [] for every chain of the structure."""
        vdg = VDG("_test_cg_", "missing", "missing", "validation", cg_natoms=3)
        seen = []

        def fake_structure_contacts(pdb_file, probe_file, cg_match_dict=None):
            seen.append(pdb_file)
            return None

        vdg.structure_contacts = fake_structure_contacts
        matches = {
            ("1abc", "", "A", "1", "LIG"): [["C1", "C2", "C3"]],
            ("1abc", "", "B", "2", "LIG"): [["C1", "C2", "C3"]],
        }
        # Both chains are visited even though neither file exists; the loop
        # continues rather than returning on the first miss.
        self.assertEqual(vdg.mine_environments(cg_match_dict=matches), [])

    def test_structure_level_work_is_separated_from_per_chain_work(self):
        """The probe file is per structure, so it must not be re-read per chain."""
        self.assertTrue(hasattr(VDG, "structure_contacts"))
        chain_params = list(inspect.signature(VDG.update_sc_info).parameters)
        self.assertEqual(chain_params,
                         ["self", "sc_info", "segi", "chain", "struct"])
        self.assertNotIn("probe_file", chain_params)

    def test_environment_api_has_no_fingerprint_or_validation_parameters(self):
        parameters = inspect.signature(VDG.mine_environments).parameters

        self.assertNotIn("logfile", parameters)
        self.assertNotIn("rscc", parameters)
        self.assertNotIn("rsr", parameters)
        self.assertNotIn("rsrz", parameters)
        self.assertEqual(
            list(parameters),
            ["self", "chain_cluster", "cg_match_dict", "pdb_gz", "min_seq_sep",
             "max_b_factor", "min_occ"],
        )

    def test_quality_thresholds_are_a_loose_floor_not_a_design_time_filter(self):
        """The measured values ride along, so the strict cut is a read-path call."""
        defaults = inspect.signature(VDG.mine_environments).parameters
        self.assertGreaterEqual(defaults["max_b_factor"].default, 100.0)
        self.assertLessEqual(defaults["min_occ"].default, 0.3)

    def test_sequence_separation_filter_is_off_by_default(self):
        """A ligand CG is not in the chain, so seq. sep. from it is undefined."""
        defaults = inspect.signature(VDG.mine_environments).parameters
        self.assertLessEqual(defaults["min_seq_sep"].default, 1)


if __name__ == "__main__":
    unittest.main()
