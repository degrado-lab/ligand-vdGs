"""Contract tests for the lean external vdG-miner environment API."""
import inspect
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path


VDG_MINER = Path(__file__).resolve().parents[1] / "external" / "vdG-miner" / "vdg_miner"
if str(VDG_MINER / "vdg") not in sys.path:
    sys.path.insert(0, str(VDG_MINER / "vdg"))

from vdg import VDG  # noqa: E402
from ligand_vdgs.functions import sasa
from ligand_vdgs.functions.sasa import radius_of  # noqa: E402


class EnvironmentApiTests(unittest.TestCase):
    def test_constructor_does_not_build_legacy_fingerprint_features(self):
        vdg = VDG("_test_cg_", "pdb", "validation", cg_natoms=3)

        self.assertFalse(hasattr(vdg, "fingerprint_cols"))
        self.assertFalse(hasattr(vdg, "contact_cols"))
        self.assertFalse(hasattr(vdg, "ABPLE_cols"))
        self.assertFalse(hasattr(vdg, "relpos_cols"))

    def test_environment_api_requires_a_mining_source(self):
        vdg = VDG("_test_cg_", "pdb", "validation", cg_natoms=3)

        with self.assertRaisesRegex(ValueError, "chain_cluster or cg_match_dict"):
            vdg.mine_environments()

    def test_environment_api_always_returns_a_list(self):
        vdg = VDG("_test_cg_", "missing", "validation", cg_natoms=3)
        matches = {("1abc", "", "A", "1", "LIG"): [["C1", "C2", "C3"]]}

        self.assertEqual(vdg.mine_environments(cg_match_dict=matches), [])

    def test_a_missing_chain_file_does_not_discard_the_whole_structure(self):
        """One unreadable chain used to return [] for every chain of the structure."""
        vdg = VDG("_test_cg_", "missing", "validation", cg_natoms=3)
        seen = []

        def fake_structure_contacts(pdb_file, cg_match_dict=None):
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
        """The gate runs over the whole structure, not per chain."""
        self.assertTrue(hasattr(VDG, "structure_contacts"))
        chain_params = list(inspect.signature(VDG.update_sc_info).parameters)
        self.assertEqual(chain_params,
                         ["self", "sc_info", "segi", "chain", "struct"])
        self.assertNotIn("probe_file", chain_params)

    def test_nothing_in_the_api_takes_probe_output_any_more(self):
        self.assertNotIn("probe_file",
                         inspect.signature(VDG.structure_contacts).parameters)
        self.assertNotIn("probe_dir",
                         inspect.signature(VDG.__init__).parameters)
        vdg = VDG("_test_cg_", "pdb", "validation", cg_natoms=3)
        self.assertFalse(hasattr(vdg, "probe_dir"))

    def test_the_membership_threshold_is_a_recorded_parameter(self):
        """theta is calibration output, so it must not be a hidden literal."""
        vdg = VDG("_test_cg_", "pdb", "validation", cg_natoms=3)
        self.assertEqual(vdg.min_contact_area, sasa.MIN_CONTACT_AREA)
        self.assertEqual(
            VDG("_test_cg_", "pdb", "validation", cg_natoms=3,
                min_contact_area=2.5).min_contact_area, 2.5)

    def test_membership_thresholds_buried_plus_shared_not_buried_alone(self):
        """Criterion B (DR-3). Exclusive area alone has a 0.9485 recall ceiling.

        The discriminating case: a residue whose entire occluded patch is shared, so
        `buried_area` is 0 and only the shared term can admit it. A gate written as
        `buried_area > theta` drops it at every theta, which is exactly the 3,311
        Probe pairs the sweep found unreachable.
        """
        self.assertIs(VDG.structure_contacts.__globals__["sasa"].contact_area,
                      sasa.contact_area)
        c = sasa.ResidueContact(0.0, 0.0, 0, 1, 4.0, shared_area=9.0)
        self.assertEqual(sasa.contact_area(c), 9.0)
        self.assertGreater(sasa.contact_area(c), sasa.MIN_CONTACT_AREA,
                           "a fully shared patch of 9 A^2 must clear theta")
        # And the exclusive rule would still drop it, which is the point.
        self.assertEqual(c.buried_area, 0.0)

    def test_theta_is_meaningful_at_the_frozen_point_count(self):
        """The invariant, not the literal.

        theta is being re-picked (the pi test at 1,000 structures retired 5.0, and
        the whole-ligand sweep it came from is the wrong scale for a fragment), so
        pinning the number here would just be churn that hides the real constraint:
        theta has to be frozen against N_SPHERE_POINTS, and it has to be several
        area quanta or it is measuring discretisation rather than chemistry.
        """
        self.assertEqual(sasa.N_SPHERE_POINTS, 512,
                         "theta is quantised by the point count; changing one "
                         "without the other invalidates the calibration")
        # theta is 0.0: the pre-registered read of the mirror-scale divergence
        # test, where every band diverged so no flat region exists below the
        # lowest one. Pinning a literal here would be churn; what must not
        # silently change is the RATIONALE, so that is what is pinned.
        self.assertGreaterEqual(sasa.MIN_CONTACT_AREA, 0.0,
                                "a negative theta is not a threshold")
        if sasa.MIN_CONTACT_AREA > 0.0:
            # Any NON-zero theta is a claim that some positive area is too small
            # to be contact. Expressed as an area that claim is element-dependent:
            # the point area is 4*pi*(r+probe)^2/n_points, so a theta meant as a
            # one-point floor holds for C/N/O/F/Cl and fails for P, Br and I. A
            # theta between zero and the LARGEST element's point area admits
            # single-point contacts for the heavy halogens while excluding them
            # for carbon -- an inconsistency no calibration intended.
            largest = max(radius_of(e) for e in ('C', 'N', 'O', 'S', 'P', 'I'))
            point_area = (4 * 3.141592653589793 * (largest + sasa.PROBE_RADIUS) ** 2
                          / sasa.N_SPHERE_POINTS)
            self.assertGreaterEqual(
                sasa.MIN_CONTACT_AREA, point_area,
                "a positive theta below the largest element's point area "
                "({:.4f} A^2) excludes single-point contacts for light elements "
                "but admits them for P/Br/I; use an n_points floor if a "
                "point-count guarantee is what is wanted".format(point_area))

    def test_environment_api_has_no_fingerprint_or_validation_parameters(self):
        parameters = inspect.signature(VDG.mine_environments).parameters

        self.assertNotIn("logfile", parameters)
        self.assertNotIn("rscc", parameters)
        self.assertNotIn("rsr", parameters)
        self.assertNotIn("rsrz", parameters)
        self.assertEqual(
            list(parameters),
            ["self", "chain_cluster", "cg_match_dict", "pdb_gz", "min_seq_sep",
             "max_b_factor", "min_occ", "include_non_aa_partners"],
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


def _atom(serial, name, resname, chain, resnum, xyz, element, hetatm=False):
    record = "HETATM" if hetatm else "ATOM  "
    # Chain occupies columns 21-22 here on purpose: the database has two-character
    # chain IDs and the gate must not care.
    return ("{rec}{serial:>5} {name:<4}{alt:1}{resname:>3}{chain:>2}"
            "{resnum:>4}{icode:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
            "{occ:6.2f}{b:6.2f}          {el:>2}\n").format(
                rec=record, serial=serial, name=name, alt=" ", resname=resname,
                chain=chain, resnum=resnum, icode=" ", x=xyz[0], y=xyz[1],
                z=xyz[2], occ=1.0, b=20.0, el=element)


class GateBehaviourTests(unittest.TestCase):
    """Acceptance cases from contact-swap-plan.md, on hand-built structures.

    Written as whole PDB files and mined end to end, not as calls to the gate
    function, because the failures that matter here are in the wiring: which atoms
    are handed to the gate, which residues survive the AA filter, and whether the
    per-slot strengths line up with ``env[1:]``.
    """

    LIG_CHAIN = "A"

    def _mine(self, atoms, cg_names, stem="1tst", min_contact_area=0.0, **kwargs):
        """Mine a hand-built structure.

        `min_contact_area` defaults to 0 rather than the calibrated theta: these
        cases test the *wiring* -- which atoms reach the gate, which residues survive
        the filters, whether the per-slot lists line up -- and hand-placed atoms are
        not calibrated geometry. Tests that are about theta set it explicitly.
        """
        tmp = tempfile.mkdtemp()
        shard = os.path.join(tmp, stem[1:3].lower())
        os.makedirs(shard)
        with open(os.path.join(shard, stem + ".pdb"), "w") as handle:
            handle.writelines(atoms)
            handle.write("END\n")
        matches = {(stem, "", self.LIG_CHAIN, "900", "LIG"): [list(cg_names)]}
        vdg = VDG("_test_cg_", tmp, "validation", cg_natoms=len(cg_names),
                  min_contact_area=min_contact_area)
        try:
            return vdg.mine_environments(cg_match_dict=matches, **kwargs)
        finally:
            shutil.rmtree(tmp)

    def _ligand(self, extra=()):
        # Three-atom CG at the origin, plus whatever non-CG ligand atoms a case
        # needs. The ligand is HETATM resname LIG, chain A, resnum 900.
        atoms = [
            _atom(1, "C1", "LIG", self.LIG_CHAIN, 900, (0.0, 0.0, 0.0), "C", True),
            _atom(2, "O1", "LIG", self.LIG_CHAIN, 900, (1.4, 0.0, 0.0), "O", True),
            _atom(3, "N1", "LIG", self.LIG_CHAIN, 900, (-1.4, 0.0, 0.0), "N", True),
        ]
        atoms.extend(extra)
        return atoms

    def _residue(self, serial, resnum, xyz, resname="ALA", chain="A",
                 name="CB", element="C"):
        return _atom(serial, name, resname, chain, resnum, xyz, element)

    def test_a_residue_touching_the_cg_becomes_a_vdm_with_its_strengths(self):
        atoms = self._ligand()
        atoms.append(self._residue(10, 1, (0.0, 3.6, 0.0)))
        atoms.append(self._residue(11, 2, (1.4, -3.6, 0.0), resname="SER",
                                   name="OG", element="O"))
        envs = self._mine(atoms, ["C1", "O1", "N1"])
        self.assertEqual(len(envs), 1)
        env = envs[0]
        self.assertEqual(len(env["env"]) - 1, len(env["buried_area"]))
        self.assertEqual(len(env["buried_area"]), len(env["n_atom_pairs"]))
        self.assertEqual(len(env["buried_area"]), len(env["min_heavy_dist"]))
        self.assertTrue(all(a > 0 for a in env["buried_area"]))
        self.assertTrue(all(d < 6.5 for d in env["min_heavy_dist"]))

    def test_strengths_stay_aligned_when_the_neighbour_order_is_reversed(self):
        """env[1:] is sorted by resindex, so a lookup by key is the only safe one.

        Two residues at deliberately different distances; swapping which one is
        written first must swap the values with them, not leave them in file order.
        """
        near = (0.0, 3.4, 0.0)
        far = (0.0, -5.2, 0.0)
        first = self._ligand()
        first.append(self._residue(10, 1, near))
        first.append(self._residue(11, 2, far))
        second = self._ligand()
        second.append(self._residue(10, 1, far))
        second.append(self._residue(11, 2, near))
        a = self._mine(first, ["C1", "O1", "N1"])[0]
        b = self._mine(second, ["C1", "O1", "N1"])[0]
        self.assertAlmostEqual(a["min_heavy_dist"][0], b["min_heavy_dist"][1],
                               places=3)
        self.assertAlmostEqual(a["min_heavy_dist"][1], b["min_heavy_dist"][0],
                               places=3)
        self.assertGreater(a["buried_area"][0], a["buried_area"][1])
        self.assertGreater(b["buried_area"][1], b["buried_area"][0])

    def test_a_residue_reaching_only_a_non_cg_ligand_atom_is_dropped(self):
        """Case (b): 4.4 A to the ligand, 9 A to the CG. Buries no CG surface."""
        extra = [_atom(4, "C9", "LIG", self.LIG_CHAIN, 900, (9.0, 0.0, 0.0),
                       "C", True)]
        atoms = self._ligand(extra)
        atoms.append(self._residue(10, 1, (13.4, 0.0, 0.0)))
        atoms.append(self._residue(11, 2, (0.0, 3.6, 0.0)))
        envs = self._mine(atoms, ["C1", "O1", "N1"])
        self.assertEqual(len(envs), 1)
        resnums = [e[3] for e in envs[0]["env"][1:]]
        self.assertEqual(resnums, [2], "only the residue touching the CG survives")

    def test_a_residue_hidden_behind_another_is_theta_s_job_not_the_gate_s(self):
        """Case (b2), restated for criterion B.

        Under exclusive area a shadowed residue was dropped absolutely: it buries
        nothing, so no theta admitted it. That is precisely the blindness criterion B
        removes -- it is also why 3,311 real Probe contacts were unreachable -- so a
        shadowed residue now gets its 1/k share and the threshold decides. Both ends
        are asserted, because "dropped" and "kept" are each wrong on their own:
        residue 2 grazes the far side and earns 0.49 A^2, which theta = 5 rejects.
        """
        atoms = self._ligand()
        atoms.append(self._residue(10, 1, (0.0, 3.3, 0.0)))
        atoms.append(self._residue(11, 2, (0.0, 6.2, 0.0)))

        loose = self._mine(atoms, ["C1", "O1", "N1"], min_contact_area=0.0)
        resnums = [e[3] for e in loose[0]["env"][1:]]
        self.assertEqual(resnums, [1, 2])
        idx = resnums.index(2)
        self.assertEqual(loose[0]["buried_area"][idx], 0.0,
                         "residue 2 is fully behind residue 1: nothing is exclusive")
        self.assertGreater(loose[0]["shared_area"][idx], 0.0,
                           "but it does occlude, so criterion B can see it at all")

        # At the calibrated theta (0.0) this residue IS admitted: 0 is the
        # pre-registered read of the mirror-scale divergence test, where the
        # faintest band was divergent like every other. The environment records
        # buried_area 0 with a positive shared_area, so a reader that wants it gone
        # can filter on either at query time -- which is the argument for 0 being
        # the recoverable direction.
        calibrated = self._mine(atoms, ["C1", "O1", "N1"],
                                min_contact_area=sasa.MIN_CONTACT_AREA)
        self.assertIn(2, [e[3] for e in calibrated[0]["env"][1:]],
                      "theta 0 admits a contact whose patch is entirely shared")

        # And it is still theta's job, not the gate's: raise theta above this
        # residue's shared area and it goes. Both ends are asserted because either
        # alone would pass for the wrong reason -- "kept" could mean the gate
        # ignores theta, "dropped" could mean shared patches are invisible.
        shared = loose[0]["shared_area"][idx]
        tight = self._mine(atoms, ["C1", "O1", "N1"],
                           min_contact_area=shared * 1.01)
        self.assertEqual([e[3] for e in tight[0]["env"][1:]], [1],
                         "above its shared area the grazing shadow is not a contact")

    def test_a_two_character_chain_id_structure_still_yields_environments(self):
        """Case (e): the old probe/PDB hash dropped 670 of these outright.

        `cg.py` and ProDy both read the chain from column 22, so both see "4" and
        the ligand key still resolves. What used to kill these files was the probe
        line hash over columns 13-26, which included column 21 and so never matched
        the PDB line. Nothing reads probe output now, so the structure survives.
        The residue-level collision two-character chains also cause (see
        tests/test_chain_ids.py) is a separate defect the repair pass fixes; this
        only asserts the structure is not dropped wholesale.
        """
        atoms = [
            _atom(1, "C1", "LIG", "A4", 900, (0.0, 0.0, 0.0), "C", True),
            _atom(2, "O1", "LIG", "A4", 900, (1.4, 0.0, 0.0), "O", True),
            _atom(3, "N1", "LIG", "A4", 900, (-1.4, 0.0, 0.0), "N", True),
            _atom(10, "CB", "ALA", "A4", 1, (0.0, 3.6, 0.0), "C"),
            _atom(11, "OG", "SER", "A4", 2, (1.4, -3.6, 0.0), "O"),
        ]
        stem = "8tst"
        tmp = tempfile.mkdtemp()
        shard = os.path.join(tmp, stem[1:3].lower())
        os.makedirs(shard)
        with open(os.path.join(shard, stem + ".pdb"), "w") as handle:
            handle.writelines(atoms)
            handle.write("END\n")
        matches = {(stem, "", "4", "900", "LIG"): [["C1", "O1", "N1"]]}
        vdg = VDG("_test_cg_", tmp, "validation", cg_natoms=3)
        try:
            envs = vdg.mine_environments(cg_match_dict=matches)
        finally:
            shutil.rmtree(tmp)
        self.assertEqual(len(envs), 1)
        self.assertEqual(len(envs[0]["env"]) - 1, 2)

    def test_a_metal_is_measured_but_not_emitted_until_the_switch_is_flipped(self):
        """Case (c): the gate sees it; the writer contract decides whether it ships."""
        atoms = self._ligand()
        atoms.append(_atom(10, "ZN", "ZN", "A", 500, (0.0, 2.4, 0.0), "ZN", True))
        atoms.append(self._residue(11, 1, (0.0, -3.6, 0.0)))
        default = self._mine(atoms, ["C1", "O1", "N1"])
        self.assertEqual([e[3] for e in default[0]["env"][1:]], [1])
        opened = self._mine(atoms, ["C1", "O1", "N1"],
                            include_non_aa_partners=True)
        self.assertIn(500, [e[3] for e in opened[0]["env"][1:]])

    def test_a_water_with_only_one_leg_is_not_a_bridge(self):
        """Case (d): a water packed against the CG but donating to nothing."""
        atoms = self._ligand()
        atoms.append(_atom(10, "O", "HOH", "A", 700, (0.0, 3.0, 0.0), "O", True))
        atoms.append(self._residue(11, 1, (0.0, 9.0, 0.0)))
        envs = self._mine(atoms, ["C1", "O1", "N1"])
        self.assertEqual(envs, [], "a lone water is not an environment")

    def test_a_bridging_water_brings_its_partner_in(self):
        atoms = self._ligand()
        atoms.append(_atom(10, "O", "HOH", "A", 700, (-1.4, 2.8, 0.0), "O", True))
        atoms.append(self._residue(11, 1, (-1.4, 5.6, 0.0), resname="SER",
                                   name="OG", element="O"))
        atoms.append(self._residue(12, 2, (1.4, -3.6, 0.0), resname="SER",
                                   name="OG", element="O"))
        envs = self._mine(atoms, ["C1", "O1", "N1"])
        self.assertEqual(len(envs), 1)
        resnums = [e[3] for e in envs[0]["env"][1:]]
        self.assertIn(1, resnums, "the water-bridged residue must be admitted")
        idx = resnums.index(1)
        self.assertEqual(envs[0]["buried_area"][idx], 0.0,
                         "it buries no CG surface; only the water does")
        self.assertGreater(envs[0]["min_heavy_dist"][idx], 4.5)

    def test_the_water_itself_is_never_a_vdm(self):
        atoms = self._ligand()
        atoms.append(_atom(10, "O", "HOH", "A", 700, (-1.4, 2.8, 0.0), "O", True))
        atoms.append(self._residue(11, 1, (-1.4, 5.6, 0.0), resname="SER",
                                   name="OG", element="O"))
        envs = self._mine(atoms, ["C1", "O1", "N1"])
        for env in envs:
            self.assertNotIn(700, [e[3] for e in env["env"][1:]])


if __name__ == "__main__":
    unittest.main()
