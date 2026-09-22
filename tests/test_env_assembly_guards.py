"""Environment-assembly guards whose failures become silent downstream.

The cases cover duplicate CG atom names, vdMs with no heavy atoms, and CGs whose
OpenBabel-resolved atoms violate connectivity or the provisional SMARTS-edge
geometry envelope. Residue membership remains the upstream SASA gate's decision.
"""
import os
import tempfile
import unittest

import numpy as np
import prody as pr
from rdkit import Chem

from ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs import (
    _cg_bond_components, _cg_bond_length_violations, _cg_smarts_bond_class,
    _get_atomgroup_for_env)
from tests.vacuity import assert_discriminates

BIOUNIT = "1abc"
CG_NAME = "_test_cg_"          # must not be a key in vdg_miner constants.cg_atoms
CG_ATOM_NAMES = ["C1", "C2", "C3"]

def _atomgroup(duplicate_name=None, vdm_offset=3.0, tail_atom=False,
               split_cg=False, blank_element_h=False, second_vdm=False,
               only_blank_element_h=False, bridged_wrong_bond=False,
               collapsed_bond=False):
    """Build a synthetic three-atom ligand and nearby protein environment."""
    # Kept off the origin: align_coords_sanity_check rejects a zero-norm row, so a
    # CG atom at (0, 0, 0) would fail the frame check for reasons unrelated to this
    # test.
    names = list(CG_ATOM_NAMES)
    coords = [(10.0, 10.0, 10.0), (11.5, 10.0, 10.0), (10.75, 11.3, 10.0)]
    if split_cg:
        coords[2] = (10.75, 25.0, 10.0)
    if bridged_wrong_bond:
        # C3 is 3.8 A from its real SMARTS neighbour C1, but 2.3 A from C2.
        # The old all-pairs connectivity test therefore accepted this row.
        coords[2] = (13.8, 10.0, 10.0)
    if collapsed_bond:
        coords[1] = (10.3, 10.0, 10.0)
    resnames, resnums, elements = ["LIG"] * 3, [900] * 3, ["C"] * 3

    if duplicate_name is not None:
        names.append(duplicate_name)
        coords.append((30.0, 30.0, 30.0))   # far away, as in the real defect
        resnames.append("LIG"); resnums.append(900); elements.append("C")

    if tail_atom:
        names.append("C9")
        coords.append((10.0 + vdm_offset - 4.0, 10.0, 10.0))
        resnames.append("LIG"); resnums.append(900); elements.append("C")

    # Protein residue within 5 A of the ligand, so the environment selection keeps it.
    x = 10.0 + vdm_offset
    if not only_blank_element_h:
        names += ["N", "CA", "C"]
        coords += [(x, 10.5, 10.0), (x + 1.0, 11.2, 10.0), (x + 2.0, 10.8, 10.0)]
        resnames += ["ALA"] * 3; resnums += [10] * 3; elements += ["N", "C", "C"]

    if blank_element_h or only_blank_element_h:
        # 2 A from CG atom C1 -- close enough that a name-blind filter would
        # treat it as a heavy contact -- but a hydrogen.
        names.append("HB1")
        coords.append((10.0, 12.0, 10.0))
        resnames.append("ALA"); resnums.append(10); elements.append("")

    if second_vdm:
        # GLY 11, genuinely in contact with the CG. Present so the far residue's
        # rejection can be shown not to take a good slot down with it.
        names += ["N", "CA", "C"]
        coords += [(10.2, 12.8, 10.0), (11.0, 13.4, 10.0), (12.0, 13.0, 10.0)]
        resnames += ["GLY"] * 3; resnums += [11] * 3; elements += ["N", "C", "C"]

    ag = pr.AtomGroup("synthetic")
    ag.setCoords(np.array(coords, dtype=float))
    ag.setNames(names)
    ag.setResnames(resnames)
    ag.setResnums(np.array(resnums, dtype=int))
    ag.setChids(["A"] * len(names))
    ag.setSegnames([""] * len(names))
    ag.setElements(elements)
    ag.setOccupancies(np.ones(len(names)))
    ag.setAltlocs([" "] * len(names))
    return ag

def _vdm_resnums(atomgroup):
    """Resnums the assembled atomgroup actually offers as vdM slots.

    get_vdm_res_features selects vdMs by the VDM_OCC occupancy marker alone, so
    this is the same question it asks -- a residue left unmarked is not a slot,
    even though its atoms are still in the group.
    """
    sel = atomgroup.select("occupancy > 1.5 and occupancy < 2.5")
    return set() if sel is None else set(int(r) for r in sel.getResnums())

def _run(row_rejections=None, **kwargs):
    environment = [(BIOUNIT, "", "A", 900, 1), (BIOUNIT, "", "A", 10, 1)]
    if kwargs.get("second_vdm", False):
        environment.append((BIOUNIT, "", "A", 11, 1))
    pdb_dir = tempfile.gettempdir()
    with tempfile.NamedTemporaryFile("w", suffix=".log", delete=False) as fh:
        logfile = fh.name
    try:
        return _get_atomgroup_for_env(
            environment, pdb_dir, CG_NAME,
            {(BIOUNIT, "", "A", "900", "LIG"): [CG_ATOM_NAMES]},
            align_atoms=[0, 1, 2], logfile=logfile,
            pdb_cache={os.path.join(pdb_dir, BIOUNIT[1:3].lower(), BIOUNIT + ".pdb"):
                       _atomgroup(**kwargs)},
            cg_bonds=((0, 1, 'single_or_aromatic'),
                      (0, 2, 'single_or_aromatic')),
            expected_elements=('C', 'C', 'C'), row_rejections=row_rejections)
    finally:
        os.unlink(logfile)

class DuplicateAtomNameTests(unittest.TestCase):
    def test_wellformed_residue_is_accepted(self):
        # Guards against the rejection being so broad it drops ordinary residues.
        self.assertIsNotNone(_run())

    def test_duplicate_cg_atom_name_is_rejected(self):
        # C1 is a CG atom, so the ambiguity reaches the CG coordinates.
        self.assertIsNone(_run(duplicate_name="C1"))

    def test_duplicate_non_cg_atom_name_is_tolerated(self):
        # A duplicate elsewhere in the residue cannot corrupt the CG, so it is
        # not grounds to discard the environment.
        self.assertIsNotNone(_run(duplicate_name="C9"))

class CgBondPlausibilityTests(unittest.TestCase):
    """Every SMARTS edge must have chemically plausible deposited geometry."""

    def test_bonded_cg_is_accepted(self):
        self.assertIsNotNone(_run())

    def test_disconnected_cg_is_rejected(self):
        # C3 sits 13 A from the other two: a plausible SMARTS match cannot look
        # like this, so the atom set came from a bad perception/name resolution.
        self.assertIsNone(_run(split_cg=True))

    def test_wrong_geminal_edge_is_rejected(self):
        # Falsifier: the buggy component check accepts C1--C2--C3 through the
        # non-SMARTS C2--C3 proximity edge.
        self.assertEqual(_cg_bond_components(
            _atomgroup(bridged_wrong_bond=True).getCoords()[:3]), [[0, 1, 2]])
        self.assertIsNone(_run(bridged_wrong_bond=True))

    def test_collapsed_smarts_bond_is_rejected(self):
        self.assertEqual(_cg_bond_components(
            _atomgroup(collapsed_bond=True).getCoords()[:3]), [[0, 1, 2]])
        self.assertIsNone(_run(collapsed_bond=True))

    def test_radius_gate_is_non_vacuous_at_both_bounds(self):
        def accepted(distance):
            return not _cg_bond_length_violations(
                np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]]),
                ('C', 'C'), ((0, 1, 'single'),))
        assert_discriminates(
            accepted, [1.50], [0.30, 2.30], 'SMARTS C-C bond-length gate')

    def test_smarts_bond_class_is_preserved(self):
        def bond_class(smarts):
            mol = Chem.MolFromSmarts(smarts)
            self.assertIsNotNone(mol)
            self.assertEqual(mol.GetNumBonds(), 1)
            return _cg_smarts_bond_class(mol.GetBondWithIdx(0))

        self.assertEqual(bond_class('[C]-[N]'), 'single')
        self.assertEqual(bond_class('[C]=[N]'), 'double')
        self.assertEqual(bond_class('[C]#[N]'), 'triple')
        self.assertEqual(bond_class('[c]:[n]'), 'aromatic')
        self.assertEqual(bond_class('[C][N]'), 'single_or_aromatic')
        self.assertEqual(bond_class('[C]~[N]'), 'ambiguous')

    def test_ambiguous_bond_uses_conservative_fallback(self):
        def accepted(distance):
            return not _cg_bond_length_violations(
                np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]]),
                ('C', 'C'), ((0, 1, 'ambiguous'),))
        assert_discriminates(
            accepted, [1.50], [0.30, 2.30],
            'ambiguous SMARTS bond uses broad covalent envelope')

    def test_multiple_bad_edges_count_as_one_rejected_row(self):
        rejections = {}
        self.assertIsNone(_run(split_cg=True, collapsed_bond=True,
                               row_rejections=rejections))
        self.assertEqual(rejections, {'cg_bond_geometry': 1})

class NoWriterSideDistanceGateTests(unittest.TestCase):
    """Membership is the SASA gate's decision; the writer must not re-litigate it.

    The removed guard was a 4.5 A heavy-atom cutoff. Under buried SASA a residue
    can bury CG surface out to ~6.5 A, and real members sit in the 4.5-6.5 A band
    (Met-SD off a ring face, ~15% of members in a 12-structure sweep), so every
    case below that used to be dropped must now survive.
    """

    def test_contacting_vdm_is_accepted(self):
        ag = _run(vdm_offset=3.0)
        self.assertIsNotNone(ag)
        self.assertEqual(_vdm_resnums(ag), {10})

    def test_vdm_beyond_the_old_cutoff_is_kept(self):
        # The discriminating case: 5.5 A from the CG is past the deleted 4.5 A
        # cutoff and inside the SASA gate's reach. This assertion is the inverse
        # of the one it replaces.
        ag = _run(vdm_offset=5.5, tail_atom=True)
        self.assertIsNotNone(ag)
        self.assertEqual(_vdm_resnums(ag), {10})

    def test_vdm_far_from_the_cg_is_still_kept(self):
        # Even a residue the gate would never admit is kept here: the writer has
        # no contact information, and inventing one would re-create the guard.
        # Whether this residue is a member was settled upstream.
        ag = _run(vdm_offset=8.0, tail_atom=True)
        self.assertIsNotNone(ag)
        self.assertEqual(_vdm_resnums(ag), {10})

    def test_an_all_hydrogen_residue_is_still_dropped(self):
        # The one residue-level check that survives, and the reason it reads the
        # atom name as well as the element: with a blank element column a plain
        # `not element H D` keeps every hydrogen, so this residue would look as
        # though it had heavy atoms for the gate to have measured.
        ag = _run(vdm_offset=8.0, tail_atom=True, only_blank_element_h=True)
        self.assertEqual(_vdm_resnums(ag) if ag is not None else set(), set())

    def test_dropping_a_residue_does_not_drop_the_environment(self):
        # Subset sizes 1 and 2 are materialized from one shared reconstruction,
        # so discarding the environment over one bad residue would destroy
        # GLY 11's size-1 vdG as well.
        ag = _run(vdm_offset=8.0, tail_atom=True, second_vdm=True,
                  only_blank_element_h=True)
        self.assertIsNotNone(ag)
        self.assertEqual(_vdm_resnums(ag), {11})

if __name__ == "__main__":
    unittest.main()
