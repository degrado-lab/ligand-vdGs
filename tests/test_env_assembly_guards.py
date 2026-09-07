"""The two environment-assembly guards in _get_atomgroup_for_env.

Both exist because a vdG can come out looking well formed while being wrong, so
nothing downstream can detect the failure.

1. A CG atom name matching two atoms in one residue. Some prepared parent PDBs
   carry a residue with duplicated heavy-atom names (2y1x SAH A:1001 has two named
   N, 51 atoms against 46 in its other three copies, the spurious N 53 A from CA).
   The spurious atom is really THR D:478's backbone N: prepwizard cannot build an
   amino acid modelled with an N and no CA, so it re-emits the orphan under a
   ligand's resname/chain/resnum. _prep_filters.drop_prepwizard_hazard_residues (run from
   both s01 and s02) now removes the cause, but prepwizard also renames HET residues to a standard
   amino acid (CYT -> CYS in 5buv), which no input filter can prevent, so this guard
   is still load-bearing.
   Disambiguating by proximity -- what fingerprint_helpers._pick_atom_by_com did --
   can pick the wrong atom and write coordinates that contradict the record's own
   recorded identity.

2. A vdM that contacts the ligand but not the CG. The environment comes from probe
   contacts against the whole ligand residue, so without a check a residue touching
   one end of a large cofactor becomes a vdM of a CG matched at the other end.

   The heavy-atom side of that check reads the element *and* the atom name, because
   a blank element column (common in older and hand-edited PDBs) makes a plain
   `not element H D` selection keep every hydrogen, and an H-only near approach then
   passes a cutoff meant for heavy atoms.

3. A CG whose atoms are not one bonded component. The atom set comes from
   OpenBabel's perception of the parent, which nothing else checks against geometry;
   a SMARTS match is a connected subgraph, so atoms tens of A apart mean the match
   was resolved onto the wrong atoms.
"""
import os
import tempfile
import unittest

import numpy as np
import prody as pr

from ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs import _get_atomgroup_for_env

BIOUNIT = "1abc"
CG_NAME = "_test_cg_"          # must not be a key in vdg_miner constants.cg_atoms
CG_ATOM_NAMES = ["C1", "C2", "C3"]


def _atomgroup(duplicate_name=None, vdm_offset=3.0, tail_atom=False,
               split_cg=False, blank_element_h=False, second_vdm=False):
    """Ligand LIG A:900 (3 CG atoms) beside protein ALA A:10, close enough to pair.

    With duplicate_name set, a fourth ligand atom is added carrying that name, so
    the name matches two atoms -- the malformed case.

    vdm_offset slides ALA away along x. tail_atom adds a non-CG ligand atom near
    ALA, so the ligand *residue* still falls inside the 5 A environment selection
    while the CG atoms themselves do not contact it -- the large-cofactor case
    (a residue touching one end of FAD, a CG matched at the other).

    split_cg moves the third CG atom out of bonding range of the other two, the
    shape a CG atom resolved onto the wrong atom produces. blank_element_h gives
    ALA a hydrogen close to the CG whose element column is empty, so a
    name-blind heavy-atom filter would count it as a contact.
    """
    # Kept off the origin: align_coords_sanity_check rejects a zero-norm row, so a
    # CG atom at (0, 0, 0) would fail the frame check for reasons unrelated to this
    # test.
    names = list(CG_ATOM_NAMES)
    coords = [(10.0, 10.0, 10.0), (11.5, 10.0, 10.0), (10.75, 11.3, 10.0)]
    if split_cg:
        coords[2] = (10.75, 25.0, 10.0)
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
    names += ["N", "CA", "C"]
    coords += [(x, 10.5, 10.0), (x + 1.0, 11.2, 10.0), (x + 2.0, 10.8, 10.0)]
    resnames += ["ALA"] * 3; resnums += [10] * 3; elements += ["N", "C", "C"]

    if blank_element_h:
        # 2 A from CG atom C1, well inside the 4.5 A cutoff, but a hydrogen.
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


def _run(**kwargs):
    second = kwargs.get("second_vdm", False)
    environment = [(BIOUNIT, "", "A", 900, 1), (BIOUNIT, "", "A", 10, 1)]
    if second:
        environment.append((BIOUNIT, "", "A", 11, 1))
    pdb_dir = tempfile.gettempdir()
    pdb_file = os.path.join(pdb_dir, BIOUNIT[1:3].lower(), BIOUNIT + ".pdb")
    cg_match_dict = {(BIOUNIT, "", "A", "900", "LIG"): [CG_ATOM_NAMES]}
    with tempfile.NamedTemporaryFile("w", suffix=".log", delete=False) as fh:
        logfile = fh.name
    try:
        return _get_atomgroup_for_env(
            environment, pdb_dir, CG_NAME, cg_match_dict,
            align_atoms=[0, 1, 2], logfile=logfile,
            pdb_cache={pdb_file: _atomgroup(**kwargs)})
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
    """The mined CG atoms must form one bonded component."""

    def test_bonded_cg_is_accepted(self):
        self.assertIsNotNone(_run())

    def test_disconnected_cg_is_rejected(self):
        # C3 sits 13 A from the other two: a plausible SMARTS match cannot look
        # like this, so the atom set came from a bad perception/name resolution.
        self.assertIsNone(_run(split_cg=True))


class CgContactCutoffTests(unittest.TestCase):
    """A vdM must contact the CG, not merely the ligand it was cut from."""

    def test_contacting_vdm_is_accepted(self):
        ag = _run(vdm_offset=3.0)
        self.assertIsNotNone(ag)
        self.assertEqual(_vdm_resnums(ag), {10})

    def test_vdm_far_from_cg_is_not_a_slot(self):
        # The ligand still reaches ALA through its C9 tail, so the residue passes
        # the 5 A environment selection, but no CG atom is within the cutoff --
        # the large-cofactor case this guard exists for. The residue is dropped
        # from the vdG, which for this single-vdM environment leaves no slots.
        ag = _run(vdm_offset=8.0, tail_atom=True)
        self.assertEqual(_vdm_resnums(ag) if ag is not None else set(), set())

    def test_hydrogen_with_blank_element_does_not_count_as_contact(self):
        # Same rejected geometry, plus an ALA hydrogen 2 A from the CG carrying an
        # empty element column. `not element H D` keeps that atom, which would turn
        # a non-contact into a 2 A "contact" and admit the slot.
        ag = _run(vdm_offset=8.0, tail_atom=True, blank_element_h=True)
        self.assertEqual(_vdm_resnums(ag) if ag is not None else set(), set())

    def test_far_vdm_does_not_take_a_contacting_one_with_it(self):
        # The reason the cutoff drops a residue rather than the environment.
        # Subset sizes 1 and 2 are materialized from one shared reconstruction,
        # so discarding the environment here would destroy GLY 11's size-1 vdG
        # as well -- a slot that is genuinely in contact and belongs in the
        # library.
        ag = _run(vdm_offset=8.0, tail_atom=True, second_vdm=True)
        self.assertIsNotNone(ag)
        self.assertEqual(_vdm_resnums(ag), {11})


if __name__ == "__main__":
    unittest.main()
