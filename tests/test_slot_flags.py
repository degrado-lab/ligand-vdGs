"""Slot labels and provenance flags emitted by reorder_vdg_subset.

Two things are under test and they carry different information on purpose:

  * the *label* says which moiety contacts the CG, and is the matchability key
    (bucket file name, and what aa_perm_indices may permute). GLY and PRO
    backbones get their own labels because they are not interchangeable with an
    ordinary backbone -- GLY has no CB to clash, PRO has no amide N-H to donate.
  * the *flag* says what the label does not: which moiety of the canonical
    residue is nearest the CG, plus whether the residue carries non-canonical
    atoms at all. On an 'X' slot that is the only record of what the real
    residue was doing.
"""

import unittest

import numpy as np
from prody import AtomGroup

from ligand_vdgs.functions.align_and_cluster import reorder_vdg_subset
from ligand_vdgs.functions import vdg_struct_utils as su


# A CG sitting on the +x axis, well clear of the backbone built below.
CG_COORDS = np.array([[6.0, 0.0, 0.0], [7.2, 0.0, 0.0]], dtype=np.float32)

# N, CA, C, O of a residue placed near the origin: ~4 A from the CG, so a
# backbone contact but no sidechain contact unless a sidechain reaches out.
_BACKBONE = {
    'N':  (0.0, 1.4, 0.0),
    'CA': (1.5, 1.4, 0.0),
    'C':  (2.2, 0.1, 0.0),
    'O':  (1.6, -1.0, 0.0),
}

# A sidechain reaching toward the CG, and one pointing away from it.
_REACHING = [('CB', (3.0, 1.0, 0.0)), ('CG', (4.4, 0.6, 0.0))]
_REMOTE = [('CB', (-2.0, 3.0, 0.0)), ('CG', (-3.2, 3.4, 0.0))]


def build_residue(resname, extra_atoms=()):
    """One-residue AtomGroup: canonical backbone plus whatever `extra_atoms` adds.

    extra_atoms: iterable of (atom_name, (x, y, z)).
    """
    names, coords = [], []
    for name, xyz in _BACKBONE.items():
        names.append(name)
        coords.append(xyz)
    for name, xyz in extra_atoms:
        names.append(name)
        coords.append(xyz)

    ag = AtomGroup('test')
    ag.setCoords(np.array(coords, dtype=np.float32))
    ag.setNames(names)
    ag.setResnames([resname] * len(names))
    ag.setResnums([1] * len(names))
    ag.setChids(['A'] * len(names))
    ag.setSegnames([''] * len(names))
    ag.setElements([n[0] for n in names])
    ag.setOccupancies([1.0] * len(names))
    return ag


def classify(resname, extra_atoms=()):
    """Return (label, slot_flag) for a single-residue vdG against CG_COORDS."""
    ag = build_residue(resname, extra_atoms)
    bb = np.array([_BACKBONE['N'], _BACKBONE['CA'], _BACKBONE['C']],
                  dtype=np.float32)
    scrr = ('', 'A', 1, resname)
    o = np.array(_BACKBONE['O'], dtype=np.float32)
    vdms_dict = {0: [resname, [scrr, bb, {0: ['vdm', bb[1]]}, o]]}
    aas, _bbs, _seqs, _cas, _scrr, flags, _os = reorder_vdg_subset(
        [0], vdms_dict, CG_COORDS, ag)
    return aas[0], flags[0]


class BackboneLabelTests(unittest.TestCase):
    def test_canonical_gly_gets_the_backbone_label(self):
        # No sidechain by chemistry, so the backbone is necessarily the contact.
        # Whether a *query* residue can host the result is decided on the read
        # path against its real atoms, not by giving glycine its own label.
        self.assertEqual(classify('GLY'), (su.BB_LABEL, su.SLOT_NO_SC))

    def test_proline_backbone_contact_gets_the_backbone_label(self):
        # Sidechain resolved but pointing away, so the backbone is the contact.
        # PRO's inability to donate an N-H is a read-path test, not a bucket.
        label, flag = classify('PRO', _REMOTE)
        self.assertEqual(label, su.BB_LABEL)
        self.assertEqual(flag, su.SLOT_BB_CLOSER)

    def test_ordinary_backbone_contact_is_plain_bb(self):
        label, flag = classify('LEU', _REMOTE)
        self.assertEqual(label, su.BB_LABEL)
        self.assertEqual(flag, su.SLOT_BB_CLOSER)

    def test_unresolved_sidechain_is_bb_with_no_sc_flag(self):
        # A LYS with only backbone atoms modeled. Label and flag are both
        # identical to glycine's, so neither separates "no sidechain by
        # chemistry" from "sidechain never resolved" -- only the stored
        # nr_scrr_resname does. Asserted here so the conflation stays a
        # known property rather than a surprise at analysis time.
        label, flag = classify('LYS')
        self.assertEqual(label, su.BB_LABEL)
        self.assertEqual(flag, su.SLOT_NO_SC)
        self.assertEqual((label, flag), classify('GLY'))

    def test_sidechain_contact_keeps_resname(self):
        label, flag = classify('LEU', _REACHING)
        self.assertEqual(label, 'LEU')
        self.assertEqual(flag, su.SLOT_SC)

    def test_terminal_oxygen_is_not_a_modification(self):
        label, flag = classify('GLY', [('OXT', (2.9, 0.0, 0.0))])
        self.assertEqual(label, su.BB_LABEL)
        self.assertEqual(flag, su.SLOT_NO_SC)


class ModifiedResidueTests(unittest.TestCase):
    def test_chromophore_atoms_touching_cg_are_not_a_glycine(self):
        # GFP-style: a residue named GLY carrying fused ring atoms, and those
        # atoms are what contacts the CG. Neither 'GLY' nor a backbone label
        # describes it. The flag still reports the real glycine underneath.
        label, flag = classify('GLY', [('CA2', (4.0, 0.4, 0.0)),
                                       ('CB2', (5.2, 0.2, 0.0))])
        self.assertEqual(label, su.NONCANONICAL_AA_LABEL)
        self.assertEqual(su.slot_reason(flag), su.SLOT_NO_SC)
        self.assertTrue(su.slot_is_modified(flag))

    def test_chromophore_atoms_remote_from_cg_stay_matchable(self):
        # Same modification, pointing away: the backbone is what contacts the
        # CG, so the geometry is reproducible and keeps a backbone label --
        # flagged modified so it can still be excluded downstream.
        label, flag = classify('GLY', [('CA2', (-2.0, 3.0, 0.0)),
                                       ('CB2', (-3.2, 3.4, 0.0))])
        self.assertEqual(label, su.BB_LABEL)
        self.assertEqual(su.slot_reason(flag), su.SLOT_NO_SC)
        self.assertTrue(su.slot_is_modified(flag))

    def test_modified_sidechain_remote_from_cg_keeps_resname(self):
        # An alkylated LYS whose canonical sidechain reaches the CG while the
        # added atoms do not: still an observation of a lysine sidechain.
        label, flag = classify('LYS', _REACHING + [('CX', (-4.0, 4.0, 0.0))])
        self.assertEqual(label, 'LYS')
        self.assertEqual(su.slot_reason(flag), su.SLOT_SC)
        self.assertTrue(su.slot_is_modified(flag))

    def test_selenomethionine_is_a_methionine(self):
        # Se-Met renamed to MET by structure prep, keeping SE in place of SD.
        # A phasing reagent, not a functional modification: it must not be
        # thrown into 'X' -- it was ~51% of every X slot before this case
        # existed.
        label, flag = classify('MET', [('CB', (3.0, 1.0, 0.0)),
                                       ('CG', (4.0, 0.8, 0.0)),
                                       ('SE', (4.9, 0.4, 0.0)),
                                       ('CE', (5.4, 1.8, 0.0))])
        self.assertEqual(label, 'MET')
        self.assertEqual(flag, su.SLOT_SC)
        self.assertFalse(su.slot_is_modified(flag))

    def test_nonstandard_resname_has_no_reference_atom_set(self):
        # MSE is not one of the 20, so there is nothing to compare its atoms to;
        # it is classified by geometry alone and keeps its own resname.
        label, flag = classify('MSE', [('CB', (3.0, 1.0, 0.0)),
                                       ('SE', (4.4, 0.6, 0.0))])
        self.assertEqual(label, 'MSE')
        self.assertEqual(flag, su.SLOT_SC)


class SplitResidueHeavyAtomTests(unittest.TestCase):
    def test_returns_none_for_unknown_resname(self):
        ag = build_residue('MSE', [('SE', (3.0, 1.0, 0.0))])
        self.assertEqual(su.split_residue_heavy_atoms(ag, 'MSE'),
                         (None, None, None))

    def test_canonical_atom_set_has_no_extras(self):
        ag = build_residue('GLY')
        bb, sc, extra = su.split_residue_heavy_atoms(ag, 'GLY')
        self.assertEqual(bb.shape, (4, 3))
        self.assertEqual(sc.shape, (0, 3))
        self.assertEqual(extra.shape, (0, 3))

    def test_terminal_oxygen_is_neither_backbone_nor_sidechain(self):
        # Matches ProDy's own split, so switching to name-based classification
        # does not quietly move OXT into the sidechain.
        ag = build_residue('GLY', [('OXT', (2.9, 0.0, 0.0))])
        bb, sc, extra = su.split_residue_heavy_atoms(ag, 'GLY')
        self.assertEqual(bb.shape, (4, 3))
        self.assertEqual(sc.shape, (0, 3))
        self.assertEqual(extra.shape, (0, 3))

    def test_hydrogens_are_never_reported_as_extras(self):
        ag = build_residue('GLY', [('HA2', (2.0, 2.4, 0.0))])
        _bb, _sc, extra = su.split_residue_heavy_atoms(ag, 'GLY')
        self.assertEqual(extra.shape, (0, 3))

    def test_selenium_counts_as_the_met_sidechain(self):
        ag = build_residue('MET', [('CB', (3.0, 1.0, 0.0)),
                                   ('SE', (4.9, 0.4, 0.0))])
        _bb, sc, extra = su.split_residue_heavy_atoms(ag, 'MET')
        self.assertEqual(sc.shape, (2, 3))
        self.assertEqual(extra.shape, (0, 3))


class QuerySlotLabelTests(unittest.TestCase):
    """The read-path counterpart: which labels a query residue may match under."""

    def test_every_residue_may_match_the_one_backbone_label(self):
        # Backbone geometries are pooled; whether a given one is hostable is a
        # per-vdG question answered by hit_finder_core.backbone_slots_can_host,
        # which sees the query residue's real sidechain.
        for resname in ('GLY', 'PRO', 'LEU', 'SER', 'TRP'):
            self.assertEqual(su.query_slot_labels(resname),
                             (resname, su.BB_LABEL))

    def test_every_residue_has_exactly_two_options(self):
        # The BSR enumeration is a product over these, so a third option for any
        # residue would multiply the query's bucket lookups.
        for resname in ('GLY', 'PRO', 'ALA', 'ARG', 'TRP'):
            self.assertEqual(len(su.query_slot_labels(resname)), 2)

    def test_query_never_produces_the_noncanonical_label(self):
        # 'X' buckets are deliberately unreachable: the geometry is real but its
        # chemical environment is not reproducible from a query structure.
        for resname in ('GLY', 'PRO', 'MET', 'CYS', 'LYS'):
            self.assertNotIn(su.NONCANONICAL_AA_LABEL,
                             su.query_slot_labels(resname))


class BackboneLabelHelperTests(unittest.TestCase):
    def test_label_per_residue(self):
        # One backbone label for every donor residue, glycine and proline included.
        for resname in ('GLY', 'PRO', 'LYS'):
            self.assertEqual(su.bb_label_for(resname), su.BB_LABEL)

    def test_no_label_contains_the_bucket_separator(self):
        # Bucket file names join sorted labels with '_', so a label carrying one
        # could not be split back apart.
        for label in su.BB_LABELS | {su.NONCANONICAL_AA_LABEL}:
            self.assertNotIn('_', label)

    def test_flag_fields_are_independent(self):
        for reason in (su.SLOT_SC, su.SLOT_NO_SC, su.SLOT_BB_CLOSER):
            self.assertEqual(su.slot_reason(reason), reason)
            self.assertFalse(su.slot_is_modified(reason))
            self.assertEqual(su.slot_reason(reason | su.SLOT_MODIFIED), reason)
            self.assertTrue(su.slot_is_modified(reason | su.SLOT_MODIFIED))


if __name__ == '__main__':
    unittest.main()
