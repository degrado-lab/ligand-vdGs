"""The CCD, not perception, decides ligand chemistry.

Every case here is one the toolkits get wrong or answer differently, which is
the reason rebuild-notes 6 makes the template authoritative: bond order and
protonation from coordinates are guesses, and a wrong guess changes a
fragment key without raising anything.
"""
import os
import tempfile

import pytest
from openbabel import openbabel as ob
from rdkit import Chem, RDLogger

from ligand_vdgs.functions import ccd_templates, ligand_perception as lp

RDLogger.DisableLog('rdApp.*')

pytestmark = pytest.mark.skipif(
    not os.path.isfile(ccd_templates.template_db_path()),
    reason=f'{ccd_templates.template_db_path()} not built '
           '(scripts/build_ccd_templates.py)')

# Acetate, whose CCD template is the anion: C(=O)[O-] with a methyl.
# A three-letter code the CCD does not have. Three letters matters: the PDB
# resname column is 18-20, so a four-character code shifts every field after
# it and OpenBabel then parses a different molecule -- which silently makes a
# template-vs-perception comparison compare two different structures.
NO_TEMPLATE = 'ZQ8'

ACETATE = """\
HETATM    1  C   {r} A 900       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O   {r} A 900       1.250   0.000   0.000  1.00  0.00           O
HETATM    3  OXT {r} A 900      -0.700   1.210   0.000  1.00  0.00           O
HETATM    4  CH3 {r} A 900      -0.750  -1.290   0.000  1.00  0.00           C
END
"""
# Same ligand with the methyl carbon unmodelled -- routine partial density.
ACETATE_PARTIAL = "\n".join(
    line for line in ACETATE.splitlines() if ' CH3 ' not in line) + "\n"


def _matches(obmol, smarts):
    pattern = ob.OBSmartsPattern()
    assert pattern.Init(smarts), smarts
    return len(list(pattern.GetUMapList())) if pattern.Match(obmol) else 0


def _by_name(obmol):
    out = {}
    for i in range(1, obmol.NumAtoms() + 1):
        atom = obmol.GetAtom(i)
        residue = atom.GetResidue()
        name = residue.GetAtomID(atom).strip() if residue is not None else ''
        heavy = sum(1 for n in ob.OBAtomAtomIter(atom) if n.GetAtomicNum() != 1)
        out[name] = (heavy, atom.GetImplicitHCount(), atom.GetFormalCharge())
    return out


class ParserTests:
    pass


def test_monatomic_ions_are_parsed():
    # Every monatomic ion is written as plain `_chem_comp_atom.key value` items
    # rather than a loop; a loop-only parser drops all of them silently.
    zinc = ccd_templates.get_template('ZN')
    assert zinc is not None
    assert [(a.name, a.element, a.charge) for a in zinc.atoms] == [('ZN', 'Zn', 2)]


def test_hydrogen_counts_come_from_the_template_graph():
    counts = ccd_templates.template_h_counts(ccd_templates.get_template('TYR'))
    assert counts['CA'] == 1 and counts['CB'] == 2 and counts['CG'] == 0
    assert counts['OH'] == 1
    assert 'HA' not in counts        # hydrogens are not keys, only counts


def test_unknown_component_is_none_not_an_error():
    assert ccd_templates.get_template(NO_TEMPLATE) is None


class TemplateBeatsPerception:
    pass


def test_protonation_comes_from_the_template():
    # OpenBabel reads acetate's carboxylate as a neutral acid; the CCD says
    # anion. This is cg_num_h and cg_formal_charge for every carboxylate in
    # the library.
    templated = lp.perceive_ligand_instance(ACETATE.format(r='ACT'), 'ACT')
    perceived = lp.perceive_ligand_instance(ACETATE.format(r=NO_TEMPLATE), NO_TEMPLATE)
    assert templated.provenance == lp.PERCEPTION_CCD_TEMPLATE
    assert perceived.provenance == lp.PERCEPTION_OPENBABEL
    assert _by_name(templated.obmol)['OXT'] == (1, 0, -1)
    assert _by_name(perceived.obmol)['OXT'] == (1, 1, 0), 'case stopped discriminating'


def test_an_unobserved_neighbour_does_not_truncate_degree_or_the_h_flag():
    """The case perception gets wrong by construction (rebuild-notes 6).

    With the methyl unmodelled, the carboxyl carbon has two observed
    neighbours. OpenBabel therefore reads it as `D2` and hands it a hydrogen,
    so it keys as `!H0` -- the same fragment as a genuine CH. The template
    knows the third neighbour exists.
    """
    templated = lp.perceive_ligand_instance(ACETATE_PARTIAL.format(r='ACT'), 'ACT')
    perceived = lp.perceive_ligand_instance(ACETATE_PARTIAL.format(r=NO_TEMPLATE), NO_TEMPLATE)
    assert _matches(templated.obmol, '[C;D3]') == 1
    assert _matches(templated.obmol, '[C;H0]') == 1
    # Discriminating: without the template both answers flip.
    assert _matches(perceived.obmol, '[C;D3]') == 0
    assert _matches(perceived.obmol, '[C;H0]') == 0


def test_the_unobserved_atom_is_present_but_unnamed():
    # It has to be a real graph atom for `D<n>` to count it, and it has to be
    # nameless so find_cg_matches drops any match that reaches it -- which is
    # how "a fragment with an unobserved atom is skipped" is implemented.
    templated = lp.perceive_ligand_instance(ACETATE_PARTIAL.format(r='ACT'), 'ACT')
    named = lp.pdb_atom_names(templated.obmol)
    assert templated.obmol.NumAtoms() == 4
    assert set(named.values()) == {'C', 'O', 'OXT'}
    assert len(lp.phantom_atom_indices(templated.obmol)) == 1


def test_a_fully_observed_ligand_gains_no_phantoms():
    templated = lp.perceive_ligand_instance(ACETATE.format(r='ACT'), 'ACT')
    assert lp.phantom_atom_indices(templated.obmol) == set()


def test_an_observed_atom_with_no_template_name_falls_back():
    # The only real mapping failure. A template atom that is *absent* is not
    # one -- leaving atoms and partial density are both normal.
    bad = ACETATE.format(r='ACT').replace(' CH3 ', ' QQQ ')
    assert lp.perceive_ligand_instance(bad, 'ACT').provenance == \
        lp.PERCEPTION_OPENBABEL


def test_a_covalently_attached_ligand_keeps_its_template():
    # OXT/HXT are pdbx_leaving_atom_flag = Y and are absent by design once the
    # ligand is linked. Treating that as a mapping failure would send every
    # covalent adduct down the perception path.
    linked = "\n".join(line for line in ACETATE.format(r='ACT').splitlines()
                       if ' OXT ' not in line) + "\n"
    perceived = lp.perceive_ligand_instance(linked, 'ACT')
    assert perceived.provenance == lp.PERCEPTION_CCD_TEMPLATE


class GraphSide:
    pass


def test_graph_side_uses_the_template_and_reports_it():
    graph = lp.perceive_ligand_graph('ATP')
    assert graph.provenance == lp.PERCEPTION_CCD_TEMPLATE
    assert graph.mol.GetNumAtoms() == 31          # heavy atoms only
    # Adenine: a fused 5+6 sharing two atoms.
    assert sum(1 for a in graph.mol.GetAtoms() if a.GetIsAromatic()) == 9


def test_graph_side_falls_back_to_smiles_for_a_non_ccd_ligand():
    graph = lp.perceive_ligand_graph(NO_TEMPLATE, smiles='CCO')
    assert graph.provenance == lp.PERCEPTION_SMILES
    assert graph.mol.GetNumAtoms() == 3


def test_a_template_built_graph_round_trips_through_the_key_writer():
    """Keys from a template graph must match themselves.

    Aromaticity here comes from the CCD's flags rather than RDKit's model, so
    this is the check that the key writer accepts that graph at all -- a key
    that cannot match itself is the dedup failure mode.
    """
    from ligand_vdgs.functions import Frags
    from ligand_vdgs.functions.utils import fragment_keys_equivalent
    for comp_id in ('ATP', 'TYR', 'HEM'):
        mol = lp.perceive_ligand_graph(comp_id).mol
        stripped = Frags.manually_remove_Hs(mol, 'single')
        assert stripped is not None, comp_id
        keys = Frags.get_fragments(2, stripped[0][0], 4, 5)
        assert keys, comp_id
        for key in keys:
            assert fragment_keys_equivalent(key, key), (comp_id, key)


def test_aromatic_ring_atoms_pool_across_pyridine_and_pyrrole():
    # The vocabulary decision `ring_query_for_atom` documents: aromatic atoms
    # carry no ring-size primitive, so a 5- and a 6-ring aromatic carbon are
    # the same key atom. Template-built graphs must not break that.
    from ligand_vdgs.functions import Frags
    for comp_id in ('ATP', 'TYR'):
        mol = lp.perceive_ligand_graph(comp_id).mol
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                assert Frags.ring_query_for_atom(atom) is None


def _block_from_template(comp_id, drop=()):
    """A PDB block of a component's heavy atoms at their CCD ideal coordinates.

    Geometry only has to be good enough for OpenBabel to read the file; the
    template is what decides the chemistry, which is the point of the test.
    """
    template = ccd_templates.get_template(comp_id)
    lines = []
    serial = 0
    for atom in template.atoms:
        if atom.element in ('H', 'D') or atom.name in drop:
            continue
        serial += 1
        x, y, z = atom.xyz
        # Columns matter: read_ligand_blocks parses the text by column, and a
        # name or element in the wrong place silently yields no ligand at all.
        # Name is 13-16, left-justified only when it fills all four.
        name = atom.name if len(atom.name) >= 4 else f' {atom.name:<3s}'
        lines.append(
            f'HETATM{serial:5d} {name:<4s} {comp_id:>3s} A 900    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          '
            f'{atom.element.upper():>2s}')
    return '\n'.join(lines) + '\nEND\n'


@pytest.mark.parametrize('comp_id', ['ACT', 'TYR', 'ATP', 'HEM', 'NAD'])
def test_keys_written_from_the_graph_match_the_same_ligand_as_a_block(comp_id):
    """The cross-side agreement the whole template swap is for (verify item 5).

    Keys are written from the template graph (RDKit) and matched at mining
    time against an OBMol (OpenBabel). If the two disagree on one bond order
    or one aromatic flag, the affected keys mine zero sites and nothing
    raises. Aromatics are the case the acetate tests cannot reach: an aromatic
    key bond has to survive `apply_template_to_obmol` writing it and
    `SetAromaticPerceived` keeping OpenBabel from re-deriving it.
    """
    from ligand_vdgs.functions import Frags
    from ligand_vdgs.generate_vdgs.estimate_frag_cost import add_vdg_miner_paths
    add_vdg_miner_paths()
    import cg

    mol = lp.perceive_ligand_graph(comp_id).mol
    stripped = Frags.manually_remove_Hs(mol, 'single')
    assert stripped is not None, comp_id
    keys = sorted(Frags.get_fragments(2, stripped[0][0], 4, 5))
    assert keys, comp_id

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f'{comp_id.lower()}.pdb')
        with open(path, 'w') as handle:
            handle.write(_block_from_template(comp_id))
        missed = [k for k in keys if not cg.find_cg_matches(k, path)[0]]
    assert not missed, (f'{comp_id}: {len(missed)} of {len(keys)} keys written '
                        f'from the template graph do not match the same '
                        f'ligand as a block, e.g. {missed[:3]}')


def test_at_least_one_case_actually_exercises_aromatic_bonds():
    """Guard on the test above: it proves nothing if no key is aromatic."""
    from ligand_vdgs.functions import Frags
    mol = lp.perceive_ligand_graph('ATP').mol
    stripped = Frags.manually_remove_Hs(mol, 'single')
    keys = Frags.get_fragments(2, stripped[0][0], 4, 5)
    assert any(any(c.islower() for c in k if c.isalpha()) for k in keys)


def test_a_bond_openbabel_invents_is_deleted():
    """The delete branch: perception can bond atoms the template does not.

    Acetate with the methyl carbon moved to 1.5 A from the carbonyl O, which
    OpenBabel bonds. Left in place it would make that O `D2` and the whole
    carboxylate key wrong.
    """
    moved = ACETATE.format(r='ACT').replace(
        'HETATM    4  CH3 ACT A 900      -0.750  -1.290   0.000',
        'HETATM    4  CH3 ACT A 900       1.900   0.850   0.000')
    perceived = lp.perceive_ligand_instance(moved, NO_TEMPLATE)
    templated = lp.perceive_ligand_instance(moved, 'ACT')
    assert _matches(perceived.obmol, '[O;D2]') == 1, 'case stopped discriminating'
    assert _matches(templated.obmol, '[O;D2]') == 0
    assert _matches(templated.obmol, '[O;D1]') == 2


def test_alt_atom_names_are_matched():
    # CCD alt_atom_id carries the star form of nucleotide primes and the
    # quoted `"N A"` form for heme/chlorophyll nitrogens. Matching atom_id
    # alone sent every CLA and HEC ligand to perception (4.5% of a census).
    by_name = ccd_templates.index_by_name(ccd_templates.get_template('HEC'))
    assert 'NA' in by_name and by_name['NA'].element == 'N'
    assert 'N A' in by_name
