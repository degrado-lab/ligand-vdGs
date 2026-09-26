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
from vacuity import assert_discriminates

RDLogger.DisableLog('rdApp.*')

pytestmark = pytest.mark.skipif(
    not os.path.isfile(ccd_templates.template_db_path()),
    reason=f'{ccd_templates.template_db_path()} not built '
           '(scripts/build_ccd_templates.py)')

def test_store_identity_pins_template_path_parser_and_stat():
    identity = ccd_templates.store_identity()
    path = os.path.abspath(ccd_templates.template_db_path())
    stat = os.stat(path)
    assert identity['ccd_dir'] == os.path.abspath(ccd_templates.ccd_dir())
    assert identity['template_db'] == path
    assert identity['parser_version'] == ccd_templates.TEMPLATE_PARSER_VERSION
    assert identity['template_db_size'] == stat.st_size
    assert identity['template_db_mtime_ns'] == stat.st_mtime_ns

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
    from ligand_vdgs.functions import frag_enumeration
    from ligand_vdgs.functions.utils import fragment_keys_equivalent
    for comp_id in ('ATP', 'TYR', 'HEM'):
        mol = lp.perceive_ligand_graph(comp_id).mol
        stripped = frag_enumeration.manually_remove_Hs(mol, 'single')
        assert stripped is not None, comp_id
        keys = frag_enumeration.get_fragments(2, stripped[0][0], 4, 5)
        assert keys, comp_id
        for key in keys:
            assert fragment_keys_equivalent(key, key), (comp_id, key)


def test_aromatic_ring_atoms_pool_across_pyridine_and_pyrrole():
    # The vocabulary decision `ring_query_for_atom` documents: aromatic atoms
    # carry no ring-size primitive, so a 5- and a 6-ring aromatic carbon are
    # the same key atom. Template-built graphs must not break that.
    from ligand_vdgs.functions import frag_enumeration
    for comp_id in ('ATP', 'TYR'):
        mol = lp.perceive_ligand_graph(comp_id).mol
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                assert frag_enumeration.ring_query_for_atom(atom) is None


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
    from ligand_vdgs.functions import frag_enumeration
    from ligand_vdgs.generate_vdgs.estimate_frag_cost import add_vdg_miner_paths
    add_vdg_miner_paths()
    import cg

    mol = lp.perceive_ligand_graph(comp_id).mol
    stripped = frag_enumeration.manually_remove_Hs(mol, 'single')
    assert stripped is not None, comp_id
    keys = sorted(frag_enumeration.get_fragments(2, stripped[0][0], 4, 5))
    assert keys, comp_id

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f'{comp_id.lower()}.pdb')
        with open(path, 'w') as handle:
            handle.write(_block_from_template(comp_id))
        missed = [k for k in keys if not cg.find_cg_matches(k, path)[0]]
    assert not missed, (f'{comp_id}: {len(missed)} of {len(keys)} keys written '
                        f'from the template graph do not match the same '
                        f'ligand as a block, e.g. {missed[:3]}')


def _any_key_aromatic(keys):
    """Aromaticity from the parsed query mol. The shipped `c.islower()` form fired
    on the `l` of Cl and on the lowercase ring primitive `r` -- and every real key
    here carries `r5`/`r6` -- so it could not reject an aliphatic set (audit V7).
    """
    return any(a.GetIsAromatic() for k in keys for a in Chem.MolFromSmarts(k).GetAtoms())

def test_at_least_one_case_actually_exercises_aromatic_bonds():
    """Guard on the test above: it proves nothing if no key is aromatic."""
    from ligand_vdgs.functions import frag_enumeration
    mol = lp.perceive_ligand_graph('ATP').mol
    stripped = frag_enumeration.manually_remove_Hs(mol, 'single')
    keys = tuple(frag_enumeration.get_fragments(2, stripped[0][0], 4, 5))
    assert not [k for k in keys if Chem.MolFromSmarts(k) is None], keys
    assert_discriminates(_any_key_aromatic, [keys],
                         [('[C;D4]-[Cl;D1]', '[C;r6;D3]')], 'ATP keys aromatic')


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


# Real parent-DB block (1dp9 IMD 500): the CONECT naming a missing atom is ignored, so
# HN1 (0.4 A from C5) stays unbonded and survives DeleteHydrogens under a template name.
IMD_STRAY_H = """\
HETATM 1269  N1  IMD A 500      47.634  53.196  34.362  1.00 26.80           N
HETATM 1270  C2  IMD A 500      47.257  52.778  35.618  1.00 26.69           C
HETATM 1271  N3  IMD A 500      47.433  51.426  35.682  1.00 25.85           N
HETATM 1272  C4  IMD A 500      47.916  50.990  34.483  1.00 28.80           C
HETATM 1273  C5  IMD A 500      48.043  52.099  33.653  1.00 27.12           C
HETATM 1274  HN1 IMD A 500      47.982  52.527  33.690  1.00 26.80           H
HETATM 1275  HC2 IMD A 500      46.899  53.496  36.341  1.00 26.69           H
HETATM 1276  HC4 IMD A 500      48.123  49.941  34.333  1.00 28.80           H
HETATM 1277  HC5 IMD A 500      48.383  52.211  32.634  1.00 27.12           H
CONECT 1269 1270 1273 1320 1274
CONECT 1270 1269 1271 1275
CONECT 1271 1270 1272
CONECT 1272 1271 1273 1276
CONECT 1273 1269 1272 1277
CONECT 1274 1269
CONECT 1275 1270
CONECT 1276 1272
CONECT 1277 1273
END
"""

def _raw_after_delete_hydrogens(block):
    """The pre-fix path: read, perceive, DeleteHydrogens, nothing else."""
    mol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat('pdb')
    conv.ReadString(mol, block)
    mol.PerceiveBondOrders()
    mol.DeleteHydrogens()
    return mol

def _has_no_hydrogen(mol):
    return all(mol.GetAtom(i).GetAtomicNum() != 1 for i in range(1, mol.NumAtoms() + 1))

def test_a_hydrogen_deletehydrogens_keeps_under_a_template_name_is_removed():
    # Falsifier: the kept HN1 maps to template HN1, the template bond N1-HN1 is laid
    # down, and N1 reads D3 (a pyrrole-type n with a substituent) instead of D2.
    templated = lp.perceive_ligand_instance(IMD_STRAY_H, 'IMD')
    assert templated.provenance == lp.PERCEPTION_CCD_TEMPLATE
    assert_discriminates(_has_no_hydrogen, [templated.obmol],
                         [_raw_after_delete_hydrogens(IMD_STRAY_H)], 'IMD stray HN1')
    assert _by_name(templated.obmol)['N1'] == (2, 1, 0)
    assert _matches(templated.obmol, '[n;D2]') == 2 and _matches(templated.obmol, '[n;D3]') == 0

def test_an_unbonded_hydrogen_with_an_unknown_name_does_not_block_the_template():
    # A86/CLA/FAD pattern: a misplaced H whose name the template lacks sent the whole
    # instance to OpenBabel.
    stray = ACETATE.format(r='ACT').replace(
        'END', 'HETATM    5  HQQ ACT A 900       6.000   6.000   6.000  1.00  0.00           H\nEND')
    assert not _has_no_hydrogen(_raw_after_delete_hydrogens(stray)), 'case stopped discriminating'
    perceived = lp.perceive_ligand_instance(stray, 'ACT')
    assert perceived.provenance == lp.PERCEPTION_CCD_TEMPLATE
    assert set(lp.pdb_atom_names(perceived.obmol).values()) == {'C', 'O', 'OXT', 'CH3'}

def test_alt_atom_names_are_matched():
    # CCD alt_atom_id carries the star form of nucleotide primes and the
    # quoted `"N A"` form for heme/chlorophyll nitrogens. Matching atom_id
    # alone sent every CLA and HEC ligand to perception (4.5% of a census).
    template = ccd_templates.get_template('HEC')
    assert ccd_templates.resolve_names(['NA'], template)[0].element == 'N'
    assert ccd_templates.resolve_names(['N A'], template)[0].name == 'NA'

ABU_LEGACY = """\
HETATM 1410  N   ABU A1457      77.773  33.156  47.937  1.00 30.07           N1+
HETATM 1411  CA  ABU A1457      77.466  34.521  48.359  1.00 28.85           C
HETATM 1412  CB  ABU A1457      78.726  35.420  48.123  1.00 29.69           C
HETATM 1413  CG  ABU A1457      79.157  35.943  49.499  1.00 30.77           C
HETATM 1414  CD  ABU A1457      80.458  36.693  49.494  1.00 32.13           C
HETATM 1415  OE1 ABU A1457      80.730  37.350  48.483  1.00 33.40           O
HETATM 1416  OE2 ABU A1457      81.151  36.616  50.536  1.00 33.89           O1-
HETATM 1417 H__1 ABU A1457      76.966  32.568  48.085  1.00 30.07           H
HETATM 1420  HA1 ABU A1457      77.211  34.525  49.419  1.00 28.85           H
HETATM 1422  HB1 ABU A1457      78.457  36.261  47.484  1.00 29.69           H
END
"""
ABU_LEGACY_NAMES = ['N', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2']
ABU_CURRENT_NAMES = ['N', 'CD', 'CB', 'CG', 'C', 'O', 'OXT']

def _one_to_one(table, names):
    return all(n in table for n in names) and len({table[n].name for n in names}) == len(names)

def test_a_fully_legacy_deposition_resolves_under_alt_names():
    # 4atq ABU 1457 (67 roster instances): legacy CA/CD are current CD/C, so the merged
    # table sends legacy CD to current CD and collides with legacy CA.
    template = ccd_templates.get_template('ABU')
    alt = {a.alt_name: a for a in template.atoms if a.alt_name}
    assert_discriminates(lambda t: _one_to_one(t, ABU_LEGACY_NAMES), [alt],
                         [ccd_templates._index_by_name(template)], 'ABU legacy names')
    got = dict(zip(ABU_LEGACY_NAMES, (e.name for e in ccd_templates.resolve_names(ABU_LEGACY_NAMES, template))))
    assert (got['CA'], got['CD'], got['OE2']) == ('CD', 'C', 'OXT')
    perceived = lp.perceive_ligand_instance(ABU_LEGACY, 'ABU')
    assert perceived.provenance == lp.PERCEPTION_CCD_TEMPLATE
    # Chemistry follows the deposited atom: legacy CD is the carboxyl C.
    by_name = _by_name(perceived.obmol)
    assert by_name['CD'] == (3, 0, 0) and by_name['CA'][0] == 2

def test_current_names_win_when_both_vocabularies_map():
    # An alt-first resolver would read current CD (the alpha C) as legacy CD (carboxyl).
    template = ccd_templates.get_template('ABU')
    alt = {a.alt_name: a for a in template.atoms if a.alt_name}
    partial = ['N', 'CD', 'CB', 'CG']   # carboxyl unobserved, so alt names also map
    assert _one_to_one(alt, partial) and alt['CD'].name == 'C', 'case stopped discriminating'
    assert [e.name for e in ccd_templates.resolve_names(partial, template)] == partial
    assert [e.name for e in ccd_templates.resolve_names(ABU_CURRENT_NAMES, template)] == ABU_CURRENT_NAMES

def test_a_genuine_duplicate_name_resolves_under_neither_vocabulary():
    assert ccd_templates.resolve_names(['C', 'C', 'O'], ccd_templates.get_template('ACT')) is None

FES_MISSING_S2 = """\
HETATM 9600  FE1 FES E 501       8.958  81.453  73.804  1.00 40.05          Fe
HETATM 9601  FE2 FES E 501      10.333  81.857  76.137  1.00 40.33          Fe
HETATM 9602  S1  FES E 501       8.252  81.008  75.864  1.00 39.83           S
END
"""

def test_templating_keeps_deposited_names_and_leaves_phantoms_unnamed():
    # Falsifier: the template's bond edits clear OB's chains flag; re-perception on the
    # next GetResidue renamed FE1 -> FE and named the phantom, so the miner stored
    # names absent from the file and counted matches reaching unobserved atoms.
    perceived = lp.perceive_ligand_instance(FES_MISSING_S2, 'FES')
    assert perceived.provenance == lp.PERCEPTION_CCD_TEMPLATE
    reperceived = ob.OBMol(perceived.obmol)
    reperceived.UnsetFlag(ob.OB_CHAINS_MOL)
    assert_discriminates(lambda m: sorted(lp.pdb_atom_names(m).values()) == ['FE1', 'FE2', 'S1'],
                         [perceived.obmol], [reperceived], 'FES names after templating')
    assert len(lp.phantom_atom_indices(perceived.obmol)) == 1
