"""Keys must separate the chemistry the old library pooled.

rebuild-notes 5, verification item 1: annotation must differ on phenol/diaryl
ether, aniline/diarylamine, phosphate mono-/di-ester, tert-alcohol/di-tert-alkyl
ether. Each pair has the same heavy-atom skeleton and the same ring context, so
only `D<n>` and `H0`/`!H0` can tell them apart.
"""
import pytest
from rdkit import Chem, RDLogger

from ligand_vdgs.functions import Frags

RDLogger.DisableLog('rdApp.*')


def keys_for(smiles, bond_radius=2, min_size=4, max_size=5):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    Chem.SanitizeMol(mol)
    stripped = Frags.manually_remove_Hs(mol, 'single')
    assert stripped is not None, smiles
    return set(Frags.get_fragments(bond_radius, stripped[0][0],
                                   min_size, max_size))


@pytest.mark.parametrize('name,a,b', [
    ('phenol/diaryl ether', 'Oc1ccccc1', 'c1ccccc1Oc1ccccc1'),
    ('aniline/diarylamine', 'Nc1ccccc1', 'c1ccccc1Nc1ccccc1'),
    ('phosphate mono/di-ester', 'COP(=O)(O)O', 'COP(=O)(O)OC'),
    ('tert-alcohol/di-tert-ether', 'CC(C)(C)O', 'CC(C)(C)OC(C)(C)C'),
])
def test_annotation_separates_pooled_chemistry(name, a, b):
    keys_a, keys_b = keys_for(a), keys_for(b)
    assert keys_a and keys_b, name
    assert not (keys_a & keys_b), (name, sorted(keys_a & keys_b))


def test_bridging_and_terminal_heteroatoms_get_different_degrees():
    # The 47.1% defect directly: the ether O is D2, the hydroxyl O is D1, even
    # though both are degree 1 inside the cut fragment.
    assert any(';D2]' in k for k in keys_for('c1ccccc1Oc1ccccc1'))
    assert any(';D1]' in k for k in keys_for('Oc1ccccc1'))


def test_ring_size_still_separates_furanose_from_pyranose():
    # The annotation must not swamp the ring split that was already there.
    assert not (keys_for('C1CCOC1') & keys_for('C1CCOCC1'))


def test_carbon_h_flag_is_binary_and_present_on_aromatic_carbon():
    keys = keys_for('Cc1ccc(O)cc1')  # all-carbon fragments are filtered out
    assert keys
    for key in keys:
        assert 'H1' not in key and 'H2' not in key and 'H3' not in key, key
    assert any('[c;H0]' in k for k in keys), keys       # the substituted ring C
    assert any('[c;!H0]' in k for k in keys), keys


def test_bracket_hydrogens_in_the_source_smiles_do_not_zero_the_h_flag():
    """RemoveAllHs(sanitize=False) drops graph H without crediting the heavy atom.

    Discriminating case: written with bracket H, the methyl carbon of ethanol
    read zero hydrogens and would have been keyed `H0` -- the same fragment as a
    quaternary carbon.
    """
    assert keys_for('[H]OCC[H]', bond_radius=2, min_size=3) == \
        keys_for('OCC', bond_radius=2, min_size=3)


def test_every_atom_of_a_key_is_annotated():
    # A partially annotated key is a silent mis-key, not a loose query.
    from ligand_vdgs.functions.utils import (_BRACKET_ATOM_INNER,
                                             split_bracket_annotations)
    checked = 0
    for smiles in ['Cc1ccc(O)cc1', 'COP(=O)(O)OC', 'C1CCOC1', 'CC(=O)[O-]']:
        keys = keys_for(smiles)
        # Without this the test asserts nothing when a ligand yields no
        # fragments -- which happens for real reasons (the all-carbon filter
        # drops toluene entirely).
        assert keys, smiles
        for key in keys:
            checked += 1
            n_bracket = len(_BRACKET_ATOM_INNER.findall(key))
            assert n_bracket == key.count('[') , key
            for inner in _BRACKET_ATOM_INNER.findall(key):
                _, ring, degree, h_token = split_bracket_annotations(inner, key)
                assert degree is not None or h_token is not None, (key, inner)
    assert checked >= 4, checked
