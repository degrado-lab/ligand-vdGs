"""The two readers of the fragment-key annotation vocabulary must agree.

`utils.split_bracket_annotations` reads the SMARTS *text* (it has to: the
isotope it produces must be in the string before MolFromSmarts sees it).
`utils._query_primitive` reads a *parsed* query atom's `DescribeQuery()`.
They cannot share code, so they are pinned by this test instead. If session 2
widens what `_query_primitive` accepts, this is where it fails.
"""
import os
import pickle

import pytest
from rdkit import Chem, RDLogger

from ligand_vdgs.functions.utils import (_QUERY_DEGREE_LINE, _QUERY_HCOUNT_LINE,
                                         _BRACKET_ATOM_INNER, _query_primitive,
                                         mol_from_fragment,
                                         split_bracket_annotations)

RDLogger.DisableLog('rdApp.*')

ANNOTATED_KEYS = [
    '[O;D1;!R][P;D4;!R]([O;D1;!R])([O;D1;!R])[O;D2;!R]',
    '[O;D2;!R][P;D4;!R]([O;D1;!R])([O;D1;!R])[O;D2;!R]',
    '[c;H0][O;D1;!R]',
    '[c;H0][N;D2;!R]',
    '[C;r5;!H0][C;r5;!H0][O;r5;D2][C;r5;!H0][C;r5;!H0]',
    '[C;!R;H0](=[O;!R;D1])[O-;D1;!R]',
    '[C;!R;!H0][C;!R;!H0][O;!R;D2][C;!R;!H0]',
    '[Br;D1;!R][c;H0]',
]

# Forms the writer never emits. `_query_primitive` describes them happily;
# split_bracket_annotations must refuse, because an unsatisfiable leftover
# token makes a key stop matching itself and dedup silently records a
# duplicate as new.
OUT_OF_VOCABULARY = ['[O;D1,D2]', '[O;!D1]', '[C;$(C=O);!R]', '[C;H1]', '[C;+0;!R]']


def _text_side(key):
    """(degree, h_token) per bracket atom, in written order."""
    return [split_bracket_annotations(m.group(1), key)[2:]
            for m in _BRACKET_ATOM_INNER.finditer(key)]


def _query_side(key):
    """The same, read back off the parsed query atoms."""
    mol = mol_from_fragment(key)
    assert mol is not None, key
    out = []
    bracket_atoms = [i for i, m in enumerate(_BRACKET_ATOM_INNER.finditer(key))]
    assert len(bracket_atoms) == mol.GetNumAtoms(), (
        f'{key}: every atom must be a bracket atom for this comparison')
    for atom in mol.GetAtoms():
        degree = _query_primitive(atom, _QUERY_DEGREE_LINE)
        hcount = _query_primitive(atom, _QUERY_HCOUNT_LINE)
        out.append((degree, hcount))
    return out


@pytest.mark.parametrize('key', ANNOTATED_KEYS)
def test_both_readers_see_the_same_degree_and_h_on_every_atom(key):
    text = _text_side(key)
    query = _query_side(key)
    assert len(text) == len(query), key
    # Both readers agreeing on "no primitive here" for every atom would pass
    # this test while exercising neither parser. Every key in ANNOTATED_KEYS
    # carries D or H on every atom, so anything less means the text reader
    # silently stopped recognising the vocabulary.
    assert all(d is not None or h is not None for d, h in text), key
    assert any(q is not None for pair in query for q in pair), key
    for i, ((degree, h_token), (q_degree, q_hcount)) in enumerate(zip(text, query)):
        expected_degree = None if degree is None else ((degree, False),)
        assert q_degree == expected_degree, (key, i, q_degree, degree)
        expected_h = {None: None, 'H0': ((0, False),), '!H0': ((0, True),)}[h_token]
        assert q_hcount == expected_h, (key, i, q_hcount, h_token)


@pytest.mark.parametrize('key', OUT_OF_VOCABULARY)
def test_out_of_vocabulary_forms_are_refused_by_the_text_reader(key):
    inner = _BRACKET_ATOM_INNER.search(key).group(1)
    with pytest.raises(ValueError):
        split_bracket_annotations(inner, key)


def test_every_production_key_parses_on_the_text_side():
    """The shipped dictionary must stay inside the closed grammar."""
    path = 'resources/database_frags_dict.pkl'
    if not os.path.exists(path):
        pytest.skip(f'{path} not present')
    with open(path, 'rb') as handle:
        frags_dict = pickle.load(handle)
    keys = [k for bucket in frags_dict.values() for k in bucket]
    assert keys, 'dictionary is empty; test would pass vacuously'
    bad = []
    for key in keys:
        try:
            _text_side(key)
        except ValueError as e:
            bad.append((key, str(e)))
    assert not bad, bad[:10]
