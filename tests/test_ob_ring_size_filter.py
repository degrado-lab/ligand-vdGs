"""OpenBabel's r<n> is looser than RDKit's; the miner must re-impose RDKit's.

See docs/pitfalls.md, "RDKit and OpenBabel disagree on what `r<n>` means".
"""
import os
import pickle
import sys

import pytest
from rdkit import Chem, RDLogger

from ligand_vdgs.generate_vdgs.estimate_frag_cost import add_vdg_miner_paths

RDLogger.DisableLog('rdApp.*')
add_vdg_miner_paths()
import cg  # noqa: E402
from openbabel import openbabel as ob  # noqa: E402

SPIRO = 'C1CCC2(CC1)CCCO2'  # 1-oxaspiro[4.5]decane; spiro C is r5 to RDKit, r5+r6 to OB


def _ob_mol(smiles):
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    mol = ob.OBMol()
    conv.ReadString(mol, smiles)
    return mol


@pytest.mark.parametrize('key', ['[C;r6]', '[C;r6][C;r6][C;r6]'])
def test_filter_restores_rdkit_semantics_on_a_spiro_system(key):
    """The discriminating case: unfiltered OB matches the spiro carbon, RDKit does not."""
    rdkit_matches = Chem.MolFromSmiles(SPIRO).GetSubstructMatches(
        Chem.MolFromSmarts(key), uniquify=True)
    mol = _ob_mol(SPIRO)
    pattern = ob.OBSmartsPattern()
    assert pattern.Init(key)
    pattern.Match(mol)
    raw = list(pattern.GetUMapList())
    constraints = cg.ring_size_constraints(key)
    kept = [m for m in raw if cg.match_satisfies_ring_sizes(mol, m, constraints)]

    assert len(raw) > len(rdkit_matches), 'case no longer discriminates'
    assert len(kept) == len(rdkit_matches)


def test_negated_ring_primitive_is_not_read_as_a_constraint():
    """`!r5` must stay OpenBabel's job; only a positive r<n> constrains."""
    assert cg.ring_size_constraints('[C;r6;!r3;!r4;!r5][C;!R]') == [6, None]
    assert cg.ring_size_constraints('[C;!R][O;!R]') == [None, None]


def test_constraints_align_with_openbabel_atom_order_for_every_shipped_key():
    """A misaligned constraint list would silently filter the wrong positions."""
    path = 'resources/database_frags_dict.pkl'
    if not os.path.exists(path):
        pytest.skip(f'{path} not present')
    with open(path, 'rb') as handle:
        frags_dict = pickle.load(handle)
    keys = [k for bucket in frags_dict.values() for k in bucket]
    assert keys, 'fragment dict is empty'
    for key in keys:
        pattern = cg.compile_smarts_patterns([key])[0]
        constraints = cg.ring_size_constraints(key)
        assert len(constraints) == pattern.NumAtoms(), key


def test_atom_count_mismatch_raises_rather_than_disabling_the_filter():
    with pytest.raises(ValueError, match='refusing to mine'):
        cg.check_ring_constraints('[C;r6]', ob.OBSmartsPattern(), [1, 2, 3])
