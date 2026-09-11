"""Fragment keys carry heavy-atom degree, so every matcher must run H-free.

RDKit and OpenBabel both count explicit H atoms in the SMARTS `D` primitive.
On a protonated ligand a hydroxyl O then reads `D2` and a methyl C matches
neither `D1` nor `D2`, so a `[O;D1]` key matches nothing -- with no error, in
either direction: mining records no sites and hit finding returns no hits.
"""
import os
import tempfile

import pytest
from rdkit import Chem, RDLogger

from ligand_vdgs.functions import Frags
from ligand_vdgs.generate_vdgs.estimate_frag_cost import add_vdg_miner_paths

RDLogger.DisableLog('rdApp.*')
add_vdg_miner_paths()
import cg  # noqa: E402

# Ethanol, fully protonated, as a prepwizard-style block: the smallest case
# where D is wrong iff the hydrogens are still in the graph.
LIGAND_BLOCK = """\
HETATM    5  C1  LIG B   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    6  C2  LIG B   1       1.520   0.000   0.000  1.00  0.00           C
HETATM    7  O1  LIG B   1       2.100   1.230   0.000  1.00  0.00           O
HETATM    8  H1  LIG B   1      -0.360   1.020   0.000  1.00  0.00           H
HETATM    9  H2  LIG B   1      -0.360  -0.510   0.890  1.00  0.00           H
HETATM   10  H3  LIG B   1      -0.360  -0.510  -0.890  1.00  0.00           H
HETATM   11  H4  LIG B   1       1.880  -0.530   0.890  1.00  0.00           H
HETATM   12  H5  LIG B   1       1.880  -0.530  -0.890  1.00  0.00           H
HETATM   13  H6  LIG B   1       3.060   1.210   0.000  1.00  0.00           H
END
"""
PROTEIN_BLOCK = """\
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  C   ALA A   1      11.500  11.000  10.000  1.00  0.00           C
ATOM      4  O   ALA A   1      12.500  11.000  10.000  1.00  0.00           O
"""
HYDROXYL_KEY = '[O;D1;!R][C;!R;!H0]'


@pytest.fixture
def protonated_pdb():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'test.pdb')
        with open(path, 'w') as handle:
            handle.write(PROTEIN_BLOCK + LIGAND_BLOCK)
        yield path


def test_find_cg_matches_sees_a_degree_key_on_a_protonated_ligand(protonated_pdb):
    matches, _ = cg.find_cg_matches(HYDROXYL_KEY, protonated_pdb)
    assert matches, 'degree-annotated key matched nothing on a protonated ligand'
    names = [tuple(m) for group in matches.values() for m in group]
    assert names == [('O1', 'C2')], names


def test_count_matching_structures_sees_a_degree_key_on_a_protonated_ligand(
        protonated_pdb):
    # Job sizing uses this; a miss here makes every degree-annotated fragment
    # look like it has zero support and drop below the build threshold.
    assert cg.count_matching_structures([HYDROXYL_KEY], protonated_pdb) == {0}


def test_hydrogens_are_what_break_it(protonated_pdb):
    """The case discriminates only if the same key misses with H's in the graph."""
    from openbabel import openbabel as ob
    conv = ob.OBConversion()
    conv.SetInFormat('pdb')
    mol = ob.OBMol()
    conv.ReadString(mol, LIGAND_BLOCK)
    mol.PerceiveBondOrders()
    pattern = ob.OBSmartsPattern()
    assert pattern.Init(HYDROXYL_KEY)
    assert not pattern.Match(mol), 'H-bearing mol now matches; case stopped discriminating'
    mol.DeleteHydrogens()
    assert pattern.Match(mol)


def test_query_ligand_mol_is_h_free_and_reads_degrees(protonated_pdb):
    mol = Frags.get_query_ligand_mol(protonated_pdb, 'CCO')
    assert mol is not None
    assert not [a for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    counts = {sm: len(mol.GetSubstructMatches(Chem.MolFromSmarts(sm)))
              for sm in ['[O;D1]', '[O;D2]', '[C;D1]', '[C;H0]', '[C;!H0]']}
    # `!H0` must stay satisfiable after the strip: it needs implicit H counts,
    # which a mol stripped by clearing valence flags would not have.
    assert counts == {'[O;D1]': 1, '[O;D2]': 0, '[C;D1]': 1,
                      '[C;H0]': 0, '[C;!H0]': 2}, counts


def test_hit_finder_rejects_an_h_bearing_query_mol():
    from ligand_vdgs.score_poses import hit_finder_core
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    with pytest.raises(ValueError, match='hydrogen'):
        hit_finder_core.match_library_frags_to_query(mol, 'unused', {})


def test_ring_size_constraints_still_align_with_degree_tokens():
    # A stray `r` match inside another token would shift every constraint.
    assert cg.ring_size_constraints('[O;D1;!R][C;r6;!H0][Br;D1;!R]') == [None, 6, None]
    assert cg.ring_size_constraints('[C;H0;r5][O;D2;r5]') == [5, 5]


def test_annotation_pickle_stays_aligned_with_the_matches_pickle(protonated_pdb):
    """The whole annotation lookup rests on "same keys, same list positions".

    Nothing downstream can detect a drift: the annotations would simply be
    attached to the wrong atoms of a plausible-looking CG.
    """
    matches, _, annots = cg.find_cg_matches(HYDROXYL_KEY, protonated_pdb,
                                            return_annotations=True)
    assert matches, 'case stopped discriminating'
    assert set(matches) == set(annots)
    for key, match_list in matches.items():
        assert len(match_list) == len(annots[key]), key
        for names, annot in zip(match_list, annots[key]):
            for field in ('heavy_degree', 'num_h', 'formal_charge', 'nbr_elems'):
                assert len(annot[field]) == len(names), (key, field)
    key = next(iter(matches))
    # ethanol's O1 then C2: the hydroxyl O has one heavy neighbour and one H,
    # the carbon two heavy neighbours and two H. Both are wrong if the
    # hydrogens were still in the graph.
    assert matches[key][0] == ['O1', 'C2']
    assert annots[key][0]['heavy_degree'] == [1, 2]
    assert annots[key][0]['num_h'] == [1, 2]
    # C2's heavy neighbour outside the match is C1.
    assert annots[key][0]['nbr_elems'] == ['', 'C']
    assert annots[key][0]['perception'] == 1
