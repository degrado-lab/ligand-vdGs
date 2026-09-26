"""self_recovery.build_units: strict truth join, rank pool, falsifier; hit_combos."""
from ligand_vdgs.tools import self_recovery as sr
from tests.vacuity import assert_discriminates

F, G = "[C;!R;!H0][C;!R;H0](=[O;!R;D1])[O;!R;D1]", "[C;!R;!H0][C;!R;H0]=[O;!R;D1]"
CG = frozenset({"CG", "CD", "OE1", "OE2"})
R1, R2 = frozenset({("", "A", 10)}), frozenset({("", "A", 11)})

def _truth(cid, residues=R1, is_nr=False, chain="A", frag=F, names=CG):
    return dict(frag=frag, subset_size=1, charge_sign="neg", aa_bucket="ARG", cluster_id=cid,
                num_parents=5, is_nr=is_nr, cg_chain=chain, cg_resnum=502, cg_names=names,
                residues=residues)

def _hit(cid, rmsd, residues=R1, frag=F, names=CG, label="A502", support=5):
    return dict(lig_instance=label, frag=frag, subset_size=1, charge_sign="neg", aa_bucket="ARG",
                vdg_cluster_id=cid, vdg_cluster_num_parents=support, vdg_rmsd=f"{rmsd:.4f}",
                rmsd_threshold="0.6000", cg_names=names, residues=residues)

def _lig(truth, hits, classes=None):
    return dict(instances={"A502": dict(chain="A", resnum=502), "B502": dict(chain="B", resnum=502)},
                enumerated_combos={"A502": {R1, R2}, "B502": {R1}}, truth=truth, hits=hits,
                contact_classes=classes or {})

def _units(lig):
    return sr.build_units(lig, sr._hits_frame(lig["hits"]))

def _unit(rows, residues="A10", label="A502"):
    (row,) = [r for r in rows if r["bsr"] == residues and r["lig_instance"] == label]
    return row

def test_strict_requires_own_residues_loose_does_not():
    # cluster 7 is s's (via R1), but the only hit reaching it came through R2.
    lig = _lig([_truth(7)], [_hit(7, 0.1, residues=R2)])
    row = _unit(_units(lig))
    assert (row["matched_strict"], row["matched_loose"]) == (False, True)
    lig_ok = _lig([_truth(7)], [_hit(7, 0.1)])
    assert_discriminates(lambda l: _unit(_units(l))["matched_strict"],
                         [lig_ok], [lig], "strict residue match")

def test_instance_is_part_of_the_unit():
    # truth observed on chain B's copy; a hit on chain A's copy must not count.
    lig = _lig([_truth(7, chain="B")], [_hit(7, 0.1, label="A502")])
    row = _unit(_units(lig), label="B502")
    assert row["matched_strict"] is False
    lig_b = _lig([_truth(7, chain="B")], [_hit(7, 0.1, label="B502")])
    assert _unit(_units(lig_b), label="B502")["matched_strict"] is True

def test_rank_pool_counts_only_overlapping_cg_atoms_strictly_below():
    other_site = frozenset({"C1", "C2", "O1"})
    hits = [_hit(7, 0.30), _hit(8, 0.10, frag=G, names=frozenset({"CD", "OE1", "CG"})),
            _hit(9, 0.05, frag=G, names=other_site), _hit(10, 0.30, residues=R2),
            _hit(11, 0.20)]
    row = _unit(_units(_lig([_truth(7)], hits)))
    # below 0.30 and overlapping CG: cluster 8 (other frag) and 11 (same frag); 9 is disjoint,
    # 10 ties (not below).
    assert (row["rank_rmsd_cross"], row["rank_rmsd_within"]) == (3, 2)
    assert row["n_pool_cross"] == 4 and row["n_pool_within"] == 3

def test_hit_combos_counts_per_instance_frag_residues():
    hits = [_hit(7, 0.1), _hit(8, 0.2), _hit(9, 0.1, residues=R2), _hit(7, 0.1, label="B502")]
    got = sr.hit_combos(sr._hits_frame(hits))
    assert got == {("A502", F, R1): 2, ("A502", F, R2): 1, ("B502", F, R1): 1}

def test_falsifier_flags_rep_self_not_at_zero():
    ok = _lig([_truth(7, is_nr=True)], [_hit(7, 0.0)])
    off = _lig([_truth(7, is_nr=True)], [_hit(7, 0.4)])
    missing = _lig([_truth(7, is_nr=True)], [])
    not_rep = _lig([_truth(7, is_nr=False)], [_hit(7, 0.4)])
    assert_discriminates(lambda l: _unit(_units(l))["falsifier_violation"],
                         [off, missing], [ok, not_rep], "falsifier")

def test_contact_class_lookup_by_instance_names_residues():
    classes = {("A502", CG, R1): "polar"}
    row = _unit(_units(_lig([_truth(7)], [_hit(7, 0.1)], classes)))
    assert row["contact_class"] == "polar"
    assert _unit(_units(_lig([_truth(7)], [_hit(7, 0.1)])))["contact_class"] == "NA"

def test_parse_bsr_combo_matches_truth_residue_keys():
    # hit strings are seg:chain:resnum; truth keys come from padded npz strings.
    assert sr.parse_bsr_combo(":A:10;:B:7") == frozenset({("", "A", 10), ("", "B", 7)})
    assert sr._res_key(" ", "A ", "10") == ("", "A", 10)
