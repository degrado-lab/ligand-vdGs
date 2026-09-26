"""W5 overlap audit instruments (tools/overlap_audit.py)."""
import csv

import numpy as np
import pytest

from ligand_vdgs.functions.frag_enumeration import group_lig_sites_by_overlap
from ligand_vdgs.tools import overlap_audit as oa
from tests.vacuity import assert_discriminates

def test_transitive_merge_is_detected_on_the_grouping_it_audits():
    # 4-atom sites along a chain: each overlaps the next by 2/4 (= the 50% link), the ends are disjoint.
    chain = [(None, None, s) for s in [(0, 1, 2, 3), (2, 3, 4, 5), (4, 5, 6, 7)]]
    groups = group_lig_sites_by_overlap(chain)
    assert len(groups) == 1  # the open bug: single linkage merges disjoint ends
    merged = [set(s) for _, _, s in groups[0]]
    assert oa.group_stats(merged) == (True, True)
    automorph = [{0, 1, 2, 3}, {0, 1, 2, 3}]  # one site under two CG orders: a legitimate group
    assert_discriminates(lambda g: oa.group_stats(g)[0], [merged], [automorph, [{0, 1, 2, 3}, {1, 2, 3, 4}]],
                         "disjoint-pair flag")

def test_relation_orients_dropped_vs_absorber():
    big, small = frozenset(range(6)), frozenset(range(4))
    assert oa.relation(big, small) == "dropped_superset"
    assert oa.relation(small, big) == "dropped_subset"
    assert oa.relation(big, frozenset(range(3, 9))) == "partial"

def _write_hits(d, recs):
    d.mkdir()
    np.savez(d / "hits.npz", frag=np.asarray([r[0] for r in recs]), subset_size=np.full(len(recs), 2, np.int8),
             match_mode=np.full(len(recs), "joint"), aa_bucket=np.full(len(recs), "ALA_GLY"),
             vdg_rmsd=np.asarray([r[2] for r in recs], np.float32), cluster_num_parents=np.ones(len(recs), np.int32),
             q_lig_atom_idx=np.asarray([list(r[3]) + [-1] * (6 - len(r[3])) for r in recs], np.int32))
    with open(d / "hits.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["fragment", "aa_bucket", "subset_size", "vdg_rmsd", "bsr_combo"])
        w.writerows([r[0], "ALA_GLY", 2, f"{r[2]:.4f}", r[1]] for r in recs)

def test_dedup_replay_finds_absorber_across_bsr_token_order(tmp_path):
    # The small fragment wins on RMSD and absorbs the 6-atom fragment written with the BSR tokens swapped.
    _write_hits(tmp_path / "s1", [("small", "A_1;B_2", 0.2, (0, 1, 2, 3)), ("big", "B_2;A_1", 0.4, range(6)),
                                  ("far", "A_1;C_3", 0.3, (0, 1, 2, 3))])
    out = tmp_path / "dedup.tsv"
    oa.main(["dedup", str(tmp_path / "s1"), "--out", str(out)])
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [(r["frag"], r["abs_frag"], r["relation"]) for r in rows] == [("big", "small", "dropped_superset")]

def test_dedup_refuses_placement_mode_and_misaligned_rows(tmp_path):
    _write_hits(tmp_path / "bb", [("x", "A_1", 0.2, (0, 1, 2))])
    with np.load(tmp_path / "bb" / "hits.npz") as d:
        arrs = dict(d)
    np.savez(tmp_path / "bb" / "hits.npz", **(arrs | dict(match_mode=np.asarray(["bb"]))))
    with pytest.raises(ValueError, match="truth-chosen"):
        oa.load_joint_hits(str(tmp_path / "bb"))
    _write_hits(tmp_path / "mis", [("x", "A_1", 0.2, (0, 1, 2)), ("y", "A_1", 0.3, (0, 1, 2))])
    lines = open(tmp_path / "mis" / "hits.tsv").read().splitlines()
    open(tmp_path / "mis" / "hits.tsv", "w").write("\n".join([lines[0], lines[2], lines[1]]) + "\n")
    with pytest.raises(ValueError, match="diverge at row 0"):
        oa.load_joint_hits(str(tmp_path / "mis"))

def test_shared_fraction_nested_requires_containment():
    key = ("1abc", "", "A", 1, (("", "A", 10),))
    a = {key: [frozenset({"C1", "C2", "O1"})], ("1abc", "", "A", 1, (("", "A", 11),)): [frozenset({"C1", "C2", "O1"})]}
    b = {key: [frozenset({"C1", "C2", "O1", "N1"})]}
    assert oa.shared_fraction(a, b, nested=True) == (0.5, 2)  # other vdM residue is not shared evidence
    partial_b = {key: [frozenset({"C2", "N1"})]}
    assert oa.shared_fraction(a, partial_b, nested=True) == (0.0, 2)
    assert oa.shared_fraction(a, partial_b, nested=False) == (0.5, 2)
