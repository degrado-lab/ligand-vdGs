"""W4 score_eval: contract refusal, M4 invariants, tie-fair endpoint, nesting, freeze and split gates."""
import itertools
import json
import os
from types import SimpleNamespace

import numpy as np
import pytest

from ligand_vdgs.tools import score_eval as se
from tests.vacuity import assert_discriminates

def contract_hits(n=4, m=3):
    rng = np.random.default_rng(0)
    return dict(frag=np.full(n, "F"), subset_size=np.full(n, 2, np.int8), match_mode=np.full(n, "bb"),
                charge_sign=np.full(n, "neut"), aa_bucket=np.full(n, "ALA_bb"), nr_idx=np.arange(n, dtype=np.int32),
                bb_rmsd=rng.random(n).astype(np.float32), vdg_rmsd=np.full(n, np.nan, np.float32),
                held_out_cg_rmsd=rng.random(n).astype(np.float32) * 4, cg_bb_dist=np.ones(n, np.float32),
                cluster_num_parents=np.arange(1, n + 1, dtype=np.int32), nr_parent_biounit=np.full(n, "1abc_1"),
                lig_instance=np.full(n, ""), min_bb_dist=np.ones(n, np.float32),
                placed_cg=rng.random((n, m, 3)).astype(np.float32), q_lig_atom_idx=np.tile(np.arange(m, dtype=np.int32), (n, 1)),
                placed_cg_element=np.full((n, m), "C"))

def _loads(tmp_path, arrays):
    p = str(tmp_path / f"h{len(os.listdir(tmp_path))}.npz")
    np.savez(p, **arrays)
    try: se.load_hits(p)
    except ValueError: return False
    return True

def test_loader_refuses_non_contract(tmp_path):
    good = contract_hits()
    bad = dict(cohort_bb_v1=dict(structure=np.full(4, "s"), fragment=good["frag"], subset_size=good["subset_size"],
                                 bb_rmsd=good["bb_rmsd"], held_out_cg_rmsd=good["held_out_cg_rmsd"],
                                 cg_bb_dist=good["cg_bb_dist"], vdg_rmsd=good["vdg_rmsd"],
                                 vdg_cluster_num_parents=good["cluster_num_parents"].astype(np.float32)),
               f4_support=dict(good, cluster_num_parents=good["cluster_num_parents"].astype(np.float32)),
               ragged=dict(good, q_lig_atom_idx=good["q_lig_atom_idx"][:, :2]))
    assert_discriminates(lambda a: _loads(tmp_path, a), [good], bad.values(), "contract loader")

def _cloud(m=40, N=4, seed=1):
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(N, 3)) * 1.5 + rng.normal(scale=0.4, size=(m, N, 3))).astype(np.float32), rng.integers(1, 9, m)

AUTOS = ((0, 1, 2, 3), (1, 0, 2, 3))

def test_m4_sigma_relabel_invariance():
    P, n = _cloud()
    Q = P.copy()
    Q[::3] = P[::3][:, list(AUTOS[1])]  # relabel every third hit by the swap automorphism
    assert np.array_equal(se.consensus_density(Q, n, AUTOS, 1.0), se.consensus_density(P, n, AUTOS, 1.0))
    # Falsifier: ignoring the automorphism group makes the relabelling visible.
    assert not np.array_equal(se.consensus_density(Q, n, AUTOS[:1], 1.0), se.consensus_density(P, n, AUTOS[:1], 1.0))

def test_m4_far_hit_self_and_single_count():
    P, n = _cloud()
    out = se.consensus_density(np.concatenate([P, P[:1] + 100.0]), np.append(n, 7), AUTOS, 1.0)
    assert np.array_equal(out[:-1], se.consensus_density(P, n, AUTOS, 1.0)) and out[-1] == 0
    assert se.consensus_density(P[:1], n[:1], AUTOS, 2.0).tolist() == [0]
    # A near-symmetric CG: both images of h' lie within r of h, yet h' votes once.
    sym = np.array([[[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, 0, 1]]], np.float32)
    assert se.consensus_density(np.concatenate([sym, sym + 0.05]), np.array([3, 5]), AUTOS, 2.0).tolist() == [5, 3]

def test_m4_brute_force():
    P, n = _cloud(m=25, seed=3)
    d = np.array([[min(np.sqrt(((P[i] - P[j][list(a)]) ** 2).sum(-1).mean()) for a in AUTOS) for j in range(25)] for i in range(25)])
    assert se.consensus_density(P, n, AUTOS, 1.0).tolist() == [(n * ((d[i] <= 1.0) & (np.arange(25) != i))).sum()
                                                               for i in range(25)]

def _brute_topk(s, p, k):
    """Average over every tie-breaking order consistent with the scores of P(>= 1 success in the top k)."""
    perms = [q for q in itertools.permutations(range(len(s))) if all(s[q[i]] >= s[q[i + 1]] for i in range(len(s) - 1))]
    return np.mean([1 - np.prod([1 - p[j] for j in q[:k]]) for q in perms])

def test_tie_fair_endpoint():
    met = se.unit_metrics(*se.candidates(np.array([1.0, 1.0, 0.0]), np.arange(3), np.array([0.5, 3.0, 0.1])))
    assert met["top1@2.0"] == 0.5 and met["top1_err"] == 1.75
    rng = np.random.default_rng(4)
    for _ in range(30):
        s, ok, frac = rng.integers(0, 3, 6).astype(float), rng.random(6) < 0.3, rng.random(6) * (rng.random(6) < 0.5)
        for k in (1, 2, 5):
            assert se.expected_topk(s, ok.astype(float), k) == pytest.approx(_brute_topk(s, ok, k)), (s, ok, k)
            assert se.expected_topk(s, frac, k) == pytest.approx(_brute_topk(s, frac, k)), (s, frac, k)

def test_candidates_max_not_sum():
    """Stage-1 candidate 7 has two hits tied at its max; its summed score (5) would outrank candidate 9 (3)."""
    s, cid, err = np.array([2.0, 2.0, 1.0, 0.0, 3.0]), np.array([7, 7, 7, 3, 9]), np.array([0.5, 3.0, 0.1, 1.0, 2.5])
    sc, p, e = se.candidates(s, cid, err)
    assert sc.tolist() == [0.0, 2.0, 3.0] and p[:, 2].tolist() == [1.0, 0.5, 0.0] and e.tolist() == [1.0, 1.75, 2.5]
    assert se.unit_metrics(sc, p, e)["top1_err"] == 2.5
    summed = np.bincount(np.unique(cid, return_inverse=True)[1], weights=s)
    assert se.unit_metrics(summed, p, e)["top1_err"] != 2.5

def test_evidence_keys_through_category_codes(tmp_path, monkeypatch):
    """Hits list GLY_bb first (1-row chunks give loader codes GLY 0, ALA 1); the sidecar table is sorted (ALA 0, GLY 1)."""
    from ligand_vdgs.functions import evidence_weights as ew
    rows = dict(frag=np.full(6, "F"), subset_size=np.full(6, 2, np.int8), charge_sign=np.full(6, "neut"),
                aa_bucket=np.repeat(["ALA_bb", "GLY_bb"], 3), nr_idx=np.tile(np.arange(3, dtype=np.int32), 2),
                first_stage_cluster_id=np.array([1, 1, 2, 1, 2, 2], np.int32), cluster_size=np.ones(6, np.int32),
                cluster_num_parents=np.array([3, 4, 5, 6, 7, 8], np.int32), n_eff=np.array([2, 2, 1, 1, 2, 2], np.int32),
                stage1_n_obs=np.ones(6, np.int32), stage1_n_entries=np.ones(6, np.int32))
    tables = {}
    sc = ew.assemble([ew.encode(rows, tables)], tables)
    ew.save_npz_atomic(side := str(tmp_path / "side.npz"), sc)
    h = contract_hits(n=4)
    h.update(aa_bucket=np.array(["GLY_bb", "ALA_bb", "GLY_bb", "ALA_bb"]), nr_idx=np.array([2, 0, 0, 2], np.int32))
    s1 = ew.stage1_ids(sc)[ew.positions(sc, h)]  # W5's own string-keyed lookup is the reference
    ents = {int(s1[0]): ["1AAA", "1BBB", "1CCC", "9ZZZ"]} | {int(x): ["1DDD"] for x in s1[1:]}
    ew.save_npz_atomic(par := str(tmp_path / "par.npz"), ew.assemble_parents(
        np.array([k for k, v in ents.items() for _ in v], np.int64), np.array([e for v in ents.values() for e in v])))
    (fam := tmp_path / "fam.txt").write_text("1AAA_1 1BBB_2\n1BBB_1 1CCC_1\n1DDD_1\n")
    np.savez(hp := str(tmp_path / "h.npz"), **h)
    np.savez(m8 := str(tmp_path / "m8.npz"), s1_id=np.sort(s1), m8=np.arange(10, 14, dtype=np.int32))
    monkeypatch.setattr(se, "_ROWS", 1)
    g, ev = se.load_hits(hp), se.load_evidence(side, par, str(fam), m8)
    assert g["cats"]["aa_bucket"].tolist() == ["GLY_bb", "ALA_bb"]
    se.attach_evidence(g, ev)
    assert g["s1"].tolist() == s1.tolist() and g["m9"].tolist() == [2, 2, 1, 1] and g["s1sup"].tolist() == [15, 7, 6, 5]
    assert g["m8"].tolist() == (np.argsort(np.argsort(s1)) + 10).tolist()  # each hit gets its own S1's M8
    # 1AAA-1BBB-1CCC chain through shared clusters (one component); unmapped 9ZZZ is its own.
    assert se.m8_unionfind(ev, g["s1"]).tolist() == [2, 1, 1, 1] and ev["m8u_entries"] == [7, 1]
    g2 = se.load_hits(hp)
    g2["nr_idx"][0] = 5
    with pytest.raises(KeyError, match="absent from the sidecar"): se.attach_evidence(g2, ev)

def test_maximal_frags():
    lab = dict(F=[(1, 2, 3)], G=[(1, 2, 3, 4)], H=[(1, 2), (7, 8)], E=[(3, 2, 1)])
    # F and E share the atom set {1,2,3} (both nested in G); H has a site outside G, so it stays.
    assert se.maximal_frags({"F", "G", "H", "E"}, lab) == {"G", "H"}
    assert se.maximal_frags({"F", "H", "E"}, lab) == {"F", "H", "E"}

def test_streamed_factorization_across_chunks(tmp_path, monkeypatch):
    """Categories first seen in later chunks, in a different order per chunk, decode exactly."""
    h = contract_hits(n=7)
    h.update(frag=np.array(["G", "F", "F", "H", "G", "F", "Hx"]), lig_instance=np.array(["B", "A", "A", "C", "B", "A", "A"]))
    np.savez_compressed(p := str(tmp_path / "h.npz"), **h)
    monkeypatch.setattr(se, "_ROWS", 3)
    g = se.load_hits(p)
    for k in ("frag", "lig_instance", "aa_bucket", "match_mode", "charge_sign"):
        assert g[k].dtype == np.int32 and g["cats"][k][g[k]].tolist() == h[k].tolist(), k
    assert np.array_equal(g["placed_cg"], h["placed_cg"]) and "placed_cg_element" not in g
    # Falsifier: chunk-local codes (no global category map) mis-decode rows past the first chunk.
    local = np.concatenate([np.unique(h["frag"][s:s + 3], return_inverse=True)[1] for s in range(0, 7, 3)])
    assert g["cats"]["frag"][local].tolist() != h["frag"].tolist()

def test_unit_groups_nan_joint_and_instances(tmp_path):
    h = contract_hits(n=6)
    h["lig_instance"] = np.array(["A", "A", "B", "B", "C", "A"])
    h["match_mode"] = np.array(["bb", "bb", "bb", "bb", "bb", "joint"])
    h["held_out_cg_rmsd"][4] = np.nan
    np.savez(p := str(tmp_path / "h.npz"), **h)
    groups, nan_units = se.unit_groups(se.load_hits(p))
    assert {k: v.tolist() for k, v in groups.items()} == {("A", "F"): [0, 1], ("B", "F"): [2, 3]}
    assert nan_units == [("C", "F")]

def test_bootstrap_weights_and_estimates():
    gidx = np.array([0, 0, 1, 2])
    Wg = se.resample_weights(gidx, 50, 0, True)
    assert (Wg[:, 0] == Wg[:, 1]).all() and (Wg.sum(1) >= 0).all()
    Ws = se.resample_weights(gidx, 50, 0, False)
    assert (Ws.sum(1) == 4).all() and not (Ws[:, 0] == Ws[:, 1]).all()
    v, si = np.array([1.0, 0.0, 1.0, 1.0, 0.0]), np.array([0, 0, 1, 2, 3])
    est, reps = se.estimates(v, si, np.ones((3, 4), int), "top1@2.0")
    assert est == 0.6 and np.allclose(reps, 0.6)
    est, reps = se.estimates(np.array([3.0, 1.0, 2.0]), np.array([0, 1, 2]), np.array([[1, 1, 1], [0, 3, 0]]), "top1_err")
    assert est == 2.0 and reps.tolist() == [2.0, 1.0]

def test_frozen_params_gate(tmp_path):
    split, tot = tmp_path / "split.tsv", tmp_path / "tot.tsv"
    split.write_text("a\n"); tot.write_text("b\n")
    good = dict(split_md5=se.md5(split), lib_totals_md5=se.md5(tot), tuned_on="dev",
                methods={st: ["M1", "M4@1.0"] for st in se.STRATA})
    def ok(p):
        path = tmp_path / "fp.json"
        path.write_text(p if isinstance(p, str) else json.dumps(p))
        try: se.frozen_params(str(path), str(split), str(tot))
        except ValueError: return False
        return True
    broken = ["{not json", dict(good, tuned_on="test"), dict(good, split_md5="0" * 32),
              {k: v for k, v in good.items() if k != "methods"}, dict(good, methods={"1": ["M1"]})]
    assert_discriminates(ok, [good], broken, "frozen params gate")

def test_split_and_test_summary_refuse_rerun(tmp_path):
    (tmp_path / "dev_test_split.tsv").write_text("x\n")
    with pytest.raises(SystemExit, match="immutable"):
        se.cmd_split(SimpleNamespace(out_dir=str(tmp_path)))
    split, tot = tmp_path / "split.tsv", tmp_path / "tot.tsv"
    split.write_text("set\tsystem_id\tpdb_id\tgroup\tclusters\tsplit\n"); tot.write_text("b\n")
    fp = tmp_path / "fp.json"
    fp.write_text(json.dumps(dict(split_md5=se.md5(split), lib_totals_md5=se.md5(tot), tuned_on="dev",
                                  methods={st: ["M1"] for st in se.STRATA})))
    (tmp_path / "t_results.tsv.DONE").write_text("{}")
    args = SimpleNamespace(split=str(split), extra_systems=[], out_dir=str(tmp_path), name="t", half="test",
                           frozen=str(fp), lib_totals=str(tot), run_dir=str(tmp_path), B=10)
    with pytest.raises(SystemExit, match="evaluated once"):
        se.cmd_summarize(args)
    fp.unlink()
    with pytest.raises(ValueError, match="unreadable"):
        se.cmd_summarize(args)

def test_m8_entity_keys(tmp_path):
    """Altloc duplicates count once, insertion codes count; a gapped chain still matches its entity;
    ambiguous or low-identity matches get their own key; the entity with most vdM residues wins."""
    atom = lambda name, alt, res, ch, num: f"ATOM  {1:5d} {name:<4}{alt}{res} {ch}{num:>5}   " + "0.000   0.000   0.000  1.00  0.00\n"
    pdb = tmp_path / "x.pdb"
    pdb.write_text(atom("CA", "A", "ALA", "A", "1 ") + atom("CA", "B", "ALA", "A", "1 ") + atom("CA", " ", "GLY", "A", "2 ")
                   + atom("CA", " ", "SER", "A", "2A") + atom("CA", " ", "TRP", "B", "5 ")
                   + "HETATM    9  C1  LIG A 900       0.000   0.000   0.000  1.00  0.00\n")
    assert se.chain_seqs(str(pdb)) == {"A": "AGS", "B": "W"}
    al = se._aligner()
    assert se.entity_identity(al, "MKAGSTLLE", "KAGLLE") == 1.0 and se.entity_identity(al, "MKAGSTLLE", "KAGWLE") == pytest.approx(5 / 6)
    rows = [dict(stem="1abc", chain="A", entity="1ABC_1", ident="1.0", second_ident="0.3"),
            dict(stem="1abc", chain="B", entity="1ABC_2", ident="1.0", second_ident="0.97"),
            dict(stem="1abc", chain="C", entity="1ABC_1", ident="0.8", second_ident="0.1")]
    ek = se.entity_keys(rows)
    assert ek == {("1abc", "A"): "1ABC_1", ("1abc", "B"): "1abc:B", ("1abc", "C"): "1abc:C"}
    ek[("2xyz", "A")], ek[("2xyz", "B")] = "2XYZ_1", "2XYZ_2"
    got = se.observation_entities(np.array(["2xyz", "2xyz", "1abc"]), np.array([["A", "B", "B"], ["B", "A", "Z"], ["A", "A", "A"]]), ek)
    assert got.tolist() == ["2XYZ_2", "2XYZ_1", "1ABC_1"]

def test_c2_sampling_fit_and_nesting():
    rng = np.random.default_rng(7)
    n = 3000
    u = dict(bb_rmsd=rng.random(n), min_bb_dist=rng.random(n) * 4, subset_size=np.full(n, 2, np.int8),
             cluster_num_parents=rng.integers(1, 50, n), m8=rng.integers(1, 9, n), m9=rng.integers(1, 5, n),
             placed_cg=rng.normal(size=(n, 3, 3)).astype(np.float32))
    err = np.where(rng.random(n) < 0.02, 1.0, 5.0)  # 2% positives: uniform sampling would keep ~40
    sc = dict(M2=u["cluster_num_parents"].astype(float), M8=u["m8"], M9=u["m9"], **{f"M4@{r}": rng.random(n) for r in se.R_GRID})
    d = se.c2_sample(u, sc, err, "k", np.random.default_rng(0))
    assert d["y"].sum() == (err <= 2).sum() and (~d["y"]).sum() == 1000 and d["w"].sum() == pytest.approx(1.0)
    assert d["w"][d["y"]].sum() == pytest.approx((err <= 2).mean())  # each class keeps its population share
    fit = se.fit_c2(d["X"][:, :4].astype(float), d["y"], d["w"])
    from sklearn.linear_model import LogisticRegression
    mu, sd, coef, b = fit
    ref = LogisticRegression(C=1.0, max_iter=1000).fit((d["X"][:, :4] - mu) / sd, d["y"], sample_weight=d["w"])
    assert np.allclose(coef, ref.coef_[0]) and b == pytest.approx(ref.intercept_[0])
    cfg = dict(r=se.R_GRID[0], C2=fit, C2_M8=fit, C2_M9=fit)
    out = se.c2_scores(u, ((0, 1, 2),), cfg)
    m4 = np.log1p(se.consensus_density(u["placed_cg"], u["cluster_num_parents"], ((0, 1, 2),), se.R_GRID[0]))
    X = np.stack([np.log1p(u["cluster_num_parents"]), -u["bb_rmsd"], u["min_bb_dist"], m4], 1)
    assert np.allclose(out["C2"], ref.decision_function((X - mu) / sd))
    # Swap test: only the support column differs between variants.
    assert not np.allclose(out["C2"], out["C2_M8"]) and np.allclose(out["C2_M8"] - out["C2"], coef[0] * (
        (np.log1p(u["m8"]) - np.log1p(u["cluster_num_parents"])) / sd[0]))

def test_c3_tiebreak_and_nested_rows():
    rows = [dict(method=f"C3@{d0}_{b0}", **{"top1@2.0": "0.5"}) for d0, b0 in se.C3_GRID]
    assert se._best_c3(rows) == "C3@2.2_inf"
    rows[4]["top1@2.0"] = "0.9"
    assert se._best_c3(rows) == f"C3@{se.C3_GRID[4][0]}_{se.C3_GRID[4][1]}"
    cfg = dict(fold_of=dict(S1=0, S2=1), pooled={"0": dict(r=0.5, c3="C3@2.2_inf"), "1": dict(r=2.0, c3="C3@3.0_0.1")})
    units = [dict(system_id=s, stratum="pooled", method=m, v=f"{s}{m}") for s in ("S1", "S2")
             for m in ("M4@0.5", "M4@2.0", "C3@2.2_inf", "C3@3.0_0.1")]
    got = {(r["system_id"], r["method"]): r["v"] for r in se.nested_rows(units, cfg)}
    assert got == {("S1", "M4*"): "S1M4@0.5", ("S2", "M4*"): "S2M4@2.0", ("S1", "C3*"): "S1C3@2.2_inf", ("S2", "C3*"): "S2C3@3.0_0.1"}
