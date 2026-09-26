"""combo_necessity.instance_rows: pruning rules keep/drop the right combos and hits."""
from ligand_vdgs.tools import combo_necessity as cn
from ligand_vdgs.tools import self_recovery as sr
from tests.vacuity import assert_discriminates

F, G = "fragF", "fragG"
A, B, C = ("", "A", 1), ("", "A", 2), ("", "A", 3)
CG = frozenset({"C1", "O1"})

def _hit(frag, residues, cid=1, rmsd=0.1):
    return dict(lig_instance="", frag=frag, subset_size=len(residues), charge_sign="neut",
                aa_bucket="X", vdg_cluster_id=cid, vdg_cluster_num_parents=3,
                vdg_rmsd=f"{rmsd:.4f}", rmsd_threshold="0.6", cg_names=CG,
                residues=frozenset(residues))

def _lig(hits, truth=()):
    E = {frozenset(x) for x in ({A}, {B}, {C}, {A, B}, {A, C}, {B, C})}
    lig = dict(resname="LIG", instances={"": dict(chain="L", resnum=1)},
               enumerated_combos={"": E}, truth=list(truth), contact_classes={}, matched_frags=[F, G])
    H = sr._hits_frame(hits)
    return dict(lig, units=sr.build_units(lig, H), hit_combos=sr.hit_combos(H))

# A productive at size 1 for F only; B productive at size 1 for G only; C never.
HITS = [_hit(F, {A}), _hit(G, {B}), _hit(F, {A, B}, cid=2), _hit(F, {A, C}, cid=3)]
DISTS = {A: 3.0, B: 4.2, C: 5.5}

def _rules(lig=None):
    return {r["rule"]: r for r in cn.instance_rows("s", lig or _lig(HITS), "", DISTS)}

def test_needed_fractions_and_invariant():
    r = _rules()["all"]
    assert (r["needed_1"], r["needed_2"]) == (2 / 3, 2 / 3)
    assert r["hit_recall"] == r["productive_recall"] == 1.0 and r["tasks_pruned"] == 0.0

def test_r1_prunes_combo_with_unproductive_residue():
    # {A,C}: C has no size-1 hit -> pruned, losing 1 of 2 size-2 hits. {A,B} kept.
    r = _rules()["R1"]
    assert r["hit_recall"] == 0.5 and r["tasks_pruned"] == 1 - 2 / 6  # {A,B} x 2 frags kept

def test_r1f_is_per_fragment():
    # For F, B is not size-1 productive (only for G) -> F's {A,B} hit is pruned by R1f, kept by R1.
    rules = _rules()
    assert (rules["R1"]["hit_recall"], rules["R1f"]["hit_recall"]) == (0.5, 0.0)

def test_r2_distance_threshold():
    rules = _rules()
    # max dist: {A,B}=4.2, {A,C}=5.5, {B,C}=5.5
    assert_discriminates(lambda name: rules[name]["hit_recall"] == 1.0,
                         ["R2@6"], ["R2@4.5", "R2@4", "R2@0"], "R2 keeps both size-2 hits")
    assert rules["R2@4.5"]["hit_recall"] == 0.5 and rules["R2@0"]["tasks_pruned"] == 1.0

def test_own_recall_counts_strict_units():
    truth = [dict(frag=F, subset_size=2, charge_sign="neut", aa_bucket="X", cluster_id=3,
                  num_parents=3, is_nr=False, cg_chain="L", cg_resnum=1, cg_names=CG,
                  residues=frozenset({A, C}))]
    rules = _rules(_lig(HITS, truth))
    assert rules["all"]["n_own2"] == 1
    assert (rules["R1"]["own_recall"], rules["R2@6"]["own_recall"]) == (0.0, 1.0)
