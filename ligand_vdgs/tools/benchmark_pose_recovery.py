"""
benchmark_pose_recovery.py

Benchmark the vdG library's ability to recover known ligand chemical group (CG)
positions from crystal structures.

Outputs (in --outdir)
---------------------
hits.tsv       joint mode only: one row per vdG hit
hits.npz       every hit, both modes; pinned contract in NPZ_FIELDS (plan_W4.md). One file per
               structure: --outdir/hits.npz, or --outdir/<name>/hits.npz with --yaml.
fragments.tsv  one row per fragment per structure per subset_scope (all, 1, 2)

hits.tsv columns
    structure, fragment, match_mode, cg_rmsd, held_out_cg_rmsd, bb_rmsd, cg_bb_dist, vdg_rmsd,
    bsr_combo, aa_bucket, charge_sign, subset_size, vdg_cluster_id, vdg_cluster_num_parents,
    vdg_index

fragments.tsv columns
    structure, fragment, match_mode, subset_scope, in_library, n_hits, and for each selection
    rule in (oracle, top_bb, top_support): <rule>_held_out_cg_rmsd, <rule>_cg_bb_dist,
    <rule>_bb_rmsd, <rule>_cluster_num_parents

Quantities (hit finder, DR-36)
------------------------------
- `held_out_cg_rmsd`: vdG CG placed by the backbone-only fit (N/CA/C), RMSD (/N_cg) to the
  nearest crystal CG labeling (`bb` rows: every site x CG automorphism; `joint` rows: the
  automorphisms of the joint fit's `q_site_idx`). The CG never enters the fit.
- `cg_rmsd`: CG-only residual of the joint bb+CG fit -- IN-SAMPLE, bounded by
  `vdg_rmsd * sqrt(n_atoms / N_cg)`. Reference only.
- `cg_bb_dist`: lever arm, mean CG-atom distance to the nearest BSR N/CA/C. Held-out error is
  expected to grow with it.
- Selection rules per fragment: `oracle` = min held_out over hits (picks by the truth; an upper
  bound on recovery), `top_bb` = min bb_rmsd, `top_support` = max vdg_cluster_num_parents
  (ties -> min bb_rmsd). Only the last two are deployable predictions.

Match mode (`--match-mode`)
---------------------------
- `joint`: bb+CG RMSD <= tau. The hit pool is selected for agreeing with the crystal CG, so
  held-out error over it is not a recovery number.
- `bb` (placement): sqrt(SSD_bb / n_atoms) <= tau, i.e. bb_rmsd <= tau * sqrt(n_atoms / N_bb),
  reading neither the query CG nor side chains (the backbone-slot gate is a virtual-CB clash + Pro N
  check on the placed CG); joint hits are a subset at equal tau. SITE-KNOWN: BSR combos are every
  subset of residues with any atom within 4.5 A of the whole crystal ligand
  (dock_utils.get_bsr_combinations; the current BSR definition, like a ligand-defined docking box),
  not the CG's own contacts; no per-CG contact filter.
- tau defaults to the hit finder's derived `normalize_rmsd(n_atoms, 'cgvdmbb')`.

Downstream analysis (pandas):

    import pandas as pd
    frags = pd.read_csv("fragments.tsv", sep="\\t")
    f = frags[frags["in_library"] & (frags["n_hits"] > 0) & (frags["subset_scope"] == "2")]
    for rule in ("oracle", "top_bb", "top_support"):
        ho = pd.to_numeric(f[f"{rule}_held_out_cg_rmsd"])
        print(rule, [(c, float((ho <= c).mean())) for c in (1.0, 1.5, 2.0)])

Usage
-----
Single structure:
    python benchmark_pose_recovery.py \\
        --pdb  crystal.pdb \\
        --smiles "CC(=O)Nc1ccc(cc1)O" \\
        --vdg-lib-dir /path/to/frag_lib \\
        --outdir results/

Multiple structures (YAML):
    python benchmark_pose_recovery.py --yaml config.yml --outdir results/

YAML format:
    structures:
      - pdb: /path/to/crystal.pdb
        smiles: "CC(=O)Nc1ccc(cc1)O"
        name: my_structure          # optional
    vdg_lib_dir: /path/to/frag_lib
    search_threshold: 1.0           # optional; Å (default: derived per fragment)
"""

import argparse
import csv
import os
import yaml
from pathlib import Path

import numpy as np
import prody
from scipy.spatial import cKDTree

from ligand_vdgs.functions.vdg_struct_utils import receptor_bb_vcb_coords
from ligand_vdgs.score_poses.hit_finder_core import MATCH_MODES, init_worker, lib_entries, score_one_model_multi_instance

def _fmt(x): return "NA" if x is None else f"{x:.4f}"

# Pinned hits.npz contract (private/project_planning/plan_W4.md; W4 consumes it). nr_idx indexes the
# (frag, subset_size, charge_sign, aa_bucket) bucket. vdg_rmsd is NaN in bb mode: that fit reads the
# crystal CG. placed_cg: vdG CG placed by the bb-only fit, atoms in the library's stored CG order;
# q_lig_atom_idx[j]: atom of ligand instance `lig_instance` (H-free mol from
# ligand_structure.get_query_ligand_mol) matched to placed atom j under the automorphism chosen by
# the crystal CG (evaluation only). placed_cg_element: library CG elements. min_bb_dist: min over placed
# atoms of the distance to the nearest receptor N/CA/C/O or virtual CB (no side chains).
NPZ_FIELDS = ("frag", "subset_size", "match_mode", "charge_sign", "aa_bucket", "nr_idx", "bb_rmsd",
              "vdg_rmsd", "held_out_cg_rmsd", "cg_bb_dist", "cluster_num_parents", "nr_parent_biounit",
              "lig_instance", "min_bb_dist", "placed_cg", "q_lig_atom_idx", "placed_cg_element")
_PADDED = dict(placed_cg=(np.nan, np.float32), q_lig_atom_idx=(-1, np.int32), placed_cg_element=("", "U2"))
RULES = ("oracle", "top_bb", "top_support")
SCOPES = ("all", "1", "2")

def _hit_arrays(records, receptor_tree):
    """frag -> dict of per-hit NPZ_FIELDS arrays, from compact records. receptor_tree: cKDTree over
    vdg_struct_utils.receptor_bb_vcb_coords."""
    out = {}
    for rec in records:
        (n, m, _), bb = rec["placed_cg"].shape, rec["match_mode"] == "bb"
        cols = dict(frag=np.full(n, rec["frag"]), subset_size=np.full(n, int(rec["subset_size"]), np.int8),
                    match_mode=np.full(n, rec["match_mode"]), charge_sign=np.asarray(rec["charge_sign"], str),
                    aa_bucket=np.full(n, rec["aa_bucket"]), nr_idx=np.asarray(rec["vdg_index"], np.int32),
                    bb_rmsd=rec["bb_rmsd"], vdg_rmsd=np.full(n, np.nan, np.float32) if bb else rec["vdg_rmsd"],
                    held_out_cg_rmsd=rec["held_out_cg_rmsd"], cg_bb_dist=rec["cg_bb_dist"],
                    cluster_num_parents=np.asarray(rec["vdg_cluster_num_parents"], np.int32),
                    nr_parent_biounit=np.asarray(rec["nr_parent_biounit"], str),
                    lig_instance=np.full(n, rec["lig_instance"]),
                    min_bb_dist=receptor_tree.query(rec["placed_cg"].reshape(-1, 3))[0].reshape(n, m).min(
                        axis=1, initial=np.inf).astype(np.float32),
                    placed_cg=rec["placed_cg"], q_lig_atom_idx=rec["q_lig_atom_idx"],
                    placed_cg_element=np.asarray(rec["placed_cg_element"], "U2"))
        d = out.setdefault(rec["frag"], {k: [] for k in NPZ_FIELDS})
        for k, v in cols.items(): d[k].append(v)
    return {f: {k: np.concatenate(v) for k, v in d.items()} for f, d in out.items()}

def write_hits_npz(path, arrays):
    """One row per hit (NPZ_FIELDS); per-atom fields padded to max_N_cg (`_PADDED`: placed_cg NaN,
    q_lig_atom_idx -1, placed_cg_element ''). Zero hits -> empty arrays."""
    hs = list(arrays.values())
    n, m = sum(len(h["bb_rmsd"]) for h in hs), max((h["placed_cg"].shape[1] for h in hs), default=0)
    padded = {k: np.full((n, m, 3) if k == "placed_cg" else (n, m), fill, dt) for k, (fill, dt) in _PADDED.items()}
    i = 0
    for h in hs:
        j, c = len(h["bb_rmsd"]), h["placed_cg"].shape[1]
        for k in _PADDED: padded[k][i:i + j, :c] = h[k]
        i += j
    empty = dict(subset_size=np.int8, nr_idx=np.int32, cluster_num_parents=np.int32)
    cols = {k: np.concatenate([h[k] for h in hs]) if hs else
            np.zeros(0, empty.get(k, np.float32 if k.endswith(("rmsd", "dist")) else "U1"))
            for k in NPZ_FIELDS if k not in _PADDED}
    np.savez_compressed(path, **cols, **padded)

def _select(h):
    """Index per rule. Oracle picks by the truth; top_bb (ties: more support) and top_support
    (ties: lower bb_rmsd) are CG-free."""
    return dict(oracle=int(np.argmin(h["held_out_cg_rmsd"])),
                top_bb=int(np.lexsort((-h["cluster_num_parents"], h["bb_rmsd"]))[0]),
                top_support=int(np.lexsort((h["bb_rmsd"], -h["cluster_num_parents"]))[0]))

def analyze_structure(pdb_path, lig_smiles, name, vdg_lib_dir, vdg_lib_entries,
                      search_threshold, match_mode="joint"):
    """Run hit-finding (metrics range over every hit). Returns (hit_rows,
    fragment_rows, hit_arrays); hit_rows only in joint mode (placement hit counts reach 1e7)."""
    print(f"\n{'='*60}")
    print(f"Structure: {name}  PDB: {pdb_path}  SMILES: {lig_smiles}  match_mode={match_mode}")
    print(f"{'='*60}")
    log_text, match_records, frags_in_lib, _, _ = score_one_model_multi_instance(
        pdbfile=os.path.basename(pdb_path), pdb_path=pdb_path, lig_smiles=lig_smiles,
        vdg_lib_dir=vdg_lib_dir, rmsd_threshold=search_threshold,
        vdg_lib_entries=vdg_lib_entries, match_mode=match_mode, compact=True)
    if log_text.strip():
        print(log_text)
    hit_rows = [dict(structure=name, fragment=rec["frag"], match_mode="joint", cg_rmsd=f"{c:.4f}",
                     held_out_cg_rmsd=f"{h:.4f}", bb_rmsd=f"{b:.4f}", cg_bb_dist=f"{lv:.4f}",
                     vdg_rmsd=f"{v:.4f}", bsr_combo=rec["bsr_combo"], aa_bucket=rec["aa_bucket"],
                     charge_sign=s, subset_size=rec["subset_size"], vdg_cluster_id=int(ci),
                     vdg_cluster_num_parents=int(npar), vdg_index=int(vi))
                for rec in match_records if rec["match_mode"] == "joint"
                for c, h, b, lv, v, s, ci, npar, vi in zip(
                    rec["in_sample_cg_rmsd"], rec["held_out_cg_rmsd"], rec["bb_rmsd"], rec["cg_bb_dist"],
                    rec["vdg_rmsd"], rec["charge_sign"], rec["vdg_cluster_id"],
                    rec["vdg_cluster_num_parents"], rec["vdg_index"])]
    arrays = _hit_arrays(match_records, cKDTree(receptor_bb_vcb_coords(prody.parsePDB(pdb_path))))
    fragment_rows = []
    for frag_smiles, in_lib in frags_in_lib.items():
        h_all = arrays.get(frag_smiles)
        for scope in SCOPES:
            h = (None if h_all is None else h_all if scope == "all" else
                 {k: v[h_all["subset_size"] == int(scope)] for k, v in h_all.items()})
            n = 0 if h is None else len(h["bb_rmsd"])
            row = dict(structure=name, fragment=frag_smiles, match_mode=match_mode,
                       subset_scope=scope, in_library=in_lib, n_hits=n)
            picks = _select(h) if n else {}
            for rule in RULES:
                i = picks.get(rule)
                row.update({f"{rule}_{k}": "NA" if i is None else
                            (int(h[k][i]) if k == "cluster_num_parents" else f"{h[k][i]:.4f}")
                            for k in ("held_out_cg_rmsd", "cg_bb_dist", "bb_rmsd", "cluster_num_parents")})
            fragment_rows.append(row)
            if scope == "all":
                status = ("NOT-IN-LIB" if not in_lib else "NO-HIT" if not n else
                          " ".join(f"{r}={row[f'{r}_held_out_cg_rmsd']}Å" for r in RULES))
                print(f"  {frag_smiles:35s}  {status}  n_hits={n}")
    return hit_rows, fragment_rows, arrays

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark vdG library CG placement vs crystal structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--yaml", metavar="CONFIG", help="YAML config (batch mode)")
    parser.add_argument("--pdb")
    parser.add_argument("--smiles")
    parser.add_argument("--name")
    parser.add_argument("--vdg-lib-dir", dest="vdg_lib_dir")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--search-threshold", type=float, default=None,
                        dest="search_threshold", metavar="Å",
                        help="tau (Å); default: derived per fragment via normalize_rmsd, as in "
                             "the hit finder.")
    parser.add_argument("--match-mode", choices=MATCH_MODES, default="joint",
                        help="joint (bb+CG) or bb (DR-36 placement mode).")
    args = parser.parse_args()

    if args.yaml:
        with open(args.yaml) as f:
            cfg = yaml.safe_load(f)
        if not isinstance(cfg, dict):
            parser.error(f"{args.yaml} is empty or not a YAML mapping")
        for key in ("structures", "vdg_lib_dir"):
            if key not in cfg:
                parser.error(f"{args.yaml} is missing required key: {key}")
        structures       = cfg["structures"]
        vdg_lib_dir      = cfg["vdg_lib_dir"]
        search_threshold = cfg.get("search_threshold")
    else:
        if not (args.pdb and args.smiles and args.vdg_lib_dir):
            parser.error("--pdb, --smiles, and --vdg-lib-dir are required without --yaml")
        structures       = [{"pdb": args.pdb, "smiles": args.smiles,
                              "name": args.name or Path(args.pdb).stem}]
        vdg_lib_dir      = args.vdg_lib_dir
        search_threshold = args.search_threshold

    if not os.path.isdir(vdg_lib_dir):
        parser.error(f"vdg_lib_dir not found: {vdg_lib_dir}")

    os.makedirs(args.outdir, exist_ok=True)
    vdg_lib_entries = lib_entries(vdg_lib_dir)
    # Eagerly populate the module-level mol cache. Direct core-library calls can
    # initialize it lazily, but doing it once here avoids work during analysis.
    init_worker(vdg_lib_entries)

    all_hit_rows      = []
    all_fragment_rows = []

    for struct_cfg in structures:
        hit_rows, fragment_rows, arrays = analyze_structure(
            pdb_path         = struct_cfg["pdb"],
            lig_smiles       = struct_cfg["smiles"],
            name             = struct_cfg.get("name", Path(struct_cfg["pdb"]).stem),
            vdg_lib_dir      = vdg_lib_dir,
            vdg_lib_entries  = vdg_lib_entries,
            search_threshold = search_threshold,
            match_mode        = args.match_mode,
        )
        all_hit_rows.extend(hit_rows)
        all_fragment_rows.extend(fragment_rows)
        name = struct_cfg.get("name", Path(struct_cfg["pdb"]).stem)
        npz_dir = args.outdir if len(structures) == 1 else os.path.join(args.outdir, name)
        os.makedirs(npz_dir, exist_ok=True)
        write_hits_npz(os.path.join(npz_dir, "hits.npz"), arrays)
        del arrays

    def _write_tsv(rows, path):
        if not rows:
            return
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"  → {path}")

    if not all_hit_rows:
        print("  No vdG hits found; hits.tsv not written.")
    if not all_fragment_rows:
        print("  No fragments analysed; fragments.tsv not written.")
    _write_tsv(all_hit_rows,      os.path.join(args.outdir, "hits.tsv"))
    _write_tsv(all_fragment_rows, os.path.join(args.outdir, "fragments.tsv"))
    print("VERDICT benchmark_pose_recovery complete", flush=True)


if __name__ == "__main__":
    main()
