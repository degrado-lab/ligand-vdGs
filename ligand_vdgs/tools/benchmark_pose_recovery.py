"""
benchmark_pose_recovery.py

Benchmark the vdG library's ability to recover known ligand chemical group (CG)
positions from crystal structures.

Outputs (in --outdir)
---------------------
hits.tsv       one row per vdG hit (all hits within --search-threshold)
fragments.tsv  one row per fragment per structure

hits.tsv columns
    structure, fragment, cg_rmsd, vdg_rmsd, bsr_combo, aa_bucket,
    subset_size, vdg_cluster_id, vdg_cluster_size, vdg_index

fragments.tsv columns
    structure, fragment, in_library, n_hits, best_cg_rmsd, best_vdg_rmsd

Downstream analysis (pandas):

    import pandas as pd
    hits  = pd.read_csv("hits.tsv",      sep="\\t")
    frags = pd.read_csv("fragments.tsv", sep="\\t")

    # Recovery rate at any threshold using the best hit per fragment:
    threshold = 1.0  # Å
    best = hits.groupby(["structure", "fragment"])["cg_rmsd"].min()
    print(f"Recovery @ {threshold} Å: {(best <= threshold).mean():.1%}")

    # CDF of best-hit cg_rmsd across all in-library fragments:
    import matplotlib.pyplot as plt
    in_lib = frags[frags["in_library"]]["best_cg_rmsd"].replace("NA", None)
    in_lib = pd.to_numeric(in_lib, errors="coerce").dropna().sort_values()
    plt.plot(in_lib.values, np.linspace(0, 1, len(in_lib)))
    plt.xlabel("best cg_rmsd (Å)"); plt.ylabel("fraction"); plt.show()

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
    search_threshold: 2.0           # optional; Å (default 2.0)
"""

import argparse
import csv
import os
import yaml
from pathlib import Path

import numpy as np

from ligand_vdgs.functions import dock_utils as dock
from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.score_poses.hit_finder_core import (
    cg_atom_order_smarts,
    init_worker,
    score_one_model,
)


def _reconstruct_R_t(rec):
    R = np.array([[float(rec[f"R{i}{j}"]) for j in range(3)] for i in range(3)],
                 dtype=np.float32)
    t = np.array([float(rec[f"t{k}"]) for k in range(3)], dtype=np.float32)
    return R, t


def _build_crystal_cg_cache(filtered_frags, vdg_lib_dir):
    """
    Return crystal CG coords keyed by (db_name, site_idx, perm_idx).

    Takes the post-mapped filtered_frags from score_one_model to ensure site_idx
    and perm_idx align with hit records (q_site_idx / q_cg_perm_idx). This avoids
    redundant PDB parsing, fragmentation, and library mapping.
    """
    cache = {}
    for db_name, grouped_sites in filtered_frags.items():
        # The mols were reordered to the library's recorded CG atom order, which
        # is not necessarily the library directory name.
        order_smarts = cg_atom_order_smarts(vdg_lib_dir, db_name)
        cache[db_name] = {}
        for site_idx, site in enumerate(grouped_sites):
            cache[db_name][site_idx] = {}
            for perm_idx, (sub_mol, _perm_inds, _orig_mol_inds) in enumerate(site):
                try:
                    cache[db_name][site_idx][perm_idx] = np.asarray(
                        dock.get_query_cg_coords(sub_mol, order_smarts), np.float32)
                except Exception as e:
                    print(f"  [WARNING] CG coords unavailable for {db_name} "
                          f"site {site_idx} perm {perm_idx}: {e}")
    return cache


def compute_cg_placement_rmsd(hit_rec, crystal_cg_cache, vdg_lib_dir, bucket_cache):
    """
    Apply the hit's R/t to the library CG and return RMSD vs the crystal CG (Å).

    Caveat: the crystal CG coords are part of the target that R/t were fitted to
    (score_one_model aligns backbone + query CG jointly), so this is the CG-only
    residual of that fit and is bounded by vdg_rmsd -- not an independent
    measurement of pose recovery.
    """
    frag      = hit_rec["frag"]
    subset    = int(hit_rec["subset_size"])
    aa_bucket = hit_rec["aa_bucket"]
    vdg_idx   = int(hit_rec["vdg_index"])
    site_idx  = int(hit_rec["q_site_idx"])
    perm_idx  = int(hit_rec["q_cg_perm_idx"])

    bucket_key = (frag, subset, aa_bucket)
    if bucket_key not in bucket_cache:
        bucket_cache[bucket_key] = vdg_npz.load_vdg_bucket(
            vdg_lib_dir, frag, subset, aa_bucket)
    bucket = bucket_cache[bucket_key]
    if bucket is None or vdg_idx >= bucket["cg"].shape[0]:
        return None

    lib_cg     = bucket["cg"][vdg_idx].astype(np.float32)
    crystal_cg = crystal_cg_cache.get(frag, {}).get(site_idx, {}).get(perm_idx)
    if crystal_cg is None or lib_cg.shape != crystal_cg.shape:
        return None

    R, t = _reconstruct_R_t(hit_rec)
    diff = lib_cg @ R + t - crystal_cg
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def analyze_structure(pdb_path, lig_smiles, name, vdg_lib_dir, vdg_lib_entries,
                      search_threshold):
    """Run hit-finding and return (hit_rows, fragment_rows)."""
    print(f"\n{'='*60}")
    print(f"Structure: {name}  PDB: {pdb_path}  SMILES: {lig_smiles}")
    print(f"{'='*60}")

    log_text, match_records, frags_in_lib, _, filtered_frags = score_one_model(
        pdbfile=os.path.basename(pdb_path),
        pdb_path=pdb_path,
        lig_smiles=lig_smiles,
        vdg_lib_dir=vdg_lib_dir,
        rmsd_threshold=search_threshold,
        vdg_lib_entries=vdg_lib_entries,
        deduplicate=True,
    )
    if log_text.strip():
        print(log_text)

    try:
        crystal_cg_cache = _build_crystal_cg_cache(filtered_frags, vdg_lib_dir)
    except Exception as e:
        print(f"  [WARNING] Could not extract crystal CG coords: {e}")
        crystal_cg_cache = {}

    bucket_cache = {}
    hits_with_rmsd = [
        (rec, compute_cg_placement_rmsd(rec, crystal_cg_cache, vdg_lib_dir, bucket_cache))
        for rec in match_records
    ]

    hit_rows = [
        dict(
            structure=name,
            fragment=rec["frag"],
            cg_rmsd=f"{cg_rmsd:.4f}" if cg_rmsd is not None else "NA",
            vdg_rmsd=f"{float(rec['vdg_rmsd']):.4f}",
            bsr_combo=rec["bsr_combo"],
            aa_bucket=rec["aa_bucket"],
            subset_size=rec["subset_size"],
            vdg_cluster_id=rec["vdg_cluster_id"],
            vdg_cluster_size=rec["vdg_cluster_size"],
            vdg_index=rec["vdg_index"],
        )
        for rec, cg_rmsd in hits_with_rmsd
    ]

    frag_hits = {}
    for rec, cg_rmsd in hits_with_rmsd:
        frag_hits.setdefault(rec["frag"], []).append((rec, cg_rmsd))

    fragment_rows = []
    for frag_smiles, in_lib in frags_in_lib.items():
        frag_hit_pairs = frag_hits.get(frag_smiles, [])
        n_hits = len(frag_hit_pairs)
        if frag_hit_pairs:
            best_vdg = min(float(rec["vdg_rmsd"]) for rec, _ in frag_hit_pairs)
            cg_values = [cg for _, cg in frag_hit_pairs if cg is not None]
            best_cg = min(cg_values) if cg_values else None
        else:
            best_vdg = None
            best_cg  = None

        fragment_rows.append(dict(
            structure=name,
            fragment=frag_smiles,
            in_library=in_lib,
            n_hits=n_hits,
            best_cg_rmsd=f"{best_cg:.4f}" if best_cg is not None else "NA",
            best_vdg_rmsd=f"{best_vdg:.4f}" if best_vdg is not None else "NA",
        ))

        if not in_lib:
            status = "NOT-IN-LIB"
        elif n_hits == 0:
            status = "NO-HIT"
        elif best_cg is not None:
            status = f"cg={best_cg:.3f}Å"
        else:
            status = f"vdg={best_vdg:.3f}Å"
        print(f"  {frag_smiles:35s}  {status}  n_hits={n_hits}")

    return hit_rows, fragment_rows


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
    parser.add_argument("--search-threshold", type=float, default=2.0,
                        dest="search_threshold", metavar="Å")
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
        search_threshold = cfg.get("search_threshold", 2.0)
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
    vdg_lib_entries = set(os.listdir(vdg_lib_dir))
    # Eagerly populate the module-level mol cache. Direct core-library calls can
    # initialize it lazily, but doing it once here avoids work during analysis.
    init_worker(vdg_lib_entries)

    all_hit_rows      = []
    all_fragment_rows = []

    for struct_cfg in structures:
        hit_rows, fragment_rows = analyze_structure(
            pdb_path         = struct_cfg["pdb"],
            lig_smiles       = struct_cfg["smiles"],
            name             = struct_cfg.get("name", Path(struct_cfg["pdb"]).stem),
            vdg_lib_dir      = vdg_lib_dir,
            vdg_lib_entries  = vdg_lib_entries,
            search_threshold = search_threshold,
        )
        all_hit_rows.extend(hit_rows)
        all_fragment_rows.extend(fragment_rows)

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


if __name__ == "__main__":
    main()
