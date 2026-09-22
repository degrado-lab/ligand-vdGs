"""
benchmark_pose_recovery.py

Benchmark the vdG library's ability to recover known ligand chemical group (CG)
positions from crystal structures.

Outputs (in --outdir)
---------------------
hits.tsv       one row per vdG hit (all hits within --search-threshold)
fragments.tsv  one row per fragment per structure

hits.tsv columns
    structure, fragment, cg_rmsd, held_out_cg_rmsd, vdg_rmsd, bsr_combo,
    aa_bucket, charge_sign, subset_size, vdg_cluster_id, vdg_cluster_num_parents, vdg_index

fragments.tsv columns
    structure, fragment, in_library, n_hits, best_cg_rmsd,
    best_held_out_cg_rmsd, best_vdg_rmsd

WHAT THIS FILE DOES NOT MEASURE
-------------------------------
Do not report a "recovery rate" from these outputs. Two separate circularities:

1. `cg_rmsd` is IN-SAMPLE. `score_one_model` fits R/t over backbone + query CG
   jointly, so the crystal CG is part of the target the transform was fitted to.
   `cg_rmsd` is that fit's CG-only residual, bounded by construction at
   `vdg_rmsd * sqrt(n_atoms / N_cg)` -- NOT by `vdg_rmsd` itself, since
   `vdg_rmsd` normalizes the joint residual by n_atoms = N_bb + N_cg while
   `cg_rmsd` normalizes by N_cg alone. So `cg_rmsd` exceeding `vdg_rmsd`
   (or the search threshold) is expected, not evidence of a leak.
   `held_out_cg_rmsd` removes this circularity: it applies the BACKBONE-ONLY
   transform (Rbb/tbb), so the CG never enters the objective.

2. SELECTION is still CG-dependent, and `held_out_cg_rmsd` does NOT fix it. A
   vdG becomes a hit only if its joint backbone+CG RMSD clears the threshold, so
   the hit pool was chosen for agreeing with the crystal CG. A min over hits
   therefore answers "among library vdGs already known to match the crystal CG,
   how well does the best one match it" -- not "can the library place this CG".
   An honest recovery number needs a search that never sees the query CG; that
   does not exist yet (see DR-17).

So: `held_out_cg_rmsd` is CG placement error under backbone-only alignment,
among CG-selected hits. It does not stand in for a recovery number in EITHER
direction.

Name the quantity, not just the direction: in RMSD this column is >= an oracle's,
which in RECOVERY-RATE terms makes it a lower bound. Same fact, opposite words --
that inversion has already caused one wrong claim.

The subset argument, stated with the threshold it actually holds at. Search
thresholds a joint RMSD normalized by n_atoms = N_bb + N_cg (hit_finder_core
:865, :1048), while a backbone-only search would normalize by N_bb. So a
joint pass at tau gives

    SSD_bb^bbfit <= SSD_bb^joint <= n_atoms * tau^2
    =>  RMSD_bb^bbfit <= tau * sqrt(n_atoms / N_bb)

i.e. the joint-selected pool is a subset of what a backbone-only search returns
only at the RELAXED threshold tau * sqrt(n_atoms/N_bb), not at tau itself. That
factor is ~1.4-1.7 at subset_size 1 (3-4 backbone atoms against a 4-5 atom CG),
so it is not negligible. Do not quote the subset claim without the threshold.

Against a deployable top-k-by-support prediction it is bounded in no proven
direction at all.

Downstream analysis (pandas):

    import pandas as pd
    hits  = pd.read_csv("hits.tsv",      sep="\\t")
    frags = pd.read_csv("fragments.tsv", sep="\\t")

    # Fit quality (NOT recovery): CDF of vdg_rmsd, the joint bb+CG residual.
    import matplotlib.pyplot as plt
    vq = hits["vdg_rmsd"].sort_values()
    plt.plot(vq.values, np.linspace(0, 1, len(vq)))
    plt.xlabel("vdg_rmsd (Å) -- fit quality, not recovery")
    plt.ylabel("fraction"); plt.show()

    # CG placement error under bb-only alignment, among CG-selected hits.
    # Label any figure from this with that whole phrase; it is not recovery.
    ho = pd.to_numeric(frags[frags["in_library"]]["best_held_out_cg_rmsd"],
                       errors="coerce").dropna().sort_values()
    plt.plot(ho.values, np.linspace(0, 1, len(ho)))
    plt.xlabel("best held_out_cg_rmsd (Å)"); plt.ylabel("fraction"); plt.show()

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
    score_one_model_multi_instance,
)


def _reconstruct_R_t(rec, prefix=""):
    """prefix="" -> joint bb+CG fit; prefix="bb" -> backbone-only fit."""
    R = np.array([[float(rec[f"R{prefix}{i}{j}"]) for j in range(3)] for i in range(3)],
                 dtype=np.float32)
    t = np.array([float(rec[f"t{prefix}{k}"]) for k in range(3)], dtype=np.float32)
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
    residual of that fit, not an independent measurement of pose recovery.
    It is bounded by vdg_rmsd * sqrt(n_atoms / N_cg), NOT by vdg_rmsd itself:
    the normalizers differ (n_atoms = N_bb + N_cg vs N_cg alone).
    Use `compute_held_out_cg_rmsd` instead; this
    value is kept only as the in-sample reference, and is scored at the stored
    `q_cg_perm_idx` so it is not pointwise comparable to the held-out one.
    """
    frag      = hit_rec["frag"]
    subset    = int(hit_rec["subset_size"])
    aa_bucket = hit_rec["aa_bucket"]
    vdg_idx   = int(hit_rec["vdg_index"])
    site_idx  = int(hit_rec["q_site_idx"])
    perm_idx  = int(hit_rec["q_cg_perm_idx"])

    charge_sign = hit_rec["charge_sign"]
    bucket_key = (frag, subset, charge_sign, aa_bucket)
    if bucket_key not in bucket_cache:
        bucket_cache[bucket_key] = vdg_npz.load_vdg_bucket(
            vdg_lib_dir, frag, subset, charge_sign, aa_bucket)
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


# Rounding an exact proper rotation to 4dp reaches |det-1| ~ 1.8e-4 and
# ortho_err ~ 1.6e-4 by itself; these bounds are that envelope, not a quality bar.
_ROT_DET_TOL, _ROT_ORTHO_TOL = 3e-4, 3e-4

def _is_proper_rotation(R):
    return (abs(float(np.linalg.det(R)) - 1.0) <= _ROT_DET_TOL
            and float(np.max(np.abs(R.T @ R - np.eye(3)))) <= _ROT_ORTHO_TOL)


def _load_hit_bucket(hit_rec, vdg_lib_dir, bucket_cache):
    """Return the hit's library CG coords, or None if the bucket is unavailable."""
    frag = hit_rec["frag"]
    key = (frag, int(hit_rec["subset_size"]), hit_rec["charge_sign"],
           hit_rec["aa_bucket"])
    if key not in bucket_cache:
        bucket_cache[key] = vdg_npz.load_vdg_bucket(
            vdg_lib_dir, frag, key[1], key[2], key[3])
    bucket = bucket_cache[key]
    vdg_idx = int(hit_rec["vdg_index"])
    if bucket is None or vdg_idx >= bucket["cg"].shape[0]: return None
    return bucket["cg"][vdg_idx].astype(np.float32)


def compute_held_out_cg_rmsd(hit_rec, crystal_cg_cache, vdg_lib_dir, bucket_cache):
    """
    CG placement error (Å) with the CG HELD OUT of the superposition.

    Applies the hit's BACKBONE-ONLY transform (Rbb/tbb from score_one_model) to
    the library CG and scores it against the crystal CG. Because the CG never
    enters the Kabsch objective, this is a prediction error rather than the
    in-sample residual `compute_cg_placement_rmsd` returns.

    Minimizes over every CG permutation recorded for the hit's site -- that set
    is the CG automorphism group, so a symmetric CG (nitro, carboxylate,
    tetrazole) is not charged for a relabeling. The hit's own `q_cg_perm_idx` is
    deliberately NOT used: the joint fit chose that perm with the crystal CG in
    the objective, so reusing it would leak the held-out quantity back in.

    Only the SUPERPOSITION is held out. The correspondence it is applied to --
    aa_perm_idx (library slot <-> query residue), vdg_idx, q_site_idx -- was
    still chosen by the joint fit, with the crystal CG in the objective. A
    residual leak, not addressed here.

    NOT pointwise comparable to `cg_rmsd`. At a FIXED correspondence the joint
    fit minimizes SSD_bb + SSD_cg while this one minimizes SSD_bb alone, so
    SSD_cg(joint) <= SSD_cg(bb-only) always. That ordering is only guaranteed at
    the same perm: `cg_rmsd` is scored at the stored `q_cg_perm_idx` while this
    minimizes over the automorphism group, so this value can land BELOW
    `cg_rmsd` on a symmetric CG. That is the comparison being ill-posed, not a
    bug. Compare the two only at a fixed perm, or min over perms on both sides.

    NOT a recovery metric: the hit pool itself was selected on joint bb+CG RMSD.
    See the "WHAT THIS FILE DOES NOT MEASURE" section of the module docstring.
    """
    if "Rbb00" not in hit_rec: return None
    lib_cg = _load_hit_bucket(hit_rec, vdg_lib_dir, bucket_cache)
    if lib_cg is None: return None

    site_perms = crystal_cg_cache.get(hit_rec["frag"], {}).get(int(hit_rec["q_site_idx"]))
    if not site_perms: return None

    R, t = _reconstruct_R_t(hit_rec, prefix="bb")
    # _is_proper_rotation only catches a corrupted/malformed R (bad det or
    # non-orthogonal) -- e.g. TSV read/write corruption. It CANNOT catch a
    # genuinely rank-deficient backbone fit (collinear or near-collinear
    # points): utils.kabsch's reflection correction (_proper_kabsch_rotations)
    # always returns a proper, orthogonal rotation by construction, even when
    # the fit is non-unique, so a degenerate fit passes this guard silently.
    # Tolerance is the 4-decimal STORAGE rounding envelope: rounding an exact
    # proper rotation to 4dp reaches |det-1| ~ 1.8e-4 on its own (measured over
    # 20k random rotations), so a tighter bound would reject valid fits.
    if not _is_proper_rotation(R): return None
    placed = lib_cg @ R + t
    best = None
    for crystal_cg in site_perms.values():
        if crystal_cg is None or crystal_cg.shape != placed.shape: continue
        diff = placed - crystal_cg
        rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
        if best is None or rmsd < best: best = rmsd
    return best


def analyze_structure(pdb_path, lig_smiles, name, vdg_lib_dir, vdg_lib_entries,
                      search_threshold):
    """Run hit-finding and return (hit_rows, fragment_rows)."""
    print(f"\n{'='*60}")
    print(f"Structure: {name}  PDB: {pdb_path}  SMILES: {lig_smiles}")
    print(f"{'='*60}")

    log_text, match_records, frags_in_lib, _, filtered_frags = score_one_model_multi_instance(
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
        (rec,
         compute_cg_placement_rmsd(rec, crystal_cg_cache, vdg_lib_dir, bucket_cache),
         compute_held_out_cg_rmsd(rec, crystal_cg_cache, vdg_lib_dir, bucket_cache))
        for rec in match_records
    ]

    hit_rows = [
        dict(
            structure=name,
            fragment=rec["frag"],
            cg_rmsd=f"{cg_rmsd:.4f}" if cg_rmsd is not None else "NA",
            held_out_cg_rmsd=f"{ho:.4f}" if ho is not None else "NA",
            vdg_rmsd=f"{float(rec['vdg_rmsd']):.4f}",
            bsr_combo=rec["bsr_combo"],
            aa_bucket=rec["aa_bucket"],
            charge_sign=rec["charge_sign"],
            subset_size=rec["subset_size"],
            vdg_cluster_id=rec["vdg_cluster_id"],
            vdg_cluster_num_parents=rec["vdg_cluster_num_parents"],
            vdg_index=rec["vdg_index"],
        )
        for rec, cg_rmsd, ho in hits_with_rmsd
    ]

    frag_hits = {}
    for rec, cg_rmsd, ho in hits_with_rmsd:
        frag_hits.setdefault(rec["frag"], []).append((rec, cg_rmsd, ho))

    fragment_rows = []
    for frag_smiles, in_lib in frags_in_lib.items():
        frag_hit_pairs = frag_hits.get(frag_smiles, [])
        n_hits = len(frag_hit_pairs)
        if frag_hit_pairs:
            best_vdg = min(float(rec["vdg_rmsd"]) for rec, _, _ in frag_hit_pairs)
            cg_values = [cg for _, cg, _ in frag_hit_pairs if cg is not None]
            ho_values = [ho for _, _, ho in frag_hit_pairs if ho is not None]
            best_cg = min(cg_values) if cg_values else None
            best_ho = min(ho_values) if ho_values else None
        else:
            best_vdg = None
            best_cg  = None
            best_ho  = None

        fragment_rows.append(dict(
            structure=name,
            fragment=frag_smiles,
            in_library=in_lib,
            n_hits=n_hits,
            best_cg_rmsd=f"{best_cg:.4f}" if best_cg is not None else "NA",
            best_held_out_cg_rmsd=f"{best_ho:.4f}" if best_ho is not None else "NA",
            best_vdg_rmsd=f"{best_vdg:.4f}" if best_vdg is not None else "NA",
        ))

        if not in_lib:
            status = "NOT-IN-LIB"
        elif n_hits == 0:
            status = "NO-HIT"
        elif best_ho is not None:
            # Held-out is the honest per-fragment number; cg= is in-sample.
            status = f"held-out cg={best_ho:.3f}Å"
        elif best_cg is not None:
            status = f"cg(in-sample)={best_cg:.3f}Å"
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
