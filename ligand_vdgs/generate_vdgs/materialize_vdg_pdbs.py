"""Materialize centroid or member PDBs from current-format ``nr_vdgs`` buckets.

Centroids share one frame; members are re-derived below
``<subset>/<charge_sign>/<bucket>/clus<id>_size<size>``.
"""

import argparse
import os

import numpy as np

from ligand_vdgs.functions import utils
from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions.Frags import check_vdg_job_status
from ligand_vdgs.functions.vdg_pdb_io import (
    NameRegistry, format_residue_tag, fresh_dir, sanitize, strip_known_exts,
    write_pdb_gz)

DEFAULT_MIN_CLUSTER_SIZE = 10
DEFAULT_REPS = 20
DEFAULT_ALIGN_CG_WEIGHT = 0.99
DEFAULT_SEED = 0

PARENT_READING_FLAGS = ("--include-full-ligand", "--carbonyl", "--sidechain",
                        "--members")

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--cg-nr-vdgs-root", required=True, help="Root of nr_vdgs.")
    p.add_argument("-o", "--out-dir", required=True, help="Output directory.")
    p.add_argument("--include-full-ligand", action="store_true", help="Include full ligand.")
    p.add_argument("--carbonyl", action="store_true",
                   help="Also write each vdM's backbone carbonyl O, re-derived from "
                        "the parent PDB. Display only: hit-finding RMSD uses the "
                        "stored N/CA/C.")
    p.add_argument("--sidechain", action=argparse.BooleanOptionalAction, default=True,
                   help="Write each vdM's sidechain heavy atoms, re-derived from the "
                        "parent PDB. On by default, but falls back to CG + bb-only.")
    p.add_argument("-P", "--pdb-dir", default=None,
                   help="RCSB-style parent database to resolve parent structures against, "
                        "overriding the build-time directory recorded in the buckets. "
                        f"Only used by {', '.join(PARENT_READING_FLAGS)}.")

    sel = p.add_argument_group("selection (both modes)")
    sel.add_argument("--members", action="store_true",
        help="Write every member of each selected cluster instead of just its "
             "centroid. Needs exactly one --aa-buckets and one --subset-sizes value, "
             "and a reachable parent PDB database: members store no coordinates.")
    sel.add_argument("--aa-buckets", nargs="+", default=None,
        help="AA bucket label(s) to materialize. Default: all buckets.")
    sel.add_argument("--subset-sizes", type=int, nargs="+", default=None,
        help="Subset size(s) to walk. Default: all present under the nr_vdgs root.")
    sel.add_argument("-n", "--top-clusters", type=int, default=None,
        help="Per bucket, take only the N largest clusters. Default: all.")
    sel.add_argument("-k", "--clusters", type=int, nargs="+", default=None,
        help="Take these cluster ID(s) instead of ranking by size. Ids are unique "
             "within a charge-sign bucket and apply to every matching sign. Needs "
             "exactly one --aa-buckets and one --subset-sizes value. Overrides "
             "--min-cluster-size.")

    g1 = p.add_argument_group("centroids mode only (the default mode)")
    g1.add_argument("--min-cluster-size", type=int, default=None,
        help="Skip clusters with fewer than this many members (inclusive). "
             "Default: 10; pass 1 for every cluster.")
    g1.add_argument("--max-files", type=int, default=None, help="Cap on files.")

    sel.add_argument("-w", "--align-cg-weight", type=float, default=None,
        help="Weight on the CG block, the rest going to the vdM backbone, when "
             "fitting onto a reference. Default: 0.99.")

    g2 = p.add_argument_group("members mode only (--members)")
    g2.add_argument("-m", "--reps", type=int, default=None,
        help="Max structures to write per cluster, counting the centroid: the "
             "centroid plus REPS-1 members sampled at random. 0 = all. Default: 20.")
    g2.add_argument("--seed", type=int, default=None,
        help="Seed for the --reps sample, so a rerun reproduces it. Default: 0.")
    return p.parse_args()

def frag_name_from_nr_root(nr_root):
    return sanitize(os.path.basename(os.path.dirname(nr_root)))

def warn_if_job_incomplete(nr_root):
    frag_dir = os.path.dirname(os.path.normpath(nr_root))
    frag_name = os.path.basename(frag_dir)
    vdg_lib_dir = os.path.dirname(frag_dir)
    if not check_vdg_job_status(frag_name, vdg_lib_dir):
        print(f"[WARNING] Fragment {frag_name!r} in {vdg_lib_dir} has no completed "
              f"vdG-generation job (no 'Job completed.' in its log); its nr_vdgs/ may "
              f"be a partial build. Materializing anyway.")

def _get_weights(num_cg_atoms, num_bb_atoms, align_cg_weight):
    return np.concatenate([
        np.full(num_cg_atoms, align_cg_weight / max(num_cg_atoms, 1)),
        np.full(num_bb_atoms, max(1.0 - align_cg_weight, 0.0) / max(num_bb_atoms, 1))])

def _kabsch_weighted_rowvec(mobile, target, weights=None):
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    if mobile.shape != target.shape:
        raise ValueError(f"Kabsch shape mismatch: {mobile.shape} vs {target.shape}")
    if mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError(f"Kabsch expects an (N, 3) array, got {mobile.shape}")
    if weights is None:
        weights = np.ones(mobile.shape[0], dtype=float)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    if weights.shape[0] != mobile.shape[0]:
        raise ValueError(f"Weight length mismatch: {weights.shape[0]} vs {mobile.shape[0]}")
    if not np.any(weights > 0):
        raise ValueError("At least one alignment weight must be positive.")
    weights = weights / weights.sum()
    mob_cent = np.sum(mobile * weights[:, None], axis=0)
    tar_cent = np.sum(target * weights[:, None], axis=0)
    R = utils._proper_kabsch_rotations(
        (mobile - mob_cent).T @ ((target - tar_cent) * weights[:, None]))
    return R, tar_cent - mob_cent @ R

def _best_perm_rigid(cg_coords, bb_flat, target_flat, weights=None,
    cg_automorphisms=None):
    cg_coords = np.asarray(cg_coords, dtype=float)
    target_flat = np.asarray(target_flat, dtype=float)
    if weights is None:
        weights = np.ones(target_flat.shape[0], dtype=float)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    weights = weights / weights.sum()

    best = None
    for perm in (cg_automorphisms or [tuple(range(cg_coords.shape[0]))]):
        mobile = np.vstack([cg_coords[np.asarray(perm, dtype=np.intp)]] +
                           ([] if bb_flat is None else
                            [np.asarray(bb_flat, dtype=float).reshape(-1, 3)]))
        R, t = _kabsch_weighted_rowvec(mobile, target_flat, weights)
        ssd = float(np.sum(weights[:, None] * (mobile @ R + t - target_flat) ** 2))
        if best is None or ssd < best[0]:
            best = (ssd, R, t)
    return best[1], best[2]

def _load_cg_automorphisms(nr_root):
    return vdg_npz.load_cg_symmetry(os.path.dirname(nr_root), frag_name=None)[1]

def align_first_to_reference(ag, cg_coords, cg_vdmbb_coords):
    cg_coords = np.asarray(cg_coords, dtype=float)
    if cg_coords.shape[0] < 3:
        raise ValueError("Need at least 3 CG atoms for first-centroid alignment.")
    R, t = _kabsch_weighted_rowvec(
        cg_coords[:3], [[0.0, 0.0, 0.0], [-1.0, 0.0, 1.0], [1.0, -1.0, 0.0]])
    ag.setCoords(ag.getCoords() @ R + t)
    return ag, np.asarray(cg_vdmbb_coords, dtype=float) @ R + t

def _apply_rigid(ag, cg_vdmbb_flat, R, t):
    ag.setCoords(ag.getCoords() @ R + t)
    return ag, np.asarray(cg_vdmbb_flat, dtype=float) @ R + t

def align_cg_to_reference(ag, cg_coords, cg_vdmbb_flat, ref_cg_coords,
                          cg_automorphisms=None):
    R, t = _best_perm_rigid(cg_coords, None, ref_cg_coords,
                            cg_automorphisms=cg_automorphisms)
    return _apply_rigid(ag, cg_vdmbb_flat, R, t)

def align_to_bucket_reference(ag, cg_coords, vdm_bb_coords, cg_vdmbb_flat,
                              ref_flat, weights, cg_automorphisms=None):
    R, t = _best_perm_rigid(cg_coords, vdm_bb_coords, ref_flat, weights,
                            cg_automorphisms)
    return _apply_rigid(ag, cg_vdmbb_flat, R, t)

def rank_clusters(cluster_size, min_cluster_size=1, top_n=None):
    cluster_size = np.asarray(cluster_size)
    order = np.argsort(-cluster_size, kind="stable")
    if min_cluster_size > 1:
        order = order[cluster_size[order] >= min_cluster_size]
    if top_n is not None:
        order = order[:top_n]
    return order

def align_member_to_centroid(ag, cg_coords, vdm_bb_coords, target_cg_vdmbb_flat,
    weights, cg_automorphisms):
    R, t = _best_perm_rigid(cg_coords, vdm_bb_coords, target_cg_vdmbb_flat,
                            weights, cg_automorphisms)
    ag.setCoords(ag.getCoords() @ R + t)
    return ag

def make_base_tag(kw):
    parts = [sanitize(strip_known_exts(kw["parent_pdb_path"]))]
    parts += [format_residue_tag(resname, seg, chain, resnum)
              for resname, seg, chain, resnum in zip(
                  kw["scrr_resname"], kw["scrr_seg"], kw["scrr_chain"], kw["scrr_resnum"])]
    parts.append(format_residue_tag(kw["cg_resname"], kw["cg_seg"], kw["cg_chain"],
                                    kw["cg_resnum"]))
    return "_".join(parts)

def _bucket_paths(nr_root, subset_name, aa_buckets=None):
    for charge_sign in vdg_npz.CHARGE_SIGNS:
        charge_dir = os.path.join(nr_root, subset_name, charge_sign)
        if not os.path.isdir(charge_dir):
            continue
        for fname in sorted(os.listdir(charge_dir)):
            if not fname.endswith(".npz"):
                continue
            aa_label = os.path.splitext(fname)[0]
            if aa_buckets is None or aa_label in aa_buckets:
                yield charge_sign, aa_label, os.path.join(charge_dir, fname)

def materialize_nr_vdgs(nr_root, out_root, max_files=None,
    aa_buckets=None, subset_sizes=None, top_clusters=None, cluster_ids=None,
    min_cluster_size=DEFAULT_MIN_CLUSTER_SIZE, include_full_ligand=False,
    include_backbone_carbonyl=False, include_sidechain=True, pdb_dir=None,
    align_cg_weight=DEFAULT_ALIGN_CG_WEIGHT):
    nr_root = os.path.abspath(nr_root)
    out_root = os.path.abspath(out_root)
    warn_if_job_incomplete(nr_root)
    cg_automorphisms = _load_cg_automorphisms(nr_root)
    extras = [flag for flag, on in (("--include-full-ligand", include_full_ligand),
                                    ("--carbonyl", include_backbone_carbonyl),
                                    ("--sidechain", include_sidechain)) if on]
    if extras and not vdg_npz.parent_extras_available(extras, pdb_dir=pdb_dir):
        include_full_ligand = include_backbone_carbonyl = include_sidechain = False
        extras = []
    made_dirs, names = set(), NameRegistry()
    fresh_dir(out_root, made_dirs)
    aa_buckets = set(aa_buckets) if aa_buckets is not None else None
    subset_sizes = {str(s) for s in subset_sizes} if subset_sizes is not None else None
    total_written = 0
    global_ref_cg = None

    for subset_name in sorted(os.listdir(nr_root)):
        if not os.path.isdir(os.path.join(nr_root, subset_name)):
            continue
        if subset_sizes is not None and subset_name not in subset_sizes:
            continue
        for charge_sign, aa_label, npz_path in _bucket_paths(
            nr_root, subset_name, aa_buckets):
            data = vdg_npz.load_bucket_npz(npz_path)
            if data is None:
                print(f"[WARNING] Skipping unreadable bucket {npz_path}.")
                continue
            if extras:
                if not vdg_npz.parent_extras_available(extras, data, pdb_dir=pdb_dir):
                    include_full_ligand = False
                    include_backbone_carbonyl = include_sidechain = False
                extras = []
            clus_id, clus_size = data["cluster_id"], data["cluster_size"]
            n_cg = data["nr_cg_coords"].shape[1]
            bucket_ref = None

            if cluster_ids is None:
                order = rank_clusters(clus_size, min_cluster_size=min_cluster_size,
                                      top_n=top_clusters)
            else:
                wanted = set(cluster_ids)
                order = [i for i in rank_clusters(clus_size)
                         if int(clus_id[i]) in wanted]
                missing = wanted - {int(clus_id[i]) for i in order}
                if missing:
                    print(f"[WARNING] Cluster id(s) {sorted(missing)} not found in "
                          f"{npz_path}; skipping.")
            for idx in order:
                idx = int(idx)
                kw = vdg_npz.nr_build_kwargs(
                    data, idx, include_full_ligand, pdb_dir=pdb_dir,
                    include_backbone_carbonyl=include_backbone_carbonyl,
                    include_sidechain=include_sidechain)
                ag, cg_vdmbb_coords = vdg_npz.build_vdg_atomgroup_from_npz(**kw)
                if global_ref_cg is None:
                    ag, cg_vdmbb_coords = align_first_to_reference(
                        ag, kw["cg_coords"], cg_vdmbb_coords)
                    global_ref_cg = ag.getCoords()[:n_cg].copy()
                elif bucket_ref is None:
                    ag, cg_vdmbb_coords = align_cg_to_reference(
                        ag, kw["cg_coords"], cg_vdmbb_coords, global_ref_cg,
                        cg_automorphisms)
                else:
                    ag, cg_vdmbb_coords = align_to_bucket_reference(
                        ag, kw["cg_coords"], kw["vdm_bb_coords"], cg_vdmbb_coords,
                        bucket_ref[0], bucket_ref[1], cg_automorphisms)
                if bucket_ref is None:
                    bucket_ref = (
                        cg_vdmbb_coords.copy(),
                        _get_weights(n_cg, kw["vdm_bb_coords"].shape[0] * 3,
                                     align_cg_weight))
                write_pdb_gz(ag, names.claim(
                    f"{frag_name_from_nr_root(nr_root)}_{charge_sign}_"
                    f"{sanitize(aa_label)}_clus{int(clus_id[idx])}_"
                    f"size{int(clus_size[idx])}_{make_base_tag(kw)}", out_root))
                total_written += 1
                if max_files is not None and total_written >= max_files:
                    return total_written

    return total_written

def _member_build_kwargs(mem, cg_coords, vdm_bb_coords, include_full_ligand=False,
                         include_backbone_carbonyl=False, include_sidechain=True):
    kw = {k: mem[k] for k in ("cg_names", "cg_elements", "cg_seg", "cg_chain",
                              "cg_resnum", "cg_resname", "scrr_seg", "scrr_chain",
                              "scrr_resnum", "scrr_resname")}
    kw.update(cg_coords=cg_coords, vdm_bb_coords=vdm_bb_coords,
              parent_pdb_path=mem["pdbpath"],
              include_full_ligand=include_full_ligand,
              include_backbone_carbonyl=include_backbone_carbonyl,
              include_sidechain=include_sidechain)
    return kw

def materialize_cluster_members(nr_root, out_root, subset_size, aa_bucket, cluster_ids=None,
    top_clusters=None, reps=None, seed=0, align_cg_weight=DEFAULT_ALIGN_CG_WEIGHT,
    include_full_ligand=False, include_backbone_carbonyl=False,
    include_sidechain=True, pdb_dir=None):
    nr_root = os.path.abspath(nr_root)
    out_root = os.path.abspath(out_root)
    warn_if_job_incomplete(nr_root)
    bucket_paths = list(_bucket_paths(nr_root, str(subset_size), {aa_bucket}))
    if not bucket_paths:
        raise FileNotFoundError(
            f"No {aa_bucket}.npz under {nr_root}/{subset_size}/<charge_sign>/")

    if (cluster_ids is None) == (top_clusters is None):
        raise ValueError("Pass exactly one of cluster_ids or top_clusters.")
    state = dict(
        frag_name=sanitize(os.path.basename(os.path.dirname(nr_root))),
        cg_automorphisms=_load_cg_automorphisms(nr_root), made_dirs=set(),
        rng=np.random.default_rng(seed), global_ref_cg=None, names=NameRegistry())
    total_written, total_skipped = 0, 0
    for charge_sign, _, npz_path in bucket_paths:
        data = vdg_npz.load_bucket_npz(npz_path)
        if data is None:
            raise ValueError(f"Bucket npz is unreadable or corrupt: {npz_path}")
        vdg_npz.require_parent_pdb_dir(data, pdb_dir=pdb_dir)
        selected_ids = cluster_ids
        if top_clusters is not None:
            selected_ids = [int(data["cluster_id"][i]) for i in
                            rank_clusters(data["cluster_size"], top_n=top_clusters)]
        out_aa_dir = os.path.join(out_root, sanitize(str(subset_size)), charge_sign,
                                  sanitize(aa_bucket))
        fresh_dir(out_aa_dir, state["made_dirs"])
        for cid in selected_ids:
            written, skipped, state["global_ref_cg"] = _materialize_one_cluster(
                data, npz_path, cid, out_aa_dir, state["frag_name"], state["names"],
                state["rng"], reps, state["global_ref_cg"], state["cg_automorphisms"],
                align_cg_weight,
                include_full_ligand, include_backbone_carbonyl, include_sidechain,
                pdb_dir)
            total_written += written
            total_skipped += skipped
    return total_written, total_skipped

def _materialize_one_cluster(data, npz_path, cid, out_aa_dir, frag_name, names,
    rng, reps, global_ref_cg, cg_automorphisms, align_cg_weight,
    include_full_ligand, include_backbone_carbonyl, include_sidechain, pdb_dir):
    matches = np.nonzero(data["cluster_id"] == cid)[0]
    if matches.size == 0:
        print(f"[WARNING] Cluster {cid} not found in {npz_path}; skipping.")
        return 0, 0, global_ref_cg
    idx = int(matches[0])
    cluster_size = int(data["cluster_size"][idx])
    kw = vdg_npz.nr_build_kwargs(
        data, idx, include_full_ligand, pdb_dir=pdb_dir,
        include_backbone_carbonyl=include_backbone_carbonyl,
        include_sidechain=include_sidechain)
    ag_cent, cg_vdmbb_cent = vdg_npz.build_vdg_atomgroup_from_npz(**kw)
    ag_cent, cg_vdmbb_cent = align_first_to_reference(
        ag_cent, kw["cg_coords"], cg_vdmbb_cent)
    n_cg = kw["cg_coords"].shape[0]
    if global_ref_cg is None:
        global_ref_cg = ag_cent.getCoords()[:n_cg].copy()
    else:
        R, t = _best_perm_rigid(ag_cent.getCoords()[:n_cg], None, global_ref_cg,
                                cg_automorphisms=cg_automorphisms)
        ag_cent, cg_vdmbb_cent = _apply_rigid(ag_cent, cg_vdmbb_cent, R, t)
    member_state = dict(
        weights=_get_weights(n_cg, kw["vdm_bb_coords"].shape[0] * 3,
                             align_cg_weight), skipped=0)
    out_dir_clus = os.path.join(out_aa_dir, f"clus{cid}_size{cluster_size}")
    os.makedirs(out_dir_clus, exist_ok=True)
    write_pdb_gz(ag_cent, names.claim(
        f"{frag_name}_NR_clus{cid}_size{cluster_size}_{make_base_tag(kw)}",
        out_dir_clus))
    mem_rows = vdg_npz.cluster_member_indices(data, cid)
    n_members_total = len(mem_rows)
    if reps and n_members_total > reps - 1:
        mem_rows = mem_rows[np.sort(
            rng.choice(n_members_total, size=reps - 1, replace=False))]
    members = vdg_npz.load_cluster_members(data, cid, pdb_dir=pdb_dir,
                                           indices=mem_rows)
    for mem in members:
        mem_struct = vdg_npz.parse_pdb_or_none(mem["pdbpath"], "this member is skipped")
        cg_coords = vdm_bb_coords = None
        if mem_struct is not None:
            cg_coords, vdm_bb_coords = vdg_npz.rederive_member_coords(
                mem["pdbpath"], mem["cg_seg"], mem["cg_chain"], mem["cg_resnum"],
                mem["cg_names"], mem["scrr_seg"], mem["scrr_chain"],
                mem["scrr_resnum"], parsed_pdb=mem_struct)
        if cg_coords is None:
            member_state["skipped"] += 1
            continue
        mem_kw = _member_build_kwargs(mem, cg_coords, vdm_bb_coords,
            include_full_ligand, include_backbone_carbonyl, include_sidechain)
        ag, _ = vdg_npz.build_vdg_atomgroup_from_npz(**mem_kw, parent_struct=mem_struct)
        ag = align_member_to_centroid(ag, cg_coords, vdm_bb_coords, cg_vdmbb_cent,
            member_state["weights"], cg_automorphisms)
        write_pdb_gz(ag, names.claim(f"{frag_name}_{make_base_tag(mem_kw)}",
                                     out_dir_clus))
    n_expected = 1 + n_members_total
    if n_expected != cluster_size:
        print(f"[WARNING] Cluster {cid}: cluster_size={cluster_size} but the npz "
              f"holds {n_expected} record(s) for it (1 centroid + "
              f"{n_members_total} non-centroid).")
    if member_state["skipped"]:
        print(f"Cluster {cid}: {member_state['skipped']} member(s) could not be re-derived "
              "from their parent PDB and were skipped.")
    return (1 + len(members) - member_state["skipped"], member_state["skipped"],
            global_ref_cg)

def _one(values, flag):
    if values is None or len(values) != 1:
        raise ValueError(f"{flag} must name exactly one value here (got "
                         f"{'nothing' if values is None else f'{len(values)} values'}).")
    return values[0]

def _check_pdb_dir(args):
    if args.pdb_dir is None:
        return
    if not (args.include_full_ligand or args.carbonyl or args.sidechain
            or args.members):
        raise ValueError(
            f"--pdb-dir has no effect without one of {', '.join(PARENT_READING_FLAGS)}; "
            "every other mode is served from the bucket npz alone.")
    if not os.path.isdir(args.pdb_dir):
        raise ValueError(f"--pdb-dir is not a directory: {args.pdb_dir}")

def _reject_other_mode_flags(args):
    bad, mode, other = ((("--reps", args.reps), ("--seed", args.seed)),
                        "centroids mode", "--members") if not args.members else (
                        (("--min-cluster-size", args.min_cluster_size),
                         ("--max-files", args.max_files)),
                        "members mode (--members)", "the default mode")
    passed = [flag for flag, value in bad if value is not None]
    if passed:
        raise ValueError(f"{', '.join(passed)} has no effect in {mode}; "
                         f"it only applies to {other}.")

def main():
    args = parse_args()
    if args.clusters is not None and args.top_clusters is not None:
        raise ValueError("Pass either --clusters or --top-clusters, not both.")
    _reject_other_mode_flags(args)
    _check_pdb_dir(args)
    if args.reps is not None and args.reps < 0:
        raise ValueError("--reps must be 0 (all) or a positive count.")
    if args.max_files is not None and args.max_files < 1:
        raise ValueError("--max-files must be a positive count (omit it for no cap).")
    if args.align_cg_weight is not None and not 0.0 <= args.align_cg_weight <= 1.0:
        raise ValueError("--align-cg-weight must be between 0 and 1.")

    if args.members:
        if args.clusters is None and args.top_clusters is None:
            raise ValueError("--members needs --clusters or --top-clusters to say "
                             "which cluster(s) to write the members of.")
        n, skipped = materialize_cluster_members(
            nr_root=args.cg_nr_vdgs_root,
            out_root=args.out_dir,
            subset_size=_one(args.subset_sizes, "--subset-sizes"),
            aa_bucket=_one(args.aa_buckets, "--aa-buckets"),
            cluster_ids=args.clusters,
            top_clusters=args.top_clusters,
            reps=DEFAULT_REPS if args.reps is None else args.reps,
            seed=DEFAULT_SEED if args.seed is None else args.seed,
            align_cg_weight=(DEFAULT_ALIGN_CG_WEIGHT if args.align_cg_weight is None
                             else args.align_cg_weight),
            include_full_ligand=args.include_full_ligand,
            include_backbone_carbonyl=args.carbonyl,
            include_sidechain=args.sidechain,
            pdb_dir=args.pdb_dir)
        print(f"Wrote {n} PDB file(s)." + (
            f" {skipped} member(s) skipped; see the [WARNING] line(s) above for the "
            "atom and parent PDB in each case." if skipped else ""))
        return

    if args.clusters is not None:
        _one(args.subset_sizes, "--subset-sizes")
        _one(args.aa_buckets, "--aa-buckets")
        if args.min_cluster_size is not None:
            raise ValueError("--min-cluster-size has no effect alongside --clusters, "
                             "which names the clusters to write outright.")

    n = materialize_nr_vdgs(
        nr_root=args.cg_nr_vdgs_root,
        out_root=args.out_dir,
        max_files=args.max_files,
        aa_buckets=args.aa_buckets,
        subset_sizes=args.subset_sizes,
        top_clusters=args.top_clusters,
        cluster_ids=args.clusters,
        min_cluster_size=(DEFAULT_MIN_CLUSTER_SIZE if args.min_cluster_size is None
                          else args.min_cluster_size),
        include_full_ligand=args.include_full_ligand,
        include_backbone_carbonyl=args.carbonyl,
        include_sidechain=args.sidechain,
        pdb_dir=args.pdb_dir,
        align_cg_weight=(DEFAULT_ALIGN_CG_WEIGHT if args.align_cg_weight is None
                         else args.align_cg_weight))
    print(f"Wrote {n} PDB file(s).")

if __name__ == "__main__":
    vdg_npz.run_cli(main)
