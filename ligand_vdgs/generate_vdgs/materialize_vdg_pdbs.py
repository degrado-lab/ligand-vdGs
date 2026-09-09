"""
Materialize vdG PDBs (CG + vdM backbone + sidechain) from nr_vdgs .npz files.

Centroids mode (default): one PDB per cluster, largest first, flat in --out-dir.

    python materialize_vdg_pdbs.py -c <vdglib>/<frag>/nr_vdgs -o <out> \
        --aa-buckets ASP_bb GLU_bb --top-clusters 25 --min-cluster-size 10

    <out>/<frag>_<AA_BUCKET>_clus<id>_size<size>_<source>_<vdm tags>_<lig tag>.pdb.gz

Members mode (--members): every member of each selected cluster, one dir per cluster.
Members store no coordinates and are re-derived from their parent PDB
(vdg_npz_utils.rederive_member_coords), so a cluster of size N yields 1 + (N-1) files.

    python materialize_vdg_pdbs.py -c <vdglib>/<frag>/nr_vdgs -o <out> \
        --members --subset-sizes 2 --aa-buckets ASP_bb --top-clusters 5 --reps 20

    <out>/<subset_size>/<AA_BUCKET>/clus<id>_size<size>/
        <frag>_NR_clus<id>_size<size>_<source>_<vdm tags>_<lig tag>.pdb.gz
        <frag>_<source>_<vdm tags>_<lig tag>.pdb.gz   (one per member)

Fits are on CG + vdM backbone N/CA/C (--align-cg-weight splits the two) against a reference
from the same AA bucket; the hop between buckets is on the CG alone, over the CG
automorphism group, because backbone slots correspond only within a bucket.

Both modes select by rank (--top-clusters) or by id (--clusters). Ids are unique only
within one bucket npz, so --clusters and --members each need exactly one --aa-buckets
and one --subset-sizes value.

Naming (functions/vdg_pdb_io.py, shared with tools/write_vdg_hit_pdbs.py): '_'-joined
fields; a residue tag is seg_chain_resnum_resname, empty segment kept ('_A_370_ASP');
1 + subset_size tags, vdM slots in slot order then the ligand. The head is NOT
positionally parseable (<frag> and <AA_BUCKET> each span a variable number of fields)
-- parse from the right: strip the extension and any '~<n>' suffix, take the last
4 * (subset_size + 1) fields, cut into tags of 4.

Parent structures resolve against the parent PDB database the library was mined from: 
the path recorded at build time, unless $PARENT_PDBS_DIR or --pdb-dir overrides it.
--sidechain is on by default. Regardless of the "extra" display atoms, every fit
stays on CG + vdM N/CA/C, so with no database reachable they are dropped with a
warning and the run still writes backbone-only vdGs from the npz. --members has no
such fallback: a member stores no coordinates, so it fails outright. 
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

# Defaults live here, not in argparse, so main() can tell "not passed" from "passed
# the default" and reject the other mode's flags instead of ignoring them.
DEFAULT_MIN_CLUSTER_SIZE = 10
DEFAULT_REPS = 20
DEFAULT_ALIGN_CG_WEIGHT = 0.99
DEFAULT_SEED = 0

# Modes that read a parent PDB; everything else is served from the npz. Sidechains
# are on by default, so every run reads parents unless --no-sidechain turns them off.
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
                   help="RCSB-style PDB mirror to resolve parent structures against, "
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
             "only within one bucket npz, so needs exactly one --aa-buckets and one "
             "--subset-sizes value. Overrides --min-cluster-size.")

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
    """Filename prefix: the fragment dir enclosing ``nr_vdgs``. A display label, not
    necessarily the CG's SMARTS. Pass an abspath'd root, else a trailing slash yields
    "nr_vdgs"."""
    return sanitize(os.path.basename(os.path.dirname(nr_root)))


def warn_if_job_incomplete(nr_root):
    """Warn (not raise -- a partial library is still worth looking at) if the job never
    completed. Both entry points bypass ``load_vdg_bucket``, where this check normally
    lives; the log is named after the on-disk dir, so use the raw basename."""
    frag_dir = os.path.dirname(os.path.normpath(nr_root))
    frag_name = os.path.basename(frag_dir)
    vdg_lib_dir = os.path.dirname(frag_dir)
    if not check_vdg_job_status(frag_name, vdg_lib_dir):
        print(f"[WARNING] Fragment {frag_name!r} in {vdg_lib_dir} has no completed "
              f"vdG-generation job (no 'Job completed.' in its log); its nr_vdgs/ may "
              f"be a partial build. Materializing anyway.")


def _get_weights(num_cg_atoms, num_bb_atoms, align_cg_weight):
    cg_w = np.full(num_cg_atoms, align_cg_weight / max(num_cg_atoms, 1), dtype=float)
    bb_total = max(1.0 - align_cg_weight, 0.0)
    bb_w = np.full(num_bb_atoms, bb_total / max(num_bb_atoms, 1), dtype=float)
    return np.concatenate([cg_w, bb_w])


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
    H = (mobile - mob_cent).T @ ((target - tar_cent) * weights[:, None])
    R = utils._proper_kabsch_rotations(H)
    return R, tar_cent - mob_cent @ R


def _best_perm_rigid(cg_coords, bb_flat, target_flat, weights=None,
    cg_automorphisms=None):
    """Rigid transform onto ``target_flat``, minimised over the CG automorphism group.

    SMARTS-slot order is not a unique correspondence for a symmetric CG, so the identity
    alone can fit the wrong transform -- which carries the vdM backbone too. Only the CG
    block is permuted; ``bb_flat`` (None for a CG-only fit) keeps its row order.
    """
    cg_coords = np.asarray(cg_coords, dtype=float)
    target_flat = np.asarray(target_flat, dtype=float)
    blocks = [] if bb_flat is None else [np.asarray(bb_flat, dtype=float).reshape(-1, 3)]
    if weights is None:
        weights = np.ones(target_flat.shape[0], dtype=float)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    weights = weights / weights.sum()

    perms = cg_automorphisms if cg_automorphisms else [tuple(range(cg_coords.shape[0]))]
    best = None
    for perm in perms:
        mobile = np.vstack([cg_coords[np.asarray(perm, dtype=np.intp)]] + blocks)
        R, t = _kabsch_weighted_rowvec(mobile, target_flat, weights)
        ssd = float(np.sum(weights[:, None] * (mobile @ R + t - target_flat) ** 2))
        if best is None or ssd < best[0]:
            best = (ssd, R, t)
    return best[1], best[2]


def _load_cg_automorphisms(nr_root):
    """CG automorphism group recorded beside the buckets. Raises if the sidecar is
    absent: the identity group would superpose a symmetric CG under the wrong
    correspondence, plausibly enough to go unnoticed."""
    _, perms = vdg_npz.load_cg_symmetry(os.path.dirname(nr_root), frag_name=None)
    return perms


def align_first_to_reference(ag, cg_coords, cg_vdmbb_coords):
    ref = np.array([[0.0, 0.0, 0.0], [-1.0, 0.0, 1.0], [1.0, -1.0, 0.0]], dtype=float)
    cg_coords = np.asarray(cg_coords, dtype=float)
    if cg_coords.shape[0] < 3:
        raise ValueError("Need at least 3 CG atoms for first-centroid alignment.")
    R, t = _kabsch_weighted_rowvec(cg_coords[:3], ref)
    ag.setCoords(ag.getCoords() @ R + t)
    return ag, np.asarray(cg_vdmbb_coords, dtype=float) @ R + t


def _apply_rigid(ag, cg_vdmbb_flat, R, t):
    """Move an AtomGroup and its CG+backbone fit block by one rigid transform.

    The block is carried separately rather than sliced back out of ``ag``: with
    --include-full-ligand or --carbonyl the AtomGroup holds display atoms that are
    not part of the fit, so its row order is not the block's.
    """
    ag.setCoords(ag.getCoords() @ R + t)
    return ag, np.asarray(cg_vdmbb_flat, dtype=float) @ R + t


def align_cg_to_reference(ag, cg_coords, cg_vdmbb_flat, ref_cg_coords,
                          cg_automorphisms=None):
    """Superpose on the CG block alone, over the automorphism group.

    Used where the target's vdM backbone slots are not in correspondence with this
    record's -- across AA buckets, whose slots differ in identity and, for a different
    subset size, in count. Within a bucket, ``align_to_bucket_reference`` fits the
    backbone too.
    """
    R, t = _best_perm_rigid(cg_coords, None, ref_cg_coords,
                            cg_automorphisms=cg_automorphisms)
    return _apply_rigid(ag, cg_vdmbb_flat, R, t)


def align_to_bucket_reference(ag, cg_coords, vdm_bb_coords, cg_vdmbb_flat,
                              ref_flat, weights, cg_automorphisms=None):
    """Superpose on CG + vdM backbone, over the automorphism group.

    Valid only against a target from the same AA bucket: one bucket's records share
    a slot count and slot AA labels, so row i of the backbone block means the same
    thing in both. Only the CG block is permuted -- backbone slot order is fixed
    within a bucket (AA-duplicate permutations were expanded before clustering).
    """
    R, t = _best_perm_rigid(cg_coords, vdm_bb_coords, ref_flat, weights,
                            cg_automorphisms)
    return _apply_rigid(ag, cg_vdmbb_flat, R, t)


def rank_clusters(cluster_size, min_cluster_size=1, top_n=None):
    """Indices into a bucket's centroid arrays, largest cluster first. Ties keep npz
    order (stable sort), so reruns agree."""
    cluster_size = np.asarray(cluster_size)
    order = np.argsort(-cluster_size, kind="stable")
    if min_cluster_size > 1:
        order = order[cluster_size[order] >= min_cluster_size]
    if top_n is not None:
        order = order[:top_n]
    return order


def align_member_to_centroid(ag, cg_coords, vdm_bb_coords, target_cg_vdmbb_flat,
    weights, cg_automorphisms):
    """Align a re-derived member onto its cluster's (already-aligned) centroid. Only the
    CG block is permuted: backbone slot order is consistent within a cluster
    (AA-duplicate permutations were expanded into distinct records before clustering)."""
    R, t = _best_perm_rigid(cg_coords, vdm_bb_coords, target_cg_vdmbb_flat,
                            weights, cg_automorphisms)
    ag.setCoords(ag.getCoords() @ R + t)
    return ag


def make_base_tag(kw):
    """``<source pdb>_<vdM tags>_<ligand tag>`` from a build-kwargs dict
    (vdg_npz.nr_build_kwargs / _member_build_kwargs).

    One tag per vdM slot, undeduplicated so the count is always the subset size, then
    the CG's residue (one per record -- get_cg_atoms refuses a CG spanning more). The
    ligand tag separates vdGs from one deposition sharing a vdM residue but not a
    ligand copy or CG site.
    """
    parts = [sanitize(strip_known_exts(kw["parent_pdb_path"]))]
    parts += [format_residue_tag(resname, seg, chain, resnum)
              for resname, seg, chain, resnum in zip(
                  kw["scrr_resname"], kw["scrr_seg"], kw["scrr_chain"], kw["scrr_resnum"])]
    parts.append(format_residue_tag(kw["cg_resname"], kw["cg_seg"], kw["cg_chain"],
                                    kw["cg_resnum"]))
    return "_".join(parts)


def materialize_nr_vdgs(nr_root, out_root, max_files=None,
    aa_buckets=None, subset_sizes=None, top_clusters=None, cluster_ids=None,
    min_cluster_size=DEFAULT_MIN_CLUSTER_SIZE, include_full_ligand=False,
    include_backbone_carbonyl=False, include_sidechain=True, pdb_dir=None,
    align_cg_weight=DEFAULT_ALIGN_CG_WEIGHT):
    """One PDB per cluster (the centroid), largest first, flat in ``out_root`` and all in
    one frame. No subset-size dir or field: a bucket label has one AA token per vdM slot,
    so it already says its subset size. ``cluster_ids`` names clusters outright; ids only
    mean something within one npz, so pass a single bucket and subset size.

    Centroids are fitted on CG + vdM backbone (``align_cg_weight`` splits the two, as in
    members mode) against the first centroid of their own bucket. Each bucket's first
    centroid is placed on the run-wide CG frame by its CG alone, since backbone slots
    are only in correspondence within a bucket -- so a single-bucket run is CG+backbone
    throughout, and a multi-bucket run is CG+backbone within each bucket.
    """
    nr_root = os.path.abspath(nr_root)
    out_root = os.path.abspath(out_root)
    frag_name = frag_name_from_nr_root(nr_root)
    warn_if_job_incomplete(nr_root)
    # Sidecar first: a missing cg_symmetry.npz is fatal, and raising after fresh_dir
    # has claimed out_root leaves a directory the next attempt refuses to reuse.
    cg_automorphisms = _load_cg_automorphisms(nr_root)
    extras = [flag for flag, on in (("--include-full-ligand", include_full_ligand),
                                    ("--carbonyl", include_backbone_carbonyl),
                                    ("--sidechain", include_sidechain)) if on]
    if extras and not vdg_npz.parent_extras_available(extras, pdb_dir=pdb_dir):
        include_full_ligand = include_backbone_carbonyl = include_sidechain = False
        extras = []
    # The freshness unit is whatever shares a frame -- here the whole run -- so name
    # every bucket in one invocation; fresh_dir refuses a non-empty out_root. Partial
    # sets are still possible (--max-files, a skipped bucket): check the return value.
    made_dirs, names = set(), NameRegistry()
    fresh_dir(out_root, made_dirs)
    aa_buckets = set(aa_buckets) if aa_buckets is not None else None
    subset_sizes = {str(s) for s in subset_sizes} if subset_sizes is not None else None
    total_written = 0
    global_ref_cg = None  # cross-bucket CG frame: every bucket's first centroid lands here

    for subset_name in sorted(os.listdir(nr_root)):
        subset_dir = os.path.join(nr_root, subset_name)
        if not os.path.isdir(subset_dir):
            continue
        if subset_sizes is not None and subset_name not in subset_sizes:
            continue
        for fname in sorted(os.listdir(subset_dir)):
            if not fname.endswith(".npz"):
                continue
            aa_label = os.path.splitext(fname)[0]
            if aa_buckets is not None and aa_label not in aa_buckets:
                continue
            npz_path = os.path.join(subset_dir, fname)
            # Read eagerly so a corrupt bucket is caught by name, not mid-write. Skip
            # rather than abort: raising would force clearing out_root by hand.
            data = vdg_npz.load_bucket_npz(npz_path)
            if data is None:
                print(f"[WARNING] Skipping unreadable bucket {npz_path}.")
                continue
            if extras:  # bucket-level field, so one check covers the run
                if not vdg_npz.parent_extras_available(extras, data, pdb_dir=pdb_dir):
                    include_full_ligand = False
                    include_backbone_carbonyl = include_sidechain = False
                extras = []
            cg = data["nr_cg_coords"]
            clus_id, clus_size = data["cluster_id"], data["cluster_size"]
            n_cg = cg.shape[1]
            # Reset per bucket: the backbone target is only meaningful against records
            # that share this bucket's slot count and slot AA labels.
            bucket_ref_flat, weights = None, None

            if cluster_ids is None:
                order = rank_clusters(clus_size, min_cluster_size=min_cluster_size,
                                      top_n=top_clusters)
            else:
                # A named id outranks the size floor; ranking only fixes write order.
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
                elif bucket_ref_flat is None:
                    ag, cg_vdmbb_coords = align_cg_to_reference(
                        ag, kw["cg_coords"], cg_vdmbb_coords, global_ref_cg,
                        cg_automorphisms)
                else:
                    ag, cg_vdmbb_coords = align_to_bucket_reference(
                        ag, kw["cg_coords"], kw["vdm_bb_coords"], cg_vdmbb_coords,
                        bucket_ref_flat, weights, cg_automorphisms)
                if bucket_ref_flat is None:
                    bucket_ref_flat = cg_vdmbb_coords.copy()
                    weights = _get_weights(n_cg, kw["vdm_bb_coords"].shape[0] * 3,
                                           align_cg_weight)
                stem = (f"{frag_name}_{sanitize(aa_label)}_clus{int(clus_id[idx])}"
                        f"_size{int(clus_size[idx])}_{make_base_tag(kw)}")
                # Unique by construction, but nothing checks the path on write.
                write_pdb_gz(ag, names.claim(stem, out_root))
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
    """Members of one or more clusters, aligned to each cluster's centroid, one
    directory per cluster.

    Pass either ``cluster_ids`` or ``top_clusters=N``. ``reps`` caps structures per
    cluster **including the centroid** (None/0 = all). Members are sampled uniformly
    (``seed`` fixes the draw): npz order tracks the parent PDB list, so taking the front
    would return near-identical structures from one deposition series.
    """
    nr_root = os.path.abspath(nr_root)
    out_root = os.path.abspath(out_root)
    warn_if_job_incomplete(nr_root)
    npz_path = os.path.join(nr_root, str(subset_size), f"{aa_bucket}.npz")
    if not os.path.isfile(npz_path):
        raise FileNotFoundError(f"No such bucket npz: {npz_path}")
    # Named outright, so an unreadable bucket fails the request; unlike the centroids
    # walk, there is nothing to fall back to.
    data = vdg_npz.load_bucket_npz(npz_path)
    if data is None:
        raise ValueError(f"Bucket npz is unreadable or corrupt: {npz_path}")
    # Members carry no coordinates, so every one of them is re-derived from its parent
    # biounit -- unconditionally, unlike centroids mode. 
    vdg_npz.require_parent_pdb_dir(data, pdb_dir=pdb_dir)

    if (cluster_ids is None) == (top_clusters is None):
        raise ValueError("Pass exactly one of cluster_ids or top_clusters.")
    if top_clusters is not None:
        order = rank_clusters(data["cluster_size"], top_n=top_clusters)
        cluster_ids = [int(data["cluster_id"][i]) for i in order]

    frag_name = sanitize(os.path.basename(os.path.dirname(nr_root)))
    cg_automorphisms = _load_cg_automorphisms(nr_root)

    # Centroids share one frame across clusters, so the bucket dir -- not each
    # clus<id> subdir -- is what a rerun must not mix with.
    out_aa_dir = os.path.join(out_root, sanitize(str(subset_size)), sanitize(aa_bucket))
    made_dirs = set()
    fresh_dir(out_aa_dir, made_dirs)
    rng = np.random.default_rng(seed)

    total_written, total_skipped = 0, 0
    global_ref_cg = None  # CG block of the first centroid; frame shared by all clusters
    # Run-wide: a tag names a parent PDB and residue set, not the CG site, so it can
    # recur across clusters. Separate dirs prevent overwrites; PyMOL names by basename.
    names = NameRegistry()
    for cid in cluster_ids:
        matches = np.nonzero(data["cluster_id"] == cid)[0]
        if matches.size == 0:
            print(f"[WARNING] Cluster {cid} not found in {npz_path}; skipping.")
            continue
        idx = int(matches[0])
        cluster_size = int(data["cluster_size"][idx])

        kw = vdg_npz.nr_build_kwargs(
            data, idx, include_full_ligand, pdb_dir=pdb_dir,
            include_backbone_carbonyl=include_backbone_carbonyl,
            include_sidechain=include_sidechain)
        ag_cent, cg_vdmbb_cent = vdg_npz.build_vdg_atomgroup_from_npz(**kw)
        ag_cent, cg_vdmbb_cent = align_first_to_reference(ag_cent, kw["cg_coords"], cg_vdmbb_cent)
        n_cg = kw["cg_coords"].shape[0]

        # One frame for every cluster's centroid, so centroids copied out of different
        # cluster dirs are comparable; align_first_to_reference pins only 3 CG atoms
        # (~2 A of avoidable spread on a flexible 5-atom CG, measured on NCCCS/COCCn).
        # Members follow their own centroid, so within-cluster geometry is untouched.
        if global_ref_cg is None:
            global_ref_cg = ag_cent.getCoords()[:n_cg].copy()
        else:
            R, t = _best_perm_rigid(ag_cent.getCoords()[:n_cg], None, global_ref_cg,
                                    cg_automorphisms=cg_automorphisms)
            ag_cent, cg_vdmbb_cent = _apply_rigid(ag_cent, cg_vdmbb_cent, R, t)

        weights = _get_weights(n_cg, kw["vdm_bb_coords"].shape[0] * 3, align_cg_weight)

        out_dir_clus = os.path.join(out_aa_dir, f"clus{cid}_size{cluster_size}")
        os.makedirs(out_dir_clus, exist_ok=True)  # out_aa_dir already claimed

        # Members name only their source; the dir gives id and size. The centroid
        # repeats them so it stays self-identifying once copied out.
        write_pdb_gz(ag_cent, names.claim(
            f"{frag_name}_NR_clus{cid}_size{cluster_size}_{make_base_tag(kw)}",
            out_dir_clus))
        total_written += 1

        # Sample rows before building per-member dicts: with a small --reps and a big
        # cluster, nearly all of that work is wasted.
        mem_rows = vdg_npz.cluster_member_indices(data, cid)
        n_members_total = len(mem_rows)
        # reps counts the always-written centroid, so it buys reps-1 members.
        if reps and n_members_total > reps - 1:
            picks = rng.choice(n_members_total, size=reps - 1, replace=False)
            mem_rows = mem_rows[np.sort(picks)]
        members = vdg_npz.load_cluster_members(data, cid, pdb_dir=pdb_dir,
                                               indices=mem_rows)
        clus_skipped = 0
        for mem in members:
            # Coords and the optional extras both read the parent biounit: parse once,
            # and bail here before rederive_member_coords re-parses it.
            mem_struct = vdg_npz.parse_pdb_or_none(
                mem["pdbpath"], "this member is skipped")
            cg_coords = vdm_bb_coords = None
            if mem_struct is not None:
                cg_coords, vdm_bb_coords = vdg_npz.rederive_member_coords(
                    mem["pdbpath"], mem["cg_seg"], mem["cg_chain"], mem["cg_resnum"],
                    mem["cg_names"], mem["scrr_seg"], mem["scrr_chain"],
                    mem["scrr_resnum"], parsed_pdb=mem_struct)
            if cg_coords is None:
                total_skipped += 1
                clus_skipped += 1
                continue
            mem_kw = _member_build_kwargs(mem, cg_coords, vdm_bb_coords,
                include_full_ligand, include_backbone_carbonyl, include_sidechain)
            ag, _ = vdg_npz.build_vdg_atomgroup_from_npz(
                **mem_kw, parent_struct=mem_struct)
            ag = align_member_to_centroid(ag, cg_coords, vdm_bb_coords, cg_vdmbb_cent,
                weights, cg_automorphisms)
            write_pdb_gz(ag, names.claim(f"{frag_name}_{make_base_tag(mem_kw)}",
                                         out_dir_clus))
            total_written += 1

        n_expected = 1 + n_members_total
        if n_expected != cluster_size:
            print(f"[WARNING] Cluster {cid}: cluster_size={cluster_size} but the npz "
                  f"holds {n_expected} record(s) for it (1 centroid + "
                  f"{n_members_total} non-centroid).")
        if clus_skipped:
            print(f"Cluster {cid}: {clus_skipped} member(s) could not be re-derived "
                  "from their parent PDB and were skipped.")

    return total_written, total_skipped


def _one(values, flag):
    """The single value of a list-valued flag, or an error naming the flag."""
    if values is None or len(values) != 1:
        got = "nothing" if values is None else f"{len(values)} values"
        raise ValueError(f"{flag} must name exactly one value here (got {got}).")
    return values[0]


def _check_pdb_dir(args):
    """Reject a --pdb-dir that can do nothing, or that points nowhere. Both fail
    silently otherwise: an unused override does nothing, and a mistyped one only warns
    per parent, writing the vdGs without the atoms it was passed for."""
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
    """Error on a flag from the mode not in effect, instead of accepting and ignoring
    it: a run that looks configured and is not."""
    centroids_only = (("--min-cluster-size", args.min_cluster_size),
                      ("--max-files", args.max_files))
    members_only = (("--reps", args.reps), ("--seed", args.seed))
    bad, mode, other = ((members_only, "centroids mode", "--members")
                        if not args.members else
                        (centroids_only, "members mode (--members)", "the default mode"))
    passed = [flag for flag, value in bad if value is not None]
    if passed:
        raise ValueError(f"{', '.join(passed)} has no effect in {mode}; "
                         f"it only applies to {other}.")


def main():
    args = parse_args()
    if args.clusters is not None and args.top_clusters is not None:
        raise ValueError("Pass either --clusters or --top-clusters, not both.")
    _reject_other_mode_flags(args)  # before value checks, so the error names the mode
    _check_pdb_dir(args)
    if args.reps is not None and args.reps < 0:
        raise ValueError("--reps must be 0 (all) or a positive count.")
    # Enforced after each write, so a non-positive value reads as "no cap".
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
            pdb_dir=args.pdb_dir,
        )
        msg = f"Wrote {n} PDB file(s)."
        if skipped:
            msg += (f" {skipped} member(s) skipped; see the [WARNING] line(s) above "
                    "for the atom and parent PDB in each case.")
        print(msg)
        return

    if args.clusters is not None:
        # Same reason as members mode: ids only mean something within one npz.
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
                         else args.align_cg_weight),
    )
    print(f"Wrote {n} PDB file(s).")


if __name__ == "__main__":
    vdg_npz.run_cli(main)
