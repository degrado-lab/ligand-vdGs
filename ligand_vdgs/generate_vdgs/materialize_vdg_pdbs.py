"""Materialize centroid or member PDBs from a fragment's ``nr_vdgs`` buckets.

Centroids share one frame; members are re-derived below
``<subset>/<charge_sign>/<bucket>/clus<id>_size<size>``.
"""
import argparse
import os
import numpy as np

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions.ligand_structure import check_vdg_job_status
from ligand_vdgs.functions.vdg_materialization import (
    RES_FIELDS, best_fit as _best_fit, cap_across_buckets, move as _move,
    place_on_global as _place_on_global, rank_clusters, weights as _weights)
from ligand_vdgs.functions.vdg_pdb_io import (
    NameRegistry, format_residue_tag, fresh_dir, sanitize, strip_known_exts, write_pdb_gz)

DEFAULTS = dict(min_cluster_size=10, reps=20, align_cg_weight=0.99, seed=0)
# Parent-PDB extras: builder kwarg -> CLI flag.
EXTRA_FLAGS = {"include_full_ligand": "--include-full-ligand",
               "include_backbone_carbonyl": "--carbonyl", "include_sidechain": "--sidechain"}
PARENT_READING_FLAGS = (*EXTRA_FLAGS.values(), "--members")
MEMBER_FIELDS = ("cg_names", "cg_elements", "cg_seg", "cg_chain", "cg_resnum", "cg_resname",
                 "scrr_seg", "scrr_chain", "scrr_resnum", "scrr_resname")
REDERIVE_FIELDS = ("pdbpath", "cg_seg", "cg_chain", "cg_resnum", "cg_names",
                   "scrr_seg", "scrr_chain", "scrr_resnum")

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--fragment-dir", required=True,
                   help="A fragment's library directory (the one holding nr_vdgs/).")
    p.add_argument("-o", "--out-dir", required=True, help="Output directory.")
    p.add_argument("--include-full-ligand", action="store_true", help="Include full ligand.")
    p.add_argument("--carbonyl", action="store_true",
                   help="Also write each vdM's backbone carbonyl O, re-derived from the parent "
                        "PDB. Display only: hit-finding RMSD uses the stored N/CA/C.")
    p.add_argument("--sidechain", action=argparse.BooleanOptionalAction, default=True,
                   help="Write each vdM's sidechain heavy atoms, re-derived from the parent "
                        "PDB. On by default, but falls back to CG + bb-only.")
    p.add_argument("-P", "--pdb-dir", default=None,
                   help="RCSB-style parent database, overriding the build-time directory "
                        f"recorded in the buckets. Only used by {', '.join(PARENT_READING_FLAGS)}.")
    sel = p.add_argument_group("selection (both modes)")
    sel.add_argument("--members", action="store_true",
        help="Write every member of each selected cluster instead of just its centroid. "
             "Needs exactly one --aa-buckets and one --subset-sizes value, and a reachable "
             "parent PDB database: members store no coordinates.")
    sel.add_argument("--aa-buckets", nargs="+", default=None,
        help="AA bucket label(s) to materialize. Default: all buckets.")
    sel.add_argument("--subset-sizes", type=int, nargs="+", default=None,
        help="Subset size(s) to walk. Default: all present.")
    sel.add_argument("-n", "--top-clusters", type=int, default=None,
        help="Per (charge sign, bucket), take only the N best-supported clusters "
             "(ranked by cluster_num_parents, ties in stored order). Default: all.")
    sel.add_argument("-k", "--clusters", type=int, nargs="+", default=None,
        help="Take these cluster ID(s) instead of ranking. Ids are unique within a "
             "charge-sign bucket and apply to every matching sign. Needs exactly one "
             "--aa-buckets and one --subset-sizes value.")
    sel.add_argument("-w", "--align-cg-weight", type=float, default=None,
        help="Weight on the CG block, the rest going to the vdM backbone, when fitting "
             "onto a reference. Default: 0.99.")
    g1 = p.add_argument_group("centroids mode only (the default mode)")
    g1.add_argument("--min-cluster-size", type=int, default=None,
        help="Skip clusters with fewer than this many members (inclusive). Default: 10; "
             "pass 1 for every cluster.")
    g1.add_argument("--max-files", type=int, default=None,
        help="Global cap: after the per-bucket selection, keep the MAX_FILES best-supported "
             "clusters (by cluster_num_parents) across all charge signs, buckets, and subset "
             "sizes. Default: no cap.")
    g2 = p.add_argument_group("members mode only (--members)")
    g2.add_argument("-m", "--reps", type=int, default=None,
        help="Max structures per cluster, counting the centroid: the centroid plus REPS-1 "
             "members sampled at random. 0 = all. Default: 20.")
    g2.add_argument("--seed", type=int, default=None,
        help="Seed for the --reps sample. Default: 0.")
    return p.parse_args()

def open_fragment(frag_dir):
    """Return (nr_vdgs root, sanitized fragment name, CG automorphisms); warn on a partial build."""
    frag_dir = os.path.abspath(frag_dir)
    nr_root = os.path.join(frag_dir, "nr_vdgs")
    if not os.path.isdir(nr_root):
        raise ValueError(f"[ERROR] No nr_vdgs/ under {frag_dir!r}; -c takes a fragment directory.")
    name, lib_dir = os.path.basename(frag_dir), os.path.dirname(frag_dir)
    if not check_vdg_job_status(name, lib_dir):
        print(f"[WARNING] Fragment {name!r} in {lib_dir} has no completed vdG-generation job "
              "(no 'Job completed.' in its log); its nr_vdgs/ may be a partial build. "
              "Materializing anyway.")
    return nr_root, sanitize(name), vdg_npz.load_cg_symmetry(frag_dir, frag_name=None)[1]

def make_base_tag(kw):
    """Parent stem, then each vdM residue tag, then the ligand residue tag."""
    return "_".join([sanitize(strip_known_exts(kw["parent_pdb_path"])),
                     *(format_residue_tag(*r) for r in zip(*(kw["scrr_" + f] for f in RES_FIELDS))),
                     format_residue_tag(*(kw["cg_" + f] for f in RES_FIELDS))])

def _bucket_paths(nr_root, subset, aa_buckets=None):
    """Yield (charge_sign, bucket label, npz path) for one subset size."""
    for sign in vdg_npz.CHARGE_SIGNS:
        sign_dir = os.path.join(nr_root, subset, sign)
        if not os.path.isdir(sign_dir):
            continue
        yield from ((sign, f[:-4], os.path.join(sign_dir, f)) for f in sorted(os.listdir(sign_dir))
                    if f.endswith(".npz") and (aa_buckets is None or f[:-4] in aa_buckets))

def _usable_extras(extras, data=None, pdb_dir=None):
    """Turn every parent-read extra off if the parent database is unreachable."""
    flags = [EXTRA_FLAGS[k] for k, on in extras.items() if on]
    if flags and not vdg_npz.parent_extras_available(flags, data, pdb_dir=pdb_dir):
        return dict.fromkeys(extras, False)
    return extras

def materialize_nr_vdgs(frag_dir, out_root, max_files=None, aa_buckets=None, subset_sizes=None,
                        top_clusters=None, cluster_ids=None,
                        min_cluster_size=DEFAULTS["min_cluster_size"],
                        align_cg_weight=DEFAULTS["align_cg_weight"], pdb_dir=None, **extras):
    """Write selected cluster centroids into one flat directory; return the file count.
    ``extras`` are the EXTRA_FLAGS builder kwargs."""
    nr_root, frag_name, perms = open_fragment(frag_dir)
    out_root = os.path.abspath(out_root)
    extras = _usable_extras(extras, pdb_dir=pdb_dir)
    fresh_dir(out_root, set())
    aa_buckets = None if aa_buckets is None else set(aa_buckets)
    subsets = sorted(s for s in os.listdir(nr_root) if os.path.isdir(os.path.join(nr_root, s))
                     and (subset_sizes is None or s in {str(x) for x in subset_sizes}))
    # Pass 1 picks clusters per bucket; --max-files then caps across all buckets by
    # support, so iteration order (e.g. neg before neut) cannot decide the cut.
    picks, supports = [], []
    for subset in subsets:
        for sign, label, npz_path in _bucket_paths(nr_root, subset, aa_buckets):
            data = vdg_npz.load_bucket_npz(npz_path)
            if data is None:
                print(f"[WARNING] Skipping unreadable bucket {npz_path}.")
                continue
            ids, support = data["cluster_id"], data["cluster_num_parents"]
            if cluster_ids is None:
                order = rank_clusters(support, data["cluster_size"], min_cluster_size, top_clusters)
            else:
                order = [i for i in rank_clusters(support, data["cluster_size"])
                         if int(ids[i]) in cluster_ids]
                if missing := set(cluster_ids) - {int(ids[i]) for i in order}:
                    print(f"[WARNING] Cluster id(s) {sorted(missing)} not found in {npz_path}; "
                          "skipping.")
            picks.append((sign, label, npz_path, [int(i) for i in order]))
            supports.append([int(support[i]) for i in order])
    kept = cap_across_buckets(supports, max_files)
    names, global_ref, checked = NameRegistry(), None, False
    for (sign, label, npz_path, order), n_keep in zip(picks, kept):
        if not n_keep:
            continue
        data = vdg_npz.load_bucket_npz(npz_path)
        if not checked:
            extras, checked = _usable_extras(extras, data, pdb_dir), True
        n_cg, bucket_ref = data["nr_cg_coords"].shape[1], None
        for idx in order[:n_keep]:
            kw = vdg_npz.nr_build_kwargs(data, idx, pdb_dir=pdb_dir, **extras)
            ag, coords = vdg_npz.build_vdg_atomgroup_from_npz(**kw)
            cg, bb = kw["cg_coords"], kw["vdm_bb_coords"]
            # A bucket's first centroid is fitted on the CG alone; later ones on CG + bb.
            if bucket_ref is None:
                coords, global_ref = _place_on_global(ag, cg, coords, global_ref, perms)
                bucket_ref = coords.copy(), _weights(n_cg, np.size(bb) // 3, align_cg_weight)
            else:
                coords = _move(ag, coords, *_best_fit(cg, bb, *bucket_ref, perms))
            write_pdb_gz(ag, names.claim(
                f"{frag_name}_{sign}_{sanitize(label)}_clus{int(data['cluster_id'][idx])}_"
                f"size{int(data['cluster_size'][idx])}_{make_base_tag(kw)}", out_root))
    return sum(kept)

def materialize_cluster_members(frag_dir, out_root, subset_size, aa_bucket, cluster_ids=None,
                                top_clusters=None, reps=DEFAULTS["reps"], seed=DEFAULTS["seed"],
                                align_cg_weight=DEFAULTS["align_cg_weight"], pdb_dir=None,
                                **extras):
    """Write each selected cluster's centroid and sampled members; return (written, skipped)."""
    if (cluster_ids is None) == (top_clusters is None):
        raise ValueError("Pass exactly one of cluster_ids or top_clusters.")
    nr_root, frag_name, perms = open_fragment(frag_dir)
    bucket_paths = list(_bucket_paths(nr_root, str(subset_size), {aa_bucket}))
    if not bucket_paths:
        raise FileNotFoundError(f"No {aa_bucket}.npz under {nr_root}/{subset_size}/<charge_sign>/")
    made_dirs, names, rng, global_ref = set(), NameRegistry(), np.random.default_rng(seed), None
    written = skipped = 0
    for sign, _, npz_path in bucket_paths:
        data = vdg_npz.load_bucket_npz(npz_path)
        if data is None:
            raise ValueError(f"Bucket npz is unreadable or corrupt: {npz_path}")
        vdg_npz.require_parent_pdb_dir(data, pdb_dir=pdb_dir)
        ids = cluster_ids if top_clusters is None else [
            int(data["cluster_id"][i]) for i in
            rank_clusters(data["cluster_num_parents"], data["cluster_size"], top_n=top_clusters)]
        out_dir = os.path.join(os.path.abspath(out_root), sanitize(str(subset_size)), sign,
                               sanitize(aa_bucket))
        fresh_dir(out_dir, made_dirs)
        for cid in ids:
            w, s, global_ref = _write_cluster(data, npz_path, cid, out_dir, frag_name, names, rng,
                                              reps, global_ref, perms, align_cg_weight, pdb_dir,
                                              extras)
            written, skipped = written + w, skipped + s
    return written, skipped

def _write_cluster(data, npz_path, cid, out_dir, frag_name, names, rng, reps, global_ref, perms,
                   align_cg_weight, pdb_dir, extras):
    """Write one cluster's centroid plus up to reps-1 members fitted onto it.
    Returns (written, skipped, global CG reference)."""
    hits = np.nonzero(data["cluster_id"] == cid)[0]
    if not hits.size:
        print(f"[WARNING] Cluster {cid} not found in {npz_path}; skipping.")
        return 0, 0, global_ref
    idx = int(hits[0])
    size = int(data["cluster_size"][idx])
    kw = vdg_npz.nr_build_kwargs(data, idx, pdb_dir=pdb_dir, **extras)
    ag, cent = vdg_npz.build_vdg_atomgroup_from_npz(**kw)
    cent, global_ref = _place_on_global(ag, kw["cg_coords"], cent, global_ref, perms)
    weights = _weights(len(kw["cg_coords"]), np.size(kw["vdm_bb_coords"]) // 3, align_cg_weight)
    clus_dir = os.path.join(out_dir, f"clus{cid}_size{size}")
    os.makedirs(clus_dir, exist_ok=True)
    write_pdb_gz(ag, names.claim(f"{frag_name}_NR_clus{cid}_size{size}_{make_base_tag(kw)}",
                                 clus_dir))
    rows = vdg_npz.cluster_member_indices(data, cid)
    n_members = len(rows)
    if reps and n_members > reps - 1:
        rows = rows[np.sort(rng.choice(n_members, size=reps - 1, replace=False))]
    members = vdg_npz.load_cluster_members(data, cid, pdb_dir=pdb_dir, indices=rows)
    skipped = 0
    for mem in members:
        struct = vdg_npz.parse_pdb_or_none(mem["pdbpath"], "this member is skipped")
        cg, bb = (None, None) if struct is None else vdg_npz.rederive_member_coords(
            *(mem[k] for k in REDERIVE_FIELDS), parsed_pdb=struct)
        if cg is None:
            skipped += 1
            continue
        mem_kw = dict({k: mem[k] for k in MEMBER_FIELDS}, cg_coords=cg, vdm_bb_coords=bb,
                      parent_pdb_path=mem["pdbpath"], **extras)
        ag, _ = vdg_npz.build_vdg_atomgroup_from_npz(**mem_kw, parent_struct=struct)
        vdg_npz.apply_rigid_transform(ag, *_best_fit(cg, bb, cent, weights, perms))
        write_pdb_gz(ag, names.claim(f"{frag_name}_{make_base_tag(mem_kw)}", clus_dir))
    if 1 + n_members != size:
        print(f"[WARNING] Cluster {cid}: cluster_size={size} but the npz holds {1 + n_members} "
              f"record(s) for it (1 centroid + {n_members} non-centroid).")
    if skipped:
        print(f"Cluster {cid}: {skipped} member(s) could not be re-derived from their parent "
              "PDB and were skipped.")
    return 1 + len(members) - skipped, skipped, global_ref

def _one(values, flag):
    if values is None or len(values) != 1:
        raise ValueError(f"{flag} must name exactly one value here (got "
                         f"{'nothing' if values is None else f'{len(values)} values'}).")
    return values[0]

def _validate(args):
    """Reject flag combinations that would be silently ignored or are out of range."""
    if args.clusters is not None and args.top_clusters is not None:
        raise ValueError("Pass either --clusters or --top-clusters, not both.")
    mode, other = (("members mode (--members)", "the default mode") if args.members
                   else ("centroids mode", "--members"))
    if passed := [flag for flag, value in (
            (("--min-cluster-size", args.min_cluster_size), ("--max-files", args.max_files))
            if args.members else (("--reps", args.reps), ("--seed", args.seed)))
            if value is not None]:
        raise ValueError(f"{', '.join(passed)} has no effect in {mode}; "
                         f"it only applies to {other}.")
    if args.pdb_dir is not None:
        if not (args.include_full_ligand or args.carbonyl or args.sidechain or args.members):
            raise ValueError(f"--pdb-dir has no effect without one of "
                             f"{', '.join(PARENT_READING_FLAGS)}; every other mode is served "
                             "from the bucket npz alone.")
        if not os.path.isdir(args.pdb_dir):
            raise ValueError(f"--pdb-dir is not a directory: {args.pdb_dir}")
    if args.reps is not None and args.reps < 0:
        raise ValueError("--reps must be 0 (all) or a positive count.")
    if args.max_files is not None and args.max_files < 1:
        raise ValueError("--max-files must be a positive count (omit it for no cap).")
    if args.align_cg_weight is not None and not 0.0 <= args.align_cg_weight <= 1.0:
        raise ValueError("--align-cg-weight must be between 0 and 1.")
    if args.members and args.clusters is None and args.top_clusters is None:
        raise ValueError("--members needs --clusters or --top-clusters to say which "
                         "cluster(s) to write the members of.")
    if not args.members and args.clusters is not None:
        _one(args.subset_sizes, "--subset-sizes")
        _one(args.aa_buckets, "--aa-buckets")
        if args.min_cluster_size is not None:
            raise ValueError("--min-cluster-size has no effect alongside --clusters, which "
                             "names the clusters to write outright.")

def main():
    args = parse_args()
    _validate(args)
    for key, value in DEFAULTS.items():
        if getattr(args, key) is None:
            setattr(args, key, value)
    common = dict(frag_dir=args.fragment_dir, out_root=args.out_dir,
                  top_clusters=args.top_clusters, cluster_ids=args.clusters,
                  align_cg_weight=args.align_cg_weight, pdb_dir=args.pdb_dir,
                  include_full_ligand=args.include_full_ligand,
                  include_backbone_carbonyl=args.carbonyl, include_sidechain=args.sidechain)
    if not args.members:
        n = materialize_nr_vdgs(max_files=args.max_files, aa_buckets=args.aa_buckets,
                                subset_sizes=args.subset_sizes,
                                min_cluster_size=args.min_cluster_size, **common)
        print(f"Wrote {n} PDB file(s).")
        return
    n, skipped = materialize_cluster_members(
        subset_size=_one(args.subset_sizes, "--subset-sizes"),
        aa_bucket=_one(args.aa_buckets, "--aa-buckets"), reps=args.reps, seed=args.seed, **common)
    print(f"Wrote {n} PDB file(s)." + (
        f" {skipped} member(s) skipped; see the [WARNING] line(s) above for the atom and "
        "parent PDB in each case." if skipped else ""))

if __name__ == "__main__":
    vdg_npz.run_cli(main)
