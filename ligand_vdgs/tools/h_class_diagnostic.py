"""h_class_diagnostic.py

Per fragment key / CG-atom slot orbit (automorphism-equivalent slots pooled)
/ AA bucket: does H-class (``num_h == 0`` vs. ``> 0``) split the observations
informatively? Run after every library (re)build and whenever the fragment
key scheme changes, before deciding which ``H0``/``!H0`` key variants to pool
at read time. Full rationale: docs/database_generation_guide.md, "When to run
the H-class diagnostic".

Reads ``{nr,mem}_cg_heavy_degree`` / ``{nr,mem}_cg_num_h`` (int8, n_rows x
n_cg_atoms) written alongside the usual per-row arrays. Sentinel contract: a
NEGATIVE value means "unreadable"; 0 is a real count, so ``heavy_degree == 0``
(impossible for a connected atom) is a hard error, and ``heavy_degree < 0``
with ``num_h >= 0`` only warns. A bucket violating this aborts the run;
``--dry-run`` only runs this check, over the whole library. Optional
``nr_vdm_o_coords`` (backbone carbonyl O, never part of RMSD) lets
``--contact-atoms`` include ``O``.

Per orbit/bucket, measures (a) class balance (nr+mem observations, clusters,
parents) per H-class; (b) contact rate by class (nr rows only: CG atom within
``--cutoff`` of a stored backbone atom), with a cluster-block bootstrap CI on
the H-vs-noH difference; (c) verdict -- ``single-class`` (a class <
``--min-class-frac`` of known pairs) > ``low-support`` (< ``--min-support``
nr rows, or too few valid bootstrap reps) > ``informative``/``uninformative``
by whether the CI excludes 0.

Output: a TSV, one row per key/slot orbit/bucket/class, plus a stdout
summary. Bucket ``X`` (non-canonical) is never read.

Usage
-----
    python ligand_vdgs/tools/h_class_diagnostic.py --lib /path/to/frag_lib \\
        --out h_class.tsv [--keys keys.txt] [--buckets bb ASP GLU] \\
        [--subset-size 1] [--cutoff 3.5] [--contact-atoms N,CA,C,O]
    python ligand_vdgs/tools/h_class_diagnostic.py --lib /path/to/frag_lib --dry-run
"""

import argparse
import itertools
import os
import sys

import numpy as np

from ligand_vdgs.functions import parent_db, utils
from ligand_vdgs.functions.Frags import check_vdg_job_status
from ligand_vdgs.functions.vdg_npz_utils import (
    cg_symmetry_path, load_bucket_npz, load_cg_symmetry, load_vdg_bucket,
    make_aa_bucket, vdg_npz_path)
from ligand_vdgs.functions.vdg_struct_utils import (
    BB_LABEL, CANONICAL_HEAVY_ATOMS, NONCANONICAL_AA_LABEL)

H_CLASS_FIELDS = ("nr_cg_heavy_degree", "nr_cg_num_h",
                  "mem_cg_heavy_degree", "mem_cg_num_h")
# Order of the stored vdM backbone triplet (see build_vdg_atomgroup_from_npz).
BB_ATOM_ORDER = ("N", "CA", "C")
SLOT_LABELS = tuple(sorted(CANONICAL_HEAVY_ATOMS)) + (BB_LABEL,)
O_COORDS_FIELD = "nr_vdm_o_coords"   # optional (n_nr, num_vdms, 3)

TSV_COLUMNS = (
    "key", "cg_smarts", "subset_size", "bucket", "slot_orbit", "orbit_elements",
    "h_class", "n_slot_obs", "n_obs", "n_clusters", "n_parents", "frac_of_all_pairs",
    "heavy_degree_hist", "n_nr_pairs", "contact_rate", "contact_atoms", "cutoff",
    "diff_H_minus_noH", "ci_lo", "ci_hi", "n_boot_valid", "verdict", "reason")


class ContractError(RuntimeError):
    """A bucket does not carry the H-class fields this tool reads."""


# --- Library discovery ---

def completed_fragments(lib_dir, keys=None):
    """Fragment dir names with a completed job log and a symmetry sidecar.

    ``keys`` maps to directory names via ``utils.smiles_to_filename``;
    without it every directory at the library root is a candidate.
    """
    candidates = ([utils.smiles_to_filename(k) for k in keys] if keys is not None
                  else sorted(os.listdir(lib_dir)))
    done, skipped = [], []
    for frag in candidates:
        if not os.path.isdir(os.path.join(lib_dir, frag)):
            skipped.append((frag, "no such fragment directory"))
        elif not check_vdg_job_status(frag, lib_dir):
            skipped.append((frag, "no 'Job completed.' in its log"))
        elif not os.path.isfile(cg_symmetry_path(lib_dir, frag)):
            skipped.append((frag, "no nr_vdgs/cg_symmetry.npz"))
        else:
            done.append(frag)
    for frag, why in skipped:
        print(f"[WARNING] skipping {frag!r}: {why}", file=sys.stderr)
    return done


def bucket_names(subset_size):
    """Every bucket name a subset size can produce, X-free, sorted."""
    combos = itertools.combinations_with_replacement(SLOT_LABELS, subset_size)
    return sorted({make_aa_bucket(c) for c in combos})


def read_bucket(lib_dir, frag, subset_size, bucket):
    """Full array dict of one bucket, or None if absent/non-canonical/unreadable.

    ``load_vdg_bucket`` carries the job-status warning and corrupt-file
    handling but strips mem_ rows and the H-class fields, so the full dict is
    read separately once it has confirmed the file is good.
    """
    if NONCANONICAL_AA_LABEL in bucket.split("_"):
        return None
    if load_vdg_bucket(lib_dir, frag, subset_size, bucket) is None:
        return None
    return load_bucket_npz(vdg_npz_path(lib_dir, frag, subset_size, bucket))


# --- Contract check ---

def contract_problems(data, npz_path):
    """Human-readable problems with the H-class fields of one bucket ([] if ok)."""
    n_cg = int(data["cg_elements"].shape[0])
    n_nr, n_mem = len(data["cluster_id"]), len(data["mem_cluster_id"])
    expected = dict(zip(H_CLASS_FIELDS,
                        [(n_nr, n_cg), (n_nr, n_cg), (n_mem, n_cg), (n_mem, n_cg)]))
    problems = []
    for field, shape in expected.items():
        if field not in data:
            problems.append(f"{npz_path}: missing field {field}; expected int8 {shape}")
            continue
        arr = data[field]
        if arr.shape != shape:
            problems.append(f"{npz_path}: {field} has shape {arr.shape}, expected {shape}")
        if not np.issubdtype(arr.dtype, np.integer):
            problems.append(f"{npz_path}: {field} has dtype {arr.dtype}, expected an integer dtype")
        elif arr.dtype != np.int8:
            print(f"[WARNING] {npz_path}: {field} is {arr.dtype}, not int8; read anyway.",
                  file=sys.stderr)
    # 0 degree is impossible (every CG atom has >= 1 heavy neighbour), so it can only
    # be a writer storing 0 for "unreadable" instead of the required negative sentinel
    # (error). deg < 0 with num_h >= 0 is merely suspicious -- H count claimed for an
    # atom whose graph was unreadable -- so it only warns.
    for prefix in ("nr", "mem"):
        deg, nh = f"{prefix}_cg_heavy_degree", f"{prefix}_cg_num_h"
        if deg not in data or nh not in data or data[deg].shape != data[nh].shape:
            continue
        d, h = np.asarray(data[deg]), np.asarray(data[nh])
        n_zero = int((d == 0).sum())
        if n_zero:
            problems.append(f"{npz_path}: {deg} has {n_zero} entries == 0; "
                             "the contract requires a NEGATIVE sentinel for unreadable atoms.")
        n_rev = int(((d < 0) & (h >= 0)).sum())
        if n_rev:
            print(f"[WARNING] {npz_path}: {n_rev} entries have {deg} < 0 but {nh} >= 0.",
                  file=sys.stderr)
    return problems


def require_contract(data, npz_path):
    problems = contract_problems(data, npz_path)
    if problems:
        raise ContractError(
            "This library does not carry the per-observation H-class fields "
            f"({', '.join(H_CLASS_FIELDS)}) the diagnostic reads:\n  "
            + "\n  ".join(problems)
            + "\nRebuild the library with a writer that records them, or run "
              "--dry-run to list every affected bucket.")


def dry_run(lib_dir, frags, subset_size, buckets):
    """Validate the contract over every requested bucket; return #problem buckets."""
    n_checked = n_bad = 0
    for frag in frags:
        for bucket in buckets:
            data = read_bucket(lib_dir, frag, subset_size, bucket)
            if data is None:
                continue
            n_checked += 1
            problems = contract_problems(
                data, vdg_npz_path(lib_dir, frag, subset_size, bucket))
            if problems:
                n_bad += 1
                print("\n".join(problems))
    print(f"[dry-run] {n_checked} bucket(s) checked in {len(frags)} fragment(s); "
          f"{n_bad} lack or mis-shape the H-class fields "
          f"({', '.join(H_CLASS_FIELDS)}).")
    return n_bad


# --- Per-bucket statistics ---

def slot_orbits(automorphisms, n_cg):
    """Connected components of slot indices under the automorphism group."""
    parent = list(range(n_cg))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for perm in automorphisms:
        for i, j in enumerate(perm):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[max(ri, rj)] = min(ri, rj)
    orbits = {}
    for i in range(n_cg):
        orbits.setdefault(find(i), []).append(i)
    return [tuple(v) for _, v in sorted(orbits.items())]


def contact_atom_indices(contact_atoms):
    """Indices into BB_ATOM_ORDER for the non-'O' entries of --contact-atoms."""
    idx = []
    for name in contact_atoms:
        if name == "O":
            continue
        if name not in BB_ATOM_ORDER:
            raise ValueError(f"--contact-atoms must be drawn from {BB_ATOM_ORDER + ('O',)}; "
                             f"got {name!r}.")
        idx.append(BB_ATOM_ORDER.index(name))
    return idx


def nr_contacts(data, atom_idx, cutoff, want_o):
    """(n_nr, n_cg) bool: CG atom within cutoff of any selected stored residue atom."""
    cg = np.asarray(data["nr_cg_coords"], dtype=np.float64)          # (n, n_cg, 3)
    bb = np.asarray(data["nr_vdm_bb_coords"], dtype=np.float64)      # (n, v, 3, 3)
    res = bb[:, :, atom_idx, :].reshape(bb.shape[0], -1, 3)          # (n, v*a, 3)
    if want_o:
        if O_COORDS_FIELD not in data:
            raise ContractError(
                f"--contact-atoms includes O but this bucket has no {O_COORDS_FIELD} "
                "(n_nr, num_vdms, 3) field; rebuild with a writer that stores the "
                "backbone carbonyl O, or drop O from --contact-atoms.")
        res = np.concatenate([res, np.asarray(data[O_COORDS_FIELD], dtype=np.float64)], axis=1)
    d = np.linalg.norm(cg[:, :, None, :] - res[:, None, :, :], axis=-1)
    # An absent O is NaN; a plain min would turn the whole row NaN (-> False).
    return np.where(np.isnan(d), np.inf, d).min(axis=2) <= cutoff


def h_class_of(num_h):
    """1 = H, 0 = noH, -1 = unknown (negative sentinel)."""
    num_h = np.asarray(num_h)
    return np.where(num_h < 0, -1, (num_h > 0).astype(np.int8))


def bootstrap_diff(nH, cH, nN, cN, n_boot, rng, chunk=200):
    """Cluster-block bootstrap of rate_H - rate_noH.

    Inputs are per-nr-row (per-cluster) pair counts: H pairs, contacting H
    pairs, noH pairs, contacting noH pairs. Reps where a class vanishes are
    NaN and excluded from the percentile CI; their count is reported.
    """
    n = len(nH)
    counts = np.stack([nH, cH, nN, cN], axis=1).astype(np.float64)   # (n, 4)
    diffs = np.empty(n_boot)
    p = np.full(n, 1.0 / n)
    for start in range(0, n_boot, chunk):
        reps = min(chunk, n_boot - start)
        weights = rng.multinomial(n, p, size=reps).astype(np.float64)  # (reps, n)
        s = weights @ counts                                            # (reps, 4)
        with np.errstate(invalid="ignore", divide="ignore"):
            diffs[start:start + reps] = s[:, 1] / s[:, 0] - s[:, 3] / s[:, 2]
    valid = np.isfinite(diffs)
    if valid.sum() == 0:
        return np.nan, np.nan, 0
    lo, hi = np.percentile(diffs[valid], [2.5, 97.5])
    return float(lo), float(hi), int(valid.sum())


def _degree_hist(degrees):
    vals, cnts = np.unique(degrees, return_counts=True)
    return ";".join(f"{int(v)}:{int(c)}" for v, c in zip(vals, cnts)) or ""


def _fmt(x):
    return "" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.4f}"


def analyse_bucket(data, orbits, args, rng):
    """Rows (dicts, TSV_COLUMNS keys minus key/smarts/bucket) for one bucket."""
    n_nr, n_mem = len(data["cluster_id"]), len(data["mem_cluster_id"])
    elements = [str(e) for e in data["cg_elements"]]

    cls_nr = h_class_of(data["nr_cg_num_h"])                    # (n_nr, n_cg)
    cls_mem = h_class_of(data["mem_cg_num_h"])                  # (n_mem, n_cg)
    cls_all = np.concatenate([cls_nr, cls_mem], axis=0)
    deg_all = np.concatenate([data["nr_cg_heavy_degree"],
                              data["mem_cg_heavy_degree"]], axis=0)
    clusters_all = np.concatenate([data["cluster_id"], data["mem_cluster_id"]])
    parents_all = np.asarray([parent_db.entry_of(b) for b in
                              np.concatenate([data["nr_parent_biounit"],
                                              data["mem_parent_biounit"]])])
    contact = nr_contacts(data, args.contact_idx, args.cutoff, args.want_o)

    rows = []
    for orbit in orbits:
        orbit = list(orbit)
        # (observation, slot) pairs pooled over the orbit's slots.
        c_pairs = cls_all[:, orbit].ravel()
        d_pairs = deg_all[:, orbit].ravel()
        obs_idx = np.repeat(np.arange(n_nr + n_mem), len(orbit))
        n_pairs = c_pairs.size
        n_known = int((c_pairs >= 0).sum())

        # nr-only, for geometry: per-row counts feed the block bootstrap.
        c_nr = cls_nr[:, orbit]                                  # (n_nr, k)
        k_nr = contact[:, orbit]
        nH = (c_nr == 1).sum(axis=1); cH = ((c_nr == 1) & k_nr).sum(axis=1)
        nN = (c_nr == 0).sum(axis=1); cN = ((c_nr == 0) & k_nr).sum(axis=1)
        rate = {"H": cH.sum() / nH.sum() if nH.sum() else np.nan,
                "noH": cN.sum() / nN.sum() if nN.sum() else np.nan}
        nr_rows = {"H": int((nH > 0).sum()), "noH": int((nN > 0).sum())}

        frac = {"H": (c_pairs == 1).sum() / n_known if n_known else np.nan,
                "noH": (c_pairs == 0).sum() / n_known if n_known else np.nan}
        diff = ci_lo = ci_hi = np.nan
        n_valid = 0
        if n_known == 0 or min(frac.values()) < args.min_class_frac:
            small = min(frac, key=frac.get) if n_known else "both"
            verdict, reason = "single-class", (
                f"{small} < {args.min_class_frac:.0%} of known-class pairs")
        elif min(nr_rows.values()) < args.min_support:
            verdict, reason = "low-support", (
                f"nr rows H={nr_rows['H']} noH={nr_rows['noH']} < {args.min_support}")
        else:
            diff = rate["H"] - rate["noH"]
            ci_lo, ci_hi, n_valid = bootstrap_diff(nH, cH, nN, cN, args.n_boot, rng)
            if n_valid < args.n_boot // 2:
                verdict, reason = "low-support", f"only {n_valid} valid bootstrap reps"
            elif ci_lo > 0 or ci_hi < 0:
                verdict, reason = "informative", "95% CI excludes 0"
            else:
                verdict, reason = "uninformative", "95% CI spans 0"

        for label, code in (("H", 1), ("noH", 0), ("unknown", -1)):
            sel = c_pairs == code
            if label == "unknown" and not sel.any():
                continue
            rows.append(dict(
                slot_orbit="+".join(str(i) for i in orbit),
                orbit_elements="".join(elements[i] for i in orbit),
                h_class=label,
                n_slot_obs=int(sel.sum()),
                n_obs=int(np.unique(obs_idx[sel]).size),
                n_clusters=int(np.unique(clusters_all[obs_idx[sel]]).size),
                n_parents=int(np.unique(parents_all[obs_idx[sel]]).size),
                frac_of_all_pairs=(f"{sel.sum() / n_pairs:.4f}" if n_pairs else ""),
                heavy_degree_hist=_degree_hist(d_pairs[sel]),
                n_nr_pairs=int((c_nr == code).sum()),
                contact_rate=_fmt(rate.get(label, np.nan)),
                contact_atoms=",".join(args.contact_atoms),
                cutoff=f"{args.cutoff:g}",
                diff_H_minus_noH=_fmt(diff), ci_lo=_fmt(ci_lo), ci_hi=_fmt(ci_hi),
                n_boot_valid=n_valid, verdict=verdict, reason=reason,
            ))
    return rows


# --- Driver ---

def select_buckets(lib_dir, frag, subset_size, requested, min_bucket_obs):
    """(bucket, data) pairs to analyse for one fragment.

    Default: bb, ASP, GLU (subset size 1) plus any bucket with at least
    ``min_bucket_obs`` observations (nr + mem rows).
    """
    names = requested if requested else bucket_names(subset_size)
    always = set() if requested else {BB_LABEL, "ASP", "GLU"}
    out = []
    for bucket in names:
        data = read_bucket(lib_dir, frag, subset_size, bucket)
        if data is None:
            if requested:
                why = ("holds non-canonical contacts and is excluded from H-class "
                       "statistics" if NONCANONICAL_AA_LABEL in bucket.split("_")
                       else f"(subset size {subset_size}) not found")
                print(f"[WARNING] {frag}: bucket {bucket!r} {why}.", file=sys.stderr)
            continue
        n_obs = len(data["cluster_id"]) + len(data["mem_cluster_id"])
        if requested or bucket in always or n_obs >= min_bucket_obs:
            out.append((bucket, data))
    return out


def print_summary(rows):
    print("\nH-class diagnostic summary (contact = CG atom within cutoff of a stored "
          "vdM backbone atom; carbonyl O and sidechains are not in the npz)")
    by_key = {}
    for r in rows:
        by_key.setdefault(r["key"], []).append(r)
    for key in sorted(by_key):
        print(f"\n{key}")
        groups = {}
        for r in by_key[key]:
            groups.setdefault((r["bucket"], r["slot_orbit"]), {})[r["h_class"]] = r
        for tag in sorted(groups):
            h, n = groups[tag].get("H"), groups[tag].get("noH")
            parts = []
            for label, x in (("H", h), ("noH", n)):
                if x is None:
                    parts.append(f"{label} n=0")
                    continue
                rate = f" contact={x['contact_rate']}" if x["contact_rate"] else ""
                parts.append(f"{label} n={x['n_slot_obs']} ({x['n_clusters']} clus, "
                             f"{x['n_parents']} parents){rate}")
            ref = h or n
            ci = (f" diff={ref['diff_H_minus_noH']} [{ref['ci_lo']}, {ref['ci_hi']}]"
                  if ref["diff_H_minus_noH"] else "")
            print(f"  {ref['bucket']:<8} slot {ref['slot_orbit']:<7} "
                  f"{ref['orbit_elements']:<5} {' | '.join(parts)}{ci} "
                  f"-> {ref['verdict']} ({ref['reason']})")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lib", required=True, help="vdG library root")
    p.add_argument("--keys", help="file with one fragment key per line "
                   "(default: every completed fragment in --lib)")
    p.add_argument("--buckets", nargs="+", help="AA buckets to analyse (default: bb, "
                   "ASP, GLU plus any with >= --min-bucket-obs observations)")
    p.add_argument("--min-bucket-obs", type=int, default=0,
                   help="skip non-default buckets with fewer observations (0 = analyse "
                        "every bucket; only saves time)")
    p.add_argument("--subset-size", type=int, default=1)
    p.add_argument("--cutoff", type=float, default=3.5, help="contact cutoff, A")
    p.add_argument("--contact-atoms", default="N,CA,C",
                   help="stored vdM backbone atoms that count as contacts")
    p.add_argument("--min-class-frac", type=float, default=0.0,
                   help="class share below which a slot is labelled single-class "
                        "(0 = only when a class is absent; labels only, never skips)")
    p.add_argument("--min-support", type=int, default=5,
                   help="nr rows (clusters) per class below which the verdict is "
                        "low-support instead of a CI verdict (labels only, never skips)")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", help="output TSV (required unless --dry-run)")
    p.add_argument("--dry-run", action="store_true",
                   help="only check that the H-class fields are present")
    args = p.parse_args(argv)
    args.contact_atoms = [a.strip() for a in args.contact_atoms.split(",") if a.strip()]
    args.contact_idx = contact_atom_indices(args.contact_atoms)
    args.want_o = "O" in args.contact_atoms
    if not args.dry_run and not args.out:
        p.error("--out is required unless --dry-run")
    if args.subset_size not in (1, 2):
        p.error("--subset-size must be 1 or 2")
    return args


def main(argv=None):
    args = parse_args(argv)
    keys = None
    if args.keys:
        with open(args.keys) as fh:
            keys = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    frags = completed_fragments(args.lib, keys)
    if not frags:
        raise SystemExit("[ERROR] no completed fragments to analyse.")

    if args.dry_run:
        buckets = args.buckets or bucket_names(args.subset_size)
        raise SystemExit(1 if dry_run(args.lib, frags, args.subset_size, buckets) else 0)

    rng = np.random.default_rng(args.seed)
    rows = []
    for frag in frags:
        key = utils.filename_to_smiles(frag)
        cg_smarts, automorphisms = load_cg_symmetry(args.lib, frag)
        for bucket, data in select_buckets(args.lib, frag, args.subset_size,
                                           args.buckets, args.min_bucket_obs):
            npz_path = vdg_npz_path(args.lib, frag, args.subset_size, bucket)
            try:
                require_contract(data, npz_path)
            except ContractError as err:
                raise SystemExit(f"[ERROR] {err}")
            orbits = slot_orbits(automorphisms, len(data["cg_elements"]))
            for r in analyse_bucket(data, orbits, args, rng):
                r.update(key=key, cg_smarts=cg_smarts, subset_size=args.subset_size,
                         bucket=bucket)
                rows.append(r)

    rows.sort(key=lambda r: (r["key"], r["bucket"], r["slot_orbit"], r["h_class"]))
    with open(args.out, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in TSV_COLUMNS) + "\n")
    print_summary(rows)
    print(f"\nWrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
