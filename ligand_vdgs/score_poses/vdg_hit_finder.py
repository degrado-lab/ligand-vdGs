"""CLI wrapper for vdG hit-finding. Core logic lives in hit_finder_core.py,
which downstream tools should import directly."""

import os
import csv
import argparse
import time
import traceback
import multiprocessing as mp
from contextlib import redirect_stdout, redirect_stderr

import pandas as pd
import prody as pr

from ligand_vdgs.functions import Frags
from ligand_vdgs.functions.utils import convert_time_elapsed

from ligand_vdgs.score_poses.hit_finder_core import init_worker, process_work_item

RESULT_FIELDS = [
    "pdbfile", "struct_id", "frag", "query_frag", "subset_size", "bsr_combo",
    "aa_bucket", "vdg_index", "vdg_cluster_id", "vdg_cluster_size", "vdg_rmsd",
    "rmsd_threshold", "aa_perm_idx", "q_site_idx", "q_cg_perm_idx", "q_atom_indices",
    "R00", "R01", "R02", "R10", "R11", "R12", "R20", "R21", "R22", "t0", "t1", "t2",
]


def _collect_pool_results(work, async_results):
    """Collect pool tasks independently so one raised task does not lose the rest."""
    results = []
    for work_item, async_result in zip(work, async_results):
        try:
            results.append(async_result.get())
        except Exception:
            pdbfile, pdb_path = work_item[:2]
            error_text = (f"Pool task failed for pdbfile={pdbfile}, pdb_path={pdb_path}\n"
                          f"{traceback.format_exc()}")
            results.append((pdbfile, "", None, [], {}, error_text))
    return results


# -------------- write results ----------- #

def _sample_name_from_pdbfile(pdbfile):
    base = os.path.basename(str(pdbfile))
    for ext in (".pdb.gz", ".pdb", ".cif.gz", ".cif", ".gz"):
        if base.endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]


def _describe_rmsd_threshold(rmsd_threshold):
    if rmsd_threshold is None:
        return 'per-combo normalize_rmsd(n_atoms, "cgvdmbb")'
    return str(rmsd_threshold)


def _write_results(outdir, all_matches, rmsd_records, rmsd_threshold):
    os.makedirs(outdir, exist_ok=True)

    if rmsd_records:
        rmsd_path = os.path.join(outdir, "ligand_rmsd_summary.tsv")
        with open(rmsd_path, "w") as fh:
            fh.write("pdbfile\tligand_rmsd\n")
            for pdbfile, rmsd in sorted(rmsd_records):
                fh.write(f"{pdbfile}\t{rmsd:.3f}\n")
        print(f"Wrote ligand RMSD summary to: {rmsd_path}")

    results_path = os.path.join(outdir, "vdg_hits.tsv")
    with open(results_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=RESULT_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(all_matches)
    print(f"Wrote vdG match results to:   {results_path}")

    if not all_matches:
        return

    df = pd.crosstab(
        pd.Series([_sample_name_from_pdbfile(r["pdbfile"]) for r in all_matches], name="pdb"),
        pd.Series([r["frag"] for r in all_matches]),
    ).sort_index(axis=0).sort_index(axis=1)

    txt_path = os.path.join(outdir, "results_summary.txt")
    with open(txt_path, "w") as fh:
        fh.write(f"vdG hit counts per sample (rows) and fragment (columns); counts use "
                 f"rmsd_threshold={_describe_rmsd_threshold(rmsd_threshold)}\n\n")
        fh.write(df.to_string())
        fh.write("\n")
    print(f"Wrote summary table to:        {txt_path}")


# ------------------- main ---------------- #

def _load_ref_ligand(ref_path, lig_smiles):
    """Return the reference ligand mol, or None (a bad reference must not block the
    query models; ligand RMSD is simply omitted)."""
    parse = pr.parseCIF if ref_path.endswith((".cif", ".cif.gz")) else pr.parsePDB
    try:
        # get_query_ligand_mol reports RDKit parse failures by returning None.
        mol = Frags.get_query_ligand_mol(parse(ref_path), lig_smiles)
    except (ValueError, TypeError) as e:
        mol = None
        print(f"[WARNING] Failed to extract reference ligand from {ref_path}: {e}; "
              "ligand RMSD will not be computed.", flush=True)
    if mol is None:
        print(f"[WARNING] Failed to build RDKit reference ligand from {ref_path}; "
              "ligand RMSD will not be computed.", flush=True)
    return mol


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="vdg-hit-finder",
        description="Find vdG hits for a ligand across a directory of models.",
    )
    p.add_argument("--smiles", required=True, help="Ligand SMILES")
    p.add_argument("--query-dir", required=True, help="Directory of query PDB/CIF models")
    p.add_argument("--vdg-lib-dir", required=True, help="vdG library directory")
    p.add_argument(
        "--rmsd-threshold", type=float, default=None,
        help="Max bb+CG RMSD (Å) per hit. Default: derived from the total number "
             "of bb+CG atoms via normalize_rmsd() (same fn used during library clustering).",
    )
    p.add_argument("--ref-pdb", help="Ground-truth PDB/CIF with the same ligand; computes ligand RMSD.")
    p.add_argument("--outdir", help="Output dir (default: ./vdg-hits/<basename(query_dir)>)")
    p.add_argument("--nprocs", type=int, default=4, help="Worker processes, parallel over models (default: 4)")
    p.add_argument("--print-bsr-selection", action="store_true",
                    help="Print ProDy binding-site selection for each model.")
    p.add_argument(
        "--contact-cutoff", type=float, default=None,
        help="Skip BSR combos with no atom within this distance (Å) of any CG atom. "
             "Disabled by default; 3.8 is a reasonable value.",
    )
    p.add_argument(
        "--no-dedup", dest="deduplicate", action="store_false",
        help="Disable hit deduplication (default: hits sharing a BSR combo and at least "
             "--min-shared-atoms query-ligand atoms collapse to the best-RMSD hit).",
    )
    p.add_argument(
        "--min-shared-atoms", type=int, default=3,
        help="Shared query-ligand atom count to call two hits duplicates (default: 3).",
    )
    args = p.parse_args(argv)

    outdir = args.outdir or os.path.join(
        os.getcwd(), "vdg-hits", os.path.basename(os.path.normpath(args.query_dir)))
    nprocs = max(1, args.nprocs)
    os.makedirs(outdir, exist_ok=True)

    pdbs = [f for f in sorted(os.listdir(args.query_dir))
            if f.lower().endswith((".pdb", ".pdb.gz", ".cif", ".cif.gz"))]
    if not pdbs:
        print("No PDB/CIF files found in query_dir.")
        return

    ref_lig_mol = _load_ref_ligand(args.ref_pdb, args.smiles) if args.ref_pdb else None
    vdg_lib_entries = set(os.listdir(args.vdg_lib_dir))

    work = [
        (pdbfile, os.path.join(args.query_dir, pdbfile), args.smiles, args.vdg_lib_dir,
         args.rmsd_threshold, ref_lig_mol, bool(args.print_bsr_selection),
         args.contact_cutoff, args.deduplicate, args.min_shared_atoms)
        for pdbfile in pdbs
    ]

    log_path = os.path.join(outdir, "hit_finder_log.txt")
    start_time = time.perf_counter()
    rmsd_records, all_matches = [], []
    all_frags_in_lib = {}

    with open(log_path, "w") as log_fh, redirect_stdout(log_fh), redirect_stderr(log_fh):
        print("Config:")
        for label, value in [
            ("ligand_smiles", args.smiles),
            ("query_dir", args.query_dir),
            ("vdg_lib_dir", args.vdg_lib_dir),
            ("outdir", outdir),
            ("rmsd_threshold", f"{_describe_rmsd_threshold(args.rmsd_threshold)} "
                               f"({'derived' if args.rmsd_threshold is None else 'fixed'})"),
            ("contact_cutoff", f"{args.contact_cutoff} (None = disabled)"),
            ("deduplicate", args.deduplicate),
            ("min_shared_atoms", args.min_shared_atoms),
            ("nprocs", nprocs),
        ]:
            print(f"  {label:<20}: {value}")
        print()

        n_workers = min(len(work), nprocs)
        if n_workers > 1:
            ctx = mp.get_context("spawn")
            with ctx.Pool(processes=n_workers, initializer=init_worker,
                          initargs=(vdg_lib_entries,)) as pool:
                async_results = [pool.apply_async(process_work_item, (w,)) for w in work]
                results = _collect_pool_results(work, async_results)
        else:
            init_worker(vdg_lib_entries)
            results = [process_work_item(w) for w in work]

        failed_models = []
        for (pdbfile, log_text, lig_rmsd_value, match_records_i,
             frags_in_lib_i, error_text) in results:
            if log_text:
                log_fh.write(log_text if log_text.endswith("\n") else log_text + "\n")
                log_fh.flush()

            if error_text is not None:
                failed_models.append(pdbfile)
                print(f"[ERROR] Model failed: {pdbfile}", flush=True)
                print(error_text, file=log_fh, flush=True)
                continue

            if lig_rmsd_value is not None:
                rmsd_records.append((pdbfile, lig_rmsd_value))
            all_matches.extend(match_records_i)
            for frag_smiles, ok in frags_in_lib_i.items():
                all_frags_in_lib[frag_smiles] = all_frags_in_lib.get(frag_smiles, False) or ok

        print(f"\nModel processing summary: total={len(results)}, "
              f"completed={len(results) - len(failed_models)}, failed={len(failed_models)}",
              flush=True)
        if failed_models:
            print("Failed models:", flush=True)
            for pdbfile in failed_models:
                print(f"  - {pdbfile}", flush=True)

        if all_frags_in_lib:
            Frags.summarize_frags(all_frags_in_lib, frags_to_exclude=[],
                                  frags_to_include="all", logfile_fh=log_fh)

        h, m, s = convert_time_elapsed(time.perf_counter() - start_time)
        print(f"\nTotal job time: {h} h, {m} mins, and {round(s)} secs.", flush=True)

    _write_results(outdir, all_matches, rmsd_records, args.rmsd_threshold)


if __name__ == "__main__":
    main()
