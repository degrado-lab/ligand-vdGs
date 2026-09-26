"""CLI wrapper for vdG hit-finding. Core logic lives in hit_finder_core.py,
which downstream tools should import directly."""

import os
import argparse
import time
import traceback
import multiprocessing as mp
from contextlib import redirect_stdout, redirect_stderr

import prody as pr

from ligand_vdgs.functions import ligand_structure
from ligand_vdgs.functions.utils import convert_time_elapsed
from ligand_vdgs.score_poses.hit_finder_output import (
    RESULT_FIELDS, _describe_rmsd_threshold, write_results as _write_results)

from ligand_vdgs.score_poses.hit_finder_core import (MATCH_MODES, init_worker, lib_entries,
    process_work_item)

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

# ------------------- main ---------------- #

def _load_ref_ligand(ref_path, lig_smiles):
    """Return the reference ligand mol, or None (a bad reference must not block the
    query models; ligand RMSD is simply omitted)."""
    parse = pr.parseCIF if ref_path.endswith((".cif", ".cif.gz")) else pr.parsePDB
    try:
        # get_query_ligand_mol reports RDKit parse failures by returning None.
        mol = ligand_structure.get_query_ligand_mol(parse(ref_path), lig_smiles)
    except (ValueError, TypeError) as e:
        mol = None
        print(f"[WARNING] Failed to extract reference ligand from {ref_path}: {e}; "
              "ligand RMSD will not be computed.", flush=True)
    if mol is None:
        print(f"[WARNING] Failed to build RDKit reference ligand from {ref_path}; "
              "ligand RMSD will not be computed.", flush=True)
    return mol

def main(argv=None):
    p = argparse.ArgumentParser(prog="vdg-hit-finder",
                                description="Find vdG hits for a ligand across a directory of models.")
    p.add_argument("--smiles", required=True, help="Ligand SMILES")
    p.add_argument("--query-dir", required=True, help="Directory of query PDB/CIF models")
    p.add_argument("--vdg-lib-dir", required=True, help="vdG library directory")
    p.add_argument("--rmsd-threshold", type=float, default=None,
                    help="Max bb+CG RMSD (Å) per hit. Default: derived from the total number "
                         "of bb+CG atoms via normalize_rmsd() (same fn used during library clustering).")
    p.add_argument(
        "--match-mode", choices=MATCH_MODES, default="joint",
        help="joint: match on bb+CG RMSD (pose scoring; slot gate: virtual CB + Pro N on the query CG). bb: placement mode, match on "
             "backbone alone at tau*sqrt(n_atoms/N_bb); reads neither the query CG nor side "
             "chains (slot gate: virtual CB + Pro N on the placed CG). BSR combos still come "
             "from residues within 4.5 Å of the query ligand (default: joint).")
    p.add_argument("--ref-pdb", help="Ground-truth PDB/CIF with the same ligand; computes ligand RMSD.")
    p.add_argument("--outdir", help="Output dir (default: ./vdg-hits/<basename(query_dir)>)")
    p.add_argument("--nprocs", type=int, default=4, help="Worker processes, parallel over models (default: 4)")
    p.add_argument("--print-bsr-selection", action="store_true",
                    help="Print ProDy binding-site selection for each model.")
    p.add_argument("--contact-cutoff", type=float, default=None,
                    help="Joint mode only: skip BSR combos with no backbone N/CA/C or virtual-CB atom "
                         "within this distance (Å) of any query CG atom. Disabled by default; 3.8 is a reasonable value.")
    args = p.parse_args(argv)
    if args.match_mode == "bb" and args.contact_cutoff is not None:
        p.error("[ERROR] --contact-cutoff reads the query CG; it is joint-mode only.")

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
    vdg_lib_entries = lib_entries(args.vdg_lib_dir)

    work = [(pdbfile, os.path.join(args.query_dir, pdbfile), args.smiles, args.vdg_lib_dir,
             args.rmsd_threshold, ref_lig_mol, bool(args.print_bsr_selection),
             args.contact_cutoff, args.match_mode) for pdbfile in pdbs]

    log_path = os.path.join(outdir, "hit_finder_log.txt")
    start_time = time.perf_counter()
    rmsd_records, all_matches = [], []
    all_frags_in_lib = {}

    with open(log_path, "w") as log_fh, redirect_stdout(log_fh), redirect_stderr(log_fh):
        print("Config:")
        for label, value in (
            ("ligand_smiles", args.smiles),
            ("query_dir", args.query_dir),
            ("vdg_lib_dir", args.vdg_lib_dir),
            ("outdir", outdir),
            ("rmsd_threshold", f"{_describe_rmsd_threshold(args.rmsd_threshold)} "
                               f"({'derived' if args.rmsd_threshold is None else 'fixed'})"),
            ("contact_cutoff", f"{args.contact_cutoff} (None = disabled)"),
            ("match_mode", args.match_mode),
            ("nprocs", nprocs)):
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
            ligand_structure.summarize_frags(all_frags_in_lib, frags_to_exclude=[],
                                             frags_to_include="all", logfile_fh=log_fh)

        h, m, s = convert_time_elapsed(time.perf_counter() - start_time)
        print(f"\nTotal job time: {h} h, {m} mins, and {round(s)} secs.", flush=True)

    _write_results(outdir, all_matches, rmsd_records, args.rmsd_threshold)

if __name__ == "__main__":
    main()
