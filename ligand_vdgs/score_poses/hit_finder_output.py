"""TSV and summary-file output for the hit-finder CLI."""

import csv
import os

import pandas as pd

RESULT_FIELDS = [
    "pdbfile", "struct_id", "lig_instance", "frag", "query_frag", "subset_size", "bsr_combo",
    "aa_bucket", "charge_sign", "vdg_index", "vdg_cluster_id", "vdg_cluster_num_parents", "vdg_rmsd",
    "rmsd_threshold", "aa_perm_idx", "q_site_idx", "q_cg_perm_idx", "q_atom_indices",
    # R/t are joint-fit transforms; Rbb/tbb fit the backbone with the CG held out.
    "match_mode", "bb_rmsd", "held_out_cg_rmsd", "cg_bb_dist",
    "R00", "R01", "R02", "R10", "R11", "R12", "R20", "R21", "R22", "t0", "t1", "t2",
    "Rbb00", "Rbb01", "Rbb02", "Rbb10", "Rbb11", "Rbb12", "Rbb20", "Rbb21", "Rbb22",
    "tbb0", "tbb1", "tbb2"]

def _sample_name_from_pdbfile(pdbfile):
    base = os.path.basename(str(pdbfile))
    return next((base[:-len(ext)] for ext in (".pdb.gz", ".pdb", ".cif.gz", ".cif", ".gz")
                 if base.endswith(ext)), os.path.splitext(base)[0])

def _describe_rmsd_threshold(rmsd_threshold):
    return 'per-combo normalize_rmsd(n_atoms, "cgvdmbb")' if rmsd_threshold is None else str(rmsd_threshold)

def write_results(outdir, all_matches, rmsd_records, rmsd_threshold):
    os.makedirs(outdir, exist_ok=True)
    if rmsd_records:
        path = os.path.join(outdir, "ligand_rmsd_summary.tsv")
        with open(path, "w") as fh:
            fh.write("pdbfile\tligand_rmsd\n")
            fh.writelines(f"{pdbfile}\t{rmsd:.3f}\n" for pdbfile, rmsd in sorted(rmsd_records))
        print(f"Wrote ligand RMSD summary to: {path}")

    path = os.path.join(outdir, "vdg_hits.tsv")
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_matches)
    print(f"Wrote vdG match results to:   {path}")
    if not all_matches:
        return

    summary = pd.crosstab(
        pd.Series([_sample_name_from_pdbfile(r["pdbfile"]) for r in all_matches], name="pdb"),
        pd.Series([r["frag"] for r in all_matches])).sort_index(axis=0).sort_index(axis=1)
    path = os.path.join(outdir, "results_summary.txt")
    with open(path, "w") as fh:
        fh.write("vdG hit counts per sample (rows) and fragment (columns); counts use "
                 f"rmsd_threshold={_describe_rmsd_threshold(rmsd_threshold)}\n\n")
        fh.write(pd.crosstab(
            pd.Series([_sample_name_from_pdbfile(r["pdbfile"]) for r in all_matches], name="pdb"),
            pd.Series([r["frag"] for r in all_matches])).sort_index(axis=0).sort_index(axis=1).to_string() + "\n")
    print(f"Wrote summary table to:        {path}")
