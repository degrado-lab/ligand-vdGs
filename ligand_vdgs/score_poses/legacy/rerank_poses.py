# rerank_poses.py
"""
Given a set of query protein–ligand complexes and a vdG library:

  - Fragment each ligand into CG-like fragments (SMILES).
  - For each frag + binding-site residue combo:
      - Load vdG centroids for matching AA bucket from npz.
      - Enumerate vdG AA permutations (incl. 'bb' wildcard).
      - Enumerate ligand CG permutations per site.
      - Align vdG [backbone + CG] to query [backbone + CG] with batched Kabsch.
      - Record the best RMSD hit per (frag, BSR combo, vdG centroid).

Outputs (per config.yml):
  outdir/<yml_basename>/:
    - ligand_rmsd_summary.tsv
    - vdg_match_results.tsv
    - rerank_log.txt
"""

import sys, os, io, csv, yaml, logging, tempfile, time
import multiprocessing as mp
import numpy as np
import pandas as pd
import prody as pr
from contextlib import redirect_stdout, redirect_stderr
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign as MA

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "functions"))
import dock_utils as dock      # binding site combos, CG coords
import Frags                   # ligand fragmentation
import vdg_npz_utils as vdg_npz
from utils import kabsch, convert_time_elapsed, smiles_equiv
# kabsch:
#   - X: [M, N, 3] or [N, 3]
#   - Y: [M, N, 3]
#   -> R [M,3,3], t [M,3], ssd [M]

with open(sys.argv[1], "r") as f:
    config = yaml.safe_load(f)

rmsd_threshold     = config["rmsd_threshold"]
query_pdbs_dir     = config["query_pdbs_dir"]
solved_struct_path = config["solved_struct_path"]
lig_smiles         = config["lig_smiles"]
vdg_lib_dir        = config["vdg_lib_dir"]
base_outdir        = config["outdir"]
frags_to_exclude   = config["frags_to_exclude"]
frags_to_include   = config["frags_to_include"]
query_pdbs         = config.get("query_pdbs", "all")
overwrite_existing = config.get("overwrite_existing", False)  # kept for config, not used here
yml_basename       = os.path.splitext(os.path.basename(sys.argv[1]))[0]
outdir             = os.path.join(base_outdir, yml_basename)

# Caches shared across worker processes (set inside workers)
_SOLVED_STRUCT   = None
_ALL_BSR_COMBOS  = None
_REF_LIG_MOL     = None
EXCELLENT_MATCH_CUTOFF = 0.4 # Angstroms

def _init_rdkit_logging(stream_like):
    """Route RDKit log messages to a local buffer so mp worker logs don't interleave."""
    rdBase.LogToPythonLogger()
    logger = logging.getLogger("rdkit")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(stream_like)
    sh.setFormatter(logging.Formatter("%(levelname)s [%(name)s]: %(message)s"))
    logger.addHandler(sh)

def check_frag_in_vdg_lib(sub_smiles, frags_in_lib, vdg_lib_dir):
    """
    Decide whether a fragment SMILES has a vdG library entry.
    Handles SMILES equivalence + yml exclude list.

    Returns:
        frags_in_lib : dict[db_name -> job_complete_bool]
        db_name      : SMILES string used as vdG lib directory name
    """
    # Prefer exact SMILES dir; fall back to equivalent SMILES if present
    found_db_name = False
    if sub_smiles in os.listdir(vdg_lib_dir):
        db_name, found_db_name = sub_smiles, True
    else:
        for existing in os.listdir(vdg_lib_dir):
            if smiles_equiv(existing, sub_smiles, check_atom_order=False):
                db_name, found_db_name = existing, True
                break
    if not found_db_name:
        db_name = sub_smiles

    # Exclude list from yml (use vdG lib name)
    if Frags.check_in_exclude_list(db_name, frags_to_exclude):
        return frags_in_lib, db_name

    # Check vdG lib job status (did clus+dedup finish?)
    if db_name not in frags_in_lib:
        if db_name in os.listdir(vdg_lib_dir):
            frags_in_lib[db_name] = Frags.check_vdg_job_status(db_name, vdg_lib_dir)
        else:
            frags_in_lib[db_name] = False
    return frags_in_lib, db_name

def should_run_frag(frag_name, frags_to_include, frags_in_lib):
    """
    Return True if this frag (vdG lib name) should be sampled,
    given include list + lib status.
    """
    if frags_to_include == "all":
        pass
    elif isinstance(frags_to_include, list):
        if frag_name not in frags_to_include:
            return False
    else:
        raise ValueError("frags_to_include must be 'all' or a list of SMILES.")

    if frag_name not in frags_in_lib or not frags_in_lib[frag_name]:
        return False
    return True

def _bsr_label_to_string(bsr_combo):
    """
    Human/TSV-friendly encoding of a binding-site residue combo:
    list[(seg, chain, resnum)] -> 'seg:chain:resnum;seg:chain:resnum;...'
    """
    return ";".join(f"{seg}:{chain}:{resnum}" for (seg, chain, resnum) in bsr_combo)

def _reorder_sub_to_target_smiles(sub_mol, target_smiles):
    """
    Reorder atoms in `sub_mol` so that their element order matches `target_smiles`.

    Both are assumed to be SMILES-equivalent (checked earlier via smiles_equiv).
    We use `target_smiles` as a SMARTS pattern and map it onto `sub_mol`.
    """
    pattern = Chem.MolFromSmarts(target_smiles)
    if pattern is None:
        raise ValueError(f"Could not parse target SMILES as SMARTS: {target_smiles}")

    # We expect a full-graph match: same atom count and connectivity
    matches = sub_mol.GetSubstructMatches(pattern, uniquify=False)
    if not matches:
        raise ValueError(
            f"Could not substructure-match fragment to target SMILES.\n"
            f"target_smiles={target_smiles}")

    perm = matches[0]
    if len(perm) != sub_mol.GetNumAtoms():
        raise ValueError(
            f"Permutation length {len(perm)} != num atoms {sub_mol.GetNumAtoms()} "
            f"for target_smiles={target_smiles}")

    # RenumberAtoms: new index i gets old atom perm[i]
    return Chem.RenumberAtoms(sub_mol, perm)

def _process_one_pdb(pdbfile):
    """
    Run vdG matching for one PDB; capture:
      - ligand RMSD to solved structure
      - vdG matches (best RMSD per vdG centroid / frag / BSR combo)
      - log text for this PDB
    """
    lig_rmsd_value, failed_buckets, match_records = None, [], []
    global _SOLVED_STRUCT, _ALL_BSR_COMBOS, _REF_LIG_MOL

    # Lazily initialize solved struct + BSR combos in the worker
    if _SOLVED_STRUCT is None or _ALL_BSR_COMBOS is None:
        _SOLVED_STRUCT = pr.parsePDB(solved_struct_path)
        lig_obj = _SOLVED_STRUCT.hetatm.select(
            "not (ion or resname SEP or resname TPO or resname MSE)")
        if lig_obj is None:
            raise ValueError(
                f"ProDy could not find ligand atoms in {solved_struct_path}. "
                "You may need to reformat the lig resname in the PDB file.")
        ligname = list(set(lig_obj.getResnames()))[0]
        _ALL_BSR_COMBOS = dock.get_bsr_combinations(_SOLVED_STRUCT, ligname)

    # Build reference RDKit ligand once (for RMSD vs target structure)
    if _REF_LIG_MOL is None:
        _solved = _SOLVED_STRUCT or pr.parsePDB(solved_struct_path)
        lig_obj = _solved.hetatm.select(
            "not (ion or resname SEP or resname TPO or resname MSE)")
        if lig_obj is None:
            raise ValueError(
                f"ProDy could not find ligand atoms in {solved_struct_path}. "
                "You may need to reformat the lig resname in the PDB file.")
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as _tmp_ref_pdb:
            pr.writePDB(_tmp_ref_pdb.name, lig_obj)
            _ref_pdb_mol = Chem.MolFromPDBFile(_tmp_ref_pdb.name, removeHs=True)
        _ref_template = Chem.MolFromSmiles(lig_smiles)
        if _ref_pdb_mol is None or _ref_template is None:
            raise ValueError("Failed to build RDKit molecules for target ligand RMSD "
                "(MolFromPDBFile / MolFromSmiles).")
        _ref_assigned = AllChem.AssignBondOrdersFromTemplate(_ref_template, _ref_pdb_mol)
        (_REF_LIG_MOL, _ref_perm_inds), _ref_smiles_no_H = Frags.manually_remove_Hs(
            _ref_assigned, return_single_mol_or_perms="single")

    # Capture all stdout/stderr for this pdbfile into a buffer
    buf = io.StringIO()
    with redirect_stdout(buf), redirect_stderr(buf):
        _init_rdkit_logging(buf)

        # Fragment ligand in this PDB
        frags_in_lib = {}
        pdb_path = os.path.join(query_pdbs_dir, pdbfile)
        orig_filtered_frags, pdb_mol_reassigned = Frags.get_frags_from_pdbfile(
            pdb_path, lig_smiles, quiet=True)

        # Symmetry-aware heavy-atom ligand RMSD vs solved structure
        if pdb_mol_reassigned is not None and _REF_LIG_MOL is not None:
            try:
                if pdb_mol_reassigned.GetNumAtoms() != _REF_LIG_MOL.GetNumAtoms():
                    print(
                        f'WARNING: ligand atom count mismatch; skipping RMSD for {pdbfile}. '
                        f'query={pdb_mol_reassigned.GetNumAtoms()} '
                        f'ref={_REF_LIG_MOL.GetNumAtoms()}')
                else:
                    lig_rmsd_value = MA.CalcRMS(pdb_mol_reassigned, _REF_LIG_MOL)
                    print(
                        f'Ligand→target RMSD (heavy-atom, symmetry-corrected): '
                        f'{lig_rmsd_value:.3f} Å')
            except Exception as e:
                print(f"WARNING: failed to compute ligand RMSD for {pdbfile}: {e}")

        # ------------------------------------------------------------------ #
        # Filter fragments by vdG lib presence / job status AND
        #   normalize everything to vdG lib name (db_name).
        #
        # Result:
        #   filtered_frags[frag_name] = grouped_by_site
        #     where:
        #       - frag_name == db_name (the vdG lib SMILES)
        #       - each `sub` mol has been reordered to match `frag_name`
        # ------------------------------------------------------------------ #
        filtered_frags = {}  # key: frag_name (db_name), value: grouped_by_site
        for sub_smiles, substruct_perm_grouped_by_site in orig_filtered_frags.items():
            frags_in_lib, db_name = check_frag_in_vdg_lib(
                sub_smiles, frags_in_lib, vdg_lib_dir)
            if not frags_in_lib.get(db_name, False):
                continue

            # If the ligand fragment SMILES == vdG lib SMILES, we can reuse as-is.
            # Otherwise, reorder each sub-mol to match db_name atom order.
            if db_name == sub_smiles:
                new_grouped = substruct_perm_grouped_by_site
            else:
                new_grouped = []
                for substruct_site in substruct_perm_grouped_by_site:
                    new_site = []
                    for sub, perm_inds, orig_mol_inds in substruct_site:
                        try:
                            sub_reordered = _reorder_sub_to_target_smiles(sub, db_name)
                        except Exception as e:
                            print(f"WARNING: failed to reorder fragment {sub_smiles} "
                                  f"to vdG name {db_name}: {e}")
                            # Skip this instance; others may still be usable
                            continue
                        new_site.append((sub_reordered, perm_inds, orig_mol_inds))
                    if new_site:
                        new_grouped.append(new_site)

            if not new_grouped:
                continue

            # Merge if multiple ligand SMILES map to the same vdG name
            if db_name in filtered_frags:
                filtered_frags[db_name].extend(new_grouped)
            else:
                filtered_frags[db_name] = new_grouped

        # For each fragment that exists in vdG lib (using vdG names):
        for frag_name, substruct_perm_grouped_by_site in filtered_frags.items():
            if not should_run_frag(frag_name, frags_to_include, frags_in_lib):
                print(f"Skipping {frag_name}..........................")
                continue
            print(f"Searching for {frag_name}..........................")

            # Precompute query CG permutations grouped by site on the ligand.
            # Each site: list of (q_cg_coords, perm_inds, orig_mol_inds)
            grouped_q_cg_perms = []
            for substruct_site in substruct_perm_grouped_by_site:
                site_perms = []
                for substruct_perm in substruct_site:
                    sub, perm_inds, orig_mol_inds = substruct_perm
                    # IMPORTANT: frag_name is now the vdG lib SMILES.
                    q_cg_coords = dock.get_query_cg_coords(sub, frag_name)
                    site_perms.append(
                        (np.asarray(q_cg_coords, np.float32), perm_inds, orig_mol_inds))
                grouped_q_cg_perms.append(site_perms)

            # Iterate over all binding-site residue combinations
            for combo, bsr_combo, bsrAAs, coords in _ALL_BSR_COMBOS:
                # combo: tuple of AA identities incl. 'bb' variants
                # coords: [[N,CA,C] per residue]
                bsr_incl_bb_identities, input_bsr_bb_coords = combo, coords

                # Sort BSR AAs + backbone coords into a deterministic canonical order
                paired = list(zip(bsr_incl_bb_identities, input_bsr_bb_coords))
                paired_sorted = sorted(paired, key=lambda pair: pair[0])
                bsr_incl_bb_identities, input_bsr_bb_coords = zip(*paired_sorted)
                bsr_incl_bb_identities = list(bsr_incl_bb_identities)
                input_bsr_bb_coords = np.asarray(input_bsr_bb_coords, dtype=np.float32
                    ).reshape(-1, 3)

                subset_size = len(bsr_incl_bb_identities)
                aa_bucket = vdg_npz.make_aa_bucket(bsr_incl_bb_identities, sort=True)

                # Precompute all query bb+CG permutations for this (frag, BSR combo).
                Y_list, meta = [], []
                for site_idx, q_cg_site in enumerate(grouped_q_cg_perms):
                    for cg_perm_idx, (q_cg_coords, perm_inds, orig_mol_inds) in enumerate(
                        q_cg_site):
                        q_bb_and_cg = np.concatenate((input_bsr_bb_coords, q_cg_coords), 
                            axis=0).astype(np.float32)
                        Y_list.append(q_bb_and_cg)
                        meta.append((site_idx, cg_perm_idx))
                if not Y_list:
                    continue
                try:
                    Y = np.stack(Y_list, axis=0)  # [num_perms, n_atoms, 3]
                except ValueError as e:
                    print(f"WARNING: skipping frag {frag_name}, BSR combo "
                          f"{_bsr_label_to_string(bsr_combo)} due to CG shape mismatch: {e}")
                    continue
                n_atoms = Y.shape[1]

                # Load vdG npz bucket for this fragment / subset_size / AA bucket
                bucket = vdg_npz.load_vdg_bucket(vdg_lib_dir, frag_name, subset_size, aa_bucket, use_cache=True)
                if bucket is None:
                    continue

                # Iterate over vdG centroids in this bucket
                for vdg_centroid in vdg_npz.iter_bucket_centroids(bucket):
                    vdg_AAs = vdg_centroid["resnames"]          # per vdM, incl 'bb'
                    vdg_bb  = vdg_centroid["bb"]                # [num_vdms, 3, 3]
                    vdg_cg  = vdg_centroid["cg"]                # [n_cg, 3]
                    vdg_idx = vdg_centroid["index"]
                    clus_id = vdg_centroid["cluster_id"]

                    # Map query AA identities (incl 'bb' wildcard) onto vdG vdM positions.
                    resind_perms = dock.map_aa_identities_to_vdg_resinds(vdg_AAs, 
                        bsr_incl_bb_identities, wildcard="bb")
                    if not resind_perms:
                        continue

                    best_rmsd = None
                    best_aa_perm_idx = None
                    best_q_site_idx = None
                    best_q_cg_perm_idx = None

                    # For each AA-permutation, build vdG bb+CG coords and align
                    for aa_perm_idx, resind_perm in enumerate(resind_perms):
                        # Flatten vdG backbone in permuted order: [subset_size*3, 3]
                        vdg_perm_bb_coords = vdg_npz.permute_backbone_flat(vdg_bb, resind_perm)
                        try:
                            db_bb_and_cg = np.concatenate((vdg_perm_bb_coords, vdg_cg), 
                                axis=0).astype(np.float32)
                        except Exception: # shape mismatch; skip this perm
                            continue

                        if db_bb_and_cg.shape[0] != n_atoms:
                            # Inconsistent vdG CG size vs query CG size; skip
                            continue

                        # Batched Kabsch on backbone+CG vs all query perms
                        _, _, ssd_batch = kabsch(db_bb_and_cg, Y)  # X is [n_atoms, 3]
                        rmsd_batch = np.sqrt(ssd_batch / float(n_atoms))

                        idx_ok = np.where(rmsd_batch <= rmsd_threshold)[0]
                        if idx_ok.size:
                            i_loc = int(idx_ok[np.argmin(rmsd_batch[idx_ok])])
                            rmsd_loc = float(rmsd_batch[i_loc])
                            if best_rmsd is None or rmsd_loc < best_rmsd:
                                best_rmsd = rmsd_loc
                                best_aa_perm_idx = aa_perm_idx
                                best_q_site_idx, best_q_cg_perm_idx = meta[i_loc]
                                if best_rmsd <= EXCELLENT_MATCH_CUTOFF:
                                    break

                    # If we found any acceptable RMSD for this vdG centroid, record it
                    if best_rmsd is not None:
                        sample    = os.path.basename(pdbfile).removesuffix(".pdb")
                        struct_id = sample[:4]
                        match_records.append(
                            dict(
                                pdbfile=pdbfile,
                                sample=sample,
                                struct_id=struct_id,
                                frag=frag_name,  # vdG lib SMILES
                                subset_size=subset_size,
                                bsr_combo=_bsr_label_to_string(bsr_combo),
                                aa_bucket=aa_bucket,
                                vdg_index=vdg_idx,
                                vdg_cluster_id=clus_id,
                                vdg_rmsd=f"{float(best_rmsd):.2f}",
                                aa_perm_idx=int(best_aa_perm_idx)
                                if best_aa_perm_idx is not None
                                else "",
                                q_site_idx=int(best_q_site_idx)
                                if best_q_site_idx is not None
                                else "",
                                q_cg_perm_idx=int(best_q_cg_perm_idx)
                                if best_q_cg_perm_idx is not None
                                else "",))

    log_text = buf.getvalue()
    return (pdbfile, log_text, failed_buckets, lig_rmsd_value, match_records)

def _write_results_and_summary(all_matches, rmsd_records, outdir, print_score_table=True):
    os.makedirs(outdir, exist_ok=True)

    # (1) ligand RMSD summary
    rmsd_path = os.path.join(outdir, "ligand_rmsd_summary.tsv")
    with open(rmsd_path, "w") as fh:
        fh.write("pdbfile\tligand_rmsd\n")
        for pdbfile, rmsd in sorted(rmsd_records):
            fh.write(f"{pdbfile}\t{rmsd:.3f}\n")

    # (2) per-match vdG results
    results_path = os.path.join(outdir, "vdg_match_results.tsv")
    fieldnames = [
        "pdbfile",
        "sample",
        "struct_id",
        "frag",
        "subset_size",
        "bsr_combo",
        "aa_bucket",       # canonical AA bucket label (e.g. "ALA_bb")
        "vdg_index",       # centroid index within npz file for this bucket
        "vdg_cluster_id",  # cluster_id from clustering
        "vdg_rmsd",        # best bb+CG RMSD for this vdG vs query
        "aa_perm_idx",     # which AA→vdM permutation was best
        "q_site_idx",      # which ligand site (if multiple instances)
        "q_cg_perm_idx",   # which CG atom permutation within that site
    ]
    with open(results_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for rec in all_matches:
            w.writerow(rec)

    print(f"\nWrote ligand RMSD summary to: {rmsd_path}")
    print(f"Wrote vdG match results to:   {results_path}")

    if not print_score_table or not all_matches:
        return

    # (3) simple hit-count summary: sample×frag counts; sorted by total hits per sample
    scores = {}
    for rec in all_matches:
        s, f = rec["sample"], rec["frag"]
        scores.setdefault(s, {}).setdefault(f, 0.0)
        scores[s][f] += 1.0

    df = pd.DataFrame(scores).fillna(0).T
    if df.empty:
        return

    # Sort rows by total hits (descending); keep an "Average" row at the bottom.
    sample_totals = df.sum(axis=1)
    df = df.loc[sample_totals.sort_values(ascending=False).index]
    avg_row = df.mean().round(2)
    df.loc["Average"] = avg_row
    df = pd.concat([df.drop("Average"), df.loc[["Average"]]])

    print(f"\n--- vdG hit counts per sample (rows) and fragment (cols); "
          f"counts use rmsd_threshold={rmsd_threshold} ---")
    header = " " * 12 + "  ".join(f"{col:7}" for col in df.columns)
    print(header)
    print("-" * len(header))
    for idx, row in df.iterrows():
        if idx == "Average":
            formatted = "  ".join(f"{v:7.2f}" for v in row)
        else:
            formatted = "  ".join(f"{int(v):7d}" for v in row)
        print(f"{idx:12} {formatted}")

def main():
    print("Config variables:")
    for k in ["rmsd_threshold", "query_pdbs_dir", "solved_struct_path", "lig_smiles",
        "vdg_lib_dir", "base_outdir", "outdir", "frags_to_exclude", "frags_to_include",
        "query_pdbs", "overwrite_existing", "yml_basename",]:
        print(f"{k}: {globals()[k]}")

    # Precompute binding-site residues + BSR combos once in main
    solved_struct = pr.parsePDB(solved_struct_path)
    lig_obj = solved_struct.hetatm.select(
        "not (ion or resname SEP or resname TPO or resname MSE)")
    if lig_obj is None:
        raise ValueError(f"ProDy could not find ligand atoms in {solved_struct_path}. "
                          "You may need to reformat the lig resname in the PDB file.")
    ligname = list(set(lig_obj.getResnames()))
    assert len(ligname) == 1
    ligname = ligname[0]

    all_bsr_combos = dock.get_bsr_combinations(solved_struct, ligname)
    global _SOLVED_STRUCT, _ALL_BSR_COMBOS
    _SOLVED_STRUCT, _ALL_BSR_COMBOS = solved_struct, all_bsr_combos

    # Decide which PDBs to run
    all_pdbs = sorted(os.listdir(query_pdbs_dir))
    if query_pdbs == "all":
        pdb_list = all_pdbs
    elif isinstance(query_pdbs, list):
        pdb_list = [p for p in all_pdbs if p in query_pdbs]
    else:
        raise ValueError("`query_pdbs` must be 'all' or a list of pdb files.")
    if not pdb_list:
        print("No pdb files selected to run.")
        return
    os.makedirs(outdir, exist_ok=True)

    # Print per-fragment summary once (using first PDB) to show which frags are in vdG lib.
    first_pdb = pdb_list[0]
    prelog_buf = io.StringIO()
    with redirect_stdout(prelog_buf), redirect_stderr(prelog_buf):
        first_pdb_path = os.path.join(query_pdbs_dir, first_pdb)
        orig_filtered_frags_once, _pdbmol = Frags.get_frags_from_pdbfile(
            first_pdb_path, lig_smiles, quiet=False)
        frags_in_lib_once, dbname_keys_once = {}, {}
        for sub_smiles in orig_filtered_frags_once.keys():
            frags_in_lib_once, db_name = check_frag_in_vdg_lib(
                sub_smiles, frags_in_lib_once, vdg_lib_dir)
            dbname_keys_once[sub_smiles] = db_name
        Frags.summarize_frags(frags_in_lib_once, frags_to_exclude, frags_to_include)

    print("Fragment list and summary:")
    header_text = prelog_buf.getvalue()
    if header_text:
        print(header_text, end="" if header_text.endswith("\n") else "\n",flush=True)

    # Multiprocessing: per-PDB worker
    n_workers = min(len(pdb_list), max(1, mp.cpu_count() - 1))
    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(_process_one_pdb, pdb_list)

    # Collect ligand RMSDs + matches
    rmsd_records, all_matches = [], []
    for pdbfile, log_text, _failed, lig_rmsd_value, match_records in results:
        sep = f"\n==================== {pdbfile} ====================\n"
        sys.stdout.write(sep)
        if log_text:
            sys.stdout.write(log_text)
            if not log_text.endswith("\n"):
                sys.stdout.write("\n")
        sys.stdout.flush()
        if lig_rmsd_value is not None:
            rmsd_records.append((pdbfile, lig_rmsd_value))
        all_matches.extend(match_records)

    # Collect any bucket-level failures (currently unused, but left for symmetry)
    all_failed = set()
    for _pdbfile, _log_text, _failed, _lig_rmsd, _matches in results:
        all_failed.update(_failed)
    if all_failed:
        print("\nvdG buckets that could not be handled (unique across all jobs):")
        for f in sorted(all_failed):
            print(f" - {f}")
        print("", flush=True)

    _write_results_and_summary(all_matches, rmsd_records, outdir, print_score_table=True)

if __name__ == "__main__":
    os.makedirs(outdir, exist_ok=True)
    log_path = os.path.join(outdir, "rerank_log.txt")
    start_time = time.perf_counter()
    with open(log_path, "w") as _log_fh, redirect_stdout(_log_fh), redirect_stderr(_log_fh):
        main()
        elapsed = time.perf_counter() - start_time
        h, m, s = convert_time_elapsed(elapsed)
        print(f"\nTotal job time: {h} h, {m} mins, and {round(s)} secs.")
