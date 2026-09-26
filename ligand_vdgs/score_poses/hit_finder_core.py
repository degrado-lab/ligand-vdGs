"""
Core vdG hit-finding: score a query structure's ligand against the vdG fragment
library. The RDKit mol and bucket caches
are per process (eager via init_worker in pool workers, else lazy), so scoring
many models in one process reuses the library arrays.
"""

import io
import multiprocessing as mp
import os
import tempfile
import traceback
from contextlib import redirect_stdout, redirect_stderr

import numpy as np
import prody as pr
from ligand_vdgs.functions import ligand_structure
from ligand_vdgs.functions import dock_utils as dock
from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions.vdg_struct_utils import bsr_contact_atoms, vcb_slot_blockers
from ligand_vdgs.functions.utils import kabsch, kabsch_ssd, best_inplace_symmetry_rmsd, normalize_rmsd
from ligand_vdgs.functions.vdg_fp_utils import (fp_tolerances,
    prefilter_query_indices_pair, prefilter_query_indices_single)
from ligand_vdgs.score_poses import hit_finder_matching as _matching
from ligand_vdgs.score_poses.hit_finder_matching import (
    BUCKET_CACHE_MAX_BYTES, CACHE_MISS, _BoundedBucketCache, _cg_symmetry_for,
    _combo_worker_struct,
    _struct_id_from_pdbfile, _worker_bucket_cache, bsr_label_to_string,
    cg_atom_order_smarts, init_rdkit_logging, init_worker, lib_entries,
    match_library_frags_to_query, query_resonance_forms)
from ligand_vdgs.score_poses.hit_finder_geometry import (
    ACCEPTOR_ELEMENTS, BB_SLOT_SIDECHAIN_CLASH, PLACEMENT_CHUNK as _PLACEMENT_CHUNK,
    PRO_NH_DONOR_CUTOFF, _dist, _min_dist2,
    backbone_slots_can_host, backbone_slots_host_mask, cg_bb_dist, cg_center,
    fp_pair_from_ca_and_cgcom, fp_single_from_ca_and_cgcom,
    has_any_bsr_atom_cg_contact, held_out_placement,
    proper_rotation_mask as _proper_rotation_mask)

EXCELLENT_MATCH_CUTOFF = 0.3
MATCH_MODES = ("joint", "bb")

def _hit_record(common, bucket, nr_idx, match_mode, rmsd, aa_perm_idx, site_idx, perm_idx,
                grouped_q_cg_perms, R, t, R_bb, t_bb, bb_rmsd, held_out, lever):
    rec = dict(common, charge_sign=str(bucket["charge_signs"][nr_idx]),
               vdg_index=int(bucket["partition_indices"][nr_idx]),
               vdg_cluster_id=int(bucket["cluster_id"][nr_idx]),
               vdg_cluster_num_parents=int(bucket["cluster_num_parents"][nr_idx]),
               vdg_rmsd=f"{rmsd:.4f}", aa_perm_idx=int(aa_perm_idx), q_site_idx=int(site_idx), q_cg_perm_idx=int(perm_idx),
               q_atom_indices=";".join(str(i) for i in sorted(
                   grouped_q_cg_perms[site_idx][perm_idx][-1])),
               match_mode=match_mode, bb_rmsd=f"{bb_rmsd:.4f}",
               held_out_cg_rmsd=f"{held_out:.4f}", cg_bb_dist=f"{lever:.4f}")
    for prefix, M, v in (("", R, t), ("bb", R_bb, t_bb)):
        M, v = np.round(M, 4), np.round(v, 4)
        rec.update({f"R{prefix}{i}{j}": f"{M[i, j]:.4f}" for i in range(3) for j in range(3)})
        rec.update({f"t{prefix}{k}": f"{v[k]:.4f}" for k in range(3)})
    return rec

def _compact_record(common, bucket, match_mode, sel, vdg_rmsd, bb_rmsd, held_out, lever, placed, q_lig_idx,
                    **extra):
    """One record for all hits of a bucket; per-hit fields are arrays (hits.npz contract)."""
    return dict(common, match_mode=match_mode, compact=True, vdg_index=bucket["partition_indices"][sel],
                charge_sign=bucket["charge_signs"][sel], vdg_cluster_id=bucket["cluster_id"][sel],
                vdg_cluster_num_parents=bucket["cluster_num_parents"][sel].astype(np.int32),
                nr_parent_biounit=bucket["parent_biounit"][sel], vdg_rmsd=np.asarray(vdg_rmsd, np.float32),
                bb_rmsd=np.asarray(bb_rmsd, np.float32), held_out_cg_rmsd=np.asarray(held_out, np.float32),
                cg_bb_dist=np.asarray(lever, np.float32), placed_cg=np.asarray(placed, np.float32),
                q_lig_atom_idx=np.asarray(q_lig_idx, np.int32), placed_cg_element=bucket["cg_elements"][sel],
                **extra)

def _score_one_bsr_combo(pdbfile, struct, frag_name, query_frag, grouped_q_cg_perms,
                          combo_item, vdg_lib_dir, bucket_cache, bb_blockers_cache,
                           rmsd_threshold, contact_cutoff,
                          lig_instance_label, match_mode="joint", compact=False):
    """Hits of one fragment on one BSR combo. No mode reads side chains. match_mode='joint' matches
    on bb+CG RMSD <= tau (slot gate and contact filter on the query CG vs backbone + virtual CB);
    'bb' (placement, DR-36) matches on the backbone alone, sqrt(SSD_bb / n_atoms) <= tau, reads
    neither the query CG nor side chains (slot gate = virtual CB + Pro N on the placed CG), and
    reports the held-out CG error of the bb-only placement.
    compact: one record per bucket whose per-hit fields are arrays (`_compact_record`)."""
    if match_mode not in MATCH_MODES:
        raise ValueError(f"[ERROR] match_mode must be one of {MATCH_MODES}, got {match_mode!r}")
    if match_mode == "bb" and contact_cutoff is not None:
        raise ValueError("[ERROR] contact_cutoff reads the query CG and side chains; it is joint-mode only.")
    combo, bsr_combo, _bsr_AAs, coords = combo_item
    match_records = []
    (bsr_incl_bb_identities, input_bsr_bb_coords, bsr_combo,
     bsr_resnames) = zip(*sorted(zip(combo, coords, bsr_combo, _bsr_AAs), key=lambda x: x[0]))
    bsr_incl_bb_identities = list(bsr_incl_bb_identities)
    bsr_combo = list(bsr_combo)
    input_bsr_bb_coords = np.asarray(input_bsr_bb_coords, np.float32)
    subset_size = input_bsr_bb_coords.shape[0]
    bb_flat = input_bsr_bb_coords.reshape(-1, 3)
    N_bb = bb_flat.shape[0]
    input_bsr_ca_coords = input_bsr_bb_coords[:, 1, :]
    aa_bucket = vdg_npz.make_aa_bucket(bsr_incl_bb_identities)

    total_q_perms = sum(len(site) for site in grouped_q_cg_perms)
    if total_q_perms == 0:
        return match_records

    N_cg = grouped_q_cg_perms[0][0][0].shape[0]
    n_atoms = N_bb + N_cg
    effective_rmsd = (normalize_rmsd(n_atoms, "cgvdmbb")
                      if rmsd_threshold is None else rmsd_threshold)
    fp_tol = fp_tolerances(effective_rmsd, n_atoms, N_cg, subset_size)

    # Joint-mode contact filter: backbone + virtual CB only (no side chains in any mode).
    bsr_atom_coords = None if contact_cutoff is None else bsr_contact_atoms(input_bsr_bb_coords, bsr_resnames)
    blocker_key = (tuple(bsr_combo), tuple(bsr_incl_bb_identities))
    if blocker_key in bb_blockers_cache:
        bb_blockers = bb_blockers_cache[blocker_key]
    else:
        bb_blockers = vcb_slot_blockers(input_bsr_bb_coords, bsr_incl_bb_identities, bsr_resnames)
        bb_blockers_cache[blocker_key] = bb_blockers

    Y = np.empty((total_q_perms, n_atoms, 3), np.float32)
    meta = np.empty((total_q_perms, 2), dtype=np.int32)
    query_contact_ok = np.empty(total_q_perms, dtype=bool)
    query_bb_slot_ok = np.empty(total_q_perms, dtype=bool)

    if subset_size == 1:
        fp_q0 = np.empty(total_q_perms, np.float32)
        fp_q1 = fp_q2 = None
    else:
        fp_q0 = np.empty(total_q_perms, np.float32)
        fp_q1 = np.empty(total_q_perms, np.float32)
        fp_q2 = np.empty(total_q_perms, np.float32)

    idx = 0
    consistent = True
    for site_idx, q_cg_site in enumerate(grouped_q_cg_perms):
        for cg_perm_idx, (q_cg_coords, q_cg_com, q_acceptor,
                          perm_inds, orig_mol_inds) in enumerate(q_cg_site):
            if q_cg_coords.shape[0] != N_cg:
                consistent = False
                break

            Y[idx, :N_bb] = bb_flat
            Y[idx, N_bb:] = q_cg_coords
            meta[idx, 0] = site_idx
            meta[idx, 1] = cg_perm_idx

            if contact_cutoff is not None:
                query_contact_ok[idx] = has_any_bsr_atom_cg_contact(
                    bsr_atom_coords, q_cg_coords, cutoff=contact_cutoff)

            if match_mode == "joint":
                query_bb_slot_ok[idx] = backbone_slots_can_host(bb_blockers, q_cg_coords, q_acceptor)

            if subset_size == 1:
                fp_q0[idx] = fp_single_from_ca_and_cgcom(
                    input_bsr_ca_coords, q_cg_com)
            else:
                d1, d2, d12 = fp_pair_from_ca_and_cgcom(
                    input_bsr_ca_coords, q_cg_com)
                fp_q0[idx] = d1
                fp_q1[idx] = d2
                fp_q2[idx] = d12
            idx += 1

        if not consistent:
            break

    if not consistent or idx == 0:
        return match_records

    if contact_cutoff is not None and not np.any(query_contact_ok):
        return match_records

    if match_mode == "bb":  # placement: no query-CG gate; the slot gate acts on the placed CG
        query_bb_slot_ok[:] = True
    elif not np.any(query_bb_slot_ok):
        return match_records

    bucket_key = (vdg_lib_dir, frag_name, subset_size, aa_bucket)
    bucket = bucket_cache.get(bucket_key)
    if bucket is CACHE_MISS:
        bucket = vdg_npz.load_vdg_bucket_all_signs(vdg_lib_dir, frag_name, subset_size, aa_bucket)
        bucket_cache.put(bucket_key, bucket)
    if bucket is None:
        return match_records

    bucket_cg      = bucket["cg"]
    bucket_bb      = bucket["bb"]
    n_nr_vdgs    = bucket_cg.shape[0]

    bucket_parts = bucket["aa_bucket_parts"]
    if list(bucket_parts) != list(bsr_incl_bb_identities):
        print(f"[WARNING] ({pdbfile}) Bucket {aa_bucket} stores slot "
              f"labels {list(bucket_parts)} but was looked up for "
              f"{list(bsr_incl_bb_identities)}; skipping.", flush=True)
        return match_records
    resind_perms = vdg_npz.aa_perm_indices(bucket_parts)

    n_perms = len(resind_perms)
    perms_arr = np.asarray(resind_perms, dtype=np.intp)
    bb_ssd = kabsch_ssd(bb_flat, bucket_bb[:, perms_arr, :, :].reshape(
        n_nr_vdgs * n_perms, -1, 3)).reshape(n_nr_vdgs, n_perms)
    # Both modes match on sqrt(SSD_bb / n_atoms) <= tau, i.e. bb_rmsd = sqrt(SSD_bb / N_bb)
    # <= tau * sqrt(n_atoms / N_bb): every joint hit is also a placement hit.
    bb_rmsd_all = np.sqrt(bb_ssd / n_atoms)
    q_ok = query_bb_slot_ok & query_contact_ok if contact_cutoff is not None else query_bb_slot_ok
    if not q_ok.any():
        return match_records
    q_cgs, q_meta = Y[q_ok, N_bb:], meta[q_ok]
    q_lever = [cg_bb_dist(q, bb_flat) for q in q_cgs]
    # Query-ligand atom index per CG atom of each labeling (placed atom j <-> q_cgs[k][j]).
    q_lig_idx = np.array([grouped_q_cg_perms[s][p][3] for s, p in q_meta], np.int32)
    common = dict(pdbfile=pdbfile, struct_id=_struct_id_from_pdbfile(pdbfile),
                  lig_instance=lig_instance_label, frag=frag_name, query_frag=query_frag,
                  subset_size=subset_size, bsr_combo=bsr_label_to_string(bsr_combo),
                  aa_bucket=aa_bucket, rmsd_threshold=f"{effective_rmsd:.4f}")

    if match_mode == "bb":
        if bucket_cg.shape[1] != N_cg:
            return match_records
        matched = bb_rmsd_all <= effective_rmsd
        nr_sel = np.flatnonzero(matched.any(axis=1))
        if nr_sel.size == 0:
            return match_records
        best_p = np.where(matched, bb_rmsd_all, np.inf)[nr_sel].argmin(axis=1)
        vdg_bb_sel = bucket_bb[nr_sel[:, None], perms_arr[best_p]].reshape(nr_sel.size, N_bb, 3)
        vdg_cg_sel = bucket_cg[nr_sel]
        R_bb, t_bb, held_out, k, placed = held_out_placement(vdg_bb_sel, vdg_cg_sel, bb_flat, q_cgs)
        # Joint fit at the placement's correspondence, for vdg_rmsd and R/t.
        R, t, joint_ssd = kabsch(np.concatenate((vdg_bb_sel, vdg_cg_sel), axis=1),
                                 np.concatenate((np.broadcast_to(bb_flat, vdg_bb_sel.shape),
                                                 q_cgs[k]), axis=1))
        ok = _proper_rotation_mask(R) & _proper_rotation_mask(R_bb)
        if not ok.all():
            print(f"[WARNING] ({pdbfile}) {int((~ok).sum())} non-unitary rotation(s) in bucket "
                  f"{aa_bucket}; discarding those hits.", flush=True)
        ok &= backbone_slots_host_mask(bb_blockers, placed, np.isin(bucket["cg_elements"][nr_sel], ACCEPTOR_ELEMENTS))
        if compact:
            sel = nr_sel[ok]
            match_records.append(_compact_record(
                common, bucket, "bb", sel, np.sqrt(joint_ssd[ok] / n_atoms),
                np.sqrt(bb_ssd[sel, best_p[ok]] / N_bb), held_out[ok],
                np.asarray(q_lever, np.float32)[k[ok]], placed[ok], q_lig_idx[k[ok]]))
            return match_records
        match_records.extend(
            _hit_record(common, bucket, nr_idx, "bb", np.sqrt(joint_ssd[i] / n_atoms), best_p[i],
                        *q_meta[k[i]], grouped_q_cg_perms, R[i], t[i], R_bb[i], t_bb[i],
                        np.sqrt(bb_ssd[nr_idx, best_p[i]] / N_bb), held_out[i], q_lever[k[i]])
            for i, nr_idx in enumerate(nr_sel) if ok[i])
        return match_records

    joint_hits = []  # compact: (nr_idx, vdg_rmsd, bb_rmsd, held_out, lever, placed, q_lig_idx, in-sample cg_rmsd)
    for nr_idx in range(n_nr_vdgs):
        vdg_cg      = bucket_cg[nr_idx]
        vdg_bb      = bucket_bb[nr_idx]
        vdg_idx     = int(nr_idx)

        if vdg_cg.shape[0] != N_cg:
            continue

        vdg_cg_com = cg_center(vdg_cg)

        best_rmsd = best_aa_perm_idx = best_q_site_idx = None
        best_query = best_R = best_t = None

        for aa_perm_idx, resind_perm in enumerate(resind_perms):
            res_subset_bb = vdg_bb[np.asarray(resind_perm, dtype=np.intp)]
            res_subset_ca = res_subset_bb[:, 1, :]

            if subset_size == 1:
                fp_v0 = fp_single_from_ca_and_cgcom(res_subset_ca, vdg_cg_com)
                idxs = prefilter_query_indices_single(
                    fp_q0, fp_v0, fp_tol[0])
            else:
                fp_v0, fp_v1, fp_v2 = fp_pair_from_ca_and_cgcom(
                    res_subset_ca, vdg_cg_com)
                idxs = prefilter_query_indices_pair(
                    fp_q0, fp_q1, fp_q2,
                    fp_v0, fp_v1, fp_v2,
                    fp_tol,)

            if idxs is None or idxs.size == 0:
                continue

            if contact_cutoff is not None:
                idxs = idxs[query_contact_ok[idxs]]
                if idxs.size == 0:
                    continue

            idxs = idxs[query_bb_slot_ok[idxs]]
            if idxs.size == 0:
                continue

            if bb_rmsd_all[nr_idx, aa_perm_idx] > effective_rmsd:
                continue

            Y_sub = Y[idxs]

            db_bb_and_cg = np.concatenate(
                (res_subset_bb.reshape(-1, 3), vdg_cg), axis=0).astype(np.float32,
                                                     copy=False)
            if db_bb_and_cg.shape[0] != n_atoms:
                continue

            rmsd_batch = np.sqrt(
                kabsch_ssd(db_bb_and_cg, Y_sub) / float(n_atoms))

            idx_ok = np.where(rmsd_batch <= effective_rmsd)[0]
            if idx_ok.size == 0:
                continue

            local_idx = int(idx_ok[np.argmin(rmsd_batch[idx_ok])])
            rmsd_loc  = float(rmsd_batch[local_idx])

            if best_rmsd is None or rmsd_loc < best_rmsd:
                R_one, t_one, _ = kabsch(
                    db_bb_and_cg, Y_sub[local_idx][None, ...])
                R_cand = R_one[0]
                det_err = abs(float(np.linalg.det(R_cand)) - 1.0)
                ortho_err = float(np.max(np.abs(R_cand.T @ R_cand - np.eye(3))))
                if det_err > 1e-5 or ortho_err > 1e-4:
                    print(f"[WARNING] ({pdbfile}) Non-unitary rotation "
                          f"matrix for vdG {vdg_idx} bucket {aa_bucket} "
                          f"(det_err={det_err:.2e}, "
                          f"ortho_err={ortho_err:.2e}); discarding hit.",
                          flush=True)
                else:
                    best_rmsd          = rmsd_loc
                    best_aa_perm_idx   = aa_perm_idx
                    best_q_site_idx    = int(meta[idxs[local_idx], 0])
                    # Keep the labeling paired with the best fit; later candidates may lose.
                    best_query         = (int(meta[idxs[local_idx], 1]), Y_sub[local_idx, N_bb:])
                    best_R             = R_cand
                    best_t             = t_one[0]

                    if best_rmsd <= EXCELLENT_MATCH_CUTOFF:
                        break

        if best_rmsd is not None:
            # Joint rows: held out over the automorphisms of the joint fit's own site only.
            same_site = np.flatnonzero(q_meta[:, 0] == best_q_site_idx)
            R_bb, t_bb, held_out, k, placed = held_out_placement(
                vdg_bb[perms_arr[best_aa_perm_idx]].reshape(1, N_bb, 3), vdg_cg[None],
                bb_flat, q_cgs[same_site])
            k = same_site[k]
            if compact:
                joint_hits.append((nr_idx, best_rmsd, np.sqrt(bb_ssd[nr_idx, best_aa_perm_idx] / N_bb),
                                   held_out[0], q_lever[k[0]], placed[0], q_lig_idx[k[0]],
                                   np.sqrt(np.mean(np.sum((vdg_cg @ best_R + best_t - best_query[1]) ** 2, axis=1)))))
                continue
            match_records.append(_hit_record(
                common, bucket, nr_idx, "joint", best_rmsd, best_aa_perm_idx, best_q_site_idx,
                best_query[0], grouped_q_cg_perms, best_R, best_t, R_bb[0], t_bb[0],
                np.sqrt(bb_ssd[nr_idx, best_aa_perm_idx] / N_bb), held_out[0],
                q_lever[k[0]]))

    if compact and joint_hits:
        cols = list(zip(*joint_hits))
        match_records.append(_compact_record(common, bucket, "joint", np.asarray(cols[0], np.intp), *cols[1:7],
                                             in_sample_cg_rmsd=np.asarray(cols[7], np.float32)))
    return match_records

def _combo_worker(task):
    (pdbfile, pdb_path, frag_name, query_frag, grouped_q_cg_perms, combo_item,
     vdg_lib_dir, rmsd_threshold, contact_cutoff, lig_instance_label, match_mode, compact) = task
    local_buf = io.StringIO()
    with redirect_stdout(local_buf), redirect_stderr(local_buf), init_rdkit_logging(local_buf):
        try:
            # Scoring writes into the buffer; read it only after scoring returns.
            records = _score_one_bsr_combo(
                pdbfile, _combo_worker_struct(pdb_path), frag_name, query_frag, grouped_q_cg_perms,
                combo_item, vdg_lib_dir, _worker_bucket_cache(), {}, rmsd_threshold,
                contact_cutoff, lig_instance_label, match_mode, compact)
            return local_buf.getvalue(), records, None
        except Exception:
            return local_buf.getvalue(), [], traceback.format_exc()

def _raise_if_combo_task_errors(errors, n_tasks):
    if errors:
        raise RuntimeError(
            f"{len(errors)}/{n_tasks} combo tasks failed; first traceback:\n{errors[0]}")

def _match_record_sort_key(rec):
    if rec.get("compact"):  # one record per (frag, combo, bucket); fields are arrays
        return rec["frag"], rec["bsr_combo"], rec["aa_bucket"]
    return (rec["frag"], rec["bsr_combo"], rec["aa_bucket"], rec["charge_sign"],
            rec["vdg_index"],
            rec["aa_perm_idx"], rec["q_site_idx"], rec["q_cg_perm_idx"])

def score_one_model(
    pdbfile, pdb_path, lig_smiles, vdg_lib_dir, rmsd_threshold=None,
    ref_lig_mol=None, print_bsr_selection=False, contact_cutoff=None,
    vdg_lib_entries=None, lig_instance_label='',
    nprocs=1, match_mode="joint", compact=False):
    if vdg_lib_entries is None:
        vdg_lib_entries = lib_entries(vdg_lib_dir)

    lig_rmsd_value = None
    match_records = []
    frags_in_lib = {}
    filtered_frags = {}
    buf = io.StringIO()
    bucket_cache = _worker_bucket_cache()

    try:
        with redirect_stdout(buf), redirect_stderr(buf), init_rdkit_logging(buf):
            struct = (pr.parseCIF(pdb_path)
                      if pdb_path.endswith((".cif", ".cif.gz"))
                      else pr.parsePDB(pdb_path))

            try:
                lig_mol_noH = ligand_structure.get_query_ligand_mol(struct, lig_smiles)
            except (ValueError, TypeError) as e:
                print(f"[WARNING] ({pdbfile}) Ligand extraction failed: {e}\n"
                      f"  Skipping this structure.", flush=True)
                return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags

            if lig_mol_noH is None:
                print(f"[WARNING] ({pdbfile}) Could not build RDKit ligand mol "
                      f"from structure; skipping.", flush=True)
                return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags

            if ref_lig_mol is not None:
                if lig_mol_noH.GetNumAtoms() != ref_lig_mol.GetNumAtoms():
                    raise ValueError("Ligand atom count mismatch: "
                        f"query={lig_mol_noH.GetNumAtoms()} "
                        f"ref={ref_lig_mol.GetNumAtoms()}")
                lig_rmsd_value = best_inplace_symmetry_rmsd(ref_lig_mol, lig_mol_noH)

            try:
                _, ligname = ligand_structure.identify_ligand_selection(struct, lig_smiles)
            except ValueError as e:
                print(f"[ERROR] ({pdbfile}) {e}; if the structure has partial "
                      f"occupancies or multiple ligands, filter before running the "
                      f"hit finder.", flush=True)
                return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags

            filtered_frags, query_frag_map, frags_in_lib = match_library_frags_to_query(
                lig_mol_noH, vdg_lib_dir, vdg_lib_entries,
                frags_in_lib=frags_in_lib,
                warn=lambda m: print(f"[WARNING] ({pdbfile}) {m}", flush=True),)

            if not filtered_frags:
                return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags

            all_bsr_combos = dock.get_bsr_combinations(
                struct, ligname, quiet=not print_bsr_selection, pdbfile=pdbfile)

            bb_blockers_cache = {}

            frag_q_cg_perms = {}
            for frag_name, grouped_sites in filtered_frags.items():
                query_frag = query_frag_map.get(frag_name, frag_name)
                cg_order_smarts = cg_atom_order_smarts(vdg_lib_dir, frag_name)

                try:
                    grouped_q_cg_perms = []
                    for site in grouped_sites:
                        perms = []
                        for sub, perm_inds, orig_mol_inds in site:
                            q_cg_coords = np.asarray(dock.get_query_cg_coords(
                                sub, cg_order_smarts), np.float32)
                            q_acceptor = np.array(
                                [a.GetSymbol() in ACCEPTOR_ELEMENTS
                                 for a in sub.GetAtoms()], dtype=bool)
                            perms.append((q_cg_coords, cg_center(q_cg_coords), q_acceptor,
                                          perm_inds, orig_mol_inds))
                        if perms:
                            grouped_q_cg_perms.append(perms)
                except ValueError as e:
                    print(f"[WARNING] ({pdbfile}) Fragment {frag_name!r}: {e}\n"
                          f"  Skipping this fragment.", flush=True)
                    continue

                if not grouped_q_cg_perms:
                    continue

                frag_q_cg_perms[frag_name] = (query_frag, grouped_q_cg_perms)

            tasks = [
                (pdbfile, pdb_path, frag_name, query_frag, grouped_q_cg_perms, combo_item,
                 vdg_lib_dir, rmsd_threshold, contact_cutoff, lig_instance_label, match_mode,
                 compact)
                for frag_name, (query_frag, grouped_q_cg_perms) in frag_q_cg_perms.items()
                for combo_item in all_bsr_combos]

            if nprocs > 1 and len(tasks) > 1:
                ctx = mp.get_context("spawn")
                errors = []
                with ctx.Pool(processes=min(nprocs, len(tasks))) as pool:
                    for local_log, local_records, err in pool.imap_unordered(_combo_worker, tasks):
                        buf.write(local_log)
                        if err:
                            errors.append(err)
                        match_records.extend(local_records)
                _raise_if_combo_task_errors(errors, len(tasks))
                match_records.sort(key=_match_record_sort_key)
            else:
                for frag_name, (query_frag, grouped_q_cg_perms) in frag_q_cg_perms.items():
                    for combo_item in all_bsr_combos:
                        match_records.extend(_score_one_bsr_combo(
                            pdbfile, struct, frag_name, query_frag, grouped_q_cg_perms,
                            combo_item, vdg_lib_dir, bucket_cache,
                            bb_blockers_cache, rmsd_threshold, contact_cutoff,
                            lig_instance_label, match_mode, compact))

        return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags

    except Exception as e:
        log_text = buf.getvalue()
        tb_str = traceback.format_exc()
        raise RuntimeError(
            f"Worker failed for pdbfile={pdbfile}, pdb_path={pdb_path}\n"
            f"--- traceback ---\n{tb_str}\n"
            f"--- worker log start ---\n{log_text}\n"
            f"--- worker log end ---") from e

def score_one_model_multi_instance(pdbfile, pdb_path, lig_smiles, vdg_lib_dir, **kwargs):
    struct = (pr.parseCIF(pdb_path) if pdb_path.endswith((".cif", ".cif.gz"))
              else pr.parsePDB(pdb_path))
    try:
        ligname, instances = ligand_structure.classify_ligand_instances(struct, lig_smiles)
    except ValueError:
        instances = []
    if len(instances) <= 1:
        return score_one_model(pdbfile, pdb_path, lig_smiles, vdg_lib_dir, **kwargs)

    log_texts, match_records, frags_in_lib, filtered_frags = [], [], {}, {}
    for resindex, label in instances:
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tf:
            instance_path = tf.name
        try:
            pr.writePDB(instance_path, ligand_structure.select_struct_for_ligand_instance(
                struct, ligname, [r for r, _ in instances if r != resindex]))
            log_text, recs, this_frags_in_lib, _, this_filtered = score_one_model(
                pdbfile, instance_path, lig_smiles, vdg_lib_dir,
                lig_instance_label=label, **kwargs)
        finally:
            os.remove(instance_path)
        log_texts.append(log_text)
        match_records.extend(recs)
        frags_in_lib.update(this_frags_in_lib)
        for frag_name, sites in this_filtered.items():
            filtered_frags.setdefault(frag_name, []).extend(sites)
    return ''.join(log_texts), match_records, frags_in_lib, None, filtered_frags

def process_work_item(args):
    (pdbfile, pdb_path, lig_smiles, vdg_lib_dir, rmsd_threshold, ref_lig_mol,
     print_bsr_selection, contact_cutoff, match_mode) = args
    vdg_lib_entries = _matching._WORKER_MOL_CACHE_ENTRIES

    if vdg_lib_entries is None:
        vdg_lib_entries = lib_entries(vdg_lib_dir)
        init_worker(vdg_lib_entries)
    try:
        log_text, match_records, frags_in_lib, lig_rmsd_value, _ = score_one_model_multi_instance(
            pdbfile=pdbfile, pdb_path=pdb_path, lig_smiles=lig_smiles,
            vdg_lib_dir=vdg_lib_dir, rmsd_threshold=rmsd_threshold,
            ref_lig_mol=ref_lig_mol, print_bsr_selection=print_bsr_selection,
            contact_cutoff=contact_cutoff, vdg_lib_entries=vdg_lib_entries,
            match_mode=match_mode)
        return pdbfile, log_text, lig_rmsd_value, match_records, frags_in_lib, None
    except Exception:
        error_text = (f"Worker failed for pdbfile={pdbfile}, pdb_path={pdb_path}\n"
                      f"{traceback.format_exc()}")
        return pdbfile, "", None, [], {}, error_text
