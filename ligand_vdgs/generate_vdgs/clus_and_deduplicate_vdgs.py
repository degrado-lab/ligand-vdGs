# clus_and_deduplicate_vdgs.py
'''
Removes redundant vdGs from the vdG library.
A vdG is considered redundant if it meets both criteria after clustering:

1) Binding pose equivalence
   vdGs are clustered by RMSD of [CG + vdM backbone atoms] using a size-normalized cutoff.
   Members belong to the same cluster when their distance ≤ cutoff.

2) Local environment similarity
   Within each pose cluster, vdGs are subclustered on a composite distance:
       `distance = flank-CA RMSD + (sequence dissimilarity)/2`
   with threshold:
       `threshold = size-normalized flank-CA RMSD cutoff + (1 − seq_sim_thresh)/2`, 
   where sequence dissimilarity = 1 − (fractional sequence identity).  

---
Workflow Overview

1. **Scan & stream**  
   Parse each PDB, extract vdGs, and write AA-combo buckets to disk (JSONL gzip).
   Each record is a vdG environment (CG + vdM backbone + flanking seq/CA) written
   to a per–AA–composition bucket file in scratch.

2. **Primary clustering (pose)**  
   Within each AA bucket, cluster on RMSD of [CG + vdM backbone] using a size-scaled cutoff.

3. **Secondary subclustering (local environment)**  
   For each pose cluster, subcluster on  
       `distance = flank-CA RMSD + (sequence dissimilarity)/2`  
   using the threshold defined above.  The flank window is controlled by `--flank`.

4. **De-duplication & output**  
   Merge AA/CG-permutation duplicates, and retain cluster centroids and a few additional 
   cluster members for manual inspection. Cluster centroids are considered nr vdGs.
   PDBs can be reconstructed later from these output arrays.

5. **Cleanup**
   Remove temporary streaming directories.
''' 

import os
import sys
import argparse
import time
import tracemalloc
import shutil
import numpy as np
import prody as pr
import multiprocessing as mp
import json, gzip, hashlib
import traceback
import getpass
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import align_and_cluster as clust
import utils
from clus_helpers import (get_vdg_AA_permutations, combine_cg_and_vdmbb_coords, 
   flatten_flanking_seqs, flatten_flanking_CAs, get_cg_coords, get_vdg_subsets, 
   select_diverse_pdbIDs)

def parse_args():
   parser = argparse.ArgumentParser()
   parser.add_argument('-c', "--cg", type=str, 
                       help="SMARTS pattern of the fragment/CG/FG. Informational only.")
   parser.add_argument('-v', "--vdglib-dir", type=str, required=True,
                       help="Directory for the vdms of this CG.")
   parser.add_argument('-w', "--align-cg-weight", type=float, default=0.99, 
                       help="Weight assigned to CG atoms when superposing output vdGs "
                       "(not for clustering). Example: 0.5 assigns half of the weight " 
                       "to CG atoms and half to vdM backbone atoms. Defaults to 0.99.")
   parser.add_argument('-q', "--seq", default=0.40, type=float, 
                       help="Sequence similarity threshold (fraction between 0 and 1) "
                       "used to set the subclustering threshold. Higher similarity "
                       "indicates redundancy. Defaults to 0.40.")
   parser.add_argument('-f', "--flank", default=2, type=int, 
                       help="Number of residues on each side (+/-) of each vdM to " 
                       "include for computing sequence and backbone similarity. "
                       "Defaults to 2.")
   parser.add_argument('-s', "--symmetry-classes", nargs='+', type=int,
                       help="Integers defining CG symmetry. Ex: carboxylate CC(=O)[O-] "
                       "would be '0 1 2 2' b/c of its equivalent oxygens. If this arg " 
                       "is not provided, the CG is assumed to be asymmetric.")
   parser.add_argument('-n', "--size-subset", type=int, required=True,  
                       help="The number of residues in the vdG subset. The purpose "
                       "of this arg is for parallelization (running with different "
                       "-n values concurrently).")
   parser.add_argument('-l', "--logfile", type=str, default='log', help="Path to log.")
   parser.add_argument('-m', "--max-num-vdgs-to-clus", default=2500, type=int, 
                       help="Cap on # samples per AA-composition bucket to cluster. The cap "
                       "applies to distinct PDB IDs; a single PDB can still contribute "
                       "multiple vdGs after this filter.")
   parser.add_argument('--num-procs', type=int, default=10,
                       help="Number of AA composition buckets to run concurrently.")
   return parser.parse_args()

def _stream_root(vdglib_dir):
    # Root for AA composition streaming buckets.
    # Tries $TMPDIR, /scratch, /tmp, and vdglib_dir as fallback.
    user = os.environ.get("USER") or getpass.getuser() or "unknown"
    vdglib_tag = os.path.basename(os.path.abspath(vdglib_dir.rstrip(os.sep)))
    tmpdir_env = os.environ.get("TMPDIR")
    if tmpdir_env:
        root = os.path.join(tmpdir_env, f"vdg_stream_{vdglib_tag}")
    elif os.path.isdir("/scratch"):
        root = os.path.join("/scratch", user, f"vdg_stream_{vdglib_tag}")
    elif os.path.isdir("/tmp"):
        root = os.path.join("/tmp", user, f"vdg_stream_{vdglib_tag}")
    else:
        root = os.path.join(vdglib_dir, "stream_tmp")
    os.makedirs(root, exist_ok=True)
    return root

def _aa_tmp_dir(vdglib_dir, size_subset):
    # Write per-AA-bucket records to disk to free memory
    return os.path.join(_stream_root(vdglib_dir), str(size_subset))

def _aa_tmp_path(vdglib_dir, size_subset, aa_key):
    # AA bucket file name is based on AA key plus a short SHA1 suffix
    h = hashlib.sha1(aa_key.encode()).hexdigest()[:16]
    fname = f"{aa_key[:80]}__{h}.jsonl.gz"
    return os.path.join(_aa_tmp_dir(vdglib_dir, size_subset), fname)

def _spill_record(vdglib_dir, size_subset, reordered_aas, record):
    '''
    Append one vdG record to its AA bucket file (gzip JSONL).

    Contents of each record:
        {
          "cg_coords": [[x,y,z], ...],
          "bbcoords":  [ [[Nx,Ny,Nz],[CAx,CAy,CAz],[Cx,Cy,Cz]],  ... ],
          "flankseqs": [ ["-2AA","-1AA","vdm","+1AA","+2AA"], ... ],
          "flankCAs":  [ [[x,y,z],...], ... ],
          "pdbpath":   "vdg pdb path",
          "scrr":      [ [seg, chain, resnum, resname], ... ]
        }
    '''
    aa_key = "_".join(reordered_aas)
    path = _aa_tmp_path(vdglib_dir, size_subset, aa_key)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "at", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")

def _load_bucket(path):
    '''
    Load a whole AA bucket JSONL.gz into memory.
    Each vdG entry becomes:
        [cg_coords, bbcoords, flankseqs, flankCAs, pdbpath, scrr]
    '''
    _vdgs = []
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            _vdgs.append([np.asarray(rec["cg_coords"]),
                         [np.asarray(r) for r in rec["bbcoords"]],
                         rec["flankseqs"],
                         [np.asarray(r) for r in rec["flankCAs"]],
                         rec["pdbpath"],
                         rec["scrr"],])
    return _vdgs

def _extract_env(all_cg_coords, all_vdmbb_coords, all_flankseqs, all_flankCAs,
                 all_pdbpaths, all_scrr_cg_perm, idx):
    """
    Extract a single vdG environment (centroid or member) from the global
    arrays built after AA+CG permutations.

    Returns a dict with:
      {
        "cg_coords": (n_cg, 3),
        "vdm_bb_coords": (num_vdms, 3, 3),
        "flank_CA_coords": (num_vdms, L, 3),
        "flank_seq": (num_vdms, L) list-of-lists,
        "pdbpath": str,
        "scrr": [ [seg, chain, resnum, resname], ... ],  # one per vdM
      }
    """
    cg_coords = np.asarray(all_cg_coords[idx], dtype=np.float64)
    vdm_bb_coords = np.asarray(all_vdmbb_coords[idx], dtype=np.float64)
    flank_CAs = np.asarray(all_flankCAs[idx], dtype=np.float64)
    flank_seq = all_flankseqs[idx]
    pdbpath = all_pdbpaths[idx]
    scrrs, _cg_perm_label = all_scrr_cg_perm[idx]   # cg_perm label is informational
    return {"cg_coords": cg_coords,
            "vdm_bb_coords": vdm_bb_coords,
            "flank_CA_coords": flank_CAs,
            "flank_seq": flank_seq,
            "pdbpath": pdbpath,
            "scrr": scrrs,}

def _bucket_npz_path(vdglib_dir, size_subset, reordered_AAs):
    aa_label = "_".join(reordered_AAs)
    nr_root = os.path.join(vdglib_dir, "nr_vdgs")
    size_dir = os.path.join(nr_root, str(size_subset))
    aa_dir = os.path.join(size_dir, aa_label)
    os.makedirs(aa_dir, exist_ok=True)
    return os.path.join(aa_dir, f"{aa_label}.npz")

def _write_bucket_npz(vdglib_dir, size_subset, reordered_AAs, clusters):
    """
    clusters: list of
      {"cluster_id": int,
       "cluster_size": int,
       "first_stage_cluster_id": int,
       "second_stage_cluster_id": int,
       "centroid": env_dict,
       "members": [env_dict, ... up to 29],}
    env_dict as returned by _extract_env.
    """
    if len(clusters) == 0:
        # No clusters for this AA bucket; nothing to write.
        return

    C = len(clusters)
    first_centroid = clusters[0]["centroid"]

    n_cg = first_centroid["cg_coords"].shape[0]
    num_vdms = first_centroid["vdm_bb_coords"].shape[0]
    # Use central vdM flanking length as canonical L
    L = len(first_centroid["flank_seq"][0])

    # Member dimension: up to 29, but allow zero if cluster sizes are 1 everywhere.
    max_members = max(len(c["members"]) for c in clusters)
    M = max_members

    # Core cluster arrays
    cluster_id = np.empty((C,), dtype=int)
    cluster_size = np.empty((C,), dtype=int)
    first_stage_cluster_id = np.empty((C,), dtype=int)
    second_stage_cluster_id = np.empty((C,), dtype=int)

    centroid_cg_coords = np.empty((C, n_cg, 3), dtype=np.float64)
    centroid_vdm_bb_coords = np.empty((C, num_vdms, 3, 3), dtype=np.float64)
    centroid_flank_CA_coords = np.empty((C, num_vdms, L, 3), dtype=np.float64)
    centroid_flank_seq = np.empty((C, num_vdms, L), dtype="U4")
    centroid_parent_pdb_path = np.empty((C,), dtype="U256")
    centroid_scrr_seg = np.empty((C, num_vdms), dtype="U8")
    centroid_scrr_chain = np.empty((C, num_vdms), dtype="U2")
    centroid_scrr_resnum = np.empty((C, num_vdms), dtype=int)
    centroid_scrr_resname = np.empty((C, num_vdms), dtype="U4")

    # Member arrays (may have M = 0; that's fine for numpy)
    if M > 0:
        member_cg_coords = np.full((C, M, n_cg, 3), np.nan, dtype=np.float64)
        member_vdm_bb_coords = np.full((C, M, num_vdms, 3, 3), np.nan, dtype=np.float64)
        member_flank_CA_coords = np.full((C, M, num_vdms, L, 3), np.nan, dtype=np.float64)
        member_flank_seq = np.empty((C, M, num_vdms, L), dtype="U4")
        member_parent_pdb_path = np.empty((C, M), dtype="U256")
        member_scrr_seg = np.empty((C, M, num_vdms), dtype="U8")
        member_scrr_chain = np.empty((C, M, num_vdms), dtype="U2")
        member_scrr_resnum = np.empty((C, M, num_vdms), dtype=int)
        member_scrr_resname = np.empty((C, M, num_vdms), dtype="U4")
    else:
        member_cg_coords = None
        member_vdm_bb_coords = None
        member_flank_CA_coords = None
        member_flank_seq = None
        member_parent_pdb_path = None
        member_scrr_seg = None
        member_scrr_chain = None
        member_scrr_resnum = None
        member_scrr_resname = None

    member_count = np.zeros((C,), dtype=int)

    # Fill arrays
    for i, clus in enumerate(clusters):
        cluster_id[i] = clus["cluster_id"]
        cluster_size[i] = clus["cluster_size"]
        first_stage_cluster_id[i] = clus["first_stage_cluster_id"]
        second_stage_cluster_id[i] = clus["second_stage_cluster_id"]

        cent = clus["centroid"]
        cg = cent["cg_coords"]
        vbb = cent["vdm_bb_coords"]
        fCA = cent["flank_CA_coords"]
        fseq = cent["flank_seq"]
        pdbpath = cent["pdbpath"]
        scrr = cent["scrr"]

        centroid_cg_coords[i] = cg
        centroid_vdm_bb_coords[i] = vbb
        centroid_flank_CA_coords[i] = fCA
        centroid_flank_seq[i] = np.array(fseq, dtype="U4")
        centroid_parent_pdb_path[i] = pdbpath

        # scrr: list of [seg, chain, resnum, resname]
        segs = []
        chains = []
        resnums = []
        resnames = []
        for seg, ch, resnum, resname in scrr:
            segs.append(str(seg))
            chains.append(str(ch))
            resnums.append(int(resnum))
            resnames.append(str(resname))
        centroid_scrr_seg[i] = np.array(segs, dtype="U8")
        centroid_scrr_chain[i] = np.array(chains, dtype="U2")
        centroid_scrr_resnum[i] = np.array(resnums, dtype=int)
        centroid_scrr_resname[i] = np.array(resnames, dtype="U4")

        # Members
        members = clus["members"]
        mcount = len(members)
        member_count[i] = mcount

        if M == 0 or mcount == 0:
            continue

        for j, mem in enumerate(members):
            if j >= M:
                break
            mcg = mem["cg_coords"]
            mvbb = mem["vdm_bb_coords"]
            mfCA = mem["flank_CA_coords"]
            mfseq = mem["flank_seq"]
            mpdb = mem["pdbpath"]
            mscrr = mem["scrr"]

            member_cg_coords[i, j] = mcg
            member_vdm_bb_coords[i, j] = mvbb
            member_flank_CA_coords[i, j] = mfCA
            member_flank_seq[i, j] = np.array(mfseq, dtype="U4")
            member_parent_pdb_path[i, j] = mpdb

            m_segs = []
            m_chains = []
            m_resnums = []
            m_resnames = []
            for seg, ch, resnum, resname in mscrr:
                m_segs.append(str(seg))
                m_chains.append(str(ch))
                m_resnums.append(int(resnum))
                m_resnames.append(str(resname))
            member_scrr_seg[i, j] = np.array(m_segs, dtype="U8")
            member_scrr_chain[i, j] = np.array(m_chains, dtype="U2")
            member_scrr_resnum[i, j] = np.array(m_resnums, dtype=int)
            member_scrr_resname[i, j] = np.array(m_resnames, dtype="U4")

    npz_path = _bucket_npz_path(vdglib_dir, size_subset, reordered_AAs)
    # Save everything in a self-describing npz for downstream materialization.
    save_kwargs = dict(
        cluster_id=cluster_id,
        cluster_size=cluster_size,
        first_stage_cluster_id=first_stage_cluster_id,
        second_stage_cluster_id=second_stage_cluster_id,
        centroid_cg_coords=centroid_cg_coords,
        centroid_vdm_bb_coords=centroid_vdm_bb_coords,
        centroid_flank_CA_coords=centroid_flank_CA_coords,
        centroid_flank_seq=centroid_flank_seq,
        centroid_parent_pdb_path=centroid_parent_pdb_path,
        centroid_scrr_seg=centroid_scrr_seg,
        centroid_scrr_chain=centroid_scrr_chain,
        centroid_scrr_resnum=centroid_scrr_resnum,
        centroid_scrr_resname=centroid_scrr_resname,
        member_count=member_count,)

    if M > 0:
        save_kwargs.update(
            member_cg_coords=member_cg_coords,
            member_vdm_bb_coords=member_vdm_bb_coords,
            member_flank_CA_coords=member_flank_CA_coords,
            member_flank_seq=member_flank_seq,
            member_parent_pdb_path=member_parent_pdb_path,
            member_scrr_seg=member_scrr_seg,
            member_scrr_chain=member_scrr_chain,
            member_scrr_resnum=member_scrr_resnum,
            member_scrr_resname=member_scrr_resname,)

    np.savez(npz_path, **save_kwargs)

def _run_one_bucket_strict(args):
    """
    Run clustering for one AA-composition bucket in a separate process.

    args:
        (bucket_path, seq_sim_thresh, reordered_AAs, symmetry_classes,
         vdglib_dir, align_cg_weight, num_flanking, logfile, size_subset,
         max_num_to_clus)
    """
    (bucket_path, seq_sim_thresh, reordered_AAs, symmetry_classes,
     vdglib_dir, align_cg_weight, num_flanking, logfile, size_subset,
     max_num_to_clus) = args

    _vdgs = _load_bucket(bucket_path)

    # Filter out any vdGs with NaN or inf CG coords
    filtered_vdgs = []
    for (cg_coords, bbcoords, flankseqs, flankCAs, pdbpath, scrr) in _vdgs:
        if not np.isfinite(cg_coords).all():
            with open(logfile, 'a') as f:
                f.write(f'[WARNING] dropping vdG from {pdbpath} from bucket '
                        f'{tuple(reordered_AAs)} due to NaN or inf CG coords.\n')
            continue
        filtered_vdgs.append([cg_coords, bbcoords, flankseqs, flankCAs, pdbpath, scrr])

    _vdgs = filtered_vdgs

    # Cluster vdG subsets with identical vdM AA compositions to identify redundancy.
    # If there are > max_num_to_clus samples, select only the PDBs with highest PDB
    # ID diversity. Note that a single PDB can contribute multiple vdGs, so even
    # after applying the max_num_to_clus cap (which limits PDB IDs, not vdGs), the
    # resulting num. of vdGs may still exceed max_num_to_clus.
    if len(_vdgs) == 0:
        return ('_'.join(reordered_AAs), 0)

    if len(_vdgs) > max_num_to_clus:
        orig_num_vdgs = len(_vdgs)
        pdbIDs = [os.path.basename(z[4])[:4] for z in _vdgs]
        diverse_pdbIDs = select_diverse_pdbIDs(pdbIDs, max_num_to_clus)
        _vdgs = [z for z in _vdgs if os.path.basename(z[4])[:4] in diverse_pdbIDs]
        if orig_num_vdgs != len(_vdgs):
            with open(logfile, 'a') as file:
                file.write(
                    f'\tThere are {orig_num_vdgs} vdgs for the {tuple(reordered_AAs)} '
                    f'subset, so only {len(_vdgs)} vdgs with maximum PDB ID '
                    f'diversity were selected for clustering.\n')

    # If still empty after diversity selection, nothing to cluster.
    if len(_vdgs) == 0:
        return ('_'.join(reordered_AAs), 0)

    try:
        # Expand vdGs over AA perms (identical AAs) and CG symm perms.
        reordered_AAs_str = '_'.join(reordered_AAs)

        (all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords, all_AA_perm_flankingseqs,
         all_AA_perm_flankingCAs, all_AA_perm_pdbpaths,
         all_AA_perm_vdm_scrr) = get_vdg_AA_permutations(reordered_AAs, _vdgs)

        (all_cg_coords, all_vdmbb_coords, all_flankseqs, all_flankCAs, all_pdbpaths,
         all_scrr_cg_perm) = clust.get_vdg_AA_and_cg_perms(all_AA_perm_cg_coords,
         all_AA_perm_vdm_bbcoords, all_AA_perm_flankingseqs, all_AA_perm_flankingCAs,
         all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr, symmetry_classes,)

        # Combined coordinates + flattened flanking info for clustering.
        all_cgvdmbb_coords = combine_cg_and_vdmbb_coords(
            all_cg_coords, all_vdmbb_coords)
        all_flat_flankseqs = flatten_flanking_seqs(all_flankseqs)
        all_flat_flankCAs = flatten_flanking_CAs(all_flankCAs)

        num_cg_atoms = len(all_cg_coords[0])
        num_vdmbb_res = len(all_vdmbb_coords[0])
        num_vdmbb_atoms = num_vdmbb_res * 3
        num_cgvdmbb_atoms = num_cg_atoms + num_vdmbb_atoms

        # Stage 1: cluster on RMSD of [CG + vdM backbone]
        cgvdmbb_rmsd_cut = normalize_rmsd(num_cgvdmbb_atoms, 'cgvdmbb')
        cgvdmbb_data_to_clus = zip([all_cgvdmbb_coords], ['cgvdmbb'])

        cgvdmbb_clus_assignments, cgvdmbb_clus_centroids = clust.get_leader_clusters(
            cgvdmbb_data_to_clus, cgvdmbb_rmsd_cut, reordered_AAs_str,
            len(reordered_AAs), vdglib_dir, final_exact_medoid_pass=True,
            final_reassign_once=True,)

        # Merge degenerate vdGs (different AA/CG perms of the same PDB + SCRR)
        reassigned_cgvdmbb_clus = clust.reassign_cgvdmbb_clusters(
            cgvdmbb_clus_assignments, all_pdbpaths, all_scrr_cg_perm,)

        # Stage 2: within each pose cluster, cluster on residues flanking vdm
        clusters_out = []
        cluster_counter = 0

        for cgvdmbb_clusnum, indices_in_cgvdmbb_cluster in \
                reassigned_cgvdmbb_clus.items():

            indices_in_cgvdmbb_cluster = list(indices_in_cgvdmbb_cluster)
            if len(indices_in_cgvdmbb_cluster) == 0:
                continue

            # Build cluster-specific flattened flanking arrays for Stage 2.
            clus_flat_flankseqs = [all_flat_flankseqs[i]
                                   for i in indices_in_cgvdmbb_cluster]
            clus_flat_flankCAs = [all_flat_flankCAs[i]
                                  for i in indices_in_cgvdmbb_cluster]

            num_flank_bb_coords = len(clus_flat_flankCAs[0])
            flankbb_rmsd_cut = normalize_rmsd(num_flank_bb_coords, 'flankbb')
            flankseq_and_bb_thresh = flankbb_rmsd_cut + (1.0 - seq_sim_thresh) / 2.0

            flankseq_and_bb_to_clus = zip(
                [clus_flat_flankseqs, clus_flat_flankCAs],
                ['flankseq', 'flankbb'],)

            flankingseq_and_bb_cluster_assignments, flankingseq_and_bb_clus_centroids = \
                clust.get_leader_clusters(flankseq_and_bb_to_clus, 
                flankseq_and_bb_thresh, reordered_AAs_str, len(reordered_AAs), 
                vdglib_dir, final_exact_medoid_pass=True, final_reassign_once=True,)

            # Each Stage 2 cluster = final vdG cluster.
            for sub_clus_num, local_indices in \
                    flankingseq_and_bb_cluster_assignments.items():
                if len(local_indices) == 0:
                    continue

                cluster_counter += 1
                centroid_local_idx = flankingseq_and_bb_clus_centroids[sub_clus_num]
                cluster_size = len(local_indices)

                # Map local indices -> global indices into all_* arrays
                global_centroid_idx = indices_in_cgvdmbb_cluster[centroid_local_idx]
                member_global_indices = [
                    indices_in_cgvdmbb_cluster[li]
                    for li in local_indices
                    if li != centroid_local_idx]

                # Deduplicate in original order
                member_global_indices = list(dict.fromkeys(member_global_indices))
                # Cap at 29 stored members per cluster
                member_global_indices = member_global_indices[:29]

                centroid_env = _extract_env(all_cg_coords, all_vdmbb_coords,
                    all_flankseqs, all_flankCAs, all_pdbpaths, all_scrr_cg_perm,
                    global_centroid_idx,)

                member_envs = [_extract_env(all_cg_coords, all_vdmbb_coords,
                    all_flankseqs, all_flankCAs, all_pdbpaths, all_scrr_cg_perm,
                    g_idx,) for g_idx in member_global_indices]

                clusters_out.append({
                    "cluster_id": cluster_counter,
                    "cluster_size": cluster_size,
                    # Track both stages explicitly so we can trace evolution
                    "first_stage_cluster_id": int(cgvdmbb_clusnum),
                    "second_stage_cluster_id": int(sub_clus_num),
                    "centroid": centroid_env,
                    "members": member_envs,})

        # Write this AA bucket summary to npz
        _write_bucket_npz(vdglib_dir, size_subset, reordered_AAs, clusters_out)
        total_vdgs = sum(c["cluster_size"] for c in clusters_out)
        return ('_'.join(reordered_AAs), total_vdgs)

    except Exception as e:
        aa_label = '_'.join(reordered_AAs) if reordered_AAs else 'UNKNOWN_AAs'
        tb = traceback.format_exc()
        raise RuntimeError(
            f"[WORKER ERROR] AA bucket: {aa_label}\n"
            f"Exception: {e}\nTraceback:\n{tb}")

def normalize_rmsd(num_atoms, atoms):
    ''' Return a size-normalized RMSD threshold (Å) for the given atom set
    ('cgvdmbb' or 'flankbb'). The threshold scales linearly with the number
    of atoms between 8 and 15:
        flankbb: 0.5 Å → 1.5 Å
        cgvdmbb: 0.5 Å → 1.0 Å
    Below 8 atoms, use the minimum; above 15, use the maximum.
    '''

    if atoms == 'flankbb':
        max_threshold = 1.5
        min_threshold = 0.5
    elif atoms == 'cgvdmbb':
        max_threshold = 1.0
        min_threshold = 0.5
    else:
        raise ValueError(f"Unknown atom set for normalize_rmsd: {atoms}")

    min_atoms = 8
    max_atoms = 15

    if num_atoms < min_atoms:
        return min_threshold
    if num_atoms > max_atoms:
        return max_threshold

    scaling_factor = (num_atoms - min_atoms) / (max_atoms - min_atoms)
    threshold = min_threshold + scaling_factor * (max_threshold - min_threshold)
    return threshold

def main():
    start_time = time.time()
    tracemalloc.start()
    args = parse_args()

    CG = args.cg
    seq_sim_thresh = args.seq
    num_flanking = args.flank
    symmetry_classes = args.symmetry_classes
    size_subset = args.size_subset
    max_num_to_clus = args.max_num_vdgs_to_clus
    logfile = args.logfile
    align_cg_weight = args.align_cg_weight

    if align_cg_weight < 0 or align_cg_weight > 1:
        raise ValueError('align_cg_weight must be a float between 0 and 1.')

    vdglib_dir = args.vdglib_dir
    vdg_pdbs_dir = os.path.join(vdglib_dir, 'vdg_pdbs')
    if not os.path.isdir(vdg_pdbs_dir):
        sys.exit(f"[ERROR] Missing directory: {vdg_pdbs_dir}")
    out_dir = os.path.join(vdglib_dir, 'nr_vdgs')

    # Iterate over the PDBs and CGs that were identified as containing the SMARTS group.
    # Spill each AA bucket to disk instead of holding all in memory.
    vdg_pdbs_in_dir = sorted(os.listdir(vdg_pdbs_dir))
    for pdbname in vdg_pdbs_in_dir:
        pdbpath = os.path.join(vdg_pdbs_dir, pdbname)
        try:
            prody_obj = pr.parsePDB(pdbpath)
        except Exception:
            with open(logfile, 'a') as file:
                file.write(f"Could not parse {pdbpath} \n")
            continue

        cg_coords = get_cg_coords(prody_obj, pdbpath)
        if cg_coords is None:  # already logged reason inside get_cg_coords()
            continue

        # Define symmetry class if it's None so it's compatible with downstream
        # CG-labelling logic (informational; clustering is orientation-invariant).
        if symmetry_classes is None:
            symmetry_classes = [i for i in range(len(cg_coords))]  # global mutation intentional

        # Characterize the vdM residues (bb coords, flanking residues, pdb paths, etc.)
        vdms_dict = clust.get_vdm_res_features(prody_obj, pdbpath, num_flanking)

        # Determine the vdM combinations, up to size_subset residues.
        vdm_resinds = list(vdms_dict.keys())
        vdg_subsets = get_vdg_subsets(vdm_resinds, size_subset)

        # Iterate over subsets of residue indices
        for vdg_subset in vdg_subsets:
            # Record these features in the same order as in vdg_subset. Then, sort all
            # based on alphabetical order of the vdM AAs. If you find a bb-only vdM,
            # assign the vdM resname as "bb".
            (re_ordered_aas, re_ordered_bbcoords, re_ordered_flankingseqs,
             re_ordered_CAs, re_ordered_scrr) = clust.reorder_vdg_subset(vdg_subset,
             vdms_dict, prody_obj.select('occupancy > 2.9 and not element H'),
             prody_obj,)  # exclude H when determining bb vs sc vdM

            # AA bucket stream: spill a single vdG record to its AA-bucket file on disk.
            rec = {"cg_coords": np.asarray(cg_coords).tolist(),
                   "bbcoords": [np.asarray(x).tolist() for x in re_ordered_bbcoords],
                   "flankseqs": re_ordered_flankingseqs,
                   "flankCAs": [np.asarray(x).tolist() for x in re_ordered_CAs],
                   "pdbpath": str(pdbpath),
                   "scrr": [[str(s), str(ch), int(r), str(rn)]
                            for (s, ch, r, rn) in re_ordered_scrr],}

            _spill_record(vdglib_dir, size_subset, re_ordered_aas, rec)

    # Evaluate the collection of vdGs of size_subset and determine redundancy.
    # Each per-AA bucket is processed independently in parallel.
    stream_dir = _aa_tmp_dir(vdglib_dir, size_subset)
    if not os.path.isdir(stream_dir):
        with open(logfile, 'a') as f:
            f.write(f"\t[STREAM ERROR] No buckets found at {stream_dir}\n")
    else:
        jobs = []
        for fname in sorted(os.listdir(stream_dir)):
            if not fname.endswith(".jsonl.gz"):
                continue
            aa_key = fname.split("__")[0]
            _reordered_AAs = aa_key.split("_")
            path = os.path.join(stream_dir, fname)

            jobs.append((path, seq_sim_thresh, list(_reordered_AAs), symmetry_classes, 
                vdglib_dir, align_cg_weight, num_flanking, logfile, size_subset,
                max_num_to_clus,))

        # Clean, isolated processes
        ctx = mp.get_context("spawn")
        try:
            with ctx.Pool(processes=int(args.num_procs), maxtasksperchild=1) as pool:
                for _ in pool.imap_unordered(_run_one_bucket_strict, jobs, chunksize=1):
                    pass
        except Exception as e:
            err_text = f'[FATAL ERROR] Worker failed with exception:\n{e}\n'
            with open(logfile, 'a') as f:
                f.write(err_text)
            print(err_text, file=sys.stderr)
            sys.exit(1) # hard-fail the whole run if any worker fails

    # Remove only this run's temp stream dir so concurrent -n runs don't clobber
    # each other.
    try:
        if os.path.isdir(stream_dir):
            shutil.rmtree(stream_dir)
    except Exception as _e:
        with open(logfile, 'a') as f:
            f.write(f'\t[STREAM WARNING]: Failed to remove {stream_dir}: {_e}\n')

    # Delete parent stream_root only if empty (i.e., last -n job).
    stream_root = _stream_root(vdglib_dir)
    try:
        os.rmdir(stream_root)
    except FileNotFoundError:
        pass  # another process already deleted it
    except OSError:
        pass  # not empty; other -n jobs still running

    # Log summary
    s = time.time() - start_time
    hours = int(s // 3600)
    minutes = int((s % 3600) // 60)
    seconds = round(s % 60, 2)

    num_vdg_pdbs = len(vdg_pdbs_in_dir)
    nr_dir_of_size_subset = os.path.join(out_dir, str(size_subset))

    # Count only centroids (i.e., clusters) as "nonredundant vdGs"
    num_centroids_of_size_subset = 0

    if os.path.isdir(nr_dir_of_size_subset):
        for aa_label in os.listdir(nr_dir_of_size_subset):
            aa_dir = os.path.join(nr_dir_of_size_subset, aa_label)
            npz_path = os.path.join(aa_dir, f"{aa_label}.npz")
            if not os.path.isfile(npz_path):
                continue
            try:
                data = np.load(npz_path)
                # One centroid per cluster: cluster_id has shape (C,)
                num_centroids_of_size_subset += data["cluster_id"].shape[0]
            except Exception:
                continue

    with open(logfile, 'a') as file:
        file.write(
            f"\nCompleted clus_and_deduplicate_vdgs.py for subset size {size_subset} "
            f"in {hours} h, {minutes} mins, and {seconds} secs.\n")
        file.write(
            f"\t{num_centroids_of_size_subset} nonredun. vdgs (cluster centroids) of "
            f"subset size {size_subset} out of {num_vdg_pdbs} input vdgs.\n")

    current, peak = tracemalloc.get_traced_memory()
    peak_mem = np.round(peak / (1024 * 1024 * 1024), 2)
    if peak_mem > 20:
        print(f"Peak memory usage at end of script for clustering subset size "
              f"{size_subset}: {peak_mem} GB")
    tracemalloc.stop()

if __name__ == "__main__":
    main()
