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
   For each fingerprint environment, reconstruct the local vdG AtomGroup from the 
   parent PDB, extract vdGs (CG + vdM backbone + flanking seq/CA), and write AA-combo 
   buckets to disk (JSONL gzip). Each record is a vdG environment written to a 
   per–AA–composition bucket file in scratch.

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
import json, gzip
import traceback
import pickle
def add_paths():
    here = os.path.dirname(os.path.abspath(__file__))
    paths = [
        os.path.join(here, "..", "functions"),
        os.path.join(here, "..", "..", "external", "vdG-miner", "vdg_miner", "programs"),
        os.path.join(here, "..", "..", "external", "vdG-miner", "vdg_miner"),]
    for p in paths:
        p = os.path.abspath(p)
        if p not in sys.path:
            sys.path.append(p)
add_paths()
import align_and_cluster as clust
from clus_helpers import (get_vdg_AA_permutations, combine_cg_and_vdmbb_coords, 
   flatten_flanking_seqs, flatten_flanking_CAs, get_cg_coords, get_vdg_subsets, 
   select_diverse_pdbIDs, _stream_root, _aa_tmp_dir, _aa_tmp_path, normalize_rmsd,)
from fingerprint_helpers import (_pick_atom_by_com, log_warning,
    align_coords_sanity_check, _resolve_duplicate_ligand_occupancies)
from constants import cg_atoms

def parse_args():
   parser = argparse.ArgumentParser()
   parser.add_argument('-c', "--cg", type=str, 
                       help="SMARTS pattern of the fragment/CG/FG. Informational only.")
   parser.add_argument('-v', "--vdglib-dir", type=str, required=True,
                       help="Directory for the vdms of this CG.")
   parser.add_argument('-P', "--pdb-dir", type=str, required=True,
                       help="Parent PDB database directory (RCSB-style mirror).")
   parser.add_argument('-F', "--fingerprints-dir", type=str, required=True,
                       help="Directory containing fingerprints/ produced by "
                            "generate_fingerprints.py (the fp_out_root).")
   parser.add_argument('--cg-match-dict-pkl', type=str,
                       help="Path to the pickled CG match dictionary if the CG "
                            "is not proteinaceous (output from smarts_to_cgs.py).")
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

def _spill_record(vdglib_dir, size_subset, reordered_aas, record):
    '''
    Append one vdG record to its AA bucket file (gzip JSONL).

    Contents of each record:
        {
          "cg_coords": [[x,y,z], ...],
          "bbcoords":  [ [[Nx,Ny,Nz],[CAx,CAy,CAz],[Cx,Cy,Cz]],  ... ],
          "flankseqs": [ ["-2AA","-1AA","vdm","+1AA","+2AA"], ... ],
          "flankCAs":  [ [[x,y,z],...], ... ],
          "pdbpath":   "source identifier (e.g. biounit_envID)",
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
        "cg_coords": (n_cg, 3) float32 array,
        "vdm_bb_coords": (num_vdms, 3, 3) float32 array,
        "flank_CA_coords": (num_vdms, L, 3) float32 array,
        "flank_seq": (num_vdms, L) list-of-lists,
        "pdbpath": str,
        "scrr": [ [seg, chain, resnum, resname], ... ],  # one per vdM
      }
    """
    cg_coords = np.asarray(all_cg_coords[idx], dtype=np.float32)
    vdm_bb_coords = np.asarray(all_vdmbb_coords[idx], dtype=np.float32)
    flank_CAs = np.asarray(all_flankCAs[idx], dtype=np.float32)
    flank_seq = all_flankseqs[idx]
    pdbpath = all_pdbpaths[idx]
    scrrs, _cg_perm_label = all_scrr_cg_perm[idx]   # cg_perm label is informational
    return {
        "cg_coords": cg_coords,
        "vdm_bb_coords": vdm_bb_coords,
        "flank_CA_coords": flank_CAs,
        "flank_seq": flank_seq,
        "pdbpath": pdbpath,
        "scrr": scrrs,
    }


def _bucket_npz_path(vdglib_dir, size_subset, reordered_AAs):
    # Place each AA bucket's npz in nr_vdgs/<size_subset>/<AA_bucket>.npz
    aa_label = "_".join(reordered_AAs)
    nr_root = os.path.join(vdglib_dir, "nr_vdgs")
    size_dir = os.path.join(nr_root, str(size_subset))
    os.makedirs(size_dir, exist_ok=True)
    return os.path.join(size_dir, f"{aa_label}.npz")

def _write_bucket_npz(vdglib_dir, size_subset, reordered_AAs, clusters):
    """
    Write one AA-bucket’s clustering results to a compact .npz file.

    Stores centroid information for each final cluster:
        - cluster IDs and sizes
        - first/second stage cluster IDs
        - centroid CG coords (float32)
        - centroid vdM backbone coords (float32)
        - centroid flanking CA coords (float32)
        - centroid flanking sequence
        - centroid parent PDB path
        - centroid SCRR info (seg/chain/resnum/resname)

    Output path:
        vdglib_dir/nr_vdgs/<size_subset>/<AA_bucket>.npz
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

    # Core cluster arrays
    cluster_id = np.empty((C,), dtype=int)
    cluster_size = np.empty((C,), dtype=int)
    first_stage_cluster_id = np.empty((C,), dtype=int)
    second_stage_cluster_id = np.empty((C,), dtype=int)

    # Store coordinates as float32 to reduce file size and I/O time
    centroid_cg_coords = np.empty((C, n_cg, 3), dtype=np.float32)
    centroid_vdm_bb_coords = np.empty((C, num_vdms, 3, 3), dtype=np.float32)
    centroid_flank_CA_coords = np.empty((C, num_vdms, L, 3), dtype=np.float32)

    centroid_flank_seq = np.empty((C, num_vdms, L), dtype="U4")
    centroid_parent_pdb_path = np.empty((C,), dtype="U256")
    centroid_scrr_seg = np.empty((C, num_vdms), dtype="U8")
    centroid_scrr_chain = np.empty((C, num_vdms), dtype="U2")
    centroid_scrr_resnum = np.empty((C, num_vdms), dtype=int)
    centroid_scrr_resname = np.empty((C, num_vdms), dtype="U4")

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

        centroid_cg_coords[i] = np.asarray(cg, dtype=np.float32)
        centroid_vdm_bb_coords[i] = np.asarray(vbb, dtype=np.float32)
        centroid_flank_CA_coords[i] = np.asarray(fCA, dtype=np.float32)
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

    npz_path = _bucket_npz_path(vdglib_dir, size_subset, reordered_AAs)
    # Save everything in a self-describing npz for downstream materialization.
    np.savez(
        npz_path,
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
    )

def _run_one_bucket_strict(args):
    """
    Run clustering for one AA-composition bucket in a separate process.

    args:
        (bucket_path, seq_sim_thresh, reordered_AAs, symmetry_classes,
         vdglib_dir, num_flanking, logfile, size_subset,
         max_num_to_clus)
    """
    (bucket_path, seq_sim_thresh, reordered_AAs, symmetry_classes,
     vdglib_dir, num_flanking, logfile, size_subset,
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

    # Nothing left after filtering
    if len(_vdgs) == 0:
        return ('_'.join(reordered_AAs), 0)

    # Cluster vdG subsets with identical vdM AA compositions to identify redundancy. 
    # If there are > max_num_to_clus samples, select only the PDBs with highest PDB ID
    # diversity. Note that a single PDB can contribute multiple vdGs, so even after  
    # applying the max_num_to_clus cap (which limits PDBs, not vdGs), the resulting num. 
    # of vdGs may still exceed max_num_to_clus.
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

    try:
        # Expand vdGs over AA perms (identical AAs) and CG symm perms.
        reordered_AAs_str = '_'.join(reordered_AAs)

        (all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords,
         all_AA_perm_flankingseqs, all_AA_perm_flankingCAs,
         all_AA_perm_pdbpaths,
         all_AA_perm_vdm_scrr) = get_vdg_AA_permutations(reordered_AAs, _vdgs)

        (all_cg_coords, all_vdmbb_coords, all_flankseqs, all_flankCAs,
         all_pdbpaths, all_scrr_cg_perm) = clust.get_vdg_AA_and_cg_perms(
            all_AA_perm_cg_coords,
            all_AA_perm_vdm_bbcoords,
            all_AA_perm_flankingseqs,
            all_AA_perm_flankingCAs,
            all_AA_perm_pdbpaths,
            all_AA_perm_vdm_scrr,
            symmetry_classes,)

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
            cgvdmbb_data_to_clus,
            cgvdmbb_rmsd_cut,
            reordered_AAs_str,
            len(reordered_AAs),
            vdglib_dir,
            final_exact_medoid_pass=True,
            final_reassign_once=True,)

        # Merge degenerate vdGs (different AA/CG perms of the same PDB + SCRR)
        reassigned_cgvdmbb_clus = clust.reassign_cgvdmbb_clusters(
            cgvdmbb_clus_assignments, all_pdbpaths, all_scrr_cg_perm)

        # Stage 2: within each pose cluster, cluster on residues flanking vdM
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
            flankseq_and_bb_thresh = (flankbb_rmsd_cut + (1.0 - seq_sim_thresh) / 2.0)

            flankseq_and_bb_to_clus = zip(
                [clus_flat_flankseqs, clus_flat_flankCAs],
                ['flankseq', 'flankbb'],)

            (flankingseq_and_bb_cluster_assignments,
             flankingseq_and_bb_clus_centroids) = clust.get_leader_clusters(
                flankseq_and_bb_to_clus,
                flankseq_and_bb_thresh,
                reordered_AAs_str,
                len(reordered_AAs),
                vdglib_dir,
                final_exact_medoid_pass=True,
                final_reassign_once=True,
            )

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

                centroid_env = _extract_env(
                    all_cg_coords, all_vdmbb_coords,
                    all_flankseqs, all_flankCAs, all_pdbpaths,
                    all_scrr_cg_perm, global_centroid_idx,
                )

                clusters_out.append({
                    "cluster_id": cluster_counter,
                    "cluster_size": cluster_size,
                    "first_stage_cluster_id": int(cgvdmbb_clusnum),
                    "second_stage_cluster_id": int(sub_clus_num),
                    "centroid": centroid_env,
                })

        # Write this AA bucket summary to npz (centroids only)
        _write_bucket_npz(vdglib_dir, size_subset, reordered_AAs, clusters_out)
        total_vdgs = sum(c["cluster_size"] for c in clusters_out)
        return ('_'.join(reordered_AAs), total_vdgs)

    except Exception as e:
        aa_label = '_'.join(reordered_AAs) if reordered_AAs else 'UNKNOWN_AAs'
        tb = traceback.format_exc()
        raise RuntimeError(
            f"[WORKER ERROR] AA bucket: {aa_label}\n"
            f"Exception: {e}\nTraceback:\n{tb}")

def _get_atomgroup_for_env(environment, pdb_dir, cg, cg_match_dict,
                           align_atoms, logfile):
    #Reconstruct the local environment AtomGroup for one fingerprint environment.
    biounit = environment[0][0]
    middle_two = biounit[1:3].lower()
    pdb_file = os.path.join(pdb_dir, middle_two, biounit + '.pdb')

    whole_struct = pr.parsePDB(pdb_file)

    scrs = [(tup[1], tup[2], '`{}`'.format(tup[3])) if tup[3] < 0 else 
            (tup[1], tup[2], tup[3]) for tup in environment]
    selstr_template = '(segment {} and chain {} and resnum {})'
    selstr_template_noseg = '(chain {} and resnum {})'
    selstrs = [selstr_template.format(*scr) if len(scr[0]) else
               selstr_template_noseg.format(*scr[1:]) for scr in scrs]
    sel = whole_struct.select(
        'same residue as within 5 of ({})'.format(' or '.join(selstrs[1:])))
    if sel is None:  # neighborhood selection empty; skip environment
        return None
    struct = sel.toAtomGroup()
    resnames = []
    align_coords = np.zeros((3, 3))

    # Track how many distinct residues we actually map each environment SCR to.
    for i, (scr, selstr) in enumerate(zip(scrs, selstrs)):
        try:
            substruct = struct.select(selstr)
            if substruct is None:  # selection empty; skip env
                return None
            resnames.append(substruct.getResnames()[0])

            # Count unique residue indices for this SCR selection to detect duplication.
            unique_res_indices = np.unique(substruct.getResindices())
            if len(unique_res_indices) != 1:
                log_warning(
                    f'[WARNING] Ambiguous residue selection in {biounit} chain {scr[1]} '
                    f'resnum {scr[2]} (maps to >1 residue); skipping environment.\n',
                    logfile,)
                return None

            if i == 0:
                if cg in cg_atoms.keys():
                    atom_names_list = cg_atoms[cg][resnames[0]]
                else:
                    key = (biounit, scrs[0][0], scrs[0][1],
                           str(scrs[0][2]), resnames[0])

                    match_list = cg_match_dict.get(key)
                    match_idx = environment[0][4] - 1  # 1-based index in env --> 0-based

                    if match_list is None:  # no CG match; possibly missing density; skip
                        return None

                    if not (0 <= match_idx < len(match_list)):  # out of range; possibly 
                                                                 # obabel issue; skip
                        return None
                    atom_names_list = match_list[match_idx]
 
                # Two-pass selection for CG atoms with ambiguity (a PDB with >2 atoms 
                # of the same CG atom name in the same residue)
                cg_atom_selstrs = ['name ' + atom_name for atom_name in atom_names_list]

                # First pass: collect selections and record non-ambiguous atoms
                sel_list = []
                unambig_atoms = []  # atoms with a single unique match

                for j, cg_selstr in enumerate(cg_atom_selstrs):
                    atom_sel = substruct.select(cg_selstr)

                    if atom_sel is None or atom_sel.numAtoms() == 0:  # no atoms; skip
                        return None

                    sel_list.append(atom_sel)
                    if atom_sel.numAtoms() == 1:
                        unambig_atoms.append(atom_sel[0])

                # Compute COM over all unambiguous CG atoms in case they're needed for 
                # disambiguation of atom names belonging to >1 atom
                com = None
                if len(unambig_atoms) > 0:
                    coords_list = []
                    for a in unambig_atoms:
                        c = np.asarray(a.getCoords())
                        # ProDy may return shape (3,) or (1, 3); normalize
                        c = c[0] if c.ndim == 2 else c
                        coords_list.append(c)
                    com = np.mean(coords_list, axis=0)

                for j, atom_sel in enumerate(sel_list):
                    atom_name = atom_names_list[j]
                    candidates = [atom for atom in atom_sel]
                    resname = resnames[0]
                    chain = scrs[0][1]
                    resnum = scrs[0][2]

                    chosen_atom = _pick_atom_by_com(
                        candidates, com, biounit, resname, chain, resnum, atom_name)

                    # If COM check failed (e.g., best_dist > 8 Å), skip this environment.
                    if chosen_atom is None:
                        return None

                    # Set the CG-encoding occupancy for this chosen atom
                    chosen_atom.setOccupancy(3.0 + j * 0.1)

                    # Fill alignment coordinates for the chosen CG atoms
                    if j in align_atoms:
                        c = np.asarray(chosen_atom.getCoords())
                        c = c[0] if c.ndim == 2 else c
                        align_coords[align_atoms.index(j)] = c

            else:
                # Non-CG residues: mark them differently
                substruct.setOccupancies(2.0)

        except Exception:
            # Environment skipped due to exception in selection processing.
            return None

    # Build local frame from align_coords
    if not align_coords_sanity_check(align_coords):  # returns T or F
        # Environment skipped: degenerate local frame.
        return None

    d01 = align_coords[0] - align_coords[1]
    d21 = align_coords[2] - align_coords[1]
    e01 = d01 / np.linalg.norm(d01)
    e21 = d21 / np.linalg.norm(d21)
    e1 = (e01 + e21) / np.linalg.norm(e01 + e21)
    e3 = np.cross(e01, e21) / np.linalg.norm(np.cross(e01, e21))
    e2 = np.cross(e3, e1)
    R = np.array([e1, e2, e3])
    t = align_coords[1]
    coords_transformed = np.dot(struct.getCoords() - t, R.T)
    struct.setCoords(coords_transformed)
    return struct

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
    pdb_dir = args.pdb_dir
    fingerprints_root = args.fingerprints_dir
    cg_match_dict_pkl = args.cg_match_dict_pkl
    vdglib_dir = args.vdglib_dir

    # Load CG match dict for non-proteinaceous CGs (if applicable).
    if cg_match_dict_pkl is not None:
        with open(cg_match_dict_pkl, 'rb') as f:
            cg_match_dict = pickle.load(f)
    else:
        cg_match_dict = {}

    # Deterministic align-atom order (as in fingerprints_to_pdbs.py)
    align_atoms = [1, 0, 2]  # arbitrary, b/c they'll be re-aligned in clustering

    # Scan fingerprints: reconstruct vdG environments and stream to AA buckets.
    fingerprints_dir = os.path.join(fingerprints_root, "fingerprints")
    if not os.path.isdir(fingerprints_dir):
        sys.exit(f"[ERROR] Missing fingerprints dir: {fingerprints_dir}")

    # Build a deterministic, sorted list of input fingerprint files.
    all_fp_files = []
    for subdir in sorted(os.listdir(fingerprints_dir)):
        if '.txt' in subdir:
            continue
        full = os.path.join(fingerprints_dir, subdir)
        if not os.path.isdir(full):
            continue
        for file in sorted(os.listdir(full)):
            if file.endswith('_fingerprints.npy'):
                all_fp_files.append((subdir, file))

    # Iterate over all fingerprint sets; for each environment, reconstruct the local
    # AtomGroup, extract vdGs of size `size_subset`, and stream them to AA buckets.
    for subdir, file in all_fp_files:
        env_path = os.path.join(
            fingerprints_dir, subdir,
            file.replace('_fingerprints.npy', '_environments.txt'))
        if not os.path.isfile(env_path):
            continue

        with open(env_path, 'r') as f_env:
            for env_idx, line in enumerate(f_env):
                environment = eval(line.strip())

                # This is the identifier we used previously for vdg_pdbs filenames.
                pdb_name = '_'.join([str(el) for el in environment[0]])
                pdb_label = pdb_name  # used as "pdbpath" string for downstream code

                # Reconstruct local environment AtomGroup aligned in CG frame.
                atomgroup = _get_atomgroup_for_env(
                    environment, pdb_dir, CG, cg_match_dict, align_atoms, logfile)
                if atomgroup is None:
                    continue

                # Resolve duplicate ligand occupancies using COM; skip env if it
                # cannot be resolved into a clean CG encoding.
                atomgroup = _resolve_duplicate_ligand_occupancies(atomgroup, pdb_label)
                if atomgroup is None:
                    continue

                # Extract CG coordinates from the environment AtomGroup.
                cg_coords = get_cg_coords(atomgroup, pdb_label)
                if cg_coords is None:  # already logged reason inside get_cg_coords()
                    continue

                # Define symmetry class if it's None so it's compatible with downstream
                # CG-labelling logic (informational; clustering is orientation-invariant).
                if symmetry_classes is None:
                    symmetry_classes = [i for i in range(len(cg_coords))]  # global mutation intentional

                # Characterize the vdM residues (bb coords, flanking residues, etc.)
                vdms_dict = clust.get_vdm_res_features(atomgroup, pdb_label, num_flanking)

                # Determine the vdM combinations, up to size_subset residues.
                vdm_resinds = list(vdms_dict.keys())
                vdg_subsets = get_vdg_subsets(vdm_resinds, size_subset)

                # Iterate over subsets of residue indices
                for vdg_subset in vdg_subsets:
                    # Record these features in the same order as in vdg_subset. Then, sort all
                    # based on alphabetical order of the vdM AAs. If you find a bb-only vdM,
                    # assign the vdM resname as "bb".
                    (re_ordered_aas, re_ordered_bbcoords, re_ordered_flankingseqs,
                     re_ordered_CAs, re_ordered_scrr) = clust.reorder_vdg_subset(
                        vdg_subset,
                        vdms_dict,
                        atomgroup.select('occupancy > 2.9 and not element H'),
                        atomgroup,)  # exclude H when determining bb vs sc vdM

                    # AA bucket stream: spill a single vdG record to its AA-bucket file on disk.
                    rec = {
                        "cg_coords": np.asarray(cg_coords).tolist(),
                        "bbcoords": [np.asarray(x).tolist() for x in re_ordered_bbcoords],
                        "flankseqs": re_ordered_flankingseqs,
                        "flankCAs": [np.asarray(x).tolist() for x in re_ordered_CAs],
                        "pdbpath": str(pdb_label),
                        "scrr": [[str(s), str(ch), int(r), str(rn)]
                                 for (s, ch, r, rn) in re_ordered_scrr],}

                    _spill_record(vdglib_dir, size_subset, re_ordered_aas, rec)

    # Evaluate the collection of vdGs of size_subset and determine redundancy.
    # Each per-AA bucket is processed independently in parallel.
    stream_dir = _aa_tmp_dir(vdglib_dir, size_subset) # $TMPDIR or /scratch
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
                vdglib_dir, num_flanking, logfile, size_subset,
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
    stream_root = _stream_root(vdglib_dir) # $TMPDIR or /scratch
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

    # "nonredundant vdGs"  = number of clusters = len(cluster_id)
    # "input environments" = sum of cluster_size over all clusters
    out_dir = os.path.join(vdglib_dir, "nr_vdgs")
    nr_dir_of_size_subset = os.path.join(out_dir, str(size_subset))

    num_centroids_of_size_subset = 0   # number of final clusters
    num_input_envs_of_size_subset = 0  # total members across all clusters

    if os.path.isdir(nr_dir_of_size_subset):
        for fname in os.listdir(nr_dir_of_size_subset):
            if not fname.endswith(".npz"):
                continue

            npz_path = os.path.join(nr_dir_of_size_subset, fname)
            if not os.path.isfile(npz_path):
                continue

            try:
                data = np.load(npz_path)
                # One centroid per cluster
                cluster_ids = data["cluster_id"]
                cluster_sizes = data["cluster_size"]

                num_centroids_of_size_subset += int(cluster_ids.shape[0])
                # cluster_size is the number of environments in each cluster
                num_input_envs_of_size_subset += int(cluster_sizes.sum())
            except Exception:
                # If anything is wrong with this npz, just skip it and continue.
                with open(logfile, 'a') as f:
                    f.write(f'\t[WARNING] Failed to read npz file {npz_path}; skipping.\n')
                continue

    with open(logfile, 'a') as file:
        file.write(
            f"\nCompleted clus_and_deduplicate_vdgs.py for subset size {size_subset} "
            f"in {hours} h, {minutes} mins, and {seconds} secs.\n")
        file.write(
            f"\t{num_centroids_of_size_subset} nonredun. vdgs of "
            f"subset size {size_subset} from {num_input_envs_of_size_subset} "
            f"input environments.\n")

    current, peak = tracemalloc.get_traced_memory()
    peak_mem = np.round(peak / (1024 * 1024 * 1024), 2)
    if peak_mem > 20:
        print(f"Peak memory usage at end of script for clustering subset size "
              f"{size_subset}: {peak_mem} GB")
    tracemalloc.stop()

if __name__ == "__main__":
    main()
