# clus_helpers.py

import os
import numpy as np
from itertools import combinations
import getpass
from ligand_vdgs.functions.vdg_struct_utils import FLANK_UNCOMPARABLE

# Field layout of the per-vdG record lists passed to unpack_vdg_records.
# The producer (clus_and_deduplicate_vdgs._load_bucket) builds records in this
# order; both sides index through these names so adding a field can't silently
# shift the positional reads below.
VDG_FIELDS = ("cg_coords", "bbcoords", "flankseqs", "flankCAs", "pdbpath", "scrr",
              "cg_names", "cg_elements", "cg_seg", "cg_chain", "cg_resnum",
              "cg_resname", "slot_flags",
              "cg_max_b", "cg_min_occ", "vdm_max_b", "vdm_min_occ")
(F_CG_COORDS, F_BBCOORDS, F_FLANKSEQS, F_FLANKCAS, F_PDBPATH, F_SCRR,
 F_CG_NAMES, F_CG_ELEMENTS, F_CG_SEG, F_CG_CHAIN, F_CG_RESNUM,
 F_CG_RESNAME, F_SLOT_FLAGS,
 F_CG_MAX_B, F_CG_MIN_OCC, F_VDM_MAX_B, F_VDM_MIN_OCC) = range(len(VDG_FIELDS))

def unpack_vdg_records(_vdgs):
    """Transpose per-vdG records into one list per VDG_FIELDS column.

    Each ``_vdg`` is a list whose field order is ``VDG_FIELDS`` (indexed here
    through the ``F_*`` constants, not literal positions):

        [cg_coords, bbcoords, flankseqs, flankCAs, pdbpath, scrr,
         cg_names, cg_elements, cg_seg, cg_chain, cg_resnum, cg_resname,
         slot_flags, cg_max_b, cg_min_occ, vdm_max_b, vdm_min_occ]

    ``slot_flags`` is the per-vdM-slot provenance code (``SLOT_*`` in
    vdg_struct_utils), stored in the same slot order as ``scrr``. The four
    trailing quality fields are returned as one ``out_quality`` tuple per record,
    so the output has one column fewer than ``VDG_FIELDS`` has fields.

    Interchangeable same-label vdM slots are **not** expanded into duplicate
    records here. Stage-1 clustering minimises over slot orderings inside its
    distance instead (``vdg_fp_utils.build_perm_group``), which keeps one
    physical environment per record -- so a vdG cannot land in two clusters, and
    no post-hoc union is needed to put the copies back together.

    CG atoms remain in the SMARTS-slot order assigned upstream. Every per-vdG
    sequence in the returned tuple is a fresh object, so no two output records
    share a mutable list.
    """
    for i, _vdg in enumerate(_vdgs):
        if len(_vdg) != len(VDG_FIELDS):
            raise ValueError(
                f"vdG record {i} has {len(_vdg)} fields, expected "
                f"{len(VDG_FIELDS)} ({', '.join(VDG_FIELDS)})")

    total = len(_vdgs)
    out_coords   = [None] * total
    out_bb       = [None] * total
    out_seqs     = [None] * total
    out_CAs      = [None] * total
    out_pdbs     = [None] * total
    out_scrr     = [None] * total
    out_names    = [None] * total
    out_elements = [None] * total
    out_seg      = [None] * total
    out_chain    = [None] * total
    out_resnum   = [None] * total
    out_resname  = [None] * total
    out_flags    = [None] * total
    out_quality  = [None] * total

    for idx, _vdg in enumerate(_vdgs):
        out_coords[idx]   = np.array(_vdg[F_CG_COORDS], dtype=np.float32)
        out_bb[idx]       = list(_vdg[F_BBCOORDS])
        out_seqs[idx]     = list(_vdg[F_FLANKSEQS])
        out_CAs[idx]      = list(_vdg[F_FLANKCAS])
        out_pdbs[idx]     = _vdg[F_PDBPATH]
        out_scrr[idx]     = list(_vdg[F_SCRR])
        out_names[idx]    = list(_vdg[F_CG_NAMES])
        out_elements[idx] = list(_vdg[F_CG_ELEMENTS])
        # Scalars: one residue per CG, so a slot ordering cannot change them.
        out_seg[idx]      = _vdg[F_CG_SEG]
        out_chain[idx]    = _vdg[F_CG_CHAIN]
        out_resnum[idx]   = _vdg[F_CG_RESNUM]
        out_resname[idx]  = _vdg[F_CG_RESNAME]
        out_flags[idx]    = list(_vdg[F_SLOT_FLAGS])
        # Measured B-factor/occupancy over the atoms that enter this vdG, kept
        # per record so a stricter cut can be applied without re-mining.
        out_quality[idx]  = (float(_vdg[F_CG_MAX_B]), float(_vdg[F_CG_MIN_OCC]),
                             float(_vdg[F_VDM_MAX_B]), float(_vdg[F_VDM_MIN_OCC]))

    return (out_coords, out_bb, out_seqs, out_CAs, out_pdbs,
            out_scrr, out_names, out_elements, out_seg, out_chain,
            out_resnum, out_resname, out_flags, out_quality)


def combine_cg_and_vdmbb_coords(all_cg, all_vdmbb):
    """Stack each vdG's CG coords and flattened vdM backbone into one (n, n_cg+n_bb, 3) array.

    Atom counts are taken from record 0 and every other record must match it; a
    ragged input is reported here rather than as a numpy broadcast error from the
    assignment below.
    """
    if not all_cg:
        return []
    if len(all_vdmbb) != len(all_cg):
        raise ValueError(
            f"combine_cg_and_vdmbb_coords: {len(all_cg)} CG records vs "
            f"{len(all_vdmbb)} vdM-backbone records")
    n = len(all_cg)
    n_cg = len(all_cg[0])
    n_vdms = len(all_vdmbb[0])
    n_bb = n_vdms * 3  # num_vdms × 3 backbone atoms (N, CA, C)
    out = np.empty((n, n_cg + n_bb, 3), dtype=np.float32)
    for i, (_cg, _vdmbb) in enumerate(zip(all_cg, all_vdmbb)):
        if len(_cg) != n_cg or len(_vdmbb) != n_vdms:
            raise ValueError(
                f"combine_cg_and_vdmbb_coords: record {i} has {len(_cg)} CG atoms "
                f"and {len(_vdmbb)} vdMs, expected {n_cg} and {n_vdms} (from record 0)")
        out[i, :n_cg] = np.asarray(_cg, dtype=np.float32)
        out[i, n_cg:] = np.asarray([atom for res in _vdmbb for atom in res], dtype=np.float32)
    return out

def flatten_flanking_CAs(cgvdmbb_clus_flankingCAs):
    cgvdmbb_clus_flat_flankCAs = []
    for vdg_flankingCAs in cgvdmbb_clus_flankingCAs:
        # `vdg_flankingCAs` is a list of vdm residues, so flatten it
        flat_flanking_CAs = []
        for res in vdg_flankingCAs:
            for CA_coord in res:
                flat_flanking_CAs.append(np.asarray(CA_coord, dtype=np.float32))
        cgvdmbb_clus_flat_flankCAs.append(flat_flanking_CAs)
    return cgvdmbb_clus_flat_flankCAs

def flatten_flanking_seqs(flankingCAs_clus_flankingseqs):
    flattened_flankingseqs_for_vdgs_in_flankingCA_clus = []
    for vdg_flankingseq in flankingCAs_clus_flankingseqs:
        flat_vdg_flankingseq = []
        for vdm_res in vdg_flankingseq:
            flat_vdg_flankingseq += vdm_res
        flattened_flankingseqs_for_vdgs_in_flankingCA_clus.append(
            flat_vdg_flankingseq)
    return flattened_flankingseqs_for_vdgs_in_flankingCA_clus

def calc_seq_similarity(list1, list2, missing_similarity=0.0):
    """Return coverage-adjusted percent identity for flanking residues.

    ``missing_similarity`` is a percent-identity prior (0..100) assigned to
    positions that are uncomparable in either sequence -- FLANK_MISSING (no
    readable residue) or FLANK_CHAIN_BREAK (past a break, so not a neighbour at
    all). The two carry the same prior here, since neither supplies a residue to
    match, but they are distinct symbols so analysis downstream can separate
    "we could not read it" from "the chain ends here". Central ``'vdm'``
    positions are excluded from both the observed and expected counts.

    Any other token is compared as a residue name, including
    NONCANONICAL_AA_LABEL (``'X'``) -- that is a vdM *slot* label and does not
    appear in flanking sequences, but if it ever does it means a real residue
    with non-canonical atoms, not absent data.
    """
    if not 0.0 <= missing_similarity <= 100.0:
        raise ValueError("missing_similarity must be between 0 and 100")
    list1 = [i for i in list1 if i != 'vdm']
    list2 = [i for i in list2 if i != 'vdm']
    if len(list1) != len(list2):
        raise ValueError(
            f"Sequence similarity got mismatched lengths: {len(list1)} vs {len(list2)}")
    expected = len(list1)
    if expected == 0:
        return float(missing_similarity)
    pairs = [(a, b) for a, b in zip(list1, list2)
             if a not in FLANK_UNCOMPARABLE and b not in FLANK_UNCOMPARABLE]
    matches = sum(1 for a, b in pairs if a == b)
    missing = expected - len(pairs)
    effective_matches = matches + missing * (missing_similarity / 100.0)
    return (effective_matches / expected) * 100.0

def get_vdg_subsets_target_size(input_list, target_size):
    # Generate combos of residues containing `target_size` elements.
    if len(input_list) < target_size:
        return []  
    return list(combinations(input_list, target_size)) 

def select_diverse_pdbIDs(strings, k): # k = max num of distinct PDB IDs to select
   # Greedy max–min Hamming selection over distinct PDB IDs to promote dataset diversity.
   # `strings` may repeat a PDB ID (multiple vdGs from one structure); dedupe first so
   # the k selected indices correspond to k distinct IDs, not k list entries.
   unique_strings = list(dict.fromkeys(strings))
   k = min(k, len(unique_strings))
   ascii_array = strings_to_ascii_array(unique_strings)
   selected_indices = select_diverse_subset_greedy(ascii_array, k)
   return [unique_strings[i] for i in selected_indices]

def strings_to_ascii_array(strings):
   if not strings:
      return np.empty((0, 0), dtype=np.int32)
   max_len = max(len(s) for s in strings)
   return np.array([[ord(c) for c in s.ljust(max_len, '\x00')] for s in strings], dtype=np.int32)

def update_min_dists(data, selected_idx, min_dists):
   # Vectorized over rows. Already-selected rows are updated too (cheaper than
   # masking) and are harmless: select_diverse_subset_greedy masks them out of the
   # argmax, so their min_dists entry is never read.
   dist = (data != data[selected_idx]).sum(axis=1).astype(min_dists.dtype)
   np.minimum(min_dists, dist, out=min_dists)

def select_diverse_subset_greedy(data, k):
   n = data.shape[0]
   if n == 0:
      raise ValueError("data is empty")
   if k <= 0:
      raise ValueError(f"k must be positive, got {k}")
   if k > n:
      raise ValueError(f"k={k} exceeds data size n={n}")
   selected = [0]  # start with first point
   selected_mask = np.zeros(n, dtype=np.uint8)
   selected_mask[0] = 1
   # Sentinel: max possible Hamming distance (data.shape[1]) + 1, avoids wraparound
   sentinel = data.shape[1] + 1
   min_dists = np.full(n, sentinel, dtype=np.int32)

   # Initial distance pass: vectorized hamming from data[0] to all others.
   min_dists[1:] = (data[1:] != data[0]).sum(axis=1).astype(np.int32)

   for _ in range(1, min(n, k)):
      # Select unselected point with largest min-distance to any selected point.
      # Mask selected entries to a value smaller than any reachable distance.
      masked = np.where(selected_mask, np.int32(-1), min_dists.astype(np.int32))
      max_idx = int(np.argmax(masked))

      if masked[max_idx] < 0:  # all points are selected
         break

      selected.append(max_idx)
      selected_mask[max_idx] = 1

      # Update min distances to all unselected points from newly selected point
      update_min_dists(data, max_idx, min_dists)

   return selected

def _stream_root(vdglib_dir, create=False):
    # Root for AA composition streaming buckets. Pure path computation unless
    # create=True: callers that only need the path (e.g. to clean it up) must not
    # have to make the directory as a side effect of asking for it.
    # Tries $TMPDIR, /scratch, /tmp, and vdglib_dir as fallback.
    # Includes SGE job ID for isolation of concurrent runs by the same user.
    user = os.environ.get("USER") or getpass.getuser() or "unknown"
    vdglib_tag = os.path.basename(os.path.abspath(vdglib_dir.rstrip(os.sep)))
    job_id = os.environ.get("JOB_ID", "")
    if job_id:
        # SGE job runs: use job_id (and optionally task_id) for per-invocation isolation
        task_id = os.environ.get("SGE_TASK_ID", "")
        job_suffix = f"{job_id}_{task_id}" if task_id else job_id
    else:
        # Local runs: use PID as fallback
        job_suffix = str(os.getpid())

    tmpdir_env = os.environ.get("TMPDIR")
    if tmpdir_env:
        root = os.path.join(tmpdir_env, user, f"vdg_stream_{vdglib_tag}_{job_suffix}")
    elif os.path.isdir("/scratch"):
        root = os.path.join("/scratch", user, f"vdg_stream_{vdglib_tag}_{job_suffix}")
    elif os.path.isdir("/tmp"):
        root = os.path.join("/tmp", user, f"vdg_stream_{vdglib_tag}_{job_suffix}")
    else:
        root = os.path.join(vdglib_dir, user, f"stream_tmp_{job_suffix}")
    if create:
        os.makedirs(root, exist_ok=True)
    return root

def _aa_tmp_dir(vdglib_dir, size_subset, create=False):
    # Write per-AA-bucket records to disk to free memory
    d = os.path.join(_stream_root(vdglib_dir, create=create), str(size_subset))
    if create:
        os.makedirs(d, exist_ok=True)
    return d
