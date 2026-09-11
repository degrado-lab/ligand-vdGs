# clus_helpers.py

import os
import numpy as np
from itertools import combinations
import getpass
from ligand_vdgs.functions.vdg_struct_utils import (FLANK_UNCOMPARABLE,
    is_valid_backbone_coords)

# One AA-composition bucket is a set of equal-length columns, one row per vdG.
# The Stage-1 array lives in its own .npy so clustering workers can mmap it;
# everything else is one npz. Shards flushed by streaming workers hold all
# columns in one npz and are concatenated at merge.
STAGE1_FILE = "cgvdmbb.npy"
COLUMNS_FILE = "columns.npz"
COLUMN_DTYPES = {
    "cgvdmbb": np.float32,      # (N, n_cg + 3*n_res, 3): CG | N,CA,C per slot
    "flank_ca": np.float32,     # (N, n_res*(2*flank+1), 3), NaN where unreadable
    "flank_seq": "U4",          # (N, n_res*(2*flank+1)) resname / '-' / '!' / 'vdm'
    "biounit": "U32",           # (N,) parent stem (parent_db)
    "scrr_seg": "U8", "scrr_chain": "U2", "scrr_resnum": np.int32, "scrr_resname": "U4",
    "cg_names": "U4", "cg_elements": "U2",
    "cg_seg": "U8", "cg_chain": "U2", "cg_resnum": np.int32, "cg_resname": "U4",
    "slot_flags": np.int8,      # (N, n_res) SLOT_* codes
    "quality": np.float32,      # (N, 4) cg_max_b, cg_min_occ, vdm_max_b, vdm_min_occ
    "vdm_o": np.float32,        # (N, n_res, 3) backbone carbonyl O per slot, NaN when
                                # absent. Not a Stage-1 atom: never in cgvdmbb / RMSD.
}


def stage1_record_is_complete(record, n_cg, n_res):
    """Whether a streamed record has every mandatory finite Stage-1 atom."""
    cg = np.asarray(record["cg_coords"], dtype=np.float32)
    bb = np.asarray(record["bbcoords"], dtype=np.float32)
    if cg.shape != (n_cg, 3) or not np.isfinite(cg).all():
        return False
    if bb.shape != (n_res, 3, 3):
        return False
    return all(is_valid_backbone_coords(res) for res in bb)


def records_to_columns(records):
    """Columns from the per-vdG dicts a streaming worker buffers.

    A record carries ``cg_coords`` (n_cg, 3), ``bbcoords`` (n_res, 3, 3),
    ``flankseqs``/``flankCAs`` (n_res lists of 2*flank+1 tokens / coords),
    ``biounit``, ``scrr`` (n_res of [seg, chain, resnum, resname]),
    ``cg_names``, ``cg_elements``, ``cg_seg``, ``cg_chain``, ``cg_resnum``,
    ``cg_resname``, ``slot_flags`` (n_res), ``quality`` (4) and ``bbo`` (n_res
    carbonyl-O coords, NaN where absent). Every record in a
    bucket has the same shape; a ragged one raises here rather than surfacing as
    a broadcast error.
    """
    n = len(records)
    r0 = records[0]
    n_cg, n_res = len(r0["cg_coords"]), len(r0["bbcoords"])
    n_flank = n_res * len(r0["flankseqs"][0])
    cols = {
        "cgvdmbb": np.empty((n, n_cg + 3 * n_res, 3), dtype=np.float32),
        "flank_ca": np.empty((n, n_flank, 3), dtype=np.float32),
        "flank_seq": np.empty((n, n_flank), dtype="U4"),
        "biounit": np.empty(n, dtype="U32"),
        "scrr_seg": np.empty((n, n_res), dtype="U8"),
        "scrr_chain": np.empty((n, n_res), dtype="U2"),
        "scrr_resnum": np.empty((n, n_res), dtype=np.int32),
        "scrr_resname": np.empty((n, n_res), dtype="U4"),
        "cg_names": np.empty((n, n_cg), dtype="U4"),
        "cg_elements": np.empty((n, n_cg), dtype="U2"),
        "cg_seg": np.empty(n, dtype="U8"),
        "cg_chain": np.empty(n, dtype="U2"),
        "cg_resnum": np.empty(n, dtype=np.int32),
        "cg_resname": np.empty(n, dtype="U4"),
        "slot_flags": np.empty((n, n_res), dtype=np.int8),
        "quality": np.empty((n, 4), dtype=np.float32),
        "vdm_o": np.empty((n, n_res, 3), dtype=np.float32),
    }
    longest = {key: 0 for key, arr in cols.items() if arr.dtype.kind == "U"}
    for k, rec in enumerate(records):
        try:
            cols["cgvdmbb"][k, :n_cg] = rec["cg_coords"]
            cols["cgvdmbb"][k, n_cg:] = np.asarray(
                rec["bbcoords"], dtype=np.float32).reshape(3 * n_res, 3)
            cols["flank_ca"][k] = np.asarray(
                rec["flankCAs"], dtype=np.float32).reshape(n_flank, 3)
            tokens = [tok for res in rec["flankseqs"] for tok in res]
            cols["flank_seq"][k] = tokens
            longest["flank_seq"] = max(longest["flank_seq"], *map(len, tokens))
            for j, (seg, chain, resnum, resname) in enumerate(rec["scrr"]):
                cols["scrr_seg"][k, j] = seg
                cols["scrr_chain"][k, j] = chain
                cols["scrr_resnum"][k, j] = resnum
                cols["scrr_resname"][k, j] = resname
                longest["scrr_seg"] = max(longest["scrr_seg"], len(seg))
                longest["scrr_chain"] = max(longest["scrr_chain"], len(chain))
                longest["scrr_resname"] = max(longest["scrr_resname"], len(resname))
            cols["cg_names"][k] = rec["cg_names"]
            cols["cg_elements"][k] = rec["cg_elements"]
            longest["cg_names"] = max(longest["cg_names"], *map(len, rec["cg_names"]))
            longest["cg_elements"] = max(longest["cg_elements"],
                                         *map(len, rec["cg_elements"]))
            cols["slot_flags"][k] = rec["slot_flags"]
            cols["quality"][k] = rec["quality"]
            cols["vdm_o"][k] = np.asarray(rec["bbo"], dtype=np.float32).reshape(n_res, 3)
        except ValueError as exc:
            raise ValueError(f"vdG record {k} does not match record 0's shape "
                             f"(n_cg={n_cg}, n_res={n_res}, flank={n_flank}): {exc}")
        cols["cg_resnum"][k] = rec["cg_resnum"]
        for key in ("biounit", "cg_seg", "cg_chain", "cg_resname"):
            cols[key][k] = rec[key]
            longest[key] = max(longest[key], len(rec[key]))
    _check_widths(cols, longest)
    return cols


def _check_widths(cols, longest):
    """Fixed-width unicode assignment truncates silently, and a clipped biounit
    stem or label cannot be told from a real one afterwards, so the raw string
    lengths are checked against the dtype width."""
    for key, length in longest.items():
        width = cols[key].dtype.itemsize // 4
        if length > width:
            raise ValueError(f"column {key!r} holds a {length}-char value at "
                             f"dtype width {width}; widen COLUMN_DTYPES")


def concat_columns(parts):
    return {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}


def _replace_into(path, write):
    tmp = f"{path}.{os.getpid()}.tmp"
    try:
        with open(tmp, "wb") as handle:
            write(handle)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def save_shard(path, cols):
    """One npz holding every column; written whole or not at all."""
    _replace_into(path, lambda handle: np.savez(handle, **cols))


def load_shard(path):
    with np.load(path) as data:
        return {key: data[key] for key in data.files}


def save_columns(bucket_dir, cols):
    """The merged bucket: STAGE1_FILE as a bare .npy, the rest in COLUMNS_FILE."""
    os.makedirs(bucket_dir, exist_ok=True)
    _replace_into(os.path.join(bucket_dir, STAGE1_FILE),
                  lambda handle: np.save(handle, cols["cgvdmbb"]))
    rest = {key: arr for key, arr in cols.items() if key != "cgvdmbb"}
    _replace_into(os.path.join(bucket_dir, COLUMNS_FILE),
                  lambda handle: np.savez(handle, **rest))


def load_stage1(bucket_dir, mmap=True):
    return np.load(os.path.join(bucket_dir, STAGE1_FILE),
                   mmap_mode="r" if mmap else None)


def load_columns(bucket_dir, mmap_stage1=False):
    cols = load_shard(os.path.join(bucket_dir, COLUMNS_FILE))
    cols["cgvdmbb"] = load_stage1(bucket_dir, mmap=mmap_stage1)
    return cols


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
