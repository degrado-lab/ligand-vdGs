# Removes redundant vdGs via two-stage clustering:
# 1) Pose clustering: RMSD of [CG + vdM backbone] with size-normalized cutoff
# 2) Environment subclustering: flank-CA RMSD + sequence dissimilarity with threshold
# Workflow: scan environments → stream to AA buckets → cluster → write nr_vdgs npz

import os
import re
import sys
import argparse
import time
import numpy as np
import multiprocessing as mp
import concurrent.futures
# Explicit: concurrent.futures exposes only the public class names through its
# lazy __getattr__ on 3.10, so 'concurrent.futures.process' is not an attribute
# and referencing it inside an except clause would raise AttributeError exactly
# when a worker has died.
from concurrent.futures.process import BrokenProcessPool
import json
import pickle
import traceback
import shutil
import hashlib
import fcntl
import errno
import itertools
from collections import OrderedDict

# Distinct exit code for "the run completed but streamed nothing". Not 1:
# the wrapper treats it as a warning rather than a crash -- the job exits 0,
# but the fragment log never gets its 'Job completed.' line, so
# check_vdg_job_status still reports it as unfinished and the post-build
# sweep surfaces it instead of it reading as a fragment with no vdGs.
EXIT_NO_VDGS = 3

class _LRUCache:
    """LRU bounded by entry count AND, optionally, by total weight.

    Entry count alone is not a memory bound here: a parsed biounit ranges from a
    few hundred KB to hundreds of MB, so `maxsize` entries of the large kind can
    OOM a worker. `weigh` returns the per-entry cost (atom count) and evictions
    continue until both bounds hold."""
    def __init__(self, maxsize=256, max_weight=None, weigh=None):
        self._cache = OrderedDict()
        self._weights = {}
        self._maxsize = maxsize
        self._max_weight = max_weight
        self._weigh = weigh
        self._total_weight = 0
    def __contains__(self, key):
        return key in self._cache
    def __getitem__(self, key):
        self._cache.move_to_end(key)
        return self._cache[key]
    def __setitem__(self, key, value):
        if key in self._cache:
            self._cache.move_to_end(key)
            self._total_weight -= self._weights.pop(key, 0)
        self._cache[key] = value
        w = self._weigh(value) if self._weigh is not None else 0
        self._weights[key] = w
        self._total_weight += w
        while self._cache and (
                len(self._cache) > self._maxsize
                or (self._max_weight is not None
                    and self._total_weight > self._max_weight
                    and len(self._cache) > 1)):
            evicted, _ = self._cache.popitem(last=False)
            self._total_weight -= self._weights.pop(evicted, 0)

def _log_write(logfile, msg):
    """Append msg to logfile with an exclusive flock to prevent interleaving from concurrent workers."""
    with open(logfile, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            f.write(msg)
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)

# Per-process tally for _log_warn_capped. Workers are spawned, so this is
# per-worker state and needs no locking.
_WARN_COUNTS = {}
_WARN_CAP = 10

def _log_warn_capped(logfile, key, msg, cap=_WARN_CAP):
    """Log a repeating per-environment warning, but only the first `cap` per worker.

    These fire once per environment, and every one takes an exclusive flock on a
    single logfile on NFS. Under a systematic failure -- wrong --pdb-dir, a
    mismatched match dict, a mirror that stopped resolving -- that is one
    serialized network write per environment across every worker, which turns a
    fast failure into a multi-GB log and a run that never finishes. The
    suppressed warnings are not lost information: _stream_one_chunk returns the
    per-key deltas of _WARN_COUNTS and they are logged, per key, as the
    "Warnings during streaming" line at the end of streaming.

    That line is what makes the suppression safe, and it is deliberately kept
    apart from the skip counters. These are two different quantities: a
    residue-level warning (vdm_not_in_contact, vdm_all_hydrogen) drops one slot
    and leaves the environment standing, so it increments nothing in `skips`,
    while the environment-level reasons (duplicate_cg_atom_name,
    cg_not_bonded_component, pdb_parse_failed, pdb_returned_none) all collapse
    into the single `no_atomgroup` bucket. Summing the two sets would double
    count the latter and invent environment skips for the former.
    """
    n = _WARN_COUNTS.get(key, 0) + 1
    _WARN_COUNTS[key] = n
    if n <= cap:
        _log_write(logfile, msg)
    elif n == cap + 1:
        _log_write(logfile,
            f"[WARNING] further '{key}' warnings from this worker are suppressed "
            f"after {cap}; the end-of-streaming counters carry the true total.\n")

def add_vdg_miner_paths():
    # The vdG-miner submodule is not a package, so it still needs sys.path.
    here = os.path.dirname(os.path.abspath(__file__))
    paths = [
        os.path.join(here, "..", "..", "external", "vdG-miner", "vdg_miner", "programs"),
        os.path.join(here, "..", "..", "external", "vdG-miner", "vdg_miner")]
    for p in paths:
        p = os.path.abspath(p)
        if p not in sys.path:
            sys.path.append(p)

add_vdg_miner_paths()
from fingerprint_helpers import (align_coords_sanity_check,
    _resolve_duplicate_ligand_occupancies,)
from constants import cg_atoms

from ligand_vdgs.functions import align_and_cluster as clust
from ligand_vdgs.functions.clus_helpers import (unpack_vdg_records,
    combine_cg_and_vdmbb_coords, flatten_flanking_seqs, flatten_flanking_CAs,
    get_vdg_subsets_target_size, select_diverse_pdbIDs, _aa_tmp_dir,
    _stream_root, VDG_FIELDS,)
from ligand_vdgs.functions.vdg_struct_utils import (VDM_OCC, cg_slot_occupancy,
    get_cg_atoms, is_valid_backbone_coords, is_hydrogen)
from ligand_vdgs.functions.vdg_npz_utils import (write_cg_symmetry, name_selstr,
    parse_pdb_with_retry)
from ligand_vdgs.functions.vdg_fp_utils import build_perm_group, slot_orders
from ligand_vdgs.functions.dock_utils import cg_element_symbols
from ligand_vdgs.functions.compute_profile import ComputeProfile
from ligand_vdgs.functions.utils import (convert_time_elapsed, normalize_rmsd,
    identify_mol_automorphisms, mol_from_fragment, validate_atom_permutations,
    _int_or_none,)

# Max distance from a CG atom to the nearest heavy atom of a vdM residue. Numerically
# the same 4.5 A that align_and_cluster.reorder_vdg_subset uses, but that one only
# *labels* a slot (sidechain vs backbone moiety) and never rejects it, so this is the
# only place a non-contacting slot is actually dropped -- don't expect a second check
# downstream.
CG_VDM_CONTACT_CUTOFF = 4.5

# Longest distance treated as a plausible covalent bond between two CG atoms. Set
# above every bond a drug-like fragment can contain (C-I is 2.14 A, S-S 2.05 A) plus
# room for refinement error, because the check is meant to catch gross perception
# failures -- an atom picked from the wrong residue is tens of A away -- not to
# referee bond geometry.
MAX_CG_BOND_DIST = 2.5


def _cg_bond_components(coords, cutoff=MAX_CG_BOND_DIST):
    """Connected components of the CG atoms under a covalent-bond cutoff.

    A SMARTS match is a connected subgraph, so its atoms must also be connected in
    space. Checked by connectivity rather than by the largest nearest-neighbour
    distance alone, which cannot see two well-formed halves sitting far apart (the
    shape a mis-resolved atom name produces).
    """
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if n < 2:
        return [list(range(n))]
    adj = (((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
           <= cutoff ** 2)
    comps, unseen = [], set(range(n))
    while unseen:
        root = min(unseen)
        seen, stack = {root}, [root]
        while stack:
            for j in np.flatnonzero(adj[stack.pop()]):
                j = int(j)
                if j not in seen:
                    seen.add(j)
                    stack.append(j)
        comps.append(sorted(seen))
        unseen -= seen
    return comps


# One source of truth for both the pickled record keys and the positional order
# unpack_vdg_records reads them back in.
_REQUIRED_VDG_KEYS = set(VDG_FIELDS)

# Streaming chunks per worker process. >1 so the executor can rebalance across
# workers; small enough that the per-chunk dispatch cost stays a rounding error
# against a chunk's parse cost.
STREAM_CHUNKS_PER_WORKER = 4


def _has_complete_stage1_coords(vdg_data, expected_n_cg, expected_num_vdms):
    """Return whether a serialized vdG has all mandatory finite Stage-1 atoms."""
    try:
        cg_coords = np.asarray(vdg_data[0], dtype=np.float32)
        bb_coords = np.asarray(vdg_data[1], dtype=np.float32)
    except (TypeError, ValueError):
        return False
    if cg_coords.shape != (expected_n_cg, 3) or not np.isfinite(cg_coords).all():
        return False
    if bb_coords.shape != (expected_num_vdms, 3, 3):
        return False
    return all(is_valid_backbone_coords(bb) for bb in bb_coords)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cg", type=str, required=True,
        help="Output label for the fragment/CG/FG.")
    parser.add_argument("--cg-smarts", type=str, required=True,
        help="Hydrogen-free fragment SMARTS from which exact CG automorphisms "
        "are derived.",)
    parser.add_argument("-v", "--vdglib-dir", type=str, required=True,
        help="Directory for the vdms of this CG.")
    parser.add_argument("-P", "--pdb-dir", type=str, required=True,
        help="Parent PDB database directory (RCSB-style mirror).")
    parser.add_argument("-E", "--environments-dir", type=str, required=True,
        help="Directory containing environments/ produced by "
        "generate_environments.py (the environment staging root).")
    parser.add_argument("--cg-match-dict-pkl", type=str,
        help="Path to the pickled CG match dictionary if the CG "
        "is not proteinaceous (output from smarts_to_cgs.py).")
    parser.add_argument("-q", "--seq", default=0.40, type=float,
        help="Sequence similarity threshold (fraction between 0 and 1) "
        "used to set the subclustering threshold. Higher similarity "
        "indicates redundancy. Defaults to 0.40.",)
    parser.add_argument("-f", "--flank", default=2, type=int,
        help="Number of residues on each side (+/-) of each vdM to "
        "include for computing sequence and backbone similarity. "
        "Defaults to 2.",)
    parser.add_argument("--subset-sizes", type=int, nargs="+", required=True,
        choices=(1, 2),
        help="Residue counts to generate in one shared environment pass.",)
    parser.add_argument("--profile-json", type=str, default=None,
        help="Optional path for a compute-profile JSON sidecar. Omitted by "
        "default; the wrapper supplies it unless --no-profile-compute.")
    parser.add_argument("-l", "--logfile", type=str, default="log", 
        help="Path to log.")
    parser.add_argument("-m", "--max-num-vdgs-to-clus", default=None, type=_int_or_none,
        help="DEBUGGING ONLY -- never set this for a production library build; "
        "a capped run silently produces an incomplete library. Cap on # samples "
        "per AA-composition bucket to cluster. The cap applies to distinct PDB "
        "IDs; a single PDB can still contribute multiple vdGs after this filter. "
        "Default: no limit.",)
    parser.add_argument("-p", "--num-procs", default=10, type=int,
        help="Number of AA composition buckets to run concurrently.",)
    parser.add_argument("--pdb-cache-size", default=256, type=int,
        help="Max number of parsed PDB structures to keep in memory during "
        "streaming. This is an entry cap only -- a parsed biounit ranges from "
        "well under 1 MB to hundreds of MB, so memory is bounded independently "
        "by a per-worker atom budget (_PDB_CACHE_ATOM_BUDGET, split across "
        "--num-procs). Set to 0 to disable caching.",)
    args = parser.parse_args()
    # 0 (or negative) would make n_workers 0, skip streaming entirely, and write an
    # empty library with a zero exit -- indistinguishable from a fragment with no vdGs.
    if args.num_procs < 1:
        parser.error(f"--num-procs must be at least 1, got {args.num_procs}")
    return args

def _pdb_id_from_path(pdbpath):
    """Recover the bare PDB ID from a `rec["pdbpath"]` value, i.e. the path built as
    `pdb_dir/<biounit[1:3].lower()>/<biounit>.pdb` (see `_stream_one_chunk`), where
    `biounit` is `<pdbid>` or `<pdbid>_<n>`. Keep in sync with that construction."""
    return os.path.splitext(os.path.basename(pdbpath))[0].split("_")[0]

# Well under the 255-byte filename limit on ext4/xfs, with room for the hash,
# the flush index and the extension.
_MAX_BUCKET_FNAME_BYTES = 200


def _bucket_fname(aa_key, flush_idx=None):
    """Name of the on-disk file holding one AA-composition bucket.

    The aa_key is written **whole**, never truncated: it is the only record of the
    bucket's composition, and both readers (`_merge_worker_dirs` and the bucket
    dispatch in `main`) recover it by parsing this name back with
    `_aa_key_from_bucket_fname`. A truncated key would merge distinct buckets into
    one file and mislabel their residues, silently. The sha1 is kept for name stability
    and legibility; it is not what makes the name unique.
    """
    h = hashlib.sha1(aa_key.encode()).hexdigest()[:16]
    stem = f"{aa_key}__{h}" if flush_idx is None else f"{aa_key}__{h}__{flush_idx:04d}"
    fname = f"{stem}.pkl"
    if len(fname.encode()) > _MAX_BUCKET_FNAME_BYTES:
        raise ValueError(
            f"AA bucket key {aa_key!r} ({len(aa_key)} chars) does not fit in a "
            f"{_MAX_BUCKET_FNAME_BYTES}-byte filename. Labels are <=4 chars, so this "
            "needs a subset size around 40 -- far past MAX_SUBSET_SIZE. Store the key "
            "inside the pickle instead of shortening it here; the readers parse it "
            "back out of the name.")
    return fname


_BUCKET_FNAME_RE = re.compile(r"^(?P<key>[^_].*?)__[0-9a-f]{16}(?:__\d{4})?\.pkl$")


def _aa_key_from_bucket_fname(fname):
    """Recover the exact aa_key from a name written by `_bucket_fname`, or None.

    None rather than a raise: the callers walk directories, and a file that is not
    one of ours (an editor backup, a stray from a killed run) must not abort a merge
    over every bucket in the library. 
    """
    match = _BUCKET_FNAME_RE.match(fname)
    return match.group("key") if match else None

def _dump_pickle_atomic(obj, path):
    """Pickle to a temp name in the same directory, then os.replace onto `path`.

    Same contract as the npz writer: a worker killed part-way through the dump
    (the OOM killer, an SGE h_rt kill, a full scratch quota) must never leave a
    truncated file sitting on a name the readers accept, because every reader
    trusts a well-formed name and would fail at unpickling -- after streaming has
    already finished. The temp name deliberately does not end in '.pkl', which is
    what `_merge_worker_dirs` and `_bucket_jobs` walk for.
    """
    tmp = f"{path}.{os.getpid()}.tmp"
    try:
        with open(tmp, "wb") as fh:
            pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)
    except BaseException:
        # BaseException, not Exception: a queue SIGTERM arrives as SystemExit,
        # which is one of the cases this cleanup exists for.
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def _load_bucket(path):
    with open(path, "rb") as f:
        records = pickle.load(f)
    _vdgs = []
    for rec in records:
        if not isinstance(rec, dict):
            raise ValueError(f"Bucket pickle corrupted: record is not a dict ({type(rec)})")
        if not _REQUIRED_VDG_KEYS.issubset(rec.keys()):
            raise ValueError(f"Bucket pickle corrupted: missing keys {_REQUIRED_VDG_KEYS - set(rec.keys())}")
        converted = dict(rec)
        converted["cg_coords"] = np.asarray(rec["cg_coords"], dtype=np.float32)
        converted["bbcoords"] = [np.asarray(r, dtype=np.float32) for r in rec["bbcoords"]]
        converted["flankCAs"] = [np.asarray(r, dtype=np.float32) for r in rec["flankCAs"]]
        _vdgs.append([converted[k] for k in VDG_FIELDS])
    return _vdgs

def _extract_env(all_cg_coords, all_vdmbb_coords,
    all_pdbpaths, all_scrr, all_cg_names, all_cg_elements, all_cg_seg, all_cg_chain,
    all_cg_resnum, all_cg_resname, all_slot_flags, all_quality, idx):
    """The representative record for one cluster, in the shape _write_bucket_npz reads.

    Flanking sequence and flank CA coords are deliberately absent: they drive Stage-2
    subclustering and nothing else, and no npz field holds them.
    """
    return {
        "cg_coords": np.asarray(all_cg_coords[idx], dtype=np.float32),
        "vdm_bb_coords": np.asarray(all_vdmbb_coords[idx], dtype=np.float32),
        "pdbpath": all_pdbpaths[idx], "scrr": all_scrr[idx],
        "slot_flags": list(all_slot_flags[idx]),
        "cg_names": list(all_cg_names[idx]), "cg_elements": list(all_cg_elements[idx]),
        "cg_seg": all_cg_seg[idx], "cg_chain": all_cg_chain[idx],
        "cg_resnum": all_cg_resnum[idx], "cg_resname": all_cg_resname[idx],
        "quality": all_quality[idx]}

def _extract_members(all_pdbpaths, all_scrr, all_cg_names, all_cg_seg,
    all_cg_chain, all_cg_resnum, all_cg_resname, all_quality,
    member_idxs, nr_idx):
    """The cluster's members *excluding* its medoid.

    The nr vdG's identity is already in the nr_ arrays, in full and then
    some (it carries coordinates and a slot flag the member rows do not), so
    emitting it here too would store the same vdG twice.
    """
    return [_extract_member_identity(all_pdbpaths, all_scrr, all_cg_names,
                all_cg_seg, all_cg_chain, all_cg_resnum, all_cg_resname,
                all_quality, idx)
            for idx in member_idxs if idx != nr_idx]

def _bucket_npz_path(vdglib_dir, size_subset, reordered_AAs):
    aa_label = "_".join(reordered_AAs)
    size_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset))
    os.makedirs(size_dir, exist_ok=True)
    return os.path.join(size_dir, f"{aa_label}.npz")

def _write_failed_marker(vdglib_dir, size_subset, aa_label, err_text):
    """Record a bucket that failed to cluster, so a partially built nr_vdgs/ is
    never mistaken for a finished library."""
    size_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset))
    os.makedirs(size_dir, exist_ok=True)
    with open(os.path.join(size_dir, f"{aa_label}.FAILED"), "w") as fh:
        fh.write(err_text + "\n")

def _clear_bucket_outputs(vdglib_dir, size_subset, aa_label):
    """Remove any previous run's npz and .FAILED marker for one bucket.

    Called at bucket entry so the outputs on disk always describe the current
    invocation: a stale npz would be counted by _subset_output_counts even when
    this run's bucket failed, and a stale .FAILED would condemn a library this
    run rebuilt successfully."""
    size_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset))
    for name in (f"{aa_label}.npz", f"{aa_label}.FAILED"):
        try:
            os.remove(os.path.join(size_dir, name))
        except FileNotFoundError:
            pass

def _extract_member_identity(all_pdbpaths, all_scrr, all_cg_names,
    all_cg_seg, all_cg_chain, all_cg_resnum, all_cg_resname, all_quality, idx):
    return {
        "pdbpath": all_pdbpaths[idx],
        "scrr": all_scrr[idx],
        "cg_names": list(all_cg_names[idx]),
        "cg_seg": all_cg_seg[idx], "cg_chain": all_cg_chain[idx],
        "cg_resnum": all_cg_resnum[idx], "cg_resname": all_cg_resname[idx],
        "quality": all_quality[idx],}

def _renumber_clusters_by_size(clusters):
    """Sort a bucket's clusters largest-first and renumber cluster_id 1..N.
    Mutates and returns ``clusters``; _write_bucket_npz writes nr vdGs in list
    order and stamps mem_cluster_id from the same dicts, so reordering here
    keeps the nr rows and the member records consistent.
    """
    clusters.sort(key=lambda clus: -clus["cluster_size"])
    for new_id, clus in enumerate(clusters, start=1):
        clus["cluster_id"] = new_id
    return clusters


def _biounit_of(pdbpath):
    """Biounit stem of a `rec["pdbpath"]`, i.e. the `<biounit>` in
    `pdb_dir/<biounit[1:3].lower()>/<biounit>.pdb`. Unlike `_pdb_id_from_path`
    this keeps any `_<n>` assembly suffix, because it must round-trip back to a
    filename."""
    return os.path.splitext(os.path.basename(pdbpath))[0]


def _parent_dir_of(pdbpath):
    """The `pdb_dir` a `rec["pdbpath"]` was built from (strip `<mid2>/<file>`)."""
    return os.path.dirname(os.path.dirname(pdbpath))


# Order of the per-record quality tuple carried on every vdG record. Written for
# both the nr_ and mem_ row sets, so it lives here rather than being spelled out
# at each of the four sites that allocate or fill it.
_QUALITY_FIELDS = ("cg_max_b", "cg_min_occ", "vdm_max_b", "vdm_min_occ")


def _quality_arrays(prefix, n):
    return {f"{prefix}_{key}": np.empty(n, dtype=np.float32)
            for key in _QUALITY_FIELDS}


def _store_quality(arrays, prefix, row, quality):
    """Fill one row's quality columns. strict=True because the arrays are
    np.empty: a short tuple would otherwise leave uninitialized floats in the
    library with nothing to signal it."""
    for key, value in zip(_QUALITY_FIELDS, quality, strict=True):
        arrays[f"{prefix}_{key}"][row] = value


def _write_bucket_npz(vdglib_dir, size_subset, reordered_AAs, clusters):
    if not clusters: return
    C = len(clusters)
    nr0 = clusters[0]["nr"]
    n_cg, num_vdms = nr0["cg_coords"].shape[0], nr0["vdm_bb_coords"].shape[0]
    parent_pdb_dir = _parent_dir_of(nr0["pdbpath"])

    arrays = {
        # Bucket-level (not per nr vdG): the label of each residue slot. U4 fits
        # every label -- the 20 resnames, 'bb', and 'X'. A longer label would
        # silently truncate here, so widen this alongside adding one.
        "aa_bucket_parts": np.asarray(list(reordered_AAs), dtype="U4"),
        "cluster_id": np.empty(C, dtype=np.int32),
        "cluster_size": np.empty(C, dtype=np.int32),
        # Exact pose radius: the greatest symmetry-aware RMSD from any member to
        # the stored nr vdG. Sphere exclusion bounds this by the Stage-1 cutoff
        # for a whole Stage-1 cluster, but Stage 2 splits those and re-picks a
        # row, so the stored value is measured, not assumed. It is what a
        # consumer needs to reason about recall: a query within `tau` of some
        # discarded member is within `tau + radius` of the nr vdG that replaced it.
        "cluster_pose_radius": np.empty(C, dtype=np.float32),
        # Distinct parent structures behind the cluster, counting the nr vdG.
        # `cluster_size` counts observations, which NCS copies and homologous
        # entries inflate; this is the support figure to do statistics on.
        "cluster_num_parents": np.empty(C, dtype=np.int32),
        "first_stage_cluster_id": np.empty(C, dtype=np.int32),
        "second_stage_cluster_id": np.empty(C, dtype=np.int32),
        "nr_cg_coords": np.empty((C, n_cg, 3), dtype=np.float32),
        "nr_vdm_bb_coords": np.empty((C, num_vdms, 3, 3), dtype=np.float32),
        # Parent structures are stored as the biounit stem ("1f8s", or "1f8s_1"
        # when the mirror carries assembly suffixes), not the absolute path:
        # the path is identical for every record in a library, so storing it in
        # full made two arrays ~80% of the bytes np.load pulls into memory.
        # `parent_pdb_dir` below records the directory once per file so a bucket
        # stays self-describing; readers rebuild the path with
        # vdg_npz_utils.resolve_parent_pdb_path, which takes an override so a
        # relocated PDB mirror (e.g. a collaborator's copy) still resolves.
        "parent_pdb_dir": np.asarray(str(parent_pdb_dir or ""), dtype="U512"),
        "nr_parent_biounit": np.empty(C, dtype="U16"),
        "nr_scrr_seg": np.empty((C, num_vdms), dtype="U8"),
        "nr_scrr_chain": np.empty((C, num_vdms), dtype="U2"),
        "nr_scrr_resnum": np.empty((C, num_vdms), dtype=np.int32),
        "nr_scrr_resname": np.empty((C, num_vdms), dtype="U4"),
        # Atom names are per CG atom: two records can be the same chemical group
        # in differently-named ligands. The CG's *residue* is a single value per
        # record (get_cg_atoms enforces one residue per CG), and its *elements*
        # are one value per bucket because every record retains the fragment's
        # SMARTS-slot order.
        # Both were previously stored per atom per record and were pure padding.
        "nr_cg_names": np.empty((C, n_cg), dtype="U4"),
        "cg_elements": np.empty(n_cg, dtype="U2"),
        "nr_cg_seg": np.empty(C, dtype="U8"),
        "nr_cg_chain": np.empty(C, dtype="U2"),
        "nr_cg_resnum": np.empty(C, dtype=np.int32),
        "nr_cg_resname": np.empty(C, dtype="U4"),
        # Per-vdM-slot provenance (SLOT_* in vdg_struct_utils): which moiety of
        # the canonical residue is nearest the CG, plus a bit for whether the
        # residue carries non-canonical atoms. Disjoint from the bucket label by
        # construction, so nothing here restates the file name; read it together
        # with nr_scrr_resname (SLOT_NO_SC is chemistry on a glycine and
        # missing density anywhere else).
        # CAVEAT: this is the medoid's flag only, not every member's. Clustering
        # never consults it and flank-sequence similarity tolerates mismatches,
        # so a cluster can mix e.g. SLOT_NO_SC and SLOT_BB_CLOSER members;
        # weighting this code by cluster_size is invalid.
        "nr_slot_flag": np.empty((C, num_vdms), dtype=np.int8),
        # Measured on the atoms that enter the vdG -- CG atoms, and the
        # contacting residues' heavy atoms -- not on each residue's first atom.
        # Stored rather than only filtered on so a stricter cut is a read-path
        # decision instead of a reason to re-mine the PDB.
        **_quality_arrays("nr", C),}

    # A "member" here is a clustered observation that is *not* the nr vdG; the nr
    # vdG's own identity lives in the nr_ arrays.
    #
    # `mem_`, not `member_`: these rows are every member of the cluster
    # EXCEPT the nr vdG itself, whose identity lives in the nr_ arrays. Each vdG is
    # therefore stored exactly once. The name is deliberately not `member_` --
    # `cluster_size` counts the medoid, so a reader treating these as "the
    # members" would be off by one per cluster with nothing to signal it.
    M = sum(len(clus["members"]) for clus in clusters)
    expected_M = sum(clus["cluster_size"] for clus in clusters) - C
    if M != expected_M:
        raise ValueError(f"mem_ row count {M} != summed cluster_size "
            f"minus one per cluster ({expected_M}); cluster membership was lost.")
    member_arrays = {
        "mem_cluster_id": np.empty(M, dtype=np.int32),
        "mem_parent_biounit": np.empty(M, dtype="U16"),
        "mem_scrr_seg": np.empty((M, num_vdms), dtype="U8"),
        "mem_scrr_chain": np.empty((M, num_vdms), dtype="U2"),
        "mem_scrr_resnum": np.empty((M, num_vdms), dtype=np.int32),
        "mem_scrr_resname": np.empty((M, num_vdms), dtype="U4"),
        "mem_cg_names": np.empty((M, n_cg), dtype="U4"),
        "mem_cg_seg": np.empty(M, dtype="U8"),
        "mem_cg_chain": np.empty(M, dtype="U2"),
        "mem_cg_resnum": np.empty(M, dtype=np.int32),
        "mem_cg_resname": np.empty(M, dtype="U4"),
        **_quality_arrays("mem", M),}
    m = 0
    for clus in clusters:
        for mem in clus["members"]:
            member_arrays["mem_cluster_id"][m] = clus["cluster_id"]
            member_arrays["mem_parent_biounit"][m] = _biounit_of(mem["pdbpath"])
            for key in ["cg_names", "cg_seg", "cg_chain", "cg_resnum", "cg_resname"]:
                member_arrays[f"mem_{key}"][m] = np.asarray(
                    mem[key], dtype=member_arrays[f"mem_{key}"].dtype)
            segs, chains, resnums, resnames = zip(
                *[(str(s), str(c), int(r), str(n)) for s, c, r, n in mem["scrr"]])
            member_arrays["mem_scrr_seg"][m] = np.array(segs, dtype="U8")
            member_arrays["mem_scrr_chain"][m] = np.array(chains, dtype="U2")
            member_arrays["mem_scrr_resnum"][m] = np.array(resnums, dtype=np.int32)
            member_arrays["mem_scrr_resname"][m] = np.array(resnames, dtype="U4")
            _store_quality(member_arrays, "mem", m, mem["quality"])
            m += 1
    arrays.update(member_arrays)

    # Resolved position by position across the bucket's nr rows, not taken from
    # row 0: a blank PDB element column is common, and the streaming gate
    # (_stream_worker) already admits records whose blanks it could not check.
    # Comparing row 0 verbatim would then reject a later row that merely fills a
    # blank row 0 left empty. Normalized the same way the gate normalizes, so the
    # stored row is comparable to RDKit symbols.
    resolved_elements = [""] * n_cg
    for i, clus in enumerate(clusters):
        nr = clus["nr"]
        arrays["cluster_id"][i] = clus["cluster_id"]
        arrays["cluster_size"][i] = clus["cluster_size"]
        arrays["cluster_pose_radius"][i] = clus["cluster_pose_radius"]
        # Distinct PDB *entries*, not biounit stems: `1abc_1` and `1abc_2` are two
        # assemblies of one deposition, so counting stems would readmit exactly
        # the NCS/homolog inflation this field exists to strip out of
        # cluster_size.
        arrays["cluster_num_parents"][i] = len({_pdb_id_from_path(nr["pdbpath"])}
            | {_pdb_id_from_path(mem["pdbpath"]) for mem in clus["members"]})
        arrays["first_stage_cluster_id"][i] = clus["first_stage_cluster_id"]
        arrays["second_stage_cluster_id"][i] = clus["second_stage_cluster_id"]
        arrays["nr_cg_coords"][i] = np.asarray(nr["cg_coords"], dtype=np.float32)
        arrays["nr_vdm_bb_coords"][i] = np.asarray(nr["vdm_bb_coords"], dtype=np.float32)
        arrays["nr_parent_biounit"][i] = _biounit_of(nr["pdbpath"])
        _store_quality(arrays, "nr", i, nr["quality"])
        for key in ["cg_names", "cg_seg", "cg_chain", "cg_resnum", "cg_resname"]:
            arrays[f"nr_{key}"][i] = np.asarray(nr[key], dtype=arrays[f"nr_{key}"].dtype)
        # One elements row per bucket: assert rather than assume, since a
        # disagreement would mean the automorphism group does not preserve
        # elements and the whole library's atom correspondence is suspect.
        # Only positions where both rows are non-blank can disagree -- a blank is
        # an absent annotation, not the claim "no element here".
        elements = [str(e).strip().capitalize() for e in nr["cg_elements"]]
        conflict = [(k, resolved_elements[k], e) for k, e in enumerate(elements)
                    if e and resolved_elements[k] and resolved_elements[k] != e]
        if conflict:
            raise ValueError(
                f"CG element sequence differs between nr vdGs in bucket "
                f"{'_'.join(reordered_AAs)}: {conflict}. "
                f"CG automorphisms must preserve element.")
        for k, e in enumerate(elements):
            if e and not resolved_elements[k]:
                resolved_elements[k] = e
        arrays["nr_slot_flag"][i] = np.asarray(nr["slot_flags"], dtype=np.int8)
        segs, chains, resnums, resnames = zip(*[(str(s), str(c), int(r), str(n)) for s, c, r, n in nr["scrr"]])
        arrays["nr_scrr_seg"][i] = np.array(segs, dtype="U8")
        arrays["nr_scrr_chain"][i] = np.array(chains, dtype="U2")
        arrays["nr_scrr_resnum"][i] = np.array(resnums, dtype=np.int32)
        arrays["nr_scrr_resname"][i] = np.array(resnames, dtype="U4")

    arrays["cg_elements"][:] = np.asarray(resolved_elements, dtype="U2")

    # Fixed-width assignment truncates silently in numpy. A clipped biounit stem
    # or mirror path makes resolve_parent_pdb_path fail for every row in the
    # bucket, and nothing downstream can tell a clipped value from a real one.
    for _name in ("nr_parent_biounit", "mem_parent_biounit"):
        _arr = arrays.get(_name)
        if _arr is None or _arr.size == 0:
            continue
        _width = _arr.dtype.itemsize // 4
        _clipped = [v for v in np.unique(_arr) if len(v) >= _width]
        if _clipped:
            raise ValueError(
                f"{_name} values reach the {_width}-char dtype width and may be "
                f"truncated (e.g. {_clipped[:3]}); widen the dtype.")
    if len(str(parent_pdb_dir or "")) >= arrays["parent_pdb_dir"].dtype.itemsize // 4:
        raise ValueError(
            f"parent_pdb_dir {parent_pdb_dir!r} reaches the dtype width and may be "
            f"truncated; widen the dtype.")

    # Written via tmp + os.replace, never straight onto the final name: an SGE
    # kill or a quota hit part-way through savez_compressed would otherwise leave
    # a truncated npz that both readers treat as a warning and skip. 
    path = _bucket_npz_path(vdglib_dir, size_subset, reordered_AAs)
    # Deliberately does NOT end in '.npz': a SIGKILL is uncatchable, so the
    # except below cannot run, and a leftover *.npz in the bucket dir is loaded
    # as a real bucket by every reader that walks for '*.npz' -- with an AA label
    # parsed out of the temp name. Written through a file handle because
    # np.savez_compressed appends '.npz' to a *path* that lacks it, but leaves a
    # file object's name alone.
    tmp = f"{path}.{os.getpid()}.tmp"
    try:
        with open(tmp, 'wb') as handle:
            np.savez_compressed(handle, **arrays)
        os.replace(tmp, path)
    except BaseException:
        # BaseException, not Exception: SIGTERM from the queue arrives as
        # SystemExit/KeyboardInterrupt, which is the case this guard exists for.
        if os.path.exists(tmp):
            os.remove(tmp)
        raise

def _abort_profile(profile, profile_json, prefix, reason):
    """Record peak RSS and write the profile sidecar on a fatal path.

    The counters that diagnose an OOM -- rss_peak_child_mb above all -- are
    otherwise recorded only on the success path, so the run whose memory
    footprint you actually need to see is precisely the one that leaves no
    sidecar behind. `aborted_at` marks the file as a crash profile: its phase
    timings are partial and its counters stop at the failure, so it must not be
    read as a completed run.

    Safe on every path: both calls no-op when profiling is disabled. RSS is
    meaningful here only because each caller sits downstream of the executor's
    shutdown, so the dead workers have been reaped and ru_maxrss includes them.
    """
    profile.set("aborted_at", reason)
    profile.record_peak_rss(prefix)
    profile.write(profile_json)


def _run_one_bucket_strict(args):
    """
    Run clustering for each AA-composition bucket in a separate process.
    """
    (bucket_path, seq_sim_thresh, reordered_AAs, cg_automorphisms, vdglib_dir,
     logfile, size_subset, max_num_to_clus) = args

    try:
        # A stale npz and a stale .FAILED from an earlier invocation are both
        # cleared before this bucket does anything else: otherwise a bucket that
        # fails now leaves the previous run's npz beside its new .FAILED (and
        # _subset_output_counts counts it), and a bucket that succeeds now
        # leaves the previous run's .FAILED marking a library that is fine.
        _clear_bucket_outputs(vdglib_dir, size_subset, "_".join(reordered_AAs))
        _vdgs = _load_bucket(bucket_path)
        cg_automorphisms = validate_atom_permutations(cg_automorphisms)
        expected_n_cg = len(cg_automorphisms[0])
        # Dropped records are counted and sampled in the log: a systematic
        # failure (bad --cg-smarts, a truncated stage-1 stage) otherwise shows
        # up only as an empty bucket with nothing explaining it.
        _kept, _dropped_paths = [], []
        for vdg_data in _vdgs:
            if _has_complete_stage1_coords(
                    vdg_data, expected_n_cg, expected_num_vdms=size_subset):
                _kept.append(vdg_data)
            elif len(_dropped_paths) < 5:
                try:
                    _dropped_paths.append(str(vdg_data[4]))
                except (IndexError, TypeError):
                    _dropped_paths.append("<unreadable record>")
        _n_dropped = len(_vdgs) - len(_kept)
        if _n_dropped:
            _log_write(logfile,
                f"\t[WARNING] AA bucket {'_'.join(reordered_AAs)}: dropped "
                f"{_n_dropped}/{len(_vdgs)} vdGs with missing or non-finite "
                f"Stage-1 coords (e.g. {_dropped_paths}).\n")
        _vdgs = _kept

        if not _vdgs:
            return ("_".join(reordered_AAs), 0, {})

        # Cap on distinct PDB IDs (single PDB can still contribute multiple vdGs)
        if max_num_to_clus is not None and len(_vdgs) > max_num_to_clus:
            orig_num = len(_vdgs)
            vdg_pdb_ids = [_pdb_id_from_path(z[4]) for z in _vdgs]
            selected_ids = set(select_diverse_pdbIDs(vdg_pdb_ids, max_num_to_clus))
            _vdgs = [z for z, pid in zip(_vdgs, vdg_pdb_ids) if pid in selected_ids]
            if orig_num != len(_vdgs):
                _log_write(logfile, f"\t{orig_num} vdGs → {len(_vdgs)} vdGs from {len(selected_ids)} diverse PDB IDs for {tuple(reordered_AAs)}.\n")

        reordered_AAs_str = "_".join(reordered_AAs)
        (all_cg_coords, all_vdmbb_coords, all_flankseqs, all_flankCAs, all_pdbpaths,
         all_scrr, all_cg_names, all_cg_elements,
         all_cg_seg, all_cg_chain, all_cg_resnum, all_cg_resname,
         all_slot_flags, all_quality
         ) = unpack_vdg_records(_vdgs)
        if not all_cg_coords:
            _log_write(logfile, f"[WARNING] AA bucket {reordered_AAs_str} has no CG environments; skipping.\n")
            return (reordered_AAs_str, 0, {})
        all_cgvdmbb_coords = combine_cg_and_vdmbb_coords(all_cg_coords, all_vdmbb_coords)
        all_flat_flankseqs = flatten_flanking_seqs(all_flankseqs)
        all_flat_flankCAs = flatten_flanking_CAs(all_flankCAs)

        n_cg_atoms = len(all_cg_coords[0])
        num_cgvdmbb_atoms = n_cg_atoms + len(all_vdmbb_coords[0]) * 3
        cgvdmbb_rmsd_cut = normalize_rmsd(num_cgvdmbb_atoms, "cgvdmbb")
        # The vdG's full symmetry group: CG automorphisms crossed with orderings
        # of interchangeable same-label vdM slots. Both are minimised over inside
        # the Stage-1 distance, so one record is one physical environment.
        perm_group = build_perm_group(cg_automorphisms, n_cg_atoms, reordered_AAs)
        bucket_stats = {"records": len(all_cgvdmbb_coords),
                        "perm_group": len(perm_group)}
        _stage1_t0 = time.time()
        stage1_clusters, _stage1_reps, _stage1_radii = clust.get_butina_clusters(
            all_cgvdmbb_coords, cgvdmbb_rmsd_cut, n_cg_atoms,
            cg_symm_perms=cg_automorphisms, aa_bucket_parts=reordered_AAs,
            counters=bucket_stats)
        bucket_stats["stage1_wall_s"] = time.time() - _stage1_t0
        bucket_stats["stage1_clusters"] = len(stage1_clusters)
        # Rank Stage-1 clusters by size so cluster 1 is the most populated one.
        reassigned_cgvdmbb_clus = {
            new_num: members for new_num, (_old, members) in enumerate(
                sorted(stage1_clusters.items(),
                       key=lambda kv: (-len(kv[1]), kv[0])), start=1)}
        clust.clear_caches()

        # Stage 2 is timed separately from Stage 1 because the two scale on
        # different things and only Stage 1 was instrumented: the profile could
        # show Stage 1's cost but left the rest of clustering as an unattributed
        # residual, so there was no way to tell what share Stage 2 held. That
        # matters now that Stage 2 minimises over same-label slot orderings,
        # which costs ~2.5x on the same-label buckets (~9% of size-2 buckets).
        clusters_out, cluster_counter = [], 0
        _stage2_t0 = time.time()
        _stage2_same_label_s = 0.0
        _n_slot_orders = len(slot_orders(reordered_AAs))
        for cgvdmbb_clusnum, stage1_idxs in reassigned_cgvdmbb_clus.items():
            stage1_idxs = list(stage1_idxs)
            if not stage1_idxs: continue
            if len(stage1_idxs) == 1:
                cluster_counter += 1
                clusters_out.append({
                    "cluster_id": cluster_counter, "cluster_size": 1,
                    "cluster_pose_radius": 0.0,
                    "first_stage_cluster_id": int(cgvdmbb_clusnum), "second_stage_cluster_id": 1,
                    "nr": _extract_env(all_cg_coords, all_vdmbb_coords,
                        all_pdbpaths, all_scrr, all_cg_names, all_cg_elements, all_cg_seg, all_cg_chain,
                        all_cg_resnum, all_cg_resname, all_slot_flags, all_quality, stage1_idxs[0]),
                    "members": _extract_members(all_pdbpaths, all_scrr, all_cg_names,
                        all_cg_seg, all_cg_chain, all_cg_resnum, all_cg_resname,
                        all_quality, stage1_idxs, stage1_idxs[0])})
                continue
            clus_flat_seqs = [all_flat_flankseqs[i] for i in stage1_idxs]
            clus_flat_cas = [all_flat_flankCAs[i] for i in stage1_idxs]
            # Missing flanks must not tighten the geometric criterion: the cutoff
            # is based on the expected flattened size, while pair RMSDs use only
            # coordinate rows that are finite in both structures.
            expected_flank_atoms = len(clus_flat_cas[0])
            flankbb_cut = normalize_rmsd(expected_flank_atoms, "flankbb")
            # get_leader_clusters below mixes RMSD and sequence dissimilarity as
            # seq_dissim * seq_weight + RMSD, with seq_weight=0.5 by default, so
            # the threshold must live on that same combined scale: the flank RMSD
            # cutoff plus half the allowed sequence-dissimilarity budget.
            thresh = flankbb_cut + (1.0 - seq_sim_thresh) / 2.0
            # Same-label slot orderings, the ones Stage 1 already quotients out
            # via build_perm_group. Without them Stage 2 compares slot 1's flank
            # to slot 1's flank positionally and splits a swapped pair into two
            # subgroups, halving the cluster_size/cluster_num_parents of one real
            # mode. Mixed-label buckets yield only the identity and cost nothing.
            _t0 = time.time()
            stage2_assigns = clust.get_leader_clusters(
                zip([clus_flat_seqs, clus_flat_cas], ["flankseq", "flankbb"]), thresh,
                missing_seq_similarity=seq_sim_thresh,
                final_exact_medoid_pass=True, final_reassign_once=True,
                slot_orders=slot_orders(reordered_AAs))
            if _n_slot_orders > 1:
                # Attributed separately so the cost of the permutation minimisation
                # is readable on its own, not buried in the Stage-2 total.
                _stage2_same_label_s += time.time() - _t0
            stage2_before = len(clusters_out)
            for sub_num, local_idxs in stage2_assigns.items():
                if not local_idxs: continue
                cluster_counter += 1
                global_member_idxs = [stage1_idxs[i] for i in local_idxs]
                # Stage 2 partitions on flanking context, so its own medoid is
                # chosen on the wrong metric: the stored row has to be central in
                # *pose*, or the radius recorded next to it means nothing.
                global_cent_idx, subgroup_radius = clust.pose_minimax_prototype(
                    all_cgvdmbb_coords, global_member_idxs, n_cg_atoms, perm_group)
                clusters_out.append({
                    "cluster_id": cluster_counter, "cluster_size": len(local_idxs),
                    "cluster_pose_radius": subgroup_radius,
                    "first_stage_cluster_id": int(cgvdmbb_clusnum), "second_stage_cluster_id": int(sub_num),
                    "nr": _extract_env(all_cg_coords, all_vdmbb_coords,
                        all_pdbpaths, all_scrr, all_cg_names, all_cg_elements, all_cg_seg, all_cg_chain,
                        all_cg_resnum, all_cg_resname, all_slot_flags, all_quality, global_cent_idx),
                    "members": _extract_members(all_pdbpaths, all_scrr, all_cg_names,
                        all_cg_seg, all_cg_chain, all_cg_resnum, all_cg_resname,
                        all_quality, global_member_idxs, global_cent_idx)})
            if len(clusters_out) == stage2_before:
                _log_write(logfile, f"[WARNING] Stage 1 cluster {cgvdmbb_clusnum} in {reordered_AAs_str} lost {len(stage1_idxs)} envs in Stage 2.\n")

        # cluster_counter above is provisional: ids are reassigned by final
        # cluster size so that cluster 1 is always the most populated one.
        bucket_stats["stage2_wall_s"] = time.time() - _stage2_t0
        bucket_stats["stage2_same_label_wall_s"] = _stage2_same_label_s
        bucket_stats["stage2_slot_orders"] = _n_slot_orders

        _renumber_clusters_by_size(clusters_out)

        # Write this AA bucket summary to npz (per-cluster nr vdGs, plus
        # lightweight per-member identity so cluster members can be
        # re-materialized on demand -- see _write_bucket_npz).
        _write_bucket_npz(vdglib_dir, size_subset, reordered_AAs, clusters_out)
        total_vdgs = sum(c["cluster_size"] for c in clusters_out)
        return ("_".join(reordered_AAs), total_vdgs, bucket_stats)

    except Exception as e:
        # One bad bucket must not abort the pool: every nr_vdgs/*.npz already
        # written would stay on disk and the fragment dir would look like a
        # finished library. Drop a .FAILED marker beside the buckets, keep
        # going, and let main() exit non-zero at the end.
        aa_label = "_".join(reordered_AAs) if reordered_AAs else "UNKNOWN_AAs"
        tb = traceback.format_exc()
        err_text = (f"[WORKER ERROR] AA bucket: {aa_label} "
                    f"(subset size {size_subset})\n"
                    f"Exception: {e}\nTraceback:\n{tb}")
        _log_write(logfile, err_text + "\n")
        _write_failed_marker(vdglib_dir, size_subset, aa_label, err_text)
        return (aa_label, None, None)

def _parse_pdb_with_retry(pdb_file, attempts=3, delay=0.5):
    """pr.parsePDB with retries -- worker chunks are file-level while biounits are
    per-environment, so two workers occasionally parse the same PDB concurrently;
    this absorbs the transient read failures that causes on Wynton's NFS.

    Mechanism lives in vdg_npz_utils so the read path here and the one the
    library's readers use cannot drift apart."""
    return parse_pdb_with_retry(pdb_file, attempts=attempts, delay=delay)


# Errnos that mean "the filesystem was busy or the connection blipped", i.e. the
# transient class _parse_pdb_with_retry exists for. A failure carrying one of
# these is never cached as permanent: doing so would convert a recoverable blip
# into a skip of every remaining environment in that biounit, which is exactly
# what the retry was written to prevent.
_TRANSIENT_ERRNOS = frozenset(
    e for e in (errno.EAGAIN, errno.EBUSY, errno.EINTR, errno.EIO, errno.EMFILE,
                errno.ENFILE, errno.ENOMEM, errno.ESTALE, errno.ETIMEDOUT,
                getattr(errno, "EREMOTEIO", None))
    if e is not None)

# Cached in place of a parse that failed for a reason retrying cannot fix (the
# file is missing, or it opened cleanly and ProDy could not read it). Without
# this, every environment referencing a genuinely corrupt biounit re-runs the
# full retry/backoff sequence -- seconds of sleeping per environment, repeated
# for as many environments as that biounit contributed.
_PARSE_FAILED = object()


def _is_transient_parse_failure(exc):
    """Whether a parse failure is worth retrying on a later environment.

    OSError is split on errno: ENOENT/EACCES describe the file itself and will
    not change, while the errnos above describe the filesystem. Anything that is
    not an OSError got far enough to open the file and fail on its *contents*
    (ProDy's own parse errors), which is permanent too.
    """
    if isinstance(exc, OSError):
        return exc.errno in _TRANSIENT_ERRNOS
    return False


def _get_atomgroup_for_env(
    environment, pdb_dir, cg, cg_match_dict, align_atoms, logfile,
    pdb_cache=None):
    biounit = environment[0][0]
    middle_two = biounit[1:3].lower()
    pdb_file = os.path.join(pdb_dir, middle_two, biounit + ".pdb")

    if pdb_cache is not None and pdb_file in pdb_cache:
        whole_struct = pdb_cache[pdb_file]
        if whole_struct is _PARSE_FAILED:
            return None
    else:
        try:
            whole_struct = _parse_pdb_with_retry(pdb_file)
        except Exception as e:
            if pdb_cache is not None and not _is_transient_parse_failure(e):
                pdb_cache[pdb_file] = _PARSE_FAILED
            _log_warn_capped(logfile, 'pdb_parse_failed',
                f"[WARNING] Failed to parse PDB {biounit}: {e}; skipping.\n")
            return None
        if pdb_cache is not None:
            pdb_cache[pdb_file] = whole_struct

    if whole_struct is None:
        _log_warn_capped(logfile, 'pdb_returned_none',
            f"[WARNING] PDB {biounit} returned None; skipping.\n")
        return None

    # Two forms of the resnum are needed and must not be confused: ProDy selection
    # strings require negative resnums backquoted, but cg_match_dict is keyed on the
    # raw PDB resnum column ("-5", no backticks -- see find_cg_matches in
    # external/vdG-miner/vdg_miner/vdg/cg.py). Using the backquoted form in the key
    # silently drops every negative-resnum ligand as "no CG match".
    raw_resnums = [tup[3] for tup in environment]
    scrs = [(tup[1], tup[2], f"`{tup[3]}`" if tup[3] < 0 else tup[3]) for tup in environment]
    selstrs = [
        f"(segment {scr[0]} and chain {scr[1]} and resnum {scr[2]})" if scr[0]
        else f"(chain {scr[1]} and resnum {scr[2]})"
        for scr in scrs]
    if len(selstrs) < 2:
        _log_warn_capped(logfile, 'env_too_few_entries',
            f"[WARNING] {biounit}: environment has <2 entries; skipping.\n")
        return None
    try:
        sel = whole_struct.select(
            "same residue as within 5 of ({})".format(" or ".join(selstrs[1:])))
    except Exception:
        _log_warn_capped(logfile, 'prody_selection_failed',
            f"[WARNING] {biounit}: ProDy selection failed; skipping.\n")
        return None
    if sel is None:
        return None
    struct = sel.toAtomGroup()
    resnames = []
    align_coords = np.zeros((3, 3))
    cg_atom_coords = None  # set when i == 0; read by the vdM branch below

    for i, (scr, selstr) in enumerate(zip(scrs, selstrs)):
        try:
            substruct = struct.select(selstr)
            if substruct is None:  # selection empty; skip env
                return None
            resnames.append(substruct.getResnames()[0])

            unique_res_indices = np.unique(substruct.getResindices())
            if len(unique_res_indices) != 1:
                _log_warn_capped(logfile, 'ambiguous_residue',
                    f"[WARNING] {biounit} chain {scr[1]} resnum {scr[2]}: ambiguous residue; skipping.\n")
                return None

            if i == 0:
                if cg in cg_atoms.keys():
                    atom_names_list = cg_atoms[cg][resnames[0]]
                else:
                    key = (biounit, scrs[0][0], scrs[0][1],
                           str(raw_resnums[0]), resnames[0])

                    match_list = cg_match_dict.get(key)
                    match_idx = environment[0][4] - 1  # 1-based index in env --> 0-based

                    if match_list is None:  # no CG match; possibly missing density; skip
                        _log_warn_capped(logfile, 'no_cg_match',
                            f"[WARNING] {biounit} chain {scr[1]} resnum {scr[2]}: "
                            f"no CG match for key {key} in cg_match_dict; skipping.\n")
                        return None

                    if not (0 <= match_idx < len(match_list)):  # out of range; possibly
                                                                # obabel issue; skip
                        _log_warn_capped(logfile, 'match_idx_out_of_range',
                            f"[WARNING] {biounit} chain {scr[1]} resnum {scr[2]}: "
                            f"match_idx {match_idx} out of range for match_list of "
                            f"length {len(match_list)}; skipping.\n")
                        return None
                    atom_names_list = match_list[match_idx]

                cg_atom_selstrs = ["name " + name_selstr(atom_name)
                                   for atom_name in atom_names_list]

                # An atom name must identify exactly one atom in the residue. ProDy's
                # default parse keeps only altloc 'A' and blank, so a name matching
                # two atoms here is not an altloc pair -- the residue is malformed.
                #
                # Origin: the parent PDB has an amino acid modelled with its backbone
                # N alone. Prepwizard cannot build an amino acid from an N with no CA,
                # so it re-emits that orphan N as an ammonium ion carrying a *ligand's*
                # resname/chain/resnum. 
                #
                # Resolving that by proximity is a guess, and guessing wrong writes a
                # vdG whose coordinates contradict its own recorded identity. 
                chosen_atoms = []
                for cg_selstr, atom_name in zip(cg_atom_selstrs, atom_names_list):
                    atom_sel = substruct.select(cg_selstr)
                    if atom_sel is None or atom_sel.numAtoms() == 0:
                        return None
                    if atom_sel.numAtoms() > 1:
                        _log_warn_capped(logfile, 'duplicate_cg_atom_name',
                            f"[WARNING] {biounit} (chain {scrs[0][1]}, resnum "
                            f"{scrs[0][2]}, {resnames[0]}): {atom_sel.numAtoms()} atoms "
                            f"named {atom_name} in one residue; malformed parent "
                            f"residue, skipping environment.\n")
                        return None
                    chosen_atoms.append(atom_sel[0])

                cg_atom_coords = np.asarray(
                    [np.reshape(a.getCoords(), 3) for a in chosen_atoms], dtype=float)

                # The atom *set* comes from OpenBabel's perception of the parent PDB
                # (cg_match_dict), and nothing upstream checks it against geometry.
                # A SMARTS match is a connected subgraph by construction, so its atoms
                # must form one bonded component in space too; when perception invents
                # a bond (disordered/clashing ligand, two copies overlaid) the match
                # can span atoms that are nowhere near each other. Element identity is
                # checked separately, per record, in _stream_one_chunk.
                components = _cg_bond_components(cg_atom_coords)
                if len(components) > 1:
                    groups = ' | '.join(
                        '+'.join(atom_names_list[k] for k in comp)
                        for comp in components)
                    gap = min(
                        float(np.linalg.norm(cg_atom_coords[k] - cg_atom_coords[l]))
                        for a, b in itertools.combinations(components, 2)
                        for k in a for l in b)
                    _log_warn_capped(logfile, 'cg_not_bonded_component',
                        f"[WARNING] {biounit} (chain {scrs[0][1]}, resnum "
                        f"{scrs[0][2]}, {resnames[0]}): CG atoms {atom_names_list} "
                        f"split into {len(components)} bonded components ({groups}); "
                        f"closest approach between components {gap:.2f} A > "
                        f"{MAX_CG_BOND_DIST} A; implausible perception, skipping.\n")
                    return None

                for j, chosen_atom in enumerate(chosen_atoms):
                    # Set the CG-encoding occupancy for this chosen atom
                    chosen_atom.setOccupancy(cg_slot_occupancy(j))

                    # Fill alignment coordinates for the chosen CG atoms
                    if j in align_atoms:
                        c = np.asarray(chosen_atom.getCoords())
                        c = c[0] if c.ndim == 2 else c
                        align_coords[align_atoms.index(j)] = c

            else:
                # A vdM must contact the CG on *heavy* atoms, not merely be a probe
                # neighbour of it.
                _names = substruct.getNames()
                _els = substruct.getElements()
                if _els is None:  # no element column at all; names carry it
                    _els = [''] * len(_names)
                _heavy_mask = np.array(
                    [not is_hydrogen(n, e) for n, e in zip(_names, _els)], dtype=bool)
                if not _heavy_mask.any():
                    # Same reasoning as the contact cutoff below: drop the residue,
                    # not the environment. 
                    _log_warn_capped(logfile, 'vdm_all_hydrogen',
                        f"[WARNING] {biounit} (chain {scr[1]}, resnum {scr[2]}): "
                        f"no heavy atoms in this residue; dropping it from the "
                        f"environment.\n")
                    continue
                vdm_heavy_coords = np.asarray(substruct.getCoords())[_heavy_mask]
                # i == 0 sets cg_atom_coords or returns, so it is always set here; a
                # None would silently disable the contact filter for every vdM slot.
                assert cg_atom_coords is not None, "CG coords missing at the vdM branch"
                d = np.sqrt(((cg_atom_coords[:, None, :]
                              - vdm_heavy_coords[None, :, :]) ** 2).sum(-1)).min()
                if d > CG_VDM_CONTACT_CUTOFF:
                    _log_warn_capped(logfile, 'vdm_not_in_contact',
                        f"[WARNING] {biounit} (chain {scr[1]}, resnum {scr[2]}): "
                        f"nearest heavy atom is {d:.1f} A from the CG (cutoff "
                        f"{CG_VDM_CONTACT_CUTOFF} A); not a CG contact, "
                        f"dropping this residue from the environment.\n")
                    continue
                # Non-CG residues: mark them differently
                substruct.setOccupancies(VDM_OCC)

        except Exception:
            tb = traceback.format_exc()
            _log_warn_capped(logfile, 'residue_selection_exception',
                f"[WARNING] {biounit} (chain {scr[1]}, resnum {scr[2]}) skipped: "
                f"exception during residue selection:\n{tb}\n")
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

_FLUSH_RECORDS_THRESHOLD = 500_000

# Total parsed atoms a worker may hold in its PDB cache. ~40 M atoms is roughly
# 1 GB of ProDy coordinate+label arrays; divided across workers in main().
_PDB_CACHE_ATOM_BUDGET = 40_000_000


def _struct_atom_count(struct):
    """Cache weight of a parsed structure; sentinels (parse failures) weigh 0."""
    try:
        return int(struct.numAtoms())
    except Exception:
        return 0


def _flush_buckets_to_disk(bucket_records, worker_tmp_dir, flush_idx):
    """Write bucket records to disk with flush-index suffix."""
    for aa_key, records in bucket_records.items():
        fname = _bucket_fname(aa_key, flush_idx=flush_idx)
        _dump_pickle_atomic(records, os.path.join(worker_tmp_dir, fname))


def _stream_one_chunk(args):
    (environment_files_chunk, environments_dir, pdb_dir, CG, cg_match_dict_pkl,
     align_atoms, logfile, expected_n_cg, expected_element_seq,
     num_flanking, worker_tmp_root, pdb_cache_limits, subset_sizes,
     flush_threshold) = args

    if cg_match_dict_pkl:
        with open(cg_match_dict_pkl, "rb") as _fh:
            cg_match_dict = pickle.load(_fh)
    else:
        cg_match_dict = {}
    pdb_cache_size, pdb_cache_atom_budget = pdb_cache_limits
    pdb_cache = (_LRUCache(maxsize=pdb_cache_size,
                           max_weight=pdb_cache_atom_budget,
                           weigh=_struct_atom_count)
                 if pdb_cache_size > 0 else None)
    bucket_records = {size: {} for size in subset_sizes}
    worker_tmp_dirs = {
        size: os.path.join(worker_tmp_root, str(size)) for size in subset_sizes}
    for worker_tmp_dir in worker_tmp_dirs.values():
        os.makedirs(worker_tmp_dir, exist_ok=True)
    total_in_memory, flush_idx = 0, 0
    # Every silent drop below is tallied so a systematic --cg-smarts/--cg-match-dict
    # mismatch shows up as a reason in the log instead of an empty library.
    skips = {"no_atomgroup": 0, "dup_occupancy": 0, "no_cg_atoms": 0,
             "bad_cg_coords": 0, "cg_elements_mismatch": 0, "unparsable_line": 0}
    # Deltas, not absolutes: _WARN_COUNTS is module state and the executor reuses
    # a worker across chunks, so the totals at exit include earlier chunks' counts.
    warn_counts_at_entry = dict(_WARN_COUNTS)

    for subdir, file in environment_files_chunk:
        env_path = os.path.join(environments_dir, subdir, file)
        if not os.path.isfile(env_path):
            continue

        with open(env_path, "r") as f_env:
            for line in f_env:
                try:
                    record = json.loads(line)
                except Exception:
                    skips["unparsable_line"] += 1
                    continue

                # generate_environments.py emits one dict per environment: the
                # residue tuples plus the quality actually measured on the atoms
                # that enter the vdG.
                environment = record["env"]
                env_quality = (record["cg_max_b"], record["cg_min_occ"],
                               record["vdm_max_b"], record["vdm_min_occ"])

                pdb_label = "_".join([str(el) for el in environment[0]])
                _biounit = environment[0][0]
                # Format consumed by _pdb_id_from_path -- keep the two in sync.
                _actual_pdb_path = os.path.join(
                    pdb_dir, _biounit[1:3].lower(), _biounit + ".pdb")

                atomgroup = _get_atomgroup_for_env(
                    environment, pdb_dir, CG, cg_match_dict, align_atoms, logfile,
                    pdb_cache=pdb_cache)
                if atomgroup is None:
                    skips["no_atomgroup"] += 1
                    continue

                atomgroup = _resolve_duplicate_ligand_occupancies(atomgroup, pdb_label)
                if atomgroup is None:
                    skips["dup_occupancy"] += 1
                    continue

                _cg_all = get_cg_atoms(atomgroup, pdb_label)
                if _cg_all is None:
                    skips["no_cg_atoms"] += 1
                    continue
                # The last four are per-vdG scalars, not per-atom lists: a CG
                # lies within one ligand residue (enforced by get_cg_atoms).
                (cg_coords, cg_names, cg_elements,
                 cg_seg, cg_chain, cg_resnum, cg_resname) = _cg_all

                if (cg_coords.shape != (expected_n_cg, 3)
                        or not np.isfinite(cg_coords).all()):
                    skips["bad_cg_coords"] += 1
                    continue

                # The CG atom set comes from OpenBabel's perception of the parent PDB
                # (cg_match_dict); nothing else checks it against the SMARTS. Validate
                # per record and drop the record, rather than letting the disagreement
                # surface as a raise at npz-write time.
                #
                # Compared position by position, not as a multiset: the record is
                # stored in SMARTS-slot order and every downstream comparison indexes
                # the automorphism group by that order, so the right elements in the
                # wrong order is a real defect, not a relabeling. Every automorphism
                # preserves atomic number (asserted in main), so no permutation of the
                # group can change the expected sequence.
                #
                # A blank PDB element column is common, so blanks are treated as
                # wildcards here instead of failing every record in such a file.
                # Slots the SMARTS does not pin (None) are unconstrained; there is
                # deliberately no fallback reference taken from the first record of
                # the chunk, which would make acceptance depend on --num-procs.
                _elements = tuple(str(e).strip().capitalize() for e in cg_elements)
                _bad_elements = any(
                    a and b and a != b
                    for a, b in zip(_elements, expected_element_seq))
                if _bad_elements:
                    skips["cg_elements_mismatch"] += 1
                    if skips["cg_elements_mismatch"] <= 10:
                        _log_write(logfile,
                            f"[WARNING] {pdb_label}: CG elements {list(_elements)} "
                            f"disagree with the SMARTS elements "
                            f"{list(expected_element_seq)}; skipping.\n")
                    continue

                vdms_dict = clust.get_vdm_res_features(atomgroup, pdb_label, num_flanking)
                vdm_resinds = list(vdms_dict.keys())

                # is_hydrogen, not a bare element test: the element column is
                # blank in older and hand-edited PDBs, where every CG hydrogen
                # would otherwise count as heavy and skew the bb-vs-sc slot
                # labeling this set drives.
                _heavy_mask = np.array(
                    [not is_hydrogen(n, e) for n, e in zip(cg_names, cg_elements)],
                    dtype=bool)
                if _heavy_mask.any():
                    _cg_heavy = cg_coords[_heavy_mask]
                else:
                    _log_warn_capped(logfile, 'cg_all_hydrogen_fallback',
                        f"[WARNING] {pdb_label}: CG has no heavy atoms (all H/D); using full cg_coords as fallback.\n")
                    _cg_heavy = cg_coords
                for size_subset in subset_sizes:
                    vdg_subsets = get_vdg_subsets_target_size(
                        vdm_resinds, size_subset)
                    for vdg_subset in vdg_subsets:
                        try:
                            (re_ordered_aas, re_ordered_bbcoords,
                             re_ordered_flankingseqs, re_ordered_CAs,
                             re_ordered_scrr, re_ordered_slot_flags
                             ) = clust.reorder_vdg_subset(
                                vdg_subset, vdms_dict, _cg_heavy, atomgroup)
                        except ValueError as _e:
                            _log_write(logfile,
                                f'[WARNING] ({pdb_label}) {_e}; skipping vdG subset.\n')
                            continue
                        rec = {
                            # Preserve SMARTS-slot order; comparisons enumerate
                            # the CG automorphism group.
                            "cg_coords": np.asarray(cg_coords, dtype=np.float32),
                            "bbcoords": [np.asarray(x, dtype=np.float32) for x in re_ordered_bbcoords],
                            "flankseqs": re_ordered_flankingseqs,
                            "flankCAs": [np.asarray(x, dtype=np.float32) for x in re_ordered_CAs],
                            "pdbpath": _actual_pdb_path,
                            "scrr": [[str(s), str(ch), int(r), str(rn)] for (s, ch, r, rn) in re_ordered_scrr],
                            "cg_names": list(cg_names), "cg_elements": list(cg_elements),
                            "cg_seg": str(cg_seg), "cg_chain": str(cg_chain),
                            "cg_resnum": int(cg_resnum), "cg_resname": str(cg_resname),
                            "slot_flags": [int(f) for f in re_ordered_slot_flags],
                            "cg_max_b": env_quality[0], "cg_min_occ": env_quality[1],
                            "vdm_max_b": env_quality[2], "vdm_min_occ": env_quality[3],
                        }
                        aa_key = "_".join(re_ordered_aas)
                        size_records = bucket_records[size_subset]
                        size_records.setdefault(aa_key, []).append(rec)
                        total_in_memory += 1
                        if total_in_memory >= flush_threshold:
                            for size, records in bucket_records.items():
                                _flush_buckets_to_disk(
                                    records, worker_tmp_dirs[size], flush_idx)
                            bucket_records = {size: {} for size in subset_sizes}
                            flush_idx += 1
                            total_in_memory = 0

    # Write any remaining records as the final pickle for each bucket.
    for size, size_records in bucket_records.items():
        for aa_key, records in size_records.items():
            _dump_pickle_atomic(
                records, os.path.join(worker_tmp_dirs[size], _bucket_fname(aa_key)))

    warns = {key: count - warn_counts_at_entry.get(key, 0)
             for key, count in _WARN_COUNTS.items()
             if count - warn_counts_at_entry.get(key, 0) > 0}
    return worker_tmp_root, skips, warns


def _merge_worker_dirs(worker_dirs, final_dir):
    """Merge per-worker pickle files by aa_key.

    Downstream clustering (`_run_one_bucket_strict`) needs one complete file per
    aa_key, so buckets can't be re-split across files the way workers flush them.
    Instead, merge one aa_key at a time: peak memory is bounded by the largest
    single bucket rather than the sum of every bucket in the library.
    """
    os.makedirs(final_dir, exist_ok=True)
    files_by_key = {}  # aa_key -> list of source paths
    for d in worker_dirs:
        if not os.path.isdir(d):
            continue
        for fname in os.listdir(d):
            if not fname.endswith(".pkl"):
                continue
            aa_key = _aa_key_from_bucket_fname(fname)
            if aa_key is None:
                print(f"[WARNING] Ignoring unrecognized file in worker dir: "
                      f"{os.path.join(d, fname)}", file=sys.stderr)
                continue
            files_by_key.setdefault(aa_key, []).append(os.path.join(d, fname))

    for aa_key, paths in files_by_key.items():
        records = []
        for path in paths:
            try:
                with open(path, "rb") as f:
                    records.extend(pickle.load(f))
            except Exception as exc:
                # Shards are written through _dump_pickle_atomic, so a shard that
                # exists is a shard that was written whole; reaching here means
                # real corruption (bad block, truncated NFS write) rather than an
                # interrupted worker. Raise rather than skip. 
                raise RuntimeError(
                    f"Unreadable stream shard {path} for AA bucket {aa_key}: "
                    f"{type(exc).__name__}: {exc}") from exc
        _dump_pickle_atomic(records, os.path.join(final_dir, _bucket_fname(aa_key)))


def _bucket_jobs(stream_dirs, subset_sizes, seq_sim_thresh, cg_automorphisms,
                 vdglib_dir, logfile, max_num_to_clus):
    """Build clustering jobs for every requested subset size, largest first.

    Bucket cost is superlinear in record count and the distribution is very
    lopsided (`bb_bb` and `ARG` run ~500x the median). Queued in name order, one
    of those can be handed out last and set the wall time on its own; longest-
    processing-time-first puts them in while the pool is still empty. Pickle
    size is the proxy for record count -- exact, since every record is the same
    fixed-width tuple.
    """
    jobs = []
    for size_subset in subset_sizes:
        stream_dir = stream_dirs[size_subset]
        if not os.path.isdir(stream_dir):
            _log_write(logfile,
                f"\t[STREAM ERROR] No buckets found at {stream_dir}\n")
            continue
        for fname in sorted(os.listdir(stream_dir)):
            if not fname.endswith(".pkl"):
                continue
            aa_key = _aa_key_from_bucket_fname(fname)
            if aa_key is None:
                _log_write(logfile,
                    f"\t[STREAM WARNING] Ignoring unrecognized bucket file {fname}\n")
                continue
            path = os.path.join(stream_dir, fname)
            # Size captured here so the sort below does not re-stat every bucket
            # O(n log n) times on NFS.
            jobs.append((os.path.getsize(path), (
                path, seq_sim_thresh,
                aa_key.split("_"), cg_automorphisms, vdglib_dir, logfile,
                size_subset, max_num_to_clus)))
    jobs.sort(key=lambda sized_job: sized_job[0], reverse=True)
    return [job for _, job in jobs]


def _subset_output_counts(vdglib_dir, size_subset, logfile):
    """Return (cluster count, input-vdG count, unreadable AA labels) for one
    subset directory.

    An unreadable npz is returned, not just logged: every reader of the library
    skips it with a warning, so the bucket is simply absent from the library, and
    a run that exits 0 over it hands back something that looks complete. The
    caller turns the labels into a non-zero exit.
    """
    nr_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset))
    num_nr_vdgs = 0
    num_inputs = 0
    unreadable = []
    if not os.path.isdir(nr_dir):
        return num_nr_vdgs, num_inputs, unreadable
    for fname in sorted(os.listdir(nr_dir)):
        if not fname.endswith(".npz"):
            continue
        npz_path = os.path.join(nr_dir, fname)
        try:
            with np.load(npz_path) as data:
                num_nr_vdgs += int(data["cluster_id"].shape[0])
                num_inputs += int(data["cluster_size"].sum())
        except Exception as exc:
            _log_write(logfile,
                f"\t[ERROR] Unreadable npz {npz_path} ({type(exc).__name__}: "
                f"{exc}); this AA bucket is missing from the library.\n")
            unreadable.append(fname[:-len(".npz")])
    return num_nr_vdgs, num_inputs, unreadable


def _scan_failed_markers(vdglib_dir, size_subset):
    """AA labels carrying a .FAILED marker in one subset directory.

    Read at the end of the run rather than trusting this run's in-memory tally:
    a marker is the only trace left by a bucket whose process is gone, and one
    left by an earlier run over the same directory (a bucket this run never
    re-ran, so _clear_bucket_outputs never touched it) condemns the library just
    as much as one written a minute ago.
    """
    nr_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset))
    if not os.path.isdir(nr_dir):
        return []
    return sorted(fname[:-len(".FAILED")] for fname in os.listdir(nr_dir)
                  if fname.endswith(".FAILED"))


def _preexisting_bucket_outputs(vdglib_dir, subset_sizes):
    """(size, filename) for bucket outputs already present in the requested
    subset directories."""
    found = []
    for size in subset_sizes:
        nr_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size))
        if not os.path.isdir(nr_dir):
            continue
        found.extend((size, fname) for fname in sorted(os.listdir(nr_dir))
                     if fname.endswith(".npz") or fname.endswith(".FAILED"))
    return found


def main():
    start_time = time.time()
    args = parse_args()
    profile = ComputeProfile(enabled=bool(args.profile_json), cg=args.cg,
                             subset_sizes=sorted(set(args.subset_sizes)))
    CG = args.cg
    seq_sim_thresh, num_flanking = args.seq, args.flank
    symmetry_mol = mol_from_fragment(args.cg_smarts)
    if symmetry_mol is None:
        raise ValueError(
            f"Could not parse --cg-smarts as SMARTS: {args.cg_smarts!r}")
    cg_automorphisms = identify_mol_automorphisms(symmetry_mol)
    # Validate automorphisms early, before streaming
    try:
        cg_automorphisms = validate_atom_permutations(cg_automorphisms)
    except ValueError as exc:
        raise ValueError(f"[ERROR] Invalid CG automorphisms: {exc}")
    if len(cg_automorphisms[0]) < 3:
        sys.exit(
            f"[ERROR] CG '{CG}' has {len(cg_automorphisms[0])} atom(s); the "
            "local alignment frame built in _get_atomgroup_for_env requires "
            "at least 3 CG atoms (indexed by align_atoms below), or every "
            "environment is silently dropped with no diagnostic.")
    # Element sequence, in SMARTS atom order, that the mined CG must reproduce
    # position by position. A multiset would accept a CG whose atoms are the right
    # elements in the wrong order, and order is load-bearing: every record is stored
    # in SMARTS-slot order and compared through the automorphism group, which is
    # indexed by that order. Per position, not all-or-nothing: an OR/negation query
    # atom nulls out its own slot (None from cg_element_symbols) and leaves every
    # pinned slot checked, since a generic atom means the element there genuinely
    # varies and says nothing about its neighbours.
    expected_element_seq = cg_element_symbols(args.cg_smarts)
    # Comparing a mined CG to the SMARTS position by position is only equivalent to
    # comparing it up to symmetry because every automorphism preserves atomic number
    # (identify_mol_automorphisms normalizes charge and H within resonance groups,
    # never element). Assert it rather than assume it: if that ever loosens, the
    # element check silently becomes stricter than the group it is meant to match.
    # A perm mapping a pinned slot onto a generic one trips the same comparison
    # (None != 'C'), which is correct: symmetry-equivalent matches could then be
    # stored in orders the pinned slot's check disagrees about.
    for _perm in cg_automorphisms:
        if any(expected_element_seq[i] != expected_element_seq[j]
               for i, j in enumerate(_perm)):
            sys.exit(
                f"[ERROR] CG automorphism {_perm} does not preserve the element "
                f"sequence {expected_element_seq} of --cg-smarts (None = element "
                "not pinned by the pattern). The stored atom order and the "
                "symmetry group would disagree.")
    subset_sizes = tuple(dict.fromkeys(args.subset_sizes))
    max_num_to_clus, logfile = args.max_num_vdgs_to_clus, args.logfile
    pdb_dir, environments_root = args.pdb_dir, args.environments_dir
    cg_match_dict_pkl, vdglib_dir = args.cg_match_dict_pkl, args.vdglib_dir
    pdb_cache_size = args.pdb_cache_size

    # CG is non-proteinaceous (not in vdG-miner's cg_atoms table) and no match
    # dict was given: every environment's `cg_match_dict.get(key)` lookup in
    # _get_atomgroup_for_env would return None, silently dropping the whole run
    # (streams to 0 vdGs with no error). Fail fast instead of burning compute
    # on a run that can only produce an empty library.
    if CG not in cg_atoms.keys() and not cg_match_dict_pkl:
        sys.exit(
            f"[ERROR] CG '{CG}' is not in vdG-miner's proteinaceous cg_atoms "
            "table and --cg-match-dict-pkl was not provided. Every environment "
            "would be silently skipped. Pass --cg-match-dict-pkl (produced by "
            "smarts_to_cgs.py).")

    _preexisting = _preexisting_bucket_outputs(vdglib_dir, subset_sizes)
    if _preexisting:
        _examples = ", ".join(f"nr_vdgs/{size}/{fname}"
                              for size, fname in _preexisting[:3])
        sys.exit(
            f"[ERROR] {vdglib_dir} already contains bucket outputs for the "
            f"requested subset size(s) ({len(_preexisting)} file(s), e.g. "
            f"{_examples}). Buckets this run does not reproduce would be left "
            "behind and read as part of the library. Remove nr_vdgs/ for these "
            "sizes (or the whole fragment directory) and re-run.")

    # CG atom indices used to build the local alignment frame in
    # _get_atomgroup_for_env; the order is arbitrary (re-aligned in clustering)
    # but requires >=3 CG atoms, enforced above.
    align_atoms = [1, 0, 2]

    environments_dir = os.path.join(environments_root, "environments")
    if not os.path.isdir(environments_dir):
        sys.exit(f"[ERROR] Missing environments dir: {environments_dir}")

    all_environment_files = []
    for subdir in sorted(os.listdir(environments_dir)):
        full = os.path.join(environments_dir, subdir)
        if not os.path.isdir(full):
            continue
        for file in sorted(os.listdir(full)):
            if file.endswith(".jsonl"):
                all_environment_files.append((subdir, file))
    # No shards means generate_environments ran (the wrapper uses check=True, so a
    # crash would already have failed) and produced nothing. Same semantic as zero
    # streamed records below -- this fragment has no vdGs -- so it takes the same
    # exit: a warning, no 'Job completed.', but not a job failure.
    if not all_environment_files:
        err_text = (f"[WARNING] No .jsonl environment shards under {environments_dir}; "
                    "generate_environments.py produced no environments for this CG. "
                    "nr_vdgs/ is EMPTY and this fragment is NOT complete, so this log "
                    "deliberately carries no completion marker.\n")
        _log_write(logfile, err_text)
        print(err_text, file=sys.stderr)
        sys.exit(EXIT_NO_VDGS)

    # Record the SMARTS and the group beside the buckets. The hit finder reads
    # these instead of re-deriving them from the library directory name, which is
    # only the -c label and need not be the SMARTS. Written only after every
    # validation above passes, so a rejected run leaves no cg_symmetry.npz behind
    # in an otherwise-empty fragment directory.
    write_cg_symmetry(args.vdglib_dir, args.cg_smarts, cg_automorphisms)

    # Each environment is decoded and reconstructed once. Records derived for
    # every requested size are written to separate AA-composition directories.
    stream_dirs = {
        size: _aa_tmp_dir(vdglib_dir, size) for size in subset_sizes}
    stream_root = _stream_root(vdglib_dir)
    workers_root = os.path.join(stream_root, "workers")

    n_workers = min(len(all_environment_files), int(args.num_procs))
    if n_workers > 0 and all_environment_files:
        n_chunks = min(len(all_environment_files),
                       n_workers * STREAM_CHUNKS_PER_WORKER)
        chunk_size = max(1, (len(all_environment_files) + n_chunks - 1) // n_chunks)
        chunks = [all_environment_files[i:i + chunk_size]
                  for i in range(0, len(all_environment_files), chunk_size)]
        worker_dirs = [
            os.path.join(workers_root, f"worker_{i}")
            for i in range(len(chunks))]

        # Both budgets below are per-worker but the run pays for them n_workers
        # times over, so they are divided by the worker count rather than being
        # fixed per process (20 slots x 500k records blows past the documented
        # 24.4 G/job ceiling on its own).
        flush_threshold = max(25_000, _FLUSH_RECORDS_THRESHOLD // n_workers)
        pdb_cache_atom_budget = max(
            1, _PDB_CACHE_ATOM_BUDGET // n_workers) if pdb_cache_size > 0 else 0
        profile.set("stream_flush_threshold", flush_threshold)
        chunk_args = [
            (chunk, environments_dir, pdb_dir, CG, cg_match_dict_pkl,
             align_atoms, logfile, len(cg_automorphisms[0]), expected_element_seq,
             num_flanking, wdir, (pdb_cache_size, pdb_cache_atom_budget),
             subset_sizes, flush_threshold)
            for chunk, wdir in zip(chunks, worker_dirs)]

        profile.set("environment_shards", len(all_environment_files))
        profile.set("stream_workers", n_workers)
        profile.set("stream_chunks", len(chunks))
        stream_ctx = mp.get_context("spawn")
        stream_skips = {}
        # Kept out of stream_skips on purpose: these count warnings, not skipped
        # environments, and the EXIT_NO_VDGS diagnosis below reasons about
        # stream_skips as a complete account of why records were lost.
        stream_warns = {}
        try:
          with profile.phase("stream_environments"):
            # ProcessPoolExecutor, not mp.Pool: if the OOM killer takes a worker,
            # mp.Pool silently respawns it and never reports the lost task, so
            # imap_unordered blocks forever and the job burns its whole h_rt
            # before dying with a partial library (CPython bpo-22393). The
            # executor marks the pool broken and raises instead. See the
            # max_tasks_per_child note on the clustering pool below.
            pool = concurrent.futures.ProcessPoolExecutor(
                max_workers=n_workers, mp_context=stream_ctx)
            try:
                futures = [pool.submit(_stream_one_chunk, ca) for ca in chunk_args]
                for fut in concurrent.futures.as_completed(futures):
                    _, chunk_skips, chunk_warns = fut.result()
                    for reason, count in chunk_skips.items():
                        stream_skips[reason] = stream_skips.get(reason, 0) + count
                    for reason, count in chunk_warns.items():
                        stream_warns[reason] = stream_warns.get(reason, 0) + count
            finally:
                # cancel_futures matters: the executor's __exit__ would otherwise
                # drain every queued chunk before the exception surfaced, so a
                # systematic failure would still burn the whole h_rt. mp.Pool's
                # __exit__ called terminate() and stopped immediately; this is the
                # nearest equivalent (in-flight tasks still finish).
                pool.shutdown(wait=True, cancel_futures=True)
        except BrokenProcessPool as e:
            err_text = (
                "[FATAL ERROR] A streaming worker died without returning a result "
                f"({e}). The usual cause is the OOM killer: each worker holds a "
                "parsed biounit plus up to its flush threshold of buffered records. "
                "Re-run with fewer --num-procs or a larger -l mem_free (which is "
                "PER SLOT under -pe smp).\n")
            _log_write(logfile, err_text)
            print(err_text, file=sys.stderr)
            _abort_profile(profile, args.profile_json, "stream.", "stream_worker_died")
            shutil.rmtree(stream_root, ignore_errors=True)
            sys.exit(1)
        except Exception as e:
            err_text = f"[FATAL ERROR] Streaming worker failed with exception:\n{e}\n"
            _log_write(logfile, err_text)
            print(err_text, file=sys.stderr)
            _abort_profile(profile, args.profile_json, "stream.", "stream_worker_exception")
            shutil.rmtree(stream_root, ignore_errors=True)
            sys.exit(1)

        # Without this, a systematic --cg-smarts / --cg-match-dict mismatch
        # produces an empty library and a log with nothing in it.
        if any(stream_skips.values()):
            _log_write(logfile, "[INFO] Environments skipped during streaming: "
                + ", ".join(f"{r}={n}" for r, n in sorted(stream_skips.items()) if n)
                + "\n")
        # The true totals _log_warn_capped promises, per key, including the ones
        # it suppressed after the cap. Reported separately from the skip counts:
        # a residue-level warning drops one slot, not the environment, and the
        # environment-level ones are already inside no_atomgroup above.
        if any(stream_warns.values()):
            _log_write(logfile, "[INFO] Warnings during streaming (incl. suppressed): "
                + ", ".join(f"{r}={n}" for r, n in sorted(stream_warns.items()) if n)
                + "\n")
        for reason, count in stream_skips.items():
            profile.set(f"stream_skipped_{reason}", count)
        for reason, count in stream_warns.items():
            profile.set(f"stream_warned_{reason}", count)
        # Did any worker write a single record? Shards exist on disk only when a
        # record was appended, so this is the cheapest true test of "did streaming
        # produce anything" without re-reading them.
        _streamed_any = any(
            fname.endswith(".pkl")
            for wdir in worker_dirs
            for _r, _d, _files in os.walk(wdir)
            for fname in _files)
        if not _streamed_any:
            # Nothing streamed at all. This is never a normal outcome for a
            # fragment that cleared the selection threshold, and every cause is
            # a configuration or I/O problem rather than a property of the
            # chemistry: a wrong -P/--pdb-dir, an unreadable mirror, a mismatched
            # --cg-match-dict-pkl, or every environment failing to build.
            if stream_skips.get("cg_elements_mismatch") and not any(
                    n for r, n in stream_skips.items() if r != "cg_elements_mismatch"):
                cause = ("every streamed environment was rejected by the CG element "
                         f"check ({stream_skips['cg_elements_mismatch']} records), so "
                         "--cg-smarts and --cg-match-dict-pkl most likely describe "
                         "different chemical groups")
            else:
                reasons = ", ".join(f"{r}={n}" for r, n in sorted(stream_skips.items()) if n)
                cause = ("no environment produced a single record"
                         + (f" (skip reasons: {reasons})" if reasons
                            else "; the shards were read but every environment was "
                                 "discarded without recording a reason, which points at "
                                 "-P/--pdb-dir not matching the mined parents"))
            # NB: this text must never contain the completion marker itself --
            # Frags.check_vdg_job_status does an unanchored substring test over the
            # whole log, so quoting the phrase here would make the fragment read as
            # finished, which is the exact failure this branch exists to prevent.
            err_text = (f"[WARNING] Streamed 0 vdG records for CG {CG}: {cause}. "
                        f"nr_vdgs/ is EMPTY and this fragment is NOT complete, so this "
                        f"log deliberately carries no completion marker.\n")
            _log_write(logfile, err_text)
            print(err_text, file=sys.stderr)
            profile.write(args.profile_json)
            shutil.rmtree(stream_root, ignore_errors=True)
            sys.exit(EXIT_NO_VDGS)
        profile.record_peak_rss("stream.")
        worker_shards_mb = profile.record_dir_size(
            "scratch_worker_shards_mb", workers_root)

        try:
          with profile.phase("merge_worker_shards"):
            for size_subset in subset_sizes:
                _merge_worker_dirs(
                    [os.path.join(wdir, str(size_subset)) for wdir in worker_dirs],
                    stream_dirs[size_subset])
            # Merged buckets and the per-worker shards they were built from
            # coexist until this rmtree, so scratch has to hold both at once.
            merged_mb = sum(profile.record_dir_size(
                f"scratch_merged_buckets_size_{size}_mb", stream_dirs[size])
                for size in subset_sizes)
            profile.max("scratch_peak_mb", worker_shards_mb + merged_mb)
            shutil.rmtree(workers_root, ignore_errors=True)
        except Exception as e:
            # Streaming has finished by now, so an uncaught raise here would abort
            # with the whole scratch tree (worker shards plus merged buckets, the
            # peak) still on disk -- on Wynton that is what fills /scratch for the
            # next job. Report it and clean up on the way out.
            err_text = (f"[FATAL ERROR] Failed to merge stream shards into AA "
                        f"buckets:\n{e}\n")
            _log_write(logfile, err_text)
            print(err_text, file=sys.stderr)
            _abort_profile(profile, args.profile_json, "merge.", "merge_worker_shards")
            shutil.rmtree(stream_root, ignore_errors=True)
            sys.exit(1)

    # AA buckets from all subset sizes share one clustering pool.
    jobs = _bucket_jobs(
        stream_dirs, subset_sizes, seq_sim_thresh, cg_automorphisms,
        vdglib_dir, logfile, max_num_to_clus)
    profile.set("buckets_queued", len(jobs))
    if jobs:
        # jobs is sorted largest-first, so job[0] is the bucket that sets both
        # the clustering tail and per-worker peak RSS. Sizing mem_free means
        # multiplying this by --num-procs, not by one worker.
        profile.set("largest_bucket_mb",
                    round(os.path.getsize(jobs[0][0]) / (1024 * 1024), 1))
        profile.set("largest_bucket_aa_key", "_".join(jobs[0][2]))
        profile.set("largest_bucket_subset_size", jobs[0][6])
    failed_buckets = []
    if jobs:
        ctx = mp.get_context("spawn")
        try:
          with profile.phase("cluster"):
            # ProcessPoolExecutor rather than mp.Pool: a worker killed by the OOM
            # killer (the common failure on the largest buckets) leaves mp.Pool
            # hanging forever on imap_unordered -- it respawns the worker but never
            # reports the lost task -- so the job burns its full h_rt and exits with
            # a partial nr_vdgs/ and no traceback. The executor raises
            # BrokenProcessPool instead.
            #
            # TRADE-OFF, and it is real: this drops the old maxtasksperchild=1,
            # which tore each worker down after its bucket so the Python-arena
            # memory behind adj/adj_d went back to the OS. The equivalent
            # (max_tasks_per_child) needs Python 3.11 and this env is 3.10.18, so
            # workers now persist and each retains the high-water mark of the
            # largest bucket it has handled. Jobs are dispatched largest-first, so
            # the peak is reached early and reused rather than growing; watch
            # cluster.rss_peak_child_mb. On 3.11+, add max_tasks_per_child=1 here
            # and this paragraph goes away.
            pool = concurrent.futures.ProcessPoolExecutor(
                max_workers=int(args.num_procs), mp_context=ctx)
            try:
                futures = [pool.submit(_run_one_bucket_strict, job) for job in jobs]
                for fut in concurrent.futures.as_completed(futures):
                    result = fut.result()
                    # total_vdgs is None only for a bucket the worker caught an
                    # exception on; it already wrote its .FAILED marker.
                    if result[1] is None:
                        failed_buckets.append(result[0])
                        continue
                    if profile.enabled and len(result) > 2:
                        stats = dict(result[2] or {})
                        # perm_group is a group *size*, not a tally. Summing it
                        # across buckets reports a number no bucket ever had; every
                        # other counter here is genuinely additive.
                        perm_group = stats.pop("perm_group", None)
                        if perm_group is not None:
                            profile.max("stage1.perm_group_max", perm_group)
                        # Stage-2 counters carry their own prefix; merging them
                        # under "stage1." would file them as Stage-1 cost, which
                        # is the confusion this split exists to remove.
                        stage2 = {k: stats.pop(k) for k in list(stats)
                                  if k.startswith("stage2_")}
                        n_orders = stage2.pop("stage2_slot_orders", None)
                        if n_orders is not None:
                            profile.max("stage2.slot_orders_max", n_orders)
                            if n_orders > 1:
                                profile.add("stage2.same_label_buckets")
                        profile.merge({k[len("stage2_"):]: v
                                       for k, v in stage2.items()}, prefix="stage2.")
                        profile.merge(stats, prefix="stage1.")
                        profile.add("stage1.buckets")
                        profile.add("stage2.buckets")
            finally:
                # cancel_futures matters: the executor's __exit__ would otherwise
                # run every queued bucket to completion before the exception
                # surfaced, so a systematic failure would still burn the whole
                # h_rt. mp.Pool's __exit__ called terminate() and stopped at once;
                # this is the nearest equivalent (in-flight buckets still finish).
                pool.shutdown(wait=True, cancel_futures=True)
        except BrokenProcessPool as e:
            err_text = (
                "[FATAL ERROR] A clustering worker died without returning a result "
                f"({e}). The usual cause is the OOM killer on the largest AA bucket: "
                "Stage 1 builds Python adjacency lists at roughly 60 B/edge, so a "
                "dense bucket can exceed the per-slot mem_free on its own. Re-run "
                "with fewer --num-procs or a larger -l mem_free (PER SLOT under "
                "-pe smp); the largest_bucket_* counters in the profile say which "
                "bucket to size for.\n")
            _log_write(logfile, err_text)
            print(err_text, file=sys.stderr)
            _abort_profile(profile, args.profile_json, "cluster.", "cluster_worker_died")
            shutil.rmtree(stream_root, ignore_errors=True)
            sys.exit(1)
        except Exception as e:
            err_text = f"[FATAL ERROR] Worker failed with exception:\n{e}\n"
            _log_write(logfile, err_text)
            print(err_text, file=sys.stderr)
            _abort_profile(profile, args.profile_json, "cluster.", "cluster_worker_exception")
            shutil.rmtree(stream_root, ignore_errors=True)
            sys.exit(1)  # hard-fail the whole run if the pool itself dies
        profile.record_peak_rss("cluster.")
    shutil.rmtree(stream_root, ignore_errors=True)

    hours, minutes, seconds = convert_time_elapsed(time.time() - start_time)
    _log_write(logfile,
        f"\nCompleted clus_and_deduplicate_vdgs.py for subset sizes "
        f"{', '.join(map(str, subset_sizes))} "
        f"in {hours} h, {minutes} mins, and {seconds} secs.\n")
    bad_buckets = set()
    for size_subset in subset_sizes:
        num_nr_vdgs, num_inputs, unreadable = _subset_output_counts(
            vdglib_dir, size_subset, logfile)
        _log_write(logfile,
            f"\t{num_nr_vdgs} nonredun. vdgs of subset size {size_subset} "
            f"from {num_inputs} input vdGs.\n")
        profile.set(f"nr_vdgs_size_{size_subset}", num_nr_vdgs)
        profile.set(f"inputs_size_{size_subset}", num_inputs)
        bad_buckets.update(f"{size_subset}/{label} (unreadable npz)"
                           for label in unreadable)
        bad_buckets.update(f"{size_subset}/{label} (.FAILED)"
                           for label in _scan_failed_markers(vdglib_dir, size_subset))
    # A worker failure normally shows up above as its own marker; it is added here
    # too in case writing the marker is what failed.
    _marked = set(bad_buckets)
    bad_buckets.update(
        f"{label} (worker error)" for label in failed_buckets
        if not any(b.startswith(f"{size}/{label} ")
                   for size in subset_sizes for b in _marked))
    profile.set("buckets_failed", len(bad_buckets))
    profile.write(args.profile_json)

    if bad_buckets:
        # Exit non-zero so the wrapper's check=True fires and the fragment log
        # never gets its 'Job completed.' line -- an incomplete library must not
        # read as finished to check_vdg_job_status.
        msg = (f"[ERROR] {len(bad_buckets)} AA bucket(s) are missing or failed to "
               f"cluster; nr_vdgs/ is INCOMPLETE: "
               f"{', '.join(sorted(bad_buckets))}\n")
        _log_write(logfile, msg)
        print(msg, file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
