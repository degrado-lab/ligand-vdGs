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
from collections import OrderedDict, deque, namedtuple

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
    mismatched match dict, a parent db that stopped resolving -- that is one
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
from ligand_vdgs.functions import parent_db
from ligand_vdgs.functions.clus_helpers import (
    get_vdg_subsets_target_size, select_diverse_pdbIDs, _aa_tmp_dir,
    _stream_root, stage1_record_is_complete, records_to_columns, concat_columns,
    save_shard, load_shard, save_columns, load_columns, load_stage1, COLUMNS_FILE)
from ligand_vdgs.functions.vdg_struct_utils import (VDM_OCC, cg_slot_occupancy,
    get_cg_atoms, is_hydrogen)
from ligand_vdgs.functions.vdg_npz_utils import (write_cg_symmetry, name_selstr,
    parse_pdb_with_retry)
from ligand_vdgs.functions.vdg_fp_utils import build_perm_group, slot_orders
from ligand_vdgs.functions.align_and_cluster import butina_partition, neighbor_csr, stage1_edges
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


# Streaming chunks per worker process. >1 so the executor can rebalance across
# workers; small enough that the per-chunk dispatch cost stays a rounding error
# against a chunk's parse cost.
STREAM_CHUNKS_PER_WORKER = 4

# Buckets with at least this many records are clustered as phased tasks --
# Stage-1 row blocks, graph assembly, Stage-2 cluster blocks, write -- so the
# pool stays busy through the tail instead of one worker holding the largest
# bucket alone. Smaller ones run start to finish in one task.
_SPLIT_MIN_RECORDS = 8_000
# Blocks per split bucket, as a multiple of --num-procs: enough granularity for
# the pool to absorb the data-dependent cost of the exact fits and of Stage 2.
_BLOCKS_PER_PROC = 4


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
        help="Worker processes for streaming and clustering.",)
    parser.add_argument("--keep-stream-dir", action="store_true",
        help="Leave the streamed bucket columns on scratch after the run "
        "instead of deleting them, and print their location.")
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

# Well under the 255-byte filename limit on ext4/xfs, with room for the hash,
# the flush index and the extension.
_MAX_BUCKET_FNAME_BYTES = 200


def _bucket_fname(aa_key, flush_idx=None):
    """Name of the on-disk file holding one AA-composition bucket.

    The aa_key is written **whole**, never truncated: it is the only record of the
    shard's composition, and `_merge_worker_dirs` recovers it by parsing this name
    back with `_aa_key_from_bucket_fname`. A truncated key would merge distinct
    buckets into one file and mislabel their residues, silently. The sha1 is kept
    for name stability and legibility; it is not what makes the name unique.
    """
    h = hashlib.sha1(aa_key.encode()).hexdigest()[:16]
    stem = f"{aa_key}__{h}" if flush_idx is None else f"{aa_key}__{h}__{flush_idx:04d}"
    fname = f"{stem}.npz"
    if len(fname.encode()) > _MAX_BUCKET_FNAME_BYTES:
        raise ValueError(
            f"AA bucket key {aa_key!r} ({len(aa_key)} chars) does not fit in a "
            f"{_MAX_BUCKET_FNAME_BYTES}-byte filename. Labels are <=4 chars, so this "
            "needs a subset size around 40 -- far past MAX_SUBSET_SIZE. Store the key "
            "inside the shard instead of shortening it here; the reader parses it "
            "back out of the name.")
    return fname


_BUCKET_FNAME_RE = re.compile(r"^(?P<key>[^_].*?)__[0-9a-f]{16}(?:__\d{4})?\.npz$")


def _aa_key_from_bucket_fname(fname):
    """Recover the exact aa_key from a name written by `_bucket_fname`, or None.

    None rather than a raise: the callers walk directories, and a file that is not
    one of ours (an editor backup, a stray from a killed run) must not abort a merge
    over every bucket in the library. 
    """
    match = _BUCKET_FNAME_RE.match(fname)
    return match.group("key") if match else None

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

# ---- bucket output ----------------------------------------------------------

# One Stage-2 subgroup: its Stage-1 cluster rank, its number within that
# cluster, the row stored as the nr vdG, every member row (nr included), and
# the exact pose radius about the nr row.
Subgroup = namedtuple("Subgroup", "stage1_id stage2_id nr_idx member_idxs radius")

_QUALITY_FIELDS = ("cg_max_b", "cg_min_occ", "vdm_max_b", "vdm_min_occ")


def _write_bucket_npz(vdglib_dir, size_subset, reordered_AAs, cols, subgroups,
                      parent_pdb_dir):
    """Write one AA bucket: nr_ rows gathered from `cols` by each subgroup's nr
    index, mem_ rows for every other member, largest subgroup first."""
    if not subgroups:
        return
    subgroups = sorted(subgroups,
                       key=lambda g: (-len(g.member_idxs), g.stage1_id, g.stage2_id))
    C = len(subgroups)
    n_cg = cols["cg_names"].shape[1]
    num_vdms = cols["scrr_seg"].shape[1]
    nr = np.asarray([g.nr_idx for g in subgroups], dtype=np.intp)
    sizes = np.asarray([len(g.member_idxs) for g in subgroups], dtype=np.int32)
    mem = np.concatenate([g.member_idxs[g.member_idxs != g.nr_idx] for g in subgroups]
                         + [np.empty(0, dtype=np.int32)]).astype(np.intp)
    if mem.size != int(sizes.sum()) - C:
        raise ValueError(f"mem_ row count {mem.size} != summed cluster_size minus "
                         f"one per cluster ({int(sizes.sum()) - C}); cluster "
                         "membership was lost.")
    # Distinct PDB *entries*, not biounit stems: `1abc_1` and `1abc_2` are two
    # assemblies of one deposition, so counting stems would readmit exactly the
    # NCS/homolog inflation this field exists to strip out of cluster_size.
    entries = np.asarray([parent_db.entry_of(b) for b in cols["biounit"]])
    num_parents = np.asarray([np.unique(entries[g.member_idxs]).size
                              for g in subgroups], dtype=np.int32)
    cgvdmbb = np.asarray(cols["cgvdmbb"][nr], dtype=np.float32)

    # Resolved position by position across the bucket's nr rows: a blank PDB
    # element column is common, and streaming already admits records whose
    # blanks it could not check. Only positions where two rows are non-blank can
    # disagree, and a disagreement means the automorphism group does not
    # preserve elements, so the library's atom correspondence is suspect.
    elements = np.char.capitalize(np.char.strip(cols["cg_elements"][nr]))
    resolved = np.array([""] * n_cg, dtype="U2")
    for k in range(n_cg):
        seen = {e for e in elements[:, k] if e}
        if len(seen) > 1:
            raise ValueError(
                f"CG element sequence differs between nr vdGs in bucket "
                f"{'_'.join(reordered_AAs)} at atom {k}: {sorted(seen)}. "
                "CG automorphisms must preserve element.")
        if seen:
            resolved[k] = seen.pop()

    arrays = {
        # Bucket-level (not per nr vdG): the label of each residue slot. U4 fits
        # every label -- the 20 resnames, 'bb', and 'X'.
        "aa_bucket_parts": np.asarray(list(reordered_AAs), dtype="U4"),
        "cluster_id": np.arange(1, C + 1, dtype=np.int32),
        "cluster_size": sizes,
        # Exact pose radius: the greatest symmetry-aware RMSD from any member to
        # the stored nr vdG. Sphere exclusion bounds this by the Stage-1 cutoff
        # for a whole Stage-1 cluster, but Stage 2 splits those and re-picks a
        # row, so the stored value is measured, not assumed. A query within
        # `tau` of some discarded member is within `tau + radius` of the nr vdG.
        "cluster_pose_radius": np.asarray([g.radius for g in subgroups], dtype=np.float32),
        # Distinct parent structures behind the cluster, counting the nr vdG.
        # `cluster_size` counts observations, which NCS copies and homologous
        # entries inflate; this is the support figure to do statistics on.
        "cluster_num_parents": num_parents,
        "first_stage_cluster_id": np.asarray([g.stage1_id for g in subgroups], dtype=np.int32),
        "second_stage_cluster_id": np.asarray([g.stage2_id for g in subgroups], dtype=np.int32),
        "nr_cg_coords": cgvdmbb[:, :n_cg],
        "nr_vdm_bb_coords": cgvdmbb[:, n_cg:].reshape(C, num_vdms, 3, 3),
        # Backbone carbonyl O per vdM slot, same frame as nr_vdm_bb_coords (the
        # environment frame, not the parent's), NaN where the residue has no O. Stored for contact statistics
        # (tools/h_class_diagnostic.py --contact-atoms N,CA,C,O); it is not a
        # Stage-1 atom and must never be appended to the RMSD input.
        "nr_vdm_o_coords": np.asarray(cols["vdm_o"][nr], dtype=np.float32),
        # Parents are stored as the biounit stem plus the database directory once
        # per file; readers rebuild the path with vdg_npz_utils.resolve_parent_pdb_path,
        # which takes an override so a relocated database still resolves.
        "parent_pdb_dir": np.asarray(str(parent_pdb_dir or ""), dtype="U512"),
        "nr_parent_biounit": cols["biounit"][nr],
        "nr_scrr_seg": cols["scrr_seg"][nr],
        "nr_scrr_chain": cols["scrr_chain"][nr],
        "nr_scrr_resnum": cols["scrr_resnum"][nr],
        "nr_scrr_resname": cols["scrr_resname"][nr],
        # Atom names are per CG atom: two records can be the same chemical group
        # in differently-named ligands. The CG's residue is one value per record
        # and its elements one row per bucket, since every record keeps the
        # fragment's SMARTS-slot order.
        "nr_cg_names": cols["cg_names"][nr],
        "cg_elements": resolved,
        "nr_cg_seg": cols["cg_seg"][nr],
        "nr_cg_chain": cols["cg_chain"][nr],
        "nr_cg_resnum": cols["cg_resnum"][nr],
        "nr_cg_resname": cols["cg_resname"][nr],
        # Per-vdM-slot provenance (SLOT_* in vdg_struct_utils) of the nr row
        # only: clustering never consults it and flank-sequence similarity
        # tolerates mismatches, so a cluster can mix e.g. SLOT_NO_SC and
        # SLOT_BB_CLOSER members; weighting this code by cluster_size is invalid.
        "nr_slot_flag": cols["slot_flags"][nr],
        # A "mem_" row is a clustered observation that is *not* the nr vdG, whose
        # identity lives in the nr_ arrays; each vdG is stored exactly once.
        # `cluster_size` counts the nr row, so a reader treating these as "the
        # members" would be off by one per cluster.
        "mem_cluster_id": np.repeat(np.arange(1, C + 1, dtype=np.int32), sizes - 1),
        "mem_parent_biounit": cols["biounit"][mem],
        "mem_scrr_seg": cols["scrr_seg"][mem],
        "mem_scrr_chain": cols["scrr_chain"][mem],
        "mem_scrr_resnum": cols["scrr_resnum"][mem],
        "mem_scrr_resname": cols["scrr_resname"][mem],
        "mem_cg_names": cols["cg_names"][mem],
        "mem_cg_seg": cols["cg_seg"][mem],
        "mem_cg_chain": cols["cg_chain"][mem],
        "mem_cg_resnum": cols["cg_resnum"][mem],
        "mem_cg_resname": cols["cg_resname"][mem],
    }
    # Measured on the atoms that enter the vdG -- CG atoms and the contacting
    # residues' heavy atoms. Stored so a stricter cut is a read-path decision.
    for k, key in enumerate(_QUALITY_FIELDS):
        arrays[f"nr_{key}"] = cols["quality"][nr, k]
        arrays[f"mem_{key}"] = cols["quality"][mem, k]

    if len(str(parent_pdb_dir or "")) >= arrays["parent_pdb_dir"].dtype.itemsize // 4:
        raise ValueError(
            f"parent_pdb_dir {parent_pdb_dir!r} reaches the dtype width and may be "
            f"truncated; widen the dtype.")

    # Written via tmp + os.replace, never straight onto the final name: an SGE
    # kill or a quota hit part-way through savez_compressed would otherwise leave
    # a truncated npz that both readers treat as a warning and skip. The temp
    # name deliberately does not end in '.npz': a SIGKILL is uncatchable, and a
    # leftover *.npz in the bucket dir is loaded as a real bucket by every reader
    # that walks for '*.npz'. Written through a file handle because
    # np.savez_compressed appends '.npz' to a *path* that lacks it.
    path = _bucket_npz_path(vdglib_dir, size_subset, reordered_AAs)
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



# ---- clustering tasks -------------------------------------------------------

# Static per-bucket facts every task needs: `key` is (size_subset, aa_label),
# `aa_parts` the slot labels, `n` the record count, `bucket_dir` the merged
# columns on scratch.
Bucket = namedtuple("Bucket", "key size_subset aa_parts n bucket_dir")
# Static per-run facts.
Run = namedtuple("Run", "cg_automorphisms seq_sim_thresh vdglib_dir logfile parent_pdb_dir")


def _bucket_geometry(run, bucket, cols):
    n_cg = cols["cg_names"].shape[1]
    n_total = cols["cgvdmbb"].shape[1]
    cutoff = normalize_rmsd(n_total, "cgvdmbb")
    # The vdG's full symmetry group: CG automorphisms crossed with orderings
    # of interchangeable same-label vdM slots. Both are minimised over inside
    # the Stage-1 distance, so one record is one physical environment.
    perm_group = build_perm_group(run.cg_automorphisms, n_cg, bucket.aa_parts)
    return n_cg, cutoff, perm_group


def _rank_clusters(clusters):
    """Largest first, ties by formation order, so cluster 1 is the most populated."""
    return sorted(clusters, key=lambda members: -len(members))


def _row_blocks(n, n_blocks):
    """Row ranges over `range(n - 1)` with equal shares of the candidate-pair
    triangle (row i has n - 1 - i candidates)."""
    cum = np.cumsum(np.arange(n - 1, 0, -1), dtype=np.int64)
    targets = cum[-1] * np.arange(1, n_blocks) / n_blocks
    cuts = np.searchsorted(cum, targets) + 1
    bounds = np.unique(np.concatenate([[0], cuts, [n - 1]]))
    return [(int(a), int(b)) for a, b in zip(bounds[:-1], bounds[1:])]


def _stage2_blocks(sizes, n_blocks):
    """Stage-1 cluster ids (1-based ranks) grouped into blocks of similar Stage-2
    cost, which is ~quadratic in cluster size; singletons cost nothing and share
    one block."""
    sizes = np.asarray(sizes)
    ids = np.flatnonzero(sizes > 1) + 1
    order = ids[np.argsort(-sizes[ids - 1], kind="stable")]
    blocks = [[] for _ in range(max(1, min(n_blocks, order.size)))]
    load = [0] * len(blocks)
    for cid in order.tolist():
        k = int(np.argmin(load))
        blocks[k].append(cid)
        load[k] += int(sizes[cid - 1]) ** 2
    singletons = (np.flatnonzero(sizes == 1) + 1).tolist()
    if singletons:
        blocks.append(singletons)
    return [b for b in blocks if b]


def _stage2_for_clusters(run, bucket, cols, clusters, cluster_ids, perm_group, n_cg):
    """Partition each listed Stage-1 cluster by flanking context and pick each
    subgroup's pose prototype. Returns (subgroups, stats)."""
    data = np.asarray(cols["cgvdmbb"], dtype=np.float32)
    flank_seq, flank_ca = cols["flank_seq"], cols["flank_ca"]
    orders = slot_orders(bucket.aa_parts)
    n_flank = flank_ca.shape[1]
    # get_leader_clusters mixes RMSD and sequence dissimilarity as
    # seq_dissim * seq_weight + RMSD with seq_weight=0.5, so the threshold lives
    # on that combined scale: the flank RMSD cutoff plus half the allowed
    # sequence-dissimilarity budget. Missing flanks must not tighten the
    # geometric criterion: the cutoff is based on the expected flattened size,
    # while pair RMSDs use only coordinate rows finite in both structures.
    thresh = normalize_rmsd(n_flank, "flankbb") + (1.0 - run.seq_sim_thresh) / 2.0
    subgroups = []
    stats = {"stage2_wall_s": 0.0, "stage2_same_label_wall_s": 0.0,
             "stage2_slot_orders": len(orders),
             "stage2_largest_cluster": 0, "stage2_largest_cluster_s": 0.0}
    t_all = time.time()
    for cid in cluster_ids:
        members = np.asarray(clusters[cid - 1], dtype=np.int32)
        if members.size == 1:
            subgroups.append(Subgroup(cid, 1, int(members[0]), members, 0.0))
            continue
        t0 = time.time()
        # Same-label slot orderings, the ones Stage 1 already quotients out via
        # build_perm_group. Without them Stage 2 compares slot 1's flank to slot
        # 1's flank positionally and splits a swapped pair into two subgroups,
        # halving the cluster_size/cluster_num_parents of one real mode.
        assigns = clust.get_leader_clusters(
            zip([flank_seq[members].tolist(), list(flank_ca[members])],
                ["flankseq", "flankbb"]),
            thresh, missing_seq_similarity=run.seq_sim_thresh,
            final_exact_medoid_pass=True, final_reassign_once=True,
            slot_orders=orders)
        clust.clear_caches()
        n_before = len(subgroups)
        for sub_num, local in assigns.items():
            if not local:
                continue
            global_idxs = members[np.asarray(local, dtype=np.intp)]
            # Stage 2 partitions on flanking context, so its own medoid is chosen
            # on the wrong metric: the stored row has to be central in *pose*,
            # or the radius recorded next to it means nothing.
            centre, radius = clust.pose_minimax_prototype(
                data, global_idxs, n_cg, perm_group)
            subgroups.append(Subgroup(cid, int(sub_num), int(centre), global_idxs, radius))
        if len(subgroups) == n_before:
            _log_write(run.logfile,
                f"[WARNING] Stage 1 cluster {cid} in {bucket.key[1]} lost "
                f"{members.size} envs in Stage 2.\n")
        dt = time.time() - t0
        if len(orders) > 1:
            stats["stage2_same_label_wall_s"] += dt
        if members.size > stats["stage2_largest_cluster"]:
            stats["stage2_largest_cluster"] = int(members.size)
            stats["stage2_largest_cluster_s"] = dt
    stats["stage2_wall_s"] = time.time() - t_all
    return subgroups, stats


def _task_small(run, bucket, _payload):
    """A whole bucket start to finish in one process."""
    cols = load_columns(bucket.bucket_dir)
    n_cg, cutoff, perm_group = _bucket_geometry(run, bucket, cols)
    stats = {"records": bucket.n, "perm_group": len(perm_group)}
    t0 = time.time()
    clusters = _rank_clusters(clust.get_butina_clusters(
        cols["cgvdmbb"], cutoff, n_cg, perm_group, counters=stats))
    stats["stage1_wall_s"] = time.time() - t0
    stats["stage1_clusters"] = len(clusters)
    subgroups, s2 = _stage2_for_clusters(
        run, bucket, cols, clusters, range(1, len(clusters) + 1), perm_group, n_cg)
    stats.update(s2)
    t0 = time.time()
    _write_bucket_npz(run.vdglib_dir, bucket.size_subset, bucket.aa_parts, cols,
                      subgroups, run.parent_pdb_dir)
    stats["write_s"] = time.time() - t0
    return sum(len(g.member_idxs) for g in subgroups), stats


def _edges_path(bucket, k):
    return os.path.join(bucket.bucket_dir, f"edges_{k:04d}.npy")


def _task_stage1_block(run, bucket, payload):
    """Within-cutoff pairs for one block of rows, written to scratch."""
    k, start, stop = payload
    cols = {"cgvdmbb": np.asarray(load_stage1(bucket.bucket_dir, mmap=True)),
            "cg_names": load_shard(os.path.join(bucket.bucket_dir, COLUMNS_FILE))["cg_names"]}
    n_cg, cutoff, perm_group = _bucket_geometry(run, bucket, cols)
    stats = {}
    t0 = time.time()
    qi, qj = stage1_edges(cols["cgvdmbb"], cutoff, n_cg, perm_group,
                          row_range=(start, stop), counters=stats)
    stats["stage1_wall_s"] = time.time() - t0
    path = _edges_path(bucket, k)
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "wb") as handle:
        np.save(handle, np.stack([qi, qj]))
    os.replace(tmp, path)
    return stats


def _partition_path(bucket):
    return os.path.join(bucket.bucket_dir, "partition.npz")


def _task_stage1_partition(run, bucket, n_blocks):
    """Assemble the row blocks' edges into the neighbour graph and partition it.
    The partition goes to scratch (largest cluster first); the sizes come back."""
    t0 = time.time()
    parts = [np.load(_edges_path(bucket, k)) for k in range(n_blocks)]
    qi = np.concatenate([part[0] for part in parts])
    qj = np.concatenate([part[1] for part in parts])
    del parts
    clusters = _rank_clusters(butina_partition(*neighbor_csr(bucket.n, qi, qj)))
    del qi, qj
    sizes = np.asarray([len(c) for c in clusters], dtype=np.int32)
    indptr = np.zeros(sizes.size + 1, dtype=np.int64)
    np.cumsum(sizes, out=indptr[1:])
    order = np.concatenate(clusters) if clusters else np.empty(0, dtype=np.int32)
    path = _partition_path(bucket)
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "wb") as handle:
        np.savez(handle, order=order, indptr=indptr)
    os.replace(tmp, path)
    for k in range(n_blocks):
        os.remove(_edges_path(bucket, k))
    return sizes, {"stage1_partition_s": time.time() - t0,
                   "stage1_clusters": int(sizes.size)}


def _load_partition(bucket):
    with np.load(_partition_path(bucket)) as data:
        order, indptr = data["order"], data["indptr"]
    return [order[indptr[c]:indptr[c + 1]] for c in range(indptr.size - 1)]


def _task_stage2_block(run, bucket, cluster_ids):
    cols = load_columns(bucket.bucket_dir)
    n_cg, _cutoff, perm_group = _bucket_geometry(run, bucket, cols)
    clusters = _load_partition(bucket)
    return _stage2_for_clusters(run, bucket, cols, clusters, cluster_ids, perm_group, n_cg)


def _task_write(run, bucket, subgroups):
    t0 = time.time()
    cols = load_columns(bucket.bucket_dir)
    _write_bucket_npz(run.vdglib_dir, bucket.size_subset, bucket.aa_parts, cols,
                      subgroups, run.parent_pdb_dir)
    os.remove(_partition_path(bucket))
    return sum(len(g.member_idxs) for g in subgroups), {"write_s": time.time() - t0}


_TASKS = {"small": _task_small, "stage1_block": _task_stage1_block,
          "stage1_partition": _task_stage1_partition,
          "stage2_block": _task_stage2_block, "write": _task_write}


def _run_task(kind, run, bucket, payload):
    """Pool entry point. One bad bucket must not abort the pool: every
    nr_vdgs/*.npz already written would stay on disk and the fragment dir would
    look like a finished library. Failures come back as text; the driver drops
    a .FAILED marker, skips the bucket's remaining tasks, and exits non-zero
    at the end."""
    try:
        return kind, bucket.key, _TASKS[kind](run, bucket, payload), None
    except Exception as exc:
        return kind, bucket.key, None, (
            f"[WORKER ERROR] AA bucket: {bucket.key[1]} (subset size "
            f"{bucket.size_subset}, task {kind})\nException: {exc}\n"
            f"Traceback:\n{traceback.format_exc()}")

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
    pdb_file = parent_db.structure_path(pdb_dir, biounit)

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

# Total parsed atoms the streaming workers may hold in their PDB caches, divided
# across workers in main(). A default ProDy parse costs ~160 B/atom (coordinates
# plus the per-atom label, index and flag arrays), so this is ~1.5 GB over a
# 20-worker run; environments are chunked by file, so a worker rarely needs
# more than the structure it is on.
_PDB_CACHE_ATOM_BUDGET = 10_000_000


def _struct_atom_count(struct):
    """Cache weight of a parsed structure; sentinels (parse failures) weigh 0."""
    try:
        return int(struct.numAtoms())
    except Exception:
        return 0


def _flush_buckets_to_disk(bucket_records, worker_tmp_dir, flush_idx=None):
    """Write each bucket's buffered records as one column shard."""
    for aa_key, records in bucket_records.items():
        fname = _bucket_fname(aa_key, flush_idx=flush_idx)
        save_shard(os.path.join(worker_tmp_dir, fname), records_to_columns(records))


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
             "bad_cg_coords": 0, "cg_elements_mismatch": 0, "unparsable_line": 0,
             "incomplete_stage1": 0}
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
                             re_ordered_scrr, re_ordered_slot_flags,
                             re_ordered_bbo
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
                            "biounit": str(_biounit),
                            "scrr": [[str(s), str(ch), int(r), str(rn)] for (s, ch, r, rn) in re_ordered_scrr],
                            "cg_names": list(cg_names), "cg_elements": list(cg_elements),
                            "cg_seg": str(cg_seg), "cg_chain": str(cg_chain),
                            "cg_resnum": int(cg_resnum), "cg_resname": str(cg_resname),
                            "slot_flags": [int(f) for f in re_ordered_slot_flags],
                            "quality": tuple(float(q) for q in env_quality),
                            # Carbonyl O per slot, NaN where the parent lacks it.
                            # Kept out of cgvdmbb so it never enters an RMSD.
                            "bbo": [np.asarray(x, dtype=np.float32) for x in re_ordered_bbo],
                        }
                        # Every Stage-1 atom is mandatory; a record missing one
                        # would only fail much later, inside the clustering.
                        if not stage1_record_is_complete(rec, expected_n_cg, size_subset):
                            skips["incomplete_stage1"] += 1
                            continue
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

    for size, size_records in bucket_records.items():
        _flush_buckets_to_disk(size_records, worker_tmp_dirs[size])

    warns = {key: count - warn_counts_at_entry.get(key, 0)
             for key, count in _WARN_COUNTS.items()
             if count - warn_counts_at_entry.get(key, 0) > 0}
    return worker_tmp_root, skips, warns


def _cap_distinct_entries(cols, aa_key, max_num_to_clus, logfile):
    """DEBUG ONLY (-m): keep the records of at most `max_num_to_clus` diverse
    PDB entries; one entry can still contribute several records."""
    entries = [parent_db.entry_of(b) for b in cols["biounit"]]
    if len(set(entries)) <= max_num_to_clus:
        return cols
    keep = set(select_diverse_pdbIDs(entries, max_num_to_clus))
    rows = np.asarray([e in keep for e in entries])
    _log_write(logfile, f"\t{len(entries)} vdGs -> {int(rows.sum())} vdGs from "
                        f"{len(keep)} diverse PDB IDs for {aa_key}.\n")
    return {key: arr[rows] for key, arr in cols.items()}


def _merge_worker_dirs(worker_dirs, final_dir, max_num_to_clus, logfile):
    """Merge per-worker shards into one column set per aa_key under `final_dir`.

    One aa_key at a time, so peak memory is bounded by the largest single bucket
    rather than the sum of every bucket in the library. Returns
    ``{aa_key: record count}``.
    """
    os.makedirs(final_dir, exist_ok=True)
    files_by_key = {}  # aa_key -> list of shard paths
    for d in worker_dirs:
        if not os.path.isdir(d):
            continue
        for fname in os.listdir(d):
            if not fname.endswith(".npz"):
                continue
            aa_key = _aa_key_from_bucket_fname(fname)
            if aa_key is None:
                print(f"[WARNING] Ignoring unrecognized file in worker dir: "
                      f"{os.path.join(d, fname)}", file=sys.stderr)
                continue
            files_by_key.setdefault(aa_key, []).append(os.path.join(d, fname))

    counts = {}
    for aa_key, paths in sorted(files_by_key.items()):
        parts = []
        for path in sorted(paths):
            try:
                parts.append(load_shard(path))
            except Exception as exc:
                # Shards are written whole or not at all, so reaching here means
                # real corruption (bad block, truncated NFS write) rather than an
                # interrupted worker. Raise rather than skip.
                raise RuntimeError(
                    f"Unreadable stream shard {path} for AA bucket {aa_key}: "
                    f"{type(exc).__name__}: {exc}") from exc
        cols = concat_columns(parts)
        del parts
        if max_num_to_clus is not None:
            cols = _cap_distinct_entries(cols, aa_key, max_num_to_clus, logfile)
        save_columns(os.path.join(final_dir, aa_key), cols)
        counts[aa_key] = len(cols["biounit"])
    return counts


def _buckets(counts_by_size, stream_dirs):
    """Every bucket of every subset size, largest first.

    Bucket cost is superlinear in record count and the distribution is very
    lopsided (`bb_bb` and `ARG` run ~500x the median). Queued in name order, one
    of those can be handed out last and set the wall time on its own; largest
    first puts them in while the pool is still empty.
    """
    buckets = []
    for size_subset, counts in counts_by_size.items():
        for aa_key, n in counts.items():
            buckets.append(Bucket(
                key=(size_subset, aa_key), size_subset=size_subset,
                aa_parts=tuple(aa_key.split("_")), n=n,
                bucket_dir=os.path.join(stream_dirs[size_subset], aa_key)))
    buckets.sort(key=lambda b: (-b.n, b.key))
    return buckets


def _record_bucket_stats(profile, stats):
    """Fold one bucket's counters into the profile. perm_group and the Stage-2
    slot-order / largest-cluster figures are per-bucket sizes, not tallies, so
    they are kept as maxima; everything else is additive."""
    if not profile.enabled or not stats:
        return
    stats = dict(stats)
    perm_group = stats.pop("perm_group", None)
    if perm_group is not None:
        profile.max("stage1.perm_group_max", perm_group)
    n_orders = stats.pop("stage2_slot_orders", None)
    if n_orders is not None:
        profile.max("stage2.slot_orders_max", n_orders)
        if n_orders > 1:
            profile.add("stage2.same_label_buckets")
    largest = stats.pop("stage2_largest_cluster", None)
    largest_s = stats.pop("stage2_largest_cluster_s", None)
    if largest is not None and largest >= profile.counters.get("stage2.largest_cluster", 0):
        profile.set("stage2.largest_cluster", largest)
        profile.set("stage2.largest_cluster_s", largest_s)
    stage2 = {k[len("stage2_"):]: v for k, v in stats.items() if k.startswith("stage2_")}
    stage1 = {k: v for k, v in stats.items() if not k.startswith("stage2_")}
    profile.merge(stage2, prefix="stage2.")
    profile.merge(stage1, prefix="stage1.")


def _cluster_buckets(pool, num_procs, run, buckets, profile):
    """Drive every bucket through the pool; returns the labels that failed.

    Small buckets are one task each. A bucket of at least _SPLIT_MIN_RECORDS
    records is four phases: Stage-1 row blocks, then graph assembly and
    partition, then Stage-2 cluster blocks, then the write. Large buckets'
    blocks are queued first, and a bucket's next phase goes to the *front* of
    the queue when its current one completes; small buckets are released only
    to keep the pool fed, so the executor's own FIFO never holds more than one
    pool's worth of them ahead of a large bucket's next phase.
    """
    n_blocks = _BLOCKS_PER_PROC * num_procs
    pending = deque()
    state = {}
    for bucket in buckets:
        if bucket.n >= _SPLIT_MIN_RECORDS:
            blocks = _row_blocks(bucket.n, n_blocks)
            state[bucket.key] = {"bucket": bucket, "blocks": len(blocks),
                                 "blocks_done": 0, "stats": {"records": bucket.n},
                                 "subgroups": [], "stage2_blocks": 0,
                                 "stage2_done": 0}
            for k, (start, stop) in enumerate(blocks):
                pending.append(("stage1_block", bucket, (k, start, stop)))
        else:
            pending.append(("small", bucket, None))
    profile.set("buckets_split", len(state))

    failed, done, inflight = set(), set(), {}

    def _fail(bucket, err_text):
        failed.add(bucket.key)
        _log_write(run.logfile, err_text + "\n")
        _write_failed_marker(run.vdglib_dir, bucket.size_subset, bucket.key[1], err_text)

    def _fill():
        while pending and len(inflight) < 2 * num_procs:
            kind, bucket, payload = pending.popleft()
            if bucket.key in failed:
                continue
            fut = pool.submit(_run_task, kind, run, bucket, payload)
            inflight[fut] = (kind, bucket)

    def _push_front(tasks):
        for task in reversed(tasks):
            pending.appendleft(task)

    _fill()
    while inflight:
        finished, _ = concurrent.futures.wait(
            inflight, return_when=concurrent.futures.FIRST_COMPLETED)
        for fut in finished:
            kind, bucket = inflight.pop(fut)
            _kind, key, value, err_text = fut.result()
            if key in failed:
                continue
            if err_text is not None:
                _fail(bucket, err_text)
                continue
            if kind == "small":
                total, stats = value
                _record_bucket_stats(profile, stats)
                done.add(key)
                continue
            st = state[key]
            if kind == "stage1_block":
                st["blocks_done"] += 1
                for name, v in value.items():
                    st["stats"][name] = st["stats"].get(name, 0) + v
                if st["blocks_done"] == st["blocks"]:
                    _push_front([("stage1_partition", bucket, st["blocks"])])
            elif kind == "stage1_partition":
                sizes, stats = value
                st["stats"].update(stats)
                st["stats"]["perm_group"] = len(build_perm_group(
                    run.cg_automorphisms, len(run.cg_automorphisms[0]), bucket.aa_parts))
                blocks = _stage2_blocks(sizes, n_blocks)
                st["stage2_blocks"] = len(blocks)
                _push_front([("stage2_block", bucket, ids) for ids in blocks])
            elif kind == "stage2_block":
                subgroups, stats = value
                st["subgroups"].extend(subgroups)
                st["stage2_done"] += 1
                for name, v in stats.items():
                    if name.endswith("largest_cluster") or name.endswith("largest_cluster_s"):
                        continue
                    st["stats"][name] = st["stats"].get(name, 0) + v
                if stats["stage2_largest_cluster"] >= st["stats"].get("stage2_largest_cluster", 0):
                    st["stats"]["stage2_largest_cluster"] = stats["stage2_largest_cluster"]
                    st["stats"]["stage2_largest_cluster_s"] = stats["stage2_largest_cluster_s"]
                st["stats"]["stage2_slot_orders"] = stats["stage2_slot_orders"]
                if st["stage2_done"] == st["stage2_blocks"]:
                    _push_front([("write", bucket, st["subgroups"])])
                    st["subgroups"] = []
            elif kind == "write":
                total, stats = value
                st["stats"].update(stats)
                _record_bucket_stats(profile, st["stats"])
                done.add(key)
        _fill()
    if profile.enabled:
        profile.set("stage1.buckets", len(done))
        profile.set("stage2.buckets", len(done))
    return sorted(key[1] for key in failed)


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
    a marker is the only trace left by a bucket whose process is gone.
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
            fname.endswith(".npz")
            for wdir in worker_dirs
            for _r, _d, _files in os.walk(wdir)
            for fname in _files)
        if not _streamed_any:
            # Nothing streamed at all. This is never a normal outcome for a
            # fragment that cleared the selection threshold, and every cause is
            # a configuration or I/O problem rather than a property of the
            # chemistry: a wrong -P/--pdb-dir, an unreadable parent db, a mismatched
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
            counts_by_size = {}
            for size_subset in subset_sizes:
                counts_by_size[size_subset] = _merge_worker_dirs(
                    [os.path.join(wdir, str(size_subset)) for wdir in worker_dirs],
                    stream_dirs[size_subset], max_num_to_clus, logfile)
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
    buckets = _buckets(counts_by_size, stream_dirs)
    profile.set("buckets_queued", len(buckets))
    if buckets:
        profile.set("largest_bucket_records", buckets[0].n)
        profile.set("largest_bucket_aa_key", buckets[0].key[1])
        profile.set("largest_bucket_subset_size", buckets[0].size_subset)
    run = Run(cg_automorphisms=cg_automorphisms, seq_sim_thresh=seq_sim_thresh,
              vdglib_dir=vdglib_dir, logfile=logfile,
              parent_pdb_dir=os.path.abspath(pdb_dir))
    failed_buckets = []
    if buckets:
        ctx = mp.get_context("spawn")
        try:
          with profile.phase("cluster"):
            # ProcessPoolExecutor rather than mp.Pool: a worker killed by the OOM
            # killer leaves mp.Pool hanging forever on imap_unordered -- it
            # respawns the worker but never reports the lost task -- so the job
            # burns its full h_rt and exits with a partial nr_vdgs/ and no
            # traceback. The executor raises BrokenProcessPool instead.
            pool = concurrent.futures.ProcessPoolExecutor(
                max_workers=int(args.num_procs), mp_context=ctx)
            try:
                failed_buckets = _cluster_buckets(
                    pool, int(args.num_procs), run, buckets, profile)
            finally:
                # cancel_futures matters: the executor's __exit__ would otherwise
                # run every queued task to completion before the exception
                # surfaced, so a systematic failure would still burn the whole
                # h_rt (in-flight tasks still finish).
                pool.shutdown(wait=True, cancel_futures=True)
        except BrokenProcessPool as e:
            err_text = (
                "[FATAL ERROR] A clustering worker died without returning a result "
                f"({e}). The usual cause is the OOM killer. Per-task peaks are the "
                "Kabsch batch of a Stage-1 row block (~0.3 GB) and the graph "
                "assembly of a large bucket (~30 B per within-cutoff pair); "
                "largest_bucket_* in the profile says which bucket to look at, "
                "and -l mem_free is PER SLOT under -pe smp.\n")
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
    if args.keep_stream_dir:
        _log_write(logfile, f"Streamed bucket columns kept at {stream_root}\n")
        print(f"Streamed bucket columns kept at {stream_root}")
    else:
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
