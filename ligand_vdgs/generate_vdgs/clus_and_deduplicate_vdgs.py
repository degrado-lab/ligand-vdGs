import os
import re
import sys
import argparse
import time
import numpy as np
import multiprocessing as mp
import concurrent.futures
from concurrent.futures.process import BrokenProcessPool
import json
import pickle
import traceback
import shutil
import hashlib
import fcntl
import errno
import itertools
import datetime
import subprocess
from collections import OrderedDict, deque, namedtuple

EXIT_NO_VDGS = 3

class _LRUCache:
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
    with open(logfile, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            f.write(msg)
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)

_WARN_COUNTS = {}
_SIDECAR_BUFFER = []
_SIDECAR_FLUSH_EVERY = 1000
_LAST_SIDECAR_LINE = {}

def _sidecar_path(logfile):
    return logfile + ".warn.tsv"

def _flush_sidecar(logfile):
    """Batch sidecar writes to avoid one flock'd NFS write per warning.
    Call at the end of any loop that emits warnings so the tail isn't lost."""
    if not _SIDECAR_BUFFER:
        return
    _log_write(_sidecar_path(logfile), "\n".join(_SIDECAR_BUFFER) + "\n")
    _SIDECAR_BUFFER.clear()

def _sidecar_log(logfile, key, *fields):
    """Buffer one compact tab-separated record for the logfile's sidecar.
    First field is the category key; the rest are category-specific and
    not fixed-schema. Skips writes identical to the immediately preceding
    line for this logfile and reports whether a line was written."""
    line = "\t".join((key,) + tuple(
        "" if f is None else str(f).replace("\t", " ").replace("\n", "\\n")
        for f in fields))
    if _LAST_SIDECAR_LINE.get(logfile) == line:
        return False
    _LAST_SIDECAR_LINE[logfile] = line
    _SIDECAR_BUFFER.append(line)
    if len(_SIDECAR_BUFFER) >= _SIDECAR_FLUSH_EVERY:
        _flush_sidecar(logfile)
    return True

def _log_warn(logfile, key, *fields):
    """Buffer a warning's detail for the sidecar, and count it (for the
    end-of-streaming summary line in the main log) only when it's not a
    duplicate of the immediately preceding warning, so counts reflect
    unique residues rather than raw call counts."""
    if _sidecar_log(logfile, key, *fields):
        _WARN_COUNTS[key] = _WARN_COUNTS.get(key, 0) + 1

def add_vdg_miner_paths():
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
try:
    import cg
except Exception:
    cg = None

from ligand_vdgs.functions import align_and_cluster as clust
from ligand_vdgs.functions import ligand_perception, parent_db
from ligand_vdgs.functions.clus_helpers import (
    get_vdg_subsets_target_size, select_diverse_pdbIDs, _aa_tmp_dir,
    _stream_root, stage1_record_is_complete, records_to_columns, concat_columns,
    save_shard, load_shard, save_columns, load_columns, load_stage1, COLUMNS_FILE)
from ligand_vdgs.functions.vdg_struct_utils import (VDM_OCC, cg_slot_occupancy,
    get_cg_atoms, is_hydrogen)
from ligand_vdgs.functions import vdg_npz_utils
from ligand_vdgs.functions.vdg_npz_utils import (write_cg_symmetry, name_selstr,
    parse_pdb_with_retry)
from ligand_vdgs.functions.vdg_fp_utils import build_perm_group, slot_orders
from ligand_vdgs.functions.align_and_cluster import (
    butina_partition, neighbor_csr_from_chunks, stage1_edges)
from ligand_vdgs.functions.dock_utils import cg_element_symbols
from ligand_vdgs.functions.compute_profile import ComputeProfile
from ligand_vdgs.functions.utils import (convert_time_elapsed, normalize_rmsd,
    identify_mol_automorphisms, mol_from_fragment, validate_atom_permutations,
    _int_or_none,)

MAX_CG_BOND_DIST = 2.5

# Cordero covalent radii (A) for every element in the current fragment
# vocabulary. The broad ratios are a provisional, bond-class-independent
# fallback; class/order-specific calibration is deferred to the next refactor.
_CG_COVALENT_RADII = {
    'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39}
MIN_CG_BOND_RADIUS_RATIO = 0.70
MAX_CG_BOND_RADIUS_RATIO = 1.25

def _cg_smarts_bond_class(bond):
    """Preserve exact SMARTS bond order where possible, else its ambiguity."""
    if bond.GetIsAromatic() or str(bond.GetBondType()) == 'AROMATIC':
        return 'aromatic'
    description = ' '.join(bond.DescribeQuery().split()) if bond.HasQuery() else ''
    exact_query_orders = {
        'BondOrder 1 = val': 'single',
        'BondOrder 2 = val': 'double',
        'BondOrder 3 = val': 'triple'}
    if description in exact_query_orders:
        return exact_query_orders[description]
    if description == 'SingleOrAromaticBond 1 = val':
        return 'single_or_aromatic'
    if not bond.HasQuery():
        return {'SINGLE': 'single', 'DOUBLE': 'double', 'TRIPLE': 'triple'}.get(
            str(bond.GetBondType()), 'ambiguous')
    return 'ambiguous'

def _cg_bond_length_violations(coords, elements, bonds):
    """SMARTS edges outside provisional broad covalent-radius bounds.

    The invariant for every query edge i-j is
    ``0.70 <= distance(i, j) / (r_cov(i) + r_cov(j)) <= 1.25``.
    Bond class/order is retained for diagnostics and future calibration. The
    same broad envelope is deliberately used for ambiguous bonds; unknown
    elements remain violations rather than silently disabling the gate.
    """
    coords = np.asarray(coords, dtype=float)
    violations = []
    for i, j, bond_class in bonds:
        elem_i, elem_j = str(elements[i]), str(elements[j])
        radii = (_CG_COVALENT_RADII.get(elem_i), _CG_COVALENT_RADII.get(elem_j))
        if None in radii:
            violations.append((i, j, bond_class, float('nan'), float('nan'),
                               float('nan'), f'unsupported_element:{elem_i}-{elem_j}'))
            continue
        radius_sum = sum(radii)
        lower = MIN_CG_BOND_RADIUS_RATIO * radius_sum
        upper = MAX_CG_BOND_RADIUS_RATIO * radius_sum
        distance = float(np.linalg.norm(coords[i] - coords[j]))
        if not np.isfinite(distance) or not lower <= distance <= upper:
            violations.append((i, j, bond_class, distance, lower, upper, 'length'))
    return violations

def _cg_bond_components(coords, cutoff=MAX_CG_BOND_DIST):
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

STREAM_CHUNKS_PER_WORKER = 4

_SPLIT_MIN_RECORDS = 8_000
_BLOCKS_PER_PROC = 4
_PIVOT_ROWS_PER_TASK = 2_000

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
    if args.num_procs < 1:
        parser.error(f"--num-procs must be at least 1, got {args.num_procs}")
    return args

_MAX_BUCKET_FNAME_BYTES = 200

def _bucket_fname(aa_key, flush_idx=None):
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
    match = _BUCKET_FNAME_RE.match(fname)
    return match.group("key") if match else None

def _bucket_npz_path(vdglib_dir, size_subset, sign, reordered_AAs):
    aa_label = "_".join(reordered_AAs)
    size_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset), sign)
    os.makedirs(size_dir, exist_ok=True)
    return os.path.join(size_dir, f"{aa_label}.npz")

def _write_failed_marker(vdglib_dir, size_subset, sign, aa_label, err_text):
    size_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset), sign)
    os.makedirs(size_dir, exist_ok=True)
    with open(os.path.join(size_dir, f"{aa_label}.FAILED"), "w") as fh:
        fh.write(err_text + "\n")

Subgroup = namedtuple("Subgroup", "stage1_id stage2_id nr_idx member_idxs radius")

_QUALITY_FIELDS = ("cg_max_b", "cg_min_occ", "vdm_max_b", "vdm_min_occ")

_SCHEMA_JSON_WIDTH = 2048

_GIT_REV_CACHE = {}

def _git_rev(path):
    directory = os.path.dirname(path) or "."
    if directory in _GIT_REV_CACHE:
        return _GIT_REV_CACHE[directory]
    _GIT_REV_CACHE[directory] = _git_rev_uncached(directory)
    return _GIT_REV_CACHE[directory]

def _git_rev_uncached(directory):
    try:
        out = subprocess.run(["git", "-C", directory,
                              "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""

def _bucket_schema_provenance(parent_pdb_dir):
    return json.dumps({
        "schema_version": vdg_npz_utils.BUCKET_SCHEMA_VERSION,
        "build_date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "parent_pdb_dir": str(parent_pdb_dir or ""),
        "annotation_schema": vdg_npz_utils.ANNOTATION_SCHEMA,
        "ligand_vdgs_rev": _git_rev(__file__),
        "vdg_miner_rev": _git_rev(cg.__file__) if cg is not None else "",
    }, sort_keys=True)

def _write_bucket_npz(vdglib_dir, size_subset, sign, reordered_AAs, cols, subgroups,
                      parent_pdb_dir):
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
    entries = np.asarray([parent_db.entry_of(b) for b in cols["biounit"]])
    num_parents = np.asarray([np.unique(entries[g.member_idxs]).size
                              for g in subgroups], dtype=np.int32)
    cgvdmbb = np.asarray(cols["cgvdmbb"][nr], dtype=np.float32)

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
        "aa_bucket_parts": np.asarray(list(reordered_AAs), dtype="U4"),
        "charge_sign": np.asarray(str(sign), dtype="U11"),
        "cluster_id": np.arange(1, C + 1, dtype=np.int32),
        "cluster_size": sizes,
        "cluster_pose_radius": np.asarray([g.radius for g in subgroups], dtype=np.float32),
        "cluster_num_parents": num_parents,
        "first_stage_cluster_id": np.asarray([g.stage1_id for g in subgroups], dtype=np.int32),
        "second_stage_cluster_id": np.asarray([g.stage2_id for g in subgroups], dtype=np.int32),
        "nr_cg_coords": cgvdmbb[:, :n_cg],
        "nr_vdm_bb_coords": cgvdmbb[:, n_cg:].reshape(C, num_vdms, 3, 3),
        "nr_vdm_o_coords": np.asarray(cols["vdm_o"][nr], dtype=np.float32),
        "parent_pdb_dir": np.asarray(str(parent_pdb_dir or ""), dtype="U512"),
        "nr_parent_biounit": cols["biounit"][nr],
        "nr_scrr_seg": cols["scrr_seg"][nr],
        "nr_scrr_chain": cols["scrr_chain"][nr],
        "nr_scrr_resnum": cols["scrr_resnum"][nr],
        "nr_scrr_resname": cols["scrr_resname"][nr],
        "nr_cg_names": cols["cg_names"][nr],
        "cg_elements": resolved,
        "nr_cg_seg": cols["cg_seg"][nr],
        "nr_cg_chain": cols["cg_chain"][nr],
        "nr_cg_resnum": cols["cg_resnum"][nr],
        "nr_cg_resname": cols["cg_resname"][nr],
        "nr_slot_flag": cols["slot_flags"][nr],
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
    for key, col in (("cg_heavy_degree", "cg_heavy_degree"),
                     ("cg_num_h", "cg_num_h"),
                     ("cg_placed_h", "cg_placed_h"),
                     ("cg_formal_charge", "cg_formal_charge"),
                     ("cg_nbr_elems", "cg_nbr_elems"),
                     ("perception", "perception"),
                     ("vdm_buried_area", "vdm_buried_area"),
                     ("vdm_shared_area", "vdm_shared_area"),
                     ("vdm_n_atom_pairs", "vdm_n_atom_pairs"),
                     ("vdm_min_heavy_dist", "vdm_min_heavy_dist")):
        arrays[f"nr_{key}"] = cols[col][nr]
        arrays[f"mem_{key}"] = cols[col][mem]

    bad_degree = (cols["cg_heavy_degree"][np.concatenate([nr, mem])] == 0)
    if bad_degree.any():
        raise ValueError(
            f"bucket {'_'.join(reordered_AAs)}: {int(bad_degree.sum())} CG atoms "
            "carry heavy degree 0, which no atom of a connected fragment can "
            "have; use -1 for unreadable.")

    arrays["schema"] = np.asarray(_bucket_schema_provenance(parent_pdb_dir),
                                  dtype=f"U{_SCHEMA_JSON_WIDTH}")
    if len(str(arrays["schema"])) >= _SCHEMA_JSON_WIDTH:
        raise ValueError("bucket provenance JSON reaches the dtype width and "
                         "would be truncated; widen _SCHEMA_JSON_WIDTH.")
    for k, key in enumerate(_QUALITY_FIELDS):
        arrays[f"nr_{key}"] = cols["quality"][nr, k]
        arrays[f"mem_{key}"] = cols["quality"][mem, k]

    if len(str(parent_pdb_dir or "")) >= arrays["parent_pdb_dir"].dtype.itemsize // 4:
        raise ValueError(
            f"parent_pdb_dir {parent_pdb_dir!r} reaches the dtype width and may be "
            f"truncated; widen the dtype.")

    path = _bucket_npz_path(vdglib_dir, size_subset, sign, reordered_AAs)
    tmp = f"{path}.{os.getpid()}.tmp"
    try:
        with open(tmp, 'wb') as handle:
            np.savez_compressed(handle, **arrays)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise

def _abort_profile(profile, profile_json, prefix, reason):
    profile.set("aborted_at", reason)
    profile.record_peak_rss(prefix)
    profile.write(profile_json)

Bucket = namedtuple("Bucket", "key size_subset sign aa_parts n bucket_dir")
Run = namedtuple("Run", "cg_automorphisms seq_sim_thresh vdglib_dir logfile parent_pdb_dir")

def _bucket_geometry(run, bucket, cols):
    n_cg = cols["cg_names"].shape[1]
    n_total = cols["cgvdmbb"].shape[1]
    cutoff = normalize_rmsd(n_total, "cgvdmbb")
    perm_group = build_perm_group(run.cg_automorphisms, n_cg, bucket.aa_parts)
    return n_cg, cutoff, perm_group

def _rank_clusters(clusters):
    return sorted(clusters, key=lambda members: -len(members))

def _row_blocks(n, n_blocks):
    cum = np.cumsum(np.arange(n - 1, 0, -1), dtype=np.int64)
    targets = cum[-1] * np.arange(1, n_blocks) / n_blocks
    cuts = np.searchsorted(cum, targets) + 1
    bounds = np.unique(np.concatenate([[0], cuts, [n - 1]]))
    return [(int(a), int(b)) for a, b in zip(bounds[:-1], bounds[1:])]

def _stage2_blocks(sizes, n_blocks):
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
    data = np.asarray(cols["cgvdmbb"], dtype=np.float32)
    flank_seq, flank_ca = cols["flank_seq"], cols["flank_ca"]
    orders = slot_orders(bucket.aa_parts)
    n_flank = flank_ca.shape[1]
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
            centre, radius = clust.pose_minimax_prototype(
                data, global_idxs, n_cg, perm_group)
            subgroups.append(Subgroup(cid, int(sub_num), int(centre), global_idxs, radius))
        if len(subgroups) == n_before:
            _log_write(run.logfile,
                f"[WARNING] Stage 1 cluster {cid} in {bucket.sign}/{bucket.key[2]} lost "
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
    _write_bucket_npz(run.vdglib_dir, bucket.size_subset, bucket.sign, bucket.aa_parts,
                      cols, subgroups, run.parent_pdb_dir)
    stats["write_s"] = time.time() - t0
    return sum(len(g.member_idxs) for g in subgroups), stats

def _edges_path(bucket, k):
    return os.path.join(bucket.bucket_dir, f"edges_{k:04d}.npy")

def _pivots_path(bucket):
    return os.path.join(bucket.bucket_dir, "pivots.npy")

def _even_blocks(n, k):
    """Equal row spans; the pivot embedding costs the same per row, unlike `_row_blocks`."""
    cuts = np.linspace(0, n, min(max(k, 1), n) + 1).astype(np.int64)
    return [(int(a), int(b)) for a, b in zip(cuts[:-1], cuts[1:]) if b > a]

def _load_stage1_cols(bucket):
    with np.load(os.path.join(bucket.bucket_dir, COLUMNS_FILE)) as columns:
        return {"cgvdmbb": np.asarray(load_stage1(bucket.bucket_dir, mmap=True)),
                "cg_names": columns["cg_names"]}

def _task_stage1_pivots(run, bucket, _payload):
    """Pick the bucket's pivot rows and preallocate the (k, n) embedding on disk."""
    cols = _load_stage1_cols(bucket)
    n_cg, cutoff, perm_group = _bucket_geometry(run, bucket, cols)
    clust.assert_perm_group_closed(perm_group)
    _t0 = time.time()
    ids = clust.select_pivots(cols["cgvdmbb"], n_cg, perm_group,
                              clust.pivot_count(bucket.n), cutoff)
    np.lib.format.open_memmap(_pivots_path(bucket), mode="w+", dtype=np.float32,
                              shape=(len(ids), bucket.n))
    return ids, {"pivot_select_s": time.time() - _t0, "pivot_k": int(len(ids))}

def _task_stage1_pivot_block(run, bucket, payload):
    """Fill one disjoint column span of the shared pivot memmap."""
    ids, start, stop = payload
    cols = _load_stage1_cols(bucket)
    n_cg, _cutoff, perm_group = _bucket_geometry(run, bucket, cols)
    _t0 = time.time()
    out = np.lib.format.open_memmap(_pivots_path(bucket), mode="r+")
    clust.pivot_distances(cols["cgvdmbb"], ids, n_cg, perm_group,
                          row_range=(start, stop), out=out[:, start:stop])
    out.flush()
    return {"pivot_block_wall_s": time.time() - _t0}

def _task_stage1_block(run, bucket, payload):
    k, start, stop = payload
    cols = _load_stage1_cols(bucket)
    n_cg, cutoff, perm_group = _bucket_geometry(run, bucket, cols)
    stats = {}
    t0 = time.time()
    qi, qj = stage1_edges(cols["cgvdmbb"], cutoff, n_cg, perm_group,
                          row_range=(start, stop), counters=stats,
                          pivots=np.load(_pivots_path(bucket), mmap_mode="r"))
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
    t0 = time.time()

    def _chunks():
        for k in range(n_blocks):
            part = np.load(_edges_path(bucket, k))
            yield part[0], part[1]

    clusters = _rank_clusters(butina_partition(
        *neighbor_csr_from_chunks(bucket.n, _chunks)))
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
    os.remove(_pivots_path(bucket))
    return sizes, {"stage1_partition_s": time.time() - t0,
                   "stage1_clusters": int(sizes.size)}

def _load_partition(bucket):
    with np.load(_partition_path(bucket)) as data:
        order, indptr = data["order"], data["indptr"]
    return [order[indptr[c]:indptr[c + 1]] for c in range(indptr.size - 1)]

def _task_stage2_block(run, bucket, cluster_ids):
    with np.load(os.path.join(bucket.bucket_dir, COLUMNS_FILE)) as columns:
        cols = {name: columns[name]
                for name in ("cg_names", "flank_seq", "flank_ca")}
    cols["cgvdmbb"] = load_stage1(bucket.bucket_dir, mmap=True)
    n_cg, _cutoff, perm_group = _bucket_geometry(run, bucket, cols)
    clusters = _load_partition(bucket)
    return _stage2_for_clusters(run, bucket, cols, clusters, cluster_ids, perm_group, n_cg)

def _task_write(run, bucket, subgroups):
    t0 = time.time()
    cols = load_columns(bucket.bucket_dir)
    _write_bucket_npz(run.vdglib_dir, bucket.size_subset, bucket.sign, bucket.aa_parts,
                      cols, subgroups, run.parent_pdb_dir)
    os.remove(_partition_path(bucket))
    return sum(len(g.member_idxs) for g in subgroups), {"write_s": time.time() - t0}

_TASKS = {"small": _task_small, "stage1_pivots": _task_stage1_pivots,
          "stage1_pivot_block": _task_stage1_pivot_block,
          "stage1_block": _task_stage1_block,
          "stage1_partition": _task_stage1_partition,
          "stage2_block": _task_stage2_block, "write": _task_write}

def _run_task(kind, run, bucket, payload):
    t0 = time.time()
    try:
        return kind, bucket.key, _TASKS[kind](run, bucket, payload), None, t0, time.time()
    except Exception as exc:
        return kind, bucket.key, None, (
            f"[WORKER ERROR] AA bucket: {bucket.sign}/{bucket.key[2]} (subset size "
            f"{bucket.size_subset}, task {kind})\nException: {exc}\n"
            f"Traceback:\n{traceback.format_exc()}"), t0, time.time()

def _parse_pdb_with_retry(pdb_file, attempts=3, delay=0.5):
    return parse_pdb_with_retry(pdb_file, attempts=attempts, delay=delay)

_TRANSIENT_ERRNOS = frozenset(
    e for e in (errno.EAGAIN, errno.EBUSY, errno.EINTR, errno.EIO, errno.EMFILE,
                errno.ENFILE, errno.ENOMEM, errno.ESTALE, errno.ETIMEDOUT,
                getattr(errno, "EREMOTEIO", None))
    if e is not None)

_PARSE_FAILED = object()

def _is_transient_parse_failure(exc):
    if isinstance(exc, OSError):
        return exc.errno in _TRANSIENT_ERRNOS
    return False

def _get_atomgroup_for_env(
    environment, pdb_dir, cg, cg_match_dict, align_atoms, logfile,
    pdb_cache=None, cg_bonds=(), expected_elements=(), row_rejections=None):
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
            _log_warn(logfile, 'pdb_parse_failed', biounit, e)
            return None
        if pdb_cache is not None:
            pdb_cache[pdb_file] = whole_struct

    if whole_struct is None:
        _log_warn(logfile, 'pdb_returned_none', biounit)
        return None

    raw_resnums = [tup[3] for tup in environment]
    scrs = [(tup[1], tup[2], f"`{tup[3]}`" if tup[3] < 0 else tup[3]) for tup in environment]
    selstrs = [
        f"(segment {scr[0]} and chain {scr[1]} and resnum {scr[2]})" if scr[0]
        else f"(chain {scr[1]} and resnum {scr[2]})"
        for scr in scrs]
    if len(selstrs) < 2:
        _log_warn(logfile, 'env_too_few_entries', biounit)
        return None
    try:
        sel = whole_struct.select(
            "same residue as within 5 of ({})".format(" or ".join(selstrs[1:])))
    except Exception:
        _log_warn(logfile, 'prody_selection_failed', biounit)
        return None
    if sel is None:
        return None
    struct = sel.toAtomGroup()
    resnames = []
    align_coords = np.zeros((3, 3))
    cg_atom_coords = None

    for i, (scr, selstr) in enumerate(zip(scrs, selstrs)):
        try:
            substruct = struct.select(selstr)
            if substruct is None:
                return None
            resnames.append(substruct.getResnames()[0])

            unique_res_indices = np.unique(substruct.getResindices())
            if len(unique_res_indices) != 1:
                _log_warn(logfile, 'ambiguous_residue', biounit, scr[1], scr[2])
                return None

            if i == 0:
                if cg in cg_atoms.keys():
                    atom_names_list = cg_atoms[cg][resnames[0]]
                else:
                    key = (biounit, scrs[0][0], scrs[0][1],
                           str(raw_resnums[0]), resnames[0])

                    match_list = cg_match_dict.get(key)
                    match_idx = environment[0][4] - 1

                    if match_list is None:
                        _log_warn(logfile, 'no_cg_match', biounit, scr[1], scr[2], key)
                        return None

                    if not (0 <= match_idx < len(match_list)):
                        _log_warn(logfile, 'match_idx_out_of_range', biounit, scr[1],
                                  scr[2], match_idx, len(match_list))
                        return None
                    atom_names_list = match_list[match_idx]

                cg_atom_selstrs = ["name " + name_selstr(atom_name)
                                   for atom_name in atom_names_list]

                chosen_atoms = []
                for cg_selstr, atom_name in zip(cg_atom_selstrs, atom_names_list):
                    atom_sel = substruct.select(cg_selstr)
                    if atom_sel is None or atom_sel.numAtoms() == 0:
                        return None
                    if atom_sel.numAtoms() > 1:
                        _log_warn(logfile, 'duplicate_cg_atom_name', biounit,
                                  scrs[0][1], scrs[0][2], resnames[0], atom_name,
                                  atom_sel.numAtoms())
                        return None
                    chosen_atoms.append(atom_sel[0])

                cg_atom_coords = np.asarray(
                    [np.reshape(a.getCoords(), 3) for a in chosen_atoms], dtype=float)

                if cg_bonds:
                    observed_elements = [str(a.getElement()).strip().capitalize()
                                         for a in chosen_atoms]
                    bond_elements = [
                        expected_elements[i]
                        if i < len(expected_elements) and expected_elements[i] is not None
                        else observed
                        for i, observed in enumerate(observed_elements)]
                    violations = _cg_bond_length_violations(
                        cg_atom_coords, bond_elements, cg_bonds)
                    if violations:
                        details = ';'.join(
                            f'{atom_names_list[a]}-{atom_names_list[b]}:{bond_class}:'
                            f'{distance:.2f}:{lower:.2f}-{upper:.2f}:{reason}'
                            for a, b, bond_class, distance, lower, upper, reason
                            in violations)
                        _log_warn(
                            logfile, 'cg_bond_length_outlier', biounit,
                            scrs[0][1], scrs[0][2], resnames[0], details)
                        if row_rejections is not None:
                            row_rejections['cg_bond_geometry'] = (
                                row_rejections.get('cg_bond_geometry', 0) + 1)
                        return None

                components = _cg_bond_components(cg_atom_coords)
                if len(components) > 1:
                    groups = ' | '.join(
                        '+'.join(atom_names_list[k] for k in comp)
                        for comp in components)
                    gap = min(
                        float(np.linalg.norm(cg_atom_coords[k] - cg_atom_coords[l]))
                        for a, b in itertools.combinations(components, 2)
                        for k in a for l in b)
                    _log_warn(logfile, 'cg_not_bonded_component', biounit,
                              scrs[0][1], scrs[0][2], resnames[0], groups,
                              f"{gap:.2f}")
                    return None

                for j, chosen_atom in enumerate(chosen_atoms):
                    chosen_atom.setOccupancy(cg_slot_occupancy(j))

                    if j in align_atoms:
                        c = np.asarray(chosen_atom.getCoords())
                        c = c[0] if c.ndim == 2 else c
                        align_coords[align_atoms.index(j)] = c

            else:
                _names = substruct.getNames()
                _els = substruct.getElements()
                if _els is None:
                    _els = [''] * len(_names)
                _heavy_mask = np.array(
                    [not is_hydrogen(n, e) for n, e in zip(_names, _els)], dtype=bool)
                if not _heavy_mask.any():
                    _log_warn(logfile, 'vdm_all_hydrogen', biounit, scr[1], scr[2])
                    continue
                substruct.setOccupancies(VDM_OCC)

        except Exception:
            _log_warn(logfile, 'residue_selection_exception', biounit, scr[1], scr[2],
                      traceback.format_exc())
            return None

    if not align_coords_sanity_check(align_coords):
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

_PDB_CACHE_ATOM_BUDGET = 10_000_000

def _struct_atom_count(struct):
    try:
        return int(struct.numAtoms())
    except Exception:
        return 0

def _flush_buckets_to_disk(sign_records, worker_tmp_dir, flush_idx=None):
    for aa_key, records in sign_records.items():
        fname = _bucket_fname(aa_key, flush_idx=flush_idx)
        save_shard(os.path.join(worker_tmp_dir, fname), records_to_columns(records))

def _stream_one_chunk(args):
    (environment_files_chunk, environments_dir, pdb_dir, CG, cg_match_dict_pkl,
     align_atoms, logfile, expected_n_cg, expected_element_seq, cg_bonds,
     num_flanking, worker_tmp_root, pdb_cache_limits, subset_sizes,
     flush_threshold) = args

    if cg_match_dict_pkl:
        with open(cg_match_dict_pkl, "rb") as _fh:
            cg_match_dict = pickle.load(_fh)
        _annot_path = vdg_npz_utils.cg_annot_pkl_path(cg_match_dict_pkl)
        if not os.path.isfile(_annot_path):
            raise FileNotFoundError(
                f"{_annot_path} is missing beside {cg_match_dict_pkl}; re-run "
                "smarts_to_cgs.py, which writes the per-CG-atom annotations "
                "the bucket npz requires.")
        with open(_annot_path, "rb") as _fh:
            cg_annot_dict = pickle.load(_fh)
    else:
        cg_match_dict = {}
        cg_annot_dict = {}
    pdb_cache_size, pdb_cache_atom_budget = pdb_cache_limits
    pdb_cache = (_LRUCache(maxsize=pdb_cache_size,
                           max_weight=pdb_cache_atom_budget,
                           weigh=_struct_atom_count)
                 if pdb_cache_size > 0 else None)
    bucket_records = {(size, sign): {} for size in subset_sizes
                      for sign in vdg_npz_utils.CHARGE_SIGNS}
    worker_tmp_dirs = {
        (size, sign): os.path.join(worker_tmp_root, str(size), sign)
        for size in subset_sizes for sign in vdg_npz_utils.CHARGE_SIGNS}
    for worker_tmp_dir in worker_tmp_dirs.values():
        os.makedirs(worker_tmp_dir, exist_ok=True)
    total_in_memory, flush_idx = 0, 0
    skips = {"no_atomgroup": 0, "dup_occupancy": 0, "no_cg_atoms": 0,
             "bad_cg_coords": 0, "cg_elements_mismatch": 0, "unparsable_line": 0,
             "incomplete_stage1": 0, "cg_bond_geometry": 0}
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

                environment = record["env"]
                env_quality = (record["cg_max_b"], record["cg_min_occ"],
                               record["vdm_max_b"], record["vdm_min_occ"])

                try:
                    contact_by_res = _contact_strength_by_residue(
                        record, environment)
                except ValueError as _e:
                    raise ValueError(
                        f"{env_path}: {_e}. Contact strength is required at "
                        "build time and is never filled in.")

                pdb_label = "_".join([str(el) for el in environment[0]])
                _biounit = environment[0][0]

                geometry_rejections_before = skips['cg_bond_geometry']
                atomgroup = _get_atomgroup_for_env(
                    environment, pdb_dir, CG, cg_match_dict, align_atoms, logfile,
                    pdb_cache=pdb_cache, cg_bonds=cg_bonds,
                    expected_elements=expected_element_seq, row_rejections=skips)
                if atomgroup is None:
                    if skips['cg_bond_geometry'] == geometry_rejections_before:
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
                (cg_coords, cg_names, cg_elements,
                 cg_seg, cg_chain, cg_resnum, cg_resname) = _cg_all

                if (cg_coords.shape != (expected_n_cg, 3)
                        or not np.isfinite(cg_coords).all()):
                    skips["bad_cg_coords"] += 1
                    continue

                _elements = tuple(str(e).strip().capitalize() for e in cg_elements)
                _bad_elements = any(
                    a and b and a != b
                    for a, b in zip(_elements, expected_element_seq))
                if _bad_elements:
                    skips["cg_elements_mismatch"] += 1
                    _sidecar_log(logfile, 'cg_elements_mismatch', pdb_label,
                                 ",".join(_elements), ",".join(expected_element_seq))
                    continue

                _cg_key = (str(_biounit), str(environment[0][1]),
                           str(environment[0][2]), str(environment[0][3]),
                           str(cg_resname))
                if CG in cg_atoms:
                    _n_cg = len(cg_names)
                    _cg_annot = {
                        "heavy_degree": [ANNOT_UNREADABLE] * _n_cg,
                        "num_h": [ANNOT_UNREADABLE] * _n_cg,
                        "formal_charge": [ANNOT_UNREADABLE] * _n_cg,
                        "nbr_elems": [""] * _n_cg,
                        "perception": ligand_perception.PERCEPTION_ATOM_NAME_TABLE,
                    }
                else:
                    _cg_annot = _cg_annotations_in_slot_order(
                        cg_annot_dict, cg_match_dict, _cg_key,
                        environment[0][4] - 1, cg_names, pdb_label)
                _sign = _charge_sign(_cg_annot)
                _cg_placed_h = _cg_placed_h_in_slot_order(
                    atomgroup, cg_coords, cg_seg, cg_chain, cg_resnum)

                vdms_dict = clust.get_vdm_res_features(atomgroup, pdb_label, num_flanking)
                vdm_resinds = list(vdms_dict.keys())

                _heavy_mask = np.array(
                    [not is_hydrogen(n, e) for n, e in zip(cg_names, cg_elements)],
                    dtype=bool)
                if _heavy_mask.any():
                    _cg_heavy = cg_coords[_heavy_mask]
                else:
                    _log_warn(logfile, 'cg_all_hydrogen_fallback', pdb_label)
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
                            _log_warn(logfile, 'reorder_vdg_subset_failed', pdb_label, _e)
                            continue
                        _slot_keys = [(str(_s), str(_ch), int(_r))
                                      for (_s, _ch, _r, _rn) in re_ordered_scrr]
                        _absent = [_k for _k in _slot_keys
                                   if _k not in contact_by_res]
                        if _absent:
                            raise KeyError(
                                f"{pdb_label}: vdM slots {_absent} have no "
                                "contact strength in the environment record; "
                                "the gate and the assembled environment "
                                "disagree about membership.")
                        rec = {
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
                            "cg_heavy_degree": _cg_annot["heavy_degree"],
                            "cg_num_h": _cg_annot["num_h"],
                            "cg_placed_h": [int(x) for x in _cg_placed_h],
                            "cg_formal_charge": _cg_annot["formal_charge"],
                            "cg_nbr_elems": _cg_annot["nbr_elems"],
                            "perception": _cg_annot["perception"],
                            "buried_area": [contact_by_res[_k][0]
                                            for _k in _slot_keys],
                            "shared_area": [contact_by_res[_k][1]
                                            for _k in _slot_keys],
                            "n_atom_pairs": [contact_by_res[_k][2]
                                             for _k in _slot_keys],
                            "min_heavy_dist": [contact_by_res[_k][3]
                                               for _k in _slot_keys],
                            "bbo": [np.asarray(x, dtype=np.float32) for x in re_ordered_bbo],
                        }
                        if not stage1_record_is_complete(rec, expected_n_cg, size_subset):
                            skips["incomplete_stage1"] += 1
                            continue
                        aa_key = "_".join(re_ordered_aas)
                        sign_records = bucket_records[(size_subset, _sign)]
                        sign_records.setdefault(aa_key, []).append(rec)
                        total_in_memory += 1
                        if total_in_memory >= flush_threshold:
                            for key, records in bucket_records.items():
                                _flush_buckets_to_disk(
                                    records, worker_tmp_dirs[key], flush_idx)
                            bucket_records = {(size, sign): {} for size in subset_sizes
                                              for sign in vdg_npz_utils.CHARGE_SIGNS}
                            flush_idx += 1
                            total_in_memory = 0

    for key, sign_records in bucket_records.items():
        _flush_buckets_to_disk(sign_records, worker_tmp_dirs[key])

    warns = {key: count - warn_counts_at_entry.get(key, 0)
             for key, count in _WARN_COUNTS.items()
             if count - warn_counts_at_entry.get(key, 0) > 0}
    _flush_sidecar(logfile)
    return worker_tmp_root, skips, warns

ANNOT_UNREADABLE = -1

def _charge_sign(cg_annot):
    if any(d == ANNOT_UNREADABLE for d in cg_annot["heavy_degree"]):
        return "unreadable"
    net = sum(cg_annot["formal_charge"])
    return "pos" if net > 0 else "neg" if net < 0 else "neut"

def _cg_placed_h_in_slot_order(atomgroup, cg_coords, cg_seg, cg_chain, cg_resnum,
                                h_bond_cutoff=1.3):
    h_mask = (
        (atomgroup.getSegnames() == cg_seg) & (atomgroup.getChids() == cg_chain)
        & (atomgroup.getResnums() == cg_resnum) & np.array(
            [is_hydrogen(n, e) for n, e in zip(atomgroup.getNames(), atomgroup.getElements())],
            dtype=bool))
    if not h_mask.any():
        return np.zeros(len(cg_coords), dtype=np.int8)
    h_coords = atomgroup.getCoords()[h_mask]
    return np.array(
        [1 if (np.linalg.norm(h_coords - c, axis=1) <= h_bond_cutoff).any() else 0
         for c in cg_coords], dtype=np.int8)

_CONTACT_FIELDS = ("buried_area", "shared_area", "n_atom_pairs", "min_heavy_dist")

def _contact_strength_by_residue(record, environment):
    n_slots = len(environment) - 1
    for field in _CONTACT_FIELDS:
        if field not in record:
            raise ValueError(
                f"environment record has no {field!r}; contact strength is "
                "required at build time (rebuild-notes 3)")
        if len(record[field]) != n_slots:
            raise ValueError(
                f"{field!r} has {len(record[field])} values for {n_slots} vdM "
                "slots; the gate and the environment disagree")
    out = {}
    for j, tup in enumerate(environment[1:]):
        dist = float(record["min_heavy_dist"][j])
        if not np.isfinite(dist):
            raise ValueError(
                f"vdM slot {tuple(tup[1:4])} has non-finite min_heavy_dist "
                f"{dist}, so the gate never measured a residue it admitted")
        out[(str(tup[1]), str(tup[2]), int(tup[3]))] = (
            float(record["buried_area"][j]), float(record["shared_area"][j]),
            int(record["n_atom_pairs"][j]), dist)
    return out

def _cg_annotations_in_slot_order(cg_annot_dict, cg_match_dict, key, match_idx,
                                  cg_names, pdb_label):
    entry = (cg_annot_dict.get(key) or [None])[match_idx] \
        if cg_annot_dict.get(key) and match_idx < len(cg_annot_dict[key]) else None
    if entry is None:
        raise KeyError(
            f"{pdb_label}: no CG annotations for key {key} match {match_idx}; "
            "the annotation pickle is out of step with the matches pickle.")
    match_names = cg_match_dict[key][match_idx]
    position = {str(name): i for i, name in enumerate(match_names)}
    order = []
    for name in cg_names:
        if str(name) not in position:
            raise KeyError(
                f"{pdb_label}: CG atom {name!r} is not in the matched atom "
                f"names {list(match_names)}; annotation would be misassigned.")
        order.append(position[str(name)])
    return {
        "heavy_degree": [entry["heavy_degree"][i] for i in order],
        "num_h": [entry["num_h"][i] for i in order],
        "formal_charge": [entry["formal_charge"][i] for i in order],
        "nbr_elems": [entry["nbr_elems"][i] for i in order],
        "perception": int(entry["perception"]),
    }

def _cap_distinct_entries(cols, aa_key, max_num_to_clus, logfile):
    entries = [parent_db.entry_of(b) for b in cols["biounit"]]
    if len(set(entries)) <= max_num_to_clus:
        return cols
    keep = set(select_diverse_pdbIDs(entries, max_num_to_clus))
    rows = np.asarray([e in keep for e in entries])
    _log_write(logfile, f"\t{len(entries)} vdGs -> {int(rows.sum())} vdGs from "
                        f"{len(keep)} diverse PDB IDs for {aa_key}.\n")
    return {key: arr[rows] for key, arr in cols.items()}

def _merge_worker_dirs(worker_dirs, final_dir, max_num_to_clus, logfile):
    os.makedirs(final_dir, exist_ok=True)
    files_by_key = {}
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
    buckets = []
    for (size_subset, sign), counts in counts_by_size.items():
        for aa_key, n in counts.items():
            buckets.append(Bucket(
                key=(size_subset, sign, aa_key), size_subset=size_subset, sign=sign,
                aa_parts=tuple(aa_key.split("_")), n=n,
                bucket_dir=os.path.join(stream_dirs[(size_subset, sign)], aa_key)))
    buckets.sort(key=lambda b: (-b.n, b.key))
    return buckets

def _record_bucket_stats(profile, stats):
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

def _buckets_jsonl_path(run):
    return os.path.splitext(run.logfile)[0] + "_buckets.jsonl"

def _append_bucket_row(path, bucket, stats, split, blocks):
    row = {"aa_key": bucket.key[2], "sign": bucket.sign, "subset_size": bucket.size_subset,
           "records": bucket.n, "split": bool(split), "blocks": blocks}
    for name in ("perm_group", "stage1_wall_s", "stage1_clusters", "edges",
                 "stage1_partition_s", "stage2_wall_s", "stage2_slot_orders",
                 "stage2_largest_cluster", "stage2_largest_cluster_s", "write_s",
                 "block_wall_max_s", "block_wall_sum_s", "stage2_blocks",
                 "stage2_block_wall_max_s", "pivot", "internal_lb", "bb_lb",
                 "exact", "exact_rows", "scan_s", "exact_s",
                 "pivot_k", "pivot_select_s", "pivot_block_wall_max_s"):
        if name in stats:
            row[name] = stats[name]
    with open(path, "a") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")

def _occupancy(spans):
    events = sorted([(t, +1) for t, _ in spans] + [(t, -1) for _, t in spans])
    hist, live, prev = {}, 0, None
    for t, delta in events:
        if prev is not None and t > prev:
            hist[live] = hist.get(live, 0.0) + (t - prev)
        live += delta
        prev = t
    return hist

def _cluster_buckets(pool, num_procs, run, buckets, profile):
    n_blocks = _BLOCKS_PER_PROC * num_procs
    pending = deque()
    state = {}
    for bucket in buckets:
        if bucket.n >= _SPLIT_MIN_RECORDS:
            blocks = _row_blocks(bucket.n, n_blocks)
            state[bucket.key] = {"bucket": bucket, "blocks": len(blocks),
                                 "blocks_done": 0, "stats": {"records": bucket.n},
                                 "row_blocks": blocks, "pivot_blocks": 0,
                                 "pivot_done": 0, "subgroups": [],
                                 "stage2_blocks": 0, "stage2_done": 0}
            pending.append(("stage1_pivots", bucket, None))
        else:
            pending.append(("small", bucket, None))
    profile.set("buckets_split", len(state))

    failed, done, inflight = set(), set(), {}
    spans = []
    rows_path = _buckets_jsonl_path(run) if profile.enabled else None

    def _fail(bucket, err_text):
        failed.add(bucket.key)
        _log_write(run.logfile, err_text + "\n")
        _write_failed_marker(run.vdglib_dir, bucket.size_subset, bucket.sign,
                             bucket.key[2], err_text)

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
            _kind, key, value, err_text, t_start, t_end = fut.result()
            spans.append((t_start, t_end))
            if key in failed:
                continue
            if err_text is not None:
                _fail(bucket, err_text)
                continue
            if kind == "small":
                total, stats = value
                if rows_path:
                    _append_bucket_row(rows_path, bucket, stats, False, 1)
                _record_bucket_stats(profile, stats)
                done.add(key)
                continue
            st = state[key]
            if kind == "stage1_pivots":
                ids, stats = value
                st["stats"].update(stats)
                # Cap the fan-out so a bucket just over _SPLIT_MIN_RECORDS does not
                # pay per-task load overhead on a few hundred rows apiece.
                pivot_blocks = _even_blocks(
                    bucket.n, min(num_procs, max(1, bucket.n // _PIVOT_ROWS_PER_TASK)))
                st["pivot_blocks"] = len(pivot_blocks)
                _push_front([("stage1_pivot_block", bucket, (ids, a, b))
                             for a, b in pivot_blocks])
            elif kind == "stage1_pivot_block":
                st["pivot_done"] += 1
                st["stats"]["pivot_block_wall_max_s"] = max(
                    st["stats"].get("pivot_block_wall_max_s", 0.0),
                    value["pivot_block_wall_s"])
                if st["pivot_done"] == st["pivot_blocks"]:
                    _push_front([("stage1_block", bucket, (k, a, b))
                                 for k, (a, b) in enumerate(st["row_blocks"])])
            elif kind == "stage1_block":
                st["blocks_done"] += 1
                for name, v in value.items():
                    st["stats"][name] = st["stats"].get(name, 0) + v
                block_s = value.get("stage1_wall_s", 0.0)
                st["stats"]["block_wall_sum_s"] = (
                    st["stats"].get("block_wall_sum_s", 0.0) + block_s)
                st["stats"]["block_wall_max_s"] = max(
                    st["stats"].get("block_wall_max_s", 0.0), block_s)
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
                st["stats"]["stage2_block_wall_max_s"] = max(
                    st["stats"].get("stage2_block_wall_max_s", 0.0),
                    stats.get("stage2_wall_s", 0.0))
                if st["stage2_done"] == st["stage2_blocks"]:
                    _push_front([("write", bucket, st["subgroups"])])
                    st["subgroups"] = []
            elif kind == "write":
                total, stats = value
                st["stats"].update(stats)
                st["stats"]["stage2_blocks"] = st["stage2_blocks"]
                if rows_path:
                    _append_bucket_row(rows_path, bucket, st["stats"], True,
                                       st["blocks"])
                _record_bucket_stats(profile, st["stats"])
                done.add(key)
        _fill()
    if profile.enabled:
        profile.set("stage1.buckets", len(done))
        profile.set("stage2.buckets", len(done))
        hist = _occupancy(spans)
        profile.set("sched.task_wall_s", round(sum(e - s0 for s0, e in spans), 3))
        profile.set("sched.tasks", len(spans))
        profile.set("sched.span_s", round(
            max((e for _, e in spans), default=0.0)
            - min((s0 for s0, _ in spans), default=0.0), 3))
        profile.set("sched.max_concurrency", max(hist, default=0))
        for k in (1, 2):
            profile.set(f"sched.wall_s_at_most_{k}_busy",
                        round(sum(v for c, v in hist.items() if 0 < c <= k), 3))
        profile.set("sched.rows_jsonl", os.path.basename(_buckets_jsonl_path(run)))
    return sorted(f"{key[1]}/{key[2]}" for key in failed)

def _subset_output_counts(vdglib_dir, size_subset, logfile):
    num_nr_vdgs = num_inputs = 0
    unreadable = []
    for sign in vdg_npz_utils.CHARGE_SIGNS:
        nr_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset), sign)
        if not os.path.isdir(nr_dir):
            continue
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
                unreadable.append(f"{sign}/{fname[:-len('.npz')]}")
    return num_nr_vdgs, num_inputs, unreadable

def _scan_failed_markers(vdglib_dir, size_subset):
    found = []
    for sign in vdg_npz_utils.CHARGE_SIGNS:
        nr_dir = os.path.join(vdglib_dir, "nr_vdgs", str(size_subset), sign)
        if not os.path.isdir(nr_dir):
            continue
        found.extend(f"{sign}/{fname[:-len('.FAILED')]}" for fname in os.listdir(nr_dir)
                     if fname.endswith(".FAILED"))
    return sorted(found)

def _preexisting_bucket_outputs(vdglib_dir, subset_sizes):
    found = []
    for size in subset_sizes:
        base = os.path.join(vdglib_dir, "nr_vdgs", str(size))
        for nr_dir, label in ([(base, "")] +
                [(os.path.join(base, s), f"{s}/") for s in vdg_npz_utils.CHARGE_SIGNS]):
            if not os.path.isdir(nr_dir):
                continue
            found.extend((size, f"{label}{fname}") for fname in sorted(os.listdir(nr_dir))
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
    expected_element_seq = cg_element_symbols(args.cg_smarts)
    cg_bonds = tuple((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
                      _cg_smarts_bond_class(bond))
                     for bond in symmetry_mol.GetBonds())
    unsupported = sorted({e for e in expected_element_seq
                          if e is not None and e not in _CG_COVALENT_RADII})
    if unsupported:
        raise ValueError(
            f"[ERROR] --cg-smarts uses elements without calibrated covalent "
            f"radii: {unsupported}")
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
    if not all_environment_files:
        err_text = (f"[WARNING] No .jsonl environment shards under {environments_dir}; "
                    "generate_environments.py produced no environments for this CG. "
                    "nr_vdgs/ is EMPTY and this fragment is NOT complete, so this log "
                    "deliberately carries no completion marker.\n")
        _log_write(logfile, err_text)
        print(err_text, file=sys.stderr)
        sys.exit(EXIT_NO_VDGS)

    write_cg_symmetry(args.vdglib_dir, args.cg_smarts, cg_automorphisms)

    stream_dirs = {
        (size, sign): os.path.join(_aa_tmp_dir(vdglib_dir, size), sign)
        for size in subset_sizes for sign in vdg_npz_utils.CHARGE_SIGNS}
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

        flush_threshold = max(25_000, _FLUSH_RECORDS_THRESHOLD // n_workers)
        pdb_cache_atom_budget = max(
            1, _PDB_CACHE_ATOM_BUDGET // n_workers) if pdb_cache_size > 0 else 0
        profile.set("stream_flush_threshold", flush_threshold)
        chunk_args = [
            (chunk, environments_dir, pdb_dir, CG, cg_match_dict_pkl,
             align_atoms, logfile, len(cg_automorphisms[0]), expected_element_seq,
             cg_bonds, num_flanking, wdir,
             (pdb_cache_size, pdb_cache_atom_budget),
             subset_sizes, flush_threshold)
            for chunk, wdir in zip(chunks, worker_dirs)]

        profile.set("environment_shards", len(all_environment_files))
        profile.set("stream_workers", n_workers)
        profile.set("stream_chunks", len(chunks))
        stream_ctx = mp.get_context("spawn")
        stream_skips = {}
        stream_warns = {}
        try:
          with profile.phase("stream_environments"):
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

        if any(stream_skips.values()):
            _log_write(logfile, "[INFO] Environments skipped during streaming: "
                + ", ".join(f"{r}={n}" for r, n in sorted(stream_skips.items()) if n)
                + "\n")
        if any(stream_warns.values()):
            _log_write(logfile, "[INFO] Warnings during streaming (incl. suppressed): "
                + ", ".join(f"{r}={n}" for r, n in sorted(stream_warns.items()) if n)
                + "\n")
        for reason, count in stream_skips.items():
            profile.set(f"stream_skipped_{reason}", count)
        for reason, count in stream_warns.items():
            profile.set(f"stream_warned_{reason}", count)
        _streamed_any = any(
            fname.endswith(".npz")
            for wdir in worker_dirs
            for _r, _d, _files in os.walk(wdir)
            for fname in _files)
        if not _streamed_any:
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
            counts_by_size = {
                (size_subset, sign): _merge_worker_dirs(
                    [os.path.join(wdir, str(size_subset), sign) for wdir in worker_dirs],
                    stream_dirs[(size_subset, sign)], max_num_to_clus, logfile)
                for size_subset in subset_sizes for sign in vdg_npz_utils.CHARGE_SIGNS}
            profile.set("aa_sign_partitions_nonempty",
                        sum(1 for c in counts_by_size.values() if c))
            merged_mb = sum(profile.record_dir_size(
                f"scratch_merged_buckets_size_{size}_{sign}_mb", stream_dirs[(size, sign)])
                for size in subset_sizes for sign in vdg_npz_utils.CHARGE_SIGNS)
            profile.max("scratch_peak_mb", worker_shards_mb + merged_mb)
            shutil.rmtree(workers_root, ignore_errors=True)
        except Exception as e:
            err_text = (f"[FATAL ERROR] Failed to merge stream shards into AA "
                        f"buckets:\n{e}\n")
            _log_write(logfile, err_text)
            print(err_text, file=sys.stderr)
            _abort_profile(profile, args.profile_json, "merge.", "merge_worker_shards")
            shutil.rmtree(stream_root, ignore_errors=True)
            sys.exit(1)

    buckets = _buckets(counts_by_size, stream_dirs)
    profile.set("buckets_queued", len(buckets))
    if buckets:
        profile.set("largest_bucket_records", buckets[0].n)
        profile.set("largest_bucket_aa_key", buckets[0].key[2])
        profile.set("largest_bucket_sign", buckets[0].sign)
        profile.set("largest_bucket_subset_size", buckets[0].size_subset)
    run = Run(cg_automorphisms=cg_automorphisms, seq_sim_thresh=seq_sim_thresh,
              vdglib_dir=vdglib_dir, logfile=logfile,
              parent_pdb_dir=os.path.abspath(pdb_dir))
    failed_buckets = []
    if buckets:
        ctx = mp.get_context("spawn")
        try:
          with profile.phase("cluster"):
            pool = concurrent.futures.ProcessPoolExecutor(
                max_workers=int(args.num_procs), mp_context=ctx)
            try:
                failed_buckets = _cluster_buckets(
                    pool, int(args.num_procs), run, buckets, profile)
            finally:
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
            sys.exit(1)
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
    _marked = set(bad_buckets)
    bad_buckets.update(
        f"{label} (worker error)" for label in failed_buckets
        if not any(b.startswith(f"{size}/{label} ")
                   for size in subset_sizes for b in _marked))
    profile.set("buckets_failed", len(bad_buckets))
    profile.write(args.profile_json)

    if bad_buckets:
        msg = (f"[ERROR] {len(bad_buckets)} AA bucket(s) are missing or failed to "
               f"cluster; nr_vdgs/ is INCOMPLETE: "
               f"{', '.join(sorted(bad_buckets))}\n")
        _log_write(logfile, msg)
        print(msg, file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
