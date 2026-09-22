"""
Core vdG hit-finding: score a query structure's ligand against the vdG fragment
library. The RDKit mol and bucket caches
are per process (eager via init_worker in pool workers, else lazy), so scoring
many models in one process reuses the library arrays.
"""

import io
import logging
import multiprocessing as mp
import os
import re
import tempfile
import traceback
from collections import Counter, OrderedDict
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from functools import lru_cache

import numpy as np
import prody as pr
from rdkit import rdBase, Chem

from ligand_vdgs.functions import Frags
from ligand_vdgs.functions import dock_utils as dock
from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions import vdg_struct_utils as struct_utils
from ligand_vdgs.functions.utils import (kabsch, kabsch_ssd, best_inplace_symmetry_rmsd,
    filename_to_smiles, smiles_to_filename, normalize_rmsd, init_query_ring_info,
    fragment_keys_equivalent)
from ligand_vdgs.functions.vdg_fp_utils import (fp_tolerances,
    prefilter_query_indices_pair, prefilter_query_indices_single)

EXCELLENT_MATCH_CUTOFF = 0.3
BUCKET_CACHE_MAX_BYTES = 128 * 1024**2

def _struct_id_from_pdbfile(pdbfile):
    basename = os.path.basename(os.fspath(pdbfile))
    for extension in (".pdb.gz", ".cif.gz", ".pdb", ".cif", ".gz"):
        if basename.lower().endswith(extension):
            basename = basename[:-len(extension)]
            break

    structure, separator, pose = basename.rpartition("_")
    if separator and structure and pose.isdecimal():
        return structure
    return basename

_WORKER_MOL_CACHE: dict | None = None
_WORKER_MOL_CACHE_ENTRIES: frozenset | None = None
_WORKER_PATTERN_ELEMENTS: dict | None = None
_warned_incomplete_frags = set()

def _bucket_nbytes(bucket):
    if bucket is None:
        return 0
    return sum(
        value.nbytes
        for value in bucket.values()
        if isinstance(value, np.ndarray))

CACHE_MISS = object()

class _BoundedBucketCache:
    def __init__(self, max_bytes=BUCKET_CACHE_MAX_BYTES):
        self.max_bytes = max_bytes
        self._entries = OrderedDict()
        self._nbytes = 0

    def get(self, key):
        try:
            value, nbytes = self._entries.pop(key)
        except KeyError:
            return CACHE_MISS
        self._entries[key] = (value, nbytes)
        return value

    def put(self, key, value):
        nbytes = _bucket_nbytes(value)
        if nbytes > self.max_bytes:
            return

        if value is not None:
            for array in value.values():
                if isinstance(array, np.ndarray):
                    array.flags.writeable = False

        old = self._entries.pop(key, None)
        if old is not None:
            self._nbytes -= old[1]

        self._entries[key] = (value, nbytes)
        self._nbytes += nbytes
        while self._nbytes > self.max_bytes:
            _, (_, evicted_nbytes) = self._entries.popitem(last=False)
            self._nbytes -= evicted_nbytes

_WORKER_BUCKET_CACHE: _BoundedBucketCache | None = None

def _worker_bucket_cache():
    global _WORKER_BUCKET_CACHE
    if _WORKER_BUCKET_CACHE is None:
        _WORKER_BUCKET_CACHE = _BoundedBucketCache()
    return _WORKER_BUCKET_CACHE

_WORKER_STRUCT_CACHE: tuple | None = None

def _combo_worker_struct(pdb_path):
    global _WORKER_STRUCT_CACHE
    if _WORKER_STRUCT_CACHE is None or _WORKER_STRUCT_CACHE[0] != pdb_path:
        _WORKER_STRUCT_CACHE = (pdb_path, pr.parseCIF(pdb_path)
                                 if pdb_path.endswith((".cif", ".cif.gz"))
                                 else pr.parsePDB(pdb_path))
    return _WORKER_STRUCT_CACHE[1]

def init_worker(lib_entries):
    global _WORKER_MOL_CACHE, _WORKER_MOL_CACHE_ENTRIES, _WORKER_BUCKET_CACHE
    global _WORKER_PATTERN_ELEMENTS
    _WORKER_BUCKET_CACHE = None
    entries = frozenset(lib_entries)
    mol_cache = {}
    for entry in entries:
        mol = init_query_ring_info(Chem.MolFromSmarts(filename_to_smiles(entry)))
        if mol is not None:
            mol_cache[entry] = mol
    _WORKER_MOL_CACHE = mol_cache
    _WORKER_PATTERN_ELEMENTS = {
        name: Counter(a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum())
        for name, mol in mol_cache.items()}
    _WORKER_MOL_CACHE_ENTRIES = entries

def _pattern_element_counts(db_name):
    return _WORKER_PATTERN_ELEMENTS[db_name]

def _ensure_worker_mol_cache(lib_entries):
    entries = frozenset(lib_entries)
    if (_WORKER_MOL_CACHE is None
            or _WORKER_MOL_CACHE_ENTRIES != entries):
        init_worker(entries)

@contextmanager
def init_rdkit_logging(stream_like):
    rdBase.LogToPythonLogger()
    logger = logging.getLogger("rdkit")
    previous = logger.level, logger.propagate
    logger.setLevel(logging.INFO)
    logger.propagate = False
    h = logging.StreamHandler(stream_like)
    h.setFormatter(logging.Formatter("%(levelname)s [%(name)s]: %(message)s"))
    logger.addHandler(h)
    try:
        yield
    finally:
        h.flush()
        logger.removeHandler(h)
        h.close()
        logger.setLevel(previous[0])
        logger.propagate = previous[1]

def bsr_label_to_string(bsr_combo):
    return ";".join(f"{seg}:{chain}:{resnum}" for (seg, chain, resnum) in bsr_combo)

@lru_cache(maxsize=512)
def _cg_symmetry_for(vdg_lib_dir, db_name):
    return vdg_npz.load_cg_symmetry(vdg_lib_dir, db_name)

@lru_cache(maxsize=512)
def _cg_order_pattern(vdg_lib_dir, db_name):
    smarts = _cg_symmetry_for(vdg_lib_dir, db_name)[0]
    pattern = init_query_ring_info(Chem.MolFromSmarts(smarts))
    if pattern is None:
        raise ValueError(f"Fragment {db_name!r}: recorded cg_smarts {smarts!r} "
                         "does not parse as SMARTS")
    return pattern

def cg_atom_order_smarts(vdg_lib_dir, db_name):
    return _cg_symmetry_for(vdg_lib_dir, db_name)[0]

def _expand_site_by_cg_automorphisms(site, automorphisms):
    if len(automorphisms) <= 1:
        return site
    n_auto = len(automorphisms[0])
    expanded, seen = [], set()
    for sub, perm_inds, orig_mol_inds in site:
        if sub.GetNumAtoms() != n_auto:
            raise ValueError(
                f"CG automorphisms are over {n_auto} atoms but the query "
                f"fragment has {sub.GetNumAtoms()}; the library's recorded "
                "symmetry does not describe this fragment")
        for auto in automorphisms:
            new_perm = tuple(perm_inds[i] for i in auto)
            key = (new_perm, tuple(orig_mol_inds))
            if key in seen:
                continue
            try:
                expanded.append((Chem.RenumberAtoms(Chem.Mol(sub), list(auto)),
                                 new_perm, orig_mol_inds))
            except Exception:
                continue
            seen.add(key)
    return expanded or site

_warned_unsearchable_libs = set()

def _warn_unsearchable_entries(vdg_lib_dir, vdg_lib_entries, warn):
    if warn is None or vdg_lib_dir in _warned_unsearchable_libs:
        return
    _warned_unsearchable_libs.add(vdg_lib_dir)
    unparsed = [e for e in sorted(set(vdg_lib_entries) - set(_WORKER_MOL_CACHE))
                if os.path.isdir(os.path.join(vdg_lib_dir, e))]
    if unparsed:
        warn(f"{len(unparsed)} fragment directory name(s) in {vdg_lib_dir} do not "
             f"parse as SMARTS, so their vdGs cannot be searched for: {unparsed}")

_warned_charge_only_frags = set()
_SMARTS_CHARGE = re.compile(r"([A-Za-z])([+-]\d*)([;\]])")

@lru_cache(maxsize=None)
def _neutralized_pattern(db_name):
    neutral, n_subs = _SMARTS_CHARGE.subn(r"\1\3", filename_to_smiles(db_name))
    if not n_subs:
        return None
    pattern = init_query_ring_info(Chem.MolFromSmarts(neutral))
    return None if pattern is None else (pattern, neutral)

def _library_entry_for_key(key):
    exact = smiles_to_filename(key)
    if exact in _WORKER_MOL_CACHE:
        return exact
    for entry in sorted(_WORKER_MOL_CACHE):
        try:
            if fragment_keys_equivalent(key, filename_to_smiles(entry)):
                return entry
        except Exception:
            continue
    return None

def _warn_charge_only_miss(lig_mol, db_name, vdg_lib_dir, frags_in_lib, warn):
    if warn is None or db_name in _warned_charge_only_frags:
        return
    neutralized = _neutralized_pattern(db_name)
    if neutralized is None:
        return
    pattern, neutral_smarts = neutralized
    if not lig_mol.HasSubstructMatch(pattern):
        return
    _warned_charge_only_frags.add(db_name)
    twin = _library_entry_for_key(neutral_smarts)
    if twin is None:
        twin_usable = False
    elif twin in frags_in_lib:
        twin_usable = frags_in_lib[twin]
    else:
        twin_usable = Frags.check_vdg_job_status(twin, vdg_lib_dir)
    if twin_usable:
        warn(f"Fragment {db_name!r} is charged and did not match this ligand, which "
             f"carries that moiety as {neutral_smarts!r}. Its neutral twin is in the "
             f"library and is searched instead, so the moiety is covered, but this "
             f"key's own vdGs are not reachable for this query.")
    else:
        warn(f"Fragment {db_name!r} is charged and did not match this ligand, but its "
             f"neutral form {neutral_smarts!r} does: the ligand carries that moiety in "
             f"the other protonation state. No usable neutral twin is in the library "
             f"to fall back on, so this moiety is unsearchable for this query and its "
             f"vdGs are silently absent from the results.")

_MAX_QUERY_FRAG_MATCHES = 1000

def match_library_frags_to_query(lig_mol, vdg_lib_dir, vdg_lib_entries,
                                 frags_in_lib=None, warn=None):
    if frags_in_lib is None:
        frags_in_lib = {}
    filtered_frags = {}
    query_frag_map = {}

    _ensure_worker_mol_cache(vdg_lib_entries)
    _warn_unsearchable_entries(vdg_lib_dir, vdg_lib_entries, warn)
    Chem.GetSymmSSSR(lig_mol)
    n_graph_h = sum(1 for a in lig_mol.GetAtoms() if a.GetAtomicNum() == 1)
    if n_graph_h:
        raise ValueError(
            f'query ligand mol has {n_graph_h} hydrogen atoms in the graph; the '
            'library\'s keys carry heavy-atom degree (`D<n>`), which counts them, '
            'so matching must run on an H-free mol (see Frags.get_query_ligand_mol).')
    lig_elements = Counter(a.GetAtomicNum() for a in lig_mol.GetAtoms())
    n_lig_atoms = lig_mol.GetNumAtoms()
    for db_name in sorted(_WORKER_MOL_CACHE):
        pattern = _WORKER_MOL_CACHE[db_name]
        if pattern.GetNumAtoms() > n_lig_atoms:
            continue
        if any(count > lig_elements.get(atomic_num, 0)
               for atomic_num, count in _pattern_element_counts(db_name).items()):
            continue
        if not lig_mol.HasSubstructMatch(pattern):
            _warn_charge_only_miss(lig_mol, db_name, vdg_lib_dir,
                                   frags_in_lib, warn)
            continue

        if db_name not in frags_in_lib:
            frags_in_lib[db_name] = Frags.check_vdg_job_status(db_name, vdg_lib_dir)
            if (not frags_in_lib[db_name] and warn is not None
                    and db_name not in _warned_incomplete_frags):
                _warned_incomplete_frags.add(db_name)
                warn(f"Fragment {db_name!r} directory exists in {vdg_lib_dir} but "
                     f"its vdG-generation job has not completed (no 'Job completed.' "
                     f"in its log); treating it as absent from the library.")
        if not frags_in_lib[db_name]:
            continue

        target_smarts, automorphisms = _cg_symmetry_for(vdg_lib_dir, db_name)
        try:
            matches = lig_mol.GetSubstructMatches(
                _cg_order_pattern(vdg_lib_dir, db_name), uniquify=False,
                maxMatches=_MAX_QUERY_FRAG_MATCHES)
        except ValueError as e:
            if warn is not None:
                warn(f"{e}; skipping this fragment.")
            continue
        if not matches:
            continue
        if len(matches) >= _MAX_QUERY_FRAG_MATCHES:
            if warn is not None:
                warn(f"Fragment {db_name!r} hit the {_MAX_QUERY_FRAG_MATCHES}-match "
                     f"cap on this ligand, so its labelings would be incomplete; "
                     f"skipping it.")
            continue

        filtered_frags[db_name] = [
            _expand_site_by_cg_automorphisms(site, automorphisms)
            for site in Frags.group_lig_sites_by_overlap(
                [(Frags.submol_from_match(lig_mol, match), tuple(match),
                  tuple(sorted(match))) for match in matches])]
        query_frag_map[db_name] = target_smarts

    return filtered_frags, query_frag_map, frags_in_lib

def cg_center(cg_coords):
    return np.asarray(cg_coords, dtype=np.float32).mean(axis=0, dtype=np.float32)

def _dist(a, b):
    d = a - b
    return float(np.sqrt(np.sum(d * d, dtype=np.float32), dtype=np.float32))

def fp_single_from_ca_and_cgcom(ca_coords, cg_com):
    return _dist(ca_coords[0], cg_com)

def fp_pair_from_ca_and_cgcom(ca_coords, cg_com):
    return (_dist(ca_coords[0], cg_com), _dist(ca_coords[1], cg_com),
            _dist(ca_coords[0], ca_coords[1]))

def get_residue_all_atom_coords(struct, seg, chain, resnum):
    terms = []
    if seg not in (None, "", "_"):
        terms.append(f"segment {seg}")
    if chain not in (None, "", "_"):
        terms.append(f"chain {chain}")
    terms.append(f"resnum {resnum}")
    sel = struct.select(" and ".join(terms))
    if sel is None:
        return None

    coords = sel.getCoords()
    if coords is None or len(coords) == 0:
        return None

    try:
        elems, names = sel.getElements(), sel.getNames()
    except Exception:
        elems = names = None

    if elems is not None and names is not None:
        mask = np.array([not struct_utils.is_hydrogen(n, e)
                         for n, e in zip(names, elems)], dtype=bool)
        if not np.any(mask):
            return None
        coords = coords[mask]

    return np.asarray(coords, dtype=np.float32)

def get_bsr_all_atom_coords(struct, bsr_combo):
    chunks = []
    for seg, chain, resnum in bsr_combo:
        c = get_residue_all_atom_coords(struct, seg, chain, resnum)
        if c is not None and len(c) > 0:
            chunks.append(c)
    if not chunks:
        return None
    return np.concatenate(chunks, axis=0).astype(np.float32, copy=False)

BB_SLOT_SIDECHAIN_CLASH = 3.4
PRO_NH_DONOR_CUTOFF = 3.5

def backbone_slot_blockers(struct, bsr_combo, slot_labels, slot_resnames):
    blockers = []
    for (seg, chain, resnum), label, resname in zip(bsr_combo, slot_labels, slot_resnames):
        if label not in struct_utils.BB_LABELS:
            continue
        terms = []
        if seg not in (None, "", "_"):
            terms.append(f"segment {seg}")
        if chain not in (None, "", "_"):
            terms.append(f"chain {chain}")
        terms.append(f"resnum {resnum}")
        sel = struct.select(" and ".join(terms))
        if sel is None:
            continue
        _bb, sc, extra = struct_utils.split_residue_heavy_atoms(sel, str(resname))
        if sc is None:
            continue
        if extra is not None and len(extra):
            sc = np.concatenate((sc, extra), axis=0) if len(sc) else extra
        n_coord = None
        if str(resname) == "PRO":
            n_sel = sel.select("name N")
            if n_sel is not None and len(n_sel):
                n_coord = np.asarray(n_sel.getCoords()[0], dtype=np.float32)
        if len(sc) or n_coord is not None:
            blockers.append((np.asarray(sc, dtype=np.float32), n_coord))
    return blockers or None

def backbone_slots_can_host(blockers, cg_coords, cg_is_acceptor):
    if not blockers:
        return True
    cg_coords = np.asarray(cg_coords, dtype=np.float32)
    for sc, n_coord in blockers:
        if len(sc):
            d2 = ((sc[:, None, :] - cg_coords[None, :, :]) ** 2).sum(-1)
            if d2.min() < BB_SLOT_SIDECHAIN_CLASH ** 2:
                return False
        if n_coord is not None and cg_is_acceptor is not None and cg_is_acceptor.any():
            d2 = ((cg_coords[cg_is_acceptor] - n_coord) ** 2).sum(-1)
            if d2.min() < PRO_NH_DONOR_CUTOFF ** 2:
                return False
    return True

def has_any_bsr_atom_cg_contact(bsr_atom_coords, cg_coords, cutoff=3.8):
    if bsr_atom_coords is None or len(bsr_atom_coords) == 0:
        return False
    cg_coords = np.asarray(cg_coords, dtype=np.float32)
    if cg_coords.size == 0:
        return False
    diff = bsr_atom_coords[:, None, :] - cg_coords[None, :, :]
    return bool(np.any(np.sum(diff * diff, axis=2, dtype=np.float32)
                      <= np.float32(cutoff * cutoff)))

def deduplicate_hits(match_records, min_shared_atoms=3):
    def _parse_atom_indices(rec):
        s = rec.get("q_atom_indices", "")
        if not s:
            return None
        try:
            return frozenset(int(x) for x in s.split(";") if x)
        except ValueError:
            return None

    def _bsr_residue_set(rec):
        return frozenset(t for t in str(rec["bsr_combo"]).split(";") if t)

    kept = []
    kept_atoms = []
    kept_bsr = []

    for rec in sorted(match_records, key=lambda r: float(r["vdg_rmsd"])):
        bsr = _bsr_residue_set(rec)
        q_atoms = _parse_atom_indices(rec)

        is_dup = False
        if q_atoms is not None:
            for k_bsr, k_atoms in zip(kept_bsr, kept_atoms):
                if k_bsr != bsr:
                    continue
                if k_atoms is None:
                    continue
                if len(q_atoms & k_atoms) >= min_shared_atoms:
                    is_dup = True
                    break

        if not is_dup:
            kept.append(rec)
            kept_atoms.append(q_atoms)
            kept_bsr.append(bsr)

    return kept

_BUCKET_ROW_FIELDS = ("cluster_id", "cluster_num_parents", "cg", "bb", "resnames", "slot_flags")

def _load_vdg_bucket_all_signs(vdg_lib_dir, frag_name, subset_size, aa_bucket):
    loaded = [(sign, vdg_npz.load_vdg_bucket(vdg_lib_dir, frag_name, subset_size, sign,
                                              aa_bucket))
              for sign in vdg_npz.CHARGE_SIGNS]
    loaded = [(sign, bucket) for sign, bucket in loaded if bucket is not None]
    if not loaded:
        return None
    return dict(
        aa_bucket_parts=loaded[0][1]["aa_bucket_parts"],
        charge_signs=np.concatenate([np.full(len(bucket["cluster_id"]), sign, dtype="U11")
                                     for sign, bucket in loaded]),
        partition_indices=np.concatenate([np.arange(len(bucket["cluster_id"]), dtype=np.int32)
                                          for _sign, bucket in loaded]),
        **{field: np.concatenate([bucket[field] for _sign, bucket in loaded])
           for field in _BUCKET_ROW_FIELDS})

def _score_one_bsr_combo(pdbfile, struct, frag_name, query_frag, grouped_q_cg_perms,
                          combo_item, vdg_lib_dir, bucket_cache, bsr_atom_coords_cache,
                          bb_blockers_cache, rmsd_threshold, contact_cutoff,
                          lig_instance_label):
    combo, bsr_combo, _bsr_AAs, coords = combo_item
    match_records = []
    paired = sorted(zip(combo, coords, bsr_combo, _bsr_AAs), key=lambda x: x[0])
    (bsr_incl_bb_identities, input_bsr_bb_coords, bsr_combo,
     bsr_resnames) = zip(*paired)
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

    combo_key = tuple(bsr_combo)
    if combo_key in bsr_atom_coords_cache:
        bsr_atom_coords = bsr_atom_coords_cache[combo_key]
    else:
        bsr_atom_coords = get_bsr_all_atom_coords(struct, bsr_combo)
        bsr_atom_coords_cache[combo_key] = bsr_atom_coords
        if bsr_atom_coords is None and contact_cutoff is not None:
            print(f"[WARNING] ({pdbfile}) No non-H atoms for BSR combo "
                  f"{bsr_label_to_string(bsr_combo)}; the contact filter "
                  f"cannot pass, so this combo is skipped.", flush=True)
    if bsr_atom_coords is None and contact_cutoff is not None:
        return match_records

    blocker_key = (combo_key, tuple(bsr_incl_bb_identities))
    if blocker_key in bb_blockers_cache:
        bb_blockers = bb_blockers_cache[blocker_key]
    else:
        bb_blockers = backbone_slot_blockers(
            struct, bsr_combo, bsr_incl_bb_identities, bsr_resnames)
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

            query_bb_slot_ok[idx] = backbone_slots_can_host(
                bb_blockers, q_cg_coords, q_acceptor)

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

    if not np.any(query_bb_slot_ok):
        return match_records

    bucket_key = (vdg_lib_dir, frag_name, subset_size, aa_bucket)
    bucket = bucket_cache.get(bucket_key)
    if bucket is CACHE_MISS:
        bucket = _load_vdg_bucket_all_signs(vdg_lib_dir, frag_name, subset_size, aa_bucket)
        bucket_cache.put(bucket_key, bucket)
    if bucket is None:
        return match_records

    cluster_num_parents = bucket["cluster_num_parents"]
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
    bb_rmsd_all = np.sqrt(kabsch_ssd(
        bb_flat,
        bucket_bb[:, np.asarray(resind_perms, dtype=np.intp), :, :].reshape(
            n_nr_vdgs * n_perms, -1, 3),
    ) / n_atoms).reshape(n_nr_vdgs, n_perms)

    for nr_idx in range(n_nr_vdgs):
        vdg_cg      = bucket_cg[nr_idx]
        vdg_bb      = bucket_bb[nr_idx]
        vdg_idx     = int(nr_idx)
        clus_id     = int(bucket["cluster_id"][nr_idx])
        clus_num_parents = int(cluster_num_parents[nr_idx])

        if vdg_cg.shape[0] != N_cg:
            continue

        vdg_cg_com = cg_center(vdg_cg)

        best_rmsd = best_aa_perm_idx = best_q_site_idx = None
        best_q_cg_perm_idx = best_R = best_t = None
        best_R_bb = best_t_bb = None

        for aa_perm_idx, resind_perm in enumerate(resind_perms):
            resind_perm_arr = np.asarray(resind_perm, dtype=np.intp)

            res_subset_bb = vdg_bb[resind_perm_arr]
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

            vdg_perm_bb = res_subset_bb.reshape(-1, 3)
            if bb_rmsd_all[nr_idx, aa_perm_idx] > effective_rmsd:
                continue

            Y_sub = Y[idxs]

            db_bb_and_cg = np.concatenate(
                (vdg_perm_bb, vdg_cg), axis=0).astype(np.float32,
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
                    best_q_cg_perm_idx = int(meta[idxs[local_idx], 1])
                    best_R             = R_cand
                    best_t             = t_one[0]
                    R_bb_one, t_bb_one, _ = kabsch(vdg_perm_bb, bb_flat[None])
                    best_R_bb          = R_bb_one[0]
                    best_t_bb          = t_bb_one[0]

                    if best_rmsd <= EXCELLENT_MATCH_CUTOFF:
                        break

        if best_rmsd is not None:
            struct_id = _struct_id_from_pdbfile(pdbfile)
            R_rounded = np.round(best_R, 4)
            t_rounded = np.round(best_t, 4)

            best_orig_mol_inds = grouped_q_cg_perms[best_q_site_idx][best_q_cg_perm_idx][-1]
            q_atom_indices_str = ";".join(
                str(i) for i in sorted(best_orig_mol_inds))

            rec = dict(
                pdbfile=pdbfile,
                struct_id=struct_id,
                lig_instance=lig_instance_label,
                frag=frag_name,
                query_frag=query_frag,
                subset_size=subset_size,
                bsr_combo=bsr_label_to_string(bsr_combo),
                aa_bucket=aa_bucket,
                charge_sign=str(bucket["charge_signs"][vdg_idx]),
                vdg_index=int(bucket["partition_indices"][vdg_idx]),
                vdg_cluster_id=clus_id,
                vdg_cluster_num_parents=clus_num_parents,
                vdg_rmsd=f"{best_rmsd:.4f}",
                rmsd_threshold=f"{effective_rmsd:.4f}",
                aa_perm_idx=int(best_aa_perm_idx),
                q_site_idx=int(best_q_site_idx),
                q_cg_perm_idx=int(best_q_cg_perm_idx),
                q_atom_indices=q_atom_indices_str,)
            rec.update({f"R{i}{j}": f"{R_rounded[i, j]:.4f}"
                        for i in range(3) for j in range(3)})
            rec.update({f"t{k}": f"{t_rounded[k]:.4f}" for k in range(3)})
            Rbb_rounded = np.round(best_R_bb, 4)
            tbb_rounded = np.round(best_t_bb, 4)
            rec.update({f"Rbb{i}{j}": f"{Rbb_rounded[i, j]:.4f}"
                        for i in range(3) for j in range(3)})
            rec.update({f"tbb{k}": f"{tbb_rounded[k]:.4f}" for k in range(3)})

            match_records.append(rec)

    return match_records

def _combo_worker(task):
    (pdbfile, pdb_path, frag_name, query_frag, grouped_q_cg_perms, combo_item,
     vdg_lib_dir, rmsd_threshold, contact_cutoff, lig_instance_label) = task
    struct = _combo_worker_struct(pdb_path)
    bucket_cache = _worker_bucket_cache()
    local_buf = io.StringIO()
    with redirect_stdout(local_buf), redirect_stderr(local_buf), init_rdkit_logging(local_buf):
        try:
            match_records = _score_one_bsr_combo(
                pdbfile, struct, frag_name, query_frag, grouped_q_cg_perms,
                combo_item, vdg_lib_dir, bucket_cache, {}, {},
                rmsd_threshold, contact_cutoff, lig_instance_label)
            return local_buf.getvalue(), match_records, None
        except Exception:
            return local_buf.getvalue(), [], traceback.format_exc()

def _raise_if_combo_task_errors(errors, n_tasks):
    if errors:
        raise RuntimeError(
            f"{len(errors)}/{n_tasks} combo tasks failed; first traceback:\n{errors[0]}")

def _match_record_sort_key(rec):
    return (rec["frag"], rec["bsr_combo"], rec["aa_bucket"], rec["charge_sign"],
            rec["vdg_index"],
            rec["aa_perm_idx"], rec["q_site_idx"], rec["q_cg_perm_idx"])

def score_one_model(
    pdbfile, pdb_path, lig_smiles, vdg_lib_dir, rmsd_threshold=None,
    ref_lig_mol=None, print_bsr_selection=False, contact_cutoff=None,
    vdg_lib_entries=None, deduplicate=True, min_shared_atoms=3, lig_instance_label='',
    nprocs=1,):
    if vdg_lib_entries is None:
        vdg_lib_entries = set(os.listdir(vdg_lib_dir))

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
                lig_mol_noH = Frags.get_query_ligand_mol(struct, lig_smiles)
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
                _, ligname = Frags.identify_ligand_selection(struct, lig_smiles)
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

            bsr_atom_coords_cache = {}
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
                            q_cg_com = cg_center(q_cg_coords)
                            q_acceptor = np.array(
                                [a.GetSymbol() in ("N", "O", "S")
                                 for a in sub.GetAtoms()], dtype=bool)
                            perms.append((q_cg_coords, q_cg_com, q_acceptor,
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
                 vdg_lib_dir, rmsd_threshold, contact_cutoff, lig_instance_label)
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
                            combo_item, vdg_lib_dir, bucket_cache, bsr_atom_coords_cache,
                            bb_blockers_cache, rmsd_threshold, contact_cutoff,
                            lig_instance_label))

        if deduplicate and match_records:
            match_records = deduplicate_hits(match_records,
                                             min_shared_atoms=min_shared_atoms)

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
        ligname, instances = Frags.classify_ligand_instances(struct, lig_smiles)
    except ValueError:
        instances = []
    if len(instances) <= 1:
        return score_one_model(pdbfile, pdb_path, lig_smiles, vdg_lib_dir, **kwargs)

    log_texts, match_records, frags_in_lib, filtered_frags = [], [], {}, {}
    for resindex, label in instances:
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tf:
            instance_path = tf.name
        try:
            pr.writePDB(instance_path, Frags.select_struct_for_ligand_instance(
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
     print_bsr_selection, contact_cutoff, deduplicate, min_shared_atoms) = args
    vdg_lib_entries = _WORKER_MOL_CACHE_ENTRIES

    if vdg_lib_entries is None:
        vdg_lib_entries = frozenset(os.listdir(vdg_lib_dir))
        init_worker(vdg_lib_entries)
    try:
        log_text, match_records, frags_in_lib, lig_rmsd_value, _ = score_one_model_multi_instance(
            pdbfile=pdbfile, pdb_path=pdb_path, lig_smiles=lig_smiles,
            vdg_lib_dir=vdg_lib_dir, rmsd_threshold=rmsd_threshold,
            ref_lig_mol=ref_lig_mol, print_bsr_selection=print_bsr_selection,
            contact_cutoff=contact_cutoff, vdg_lib_entries=vdg_lib_entries,
            deduplicate=deduplicate, min_shared_atoms=min_shared_atoms)
        return pdbfile, log_text, lig_rmsd_value, match_records, frags_in_lib, None
    except Exception:
        error_text = (f"Worker failed for pdbfile={pdbfile}, pdb_path={pdb_path}\n"
                      f"{traceback.format_exc()}")
        return pdbfile, "", None, [], {}, error_text
