"""Library vocabulary matching and per-worker caches for hit finding."""

import logging
import os
import re
from collections import Counter, OrderedDict
from contextlib import contextmanager
from functools import lru_cache
from itertools import permutations, product

import numpy as np
import prody as pr
from rdkit import rdBase, Chem

from ligand_vdgs.functions import ligand_structure, vdg_npz_utils as vdg_npz
from ligand_vdgs.functions.utils import (_find_resonance_terminal_groups, filename_to_smiles,
    fragment_keys_equivalent, init_query_ring_info, smiles_to_filename)

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

def lib_entries(vdg_lib_dir):
    """Fragment directory names in `vdg_lib_dir` (skips root files such as provenance json)."""
    return frozenset(e.name for e in os.scandir(vdg_lib_dir) if e.is_dir())

def init_worker(vdg_lib_entries):
    global _WORKER_MOL_CACHE, _WORKER_MOL_CACHE_ENTRIES, _WORKER_BUCKET_CACHE
    global _WORKER_PATTERN_ELEMENTS
    _WORKER_BUCKET_CACHE = None
    entries = frozenset(vdg_lib_entries)
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

def query_resonance_forms(lig_mol):
    """lig_mol plus copies with each terminal resonance group's end-atom states (bond order to the
    center, formal charge, explicit Hs) permuted over its equivalent atoms. Bond orders come from one
    arbitrary Kekule form (P=O1P vs the library parent's P=O2P), so bond-order keys must be matched
    against every form or they miss the equivalent oxygens."""
    gid_of = {b: gid for b, (_z, gid) in _find_resonance_terminal_groups(lig_mol).items()}
    slots = []
    for bond_ids in ([b for b in gid_of if gid_of[b] == g] for g in sorted(set(gid_of.values()))):
        bonds = [lig_mol.GetBondWithIdx(i) for i in bond_ids]
        ends = [b.GetEndAtom() if b.GetEndAtom().GetDegree() == 1 else b.GetBeginAtom() for b in bonds]
        states = [(b.GetBondType(), a.GetFormalCharge(), a.GetNumExplicitHs(), a.GetNoImplicit())
                  for b, a in zip(bonds, ends)]
        perms = sorted(set(permutations(states)) - {tuple(states)})
        if perms:
            slots.append([(bond_ids, [a.GetIdx() for a in ends], p) for p in [tuple(states)] + perms])
    forms = [lig_mol]
    for choice in list(product(*slots))[1:]:
        rw = Chem.RWMol(lig_mol)
        for bond_ids, atom_ids, states in choice:
            for b, a, (btype, charge, n_h, no_imp) in zip(bond_ids, atom_ids, states):
                rw.GetBondWithIdx(b).SetBondType(btype)
                atom = rw.GetAtomWithIdx(a)
                atom.SetFormalCharge(charge), atom.SetNumExplicitHs(n_h), atom.SetNoImplicit(no_imp)
        form = rw.GetMol()
        form.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(form)
        forms.append(form)
    return forms

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
        twin_usable = ligand_structure.check_vdg_job_status(twin, vdg_lib_dir)
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
            'so matching must run on an H-free mol (see ligand_structure.get_query_ligand_mol).')
    lig_elements = Counter(a.GetAtomicNum() for a in lig_mol.GetAtoms())
    n_lig_atoms = lig_mol.GetNumAtoms()
    forms = query_resonance_forms(lig_mol)
    for db_name in sorted(_WORKER_MOL_CACHE):
        pattern = _WORKER_MOL_CACHE[db_name]
        if pattern.GetNumAtoms() > n_lig_atoms:
            continue
        if any(count > lig_elements.get(atomic_num, 0)
               for atomic_num, count in _pattern_element_counts(db_name).items()):
            continue
        if not any(form.HasSubstructMatch(pattern) for form in forms):
            _warn_charge_only_miss(lig_mol, db_name, vdg_lib_dir,
                                   frags_in_lib, warn)
            continue

        if db_name not in frags_in_lib:
            frags_in_lib[db_name] = ligand_structure.check_vdg_job_status(db_name, vdg_lib_dir)
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
            per_form = [form.GetSubstructMatches(_cg_order_pattern(vdg_lib_dir, db_name), uniquify=False,
                                                 maxMatches=_MAX_QUERY_FRAG_MATCHES) for form in forms]
        except ValueError as e:
            if warn is not None:
                warn(f"{e}; skipping this fragment.")
            continue
        # First form wins per match so each site's submol carries bond orders its key matched.
        matches = dict(reversed([(m, form) for form, ms in zip(forms, per_form) for m in ms]))
        if not matches:
            continue
        if max(map(len, per_form)) >= _MAX_QUERY_FRAG_MATCHES:
            if warn is not None:
                warn(f"Fragment {db_name!r} hit the {_MAX_QUERY_FRAG_MATCHES}-match "
                     f"cap on this ligand, so its labelings would be incomplete; "
                     f"skipping it.")
            continue

        # One site per distinct atom set (its labelings are CG automorphism orbits); overlapping
        # but different atom sets stay separate sites, so no transitive merge shares a q_site_idx.
        filtered_frags[db_name] = [
            _expand_site_by_cg_automorphisms(
                [(ligand_structure.submol_from_match(form, match), tuple(match), key)
                 for match, form in matches.items() if tuple(sorted(match)) == key], automorphisms)
            for key in dict.fromkeys(tuple(sorted(match)) for match in matches)]
        query_frag_map[db_name] = target_smarts

    return filtered_frags, query_frag_map, frags_in_lib
