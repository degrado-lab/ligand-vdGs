"""
Core vdG hit-finding: score a query structure's ligand against the vdG fragment
library. Public API: score_one_model, match_library_frags_to_query,
deduplicate_hits, init_worker, process_work_item. The RDKit mol and bucket caches
are per process (eager via init_worker in pool workers, else lazy), so scoring
many models in one process reuses the library arrays.
"""

import io
import logging
import os
import re
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


# Constants                                                           #

EXCELLENT_MATCH_CUTOFF = 0.3   # Å; early-exit from AA permutation loop
# Per *process*: a pool of N workers holds up to N times this resident.
BUCKET_CACHE_MAX_BYTES = 128 * 1024**2


# Query filename helpers                                            #

def _struct_id_from_pdbfile(pdbfile):
    """Structure-group ID from a ``<structure>_<pose>.<ext>`` query filename:
    strip the extension and a trailing numeric pose suffix (structure IDs are not
    assumed to be four characters)."""
    basename = os.path.basename(os.fspath(pdbfile))
    lowercase_basename = basename.lower()
    for extension in (".pdb.gz", ".cif.gz", ".pdb", ".cif", ".gz"):
        if lowercase_basename.endswith(extension):
            basename = basename[:-len(extension)]
            break

    structure, separator, pose = basename.rpartition("_")
    if separator and structure and pose.isdecimal():
        return structure
    return basename


# Per-worker RDKit mol cache                                          #

_WORKER_MOL_CACHE: dict | None = None
_WORKER_MOL_CACHE_ENTRIES: frozenset | None = None
# {db_name: {atomic number: count}} for the atoms each key pins down an element
# for; the prefilter in match_library_frags_to_query needs it per model, and it
# depends only on the library.
_WORKER_PATTERN_ELEMENTS: dict | None = None

# Module-level (per-process) so a repeatedly-matched incomplete fragment warns
# once per worker process rather than once per model scored.
_warned_incomplete_frags = set()


def _bucket_nbytes(bucket):
    """Return the NumPy storage used by a loaded vdG bucket."""
    if bucket is None:
        return 0
    return sum(
        value.nbytes
        for value in bucket.values()
        if isinstance(value, np.ndarray))


# Returned by _BoundedBucketCache.get for a never-stored key. Distinct from None,
# which load_vdg_bucket returns for a nonexistent bucket -- the common outcome;
# conflating them would re-stat every missing bucket on every model.
CACHE_MISS = object()


class _BoundedBucketCache:
    """Small LRU cache for loaded vdG buckets, bounded by array storage."""

    def __init__(self, max_bytes=BUCKET_CACHE_MAX_BYTES):
        self.max_bytes = max_bytes
        self._entries = OrderedDict()
        self._nbytes = 0

    def get(self, key):
        """The cached bucket (possibly None, meaning "no such bucket"), or
        CACHE_MISS if this key has not been looked up before."""
        try:
            value, nbytes = self._entries.pop(key)
        except KeyError:
            return CACHE_MISS
        self._entries[key] = (value, nbytes)
        return value

    def put(self, key, value):
        nbytes = _bucket_nbytes(value)
        if nbytes > self.max_bytes:
            # Larger than the whole budget: still returned, just not cached.
            return

        # Cached buckets outlive the model that loaded them, so freeze them: an
        # in-place write would corrupt every later model in this worker.
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
    """The process-wide bucket cache, created on first use. Per process, not per
    call: consecutive models visit largely the same buckets. Keys carry the
    library directory, so two libraries in one process cannot collide."""
    global _WORKER_BUCKET_CACHE
    if _WORKER_BUCKET_CACHE is None:
        _WORKER_BUCKET_CACHE = _BoundedBucketCache()
    return _WORKER_BUCKET_CACHE


def init_worker(lib_entries):
    """Build the per-process RDKit mol cache (library dir name -> Mol).
    multiprocessing.Pool initializer. Directory names are
    ``smiles_to_filename()``-encoded and must be decoded before parsing, or keys
    carrying '/' or '\\' (double-bond stereo) are permanently unreachable."""
    global _WORKER_MOL_CACHE, _WORKER_MOL_CACHE_ENTRIES, _WORKER_BUCKET_CACHE
    global _WORKER_PATTERN_ELEMENTS
    # Dropped with the mol cache; stale bucket keys carry vdg_lib_dir, so they
    # could only waste budget, not be misread.
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
    # Separate from cache truthiness: a library may have no parseable names.
    _WORKER_MOL_CACHE_ENTRIES = entries


def _pattern_element_counts(db_name):
    """{atomic number: count} for the elements a library key pins down."""
    return _WORKER_PATTERN_ELEMENTS[db_name]


def _ensure_worker_mol_cache(lib_entries):
    """Initialize or refresh the cache for ``lib_entries`` when necessary."""
    entries = frozenset(lib_entries)
    if (_WORKER_MOL_CACHE is None
            or _WORKER_MOL_CACHE_ENTRIES != entries):
        init_worker(entries)


# Small helpers                                                       #

@contextmanager
def init_rdkit_logging(stream_like):
    """Temporarily route RDKit log messages to ``stream_like``. The ``rdkit``
    logger is process-global, so the handler must be removed before the per-model
    buffer goes out of scope; otherwise a later model writes to it."""
    rdBase.LogToPythonLogger()
    logger = logging.getLogger("rdkit")
    previous_level = logger.level
    previous_propagate = logger.propagate
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
        logger.setLevel(previous_level)
        logger.propagate = previous_propagate


def bsr_label_to_string(bsr_combo):
    return ";".join(f"{seg}:{chain}:{resnum}" for (seg, chain, resnum) in bsr_combo)


@lru_cache(maxsize=512)
def _cg_symmetry_for(vdg_lib_dir, db_name):
    """``(target_smarts, automorphisms)`` from the ``cg_symmetry.npz`` sidecar:
    the group the library clustered under, and the atom order it indexes. The
    directory name is only the ``-c`` label and may differ, so re-deriving from
    it would give a different group; no fallback -- a missing sidecar raises."""
    return vdg_npz.load_cg_symmetry(vdg_lib_dir, db_name)


@lru_cache(maxsize=512)
def _cg_order_pattern(vdg_lib_dir, db_name):
    """The recorded ``cg_smarts`` as a query mol, ring info initialized. Never the
    directory name: only the recorded one is the atom order the stored CG coords
    and automorphisms use."""
    smarts = _cg_symmetry_for(vdg_lib_dir, db_name)[0]
    pattern = init_query_ring_info(Chem.MolFromSmarts(smarts))
    if pattern is None:
        raise ValueError(f"Fragment {db_name!r}: recorded cg_smarts {smarts!r} "
                         "does not parse as SMARTS")
    return pattern


def cg_atom_order_smarts(vdg_lib_dir, db_name):
    """The SMARTS whose atom order a fragment's stored CG coords use. Query mols
    from ``match_library_frags_to_query`` are built in it, and
    ``dock.get_query_cg_coords`` validates their element order against it."""
    return _cg_symmetry_for(vdg_lib_dir, db_name)[0]


def _expand_site_by_cg_automorphisms(site, automorphisms):
    """Add every equivalent CG labeling of each entry in one query site.

    RDKit matching is not resonance-normalized, so the query finds fewer labelings
    than the library clustered under (6 of a phosphate's 24, 1 of a carboxylate's
    2), and any it cannot generate silently never matches. Automorphisms preserve
    the element sequence, so ``get_query_cg_coords``'s guard still holds."""
    if len(automorphisms) <= 1:
        return site
    n_auto = len(automorphisms[0])
    expanded, seen = [], set()
    for sub, perm_inds, orig_mol_inds in site:
        if sub.GetNumAtoms() != n_auto:
            # Wrong-length permutations make every RenumberAtoms raise, which
            # the `or site` fallback would hide.
            raise ValueError(
                f"CG automorphisms are over {n_auto} atoms but the query "
                f"fragment has {sub.GetNumAtoms()}; the library's recorded "
                "symmetry does not describe this fragment")
        for auto in automorphisms:
            try:
                renumbered = Chem.RenumberAtoms(Chem.Mol(sub), list(auto))
            except Exception:
                continue
            # perm_inds is identity/dedup only, never coordinates; composing
            # keeps it describing the atom now in each slot.
            new_perm = tuple(perm_inds[i] for i in auto)
            key = (new_perm, tuple(orig_mol_inds))
            if key in seen:
                continue
            seen.add(key)
            expanded.append((renumbered, new_perm, orig_mol_inds))
    return expanded or site


# Libraries already warned about, so a pool worker reports an unsearchable
# fragment directory once rather than once per model.
_warned_unsearchable_libs = set()


def _warn_unsearchable_entries(vdg_lib_dir, vdg_lib_entries, warn):
    """Report fragment directories whose name RDKit cannot parse as SMARTS: the
    searched key set is exactly what parses, so they are unreachable while looking
    present. Non-directory files at the library root are not fragments."""
    if warn is None or vdg_lib_dir in _warned_unsearchable_libs:
        return
    _warned_unsearchable_libs.add(vdg_lib_dir)
    unparsed = [e for e in sorted(set(vdg_lib_entries) - set(_WORKER_MOL_CACHE))
                if os.path.isdir(os.path.join(vdg_lib_dir, e))]
    if unparsed:
        warn(f"{len(unparsed)} fragment directory name(s) in {vdg_lib_dir} do not "
             f"parse as SMARTS, so their vdGs cannot be searched for: {unparsed}")


# Charged keys already reported against this ligand, so a pool worker says it
# once per fragment rather than once per model.
_warned_charge_only_frags = set()

# Formal charge inside a SMARTS atom primitive: `[N+;!R]`, `[O-;!R]`, `[N+2;r5]`.
_SMARTS_CHARGE = re.compile(r"([A-Za-z])([+-]\d*)([;\]])")


@lru_cache(maxsize=None)
def _neutralized_pattern(db_name):
    """(pattern, SMARTS) for a charged library key with its charges dropped; None
    if uncharged or unparseable once neutralized (diagnostic only, never raises)."""
    neutral, n_subs = _SMARTS_CHARGE.subn(r"\1\3", filename_to_smiles(db_name))
    if not n_subs:
        return None
    pattern = init_query_ring_info(Chem.MolFromSmarts(neutral))
    return None if pattern is None else (pattern, neutral)


def _library_entry_for_key(key):
    """The library directory holding ``key``, or None. Not a string comparison:
    one fragment can be spelled several ways, so an equivalence scan
    (``utils.fragment_keys_equivalent``) avoids calling a present twin missing.
    Unparseable names count as missing -- they are unsearchable too. Runs at most
    once per charged fragment per process."""
    exact = smiles_to_filename(key)
    if exact in _WORKER_MOL_CACHE:
        return exact
    for entry in sorted(_WORKER_MOL_CACHE):
        try:
            if fragment_keys_equivalent(key, filename_to_smiles(entry)):
                return entry
        except Exception:
            # A key this predicate cannot compare is not a usable fallback
            # either; keep scanning rather than failing a query over it.
            continue
    return None


def _warn_charge_only_miss(lig_mol, db_name, vdg_lib_dir, frags_in_lib, warn):
    """Report a charged key that missed only because the ligand is drawn neutral.

    A neutral key is charge-loose (`[O;!R]` matches an anionic O) but not the
    reverse, and the miss is invisible -- a fragment matching nothing looks like
    one the ligand lacks. Fires only when the neutralized key *does* match, and
    reports whether a searchable neutral twin still covers the moiety."""
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
    # Don't prime frags_in_lib: the main loop warns on cache miss only.
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


# Cap on substructure matches of one library key in one query ligand.
# uniquify=False returns every labeling (24 per phosphate site), so reaching this
# means an assumption broke; the fragment is dropped loudly.
_MAX_QUERY_FRAG_MATCHES = 1000


def match_library_frags_to_query(lig_mol, vdg_lib_dir, vdg_lib_entries,
                                 frags_in_lib=None, warn=None):
    """
    Find every occurrence of the library's own fragment keys in a query ligand.

    The inverse of re-fragmenting the query: the enumeration parameters that built
    the library are recorded nowhere the query side can read, so matching the
    library's own keys is what makes it self-describing. Keys are matched against
    the *intact* ligand, the only place their ring primitives (`[C;r6]`) hold.

    Returns (filtered_frags, query_frag_map, frags_in_lib); filtered_frags is
    {library dir name: list of site groups}, whose site positions are exactly the
    q_site_idx recorded on hits.
    """
    if frags_in_lib is None:
        frags_in_lib = {}
    filtered_frags = {}       # db_name -> list of grouped_sites
    query_frag_map = {}       # db_name -> library key SMARTS (decoded dir name)

    _ensure_worker_mol_cache(vdg_lib_entries)
    _warn_unsearchable_entries(vdg_lib_dir, vdg_lib_entries, warn)
    # `r<n>` means *smallest* ring: GetSymmSSSR guarantees it, FastFindRings does
    # not. Idempotent, so free if get_query_ligand_mol already did it.
    Chem.GetSymmSSSR(lig_mol)
    # Keys carry heavy-atom degree, and RDKit counts graph H atoms in `D`: one
    # leftover explicit H makes a hydroxyl O read D2, so every `[O;D1]` key
    # misses and hit finding returns an empty result with no error.
    # get_query_ligand_mol already strips them; this catches a caller that
    # built lig_mol some other way.
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
        # Cheap reject before the matcher (thousands of keys per model). Open
        # elements ('[N,O]', atomic number 0) are uncounted, so this stays a
        # necessary but never sufficient condition.
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

        # Query atom order must come from the recorded cg_smarts, not the
        # directory name: only the former matches the stored coords/automorphisms.
        target_smarts, automorphisms = _cg_symmetry_for(vdg_lib_dir, db_name)
        try:
            order_pattern = _cg_order_pattern(vdg_lib_dir, db_name)
        except ValueError as e:
            if warn is not None:
                warn(f"{e}; skipping this fragment.")
            continue

        matches = lig_mol.GetSubstructMatches(
            order_pattern, uniquify=False, maxMatches=_MAX_QUERY_FRAG_MATCHES)
        if not matches:
            # Only reachable when cg_smarts is a stricter query than the
            # directory name it was labeled with.
            continue
        if len(matches) >= _MAX_QUERY_FRAG_MATCHES:
            if warn is not None:
                warn(f"Fragment {db_name!r} hit the {_MAX_QUERY_FRAG_MATCHES}-match "
                     f"cap on this ligand, so its labelings would be incomplete; "
                     f"skipping it.")
            continue

        # Both in query-ligand indices: perm_inds by CG slot (dedup only, never
        # coordinates), orig_mol_inds the site's atom set (grouping, q_atom_indices).
        instances = [(Frags.submol_from_match(lig_mol, match), tuple(match),
                      tuple(sorted(match))) for match in matches]
        grouped_sites = Frags.group_lig_sites_by_overlap(instances)
        filtered_frags[db_name] = [
            _expand_site_by_cg_automorphisms(site, automorphisms)
            for site in grouped_sites]
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
    d1  = _dist(ca_coords[0], cg_com)
    d2  = _dist(ca_coords[1], cg_com)
    d12 = _dist(ca_coords[0], ca_coords[1])
    return d1, d2, d12


def get_residue_all_atom_coords(struct, seg, chain, resnum):
    """All non-hydrogen atom coords for one residue, as float32 [N,3]; None if not found."""
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

    # Blank element columns are common in hand-edited files; fall back to the
    # atom-name convention, as split_residue_heavy_atoms does.
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
    """Concatenate all non-H atom coords from the residues in bsr_combo."""
    chunks = []
    for seg, chain, resnum in bsr_combo:
        c = get_residue_all_atom_coords(struct, seg, chain, resnum)
        if c is not None and len(c) > 0:
            chunks.append(c)
    if not chunks:
        return None
    return np.concatenate(chunks, axis=0).astype(np.float32, copy=False)


# A vdG stores only N/CA/C, so a backbone-labeled slot says nothing about the
# query residue's sidechain; these thresholds let one `bb` label replace a
# bbGLY/bbPRO split. Measured over a 193-fragment library: CB-to-CG distance at
# real backbone contacts has a 1st percentile of 4.47 A (3.4 A rejects 0.08% of
# known-good geometries, 46% of glycine-derived ones), and 0 of 1,531
# proline-derived backbone vdGs have a CG acceptor within 3.5 A of N.
BB_SLOT_SIDECHAIN_CLASH = 3.4   # A; query sidechain heavy atom vs CG
PRO_NH_DONOR_CUTOFF = 3.5       # A; backbone N vs a CG acceptor (N/O/S)


def backbone_slot_blockers(struct, bsr_combo, slot_labels, slot_resnames):
    """Per-slot obstacles to hosting a *backbone*-mediated vdG, or None.

    One (sidechain_coords, n_coord_or_None) entry per backbone-labeled slot;
    `n_coord` only for proline, which cannot donate its backbone N-H. Sidechain-
    labeled slots contribute nothing: their contact is the claim being tested."""
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
        if sc is None:      # resname outside the 20; no reference set to split on
            continue
        # Non-canonical atoms occlude too, though not declared sidechain.
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
    """Whether every backbone-labeled slot can host a CG at `cg_coords`."""
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
    """Return True if any BSR atom is within cutoff of any CG atom."""
    if bsr_atom_coords is None or len(bsr_atom_coords) == 0:
        return False
    cg_coords = np.asarray(cg_coords, dtype=np.float32)
    if cg_coords.size == 0:
        return False
    cutoff2 = np.float32(cutoff * cutoff)
    diff = bsr_atom_coords[:, None, :] - cg_coords[None, :, :]
    d2 = np.sum(diff * diff, axis=2, dtype=np.float32)
    return bool(np.any(d2 <= cutoff2))


# Deduplication                                                       #

def deduplicate_hits(match_records, min_shared_atoms=3):
    """
    Collapse hits describing the same protein-ligand interaction: same *set* of
    BSR residues, and query atom index sets overlapping by >= min_shared_atoms.
    The lower-vdg_rmsd record wins; records without q_atom_indices are kept as-is.

    Sets, not the bsr_combo string, which is ordered to match aa_bucket: one
    residue pair emits several orderings across its aa_bucket variants, so string
    equality would leave them uncollapsed (grouping outside this function must key
    on the set too). The predicate ignores frag/aa_bucket/nr vdG, so at most one
    hit survives per (BSR residue set, ligand site). Sorted by vdg_rmsd ascending.
    """
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

    sorted_recs = sorted(match_records, key=lambda r: float(r["vdg_rmsd"]))
    kept = []
    kept_atoms = []  # parallel list: frozenset | None
    kept_bsr = []    # parallel list: frozenset of residue tokens

    for rec in sorted_recs:
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


# Core per-model scoring                                              #

def score_one_model(
    pdbfile, pdb_path, lig_smiles, vdg_lib_dir, rmsd_threshold=None,
    ref_lig_mol=None, print_bsr_selection=False, contact_cutoff=None,
    vdg_lib_entries=None, deduplicate=True, min_shared_atoms=3,):
    """
    Score one query PDB/CIF file against the vdG library: fragment the query
    ligand, find candidate binding-site-residue (BSR) combos, and for each
    fragment x BSR combo search the matching (frag, subset_size, aa_bucket) npz
    bucket for the lowest-RMSD nr vdG via Kabsch alignment.

    Parameters
    ----------
    pdbfile : basename, used as a key in output records.
    pdb_path : full path to the query structure file.
    lig_smiles : SMILES of the query ligand (bond-order template).
    vdg_lib_dir : root directory of the vdG fragment library.
    rmsd_threshold : max bb+CG RMSD (Å) for a hit. None (default) derives it per
        combo from normalize_rmsd(n_atoms, "cgvdmbb") -- the same function used in
        library clustering, so generation and retrieval share a quality window.
    ref_lig_mol : rdkit Mol to compute ligand RMSD against (optional).
    print_bsr_selection : print the ProDy BSR selection string (debugging).
    contact_cutoff : if set, skip BSR combos with no BSR atom within this
        distance (Å) of any CG atom. Disabled by default; 3.8 is reasonable.
    vdg_lib_entries : set of vdg_lib_dir directory names; if None, read from disk.
    deduplicate : if True (default), run deduplicate_hits() before returning.
    min_shared_atoms : passed to deduplicate_hits (default: 3).

    Returns
    -------
    log_text : captured stdout/stderr.
    match_records : one dict per hit — pdbfile, struct_id, frag, query_frag,
        subset_size, bsr_combo, aa_bucket, vdg_index, vdg_cluster_id,
        vdg_cluster_size, vdg_rmsd, rmsd_threshold, aa_perm_idx, q_site_idx,
        q_cg_perm_idx, q_atom_indices, R00-R22, t0-t2 (db->query Kabsch
        transform). bsr_combo is ordered to match aa_bucket, so compare residue
        sets, not the string. aa_perm_idx indexes
        `vdg_npz_utils.aa_perm_indices(bucket["aa_bucket_parts"])`, never
        re-derived from resnames, and is applied as `vdg_bb[perm]` (perm[i] is the
        library slot superposed onto BSR position i); the group is inverse-closed,
        so the other reading transposes the correspondence silently.
    frags_in_lib : {library dir name: bool} -- False when the directory exists but
        its generation job did not complete.
    lig_rmsd_value : ligand RMSD vs ref_lig_mol, or None if not computed.
    filtered_frags : {library dir name: list of site groups}, indexed by q_site_idx.
    """
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

            het = struct.select("hetatm and not (water or ion or resname SEP or "
                                "resname TPO or resname MSE)")
            if het is None:
                print(f"[ERROR] ({pdbfile}) Could not find ligand HETATM atoms.", flush=True)
                return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags

            lignames = list(set(het.getResnames()))
            if len(lignames) != 1:
                n_heavy = lig_mol_noH.GetNumAtoms()
                matches = []
                for name in lignames:
                    sel = het.select(f"resname {name}")
                    if sel is None:
                        continue
                    resindices = np.unique(sel.getResindices())
                    sel_heavy = sel.select("not element H D")
                    n_sel_heavy = len(sel_heavy) if sel_heavy is not None else 0
                    if (len(resindices) > 0
                            and n_sel_heavy == n_heavy * len(resindices)):
                        matches.append(name)
                if len(matches) == 1:
                    lignames = matches
            if len(lignames) != 1:
                print(f"[ERROR] ({pdbfile}) Expected exactly one ligand resname, "
                      f"got {lignames}; if the structure has partial occupancies or "
                      f"multiple ligands, filter before running the hit finder.",
                      flush=True)
                return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags
            ligname = lignames[0]

            # query_frag_map keeps the query SMILES separate from the library
            # dir name; they differ under cross-library querying.
            filtered_frags, query_frag_map, frags_in_lib = match_library_frags_to_query(
                lig_mol_noH, vdg_lib_dir, vdg_lib_entries,
                frags_in_lib=frags_in_lib,
                warn=lambda m: print(f"[WARNING] ({pdbfile}) {m}", flush=True),)

            if not filtered_frags:
                return buf.getvalue(), match_records, frags_in_lib, lig_rmsd_value, filtered_frags

            all_bsr_combos = dock.get_bsr_combinations(
                struct, ligname, quiet=not print_bsr_selection, pdbfile=pdbfile)

            # Keyed on the BSR combo alone, and each costs a ProDy selection.
            # The fragment loop is outer, so without these every combo would be
            # re-selected once per library fragment.
            bsr_atom_coords_cache = {}
            bb_blockers_cache = {}

            for frag_name, grouped_sites in filtered_frags.items():
                query_frag = query_frag_map.get(frag_name, frag_name)
                cg_order_smarts = cg_atom_order_smarts(vdg_lib_dir, frag_name)

                try:
                    grouped_q_cg_perms = []
                    for site in grouped_sites:
                        perms = []
                        for sub, perm_inds, orig_mol_inds in site:
                            q_cg_coords = np.asarray(
                                dock.get_query_cg_coords(sub, cg_order_smarts), np.float32
                            )
                            q_cg_com = cg_center(q_cg_coords)
                            # CG atoms that could accept a backbone N-H, in slot
                            # order. From `sub`: the library pattern may leave an
                            # element open ('[N,O]') and does not know which atom
                            # this permutation puts in each slot.
                            q_acceptor = np.array(
                                [a.GetSymbol() in ("N", "O", "S")
                                 for a in sub.GetAtoms()], dtype=bool)
                            perms.append((q_cg_coords, q_cg_com, q_acceptor,
                                          perm_inds, orig_mol_inds))
                        if perms:
                            grouped_q_cg_perms.append(perms)
                except ValueError as e:
                    # The stored CG atom order does not describe this query
                    # fragment; per fragment, so skip rather than fail the model.
                    print(f"[WARNING] ({pdbfile}) Fragment {frag_name!r}: {e}\n"
                          f"  Skipping this fragment.", flush=True)
                    continue

                if not grouped_q_cg_perms:
                    continue

                for combo, bsr_combo, _bsr_AAs, coords in all_bsr_combos:
                    # Sort by AA *label*: aa_bucket is built from the sorted
                    # identities, and the query bb coords must sit in the bucket's
                    # aa_bucket_parts order (checked below) for vdg_bb[perm] to
                    # mean anything, so a label-independent key breaks it. The cost
                    # is that one residue pair gets different orderings across its
                    # aa_bucket variants, so cross-record comparison must use
                    # residue *sets* (see deduplicate_hits). bsr_resnames rides
                    # along: it is the real per-slot resname a backbone label
                    # hides, which backbone_slot_blockers needs.
                    paired = sorted(zip(combo, coords, bsr_combo, _bsr_AAs),
                                    key=lambda x: x[0])
                    (bsr_incl_bb_identities, input_bsr_bb_coords, bsr_combo,
                     bsr_resnames) = zip(*paired)
                    bsr_incl_bb_identities = list(bsr_incl_bb_identities)
                    bsr_combo = list(bsr_combo)
                    input_bsr_bb_coords = np.asarray(input_bsr_bb_coords, np.float32)
                    subset_size = input_bsr_bb_coords.shape[0]
                    bb_flat = input_bsr_bb_coords.reshape(-1, 3)
                    bb_batch = bb_flat[None]   # kabsch_ssd wants Y as [M, N, 3]
                    N_bb = bb_flat.shape[0]
                    input_bsr_ca_coords = input_bsr_bb_coords[:, 1, :]
                    aa_bucket = vdg_npz.make_aa_bucket(bsr_incl_bb_identities)

                    total_q_perms = sum(len(site) for site in grouped_q_cg_perms)
                    if total_q_perms == 0:
                        continue

                    first_cg = grouped_q_cg_perms[0][0][0]
                    N_cg = first_cg.shape[0]
                    n_atoms = N_bb + N_cg
                    effective_rmsd = (normalize_rmsd(n_atoms, "cgvdmbb")
                                      if rmsd_threshold is None else rmsd_threshold)
                    # Anything these reject is provably beyond effective_rmsd.
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
                        continue

                    # Obstacles to a backbone-labeled slot, fixed per combo. Not
                    # optional, unlike the contact filter: it stands in for the
                    # bbGLY/bbPRO split the library no longer carries.
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
                        continue

                    if contact_cutoff is not None and not np.any(query_contact_ok):
                        continue

                    # Every query permutation puts the CG somewhere the backbone
                    # slots cannot host it, so no nr vdG in this bucket can match.
                    if not np.any(query_bb_slot_ok):
                        continue

                    bucket_key = (vdg_lib_dir, frag_name, subset_size, aa_bucket)
                    bucket = bucket_cache.get(bucket_key)
                    if bucket is CACHE_MISS:
                        bucket = vdg_npz.load_vdg_bucket(
                            vdg_lib_dir, frag_name, subset_size, aa_bucket)
                        # Cached even when None, so a bucket this library does
                        # not have is not re-stat'd for every later model.
                        bucket_cache.put(bucket_key, bucket)
                    if bucket is None:
                        continue

                    cluster_ids    = bucket["cluster_id"]
                    cluster_sizes  = bucket["cluster_size"]
                    bucket_cg      = bucket["cg"]
                    bucket_bb      = bucket["bb"]
                    n_nr_vdgs    = bucket_cg.shape[0]

                    # Slots are interchangeable only when their *bucket labels*
                    # match: a backbone label is a role, not a residue, so
                    # per-nr-vdG resnames would invent swaps (a second ASP) or drop
                    # the vdG (a modified residue). Bucket-level, so computed once.
                    bucket_parts = bucket["aa_bucket_parts"]
                    if list(bucket_parts) != list(bsr_incl_bb_identities):
                        print(f"[WARNING] ({pdbfile}) Bucket {aa_bucket} stores slot "
                              f"labels {list(bucket_parts)} but was looked up for "
                              f"{list(bsr_incl_bb_identities)}; skipping.", flush=True)
                        continue
                    resind_perms = vdg_npz.aa_perm_indices(bucket_parts)

                    for nr_idx in range(n_nr_vdgs):
                        vdg_cg      = bucket_cg[nr_idx]
                        vdg_bb      = bucket_bb[nr_idx]
                        vdg_idx     = int(nr_idx)
                        clus_id     = int(cluster_ids[nr_idx])
                        clus_size   = int(cluster_sizes[nr_idx])

                        if vdg_cg.shape[0] != N_cg:
                            continue

                        vdg_cg_com = cg_center(vdg_cg)

                        best_rmsd = best_aa_perm_idx = best_q_site_idx = None
                        best_q_cg_perm_idx = best_R = best_t = None

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

                            # Backbone-only SSD is a lower bound on the full
                            # bb+CG SSD, so anything it rejects cannot pass.
                            vdg_perm_bb = res_subset_bb.reshape(-1, 3)
                            if (np.sqrt(kabsch_ssd(vdg_perm_bb, bb_batch)[0] / n_atoms)
                                    > effective_rmsd):
                                continue

                            Y_sub = Y[idxs]

                            db_bb_and_cg = np.concatenate(
                                (vdg_perm_bb, vdg_cg), axis=0).astype(np.float32,
                                                                     copy=False)
                            if db_bb_and_cg.shape[0] != n_atoms:
                                continue

                            # Scan on ssd alone; the rotation is recomputed below
                            # only for a candidate that beats the incumbent.
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

                                    # Stop at an excellent match to limit combinatorics; best_rmsd is not necessarily the global minimum.
                                    if best_rmsd <= EXCELLENT_MATCH_CUTOFF:
                                        break

                        if best_rmsd is not None:
                            struct_id = _struct_id_from_pdbfile(pdbfile)
                            R_rounded = np.round(best_R, 4)
                            t_rounded = np.round(best_t, 4)

                            # orig_mol_inds (last perm field): original query
                            # ligand atom indices.
                            best_orig_mol_inds = grouped_q_cg_perms[best_q_site_idx][best_q_cg_perm_idx][-1]
                            q_atom_indices_str = ";".join(
                                str(i) for i in sorted(best_orig_mol_inds))

                            rec = dict(
                                pdbfile=pdbfile,
                                struct_id=struct_id,
                                frag=frag_name,
                                query_frag=query_frag,
                                subset_size=subset_size,
                                bsr_combo=bsr_label_to_string(bsr_combo),
                                aa_bucket=aa_bucket,
                                vdg_index=vdg_idx,
                                vdg_cluster_id=clus_id,
                                vdg_cluster_size=clus_size,
                                vdg_rmsd=f"{best_rmsd:.4f}",
                                rmsd_threshold=f"{effective_rmsd:.4f}",
                                aa_perm_idx=int(best_aa_perm_idx),
                                q_site_idx=int(best_q_site_idx),
                                q_cg_perm_idx=int(best_q_cg_perm_idx),
                                q_atom_indices=q_atom_indices_str,)
                            for i in range(3):
                                for j in range(3):
                                    rec[f"R{i}{j}"] = f"{R_rounded[i, j]:.4f}"
                            for k in range(3):
                                rec[f"t{k}"] = f"{t_rounded[k]:.4f}"

                            match_records.append(rec)

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


# multiprocessing worker wrapper                                     #

def process_work_item(args):
    """Pool worker: unpacks the score_one_model arguments (pdbfile, pdb_path,
    lig_smiles, vdg_lib_dir, rmsd_threshold, ref_lig_mol, print_bsr_selection,
    contact_cutoff, deduplicate, min_shared_atoms) and returns ``(pdbfile,
    log_text, lig_rmsd_value, match_records, frags_in_lib, error_text)``; errors
    are returned, not raised, so one bad model cannot abort the batch.
    ``vdg_lib_entries`` comes from ``init_worker`` once per process rather than
    being pickled into every task."""
    (pdbfile, pdb_path, lig_smiles, vdg_lib_dir, rmsd_threshold, ref_lig_mol,
     print_bsr_selection, contact_cutoff, deduplicate, min_shared_atoms) = args
    vdg_lib_entries = _WORKER_MOL_CACHE_ENTRIES

    # Robust to direct use without an explicit init_worker call.
    if vdg_lib_entries is None:
        vdg_lib_entries = frozenset(os.listdir(vdg_lib_dir))
        init_worker(vdg_lib_entries)
    try:
        log_text, match_records, frags_in_lib, lig_rmsd_value, _ = score_one_model(
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
