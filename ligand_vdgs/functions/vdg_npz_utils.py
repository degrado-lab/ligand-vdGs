# vdg_npz_utils.py

import os
import time
import zipfile

import numpy as np
import prody as pr

from ligand_vdgs.functions import utils
from ligand_vdgs.functions.vdg_struct_utils import (
    NONCG_LIGAND_OCC, VDM_OCC, cg_slot_occupancy)
from ligand_vdgs.functions.Frags import check_vdg_job_status


def parse_pdb_with_retry(pdb_path, attempts=1, delay=0.5):
    """``pr.parsePDB`` with an optional linear backoff, raising on final failure.
    ``attempts`` > 1 is for callers that parse the same file from several
    concurrent processes, where a read can fail transiently on a network
    filesystem.
    """
    for i in range(attempts):
        try:
            return pr.parsePDB(str(pdb_path))
        except Exception:
            if i == attempts - 1:
                raise
            time.sleep(delay * (i + 1))


def parse_pdb_or_none(pdb_path, consequence, attempts=1, delay=0.5):
    """Parse a PDB, or warn and return None.

    ``consequence`` names what the caller loses when the parse fails ("no
    non-CG atoms are added"), so the warning still says which part of the
    output is missing -- the reason these call sites are worth sharing is the
    parse/except/None-check triple, not the message.
    """
    try:
        struct = parse_pdb_with_retry(pdb_path, attempts=attempts, delay=delay)
    except Exception as e:
        print(f"[WARNING] parsePDB failed for {pdb_path}: {e}; {consequence}.")
        return None
    if struct is None:
        print(f"[WARNING] parsePDB returned nothing for {pdb_path}; {consequence}.")
        return None
    return struct


# Below this ratio, rotation about the anchors' primary axis is too sensitive to
# PDB coordinate precision to extrapolate safely to unmatched ligand atoms.
_MIN_ROTATION_ANCHOR_RATIO = 1e-3


def _rotation_anchor_ratio(coords):
    """Return the secondary/primary spread of an alignment point set."""
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or not np.isfinite(coords).all():
        return 0.0
    singular_values = np.linalg.svd(coords - coords.mean(axis=0), compute_uv=False)
    if singular_values.size < 2 or singular_values[0] <= np.finfo(float).tiny:
        return 0.0
    return float(singular_values[1] / singular_values[0])


def _append_full_ligand_from_parent(all_coords, names, resnames, resnums, chids, segnames, 
    elements, occupancies, parent_pdb_path, cg_coords, cg_chain, cg_resnum, cg_names,
    pdb):
    """Add non-CG ligand atoms from parent PDB, mapped into cg_coords frame.

    ``pdb`` is the already-parsed parent (the caller parses it once and shares it
    with _append_backbone_carbonyls); ``parent_pdb_path`` is kept for messages.

    They are marked NONCG_LIGAND_OCC, which sits below the CG band (see the
    occupancy protocol in vdg_struct_utils) -- not above it, where a reader
    testing `occupancy > 2.9` would pick them up as CG atoms.
    """
    if not parent_pdb_path or pdb is None:
        return

    cg_coords = np.asarray(cg_coords, float)
    cg_chain, cg_resnum = str(cg_chain), int(cg_resnum)
    cg_names  = np.asarray(cg_names)

    # Map CG atoms in the *current* frame: (chain,resnum,name) -> coord.
    #
    # Atom names are not unique within a PDB residue in practice (2y1x's SAH
    # carries two atoms named N at a blank altloc), and this key is the only thing
    # tying a parent atom to its current-frame coordinate. Altloc variants, which
    # legitimately share a name, are not the case handled here: the parsePDB above
    # takes ProDy's default and so keeps only altloc 'A' and blank. Passing
    # altloc='all' there would make this guard fire on ordinary altlocs, which
    # should be resolved by occupancy (vdg_struct_utils._pick_best_altloc) instead. A collision would pair a parent atom
    # with the wrong CG coordinate and tilt the single fit that maps every non-CG
    # atom, so refuse the ligand rather than letting the last writer win.
    cg_atom_coords = {}
    duplicate_cg_keys = set()
    for coord, name in zip(cg_coords, cg_names):
        key = (cg_chain, cg_resnum, str(name))
        if key in cg_atom_coords:
            duplicate_cg_keys.add(key)
        cg_atom_coords[key] = np.asarray(coord, float)
    if duplicate_cg_keys:
        print(f"[WARNING] _append_full_ligand: CG atoms of {parent_pdb_path} share "
              f"names {sorted(duplicate_cg_keys)}, so parent atoms cannot be matched "
              "to them unambiguously; no non-CG atoms are added.")
        return

    cg_keys = set(cg_atom_coords.keys())

    # One residue: a vdG's CG is cut from a single ligand residue, so there is one
    # selection to make and one fit to solve.
    ch, rn = cg_chain, cg_resnum
    sel = pdb.select(f"chain {ch} and resnum {rn}")
    if sel is None:
        return

    coords_pdb = sel.getCoords()
    anames     = sel.getNames()
    rnames     = sel.getResnames()
    rnums      = sel.getResnums()
    chs        = sel.getChids()
    elems      = sel.getElements()
    segs       = sel.getSegnames()

    # Build matched pairs: parent-PDB CG atoms <-> current CG coords
    match_pdb, match_cg = [], []
    matched_keys, ambiguous_keys = set(), set()
    for c_pdb, an, rn_, ch_ in zip(coords_pdb, anames, rnums, chs):
        key = (str(ch_), int(rn_), str(an))
        if key not in cg_atom_coords:
            continue
        if key in matched_keys:  # two parent atoms with one name (see above)
            ambiguous_keys.add(key)
            continue
        matched_keys.add(key)
        match_pdb.append(np.asarray(c_pdb, float))
        match_cg.append(cg_atom_coords[key])
    if ambiguous_keys:
        print(f"[WARNING] _append_full_ligand: {parent_pdb_path} chain {ch} resnum "
              f"{rn} has repeated atom names {sorted(ambiguous_keys)}; which atom is "
              "the CG atom is undecidable, so non-CG atoms are skipped.")
        return

    if len(match_pdb) < 3:
        print(f"[WARNING] _append_full_ligand: only {len(match_pdb)} CG atom(s) matched "
              f"in {parent_pdb_path} for chain {ch} resnum {rn}; skipping non-CG atoms.")
        return

    X = np.asarray(match_pdb, float).reshape(-1, 3)  # mobile: [N, 3] uses fast path
    Y_points = np.asarray(match_cg, float).reshape(-1, 3)
    source_ratio = _rotation_anchor_ratio(X)
    target_ratio = _rotation_anchor_ratio(Y_points)
    if min(source_ratio, target_ratio) <= _MIN_ROTATION_ANCHOR_RATIO:
        print(
            f"[WARNING] _append_full_ligand: matched CG atoms are collinear or nearly "
            f"collinear in {parent_pdb_path} for chain {ch} resnum {rn} "
            f"(anchor ratios {source_ratio:.2e}, {target_ratio:.2e}); rotation about "
            "their primary axis is underdetermined, so non-CG atoms are skipped."
        )
        return

    Y = Y_points.reshape(1, -1, 3)  # target: cg_coords frame
    R_all, t_all, _ = utils.kabsch(X, Y)
    R, t = R_all[0], t_all[0]
    coords_mapped = coords_pdb @ R + t

    for c_map, an, rn_, ch_, el, sg, rn_name in zip(
            coords_mapped, anames, rnums, chs, elems, segs, rnames):
        key = (str(ch_), int(rn_), str(an))
        if key in cg_keys:  # skip CG atoms already present
            continue
        all_coords.append(c_map)
        names.append(str(an).strip())
        resnames.append(str(rn_name).strip())
        resnums.append(int(rn_))
        chids.append(str(ch_))
        segnames.append(str(sg).strip())
        elements.append((str(el).strip() or "C"))
        occupancies.append(NONCG_LIGAND_OCC)


# A backbone triplet re-selected from its own parent should superpose exactly on
# the stored one. Anything above this is a mismatched structure, not noise.
_CARBONYL_FIT_TOLERANCE = 0.1   # Angstrom


def _append_backbone_carbonyls(all_coords, names, resnames, resnums, chids, segnames,
    elements, occupancies, parent_pdb_path, vdm_bb_coords, scrr_seg, scrr_chain,
    scrr_resnum, scrr_resname, struct):
    """Add each vdM's backbone carbonyl O, re-derived from the parent PDB.

    ``struct`` is the already-parsed parent, shared with
    _append_full_ligand_from_parent so the biounit is read once per vdG.

    The npz stores N/CA/C only, deliberately: those three atoms are what the
    hit-finding RMSD is computed over, and storing O would invite the reading
    that it took part. The O is recovered here instead by superposing the
    parent residue's own N/CA/C onto the stored triplet and carrying its O
    through the same transform -- so this is a *display* atom, added at write
    time, and nothing about the library or the matching changes.
    """
    if not parent_pdb_path:
        print("[WARNING] backbone carbonyls requested without a parent PDB; skipping.")
        return
    if struct is None:
        return

    vdm_bb_coords = np.asarray(vdm_bb_coords, float)
    for v_idx in range(len(vdm_bb_coords)):
        seg, ch = scrr_seg[v_idx], scrr_chain[v_idx]
        resnum, resname = int(scrr_resnum[v_idx]), str(scrr_resname[v_idx])
        parent = np.empty((3, 3), dtype=float)
        missing = False
        for a, name in enumerate(("N", "CA", "C")):
            atom = _select_one_atom(struct, seg, ch, resnum, name)
            if atom is None:
                missing = True
                break
            parent[a] = np.asarray(atom.getCoords()).reshape(3)
        oxygen = None if missing else _select_one_atom(struct, seg, ch, resnum, "O")
        if missing or oxygen is None:
            print(f"[WARNING] carbonyl O of vdM slot {v_idx} (chain {ch}, "
                  f"resnum {resnum}) not resolvable in {parent_pdb_path}; skipping it.")
            continue
        # utils.kabsch returns R, t under Y ~= X @ R + t, so the parent frame is
        # X and the stored frame is Y. Applying R.T here would be a silent,
        # plausible-looking bug. Y must be batched, hence the [None].
        R, t, ssd = utils.kabsch(parent, vdm_bb_coords[v_idx][None])
        # The same three atoms of the same residue in two frames should fit
        # exactly. A real residual means the parent PDB is not the structure this
        # vdG came from, and the rotation would then place O somewhere arbitrary.
        rmsd = float(np.sqrt(max(float(ssd[0]), 0.0) / 3.0))
        if rmsd > _CARBONYL_FIT_TOLERANCE:
            print(f"[WARNING] vdM slot {v_idx} (chain {ch}, resnum {resnum}) backbone "
                  f"does not match {parent_pdb_path} (fit RMSD {rmsd:.2f} A); "
                  "skipping its carbonyl O.")
            continue
        o_coord = np.asarray(oxygen.getCoords()).reshape(1, 3) @ R[0] + t[0]

        all_coords.append(o_coord.reshape(3))
        names.append("O")
        resnames.append(resname)
        resnums.append(resnum)
        chids.append("" if ch in (None, "None") else str(ch))
        segnames.append("" if seg in (None, "", "None") else str(seg))
        elements.append("O")
        occupancies.append(VDM_OCC)


def build_vdg_atomgroup_from_npz(cg_coords, cg_names, cg_elements, cg_seg, cg_chain,
    cg_resnum, cg_resname, vdm_bb_coords, scrr_seg, scrr_chain, scrr_resnum, scrr_resname,
    include_full_ligand=False, parent_pdb_path=None, include_backbone_carbonyl=False,
    parent_struct=None,):
    """Materialize a vdG AtomGroup from the nr vdG arrays.

    CG plus vdM backbone N/CA/C always; optionally the rest of the ligand, and
    optionally each vdM's backbone carbonyl O. Both extras come from the parent
    PDB -- neither is stored in the npz -- and both read the *same* biounit, so
    it is parsed once here and shared (parsing dominates on a network
    filesystem). ``parent_struct`` lets a caller that has already parsed it (a
    member loop that just called ``rederive_member_coords``) skip even that.
    """
    cg_coords, vdm_bb_coords = np.asarray(cg_coords, float), np.asarray(vdm_bb_coords, float)
    n_cg, num_vdms = cg_coords.shape[0], vdm_bb_coords.shape[0]
    vdm_bb_flat = vdm_bb_coords.reshape(-1, 3)

    all_coords, names, resnames, resnums = [], [], [], []
    chids, segnames, elements, occupancies = [], [], [], []

    cg_names, cg_elements = np.asarray(cg_names), np.asarray(cg_elements)
    # The CG's residue is one value for the whole group, not one per atom.
    cg_resname, cg_chain, cg_resnum = str(cg_resname), str(cg_chain), int(cg_resnum)
    cg_seg = "" if cg_seg in (None, "", "None") else str(cg_seg)

    # CG: one occupancy per slot index, encoding CG atom order
    for i in range(n_cg):
        all_coords.append(cg_coords[i])
        names.append(str(cg_names[i]))
        resnames.append(cg_resname)
        resnums.append(cg_resnum)
        chids.append(cg_chain)
        segnames.append(cg_seg)
        elem_i = cg_elements[i]
        elements.append(str(elem_i) if elem_i not in (None, "", " ") else "C")
        occupancies.append(cg_slot_occupancy(i))

    # One parse for both extras below. Each helper still tolerates None (it just
    # adds nothing), so a bad parent PDB degrades the same way it did when they
    # parsed it themselves.
    if parent_struct is None and parent_pdb_path and (
            include_full_ligand or include_backbone_carbonyl):
        parent_struct = parse_pdb_or_none(
            parent_pdb_path, "no parent-derived atoms are added")

    # Optional: ligand atoms (non-CG) from parent PDB, mapped to the same frame
    if include_full_ligand and parent_pdb_path is not None:
        _append_full_ligand_from_parent(all_coords, names, resnames, resnums, chids,
            segnames, elements, occupancies, parent_pdb_path, cg_coords, cg_chain,
            cg_resnum, cg_names, parent_struct,)

    # vdM backbone atoms (N, CA, C)
    scrr_seg, scrr_chain = np.asarray(scrr_seg), np.asarray(scrr_chain)
    scrr_resnum, scrr_resname = np.asarray(scrr_resnum), np.asarray(scrr_resname)
    backbone_atom_names, backbone_elements = ["N", "CA", "C"], ["N", "C", "C"]

    for v_idx in range(num_vdms):
        seg = "" if scrr_seg[v_idx] in (None, "", "None") else str(scrr_seg[v_idx])
        ch, resnum, resname = str(scrr_chain[v_idx]), int(scrr_resnum[v_idx]), str(
            scrr_resname[v_idx])
        bb = vdm_bb_coords[v_idx]
        if bb.shape != (3, 3):
            raise ValueError(f"vdM bb wrong shape at {v_idx}: {bb.shape}")
        for a_name, a_elem, a_coord in zip(backbone_atom_names, backbone_elements, bb):
            all_coords.append(a_coord)
            names.append(a_name)
            resnames.append(resname)
            resnums.append(resnum)
            chids.append(ch)
            segnames.append(seg)
            elements.append(a_elem)
            occupancies.append(VDM_OCC)

    if include_backbone_carbonyl:
        _append_backbone_carbonyls(all_coords, names, resnames, resnums, chids,
            segnames, elements, occupancies, parent_pdb_path, vdm_bb_coords,
            scrr_seg, scrr_chain, scrr_resnum, scrr_resname, parent_struct)

    all_coords = np.asarray(all_coords, float)
    ag = pr.AtomGroup("vdg_nr")
    ag.setCoords(all_coords)
    ag.setNames(np.asarray(names))
    ag.setResnames(np.asarray(resnames))
    ag.setResnums(np.asarray(resnums, int))
    ag.setChids(np.asarray(chids))
    ag.setSegnames(np.asarray(segnames))
    ag.setElements(np.asarray(elements))
    ag.setOccupancies(np.asarray(occupancies, float))

    cg_vdmbb_flat = np.vstack([cg_coords, vdm_bb_flat])
    return ag, cg_vdmbb_flat

NR_FIELDS = ("cg_coords", "cg_names", "cg_seg", "cg_chain",
    "cg_resnum", "cg_resname", "vdm_bb_coords", "scrr_seg", "scrr_chain",
    "scrr_resnum", "scrr_resname")


def nr_build_kwargs(data, idx, include_full_ligand=False, pdb_dir=None,
                          include_backbone_carbonyl=False):
    """kwargs for ``build_vdg_atomgroup_from_npz`` from nr vdG ``idx`` of a bucket.

    Just unwraps the ``nr_``-prefixed columns; the scrr_* arrays come back in
    stored library slot order, which the hit writer reorders afterwards.
    """
    kwargs = {f: data[f"nr_{f}"][idx] for f in NR_FIELDS}
    # Bucket-level, not per nr vdG (see _write_bucket_npz).
    kwargs["cg_elements"] = data["cg_elements"]
    kwargs["parent_pdb_path"] = resolve_parent_pdb_path(
        data, str(data["nr_parent_biounit"][idx]), pdb_dir=pdb_dir)
    kwargs["include_full_ligand"] = include_full_ligand
    kwargs["include_backbone_carbonyl"] = include_backbone_carbonyl
    return kwargs


def resolve_parent_pdb_path(data, biounit, pdb_dir=None):
    """Absolute path of a parent structure from its biounit stem.

    Buckets store the stem ("1f8s") plus the build-time directory once per file,
    rather than a full path per record -- see the note in `_write_bucket_npz`.
    The layout mirrors how generation built the path:
    ``<pdb_dir>/<biounit[1:3].lower()>/<biounit>.pdb``.
 
    """
    biounit = str(biounit).strip()
    if not biounit:
        return ""
    if pdb_dir is None:
        pdb_dir = str(data["parent_pdb_dir"])
    if not pdb_dir:
        return ""
    return os.path.join(pdb_dir, biounit[1:3].lower(), biounit + ".pdb")


def _resnum_selstr(resnum):
    """ProDy requires backtick-quoting for negative resnums in selection strings."""
    return f"`{resnum}`" if resnum < 0 else str(resnum)


def name_selstr(name):
    """ProDy selection clause matching one literal atom name.

    Ligand atom names may contain characters the selection grammar treats as
    operators -- `N9+` parses as `N9` followed by a dangling `+` and raises
    SelectionError. Backticks make the name literal; they are a no-op for
    ordinary names, and they also suppress wildcard expansion, which is what we
    want here since these names come from match lists, not user patterns.
    """
    return f"`{name}`"


def _select_one_atom(struct, seg, chain, resnum, name):
    """Select a single named atom, preferring highest occupancy then altloc 'A'
    among ties (matching vdg_struct_utils._pick_best_altloc). Always returns an
    Atom (or None), never a multi-atom Selection.
    """
    seg = "" if seg in (None, "", "None") else str(seg)
    seg_clause = f"segment {seg} and " if seg else ""
    sel = struct.select(
        f"{seg_clause}chain {chain} and resnum {_resnum_selstr(resnum)} "
        f"and name {name_selstr(name)}")
    if sel is None or sel.numAtoms() == 0:
        return None
    if sel.numAtoms() == 1:
        return sel[0]
    atoms = list(sel)
    max_occ = max(a.getOccupancy() for a in atoms)
    top = [a for a in atoms if a.getOccupancy() == max_occ]
    if len(top) == 1:
        return top[0]
    for a in top:
        if a.getAltloc() == 'A':
            return a
    return top[0]


def rederive_member_coords(pdbpath, cg_seg, cg_chain, cg_resnum, cg_names,
    scrr_seg, scrr_chain, scrr_resnum, parsed_pdb=None):
    """Re-select a cluster member's CG and vdM-backbone atoms directly from its
    parent PDB by (seg, chain, resnum, name), returning their coordinates in the
    PDB's original frame (not yet aligned to the cluster's nr vdG).

    Coordinates are intentionally not persisted per member in the nr_vdgs npz
    (see clus_and_deduplicate_vdgs._extract_member_identity); this re-derives
    them on demand from the lightweight identity fields that are stored.

    Returns ``(cg_coords, vdm_bb_coords)`` or ``(None, None)`` if any named atom
    is missing (e.g. the PDB mirror changed, or the atom was never resolved).
    """
    struct = parsed_pdb
    if struct is None:  # short-circuit: an already-parsed structure is reused as is
        struct = parse_pdb_or_none(
            pdbpath, "this member's coordinates cannot be re-derived")
        if struct is None:
            return None, None

    n_cg = len(cg_names)
    cg_resnum = int(cg_resnum)  # one residue per CG, so these are scalars
    cg_coords = np.empty((n_cg, 3), dtype=float)
    for i in range(n_cg):
        atom = _select_one_atom(struct, cg_seg, cg_chain, cg_resnum, cg_names[i])
        if atom is None:
            print(f"[WARNING] rederive_member_coords: CG atom {cg_names[i]} "
                  f"(chain {cg_chain}, resnum {cg_resnum}) not found in {pdbpath}.")
            return None, None
        cg_coords[i] = np.asarray(atom.getCoords()).reshape(3)

    num_vdms = len(scrr_chain)
    vdm_bb_coords = np.empty((num_vdms, 3, 3), dtype=float)
    for v in range(num_vdms):
        for a, name in enumerate(("N", "CA", "C")):
            atom = _select_one_atom(struct, scrr_seg[v], scrr_chain[v], int(scrr_resnum[v]), name)
            if atom is None:
                print(f"[WARNING] rederive_member_coords: backbone atom {name} of vdM slot "
                      f"{v} (chain {scrr_chain[v]}, resnum {scrr_resnum[v]}) not found in {pdbpath}.")
                return None, None
            vdm_bb_coords[v, a] = np.asarray(atom.getCoords()).reshape(3)

    return cg_coords, vdm_bb_coords


def make_aa_bucket(aa_list):
    """
    Turn a list/tuple of slot labels (resnames, BB_LABELS, 'X') into the bucket name.

    Joins on '_', which is why no label may contain one (see BB_LABELS).

    Example:
        ['bb', 'ALA']    -> 'ALA_bb'  (sorted)
        ['GLU', 'GLU']   -> 'GLU_GLU'
        ['bb', 'SER']    -> 'SER_bb'
    """
    return "_".join(sorted(aa_list))

# What a truncated, half-written or bit-rotted bucket raises. zipfile.BadZipFile
# derives straight from Exception, not OSError, so it is not covered by any of
# the others and has to be named: a truncated file raises it out of np.load,
# while a flipped byte inside a member raises it (CRC-32) only when that array
# is actually read, i.e. from inside the dict construction below. Both points
# have to sit inside the same try.
CORRUPT_NPZ_ERRORS = (OSError, EOFError, zipfile.BadZipFile,
                      KeyError, ValueError, TypeError)


def load_bucket_npz(npz_path):
    """Every array of one bucket npz, read eagerly into a dict.

    Eager because ``np.load`` decompresses a member only when it is indexed, so
    a bit-rotted array raises out of whatever consumer happened to touch it,
    arbitrarily far from the file that is actually broken. Reading here puts
    every corruption at one place -- and closes the file, rather than leaving an
    open handle per bucket walked.

    Returns None (after a warning) if the file is missing or unreadable; the
    resulting dict is a drop-in for an ``np.load`` result for readers that only
    index it.
    """
    if not os.path.isfile(npz_path):
        return None
    try:
        with np.load(npz_path) as data:
            return {key: data[key] for key in data.files}
    except CORRUPT_NPZ_ERRORS as e:
        print(f"[WARNING] Could not load {npz_path}: {type(e).__name__}: {e}")
        return None


# (frag_name, vdg_lib_dir) -> completed?  Memoizes the *result*, not just the
# warning: the log read is an NFS stat+read per call, and a completed fragment
# used to pay it on every bucket load because only failures were remembered. A
# library is static for the life of a read-path process, so a cached answer
# cannot go stale within one.
_vdg_job_status_cache = {}
# Kept as the set of fragments already warned about, for callers/tests that
# inspect it.
_warned_incomplete_frags = set()


def _vdg_job_completed(frag_name, vdg_lib_dir):
    key = (frag_name, vdg_lib_dir)
    if key not in _vdg_job_status_cache:
        _vdg_job_status_cache[key] = check_vdg_job_status(frag_name, vdg_lib_dir)
    return _vdg_job_status_cache[key]


def load_vdg_bucket(vdg_lib_dir, frag_name, subset_size, aa_bucket):
    """Load one vdG npz bucket; returns short-keyed dict or None if missing/unreadable.

    A corrupt bucket is warned about and skipped rather than raised on: a
    library is thousands of buckets and one unreadable file should not abort a
    whole hit-finding run. Note this is the opposite call from
    ``load_cg_symmetry``, which raises -- there, degrading quietly would mean
    matching against the wrong CG atom correspondence.

    Also warns (once per fragment) if the fragment's vdG-generation job never
    completed -- the data still loads, since a bucket npz is written before
    clustering finishes and can be legitimately usable, but a caller should
    know the fragment may be a partial build.
    """
    if (frag_name, vdg_lib_dir) not in _warned_incomplete_frags and \
            not _vdg_job_completed(frag_name, vdg_lib_dir):
        _warned_incomplete_frags.add((frag_name, vdg_lib_dir))
        print(f"[WARNING] Fragment {frag_name!r} in {vdg_lib_dir} has no completed "
              f"vdG-generation job (no 'Job completed.' in its log); loading anyway.")

    npz_path = vdg_npz_path(vdg_lib_dir, frag_name, subset_size, aa_bucket)
    if not os.path.isfile(npz_path):
        return None

    try:
        with np.load(npz_path) as data:
            return dict(
                aa_bucket_parts=[str(x) for x in data["aa_bucket_parts"]],
                cluster_id=data["cluster_id"].astype(np.int32),
                cluster_size=data["cluster_size"].astype(np.int32),
                cg=data["nr_cg_coords"].astype(np.float32),
                bb=data["nr_vdm_bb_coords"].astype(np.float32),
                resnames=data["nr_scrr_resname"],  # dtype <U4/string
                # Per-slot provenance (SLOT_* in vdg_struct_utils).
                slot_flags=data["nr_slot_flag"].astype(np.int8),
            )
    except CORRUPT_NPZ_ERRORS as e:
        print(f"[WARNING] Could not load {npz_path}: {type(e).__name__}: {e}")
        return None

def cluster_member_indices(data, cluster_id):
    """Row indices of one cluster's members, for callers that sample before
    materializing the per-member dicts (see ``load_cluster_members(indices=...)``)."""
    return np.nonzero(data["mem_cluster_id"] == cluster_id)[0]


def load_cluster_members(data, cluster_id, pdb_dir=None, indices=None):
    """Per-member identity rows for one cluster, **excluding its nr vdG**.

    A "member" is a clustered observation that is not the non-redundant vdG, so

        cluster_size == 1 + len(load_cluster_members(...))

    ``indices`` restricts the result to those member rows (a subset of
    ``cluster_member_indices``), so a sampling caller does not build a dict per
    member of the whole cluster and then throw nearly all of them away.
    """
    if indices is None:
        indices = cluster_member_indices(data, cluster_id)
    if len(indices) == 0:
        return []

    members = []
    for i in indices:
        members.append(dict(
            pdbpath=resolve_parent_pdb_path(
                data, str(data["mem_parent_biounit"][i]), pdb_dir=pdb_dir),
            cg_names=data["mem_cg_names"][i], cg_elements=data["cg_elements"],
            cg_seg=data["mem_cg_seg"][i], cg_chain=data["mem_cg_chain"][i],
            cg_resnum=data["mem_cg_resnum"][i],
            cg_resname=data["mem_cg_resname"][i],
            scrr_seg=data["mem_scrr_seg"][i],
            scrr_chain=data["mem_scrr_chain"][i],
            scrr_resnum=data["mem_scrr_resnum"][i],
            scrr_resname=data["mem_scrr_resname"][i],
        ))
    return members


# The CG's SMARTS and automorphism group, recorded once per fragment library so
# consumers never have to re-derive them from the library directory name. The
# directory name is only a label (`clus_and_deduplicate_vdgs.py -c`), and it is
# free to differ from the `--cg-smarts` the group was actually computed from.
CG_SYMMETRY_FILENAME = "cg_symmetry.npz"


def cg_symmetry_path(vdg_lib_dir, frag_name=None):
    """Path to a fragment's symmetry sidecar.

    Generation passes the fragment's own output dir and omits ``frag_name``;
    hit finding passes the library root plus the fragment dir name.
    """
    parts = [vdg_lib_dir] if frag_name is None else [vdg_lib_dir, frag_name]
    return os.path.join(*parts, "nr_vdgs", CG_SYMMETRY_FILENAME)


def write_cg_symmetry(vdglib_dir, cg_smarts, cg_automorphisms):
    """Record the SMARTS and automorphism group the buckets were clustered under.

    Written atomically because the subset-size 1 and 2 runs share a fragment
    directory and race to write identical content.
    """
    perms = np.asarray([list(p) for p in cg_automorphisms], dtype=np.int32)
    if perms.ndim != 2:
        raise ValueError(f"Expected a rectangular permutation table, got {perms.shape}")
    path = cg_symmetry_path(vdglib_dir)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # np.savez_compressed appends '.npz' unless the name already ends in it.
    tmp = f"{path}.{os.getpid()}.tmp.npz"
    np.savez_compressed(tmp,
                        cg_smarts=np.asarray(str(cg_smarts)),
                        cg_automorphisms=perms)
    os.replace(tmp, path)


FRAGMENT_ALIASES_FILENAME = "fragment_aliases.tsv"


def load_fragment_aliases(vdg_lib_dir):
    """Map each collapsed protonation variant to the fragment that covers it.

    Written by the scheduler scripts when a variant group is reduced to one job.
    Missing file means no collapse was recorded, which is not an error.
    """
    path = os.path.join(vdg_lib_dir, FRAGMENT_ALIASES_FILENAME)
    aliases = {}
    if not os.path.isfile(path):
        return aliases
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            alias, _, representative = line.rstrip("\n").partition("\t")
            if representative:
                aliases[alias] = representative
    return aliases


def resolve_fragment_alias(vdg_lib_dir, smiles):
    """The fragment SMILES whose library directory actually holds `smiles`.

    A charged variant's vdGs live under its neutral twin, because the neutral
    SMARTS matches both protonation states and the two were mined as one job.
    Returns `smiles` unchanged when it is not an alias.
    """
    return load_fragment_aliases(vdg_lib_dir).get(smiles, smiles)


def load_cg_symmetry(vdg_lib_dir, frag_name):
    """Return ``(cg_smarts, automorphisms)`` as recorded at generation time.

    Every generation run writes the sidecar, so a missing or unreadable one is a
    broken library, not an old one, and raises. A corrupt permutation table
    raises out of validate_atom_permutations for the same reason.
    """
    path = cg_symmetry_path(vdg_lib_dir, frag_name)
    try:
        with np.load(path) as data:
            smarts = str(data["cg_smarts"])
            raw_perms = [tuple(int(i) for i in row)
                         for row in data["cg_automorphisms"]]
    except FileNotFoundError:
        # An absent sidecar is a broken library, and the plain error names it
        # better than any rewrite would; only corruption gets re-raised below.
        raise
    except CORRUPT_NPZ_ERRORS as e:
        raise ValueError(
            f"CG symmetry sidecar {path} is unreadable "
            f"({type(e).__name__}: {e}). Falling back to the identity group "
            "would superpose a symmetric CG under the wrong atom "
            "correspondence, so this is fatal; rebuild the fragment.") from e
    return smarts, utils.validate_atom_permutations(raw_perms)


def vdg_npz_path(vdg_lib_dir, frag_name, subset_size, aa_bucket):
    return os.path.join(vdg_lib_dir, frag_name, "nr_vdgs", str(subset_size),
        f"{aa_bucket}.npz")

def apply_rigid_transform(ag, R, t):
    coords = ag.getCoords()
    ag.setCoords(coords @ np.asarray(R, float) + np.asarray(t, float))
    return ag


def aa_perm_indices(bucket_parts):
    """
    Return all residue-index permutations that only swap slots with the same label.

    Input is `aa_bucket_parts`, NOT per-nr-vdG resnames. A backbone label is a
    role ("contacts via backbone"), so two 'bb' slots are interchangeable
    regardless of their underlying residue identities -- which is exactly why the
    per-nr-vdG resnames must not be consulted here.

    No special case is needed for any label; they all group alike.

    E.g.  ['ALA', 'ALA', 'SER'] -> [[0,1,2], [1,0,2]]
          ['ASP', 'HIS']        -> [[0,1]]
          ['bb', 'bb', 'bb']    -> all 6 permutations of [0,1,2]
          ['bb', 'SER']         -> [[0,1]]  (different labels, no swap)
    """
    return utils.group_preserving_permutations(bucket_parts)
