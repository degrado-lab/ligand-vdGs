# vdg_npz_utils.py

import json
import os
import time
import zipfile

import numpy as np
from ligand_vdgs.functions import parent_db
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
    with _append_vdm_parent_atoms).

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
_VDM_FIT_TOLERANCE = 0.1   # Angstrom

# Not sidechain: the stored triplet, the carbonyl (its own flag), and OXT.
_NON_SIDECHAIN_NAMES = frozenset(("N", "CA", "C", "O", "OXT"))


def _atom_element(atom, name):
    """Element of a parent atom, falling back to the atom name when the PDB's
    element column is blank."""
    elem = str(atom.getElement() or "").strip()
    return elem if elem else (name.lstrip("0123456789")[:1].upper() or "C")


def _sidechain_atoms(struct, seg, chain, resnum):
    """Heavy sidechain atoms of one residue, at most one per atom name.

    Resolved name by name through ``_select_one_atom`` rather than as one
    selection, so altlocs go through the same occupancy/altloc tie-break the rest
    of the pipeline uses. 
    """
    seg_s = "" if seg in (None, "", "None") else str(seg)
    seg_clause = f"segment {seg_s} and " if seg_s else ""
    sel = struct.select(f"{seg_clause}chain {chain} and "
                        f"resnum {_resnum_selstr(int(resnum))}")
    if sel is None:
        return []
    wanted = []
    for raw_name, raw_elem in zip(sel.getNames(), sel.getElements()):
        name = str(raw_name).strip()
        if name in _NON_SIDECHAIN_NAMES or name in wanted:
            continue
        elem = str(raw_elem or "").strip() or name.lstrip("0123456789")[:1]
        if elem.upper() in ("H", "D"):
            continue
        wanted.append(name)
    atoms = []
    for name in wanted:
        atom = _select_one_atom(struct, seg, chain, resnum, name)
        if atom is not None:
            atoms.append((name, atom))
    return atoms


def _append_vdm_parent_atoms(all_coords, names, resnames, resnums, chids, segnames,
    elements, occupancies, parent_pdb_path, vdm_bb_coords, scrr_seg, scrr_chain,
    scrr_resnum, scrr_resname, struct, include_carbonyl=False,
    include_sidechain=False):
    """Add each vdM's backbone carbonyl O and/or its sidechain, re-derived from
    the parent PDB.

    ``struct`` is the already-parsed parent, shared with
    _append_full_ligand_from_parent so the biounit is read once per vdG.

    The npz stores N/CA/C only, deliberately: those three atoms are what the
    hit-finding RMSD is computed over. The extra atoms are recovered here.
    """
    if not (include_carbonyl or include_sidechain):
        return
    what = " and ".join([w for w, on in (("carbonyl O", include_carbonyl),
                                         ("sidechain", include_sidechain)) if on])
    if not parent_pdb_path:
        print(f"[WARNING] vdM {what} requested without a parent PDB; skipping.")
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
        if missing:
            print(f"[WARNING] backbone of vdM slot {v_idx} (chain {ch}, resnum "
                  f"{resnum}) not resolvable in {parent_pdb_path}; skipping its "
                  f"{what}.")
            continue
        # utils.kabsch returns R, t under Y ~= X @ R + t, so the parent frame is
        # X and the stored frame is Y. Applying R.T here would be a silent,
        # plausible-looking bug. Y must be batched, hence the [None].
        R, t, ssd = utils.kabsch(parent, vdm_bb_coords[v_idx][None])
        # N/CA/C must align perfectly - otherwise, it's not the correct PDB.
        rmsd = float(np.sqrt(max(float(ssd[0]), 0.0) / 3.0))
        if rmsd > _VDM_FIT_TOLERANCE:
            print(f"[WARNING] vdM slot {v_idx} (chain {ch}, resnum {resnum}) backbone "
                  f"does not match {parent_pdb_path} (fit RMSD {rmsd:.2f} A); "
                  f"skipping its {what}.")
            continue

        extras = []
        if include_carbonyl:
            oxygen = _select_one_atom(struct, seg, ch, resnum, "O")
            if oxygen is None:
                print(f"[WARNING] carbonyl O of vdM slot {v_idx} (chain {ch}, "
                      f"resnum {resnum}) not resolvable in {parent_pdb_path}; "
                      "skipping it.")
            else:
                extras.append(("O", oxygen))
        if include_sidechain:
            sidechain = _sidechain_atoms(struct, seg, ch, resnum)
            # GLY has no sidechain; anything else that has none is a stripped or
            # mismatched residue worth naming.
            if not sidechain and resname != "GLY":
                print(f"[WARNING] no sidechain heavy atoms for vdM slot {v_idx} "
                      f"({resname} chain {ch}, resnum {resnum}) in "
                      f"{parent_pdb_path}.")
            extras.extend(sidechain)

        for name, atom in extras:
            coord = np.asarray(atom.getCoords()).reshape(1, 3) @ R[0] + t[0]
            all_coords.append(coord.reshape(3))
            names.append(name)
            resnames.append(resname)
            resnums.append(resnum)
            chids.append("" if ch in (None, "None") else str(ch))
            segnames.append("" if seg in (None, "", "None") else str(seg))
            elements.append(_atom_element(atom, name))
            occupancies.append(VDM_OCC)


def build_vdg_atomgroup_from_npz(cg_coords, cg_names, cg_elements, cg_seg, cg_chain,
    cg_resnum, cg_resname, vdm_bb_coords, scrr_seg, scrr_chain, scrr_resnum, scrr_resname,
    include_full_ligand=False, parent_pdb_path=None, include_backbone_carbonyl=False,
    include_sidechain=False, parent_struct=None,):
    """Materialize a vdG AtomGroup from the nr vdG arrays.

    CG plus vdM backbone N/CA/C always; optionally the rest of the ligand, each
    vdM's backbone carbonyl O, and each vdM's sidechain. All three extras come
    from the parent PDB -- none is stored in the npz -- so re-derive here.

    The returned fit block is CG + N/CA/C regardless of the extras: they are
    display atoms and never enter an alignment.
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
            include_full_ligand or include_backbone_carbonyl or include_sidechain):
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

    if include_backbone_carbonyl or include_sidechain:
        _append_vdm_parent_atoms(all_coords, names, resnames, resnums, chids,
            segnames, elements, occupancies, parent_pdb_path, vdm_bb_coords,
            scrr_seg, scrr_chain, scrr_resnum, scrr_resname, parent_struct,
            include_carbonyl=include_backbone_carbonyl,
            include_sidechain=include_sidechain)

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
                          include_backbone_carbonyl=False, include_sidechain=False):
    """kwargs for ``build_vdg_atomgroup_from_npz`` from nr vdG ``idx`` of a bucket.

    Just unwraps the ``nr_``-prefixed columns; the scrr_* arrays come back in
    stored library slot order, which the hit writer reorders afterwards.
    """
    kwargs = {f: data[f"nr_{f}"][idx] for f in NR_FIELDS}
    # Bucket-level, not per nr vdG (see _write_bucket_npz).
    kwargs["cg_elements"] = data["cg_elements"]
    # Every mode puts this path in the output filename; only these read the file.
    reads_parent = (include_full_ligand or include_backbone_carbonyl
                    or include_sidechain)
    kwargs["parent_pdb_path"] = resolve_parent_pdb_path(
        data, str(data["nr_parent_biounit"][idx]), pdb_dir=pdb_dir,
        for_reading=reads_parent)
    kwargs["include_full_ligand"] = include_full_ligand
    kwargs["include_backbone_carbonyl"] = include_backbone_carbonyl
    kwargs["include_sidechain"] = include_sidechain
    return kwargs


PDB_DIR_ENV_VAR = "PARENT_PDBS_DIR"

# Appended to every message about a missing or unusable parent database.
PDB_DIR_HELP = (
    f"${PDB_DIR_ENV_VAR} is the parent structure database this vdG library was mined "
    "from: two-character subdirs of uncompressed PDBs, i.e. <dir>/f8/1f8s.pdb. Set it "
    f"with\n    export {PDB_DIR_ENV_VAR}=/path/to/parent_pdbs\n"
    "(or pass --pdb-dir where a script offers it).")

class ParentPdbDirError(ValueError):
    "Unusable parent PDB database." 


# Parent-dir checks already run in this process, keyed by the directory they passed
# for. resolve_parent_pdb_path is called once per record, so the check has to be
# memoised to be affordable there.
_checked_pdb_dirs = set()


def effective_parent_pdb_dir(data=None, pdb_dir=None):
    """The parent PDB database.

    Precedence: an explicit ``pdb_dir`` argument, then ``$PARENT_PDBS_DIR``,
    then the build-time ``parent_pdb_dir`` stored in the bucket. The env var exists
    because that stored value is an absolute path on the build machine, so a
    library copied elsewhere resolves every parent to a path that does not exist.

    ``data=None`` consults only the first two, returning ``(None, ...)`` when
    neither is set -- for a caller checking before it has a bucket in hand.
    """
    if pdb_dir is not None:
        return str(pdb_dir), "pdb_dir argument"
    env_dir = os.environ.get(PDB_DIR_ENV_VAR)
    if env_dir:
        return env_dir, f"${PDB_DIR_ENV_VAR}"
    if data is None:
        return None, "not yet known"
    return str(data["parent_pdb_dir"]), "parent_pdb_dir recorded in the bucket"


def require_parent_pdb_dir(data=None, pdb_dir=None):
    """``effective_parent_pdb_dir``, raising if it does not name a real directory.

    Passing directories are remembered, so repeat calls are a set lookup.

    With ``data=None`` and nothing overriding, returns None: only the bucket can
    supply a directory, so the check defers to a later call that has one.
    """
    resolved, source = effective_parent_pdb_dir(data, pdb_dir=pdb_dir)
    if resolved is None:
        return None
    if resolved in _checked_pdb_dirs:
        return resolved
    if not resolved:
        raise ParentPdbDirError(
            "No parent PDB directory is available, so parent structures cannot be "
            f"read.\n{PDB_DIR_HELP}")
    if not os.path.isdir(resolved):
        raise ParentPdbDirError(
            f"Parent PDB directory does not exist: {resolved} (from {source}).\n"
            f"{PDB_DIR_HELP}")
    # isdir alone passes a parent database with the wrong layout (flat, gzipped, .cif),
    # which then costs one warning per record and still exits 0 having written
    # only what needs no parent. Try a few stems, failing only if none resolve:
    # one absent file is a partial database (already warned about per record),
    # all absent is the wrong layout.
    stems = [] if data is None else [str(b) for b in data["nr_parent_biounit"][:3]]
    stems = [b for b in stems if b]
    if stems:
        test_probes = [parent_db.structure_path(resolved, b) for b in stems]
        if not any(os.path.isfile(test_probe) for test_probe in test_probes):
            raise ParentPdbDirError(
                f"Parent PDB database {resolved} (from {source}) holds none of "
                f"{test_probes}. Parents are resolved as "
                "<dir>/<stem[1:3]>/<stem>.pdb -- uncompressed .pdb, in "
                "two-character subdirectories -- so a flat, gzipped or mmCIF "
                "database will not work even though the directory exists.\n"
                f"{PDB_DIR_HELP}")
        # Only a probed directory is remembered: a check made before any bucket was
        # loaded saw no stems, so it cannot stand in for the layout check.
        _checked_pdb_dirs.add(resolved)
    return resolved


# Printed once per distinct problem, so a run that calls the check twice (before the
# walk, then again once a bucket supplies the library's recorded path) says it once.
_warned_missing_parent_db = set()


def missing_parent_db_message(extras, reason=None):
    """The one message both PDB writers print when no parent database is available.

    Says what the output is instead (CG + vdM N/CA/C), which flags were dropped, and
    how to get them back. Kept here, not in either script, so the two cannot drift.
    """
    lines = [
        "No parent PDB database is available, so the vdGs written here are the CG "
        "plus the vdM backbone N/CA/C only."]
    if reason:
        lines.append(f"  Reason: {reason}")
    if extras:
        subject = "it adds" if len(extras) == 1 else "they add"
        lines.append(
            f"  {', '.join(extras)} dropped: the atoms {subject} are not stored in the "
            "library -- they are re-derived at write time from the structure each vdG "
            "was mined from.")
    lines.append(
        "  To include them, point at that structure database -- two-character subdirs "
        f"of uncompressed PDBs, i.e. <dir>/f8/1f8s.pdb -- and rerun:\n"
        f"      export {PDB_DIR_ENV_VAR}=/path/to/parent_pdbs\n"
        "  or pass it per run with --pdb-dir /path/to/parent_pdbs.")
    return "\n".join(lines)


def parent_extras_available(extras, data=None, pdb_dir=None):
    """True if ``extras`` (flag names) can be re-derived from a parent PDB database.

    False, with one ``missing_parent_db_message``, if none is reachable. Shared by
    both PDB writers so their behaviour is identical: the extras are display atoms,
    so a missing database is not a reason to write nothing -- the vdG itself, CG +
    vdM N/CA/C, is in the npz and needs no parent. Callers that cannot fall back on
    stored coordinates (the members mode of materialize_vdg_pdbs, whose members store
    none) must keep calling ``require_parent_pdb_dir`` and fail instead.
    """
    if not extras:
        return False
    try:
        require_parent_pdb_dir(data, pdb_dir=pdb_dir)
    except ParentPdbDirError as err:
        # ParentPdbDirError spells out $PARENT_PDBS_DIR itself; keep only its first
        # line as the reason, since the message below gives that guidance once.
        msg = missing_parent_db_message(extras, reason=str(err).split("\n")[0])
        if msg not in _warned_missing_parent_db:
            _warned_missing_parent_db.add(msg)
            print(f"[WARNING] {msg}")
        return False
    return True


def resolve_parent_pdb_path(data, biounit, pdb_dir=None, for_reading=True):
    """Absolute path of a parent structure from its biounit stem.

    Buckets store the stem ("1f8s") plus the build-time directory once per file,
    rather than a full path per record -- see the note in `_write_bucket_npz`.
    The layout mirrors how generation built the path:
    ``<pdb_dir>/<biounit[1:3].lower()>/<biounit>.pdb``.

    Raises ``ParentPdbDirError`` (once-per-directory check, memoised) if the database is
    missing or laid out differently, so any caller that reads parents gets the check
    and its how-to for free. Pass ``for_reading=False`` where the path is only a
    display tag.
    """
    biounit = str(biounit).strip()
    if not biounit:
        return ""  # a record with no parent stem
    if for_reading and effective_parent_pdb_dir(
            data, pdb_dir=pdb_dir)[0] not in _checked_pdb_dirs:
        require_parent_pdb_dir(data, pdb_dir=pdb_dir)
    resolved, _ = effective_parent_pdb_dir(data, pdb_dir=pdb_dir)
    if not resolved:
        return ""
    return parent_db.structure_path(resolved, biounit)


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
    is missing (e.g. the PDB database changed, or the atom was never resolved).
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
            # Checked here too, not only in load_vdg_bucket
            version = bucket_schema_version(data)
            if version != BUCKET_SCHEMA_VERSION:
                raise BucketSchemaMismatch(
                    f"{npz_path} was written pre-refactor; "
                    f"rebuild the frag lib.")
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


# Neighbour-element multiset packing. 
# Multiplicity is kept -- a carbon with two carbon neighbours outside the
# match is not the same environment as one with a single carbon.
NBR_ELEM_SLOTS = ('C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'X')
_NBR_ELEM_BITS = 3          # counts 0-7, saturating; 10 slots = 30 bits
_NBR_ELEM_MAX = (1 << _NBR_ELEM_BITS) - 1
_NBR_ELEM_INDEX = {sym: i for i, sym in enumerate(NBR_ELEM_SLOTS)}


def encode_nbr_elems(symbols):
    """Element symbols of one atom's out-of-match heavy neighbours -> uint32.

    Anything outside NBR_ELEM_SLOTS counts into the 'X' slot, so a metal or a
    selenium neighbour is recorded as present rather than dropped. Counts
    saturate at 7, which no real atom reaches (at most four heavy neighbours).
    """
    code = 0
    for symbol in symbols:
        i = _NBR_ELEM_INDEX.get(str(symbol).strip().capitalize(),
                                _NBR_ELEM_INDEX['X'])
        shift = i * _NBR_ELEM_BITS
        count = (code >> shift) & _NBR_ELEM_MAX
        if count < _NBR_ELEM_MAX:
            code += 1 << shift
    return code


def decode_nbr_elems(code):
    """uint32 -> tuple of element symbols, with multiplicity.

    Canonical in NBR_ELEM_SLOTS order, which is not the alphabetical order the
    miner's debug string used; compare decoded tuples, not strings.
    """
    out = []
    for i, symbol in enumerate(NBR_ELEM_SLOTS):
        out.extend([symbol] * ((int(code) >> (i * _NBR_ELEM_BITS)) & _NBR_ELEM_MAX))
    return tuple(out)


def cg_annot_pkl_path(cg_match_dict_pkl):
    """The per-CG-atom annotation pickle that sits beside a matches pickle.
    """
    base = cg_match_dict_pkl
    if base.endswith('.pkl'):
        base = base[:-len('.pkl')]
    return f'{base}.annot.pkl'


class BucketSchemaMismatch(Exception):
    """A bucket npz written by an incompatible version of the writer.

    Deliberately outside CORRUPT_NPZ_ERRORS, which downgrades a bad file to a
    warning so one unreadable bucket cannot abort a run over thousands.
    """


BUCKET_SCHEMA_VERSION = 2
# What the per-atom annotation columns mean, recorded alongside the version so a
# later reader can tell a schema bump from a semantics change.
ANNOTATION_SCHEMA = (
    "cg_heavy_degree/cg_num_h/cg_formal_charge int8 per CG atom (-1 unreadable, "
    "degree 0 invalid); cg_nbr_elems uint32 packed heavy-neighbour element "
    "multiset outside the match (decode_nbr_elems); "
    "perception int8 (0 ccd, 1 openbabel, 2 smiles, 3 atom-name table); "
    "vdm_buried_area / vdm_shared_area float32, vdm_n_atom_pairs int16, "
    "vdm_min_heavy_dist float32 per vdM slot")


def bucket_schema_version(data):
    """The schema version recorded in an open bucket npz.

    Version 1 is any bucket written before the provenance block refactor 
    existed, so its numbers are not comparable with a current build.
    """
    if "schema" not in getattr(data, "files", []):
        return 1
    try:
        return int(json.loads(str(data["schema"]))["schema_version"])
    except Exception:
        return 1


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
            version = bucket_schema_version(data)
            if version != BUCKET_SCHEMA_VERSION: # refuse
                raise BucketSchemaMismatch(
                    f"{npz_path} was written with bucket schema version "
                    f"{version}, but this code reads version "
                    f"{BUCKET_SCHEMA_VERSION}; rebuild the fragment.")
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
            # alias, representative[, kind]; kind is informational (see
            # extract_fragment_smiles.write_fragment_aliases).
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1]:
                aliases[fields[0]] = fields[1]
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


def run_cli(main):
    """Entry-point wrapper: report a bad parent PDB database as a message, not a
    traceback. Use as ``if __name__ == "__main__": vdg_npz_utils.run_cli(main)``."""
    try:
        main()
    except ParentPdbDirError as err:
        raise SystemExit(f"[ERROR] {err}")
