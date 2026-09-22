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
    for i in range(attempts):
        try:
            return pr.parsePDB(str(pdb_path))
        except Exception:
            if i == attempts - 1:
                raise
            time.sleep(delay * (i + 1))

def parse_pdb_or_none(pdb_path, consequence, attempts=1, delay=0.5):
    try:
        struct = parse_pdb_with_retry(pdb_path, attempts=attempts, delay=delay)
    except Exception as e:
        print(f"[WARNING] parsePDB failed for {pdb_path}: {e}; {consequence}.")
        return None
    if struct is None:
        print(f"[WARNING] parsePDB returned nothing for {pdb_path}; {consequence}.")
        return None
    return struct

_MIN_ROTATION_ANCHOR_RATIO = 1e-3

def _rotation_anchor_ratio(coords):
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
    if not parent_pdb_path or pdb is None:
        return

    cg_coords = np.asarray(cg_coords, float)
    cg_chain, cg_resnum = str(cg_chain), int(cg_resnum)
    cg_names  = np.asarray(cg_names)

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

    match_pdb, match_cg = [], []
    matched_keys, ambiguous_keys = set(), set()
    for c_pdb, an, rn_, ch_ in zip(coords_pdb, anames, rnums, chs):
        key = (str(ch_), int(rn_), str(an))
        if key not in cg_atom_coords:
            continue
        if key in matched_keys:
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

    X = np.asarray(match_pdb, float).reshape(-1, 3)
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

    Y = Y_points.reshape(1, -1, 3)
    R_all, t_all, _ = utils.kabsch(X, Y)
    R, t = R_all[0], t_all[0]
    coords_mapped = coords_pdb @ R + t

    for c_map, an, rn_, ch_, el, sg, rn_name in zip(
            coords_mapped, anames, rnums, chs, elems, segs, rnames):
        key = (str(ch_), int(rn_), str(an))
        if key in cg_keys:
            continue
        all_coords.append(c_map)
        names.append(str(an).strip())
        resnames.append(str(rn_name).strip())
        resnums.append(int(rn_))
        chids.append(str(ch_))
        segnames.append(str(sg).strip())
        elements.append((str(el).strip() or "C"))
        occupancies.append(NONCG_LIGAND_OCC)

_VDM_FIT_TOLERANCE = 0.1

_NON_SIDECHAIN_NAMES = frozenset(("N", "CA", "C", "O", "OXT"))

def _atom_element(atom, name):
    elem = str(atom.getElement() or "").strip()
    return elem if elem else (name.lstrip("0123456789")[:1].upper() or "C")

def _sidechain_atoms(struct, seg, chain, resnum):
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
        R, t, ssd = utils.kabsch(parent, vdm_bb_coords[v_idx][None])
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
    cg_coords, vdm_bb_coords = np.asarray(cg_coords, float), np.asarray(vdm_bb_coords, float)
    n_cg, num_vdms = cg_coords.shape[0], vdm_bb_coords.shape[0]
    vdm_bb_flat = vdm_bb_coords.reshape(-1, 3)

    all_coords, names, resnames, resnums = [], [], [], []
    chids, segnames, elements, occupancies = [], [], [], []

    cg_names, cg_elements = np.asarray(cg_names), np.asarray(cg_elements)
    cg_resname, cg_chain, cg_resnum = str(cg_resname), str(cg_chain), int(cg_resnum)
    cg_seg = "" if cg_seg in (None, "", "None") else str(cg_seg)

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

    if parent_struct is None and parent_pdb_path and (
            include_full_ligand or include_backbone_carbonyl or include_sidechain):
        parent_struct = parse_pdb_or_none(
            parent_pdb_path, "no parent-derived atoms are added")

    if include_full_ligand and parent_pdb_path is not None:
        _append_full_ligand_from_parent(all_coords, names, resnames, resnums, chids,
            segnames, elements, occupancies, parent_pdb_path, cg_coords, cg_chain,
            cg_resnum, cg_names, parent_struct,)

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
    kwargs = {f: data[f"nr_{f}"][idx] for f in NR_FIELDS}
    kwargs["cg_elements"] = data["cg_elements"]
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

PDB_DIR_HELP = (
    f"${PDB_DIR_ENV_VAR} is the parent structure database this vdG library was mined "
    "from: two-character subdirs of uncompressed PDBs, i.e. <dir>/f8/1f8s.pdb. Set it "
    f"with\n    export {PDB_DIR_ENV_VAR}=/path/to/parent_pdbs\n"
    "(or pass --pdb-dir where a script offers it).")

class ParentPdbDirError(ValueError):
    pass

_checked_pdb_dirs = set()

def effective_parent_pdb_dir(data=None, pdb_dir=None):
    if pdb_dir is not None:
        return str(pdb_dir), "pdb_dir argument"
    env_dir = os.environ.get(PDB_DIR_ENV_VAR)
    if env_dir:
        return env_dir, f"${PDB_DIR_ENV_VAR}"
    if data is None:
        return None, "not yet known"
    return str(data["parent_pdb_dir"]), "parent_pdb_dir recorded in the bucket"

def require_parent_pdb_dir(data=None, pdb_dir=None):
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
        _checked_pdb_dirs.add(resolved)
    return resolved

_warned_missing_parent_db = set()

def missing_parent_db_message(extras, reason=None):
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
    if not extras:
        return False
    try:
        require_parent_pdb_dir(data, pdb_dir=pdb_dir)
    except ParentPdbDirError as err:
        msg = missing_parent_db_message(extras, reason=str(err).split("\n")[0])
        if msg not in _warned_missing_parent_db:
            _warned_missing_parent_db.add(msg)
            print(f"[WARNING] {msg}")
        return False
    return True

def resolve_parent_pdb_path(data, biounit, pdb_dir=None, for_reading=True):
    biounit = str(biounit).strip()
    if not biounit:
        return ""
    if for_reading and effective_parent_pdb_dir(
            data, pdb_dir=pdb_dir)[0] not in _checked_pdb_dirs:
        require_parent_pdb_dir(data, pdb_dir=pdb_dir)
    resolved, _ = effective_parent_pdb_dir(data, pdb_dir=pdb_dir)
    if not resolved:
        return ""
    return parent_db.structure_path(resolved, biounit)

def _resnum_selstr(resnum):
    return f"`{resnum}`" if resnum < 0 else str(resnum)

def name_selstr(name):
    return f"`{name}`"

def _select_one_atom(struct, seg, chain, resnum, name):
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
    struct = parsed_pdb
    if struct is None:
        struct = parse_pdb_or_none(
            pdbpath, "this member's coordinates cannot be re-derived")
        if struct is None:
            return None, None

    n_cg = len(cg_names)
    cg_resnum = int(cg_resnum)
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
    return "_".join(sorted(aa_list))

CORRUPT_NPZ_ERRORS = (OSError, EOFError, zipfile.BadZipFile,
                      KeyError, ValueError, TypeError)

def load_bucket_npz(npz_path):
    if not os.path.isfile(npz_path):
        return None
    try:
        with np.load(npz_path) as data:
            version = bucket_schema_version(data)
            if version != BUCKET_SCHEMA_VERSION:
                raise BucketSchemaMismatch(
                    f"{npz_path} was written pre-refactor; "
                    f"rebuild the frag lib.")
            return {key: data[key] for key in data.files}
    except CORRUPT_NPZ_ERRORS as e:
        print(f"[WARNING] Could not load {npz_path}: {type(e).__name__}: {e}")
        return None

_vdg_job_status_cache = {}
_warned_incomplete_frags = set()

def _vdg_job_completed(frag_name, vdg_lib_dir):
    key = (frag_name, vdg_lib_dir)
    if key not in _vdg_job_status_cache:
        _vdg_job_status_cache[key] = check_vdg_job_status(frag_name, vdg_lib_dir)
    return _vdg_job_status_cache[key]

NBR_ELEM_SLOTS = ('C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'X')
_NBR_ELEM_BITS = 3
_NBR_ELEM_MAX = (1 << _NBR_ELEM_BITS) - 1
_NBR_ELEM_INDEX = {sym: i for i, sym in enumerate(NBR_ELEM_SLOTS)}

def encode_nbr_elems(symbols):
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
    out = []
    for i, symbol in enumerate(NBR_ELEM_SLOTS):
        out.extend([symbol] * ((int(code) >> (i * _NBR_ELEM_BITS)) & _NBR_ELEM_MAX))
    return tuple(out)

def cg_annot_pkl_path(cg_match_dict_pkl):
    base = cg_match_dict_pkl
    if base.endswith('.pkl'):
        base = base[:-len('.pkl')]
    return f'{base}.annot.pkl'

class BucketSchemaMismatch(Exception):
    pass

BUCKET_SCHEMA_VERSION = 4
ANNOTATION_SCHEMA = (
    "cg_heavy_degree/cg_num_h/cg_formal_charge int8 per CG atom (-1 unreadable, "
    "degree 0 invalid); cg_nbr_elems uint32 packed heavy-neighbour element "
    "multiset outside the match (decode_nbr_elems); "
    "perception int8 (0 ccd, 1 openbabel, 2 smiles, 3 atom-name table); "
    "vdm_buried_area / vdm_shared_area float32, vdm_n_atom_pairs int16, "
    "vdm_min_heavy_dist float32 per vdM slot; charge_sign U11 scalar, one of "
    "CHARGE_SIGNS: matches directly whichever partition "
    "the bucket's own nr_vdgs/<size>/<sign>/ directory names. cg_placed_h int8 "
    "per CG atom: 1 has a placed H within the bond-length "
    "cutoff in the parent structure, 0 none -- geometric, from whatever parent "
    "dir the library resolves to (not cg_num_h, the CCD template's nominal "
    "count).")

CHARGE_SIGNS = ('pos', 'neut', 'neg', 'unreadable')

def bucket_schema_version(data):
    if "schema" not in getattr(data, "files", []):
        return 1
    try:
        return int(json.loads(str(data["schema"]))["schema_version"])
    except Exception:
        return 1

def load_vdg_bucket(vdg_lib_dir, frag_name, subset_size, sign, aa_bucket):
    if (frag_name, vdg_lib_dir) not in _warned_incomplete_frags and \
            not _vdg_job_completed(frag_name, vdg_lib_dir):
        _warned_incomplete_frags.add((frag_name, vdg_lib_dir))
        print(f"[WARNING] Fragment {frag_name!r} in {vdg_lib_dir} has no completed "
              f"vdG-generation job (no 'Job completed.' in its log); loading anyway.")

    npz_path = vdg_npz_path(vdg_lib_dir, frag_name, subset_size, sign, aa_bucket)
    if not os.path.isfile(npz_path):
        legacy_path = os.path.join(vdg_lib_dir, frag_name, "nr_vdgs",
                                    str(subset_size), f"{aa_bucket}.npz")
        if os.path.isfile(legacy_path):
            raise BucketSchemaMismatch(f"{npz_path} does not exist. ")
        return None

    try:
        with np.load(npz_path) as data:
            version = bucket_schema_version(data)
            if version != BUCKET_SCHEMA_VERSION:
                raise BucketSchemaMismatch(
                    f"{npz_path} was written with bucket schema version "
                    f"{version}, but this code reads version "
                    f"{BUCKET_SCHEMA_VERSION}; rebuild the fragment.")
            return dict(
                aa_bucket_parts=[str(x) for x in data["aa_bucket_parts"]],
                charge_sign=str(data["charge_sign"]),
                cluster_id=data["cluster_id"].astype(np.int32),
                cluster_num_parents=data["cluster_num_parents"].astype(np.int32),
                cg=data["nr_cg_coords"].astype(np.float32),
                bb=data["nr_vdm_bb_coords"].astype(np.float32),
                resnames=data["nr_scrr_resname"],
                slot_flags=data["nr_slot_flag"].astype(np.int8),
            )
    except CORRUPT_NPZ_ERRORS as e:
        print(f"[WARNING] Could not load {npz_path}: {type(e).__name__}: {e}")
        return None

def cluster_member_indices(data, cluster_id):
    return np.nonzero(data["mem_cluster_id"] == cluster_id)[0]

def load_cluster_members(data, cluster_id, pdb_dir=None, indices=None):
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

CG_SYMMETRY_FILENAME = "cg_symmetry.npz"

def cg_symmetry_path(vdg_lib_dir, frag_name=None):
    parts = [vdg_lib_dir] if frag_name is None else [vdg_lib_dir, frag_name]
    return os.path.join(*parts, "nr_vdgs", CG_SYMMETRY_FILENAME)

def write_cg_symmetry(vdglib_dir, cg_smarts, cg_automorphisms):
    perms = np.asarray([list(p) for p in cg_automorphisms], dtype=np.int32)
    if perms.ndim != 2:
        raise ValueError(f"Expected a rectangular permutation table, got {perms.shape}")
    path = cg_symmetry_path(vdglib_dir)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.{os.getpid()}.tmp.npz"
    np.savez_compressed(tmp,
                        cg_smarts=np.asarray(str(cg_smarts)),
                        cg_automorphisms=perms)
    os.replace(tmp, path)

FRAGMENT_ALIASES_FILENAME = "fragment_aliases.tsv"

def load_fragment_aliases(vdg_lib_dir):
    path = os.path.join(vdg_lib_dir, FRAGMENT_ALIASES_FILENAME)
    aliases = {}
    if not os.path.isfile(path):
        return aliases
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1]:
                aliases[fields[0]] = fields[1]
    return aliases

def resolve_fragment_alias(vdg_lib_dir, smiles):
    return load_fragment_aliases(vdg_lib_dir).get(smiles, smiles)

def load_cg_symmetry(vdg_lib_dir, frag_name):
    path = cg_symmetry_path(vdg_lib_dir, frag_name)
    try:
        with np.load(path) as data:
            smarts = str(data["cg_smarts"])
            raw_perms = [tuple(int(i) for i in row)
                         for row in data["cg_automorphisms"]]
    except FileNotFoundError:
        raise
    except CORRUPT_NPZ_ERRORS as e:
        raise ValueError(
            f"CG symmetry sidecar {path} is unreadable "
            f"({type(e).__name__}: {e}). Falling back to the identity group "
            "would superpose a symmetric CG under the wrong atom "
            "correspondence, so this is fatal; rebuild the fragment.") from e
    return smarts, utils.validate_atom_permutations(raw_perms)

def vdg_npz_path(vdg_lib_dir, frag_name, subset_size, sign, aa_bucket):
    return os.path.join(vdg_lib_dir, frag_name, "nr_vdgs", str(subset_size), sign,
        f"{aa_bucket}.npz")

def apply_rigid_transform(ag, R, t):
    coords = ag.getCoords()
    ag.setCoords(coords @ np.asarray(R, float) + np.asarray(t, float))
    return ag

def aa_perm_indices(bucket_parts):
    return utils.group_preserving_permutations(bucket_parts)

def run_cli(main):
    try:
        main()
    except ParentPdbDirError as err:
        raise SystemExit(f"[ERROR] {err}")
