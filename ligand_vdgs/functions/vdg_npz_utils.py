# vdg_npz_utils.py

import os
import numpy as np
import prody as pr
import utils

def _append_full_ligand_from_parent(all_coords, names, resnames, resnums, chids, segnames, 
    elements, occupancies, parent_pdb_path, cg_coords, cg_chain, cg_resnum, cg_names):
    "Add non-CG ligand atoms from parent PDB, mapped into cg_coords frame, occ=4.0."
    if not parent_pdb_path:
        return
    try:
        pdb = pr.parsePDB(parent_pdb_path)
    except Exception as e:
        print(f"[WARNING] parsePDB failed for {parent_pdb_path}: {e}")
        return

    cg_coords = np.asarray(cg_coords, float)
    cg_chain  = np.asarray(cg_chain)
    cg_resnum = np.asarray(cg_resnum, int)
    cg_names  = np.asarray(cg_names)

    # Map CG atoms in the *current* frame: (chain,resnum,name) -> coord
    cg_atom_coords = {}
    for coord, ch, rn, name in zip(cg_coords, cg_chain, cg_resnum, cg_names):
        key = (str(ch), int(rn), str(name))
        cg_atom_coords[key] = np.asarray(coord, float)

    cg_keys     = set(cg_atom_coords.keys())
    cg_res_keys = sorted({(str(c), int(r)) for c, r in zip(cg_chain, cg_resnum)})

    for ch_res in cg_res_keys:
        ch, rn = ch_res
        sel = pdb.select(f"chain {ch} and resnum {rn}")
        if sel is None:
            continue

        coords_pdb = sel.getCoords()
        anames     = sel.getNames()
        rnames     = sel.getResnames()
        rnums      = sel.getResnums()
        chs        = sel.getChids()
        elems      = sel.getElements()
        segs       = sel.getSegnames() if hasattr(sel, "getSegnames") else np.array([""] * len(sel))

        # Build matched pairs: parent-PDB CG atoms ↔ current CG coords
        match_pdb, match_cg = [], []
        for c_pdb, an, rn_, ch_ in zip(coords_pdb, anames, rnums, chs):
            key = (str(ch_), int(rn_), str(an))
            if key in cg_atom_coords:
                match_pdb.append(np.asarray(c_pdb, float))
                match_cg.append(cg_atom_coords[key])

        if len(match_pdb) < 3:
            coords_mapped = coords_pdb  # too few shared atoms to fit; use as-is
        else:
            X = np.asarray(match_pdb, float).reshape(1, -1, 3)  # mobile: parent PDB
            Y = np.asarray(match_cg, float).reshape(1, -1, 3)  # target: cg_coords frame
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
            occupancies.append(4.0)

def build_vdg_atomgroup_from_npz(cg_coords, cg_names, cg_elements, cg_seg, cg_chain, 
    cg_resnum, cg_resname, vdm_bb_coords, vdm_atom_coords, vdm_atom_names, 
    vdm_atom_elements, scrr_seg, scrr_chain, scrr_resnum, scrr_resname,
    include_full_ligand=False, parent_pdb_path=None,):
    "Materialize a vdG AtomGroup from centroid arrays (CG + vdG; optional full ligand)."
    cg_coords, vdm_bb_coords = np.asarray(cg_coords, float), np.asarray(vdm_bb_coords, float)
    n_cg, num_vdms = cg_coords.shape[0], vdm_bb_coords.shape[0]
    vdm_bb_flat = vdm_bb_coords.reshape(-1, 3)

    all_coords, names, resnames, resnums = [], [], [], []
    chids, segnames, elements, occupancies = [], [], [], []

    cg_names, cg_elements = np.asarray(cg_names), np.asarray(cg_elements)
    cg_seg, cg_chain = np.asarray(cg_seg), np.asarray(cg_chain)
    cg_resnum, cg_resname = np.asarray(cg_resnum), np.asarray(cg_resname)

    # CG: occupancy 3.x
    for i in range(n_cg):
        all_coords.append(cg_coords[i])
        names.append(str(cg_names[i]))
        resnames.append(str(cg_resname[i]))
        resnums.append(int(cg_resnum[i]))
        chids.append(str(cg_chain[i]))
        seg_i = cg_seg[i]
        segnames.append("" if seg_i in (None, "", "None") else str(seg_i))
        elem_i = cg_elements[i]
        elements.append(str(elem_i) if elem_i not in (None, "", " ") else "C")
        occupancies.append(3.0 + 0.1 * i)

    # Optional: ligand atoms (non-CG) from parent PDB, mapped to same frame, occ=4.0
    if include_full_ligand and parent_pdb_path is not None:
        _append_full_ligand_from_parent(all_coords, names, resnames, resnums, chids, 
            segnames, elements, occupancies, parent_pdb_path, cg_coords, cg_chain, 
            cg_resnum, cg_names,)

    # vdM backbones + side chains: occupancy 2.0
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
            occupancies.append(2.0)

        coords_i = np.asarray(vdm_atom_coords[v_idx], float)
        names_i  = np.asarray(vdm_atom_names[v_idx])
        elems_i  = np.asarray(vdm_atom_elements[v_idx])
        if coords_i.size == 0:
            continue
        if not (coords_i.shape[0] == names_i.shape[0] == elems_i.shape[0]):
            raise ValueError(
                f"vdM atom array mismatch at {v_idx}: "
                f"coords={coords_i.shape[0]}, names={names_i.shape[0]}, elems={elems_i.shape[0]}")
        for a_coord, a_name, a_elem in zip(coords_i, names_i, elems_i):
            a_name_str = str(a_name)
            if a_name_str in ("N", "CA", "C"):
                continue
            all_coords.append(a_coord)
            names.append(a_name_str)
            resnames.append(resname)
            resnums.append(resnum)
            chids.append(ch)
            segnames.append(seg)
            elem_str = str(a_elem) if a_elem not in (None, "", " ") else "C"
            elements.append(elem_str)
            occupancies.append(2.0)

    all_coords = np.asarray(all_coords, float)
    ag = pr.AtomGroup("vdg_centroid")
    ag.setCoords(all_coords)
    ag.setNames(np.asarray(names, "U4"))
    ag.setResnames(np.asarray(resnames, "U4"))
    ag.setResnums(np.asarray(resnums, int))
    ag.setChids(np.asarray(chids, "U2"))
    ag.setSegnames(np.asarray(segnames, "U8"))
    ag.setElements(np.asarray(elements, "U2"))
    ag.setOccupancies(np.asarray(occupancies, float))

    cg_vdmbb_flat = np.vstack([cg_coords, vdm_bb_flat])
    return ag, cg_vdmbb_flat

def make_aa_bucket(aa_list, sort=True):
    """
    Turn a list/tuple of AA identities (incl. 'bb') into the canonical bucket name.

    Example:
        ['bb', 'ALA']  -> 'ALA_bb'  (sorted)
        ['GLU', 'GLU'] -> 'GLU_GLU'
    """
    aa_list = list(aa_list)
    if sort:
        aa_list = sorted(aa_list)
    return "_".join(aa_list)

def load_vdg_bucket(vdg_lib_dir, frag_name, subset_size, aa_bucket):
    """
    Load one vdG npz bucket and return a dict:

        {
          "cluster_id": [C],
          "cg":        [C, n_cg, 3]          # float32
          "bb":        [C, num_vdms, 3, 3]   # float32, N/CA/C per vdM
          "resnames":  [C, num_vdms],        # str, AA names
        }

    Returns None if the file is missing or cannot be read.
    """
    npz_path = vdg_npz_path(vdg_lib_dir, frag_name, subset_size, aa_bucket)
    if not os.path.isfile(npz_path):
        return None

    try:
        data = np.load(npz_path, allow_pickle=True)
        bucket = dict(
            cluster_id=data["cluster_id"].astype(int),
            cg=data["centroid_cg_coords"].astype(np.float32),
            bb=data["centroid_vdm_bb_coords"].astype(np.float32),
            resnames=data["centroid_scrr_resname"],  # dtype <U4/string
        )
    except Exception as e:
        print(f"[WARNING] Could not load {npz_path}: {e}")
        bucket = None

    return bucket

def vdg_npz_path(vdg_lib_dir, frag_name, subset_size, aa_bucket):
    return os.path.join(vdg_lib_dir, frag_name, "nr_vdgs", str(subset_size),
        f"{aa_bucket}.npz")

def iter_bucket_centroids(bucket):
    """
    Yield one vdG centroid at a time as a dict.

    Each item:
        { "index":       int,
          "cluster_id":  int,
          "cg":          [n_cg, 3] float32,
          "bb":          [num_vdms, 3, 3] float32,
          "resnames":    [num_vdms] (str)  }
    """
    if bucket is None:
        return
    cluster_ids = bucket["cluster_id"]
    cg         = bucket["cg"]
    bb         = bucket["bb"]
    resnames   = bucket["resnames"]
    n = cg.shape[0]
    for i in range(n):
        yield dict(index=int(i), cluster_id=int(cluster_ids[i]), cg=cg[i], bb=bb[i],
            resnames=list(resnames[i]))

def permute_backbone_flat(bb_coords, resind_perm):
    """
    Given:
        bb_coords   : [num_vdms, 3, 3] (N, CA, C per vdM)
        resind_perm : iterable of vdM indices into bb_coords

    Return:
        [subset_size * 3, 3] float32
        (N, CA, C for each vdM in permuted order, then flattened over residues)
    """
    bb_list = []
    for ridx in resind_perm:
        bb_list.extend(bb_coords[ridx])  # each is [3,3], appended in order N,CA,C
    return np.asarray(bb_list, dtype=np.float32)
