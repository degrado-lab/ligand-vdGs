import os
import numpy as np
import prody as pr
from itertools import permutations, product, combinations
from numba import njit
import hashlib

def get_vdg_AA_permutations(reordered_AAs, _vdgs):
    permuted_indices = permute_AA_duplicates(reordered_AAs)
    all_AA_cg_perm_cg_coords = []
    all_AA_cg_perm_vdm_bbcoords = []
    all_AA_cg_perm_flankingseqs = []
    all_AA_cg_perm_flankingCAs = []
    all_AA_cg_perm_pdbpaths = []
    all_AA_cg_perm_vdm_scrr = []

    # Iterate over all AA permutations of each vdg
    for _vdg in _vdgs:
        # CG coords and pdbpaths remain unchanged, but vdmbbs, flankingseqs, 
        # flankingCAs, and scrrs need to be permuted.
        for permutation in permuted_indices:
            all_AA_cg_perm_cg_coords.append(_vdg[0])
            all_AA_cg_perm_pdbpaths.append(_vdg[4])
            nonpermuted_vdmbb = _vdg[1]
            nonpermuted_flankingseqs = _vdg[2]
            nonpermuted_flankingCAs = _vdg[3]
            nonpermuted_vdm_scrr = _vdg[5]
            vdmbb_permutation = [nonpermuted_vdmbb[ix] for ix in permutation]
            flankingseqs_permutation = [nonpermuted_flankingseqs[ix] for ix in 
                                        permutation]
            flankingCAs_permutation = [nonpermuted_flankingCAs[ix] for ix in permutation]
            vdm_scrrs_permutation = [nonpermuted_vdm_scrr[ix] for ix in permutation]
            all_AA_cg_perm_vdm_bbcoords.append(vdmbb_permutation)
            all_AA_cg_perm_flankingseqs.append(flankingseqs_permutation)
            all_AA_cg_perm_flankingCAs.append(flankingCAs_permutation)
            all_AA_cg_perm_vdm_scrr.append(vdm_scrrs_permutation)

    return all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords, \
        all_AA_cg_perm_flankingseqs, all_AA_cg_perm_flankingCAs, \
        all_AA_cg_perm_pdbpaths, all_AA_cg_perm_vdm_scrr

def permute_AA_duplicates(seq):
    # Dictionary to store indices for each element in the sequence
    seen = {}
    for i, item in enumerate(seq):
        if item not in seen:
            seen[item] = []
        seen[item].append(i)

    permute_groups = [] # list of all index groups that have duplicates
    for indices in seen.values():
        if len(indices) > 1:  # only interested in elements with duplicates
            permute_groups.append(indices)

    if not permute_groups: # no duplicates, so return the original index list
        return [list(range(len(seq)))]

    # generate all possible permutations of indices within each group of duplicates
    permuted_idx_lists = []
    for perm_combination in product(*[permutations(group) for group in 
                                      permute_groups]):
        # start with the list of original indices
        permuted_idx = list(range(len(seq)))

        # flatten the product of permutations and assign them to the corresponding 
        # positions
        for group_idx, perm in zip(permute_groups, perm_combination):
            for orig_idx, new_idx in zip(group_idx, perm):
                permuted_idx[orig_idx] = new_idx

        permuted_idx_lists.append(permuted_idx)
   
    return permuted_idx_lists

def combine_cg_and_vdmbb_coords(all_AA_cg_perm_cg_coords, 
                                all_AA_cg_perm_vdm_bbcoords):
    all_cg_and_vdmbb_coords = []
    for _cg, _vdmbb_per_res in zip(all_AA_cg_perm_cg_coords, 
                                   all_AA_cg_perm_vdm_bbcoords):
        flattened_cg_coords = np.array(flatten_cg_coords(_cg))
        flattened_vdmbbs =    np.array(flatten_vdg_bbs(_vdmbb_per_res))
        cg_and_vdmbb = np.vstack([flattened_cg_coords, flattened_vdmbbs])
        all_cg_and_vdmbb_coords.append(cg_and_vdmbb)
    return all_cg_and_vdmbb_coords

def flatten_cg_coords(_cg):
    coords = []
    for atom in _cg:
       coords.append(atom)
    return coords

def flatten_vdg_bbs(_vdmbb):
    # Flatten backbones to a single list where each element is a numpy array of the 
    # x, y, z coords of each vdm backbone atom
    bbs = []
    for res in _vdmbb:
       for atom in res:
          bbs.append(atom)
    return bbs

def flatten_flanking_CAs(cgvdmbb_clus_flankingCAs):
    cgvdmbb_clus_flat_flankCAs = []
    for vdg_flankingCAs in cgvdmbb_clus_flankingCAs:
        # `vdg_flankingCAs` is a list of vdm residues, so flatten it
        flat_flanking_CAs = []
        for res in vdg_flankingCAs:
            for CA_coord in res:
                flat_flanking_CAs.append(CA_coord)
        cgvdmbb_clus_flat_flankCAs.append(flat_flanking_CAs)
    return cgvdmbb_clus_flat_flankCAs

def flatten_flanking_seqs(flankingCAs_clus_flankingseqs):
    flattened_flankingseqs_for_vdgs_in_flankingCA_clus = []
    for vdg_flankingseq in flankingCAs_clus_flankingseqs:
        flat_vdg_flankingseq = []
        for vdm_res in vdg_flankingseq:
            flat_vdg_flankingseq += vdm_res
        flattened_flankingseqs_for_vdgs_in_flankingCA_clus.append(
            flat_vdg_flankingseq)
    return flattened_flankingseqs_for_vdgs_in_flankingCA_clus

def kabsch(X, Y, chunk_size=30000):
    """
    Rotate and translate X into Y to minimize SSD (Kabsch, 1976).
    X, Y: arrays of shape [M, N, 3]
    Returns:
        R:   [M, 3, 3]
        t:   [M, 3]
        ssd: [M]
    """
    R_chunks, t_chunks, ssd_chunks = [], [], []
    M = len(X)

    for start in range(0, M, chunk_size):
        stop = min(start + chunk_size, M)
        Xc_in = X[start:stop]
        Yc_in = Y[start:stop]

        mask = np.logical_or(np.isnan(Xc_in), np.isnan(Yc_in))
        valid_atom = ~np.any(mask, axis=2)
        N_atoms = np.sum(valid_atom, axis=1, keepdims=True)
        zero_rows = (N_atoms == 0).squeeze(-1)
        safeN = np.where(N_atoms == 0, 1, N_atoms)

        X_nonan = np.where(mask, 0.0, Xc_in).astype(np.float32)
        Y_nonan = np.where(mask, 0.0, Yc_in).astype(np.float32)

        valid_coords = np.repeat(valid_atom[:, :, None], 3, axis=2)
        Xbar = (np.sum(X_nonan * valid_coords, axis=1, keepdims=True) /
                safeN[:, None, :].astype(np.float32))
        Ybar = (np.sum(Y_nonan * valid_coords, axis=1, keepdims=True) /
                safeN[:, None, :].astype(np.float32))

        Xc = X_nonan - Xbar
        Yc = Y_nonan - Ybar
        Xc[mask] = 0.0
        Yc[mask] = 0.0

        H = np.matmul(np.transpose(Xc, (0, 2, 1)), Yc)
        U, S, Vt = np.linalg.svd(H, full_matrices=False)
        d = np.sign(np.linalg.det(np.matmul(U, Vt))).astype(np.float32)

        D = np.zeros((H.shape[0], 3, 3), dtype=np.float32)
        D[:, 0, 0] = 1.0
        D[:, 1, 1] = 1.0
        D[:, 2, 2] = d

        R = np.matmul(U.astype(np.float32), np.matmul(D, Vt.astype(np.float32)))
        t = (Ybar - np.matmul(Xbar, R)).reshape(-1, 3)

        XRmY = np.matmul(Xc, R) - Yc
        ssd = np.sum(XRmY ** 2, axis=(1, 2)).astype(np.float64)

        if np.any(zero_rows):
            ssd[zero_rows] = np.inf

        R_chunks.append(R)
        t_chunks.append(t.astype(np.float32))
        ssd_chunks.append(ssd)

    R = np.concatenate(R_chunks) if R_chunks else np.empty((0, 3, 3), np.float32)
    t = np.concatenate(t_chunks) if t_chunks else np.empty((0, 3), np.float32)
    ssd = np.concatenate(ssd_chunks) if ssd_chunks else np.empty((0,), np.float64)
    return R, t, ssd

def calc_seq_similarity(list1, list2):
    # Percent identity of flanking residue names after dropping positions labeled 'vdm' 
    # or 'X'. Returns identity=0 (max dissimilarity) if all positions drop.
    list1 = [i for i in list1 if i != 'vdm']
    list2 = [i for i in list2 if i != 'vdm']
    assert len(list1) == len(list2)
    # Exclude residues if at least one of them is an "X"
    indices_to_exclude = []
    for ind in range(len(list1)):
        if list1[ind] == 'X' or list2[ind] == 'X':
            indices_to_exclude.append(ind)
    list1 = [item for idx, item in enumerate(list1) if idx not in indices_to_exclude]
    list2 = [item for idx, item in enumerate(list2) if idx not in indices_to_exclude]
    matches = sum(1 for a, b in zip(list1, list2) if a == b)

    # Calculate the percentage of matches. If all residues were dropped out, 
    # it should be treated as max dissimilarity.
    if len(list1) == 0:
        match_percentage = 0
    else:
        match_percentage = (matches / len(list1)) * 100
    return match_percentage

def get_res_iden(vdm_obj):
    seg =     list(set(vdm_obj.getSegnames()))
    chain =   list(set(vdm_obj.getChids()))
    resnum =  list(set(vdm_obj.getResnums()))
    resname = list(set(vdm_obj.getResnames()))
    assert len(seg) == 1
    assert len(chain) == 1
    assert len(resnum) == 1
    assert len(resname) == 1
    return seg[0], chain[0], resnum[0], resname[0]

def found_chain_break(flanking_seq_dict, chain_break_ind):
    # If a chain break is found, overwrite all the subsequent (or preceding) flanking
    # residues beyond the chain break as AA "X".
    flank_indices = list(flanking_seq_dict.keys())
    if chain_break_ind > 0 :
        inds_to_overwrite = [n for n in flank_indices if n >= chain_break_ind]
    elif chain_break_ind < 0:
        inds_to_overwrite = [n for n in flank_indices if n <= chain_break_ind]
   
    for overwrite_ind in inds_to_overwrite:
        flanking_seq_dict[overwrite_ind] = \
            ['X', np.array([np.nan, np.nan, np.nan])] # ['X', None]
    return flanking_seq_dict

def _pick_single_ca(ca_sel, prev_ca):
    # Return a single CA atom object from a ProDy selection, handling altlocs/dups
    if len(ca_sel) == 1:
        return ca_sel[0]

    # If we have multiple CAs (altlocs/dups), prefer the one closest to prev_ca.
    if prev_ca is None:
        return ca_sel[0]

    best_atom, best_dist = None, float('inf')
    for atom in ca_sel:
        d = pr.calcDistance(prev_ca, atom)
        if d < best_dist:
            best_atom, best_dist = atom, d
    return best_atom

def get_AA_and_CA_coords(prody_obj, current_resindex):
    '''
    Return (AA, CA_coords) for a residue index.
    If residue is missing/non-protein, returns AA='X' and CA_coords = [nan, nan, nan].
    '''
    # Build selection for this residue index
    if current_resindex < 0: 
        sel_str = f'resindex `{current_resindex}`'
    else:
        sel_str = f'resindex {current_resindex}'

    curr_resindex_obj = prody_obj.select(sel_str)
    if curr_resindex_obj is None or curr_resindex_obj.protein is None:
        return 'X', np.array([np.nan, np.nan, np.nan])

    # Resolve residue name
    curr_res_AA = list(set(curr_resindex_obj.getResnames()))
    if len(curr_res_AA) != 1:
        print(f'[WARNING] get_AA_and_CA_coords: \nResindex {current_resindex} in '
              f'{prody_obj.getTitle()} contains >1 AA: {curr_res_AA}. Defaulting to AA="X".')
        return 'X', np.array([np.nan, np.nan, np.nan])
    AA = curr_res_AA[0]

    # Select CA(s)
    CA_sel = curr_resindex_obj.select(sel_str + ' and name CA')
    if CA_sel is None or len(CA_sel) == 0:
        return 'X', np.array([np.nan, np.nan, np.nan])

    try:
        ca_atom = _pick_single_ca(CA_sel, prev_ca=None)
        CA_coords = ca_atom.getCoords()
        # Ensure shape is (3,)
        CA_coords = np.asarray(CA_coords, dtype=float).reshape(3,)
    except Exception:
        print(f'[WARNING] Could not resolve CA for resindex {current_resindex} in '
              f'{prody_obj.getTitle()}. Defaulting to AA="X".')
        return 'X', np.array([np.nan, np.nan, np.nan])

    return AA, CA_coords

def get_cg_coords(prody_obj, pdbpath):
    # The CG atoms have their occupancies set to >= 3.0, with unique values (e.g., 
    # 3.0, 3.1, 3.2, etc.) to allow a 1:1 correspondence of equivalent atoms between 
    # different ligands.
    cg = prody_obj.select('occupancy > 2.9')

    if cg is None or len(cg) == 0:
        print(f'[WARNING] get_cg_coords: no atoms with occupancy > 2.9 in {pdbpath}')
        return None

    num_atoms = len(cg)
    cg_coords = []
    for ind in range(num_atoms):
        occ = f'3.{ind}'
        atom = cg.select(f'occupancy == {occ}')
        if atom is None or len(atom) != 1:
            print(f'[WARNING] get_cg_coords: {0 if atom is None else len(atom)} atoms '
                  f'are occupancy {occ} in {pdbpath}.')
            return None
        atom_coords = atom.getCoords()[0]
        cg_coords.append(atom_coords)

    cg_coords = np.asarray(cg_coords, dtype=float)

    if  not np.isfinite(cg_coords).all(): # reject NaN or inf coords
            print(f'[WARNING] get_cg_coords: NaN or inf CG coords in {pdbpath}.')
            return None

    return cg_coords

def determine_flank_resnums(par, list_resnums, s, c, vdm_ca):
    # Determine which flanking residues to add to pr obj based on chain breaks. 
    resnums_to_include = []
    assert len(vdm_ca) == 1
    vdm_ca = vdm_ca[0]
    prev_ca = vdm_ca
    # Iterate over resnums 
    for r in list_resnums:
        res_sel = f'resnum `{r}`' if r < 0 else f'resnum {r}'
        if s == '':
            ca_sel = par.select(f'chain {c} and {res_sel} and name CA')
        else:
            ca_sel = par.select(f'segment {s} and chain {c} and {res_sel} and name CA')

        # If missing, stop (chain break).
        if ca_sel is None or len(ca_sel) == 0:
            return resnums_to_include

        # Collapse to one CA robustly
        try:
            curr_ca = _pick_single_ca(ca_sel, prev_ca)
        except Exception:
            return resnums_to_include

        # Check chain break distance to previous CA
        dist = pr.calcDistance(prev_ca, curr_ca)
        if dist > 4.5:  # chain break
            return resnums_to_include

        resnums_to_include.append(r)
        prev_ca = curr_ca

    return resnums_to_include

def get_bb_coords(obj):
    bb_coords = []
    for atom_name in ['N', 'CA', 'C']:
        atom_obj = obj.select(f'name {atom_name}')
        if atom_obj is None or len(atom_obj) == 0:
            return None # missing backbone atom; give up on this residue

        first_coord = atom_obj.getCoords()[0]

        if len(atom_obj) > 1: # if there are duplicates, compare them to the first one
            for i, a in enumerate(atom_obj):
                if i == 0:
                    continue  # skip the first one (reference)
                dist = pr.calcDistance(first_coord, a.getCoords())
                if dist > 0.2:
                    print(
                        f"[WARNING] Ambiguous coordinates for {obj.getTitle()} resnum "
                        f"{set(obj.getResnums())} atom {atom_name} are {dist} apart; "
                        f"using first occurrence, but investigate why.")
                    break

        # Scenarios: clean structure w/ 1 atom, duplicate atoms < 0.2A apart are 
        # likely benign, and duplicate atoms > 0.2A possibly harmful, so log warning. 
        bb_coords.append(first_coord)

    return bb_coords

def get_vdg_subsets(input_list, target_size):
    # Generate combos of residues containing `target_size` elements.
    if len(input_list) < target_size:
        return []  
    return list(combinations(input_list, target_size)) 

def get_clus_pdb_outpath(clusnum, clusnum_dir, pdbpath, scrr_cg_perm, clusmem_ind, 
                         is_centroid):
    pdbname = os.path.basename(pdbpath).removesuffix('.pdb')
    scrrs, perm = scrr_cg_perm
    perm = perm.removeprefix('cg_perm_')
    vdg_scrr_str = '_'.join(['_'.join([str(z) for z in v_s]) for v_s in scrrs])
    if is_centroid:
        pdb_str = f'clus{clusnum}_{pdbname}_{vdg_scrr_str}_CGperm{perm}_ix{clusmem_ind}_centroid.pdb.gz'
    else:
        pdb_str = f'clus{clusnum}_{pdbname}_{vdg_scrr_str}_CGperm{perm}_ix{clusmem_ind}.pdb.gz'

    return os.path.join(clusnum_dir, pdb_str)

def select_diverse_pdbIDs(strings, k): # k = max num to select
   # Greedy max–min Hamming selection over equal-length strings (PDB IDs)
   # to promote dataset diversity. Assumes PDB IDs are 4 characters long.
   ascii_array = strings_to_ascii_array(strings)
   selected_indices = select_diverse_subset_greedy(ascii_array, k)
   return [strings[i] for i in selected_indices]

def strings_to_ascii_array(strings):
   return np.array([[ord(c) for c in s] for s in strings], dtype=np.uint8)

@njit
def hamming(s1, s2):
   dist = 0
   for i in range(len(s1)):
      if s1[i] != s2[i]:
         dist += 1
   return dist

@njit
def update_min_dists(data, selected_idx, selected_mask, min_dists):
   n = data.shape[0]
   for i in range(n):
      if selected_mask[i] == 0:
         dist = hamming(data[selected_idx], data[i])
         if dist < min_dists[i]:
            min_dists[i] = dist

def select_diverse_subset_greedy(data, k):
   n = data.shape[0]
   selected = [0]  # start with first point
   selected_mask = np.zeros(n, dtype=np.uint8)
   selected_mask[0] = 1
   min_dists = np.full(n, 255, dtype=np.uint8)

   # Initial distance pass
   for i in range(1, n):
      min_dists[i] = hamming(data[0], data[i])

   for _ in range(1, min(n, k)):
      # Select max of min distances
      max_idx = -1
      max_val = -1
      for i in range(n):
         if selected_mask[i] == 0 and min_dists[i] > max_val:
            max_val = min_dists[i]
            max_idx = i

      selected.append(max_idx)
      selected_mask[max_idx] = 1

      # Parallel update of min distances
      update_min_dists(data, max_idx, selected_mask, min_dists)

   return selected

def _stream_root(vdglib_dir):
    # Root for AA composition streaming buckets.
    # Tries $TMPDIR, /scratch, /tmp, and vdglib_dir as fallback.
    user = os.environ.get("USER") or getpass.getuser() or "unknown"
    vdglib_tag = os.path.basename(os.path.abspath(vdglib_dir.rstrip(os.sep)))
    tmpdir_env = os.environ.get("TMPDIR")
    if tmpdir_env:
        root = os.path.join(tmpdir_env, f"vdg_stream_{vdglib_tag}")
    elif os.path.isdir("/scratch"):
        root = os.path.join("/scratch", user, f"vdg_stream_{vdglib_tag}")
    elif os.path.isdir("/tmp"):
        root = os.path.join("/tmp", user, f"vdg_stream_{vdglib_tag}")
    else:
        root = os.path.join(vdglib_dir, "stream_tmp")
    os.makedirs(root, exist_ok=True)
    return root

def _aa_tmp_dir(vdglib_dir, size_subset):
    # Write per-AA-bucket records to disk to free memory
    return os.path.join(_stream_root(vdglib_dir), str(size_subset))

def _aa_tmp_path(vdglib_dir, size_subset, aa_key):
    # AA bucket file name is based on AA key plus a short SHA1 suffix
    h = hashlib.sha1(aa_key.encode()).hexdigest()[:16]
    fname = f"{aa_key[:80]}__{h}.jsonl.gz"
    return os.path.join(_aa_tmp_dir(vdglib_dir, size_subset), fname)

def normalize_rmsd(num_atoms, atoms):
    ''' Return a size-normalized RMSD threshold (Å) for the given atom set
    ('cgvdmbb' or 'flankbb'). The threshold scales linearly with the number
    of atoms between 8 and 15:
        flankbb: 0.5 Å → 1.5 Å
        cgvdmbb: 0.5 Å → 1.0 Å
    Below 8 atoms, use the minimum; above 15, use the maximum.
    '''

    if atoms == 'flankbb':
        max_threshold = 1.5
        min_threshold = 0.5
    elif atoms == 'cgvdmbb':
        max_threshold = 1.0
        min_threshold = 0.5
    else:
        raise ValueError(f"Unknown atom set for normalize_rmsd: {atoms}")

    min_atoms = 8
    max_atoms = 15

    if num_atoms < min_atoms:
        return min_threshold
    if num_atoms > max_atoms:
        return max_threshold

    scaling_factor = (num_atoms - min_atoms) / (max_atoms - min_atoms)
    threshold = min_threshold + scaling_factor * (max_threshold - min_threshold)
    return threshold