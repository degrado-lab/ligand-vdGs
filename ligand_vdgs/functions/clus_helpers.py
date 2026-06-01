# clus_helpers.py

import os
import numpy as np
import prody as pr
from itertools import permutations, product, combinations
from numba import njit
import hashlib
import getpass

def get_vdg_AA_permutations(reordered_AAs, _vdgs):
    """
    Expand vdGs over AA permutations (for duplicated AAs).

    Each _vdg is:
        [cg_coords, bbcoords, flankseqs, flankCAs, pdbpath, scrr, vdm_heavycoords,
         cg_names, cg_elements, cg_seg, cg_chain, cg_resnum, cg_resname]
    """
    permuted_indices = permute_AA_duplicates(reordered_AAs)
    all_AA_cg_perm_cg_coords = []
    all_AA_cg_perm_vdm_bbcoords = []
    all_AA_cg_perm_flankingseqs = []
    all_AA_cg_perm_flankingCAs = []
    all_AA_cg_perm_pdbpaths = []
    all_AA_cg_perm_vdm_scrr = []
    all_AA_cg_perm_vdm_heavycoords = []
    all_AA_cg_perm_cg_names = []
    all_AA_cg_perm_cg_elements = []
    all_AA_cg_perm_cg_seg = []
    all_AA_cg_perm_cg_chain = []
    all_AA_cg_perm_cg_resnum = []
    all_AA_cg_perm_cg_resname = []

    # Iterate over all AA permutations of each vdg
    for _vdg in _vdgs:
        # CG coords and pdbpaths remain unchanged, but vdmbbs, flankingseqs, 
        # flankingCAs, scrrs, and vdm_heavycoords need to be permuted.
        for permutation in permuted_indices:
            all_AA_cg_perm_cg_coords.append(_vdg[0])
            all_AA_cg_perm_pdbpaths.append(_vdg[4])

            nonpermuted_vdmbb = _vdg[1]
            nonpermuted_flankingseqs = _vdg[2]
            nonpermuted_flankingCAs = _vdg[3]
            nonpermuted_vdm_scrr = _vdg[5]
            nonpermuted_vdm_heavy = _vdg[6]

            vdmbb_permutation = [nonpermuted_vdmbb[ix] for ix in permutation]
            flankingseqs_permutation = [nonpermuted_flankingseqs[ix] for ix in permutation]
            flankingCAs_permutation = [nonpermuted_flankingCAs[ix] for ix in permutation]
            vdm_scrrs_permutation = [nonpermuted_vdm_scrr[ix] for ix in permutation]
            vdm_heavy_permutation = [nonpermuted_vdm_heavy[ix] for ix in permutation]

            all_AA_cg_perm_vdm_bbcoords.append(vdmbb_permutation)
            all_AA_cg_perm_flankingseqs.append(flankingseqs_permutation)
            all_AA_cg_perm_flankingCAs.append(flankingCAs_permutation)
            all_AA_cg_perm_vdm_scrr.append(vdm_scrrs_permutation)
            all_AA_cg_perm_vdm_heavycoords.append(vdm_heavy_permutation)
            all_AA_cg_perm_cg_names.append(_vdg[7])
            all_AA_cg_perm_cg_elements.append(_vdg[8])
            all_AA_cg_perm_cg_seg.append(_vdg[9])
            all_AA_cg_perm_cg_chain.append(_vdg[10])
            all_AA_cg_perm_cg_resnum.append(_vdg[11])
            all_AA_cg_perm_cg_resname.append(_vdg[12])

    return (all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords, 
        all_AA_cg_perm_flankingseqs, all_AA_cg_perm_flankingCAs, all_AA_cg_perm_pdbpaths,
        all_AA_cg_perm_vdm_scrr, all_AA_cg_perm_vdm_heavycoords, all_AA_cg_perm_cg_names,
        all_AA_cg_perm_cg_elements, all_AA_cg_perm_cg_seg, all_AA_cg_perm_cg_chain,
        all_AA_cg_perm_cg_resnum, all_AA_cg_perm_cg_resname)

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

def combine_cg_and_vdmbb_coords(all_cg, all_vdmbb):
    out = []
    for _cg, _vdmbb in zip(all_cg, all_vdmbb):
        flattened_cg = np.asarray(_cg, dtype=np.float32)
        flattened_vdmbb = np.asarray([atom for res in _vdmbb for atom in res],
            dtype=np.float32,)
        out.append(np.vstack([flattened_cg, flattened_vdmbb]))
    return out

def flatten_flanking_CAs(cgvdmbb_clus_flankingCAs):
    cgvdmbb_clus_flat_flankCAs = []
    for vdg_flankingCAs in cgvdmbb_clus_flankingCAs:
        # `vdg_flankingCAs` is a list of vdm residues, so flatten it
        flat_flanking_CAs = []
        for res in vdg_flankingCAs:
            for CA_coord in res:
                flat_flanking_CAs.append(np.asarray(CA_coord, dtype=np.float32))
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
        return 'X', np.array([np.nan, np.nan, np.nan], dtype=np.float32)

    # Resolve residue name
    curr_res_AA = list(set(curr_resindex_obj.getResnames()))
    if len(curr_res_AA) != 1:
        print(f'[WARNING] get_AA_and_CA_coords: \nResindex {current_resindex} in '
              f'{prody_obj.getTitle()} contains >1 AA: {curr_res_AA}. Defaulting to AA="X".')
        return 'X', np.array([np.nan, np.nan, np.nan], dtype=np.float32)
    AA = curr_res_AA[0]

    # Select CA(s)
    CA_sel = curr_resindex_obj.select(sel_str + ' and name CA')
    if CA_sel is None or len(CA_sel) == 0:
        return 'X', np.array([np.nan, np.nan, np.nan], dtype=np.float32)

    try:
        ca_atom = _pick_single_ca(CA_sel, prev_ca=None)
        CA_coords = np.asarray(ca_atom.getCoords(), dtype=np.float32)
    except Exception:
        print(f'[WARNING] Could not resolve CA for resindex {current_resindex} in '
              f'{prody_obj.getTitle()}. Defaulting to AA="X".')
        return 'X', np.array([np.nan, np.nan, np.nan], dtype=np.float32)

    return AA, CA_coords

def get_cg_atoms(prody_obj, pdbpath):
    cg = prody_obj.select('occupancy > 2.9')
    if cg is None or len(cg) == 0:
        print(f'[WARNING] get_cg_atoms: no atoms with occupancy > 2.9 in {pdbpath}')
        return None
    num_atoms = len(cg)
    coords, names, elements, segs, chains, resnums, resnames = [], [], [], [], [], [], []
    for ind in range(num_atoms):
        occ = f'3.{ind}'
        atom = cg.select(f'occupancy == {occ}')
        if atom is None or len(atom) != 1:
            print(f'[WARNING] get_cg_atoms: {0 if atom is None else len(atom)} atoms '
                  f'are occupancy {occ} in {pdbpath}.')
            return None
        a = atom[0]

        #c = np.asarray(a.getCoords(), dtype=float)
        #if c.ndim == 2:
        #    # e.g. shape (1, 3)
        #    c = c[0]
        #if c.shape != (3,):
        #    print(f'[WARNING] get_cg_atoms: unexpected coord shape {c.shape} for '
        #          f'{pdbpath}; skipping CG.')
        #    return None

        coords.append(a.getCoords())
        names.append(a.getName())
        elements.append(a.getElement())
        segs.append(a.getSegname())
        chains.append(a.getChid())
        resnums.append(int(a.getResnum()))
        resnames.append(a.getResname())
    return (np.asarray(coords, dtype=np.float32),
            names, elements, segs, chains, resnums, resnames)

def get_cg_coords(prody_obj, pdbpath):
    out = get_cg_atoms(prody_obj, pdbpath)
    if out is None:
        return None
    coords, _, _, _, _, _, _ = out
    return coords

def get_cg_atom_metadata(prody_obj, pdbpath):
    out = get_cg_atoms(prody_obj, pdbpath)
    if out is None:
        return None
    _, names, elements, segs, chains, resnums, resnames = out
    return (names, elements, segs, chains, resnums, resnames)

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
                        f"using first occurrence.")
                    break

        # Scenarios: clean structure w/ 1 atom, duplicate atoms < 0.2A apart are 
        # likely benign, and duplicate atoms > 0.2A possibly harmful, so log warning. 
        bb_coords.append(np.asarray(first_coord, dtype=np.float32))

    return bb_coords

def get_vdg_subsets_target_size(input_list, target_size):
    # Generate combos of residues containing `target_size` elements.
    if len(input_list) < target_size:
        return []  
    return list(combinations(input_list, target_size)) 

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