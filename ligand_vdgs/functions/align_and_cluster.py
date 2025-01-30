# NOTE: this script was originally cluster_nodes.py in the `production` branch
# of vdG-miner. Unable to cleanly import it from ligand-vdGs because it's
# from the vdG-miner submodule, so it's now copied over to 
# ligand-vdGs/ligand_vdgs/functions.

import os
import sys
from itertools import permutations
import numpy as np
import prody as pr
from scipy.sparse import csr_matrix
import shutil

def is_contained(u, v):
    """Determine whether there exists a mask such that np.all(u == v[mask])
    
    Parameters
    ----------
    u : np.array [M]
        The array for which containment in v is to be checked.
    v : np.array [N]
        The array for which containment of u is to be checked. Requires N > M.

    Returns
    -------
    contained : bool
        The boolean value of np.all(u == v[mask]); in other words, True is 
        returned if u is contained in v, and False otherwise.
    """
    unique_u, counts_u = np.unique(u, return_counts=True)
    unique_v, counts_v = np.unique(v, return_counts=True)

    # match unique elements in u to those in v
    match_indices = np.isin(unique_u, unique_v)

    # return False if any element of u is not in v
    if not np.all(match_indices):
        return False
    
    # get indices where unique_u elements appear in unique_v
    indices_in_v = np.searchsorted(unique_u, unique_v)

    # check if counts in v are greater than or equal to those in u
    return np.all(counts_u <= counts_v[indices_in_v])

def permute_with_symmetry(symmetry_classes):
    """Get all possible permutation arrays that preserve symmetry classes.
    
    Parameters
    ----------
    symmetry_classes : list
        List of integers representing the symmetry classes of the elements.
        
    Returns
    -------
    valid_permutations : list
        List of valid permutations that preserve symmetry classes.
    """
    elements = list(range(len(symmetry_classes)))
    # Get all possible permutations
    all_permutations = permutations(elements)
    # Filter permutations based on symmetry classes
    valid_permutations = []
    for perm in all_permutations:
        is_valid = True
        for i, p in enumerate(perm):
            # Compare the symmetry classes of elements in the original list
            # and the permuted list at the same positions
            if symmetry_classes[elements.index(p)] != symmetry_classes[i]:
                is_valid = False
                break
        if is_valid:
            valid_permutations.append(np.array(perm))
    return valid_permutations

def greedy(adj_mat, min_cluster_size=1):
    """Greedy clustering algorithm based on an adjacency matrix.
        
        Takes an adjacency matrix as input.
        All values of adj_mat are 1 or 0:  1 if <= to rmsd_cutoff, 
        0 if > rmsd_cutoff.

        The diagonal of adj_mat should be 0.

    Parameters
    ----------
    adj_mat : scipy.sparse.csr_matrix
        Adjacency matrix of the graph.
    min_cluster_size : int, optional
        Minimum size of a cluster, by default 2.

    Returns
    -------
    all_mems : list
        List of arrays of cluster members.
    cents : list
        List of cluster centers.
    """
    if not isinstance(adj_mat, csr_matrix):
        try:
            adj_mat = csr_matrix(adj_mat)
        except:
            print('adj_mat distance matrix must be scipy csr_matrix '
                  '(or able to convert to one)')
            return

    assert adj_mat.shape[0] == adj_mat.shape[1], \
        'Distance matrix is not square.'

    all_mems = []
    cents = []
    indices = np.arange(adj_mat.shape[0])

    try:
        while adj_mat.shape[0] > 0:

            cent = adj_mat.sum(axis=1).argmax()
            row = adj_mat.getrow(cent)
            tf = ~row.toarray().astype(bool)[0]
            mems = indices[~tf]

            if len(mems) < min_cluster_size:
                [cents.append(i) for i in indices]
                [all_mems.append(np.array([i])) for i in indices]
                break

            cents.append(indices[cent])
            all_mems.append(mems)

            indices = indices[tf]
            adj_mat = adj_mat[tf][:, tf]
    except KeyboardInterrupt:
        pass

    return all_mems, cents

def kabsch(X, Y, w=None):
    """Rotate and translate X into Y to minimize the MSD between the two.
       
       Implements the SVD method by Kabsch et al. (Acta Crystallogr. 1976, 
       A32, 922).

    Parameters
    ----------
    X : np.array [M x N x 3]
        Array of M sets of mobile coordinates (N x 3) to be transformed by a 
        proper rotation to minimize mean squared displacement (MSD) from Y.
    Y : np.array [M x N x 3]
        Array of M sets of stationary coordinates relative to which to 
        transform X.
    W : np.array [N], optional
        Vector of weights for fitting.

    Returns
    -------
    R : np.array [M x 3 x 3]
        Proper rotation matrices required to transform each set of coordinates 
        in X such that its MSD with the corresponding coordinates in Y is 
        minimized.
    t : np.array [M x 3]
        Translation matrix required to transform X such that its MSD with Y 
        is minimized.
    msd : np.array [M]
        Mean squared displacement after alignment for each pair of coordinates.
    """
    N = X.shape[1]
    if w is None:
        w = np.ones((1, N, 1)) / N
    else:
        w = w.reshape((1, -1, 1)) / w.sum()
    # compute R using the Kabsch algorithm
    Xbar, Ybar = np.sum(X * w, axis=1, keepdims=True), \
                 np.sum(Y * w, axis=1, keepdims=True)
    # subtract Xbar and Ybar, then weight the resulting matrices
    Xc, Yc = np.sqrt(w) * (X - Xbar), np.sqrt(w) * (Y - Ybar)
    H = np.matmul(np.transpose(Xc, (0, 2, 1)), Yc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.matmul(U, Vt)))
    D = np.zeros((X.shape[0], 3, 3))
    D[:, 0, 0] = 1.
    D[:, 1, 1] = 1.
    D[:, 2, 2] = d
    R = np.matmul(U, np.matmul(D, Vt))
    t = (Ybar - np.matmul(Xbar, R)).reshape((-1, 3))
    # compute MSD from aligned coordinates XR
    XRmY = np.matmul(Xc, R) - Yc
    msd = np.sum(XRmY ** 2, axis=(1, 2))
    return R, t, msd

def cluster_structures(vdglib_dir, rmsd_cutoff, idxs, all_AA_permuted_pdbpaths, 
                       all_AA_permuted_vdm_scrrs, out_clus_dir, 
                       all_AA_permuted_bbatoms_to_align, symmetry_classes=None):
    """Greedily cluster the structures based on RMSD.

    Parameters
    ----------
    idxs : list, optional
        Indices of CG atoms (ordered by smarts pattern) on which to cluster. 
    symmetry_classes : list, optional
        Integers representing the symmetry classes of the CG atoms on which clustering is to 
        be performed. If provided, should have the same length as idxs. If not provided, the 
        atoms are assumed to be symmetrically inequivalent.
    all_AA_permuted_bbatoms_to_align : list
        Either: 
        - BB atoms of vdms if aligning/clustering on just cg+vdm only. 
            (Step 1 in header of cluster_and_deduplicate.py)
        - BB atoms of vdms+flanking residues if aligning/clustering on cg+vdm+flanking residues.
            (Step 2 in header of cluster_and_deduplicate.py)
    Returns
    -------
    cluster_assignments : dict
        key = cluster num; value = indices of the list elements belong to that clus.
    """

    # Process each pdb
    all_pdbs, all_structs, coords, all_vdg_scrrs = [], [], [], []
    if symmetry_classes is None:
        symmetry_classes = list(range(len(idxs)))
    assert len(symmetry_classes) == len(idxs), \
        'Length of symmetry_classes must match length of idxs.'
    perms = permute_with_symmetry(symmetry_classes)
    first_vdm_scrr = all_AA_permuted_vdm_scrrs[0]
    num_vdms = len(first_vdm_scrr)
    res_idxs = list(range(num_vdms)) 
    N = len(idxs) + \
        len(res_idxs) * 3 # number of atoms in CG plus N, CA, C for each aa

    flat_all_AA_permuted_bb_atoms_to_clus = flatten_AA_permuted_bbatoms(
        all_AA_permuted_bbatoms_to_align, all_AA_permuted_pdbpaths, all_AA_permuted_vdm_scrrs)

    for _pdb, _vdgscrr, _flat_bb_atoms_to_clus in zip(all_AA_permuted_pdbpaths, 
                        all_AA_permuted_vdm_scrrs, flat_all_AA_permuted_bb_atoms_to_clus):
        coords, all_pdbs, all_structs, all_vdg_scrrs = get_struct_data(_pdb, idxs, perms,
                _flat_bb_atoms_to_clus, coords, all_pdbs, all_structs, all_vdg_scrrs, _vdgscrr)
    
    coords = np.array(coords).reshape((len(all_pdbs), len(perms), -1, 3))
    assert len(all_pdbs) == len(all_structs) == len(coords) == len(all_vdg_scrrs)
    coords = coords.transpose((1, 0, 2, 3))
    M = coords.shape[1]
    if M == 1:
        return
    
    all_mems, cents, triu_indices, R, t, msd = align_and_cluster(M, idxs, res_idxs, 
                                                                 N, coords, rmsd_cutoff) 
    assert len(all_pdbs) == len(all_structs) == len(all_vdg_scrrs)
    
    cluster_dirname = os.path.join(vdglib_dir, out_clus_dir)
    cluster_num = 1
    num_supposed_to_be_output = 0
    
    # Process and output PDBs of the clusters. Note that they must be output to a 'temp' 
    # dir because there may be degenerate binding sites that will be deduplicated later.
    for cluster, cent in zip(all_mems, cents):
        num_supposed_to_be_output, vdmAA_str = process_AA_redun_cluster(cent, cluster, all_pdbs, 
                all_structs, all_vdg_scrrs, triu_indices, R, t, msd, first_vdm_scrr, 
                cluster_dirname, cluster_num, num_supposed_to_be_output)
        cluster_num += 1

    # Use dict to record which indices (of all_AA_permuted_...) belong to which cluster. 
    # Note that this dict is called 'temp' because it may contain redundant vdGs based on 
    # different AA permutations within the same binding site.
    temp_cluster_assignments = {} # key = clusnum, value = indices 
    clus_num = 1
    for cluster_result in all_mems:
        temp_cluster_assignments[clus_num] = cluster_result
        clus_num += 1
    assert num_supposed_to_be_output == sum([len(v) for k, v in temp_cluster_assignments.items()])

    # Now, merge clusters based on duplicate binding sites (diff permutations of AAs).
    temp_cluster_assignments = remove_redundant_AAperm_clusters(temp_cluster_assignments, 
                                cluster_dirname, num_vdms, vdmAA_str)
    
    # Redefine cluster_assignments, now that duplicates have been merged. basically, just
    # renumber the cluster numbers
    redefined_cluster_assignments = renumber_cluster_assignments(temp_cluster_assignments)

    # Output pdbs based on the deduplicated redefined_cluster_assignments.
    for redefined_clusnum, redefined_clusmems in redefined_cluster_assignments.items():
        for elem in redefined_clusmems:
            pdb_name = all_AA_permuted_pdbpaths[elem]
            vdgscrrs  = all_vdg_scrrs[elem]
            original_temp_clusnum = get_original_temp_clusnum(temp_cluster_assignments, 
                                                              elem)
            is_centroid = elem in cents
            if is_centroid:
                pdb_name = pdb_name[:-4] + '_centroid.pdb'

            redun_temp_pdbpath, vdmAA_str = get_output_pdbpath(pdb_name, first_vdm_scrr, 
                    vdgscrrs, cluster_dirname, original_temp_clusnum, tempdir=True)
            
            new_redefined_pdbpath, vdmAA_str = get_output_pdbpath(pdb_name, first_vdm_scrr,
                    vdgscrrs, cluster_dirname, redefined_clusnum, tempdir=False)
            
            # rename (move) the old pdbpath to new pdbpath instead of writePDB because
            # cl_struct is mutable and its reference may have been modified.
            if os.path.exists(redun_temp_pdbpath):
                os.rename(redun_temp_pdbpath, new_redefined_pdbpath)
            else:
                equiv_temp_pdb = find_equivalent_pdb(redun_temp_pdbpath, new_redefined_pdbpath)
                os.rename(equiv_temp_pdb, new_redefined_pdbpath)

    # Delete the 'temp' reordered_AA dir.
    delete_temp_reordered_dir(redun_temp_pdbpath)
    return redefined_cluster_assignments

def find_equivalent_pdb(redun_temp_pdbpath, new_redefined_pdbpath):
    # Find equivalent binding site where the order of identical AAs are flipped.
    
    n_spl = sorted(os.path.basename(new_redefined_pdbpath).split('_'))
    n_spl = [i for i in n_spl if 'clus' not in i]
    dirs_upstream_of_clus = os.path.dirname(redun_temp_pdbpath).split(os.sep)[:-1]
    dirs_upstream_of_clus_path = '/'+os.path.join(*dirs_upstream_of_clus)
    for _root, _dirs, _files in os.walk(dirs_upstream_of_clus_path):
        for e in _files:
            e_spl = sorted(e.split('_'))
            e_spl = [i for i in e_spl if 'clus' not in i]
            if e_spl == n_spl:
                return os.path.abspath(os.path.join(_root, e)) 
    
    print('No equivalent PDB found for: ', redun_temp_pdbpath)
    sys.exit(1)

def merge_equivalent_AAperm_clusters(clusnumA_sorted_biounit_vdmscrrs, clusnumB_sorted_biounit_vdmscrrs,
                              temp_cluster_assignments, _clusnumA, clusnumA_mems, _clusnumB, 
                              clusnumB_mems):
    # Compare two clusters and merge if equivalent
    for _a in clusnumA_sorted_biounit_vdmscrrs:
        for _b in clusnumB_sorted_biounit_vdmscrrs:
            if _a == _b: # then merge clusnumA with clusnumB and delete clusnumB
                temp_cluster_assignments[_clusnumA] = np.append(clusnumA_mems, clusnumB_mems)
                del temp_cluster_assignments[_clusnumB]
                return temp_cluster_assignments
    return temp_cluster_assignments

def get_clus_sorted_biounits_and_vdmscrrs(cluster_dirname, num_vdms, vdmAA_str, _clusnum, tempdir=False):

    if tempdir:
        _outdir = os.path.join(cluster_dirname, 'temp', str(num_vdms), 
                                                vdmAA_str, f'cluster_{_clusnum}')
    else:
        _outdir = os.path.join(cluster_dirname, str(num_vdms), 
                                                vdmAA_str, f'cluster_{_clusnum}')
    clus_pdbnames = os.listdir(_outdir)
    clus_sorted_biounits_and_vdmscrrs = []
    for _pdbname in clus_pdbnames:
        biounit, list_scrrs = derive_list_scrr_from_pdbname(_pdbname)
        clus_sorted_biounits_and_vdmscrrs.append([biounit, sorted(list_scrrs)])
    return clus_sorted_biounits_and_vdmscrrs

def derive_list_scrr_from_pdbname(_pdbname):
    y = _pdbname.removesuffix('.pdb')
    delimited = y.split('_')
    delimited = [i for i in delimited if i != 'centroid']
    biounit = delimited[1:6]
    scrr_chain = delimited[6:]
    assert len(scrr_chain) % 4 == 0
    list_scrrs = [scrr_chain[i:i + 4] for i in range(0, len(scrr_chain), 4)]
    return [biounit, list_scrrs]

def get_output_pdbpath(pdb_name, first_vdm_scrr, vdgscrrs, cluster_dirname, cluster_num, 
                       tempdir=False):
    # tempdir is for temp_cluster_assignments, which will be redefined later when duplicate
    # clusters are merged.
    num_vdms = len(first_vdm_scrr)
    base_pdb = os.path.basename(pdb_name)
    base_pdbname = base_pdb.removesuffix('.pdb')
    vdmAA_str = '_'.join([x[3] for x in first_vdm_scrr])
    vdg_scrr_str = '_'.join(['_'.join([str(z) for z in v_s]) for v_s in vdgscrrs])
    if tempdir: 
        cluster_subdirname = os.path.join(cluster_dirname, 'temp', str(num_vdms), vdmAA_str, f'cluster_{cluster_num}')
    else:
        cluster_subdirname = os.path.join(cluster_dirname, str(num_vdms), vdmAA_str, f'cluster_{cluster_num}')
    os.makedirs(cluster_subdirname, exist_ok=True)
    output_pdbpath = os.path.join(cluster_subdirname, f'clus{cluster_num}_{base_pdbname}_{vdg_scrr_str}.pdb')
    return output_pdbpath, vdmAA_str

def flatten_AA_permuted_bbatoms(all_AA_permuted_bbatoms_to_align, all_AA_permuted_pdbpaths, 
                                all_AA_permuted_vdm_scrrs):
    flat_all_AA_permuted_bb_atoms_to_clus = []
    for AA_bindingsite_perm in all_AA_permuted_bbatoms_to_align:
        flattened_AA_perm = []
        for _vdmres in AA_bindingsite_perm:
            for _vdmresbbcoord in _vdmres:
                flattened_AA_perm.append(_vdmresbbcoord)
        flat_all_AA_permuted_bb_atoms_to_clus.append(flattened_AA_perm)
    assert len(all_AA_permuted_pdbpaths) == len(flat_all_AA_permuted_bb_atoms_to_clus)
    return flat_all_AA_permuted_bb_atoms_to_clus

def get_struct_data(_pdb, idxs, perms, _flat_bb_atoms_to_clus, coords, all_pdbs, all_structs, 
                    all_vdg_scrrs, _vdgscrr):
    # get data for each pdb/struct
    struct = pr.parsePDB(_pdb)
    occs = struct.getOccupancies()
    all_coords = struct.getCoords()
    cg_idxs = np.array(
        [np.argwhere(occs == 3. + 0.1 * idx)[0][0] 
         for idx in idxs]
    )
    coords_to_add = []
    for perm in perms:
        perm_coords = np.zeros((len(idxs), 3))
        perm_coords[:len(perm)] = all_coords[cg_idxs[perm]]
        # Comment out below, which is for `production` branch of vdG-miner and does not
        # sample permutations of identical AAs in binding site.
        '''
        for j, name in enumerate(['N', 'CA', 'C']):
            mask = np.logical_and.reduce((occs > 1., names == name))
            perm_coords[len(perm)+j::3] = \
                all_coords[mask][np.array(res_idxs)]
        '''
        # The below if for ligand-vdGs, which samples perms of identical AAs in binding site.
        for _bbatom in _flat_bb_atoms_to_clus:
            perm_coords = np.vstack([perm_coords, _bbatom])
        coords_to_add.append(perm_coords)
    coords += coords_to_add
    all_pdbs.append(_pdb)
    all_structs.append(struct.copy())
    all_vdg_scrrs.append(_vdgscrr)
    return coords, all_pdbs, all_structs, all_vdg_scrrs

def get_weights(idxs, res_idxs, N):
    # define equal weights for the CG and the mainchain atoms of 
    # interacting residues to be used in the Kabsch alignments
    weights = np.zeros(N)
    weights[:len(idxs)] = 0.5 / len(idxs)
    weights[len(idxs):] = 0.5 / (len(res_idxs) * 3)
    assert np.abs(weights.sum() - 1.) < 1e-8
    return weights

def align_and_cluster(M, idxs, res_idxs, N, coords, rmsd_cutoff):
    # get weights for the Kabsch alignments
    weights = get_weights(idxs, res_idxs, N)
    
    # find minimal-RMSD alignments between all pairs of structures
    triu_indices = np.triu_indices(M, 1)
    L = len(triu_indices[0])
    R, t, msd, best_perms = np.zeros((L, 3, 3)), np.zeros((L, 3)), \
                            10000. * np.ones(L), np.zeros(L, dtype=int)
    for i in range(coords.shape[0]):
        # iterate over atom permutations to find the best alignments
        _R, _t, _msd = kabsch(coords[i][triu_indices[0]], 
                              coords[0][triu_indices[1]], 
                              weights)
        R[_msd < msd] = _R[_msd < msd]
        t[_msd < msd] = _t[_msd < msd]
        best_perms[_msd < msd] = i
        msd[_msd < msd] = _msd[_msd < msd]
    A = np.eye(M, dtype=int)
    A[triu_indices] = (msd <= rmsd_cutoff ** 2).astype(int)
    A = A + A.T
    all_mems, cents = greedy(A)
    return all_mems, cents, triu_indices, R, t, msd

def process_AA_redun_cluster(cent, cluster, all_pdbs, all_structs, all_vdg_scrrs,
        triu_indices, R, t, msd, first_vdm_scrr, cluster_dirname, cluster_num, 
        num_supposed_to_be_output):
    # Process and output PDBs of this cluster. Note that they must be output to a 'temp'
    # dir because there may be degenerate binding sites that will be deduplicated later.
    assert cent in cluster, f'Centroid {cent} not in cluster {cluster}.'
    for el in cluster:
        pdb_name = all_pdbs[el]
        cl_struct = all_structs[el].copy()
        vdgscrrs  = all_vdg_scrrs[el]
        moved_cl_struct = None # refresh `moved_cl_struct` in case cl_struct was 
                               # accidentally modified in a previous iteration (because
                               # prody obj changes directly modify the reference obj)
        # create the environment PDB file
        if el == cent:
            pdb_name = pdb_name[:-4] + '_centroid.pdb'
            moved_cl_struct = cl_struct.copy()
        else:
            el_mobile = np.logical_and(triu_indices[0] == el,
                                       triu_indices[1] == cent)
            cent_mobile = np.logical_and(triu_indices[0] == cent,
                                         triu_indices[1] == el)
            if np.any(el_mobile):
                # non-centroid was mobile in the Kabsch alignment
                idx = np.argwhere(el_mobile).flatten()[0]
                _R = R[idx]
                _t = t[idx]
                rmsd = np.sqrt(msd[idx])
            else:
                # centroid was mobile in the Kabsch alignment
                idx = np.argwhere(cent_mobile).flatten()[0]
                _R = R[idx].T
                _t = -np.dot(t[idx], _R)
                rmsd = np.sqrt(msd[idx])

            moved_cl_struct = cl_struct.copy() # deepcopy because of prody obj mutability
                                               # to reference obj when setCoords is called
            moved_cl_struct.setCoords(
                np.dot(cl_struct.getCoords(), _R) + _t)
            #print('RMSD =', rmsd)
        
        # Output to a 'temp' dir because the clusters will be re-defined later, when 
        # we take care of duplicates (i.e. different AA permutations of vdms in the same, 
        # or homologous, binding sites will end up in different clusters, therefore 
        # appearing unique, but in reality, the order of identical AAs are simply flipped.) 
        # The way to check for this later is to check for equivalent vdm scrrs, e.g. if 
        # cluster 1 contains xxx_GLY1_GLY3 and cluster2 contains xxx_GLY3_GLY1, then
        # merge those clusters.
        output_pdbpath, vdmAA_str = get_output_pdbpath(pdb_name, first_vdm_scrr, vdgscrrs, 
                                            cluster_dirname, cluster_num, tempdir=True)
        pr.writePDB(output_pdbpath, moved_cl_struct)
        num_supposed_to_be_output += 1
    return num_supposed_to_be_output, vdmAA_str

def remove_redundant_AAperm_clusters(temp_cluster_assignments, cluster_dirname, num_vdms, vdmAA_str):
    from pprint import pprint
    already_visited_temp_clusnum_A = []
    num_of_clusters = list(temp_cluster_assignments.keys())
    for _clusnumA in num_of_clusters:
        already_visited_temp_clusnum_A.append(_clusnumA)
        if _clusnumA not in temp_cluster_assignments.keys(): # bc dict dynamically changing
            continue
        clusnumA_mems = temp_cluster_assignments[_clusnumA]
        for _clusnumB in num_of_clusters:
            if _clusnumA == _clusnumB:
                continue
            elif _clusnumB in already_visited_temp_clusnum_A: # checking that _clusnumB not 
                                # already seen as_clusnumA ensures that only the upper triangle 
                                # is calculated in this all x all protocol 
                continue
            elif _clusnumB not in temp_cluster_assignments.keys(): # bc dict dynamically changing
                continue
            else:
                clusnumB_mems = temp_cluster_assignments[_clusnumB]
                clusnumA_sorted_biounit_vdmscrrs = get_clus_sorted_biounits_and_vdmscrrs(
                                            cluster_dirname, num_vdms, vdmAA_str, _clusnumA, tempdir=True)
                clusnumB_sorted_biounit_vdmscrrs = get_clus_sorted_biounits_and_vdmscrrs(
                                            cluster_dirname, num_vdms, vdmAA_str, _clusnumB, tempdir=True)
                temp_cluster_assignments = merge_equivalent_AAperm_clusters(
                        clusnumA_sorted_biounit_vdmscrrs, clusnumB_sorted_biounit_vdmscrrs,
                        temp_cluster_assignments, _clusnumA, clusnumA_mems, _clusnumB, 
                        clusnumB_mems)
    return temp_cluster_assignments

def renumber_cluster_assignments(temp_cluster_assignments):
    redefined_cluster_assignments = {}
    temp_cluster_assignments_keys = sorted(list(temp_cluster_assignments.keys()))
    redefined_clusnum = 1
    for temp_clus_assignment_key in temp_cluster_assignments_keys:
        redefined_cluster_assignments[redefined_clusnum] = temp_cluster_assignments[
                                                           temp_clus_assignment_key]
        redefined_clusnum += 1
    return redefined_cluster_assignments 

def get_original_temp_clusnum(temp_cluster_assignments, elem):
    original_temp_clusnum = [_cnum for _cnum, _cmems in temp_cluster_assignments.items() if elem in _cmems]
    assert len(original_temp_clusnum) == 1
    return original_temp_clusnum[0]

def delete_temp_reordered_dir(redun_temp_pdbpath):
    pre, spl = redun_temp_pdbpath.split('temp/')
    spl = spl.split('/')
    temp_reordered_dir = os.path.join(pre, 'temp', spl[0], spl[1])
    assert os.path.exists(temp_reordered_dir)
    shutil.rmtree(temp_reordered_dir) # will remove the duplicated pdbs