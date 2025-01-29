'''
Removes redundant vdGs from the vdG library. vdGs are redundant if they meet all criteria: 
    1. CG and vdM AA identities and binding pose (backbones within a specified RMSD threshold)
    2. similar positions of the flanking residues (CAs within a specified RMSD 
       threshold)
    3. sequence similarity of residues flanking the vdMs (specified by seq. similarity 
       threshold). however, when the clustering is done, it's based on dissimilarity (because
       the clustering is distance-based)
    4. lastly, deduplicate vdgs that are identical but are featurized differently solely 
       because of different AA permutations of the same binding site
'''

import os
import sys
import argparse
import time
from itertools import combinations, permutations, product
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import prody as pr
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
from align_and_cluster import cluster_structures

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', "--cg", type=str, 
                        help="The common name for the chemical group. Defaults to the "
                        "SMARTS pattern.")
    parser.add_argument('-i', '--cg-indices-to-align', nargs='+', type=int, 
                        help='Indices of CG atoms on which to align.')
    parser.add_argument('-s', '--symmetry-classes', nargs='+', type=int,
                        help='Integers representing the symmetry classes of the CG atoms on '
                        'which clustering is to be performed. If provided, should have the same '
                        'length as idxs. If not provided, the atoms are assumed to be '
                        'symmetrically inequivalent.')
    parser.add_argument('-v', "--vdglib-dir", type=str, required=True,
                        help="Directory for the vdms of this CG.")
    parser.add_argument('-t', "--seq", default=0.50,
                        help="Sequence similarity threshold for clustering sequences "
                        "flanking the vdM to determine redundancy. This should be a value "
                        "between 0 and 1. Values > seq will be considered redundant.")
    parser.add_argument('-f', "--flank", default=5,
                        help="Number of residues flanking the vdms for which to calculate "
                        "sequence similarity and backbone similarity.")
    parser.add_argument('-l', "--logfile", help="Path to log file.")
    parser.add_argument('-o', "--output-clusters", type=bool, 
                        help="Determine whether to output clustered pdbs or not.")
    
    return parser.parse_args()

'''
Characterize all subsets of vdMs within each vdG, up to 4 vdMs (i.e., all quadruples, 
triplets, pairs, and singles).

Each key of the `vdm_combos` dict denotes how many residues (vdMs) within the vdG are 
being considered, and each value is a subdict whose keys are the combinations of vdM 
AAs and whose values are a featurized list of the vdG subsets defined by those 
combinations of vdM AAs. The features of the vdG subsets include the CG coordinates, 
the backbone (N, CA, C) coordinates of the vdMs, the sequences flanking those vdMs 
(+/- 5 residues), the PDB from which the vdG is derived, and the PDB segments/chains/
resnums of the CG and vdMs.

Example:

vdm_combos = { 
              4: {
                    (Ala, Ala, Cys, Phe): [
                   #                          [  [CG coords], 
                   #                             [bb N-Ca-C coords of Ala1, Ala2, Cys3, Phe4], 
                                                [seq. +/- 5 of Ala1, Ala2, Cys3, Phe4], 
                                                [CA   +/- 5 of Ala1, Ala2, Cys3, Phe4], 
                                                [PDB path],
                                                [PDB seg/chain/resnum of Ala1, Ala2, Cys3, Phe4]
                                             ], 
                                             ...
                                          ], 
                    (Ala, Ala, Cys, Trp): ...
                 }, 
              3: {
                    (Ala, Ala, Cys): ...,
                 } 
             }

Note that with there should be an order-preserving one-to-one mapping of each vdM in the 
lists that enumerate the vdM characteristics. When dealing with multiple vdMs of the same 
AA identity (such as having 2 Ala's in this example), all AA permutations must be sampled 
when checking for redundancy.
'''

def main():
   start_time = time.time()
   args = parse_args()
   CG = args.cg
   seq_sim_thresh = args.seq
   num_flanking = args.flank
   symmetry_classes = args.symmetry_classes
   cg_indices_to_align = args.cg_indices_to_align
   logfile = args.logfile
   vdglib_dir = args.vdglib_dir
   vdg_pdbs_dir = os.path.join(vdglib_dir, 'vdg_pdbs')
   out_dir = os.path.join(vdglib_dir, 'nr_vdgs')
   output_clusters = args.output_clusters

   os.makedirs(out_dir, exist_ok=True)
   if len(os.listdir(out_dir)) > 0:
      print(f'The nr_vdgs output dir {out_dir} is not empty. Remove files or define a new output '
            'dir name to prevent accidental overwriting.')
      sys.exit(1)

   if output_clusters:
      has_clustered_files = check_output_cluster_dirs(vdglib_dir)
      if has_clustered_files:
         sys.exit(1)

   with open(logfile, 'a') as file:
        file.write(f"{'='*20} Starting deduplicate_reun_vdgs.py run {'='*20} \n")
   
   # Initialize vdm_combos dict to store the subsets of vdMs within a vdG
   vdm_combos = {}

   # Iterate over the PDBs and CGs that were identified as containing the SMARTS group
   vdg_pdbs_in_dir = os.listdir(vdg_pdbs_dir)
   for pdbname in vdg_pdbs_in_dir:
      ####### for testing ########
      #if '6jsf' not in pdbname and '7myu' not in pdbname:
      #   continue
      ####### for testing ########
      pdbpath = os.path.join(vdg_pdbs_dir, pdbname)
      prody_obj = pr.parsePDB(pdbpath)
      cg_coords = get_cg_coords(prody_obj)
      # Characterize the vdM residues (bb coords, flanking residues, pdb paths, etc.)
      vdms_dict = get_vdm_res_features(prody_obj, pdbpath, num_flanking)
      # Determine the vdM combinations, up to 4 residues
      vdm_resinds = list(vdms_dict.keys())
      vdg_subsets = get_vdg_subsets(vdm_resinds)
      # Iterate over subsets
      for vdg_subset in vdg_subsets:
         # Record these features in the same order as in vdg_subset. Then, sort all based on
         # alphabetical order of the vdm AAs.
         re_ordered_aas, re_ordered_bbcoords, re_ordered_flankingseqs, re_ordered_CAs, \
            re_ordered_scrr = reorder_vdg_subset(vdg_subset, vdms_dict)
         # Add to `vdm_combos` dict
         vdm_combos = add_vdgs_to_dict(vdm_combos, vdg_subset, re_ordered_aas, 
            re_ordered_bbcoords, re_ordered_flankingseqs, re_ordered_CAs, re_ordered_scrr,
            cg_coords, pdbpath)
   # Evaluate the complete collection of vdGs and determine redundancy 
   all_nr_vdgs = []
   for num_vdms_in_subset, _subsets in vdm_combos.items():
      for _reordered_AAs, _vdgs in _subsets.items():
         # vdG subsets that have identical vdm AA compositions may be redundant
         if len(_vdgs) <= 1: # automatically not reundant b/c no other vdgs have this vdm combo
            single_vdg = _vdgs[0]
            all_nr_vdgs.append(single_vdg)
         else:
            nr_vdgs = get_nr_vdgs_of_same_AA_comp(_vdgs, seq_sim_thresh, _reordered_AAs, 
                                                  vdglib_dir, cg_indices_to_align,
                                                  symmetry_classes, output_clusters)
            for n_v in nr_vdgs:
               all_nr_vdgs.append(n_v)
   
   # There are many duplicates in `nr_vdgs`, because different AA permutations of vdms in the same 
   # binding site will end up in different clusters, therefore appearing to be unique. 
   # The last step is to account for these duplicates by gathering all vdg subset residues and 
   # ordering them by scrr to determine which vdgs are identical.
   
   namingdict = {}
   for nonredun_vdg in all_nr_vdgs:
      nonredun_vdg_pdb = nonredun_vdg[4]
      pdbcode = nonredun_vdg_pdb.rstrip('/').removesuffix('.pdb').split('/')[-1]
      if pdbcode not in namingdict.keys():
         namingdict[pdbcode] = []
      scrr = nonredun_vdg[5]
      sorted_scrr = sorted(scrr)
      if sorted_scrr not in namingdict[pdbcode]:
         namingdict[pdbcode].append(sorted_scrr)

   # Print out each vdg (separated by the # of vdms in the vdg)
   for pdbcode, list_scrrs in namingdict.items():
      for vdg_scrr in list_scrrs:
         vdg_scrr_str = '_'.join(['_'.join([str(z) for z in v_s]) for v_s in vdg_scrr])
         pr_obj = isolate_vdg_subset_obj(pdbcode, vdg_scrr, vdg_pdbs_dir)
         # Write out pdb
         write_vdg_subset(pdbcode, vdg_scrr, vdg_scrr_str, pr_obj, out_dir)
         
   # Print out time elapsed
   seconds = time.time() - start_time
   hours = seconds // 3600
   minutes = (seconds % 3600) // 60
   seconds = seconds % 60
   seconds = round(seconds, 2)
    
   # Write out results to log file
   with open(logfile, 'a') as file:
      num_vdg_pdbs = len(vdg_pdbs_in_dir)
      num_nr_singles =  len(os.listdir(os.path.join(out_dir, '1')))
      num_nr_pairs =    len(os.listdir(os.path.join(out_dir, '2')))
      num_nr_triplets = len(os.listdir(os.path.join(out_dir, '3')))
      num_nr_quad =     len(os.listdir(os.path.join(out_dir, '4')))
      file.write(f'\tNumber of vdg subsets from {num_vdg_pdbs} vdgs:\n')
      file.write(f'\t\tNonredundant single vdms: {num_nr_singles}\n')
      file.write(f'\t\tNonredundant vdm pairs: {num_nr_pairs}\n')
      file.write(f'\t\tNonredundant vdm triplets: {num_nr_triplets}\n')
      file.write(f'\t\tNonredundant vdm quadruplets: {num_nr_quad}\n')
      file.write(f"Completed deduplicate_redun_vdgs.py in {hours} h, ")
      file.write(f"{minutes} mins, and {seconds} secs \n") 

def isolate_vdg_subset_obj(pdbcode, vdg_scrr, vdg_pdbs_dir):
   # Select just the specified vdms (vdg subset) for printing out
   # Initialize prody obj for printing out
   input_vdg_pdbpath = os.path.join(vdg_pdbs_dir, pdbcode + '.pdb')
   par = pr.parsePDB(input_vdg_pdbpath)
   pr_obj = par.select('occupancy > 2.8') # occ >= 3 is CG
   for vdm_scrr in vdg_scrr:
      if vdm_scrr[0]: # has a pdb segment defined
         res_obj = par.select(
            f'segment {vdm_scrr[0]} and chain {vdm_scrr[1]} and resnum {vdm_scrr[2]}')
      else: # no segment
         res_obj = par.select(
            f'chain {vdm_scrr[1]} and resnum {vdm_scrr[2]}')
      assert len(set(res_obj.getResindices())) == 1
      vdm_res_name = list(set(res_obj.getResnames()))
      assert len(vdm_res_name)  == 1
      assert vdm_res_name[0] == vdm_scrr[3]
      pr_obj += res_obj
      # Ensure the correct number of vdms in this vdg
      num_resinds = len(set(pr_obj.getResindices()))
   assert num_resinds == len(vdg_scrr) + 1 # +1 is because CG counts as a resind
   return pr_obj

def write_vdg_subset(pdbcode, vdg_scrr, vdg_scrr_str, pr_obj, out_dir):
   n_vdms_in_subset = str(len(vdg_scrr))
   outputname = f'{n_vdms_in_subset}_{pdbcode}_{vdg_scrr_str}.pdb'
   outputpath = os.path.join(out_dir, n_vdms_in_subset, outputname)
   if str(n_vdms_in_subset) not in os.listdir(out_dir):
      os.mkdir(os.path.join(out_dir, n_vdms_in_subset))
   pr.writePDB(outputpath, pr_obj)

def get_vdg_permutations(reordered_AAs, _vdgs):
   # Note that these are permutations of AAs only, not symmetric atoms
   permuted_indices = permute_duplicates(reordered_AAs)
   all_AA_permuted_cg_coords = []
   all_AA_permuted_vdm_bbcoords = []
   all_AA_permuted_flankingseqs = []
   all_AA_permuted_flankingCAs = []
   all_AA_permuted_pdbpaths = []
   all_AA_permuted_vdm_scrr = []

   # Iterate over all AA permutations of each vdg
   for _vdg in _vdgs:
      # CG coords and pdbpaths remain unchanged, but vdmbbs, flankingseqs, flankingCAs, and scrrs
      # need to be permuted.
      for permutation in permuted_indices:
         all_AA_permuted_cg_coords.append(_vdg[0])
         all_AA_permuted_pdbpaths.append(_vdg[4])
         nonpermuted_vdmbb = _vdg[1]
         nonpermuted_flankingseqs = _vdg[2]
         nonpermuted_flankingCAs = _vdg[3]
         nonpermuted_vdm_scrr = _vdg[5]
         vdmbb_permutation = [nonpermuted_vdmbb[ix] for ix in permutation]
         flankingseqs_permutation = [nonpermuted_flankingseqs[ix] for ix in permutation]
         flankingCAs_permutation = [nonpermuted_flankingCAs[ix] for ix in permutation]
         vdm_scrrs_permutation = [nonpermuted_vdm_scrr[ix] for ix in permutation]
         all_AA_permuted_vdm_bbcoords.append(vdmbb_permutation)
         all_AA_permuted_flankingseqs.append(flankingseqs_permutation)
         all_AA_permuted_flankingCAs.append(flankingCAs_permutation)
         all_AA_permuted_vdm_scrr.append(vdm_scrrs_permutation)

   return all_AA_permuted_cg_coords, all_AA_permuted_vdm_bbcoords, all_AA_permuted_flankingseqs, \
      all_AA_permuted_flankingCAs, all_AA_permuted_pdbpaths, all_AA_permuted_vdm_scrr

def elements_in_clusters(indices_of_elements_in_cluster, cg_coords, vdm_bbcoords, flankingseqs, 
   flankingCAs, pdbpaths, vdm_scrrs):
   assert len(cg_coords) == len(pdbpaths)
   clus_cg_coords = [cg_coords[index] for index in indices_of_elements_in_cluster]
   clus_vdmbb_coords = [vdm_bbcoords[index] for index in indices_of_elements_in_cluster]
   clus_flankingseqs = [flankingseqs[index] for index in indices_of_elements_in_cluster]
   clus_flankingCAs = [flankingCAs[index] for index in indices_of_elements_in_cluster]
   clus_pdbpaths = [pdbpaths[index] for index in indices_of_elements_in_cluster]
   clus_vdm_scrr = [vdm_scrrs[index] for index in indices_of_elements_in_cluster]
   return clus_cg_coords, clus_vdmbb_coords, clus_flankingseqs, \
      clus_flankingCAs, clus_pdbpaths, clus_vdm_scrr

def permute_duplicates(seq):
    # Dictionary to store indices for each element in the sequence.
    # These are for AA permutations, not permutations of symmetric CG atoms.
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
    for perm_combination in product(*[permutations(group) for group in permute_groups]):
        # start with the list of original indices
        permuted_idx = list(range(len(seq)))

        # flatten the product of permutations and assign them to the corresponding positions
        for group_idx, perm in zip(permute_groups, perm_combination):
            for orig_idx, new_idx in zip(group_idx, perm):
                permuted_idx[orig_idx] = new_idx

        permuted_idx_lists.append(permuted_idx)
    
    return permuted_idx_lists

def combine_cg_and_vdmbb_coords(all_AA_permuted_cg_coords, all_AA_permuted_vdm_bbcoords):
   all_cg_and_vdmbb_coords = []
   for _cg, _vdmbb_per_res in zip(all_AA_permuted_cg_coords, all_AA_permuted_vdm_bbcoords):
      vdg_cg_amd_vdmbb = []
      flattened_cg_coords = flatten_cg_coords(_cg)
      flattened_vdmbbs = flatten_vdg_bbs(_vdmbb_per_res)
      cg_and_vdmbb = flattened_cg_coords + flattened_vdmbbs
      cg_and_vdmbb = np.array(cg_and_vdmbb)
      all_cg_and_vdmbb_coords.append(cg_and_vdmbb)
   return all_cg_and_vdmbb_coords

def get_nr_vdgs_of_same_AA_comp(_vdgs, seq_sim_thresh, reordered_AAs, vdglib_dir, 
                                cg_indices_to_align, symmetry_classes, output_clusters):
   # vdG subsets that have identical vdm AA compositions may be redudant.
   # First, get all permutations for equivalent AAs. For example, if the vdms are 
   # [Ala1, Ala2, Glu], then we need to sample [Ala1, Ala2, Glu] and [Ala2, Ala1, Glu].
   # The easiest way to get the permutations is to label by the index of the ordered AAs
   
   all_AA_permuted_cg_coords, all_AA_permuted_vdm_bbcoords, all_AA_permuted_flankingseqs, \
      all_AA_permuted_flankingCAs, all_AA_permuted_pdbpaths, all_AA_permuted_vdm_scrr = \
      get_vdg_permutations(reordered_AAs, _vdgs) # no longer need all_AA_permuted_cg_coords and 
                                                 # all_AA_permuted_vdm_bbcoords b/c they'll be 
                                                 # taken care of by symmetry-aware 
                                                 # align_and_cluster.py
   
   # Analyze all AA permutations of the vdgs.
   # First, calc rmsd between all AA permutations.

   # Cluster these potential redundant permuted vdGs by rmsd of CG + vdm bb coords. If their 
   # rmsds are larger than the rmsd thershold, then they are automatically considered different 
   # vdGs. `cluster_assignments` is a dict where each key is the cluster number, and each
   # value is a list of the indices of the vdGs (indexing is relative position in 
   # all_AA_permuted_cg_coords, all_AA_permuted_vdm_bbcoords, etc.) that belong in that cluster.
   
   num_vdms = len(reordered_AAs)
   first_vdg_cg_coords = _vdgs[0][0]
   num_cg_atoms = len(first_vdg_cg_coords)
   num_coords_to_align = num_cg_atoms + (num_vdms * 3)
   rmsd_cutoff = normalize_rmsd(num_coords_to_align)
   print(f'rmsd for {num_coords_to_align} atoms to align: {rmsd_cutoff}')
   if cg_indices_to_align is None:
      cg_indices_to_align = list(range(num_cg_atoms))

   cluster_structures(vdglib_dir, rmsd_cutoff, cg_indices_to_align, all_AA_permuted_pdbpaths, 
                      all_AA_permuted_vdm_scrr, out_clus_dir='cg_vdmbb_clusters',
                      all_AA_permuted_bbatoms_to_align=all_AA_permuted_vdm_bbcoords,
                      symmetry_classes=symmetry_classes)

   #cgvdmbb_cluster_assignments = get_hierarchical_clusters(all_cg_and_vdmbb_coords, 
   #   None_in_coords=True) # `cluster_assignments` dict not ordered by cluster size

   #nr_vdgs = []
   #for cgvdmbb_clusnum, indices_of_elements_in_cg_vdmbb_cluster in \
   #   cgvdmbb_cluster_assignments.items():
   #   # Within these clusters (because each of these clusters are nonredundant from each
   #   # other in terms of CG and vdm bb positions), cluster again based on backbone stretches
   #   # flanking the vdms.
   #   cgvdmbb_clus_cg_coords, cgvdmbb_clus_vdmbb_coords, cgvdmbb_clus_flankingseqs, \
   #      cgvdmbb_clus_flankingCAs, cgvdmbb_clus_pdbpaths, cgvdmbb_clus_vdm_scrr = \
   #      elements_in_clusters(indices_of_elements_in_cg_vdmbb_cluster, 
   #      all_AA_permuted_cg_coords, all_AA_permuted_vdm_bbcoords, all_AA_permuted_flankingseqs,
   #      all_AA_permuted_flankingCAs, all_AA_permuted_pdbpaths, all_AA_permuted_vdm_scrr)
   #   
   #   # If there's only 1 element in the cluster, then it's automatically nonredundant b/c its
   #   # CG+vdm bb position (binding pose) doesn't match any other vdG's binding pose
   #   if len(indices_of_elements_in_cg_vdmbb_cluster) <= 1:
   #      vdg_descript = [cgvdmbb_clus_cg_coords[0], cgvdmbb_clus_vdmbb_coords[0], 
   #         cgvdmbb_clus_flankingseqs[0], cgvdmbb_clus_flankingCAs[0], cgvdmbb_clus_pdbpaths[0], 
   #         cgvdmbb_clus_vdm_scrr[0]] # these 6 elements define a vdg
   #      nr_vdgs.append(vdg_descript)
   #      continue

   #   # Otherwise, if there's >1 element in the cluster, it's possible that the vdgs within 
   #   # this cluster are redundant. Next thing to look at is rmsd of those flanking 
   #   # bb stretches. Align and cluster.
   #   flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster = []
   #   for vdg_flankingCAs in cgvdmbb_clus_flankingCAs:
   #      # `vdg_flankingCAs` is a list of vdm residues, so flatten it
   #      flat_flanking_CAs = []
   #      for res in vdg_flankingCAs:
   #         for CA_coord in res:
   #            flat_flanking_CAs.append(CA_coord)
   #      flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster.append(flat_flanking_CAs)
   #   # Cluster the vdgs in `flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster`. In order to
   #   # calculate rmsd for clustering, take care of the instances where coordinates are 
   #   # represented as "None" because the stretch of +/- 5 AAs was missing N-terminal or 
   #   # C-terminal residues. Only select residues that have CA coords in both of the vdgs whose
   #   # rmsd is being calculated.
   #   flankingbb_cluster_assignments = get_hierarchical_clusters(
   #      flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster, None_in_coords=True)
   #   for flankingCAs_clusnum, indices_of_elements_in_flankingCAs_cluster in \
   #      flankingbb_cluster_assignments.items():
   #      # For each cluster determined by rmsd between CA coords flanking the vdms, gather the
   #      # vdg features (cg coords, vdmbb coords, etc.) corresponding to the vdgs belonging in
   #      # that cluster. The vdgs are called based on their indices within the clusters defined 
   #      # in the first step, when the rmsds of cg + vdm bbs were being compared.
   #      flankingCAs_clus_cg_coords, flankingCAs_clus_vdmbb_coords, \
   #         flankingCAs_clus_flankingseqs, flankingCAs_clus_flankingCAs, \
   #         flankingCAs_clus_pdbpaths, flankingCAs_clus_vdm_scrr = elements_in_clusters(
   #         indices_of_elements_in_flankingCAs_cluster, cgvdmbb_clus_cg_coords, 
   #         cgvdmbb_clus_vdmbb_coords, cgvdmbb_clus_flankingseqs, cgvdmbb_clus_flankingCAs,
   #         cgvdmbb_clus_pdbpaths, cgvdmbb_clus_vdm_scrr)

   #      # If there's only 1 element in the cluster, then it's automatically nonredundant b/c 
   #      # its binding site backbone doesn't match any other vdG's binding site
   #      if len(indices_of_elements_in_flankingCAs_cluster) <= 1:
   #         vdg_descript = [flankingCAs_clus_cg_coords[0], flankingCAs_clus_vdmbb_coords[0],
   #            flankingCAs_clus_flankingseqs[0], flankingCAs_clus_flankingCAs[0],
   #            flankingCAs_clus_pdbpaths[0], flankingCAs_clus_vdm_scrr[0]] # these 6 define a vdg
   #         nr_vdgs.append(vdg_descript)
   #         continue

   #      # Otherwise, if there's more than 1 element in the cluster, it's possible that the vdgs
   #      # within this cluster are redundant. Next metric to cluster on is sequence dissimilarity 
   #      # of the residues +/- of the vdms.
   #      flattened_flankingseqs_for_vdgs_in_flankingCA_clus = []
   #      
   #      for vdg_flankingseq in flankingCAs_clus_flankingseqs:
   #         flat_vdg_flankingseq = []
   #         for vdm_res in vdg_flankingseq:
   #            flat_vdg_flankingseq += vdm_res
   #         flattened_flankingseqs_for_vdgs_in_flankingCA_clus.append(flat_vdg_flankingseq)
   #      seqsim_clus_assignments = get_hierarchical_clusters(
   #         flattened_flankingseqs_for_vdgs_in_flankingCA_clus, seq_sim_thresh=seq_sim_thresh)
   #      # If the vdgs belong in the same cluster based on sequence similarity (and therefore
   #      # backbones and cg+vdmbb binding pose), then the vdgs are redundant. Select only one
   #      # to add to the `nr_vdgs` list.
   #      for seqsim_clusnum, indices_of_elements_in_flankingseq_clus in \
   #         seqsim_clus_assignments.items():
   #         flankingseq_clus_cg_coords, flankingseq_clus_vdmbb_coords, \
   #            flankingseq_clus_flankingseqs, flankingseq_clus_flankingCAs, \
   #            flankingseq_clus_pdbpaths, flankingseq_clus_vdm_scrr = \
   #            elements_in_clusters(indices_of_elements_in_flankingseq_clus, \
   #            flankingCAs_clus_cg_coords, flankingCAs_clus_vdmbb_coords, \
   #            flankingCAs_clus_flankingseqs, flankingCAs_clus_flankingCAs, \
   #            flankingCAs_clus_pdbpaths, flankingCAs_clus_vdm_scrr)
   #        
   #         if len(indices_of_elements_in_flankingseq_clus) <= 1:
   #            vdg_descript = [flankingseq_clus_cg_coords[0], flankingseq_clus_vdmbb_coords[0], 
   #            flankingseq_clus_flankingseqs[0], flankingseq_clus_flankingCAs[0], 
   #            flankingseq_clus_pdbpaths[0], flankingseq_clus_vdm_scrr[0]]
   #            nr_vdgs.append(vdg_descript)
   #         
   #         else: # select only one, so use index 0 for indices_of_elements_in_flankingseq_clus
   #            flankingseq_clus_cg_coords, flankingseq_clus_vdmbb_coords, \
   #            flankingseq_clus_flankingseqs, flankingseq_clus_flankingCAs, \
   #            flankingseq_clus_pdbpaths, flankingseq_clus_vdm_scrr = \
   #            elements_in_clusters([0], \
   #            flankingCAs_clus_cg_coords, flankingCAs_clus_vdmbb_coords, \
   #            flankingCAs_clus_flankingseqs, flankingCAs_clus_flankingCAs, \
   #            flankingCAs_clus_pdbpaths, flankingCAs_clus_vdm_scrr)
   #            vdg_descript = [flankingseq_clus_cg_coords[0], flankingseq_clus_vdmbb_coords[0],
   #               flankingseq_clus_flankingseqs[0], flankingseq_clus_flankingCAs[0], 
   #               flankingseq_clus_pdbpaths[0], flankingseq_clus_vdm_scrr[0]]
   #            nr_vdgs.append(vdg_descript)
   #           
   return nr_vdgs 

def calc_seq_similarity(list1, list2):
    # Count how many residues match
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

    # Calculate the percentage of matches
    if len(list1) == 0:
        match_percentage = 0
    else:
        match_percentage = (matches / len(list1)) * 100
    return match_percentage

def get_overlapping_res_coords(list1, list2):
    list1_overlap = []
    list2_overlap = []
    
    # Loop over both lists and check if both elements at the same index are not None
    for i in range(len(list1)):
        if list1[i] is not None and list2[i] is not None:
            list1_overlap.append(list1[i])
            list2_overlap.append(list2[i])
    
    return np.array(list1_overlap), np.array(list2_overlap)

def get_hierarchical_clusters(data, None_in_coords=False, seq_sim_thresh=None):
   # Returns dict where key = cluster number and value = indices of the list elements that
   # belong in that cluster. This dict is not ordered by size.
   # `data` is either a list of coords, or a list of sequences.
   # None_in_coords is a param that distinguishes whether `data` is coordinates or sequences.
   # This distinction is important because rmsd cannot be calculated on non-coordinates, so if 
   # any residues are "None" for either of the 2 vdgs being compared, that whole residue is 
   # discounted.

   num_coords = len(data[0])
   thresh = normalize_rmsd(num_coords)
   dist_matrix = create_dist_matrix(data, None_in_coords) # dist is either rmsd or % seq sim
   condensed_dist_matrix = squareform(dist_matrix)
   Z = linkage(condensed_dist_matrix, method='single')
   if seq_sim_thresh:
      seq_dissimilarity_thresh = 1 - seq_sim_thresh
      clusters = fcluster(Z, seq_dissimilarity_thresh, criterion='distance')
   else:
      clusters = fcluster(Z, thresh, criterion='distance')
   cluster_assignments = {} # key = clusnum, value = indices of `coords` belonging to that cluster
   for all_vdgs_index, clusnum in enumerate(clusters):
      if clusnum not in cluster_assignments.keys():
         cluster_assignments[clusnum] = []
      cluster_assignments[clusnum].append(all_vdgs_index)
   
   return cluster_assignments  

def create_dist_matrix(data, None_in_coords):
   n = len(data)
   distance_matrix = np.zeros((n, n))
   for i in range(n):
      for j in range(i + 1, n):
         # Align and then calc rmsd or seq dissimilarity. If "None" is in the list of coords 
         # (which is what would happen if rmsd is being calculated on stretches of backbone
         # and terminal residues are missing), then select only the overlapping residues 
         # where a coordinate could be extracted for both elements being measured.
         data_i = data[i]
         data_j = data[j]
         if None_in_coords: # calc rmsd 
            data_i, data_j = get_overlapping_res_coords(data_i, data_j)
            moved_coords_i, transf = pr.superpose(data_i, data_j)
            rmsd = pr.calcRMSD(moved_coords_i, data_j)
            value = rmsd
         else: # calc seq similarity
            seq_sim = calc_seq_similarity(data_i, data_j)
            value = (100 - seq_sim) / 100 # clustering is distance-based, so the metric is DISsimilarity
         distance_matrix[i][j] = value 
         distance_matrix[j][i] = value   # Symmetric matrix
   return distance_matrix

def flatten_cg_coords(_cg):
   coords = []
   for atom in _cg:
      coords.append(atom)
   return coords

def flatten_vdg_bbs(_vdmbb):
   # Flatten backbones to a single list where each element is a numpy array of the x, y, z 
   # coords of each vdm backbone atom
   bbs = []
   for res in _vdmbb:
      for atom in res:
         bbs.append(atom)
   return bbs

def add_vdgs_to_dict(vdm_combos, vdg_subset, re_ordered_aas, re_ordered_bbcoords, 
                     re_ordered_flankingseqs, re_ordered_CAs, re_ordered_scrr, cg_coords, pdbpath):
   # Construct the `vdm_combos` dict by categorizing subset of vdgs based on the number of vdms in
   # each subset, as well as the AA compositions of the vdms. Store each vdg subset as a list
   # containing [CG coords, backbone N-CA-C coords of the vdms, sequences flanking the vdms,
   # CA coordinates flanking the vdms, the PDB path, PDB seg/chain/resnum of the vdms]
   num_vdms_in_subset = len(vdg_subset)
   if num_vdms_in_subset not in vdm_combos.keys():
      vdm_combos[num_vdms_in_subset] = {}
   re_ord_aas_tup = tuple(re_ordered_aas)
   if re_ord_aas_tup not in vdm_combos[num_vdms_in_subset].keys():
      vdm_combos[num_vdms_in_subset][re_ord_aas_tup] = []
   # Add the vdg to `vdm_combos`
   vdm_combos[num_vdms_in_subset][re_ord_aas_tup].append([cg_coords, re_ordered_bbcoords,
      re_ordered_flankingseqs, re_ordered_CAs, pdbpath, re_ordered_scrr])

   return vdm_combos

def reorder_vdg_subset(vdg_subset, vdms_dict):
   '''Reorders vdms alphabetically.'''
   # First, record
   aas_of_vdms_in_order = []
   bb_coords_of_vdms_in_order = []
   flankingseqs_of_vdms_in_order = []
   flanking_CA_coords_of_vdms_in_order = []
   seg_ch_res_of_vdms_in_order = []
   for _vdmresind in vdg_subset:
      vdmAA, vdm_features = vdms_dict[_vdmresind]
      aas_of_vdms_in_order.append(vdmAA)
      vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict = vdm_features
      seg_ch_res_of_vdms_in_order.append(vdm_seg_chain_resnum_resname)
      bb_coords_of_vdms_in_order.append(bb_coords)
      flankingseqs = []
      flankingCAs = []
      sorted_flank_indices = sorted(list(flanking_seq_dict.keys()))
      # Decompress the flanking AA and CA info
      for flank_ind in sorted_flank_indices:
         flank_resname, flank_ca = flanking_seq_dict[flank_ind]
         flankingseqs.append(flank_resname)
         flankingCAs.append(flank_ca)
      flankingseqs_of_vdms_in_order.append(flankingseqs)
      flanking_CA_coords_of_vdms_in_order.append(flankingCAs)
   # Then, re-order based on alphabetical order
   super_list = [aas_of_vdms_in_order, bb_coords_of_vdms_in_order, 
      flankingseqs_of_vdms_in_order, flanking_CA_coords_of_vdms_in_order, 
      seg_ch_res_of_vdms_in_order]
   return sort_vdGs_by_AA(super_list)

def sort_vdGs_by_AA(super_list):
    # Check if all sublists have equal length
    assert all(len(sublist) == len(super_list[0]) for sublist in super_list)
    # Combine the sublists into a list of tuples, where each tuple corresponds to the 
    # elements at the same index
    combined = list(zip(*super_list))
    # Sort the combined list based on the first element (from the first sublist)
    sorted_combined = sorted(combined, key=lambda x: x[0])
    # Unzip the sorted combined list back into sublists
    sorted_sublists = list(zip(*sorted_combined))
    return [list(sublist) for sublist in sorted_sublists]

def get_vdm_res_features(prody_obj, pdbpath, num_flanking):
   # Identify the vdM residues (occ == 2). To be safe, select > 1.5 and < 2.5.
   vdm_residues = prody_obj.select('(occupancy) > 1.5 and (occupancy < 2.5)')
   vdm_resinds = set(vdm_residues.getResindices())
   # Record features of the vdm residues (bb coords, flanking residues, pdb paths, etc.)
   vdms_dict = {}
   for vdm_resind in vdm_resinds:
      vdm_obj = vdm_residues.select(f'resindex {vdm_resind}')
      # BB coords
      bb_coords = get_bb_coords(vdm_obj)

      # Sequence of (contiguous) flanking residues
      flanking_seq_dict = {} # key = relative flank num (-1, +1, etc.), 
                             # value = list(CA coords, AA identity)
      # Walk up and down the flanking residues and check that they are actually
      # neighboring the vdM, and not jumped through a chain break. If there's a
      # chain break, report the AA as "X".
      # > First, store the CA coords.
      for flank_num in range(1, num_flanking + 1):
         negative_flank_num = -1 * flank_num
         positive_flank_num = flank_num
         # Adding the flank_num (neg. or pos.) gives you the residue index
         for f in [negative_flank_num, positive_flank_num]:
            # Get the AA identity and CA coords of the "current" resindex
            current_resindex = vdm_resind + f
            AA, CA_coords = get_AA_and_CA_coords(prody_obj, current_resindex)
            flanking_seq_dict[f] = [AA, CA_coords]
      # Get the CA coords of the central vdM as well
      flanking_seq_dict[0] = ['vdm', bb_coords[1]] # label as 'vdm' for easy exclusion
                                                   # when calculating seq. similarity
      # > Then, go through the dict for chain breaks. Start in the fwd direction, and
      #   then progress backward. 
      fwd_rel_inds = range(1, num_flanking + 1)
      back_rel_inds = [-1 * i for i in fwd_rel_inds]
      central_vdm_CA = flanking_seq_dict[0][1]
      for list_indices in [fwd_rel_inds, back_rel_inds]:
         for ind in list_indices:
            # Get distance between current flanking res (relative to central vdm)
            # and res prior to it.
            if ind == 1 or ind == -1:
               prev_CA = central_vdm_CA
            curr_CA = flanking_seq_dict[ind][1]
            if curr_CA is None:
               flanking_seq_dict = found_chain_break(flanking_seq_dict, ind)
               break 
            dist = pr.calcDistance(np.array(prev_CA), np.array(curr_CA))
            # If the dist is > 4.5A, then it's a chain break.
            if dist > 4.5:
               flanking_seq_dict = found_chain_break(flanking_seq_dict, ind)
               break 
            # Otherwise, continue walking.
            prev_CA = curr_CA
      # Get this vdm resind's PDB identifier (seg, chain, resnum)
      vdm_seg_chain_resnum_resname = get_res_iden(vdm_obj)
      # Store all into dict
      vdm_descript = [vdm_seg_chain_resnum_resname, bb_coords, flanking_seq_dict]
      vdm_AA = vdm_seg_chain_resnum_resname[-1]
      assert vdm_resind not in list(vdms_dict.keys())
      vdms_dict[vdm_resind] = [vdm_AA, vdm_descript]
   return vdms_dict

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
      flanking_seq_dict[overwrite_ind] = ['X', None]
   return flanking_seq_dict

def get_AA_and_CA_coords(prody_obj, current_resindex):
   if current_resindex < 0: 
      sel_str = f'resindex `{current_resindex}`'
   else:
      sel_str = f'resindex {current_resindex}'

   curr_resindex_obj = prody_obj.select(sel_str)
   if curr_resindex_obj is None:
      AA = 'X'
      CA_coords = None
   elif curr_resindex_obj.protein is None: # not a residue
      AA = 'X'
      CA_coords = None
   else:
      curr_res_AA = list(set(curr_resindex_obj.getResnames()))
      assert len(curr_res_AA) == 1
      AA = curr_res_AA[0]
      CA_obj = curr_resindex_obj.select(sel_str + ' and name CA')
      assert len(CA_obj) == 1
      CA_coords = CA_obj.getCoords()[0]
   return AA, CA_coords

def get_cg_coords(prody_obj):
   # The CG atoms have their occupancies set to >= 3.0, with unique values (e.g., 3.0, 
   # 3.1, 3.2, etc.) to allow a 1:1 correspondence of equivalent atoms between different
   # ligands.
   cg = prody_obj.select('occupancy >= 3.0')
   num_atoms = len(cg)
   cg_coords = []
   for ind in range(num_atoms):
      occ = f'3.{ind}'
      atom = cg.select(f'occupancy == {occ}')
      assert len(atom) == 1
      atom_coords = atom.getCoords()[0]
      cg_coords.append(atom_coords)

   return np.array(cg_coords)

def get_bb_coords(obj):
   bb_coords = []
   for atom in ['N', 'CA', 'C']:
      atom_obj = obj.select(f'name {atom}')
      assert len(atom_obj) == 1
      coord = atom_obj.getCoords()[0]
      bb_coords.append(coord)
   return bb_coords

def get_vdg_subsets(input_list):
    # Group by quadruples, triples, pairs, and singles.
    # Initialize an empty list to store all subsets
    all_subsets = []
    # Loop through subset sizes 1 to 4, generate combos, then add to list
    for r in range(1, 5):
        subsets = combinations(input_list, r)
        all_subsets.extend(subsets)
    return all_subsets

def normalize_rmsd(num_atoms):
    # Smoothly scale rmsd threshold by the number of atoms.
    # Lower bound 0.4A with 10 atoms and upper bound 1A with 20 atoms.
    min_atoms = 10
    min_threshold = 0.4
    max_atoms = 20
    max_threshold = 1 
    
    if num_atoms < min_atoms:
        return min_threshold  # below 5 atoms, use the minimum threshold (1.5 Å)
    if num_atoms > max_atoms:
        return max_threshold  # above 20 atoms, use the maximum threshold (2.5 Å)
    
    # Linear scaling for atoms between 5 and 20
    scaling_factor = (num_atoms - min_atoms) / (max_atoms - min_atoms)
    threshold = min_threshold + scaling_factor * (max_threshold - min_threshold)
    
    return threshold

def check_output_cluster_dirs(vdglib_dir):
   # Check if the output directories are empty. If they are not, exit the script to prevent
   # accidental overwriting of files.
   has_files = False
   for out_clus_dir in ['cg_vdmbb_clusters', 'flankingbb_clusters', 'flankingseq_clusters']:
      cluster_dirname = os.path.join(vdglib_dir, out_clus_dir)
      os.makedirs(cluster_dirname, exist_ok=True)
      for root, dirs, files in os.walk(cluster_dirname):
         if files:
            has_files = True
            print('Files already exist in the output directory. '
               'Exiting to prevent accidental overwriting. '
               'Please remove the files or specify a different output directory.')
            return has_files
   return has_files      

if __name__ == "__main__":
    main()
