'''
WARNING: incomplete.

Removes redundant vdGs from the vdG library. vdGs are redundant if they meet all 
criteria: 
    1. CG and vdM AA identities and binding pose (backbones within a specified RMSD 
       threshold)
    2. similar positions of the flanking residues (CAs within a specified RMSD 
       threshold)
    3. sequence similarity of residues flanking the vdMs (specified by seq. 
       similarity ithreshold). however, when the clustering is done, it's based on 
       dissimilarity (because the clustering is distance-based)
    4. lastly, deduplicate vdgs that are identical but are featurized differently 
       solely because of different AA permutations of the same binding site
'''

import os
import sys
import argparse
import time
from itertools import combinations, permutations, product
import numpy as np
import prody as pr
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import align_and_cluster as clust
import utils

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', "--cg", type=str, 
                        help="The common name for the chemical group. Defaults to "
                        "the SMARTS pattern.")
    parser.add_argument('-v', "--vdglib-dir", type=str, required=True,
                        help="Directory for the vdms of this CG.")
    parser.add_argument('-o', "--output-clus-pdbs", action='store_true', 
                        help="Output clustered PDBs.")
    parser.add_argument('-q', "--seq", default=0.50,
                        help="Sequence similarity threshold for clustering sequences "
                        "flanking the vdM to determine redundancy. This should be a "
                        "value between 0 and 1. Values > seq will be considered "
                        "redundant.")
    parser.add_argument('-f', "--flank", default=5,
                        help="Number of residues flanking the vdms for which to "
                        "calculate sequence similarity and backbone similarity.")
    parser.add_argument('-s', '--symmetry-classes', nargs='+', type=int,
                        help='Integers representing the symmetry classes of the CG '
                        'atoms on which clustering is to be performed. If provided, '
                        'should have the same length as idxs. If not provided, the '
                        'atoms are assumed to be symmetrically inequivalent.')
    parser.add_argument('-l', "--logfile", help="Path to log file.")
    
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
                                         [ [CG coords], 
                                           [bb N-Ca-C of Ala1, Ala2, Cys3, Phe4], 
                                           [seq. +/- 5 of Ala1, Ala2, Cys3, Phe4], 
                                           [CA   +/- 5 of Ala1, Ala2, Cys3, Phe4], 
                                           [PDB path],
                                           [PDB scrr of Ala1, Ala2, Cys3, Phe4]
                                             ], 
                                             ...
                                          ], 
                  (Ala, Ala, Cys, Trp): ...
                 }, 
              3: {
                  (Ala, Ala, Cys): ...,
                 } 
             }

Note that with there should be an order-preserving one-to-one mapping of each vdM in  
the lists that enumerate the vdM characteristics. When dealing with multiple vdMs of 
the same AA identity (such as having 2 Ala's in this example), all AA permutations 
must be sampled when checking for redundancy.
'''

def main():
   start_time = time.time()
   args = parse_args()
   CG = args.cg
   seq_sim_thresh = args.seq
   num_flanking = args.flank
   symmetry_classes = args.symmetry_classes
   logfile = args.logfile
   output_clus_pdbs = args.output_clus_pdbs
   print()
   print('output_clus_pdbs:', output_clus_pdbs)
   vdglib_dir = args.vdglib_dir
   vdg_pdbs_dir = os.path.join(vdglib_dir, 'vdg_pdbs')
   out_dir = os.path.join(vdglib_dir, 'nr_vdgs')
   utils.handle_existing_files(out_dir)

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
      cg_coords = clust.get_cg_coords(prody_obj)
      # Characterize the vdM residues (bb coords, flanking residues, pdb paths, etc.)
      vdms_dict = clust.get_vdm_res_features(prody_obj, pdbpath, num_flanking)
      # Determine the vdM combinations, up to 4 residues
      vdm_resinds = list(vdms_dict.keys())
      vdg_subsets = clust.get_vdg_subsets(vdm_resinds)
      # Iterate over subsets
      for vdg_subset in vdg_subsets:
         # Record these features in the same order as in vdg_subset. Then, sort all 
         # based on alphabetical order of the vdm AAs.
         re_ordered_aas, re_ordered_bbcoords, re_ordered_flankingseqs, \
            re_ordered_CAs, re_ordered_scrr = clust.reorder_vdg_subset(
            vdg_subset, vdms_dict)
         # Add to `vdm_combos` dict
         vdm_combos = clust.add_vdgs_to_dict(vdm_combos, vdg_subset, re_ordered_aas, 
            re_ordered_bbcoords, re_ordered_flankingseqs, re_ordered_CAs, 
            re_ordered_scrr, cg_coords, pdbpath)
   
   # Evaluate the complete collection of vdGs and determine redundancy 
   all_nr_vdgs = []
   for num_vdms_in_subset, _subsets in vdm_combos.items():
      for _reordered_AAs, _vdgs in _subsets.items():
         # vdG subsets that have identical vdm AA compositions may be redundant.
         if len(_vdgs) <= 1: # Automatically not reundant b/c no other vdgs have 
                             # this vdm combo
            single_vdg = _vdgs[0]
            all_nr_vdgs.append(single_vdg)
         else:
            nr_vdgs = get_nr_vdgs_of_same_AA_comp(_vdgs, seq_sim_thresh, 
                  _reordered_AAs, symmetry_classes, output_clus_pdbs, vdglib_dir)
            for n_v in nr_vdgs:
               all_nr_vdgs.append(n_v)
   
   # There are many duplicates in `nr_vdgs`, because different AA permutations of 
   # vdms in the same binding site will end up in different clusters, therefore 
   # appearing to be unique. The last step is to account for these duplicates by 
   # gathering all vdg subset residues and ordering them by scrr to determine which 
   # vdgs are identical.
   namingdict = {}
   for nonredun_vdg in all_nr_vdgs:
      nonredun_vdg_pdb = nonredun_vdg[4]
      pdbcode = nonredun_vdg_pdb.rstrip('/').rstrip('.pdb').split('/')[-1]
      if pdbcode not in namingdict.keys():
         namingdict[pdbcode] = []
      scrr_cg_perm = nonredun_vdg[5]
      sorted_scrr = sorted(scrr_cg_perm)
      if sorted_scrr not in namingdict[pdbcode]:
         namingdict[pdbcode].append(sorted_scrr)

   # Print out each vdg (separated by the # of vdms in the vdg)
   print('NEED TO ALIGN ON REF COORD bc of single_vdg.')
   for pdbcode, list_scrrs_cg_perm in namingdict.items():
      for vdg_scrr_cg_perm in list_scrrs_cg_perm:
         vdg_scrr_str = '_'.join(['_'.join([str(z) for z in v_s]) for v_s in 
                                  vdg_scrr_cg_perm])
         pr_obj = isolate_vdg_subset_obj(pdbcode, vdg_scrr_cg_perm, vdg_pdbs_dir)
         # Write out pdb
         write_vdg_subset(pdbcode, vdg_scrr_cg_perm, vdg_scrr_str, pr_obj, out_dir)
         
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

def isolate_vdg_subset_obj(pdbcode, vdg_scrr_cg_perm, vdg_pdbs_dir):
   # Select just the specified vdms (vdg subset) for printing out
   # Initialize prody obj for printing out
   input_vdg_pdbpath = os.path.join(vdg_pdbs_dir, pdbcode + '.pdb')
   par = pr.parsePDB(input_vdg_pdbpath)
   pr_obj = par.select('occupancy > 2.8') # occ >= 3 is CG
   for vdm_scrr in vdg_scrr_cg_perm:
      if vdm_scrr[0]: # has a pdb segment defined
         res_obj = par.select(
            f'segment {vdm_scrr[0]} and chain {vdm_scrr[1]} and resnum {vdm_scrr[2]}'
            )
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
   assert num_resinds == len(vdg_scrr_cg_perm) + 1 # +1 is because CG counts as a 
                                                   # resind
   return pr_obj

def write_vdg_subset(pdbcode, vdg_scrr_cg_perm, vdg_scrr_str, pr_obj, out_dir):
   n_vdms_in_subset = str(len(vdg_scrr_cg_perm))
   outputname = f'{n_vdms_in_subset}_{pdbcode}_{vdg_scrr_str}.pdb'
   outputpath = os.path.join(out_dir, n_vdms_in_subset, outputname)
   if str(n_vdms_in_subset) not in os.listdir(out_dir):
      os.mkdir(os.path.join(out_dir, n_vdms_in_subset))
   pr.writePDB(outputpath, pr_obj)

def get_nr_vdgs_of_same_AA_comp(_vdgs, seq_sim_thresh, reordered_AAs, 
                                symmetry_classes, output_clus_pdbs, vdglib_dir):
   # vdG subsets that have identical vdm AA compositions may be redudant.
   # First, get all AA permutations for equivalent AAs. For example, if the vdms are 
   # [Ala1, Ala2, Glu], then we need to sample [Ala1, Ala2, Glu] and [Ala2, Ala1, 
   # Glu]. The easiest way to get the permutations is to label by the index of the 
   # ordered AAs
   
   # Analyze all AA permutations of the vdgs.
   all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords, \
      all_AA_perm_flankingseqs, all_AA_perm_flankingCAs, \
      all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr = \
      clust.get_vdg_AA_permutations(reordered_AAs, _vdgs)

   # Add in permutations of symmetric CG atoms.
   all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords, \
      all_AA_cg_perm_flankingseqs, all_AA_cg_perm_flankingCAs, \
      all_AA_cg_perm_pdbpaths, all_AA_cg_perm_vdm_scrr_cg_perm = \
      clust.get_vdg_AA_and_cg_perms(all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords, 
                              all_AA_perm_flankingseqs, all_AA_perm_flankingCAs, 
                              all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr, 
                              symmetry_classes)
   
   # First, calc rmsd between all AA permutations
   all_cg_and_vdmbb_coords = clust.combine_cg_and_vdmbb_coords(
      all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords)

   # Cluster these potential redundant AA permuted vdGs by rmsd of CG + vdm bb 
   # coords. If their rmsds are larger than the rmsd thershold, then they are 
   # automatically considered different vdGs. `cluster_assignments` is a dict where 
   # each key is the cluster number, and each value is a list of the indices of the 
   # vdGs (indexing is relative position in all_AA_cg_perm_cg_coords, 
   # all_AA_cg_perm_vdm_bbcoords, etc.) that belong in that cluster.
   num_coords = len(all_cg_and_vdmbb_coords[0])
   rmsd_cut = normalize_rmsd(num_coords)
   print('Number of coords, rmsd_cut', num_coords, rmsd_cut)
   # `cgvdmbb_cluster_assignments` dict not ordered by # cluster size
   cgvdmbb_cluster_assignments = clust.get_hierarchical_clusters(
      all_cg_and_vdmbb_coords, rmsd_cut, None_in_coords=True) 

   if output_clus_pdbs:
      clusdir = os.path.join(
         vdglib_dir, 'clusters', 'temp', 'cgvdmbb', str(len(reordered_AAs)), 
         '_'.join(reordered_AAs)) # temp bc need to be reassigned
                                          # to handle degenerate binding sites
      utils.handle_existing_files(clusdir)
      clust.write_out_clusters(clusdir, cgvdmbb_cluster_assignments, 
         all_AA_cg_perm_cg_coords, all_AA_cg_perm_pdbpaths, 
         all_AA_cg_perm_vdm_scrr_cg_perm, symmetry_classes, 
         all_cg_and_vdmbb_coords)

   nr_vdgs = []
   for cgvdmbb_clusnum, indices_of_elements_in_cg_vdmbb_cluster in \
      cgvdmbb_cluster_assignments.items():
      # Within these clusters (because each of these clusters are nonredundant from 
      # each other in terms of CG and vdm bb positions), cluster again based on 
      # backbone stretches flanking the vdms.
      cgvdmbb_clus_cg_coords, cgvdmbb_clus_vdmbb_coords, cgvdmbb_clus_flankingseqs, \
         cgvdmbb_clus_flankingCAs, cgvdmbb_clus_pdbpaths, \
         cgvdmbb_clus_vdm_scrr_cg_perm = clust.elements_in_clusters(
         indices_of_elements_in_cg_vdmbb_cluster, all_AA_cg_perm_cg_coords, 
         all_AA_cg_perm_vdm_bbcoords, all_AA_cg_perm_flankingseqs, 
         all_AA_cg_perm_flankingCAs, all_AA_cg_perm_pdbpaths, 
         all_AA_cg_perm_vdm_scrr_cg_perm)
      
      # If there's only 1 element in the cluster, then it's automatically 
      # nonredundant b/c its CG+vdm bb position (binding pose) doesn't match any 
      # other vdG's binding pose
      if len(indices_of_elements_in_cg_vdmbb_cluster) <= 1:
         vdg_descript = [cgvdmbb_clus_cg_coords[0], cgvdmbb_clus_vdmbb_coords[0], 
            cgvdmbb_clus_flankingseqs[0], cgvdmbb_clus_flankingCAs[0], 
            cgvdmbb_clus_pdbpaths[0], cgvdmbb_clus_vdm_scrr_cg_perm[0]] # these 6 
                                                         # elements define a vdg
         nr_vdgs.append(vdg_descript)
         continue

      # Otherwise, if there's >1 element in the cluster, it's possible that the vdgs 
      # within this cluster are redundant. Next thing to look at is rmsd of those 
      # flanking bb stretches. Align and cluster.
      flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster = []
      for vdg_flankingCAs in cgvdmbb_clus_flankingCAs:
         # `vdg_flankingCAs` is a list of vdm residues, so flatten it
         flat_flanking_CAs = []
         for res in vdg_flankingCAs:
            for CA_coord in res:
               flat_flanking_CAs.append(CA_coord)
         flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster.append(flat_flanking_CAs)
      # Cluster the vdgs in `flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster`. In 
      # order to calculate rmsd for clustering, take care of the instances where 
      # coordinates are represented as "None" because the stretch of +/- 5 AAs was 
      # missing N-terminal or C-terminal residues. Only select residues that have CA 
      # coords in both of the vdgs whose rmsd is being calculated.
      flankingbb_cluster_assignments = clust.get_hierarchical_clusters(
         flattened_flankingCAs_for_vdgs_in_cgvdmbb_cluster, None_in_coords=True)
      for flankingCAs_clusnum, indices_of_elements_in_flankingCAs_cluster in \
         flankingbb_cluster_assignments.items():
         # For each cluster determined by rmsd between CA coords flanking the vdms, 
         # gather the vdg features (cg coords, vdmbb coords, etc.) corresponding to 
         # the vdgs belonging in that cluster. The vdgs are called based on their 
         # indices within the clusters defined in the first step, when the rmsds of 
         # cg + vdm bbs were being compared.
         flankingCAs_clus_cg_coords, flankingCAs_clus_vdmbb_coords, \
            flankingCAs_clus_flankingseqs, flankingCAs_clus_flankingCAs, \
            flankingCAs_clus_pdbpaths, flankingCAs_clus_vdm_scrr_cg_perm = \
            clust.elements_in_clusters(indices_of_elements_in_flankingCAs_cluster, 
            cgvdmbb_clus_cg_coords, cgvdmbb_clus_vdmbb_coords, 
            cgvdmbb_clus_flankingseqs, cgvdmbb_clus_flankingCAs,
            cgvdmbb_clus_pdbpaths, cgvdmbb_clus_vdm_scrr_cg_perm)

         # If there's only 1 element in the cluster, then it's automatically 
         # nonredundant b/c its binding site backbone doesn't match any other vdG's 
         # binding site
         if len(indices_of_elements_in_flankingCAs_cluster) <= 1:
            vdg_descript = [flankingCAs_clus_cg_coords[0], 
            flankingCAs_clus_vdmbb_coords[0], flankingCAs_clus_flankingseqs[0], 
            flankingCAs_clus_flankingCAs[0], flankingCAs_clus_pdbpaths[0], 
            flankingCAs_clus_vdm_scrr_cg_perm[0]] # these 6 define a vdg
            nr_vdgs.append(vdg_descript)
            continue

         # Otherwise, if there's more than 1 element in the cluster, it's possible 
         # that the vdgs within this cluster are redundant. Next metric to cluster on 
         # is sequence dissimilarity of the residues +/- of the vdms.
         flattened_flankingseqs_for_vdgs_in_flankingCA_clus = []
         
         for vdg_flankingseq in flankingCAs_clus_flankingseqs:
            flat_vdg_flankingseq = []
            for vdm_res in vdg_flankingseq:
               flat_vdg_flankingseq += vdm_res
            flattened_flankingseqs_for_vdgs_in_flankingCA_clus.append(
               flat_vdg_flankingseq)
         seqsim_clus_assignments = clust.get_hierarchical_clusters(
            flattened_flankingseqs_for_vdgs_in_flankingCA_clus, 
            seq_sim_thresh=seq_sim_thresh)
         # If the vdgs belong in the same cluster based on sequence similarity (and 
         # therefore backbones and cg+vdmbb binding pose), then the vdgs are 
         # redundant. Select only one to add to the `nr_vdgs` list.
         for seqsim_clusnum, indices_of_elements_in_flankingseq_clus in \
            seqsim_clus_assignments.items():
            flankingseq_clus_cg_coords, flankingseq_clus_vdmbb_coords, \
               flankingseq_clus_flankingseqs, flankingseq_clus_flankingCAs, \
               flankingseq_clus_pdbpaths, flankingseq_clus_vdm_scrr_cg_perm = \
               clust.elements_in_clusters(indices_of_elements_in_flankingseq_clus, \
               flankingCAs_clus_cg_coords, flankingCAs_clus_vdmbb_coords, \
               flankingCAs_clus_flankingseqs, flankingCAs_clus_flankingCAs, \
               flankingCAs_clus_pdbpaths, flankingCAs_clus_vdm_scrr_cg_perm)
           
            if len(indices_of_elements_in_flankingseq_clus) <= 1:
               vdg_descript = [flankingseq_clus_cg_coords[0], 
               flankingseq_clus_vdmbb_coords[0], flankingseq_clus_flankingseqs[0],
               flankingseq_clus_flankingCAs[0], flankingseq_clus_pdbpaths[0], 
               flankingseq_clus_vdm_scrr_cg_perm[0]]
               nr_vdgs.append(vdg_descript)
            
            else: # select only one, so use index 0 for 
                  # indices_of_elements_in_flankingseq_clus
               flankingseq_clus_cg_coords, flankingseq_clus_vdmbb_coords, \
               flankingseq_clus_flankingseqs, flankingseq_clus_flankingCAs, \
               flankingseq_clus_pdbpaths, flankingseq_clus_vdm_scrr_cg_perm = \
               clust.elements_in_clusters([0], \
               flankingCAs_clus_cg_coords, flankingCAs_clus_vdmbb_coords, \
               flankingCAs_clus_flankingseqs, flankingCAs_clus_flankingCAs, \
               flankingCAs_clus_pdbpaths, flankingCAs_clus_vdm_scrr_cg_perm)
               vdg_descript = [flankingseq_clus_cg_coords[0], 
                  flankingseq_clus_vdmbb_coords[0], flankingseq_clus_flankingseqs[0],
                  flankingseq_clus_flankingCAs[0], flankingseq_clus_pdbpaths[0], 
                  flankingseq_clus_vdm_scrr_cg_perm[0]]
               nr_vdgs.append(vdg_descript)
              
   return nr_vdgs 

def normalize_rmsd(num_atoms):
    # Smoothly scale rmsd threshold by the number of atoms.
    min_atoms = 8
    min_threshold = 0.3
    max_atoms = 20
    max_threshold = 1.3 
    
    if num_atoms < min_atoms:
        return min_threshold  # below 5 atoms, use the minimum threshold (1.5 Å)
    if num_atoms > max_atoms:
        return max_threshold  # above 20 atoms, use the maximum threshold (2.5 Å)
    
    # Linear scaling for atoms between 5 and 20
    scaling_factor = (num_atoms - min_atoms) / (max_atoms - min_atoms)
    threshold = min_threshold + scaling_factor * (max_threshold - min_threshold)
    
    return threshold

if __name__ == "__main__":
    main()
