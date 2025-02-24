'''
Removes redundant vdGs from the vdG library. vdGs are redundant if they meet all 
criteria: 
   1. Sequence similarity of residues flanking the vdMs. 
         Note that the clustering criterion is actually dissimilarity, because the 
         "distances" between cluster members must be *below* the threshold.  
   2. vdM AA identities and binding pose (CG + vdM backbone atoms within a 
      specified RMSD threshold). 
         Each AA permutation (for multiple vdMs of the same AA) is listed as a 
         distinct vdG. Each CG permutation (for symmetrically equivalent CG atoms) 
         is also listed as a distinct vdG. E.g., carboxylates CGs have 2 permutations 
         because of 2 identical oxygens. After clustering, duplicate vdGs (different 
         AA and CG permutations of the same vdG) are removed.
   3. Similar positions of the residues flanking the vdMs (rmsd of their CAs within 
      a specified threshold).
         The span of contiguous residues flanking the vdMs is defined by the 
         --flank argument; see --help.

First, vdGs are clustered on rmsd of cg + vdm bb. Then, secondary clustering occurs on 
a combined dist metric of flank bb seq dissimilarity + flank bb rmsd.
'''

import os
import sys
import argparse
import time
import shutil
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
    #parser.add_argument('-o', "--output-clus-pdbs", action='store_true', 
    #                    help="Output clustered PDBs.")
    parser.add_argument('-w', "--align-cg-weight", type=float, default=0.90, 
                        help="Fraction of weights to assign to CG atoms (collectively) "
                        "when superposing output vdGs. Not weights for clustering. "
                        "Example: 0.5 means 1/2 of weight is assigned to CG atoms and "
                        "the remaining 1/2 goes to the vdM backbone atoms.")
    parser.add_argument('-q', "--seq", default=0.40,
                        help="Sequence similarity threshold for clustering sequences "
                        "flanking the vdM to determine redundancy. This should be a "
                        "value between 0 and 1. Values > seq will be considered "
                        "redundant.")
    parser.add_argument('-f', "--flank", default=2,
                        help="Number of residues +/- the vdms for which to "
                        "calculate sequence similarity and backbone similarity.")
    parser.add_argument('-s', '--symmetry-classes', nargs='+', type=int,
                        help='Integers representing the symmetry classes of the CG '
                        'atoms on which clustering is to be performed. If provided, '
                        'should have the same length as idxs. If not provided, the '
                        'atoms are assumed to be symmetrically inequivalent.')
    parser.add_argument('-n', '--size-subset', type=int, 
                        help='Specify the number of residues in the vdG subset. '
                        'The purpose of this arg is for parallelization (running '
                        'this script with different -n values concurrently. 1, 2, '
                        '3, and 4 are recommended.')
    parser.add_argument('-l', "--logfile", help="Path to log file.")
    parser.add_argument('-p', "--print-flankbb", action='store_true', 
                        help="Include flanking bb residues when writing out PDB.")
    
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
   print_flankbb = args.print_flankbb
   size_subset = args.size_subset
   if size_subset not in [1, 2, 3, 4]:
      raise ValueError('size_subset must be an integer between 1 and 4.')
   logfile = args.logfile
   align_cg_weight = args.align_cg_weight
   if align_cg_weight < 0 or align_cg_weight > 1:
      raise ValueError('align_cg_weight must be a float between 0 and 1.')
   vdglib_dir = args.vdglib_dir
   vdg_pdbs_dir = os.path.join(vdglib_dir, 'vdg_pdbs')
   out_dir = os.path.join(vdglib_dir, 'nr_vdgs')

   with open(logfile, 'a') as file:
        file.write(f"{'='*20} Starting deduplicate_reun_vdgs.py run {'='*20} \n")
   
   # Initialize vdm_combos dict to store the subsets of vdMs within a vdG
   vdm_combos = {}

   # Iterate over the PDBs and CGs that were identified as containing the SMARTS group
   atomgroup_dict = {}
   vdg_pdbs_in_dir = os.listdir(vdg_pdbs_dir)
   for pdbname in vdg_pdbs_in_dir:
      pdbpath = os.path.join(vdg_pdbs_dir, pdbname)
      prody_obj = pr.parsePDB(pdbpath)
      if len(vdg_pdbs_in_dir) < 10000:
         atomgroup_dict[os.path.join(vdg_pdbs_dir, pdbname)] = prody_obj
      cg_coords = clust.get_cg_coords(prody_obj)
      if cg_coords is None:
         logfile.write(f'\tIncorrect number of occ >= 3 atoms in vdg_pdbs/{pdbname}.\n')
         continue
      # define symmetry class if it's None so it's compatible with pdb output naming
      if symmetry_classes is None:
         symmetry_classes = [i for i in range(len(cg_coords))]
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
   
   
   # Evaluate the complete collection of vdGs (of size_subset) and determine redundancy 
   for num_vdms_in_subset, _subsets in vdm_combos.items():
      if num_vdms_in_subset != size_subset:
         continue
      for _reordered_AAs, _vdgs in _subsets.items():
         # vdG subsets that have identical vdm AA compositions may be redundant.
         cluster_vdgs_of_same_AA_comp(_vdgs, seq_sim_thresh, _reordered_AAs, 
               symmetry_classes, vdglib_dir, align_cg_weight, num_flanking, 
               atomgroup_dict, print_flankbb, logfile)
         '''Final results: vdgs that end up in the same cg+vdmbb, flanking bb, and 
         flanking seq clusters are redundant. Select only the centroid to be output.'''
         copy_nr_to_outdir(vdglib_dir, out_dir, _reordered_AAs)

   '''
   Clean up the entire output from clus_and_deduplicate_vdgs.py.
   At each clustering stage (flankseq, flankbb, then cgvdmbb), the clustered PDBs are 
   output so that they can be loaded in for the next step. (The libraries 
   are usually way too large to store in memory.) At the end of clustering (at the cgvdmbb
   level), the clustered output for cgvdmbb is organized by: 
   CG/clusters/cgvdmbb/cgvdmbb_clusnum/AAs/flankseq_clusnum/flankbb_clusnum/cgvdmbb_clusnum.
   Therefore, it's easy to trace the flankseq, flankbb, and cgvdmbb cluster numbers, and
   earlier output (clusters/flankseq/ and clusters/flankbb) can be deleted, as well as
   the temp dir used for the different permutations of cgvdmbb.
   '''
   delete_clusterdirs(vdglib_dir, logfile, size_subset)

   # Print out time elapsed
   seconds = time.time() - start_time
   hours = round(seconds // 3600)
   minutes = (seconds % 3600) // 60
   seconds = seconds % 60
   seconds = round(seconds, 2)
    
   # Write out results to log file
   with open(logfile, 'a') as file:
      num_vdg_pdbs = len(vdg_pdbs_in_dir)
      nr_dir_of_size_subset = os.path.join(out_dir, str(size_subset))
      num_vdgs_of_size_subset = 0
      if not os.path.exists(nr_dir_of_size_subset):
         pass
      else:
         for aa_subset in os.listdir(nr_dir_of_size_subset):
            aa_subset_len = len(os.listdir(os.path.join(nr_dir_of_size_subset, aa_subset)))
            num_vdgs_of_size_subset += aa_subset_len
      
      file.write(f'\t{num_vdgs_of_size_subset} nonredun. vdgs of subset size '
                 f'{size_subset} out of {num_vdg_pdbs} vdgs.\n')
      file.write(f"Completed clus_and_deduplicate_vdgs.py in {hours} h, ")
      file.write(f"{minutes} mins, and {seconds} secs \n") 

def copy_nr_to_outdir(vdglib_dir, nr_dir, reordered_AAs):
   clustersdir = os.path.join(vdglib_dir, 'clusters', 'cgvdmbb', 
                              str(len(reordered_AAs)), '_'.join(reordered_AAs))
   nr_dir = os.path.join(nr_dir, str(len(reordered_AAs)), '_'.join(reordered_AAs))
   utils.handle_existing_files(nr_dir)
   os.makedirs(nr_dir, exist_ok=True)
   for cgvdmbbclus in os.listdir(clustersdir):
      cgvdmbbclusdir = os.path.join(clustersdir, cgvdmbbclus)
      for pdb in os.listdir(cgvdmbbclusdir):
         if 'centroid' not in pdb: 
            continue
         else:
            # copy the centroid to the output dir
            biounit = '_'.join(pdb.split('_')[:5])
            vdmscrr_list = pdb.split('_')[5:-1]
            # remove resnames for conciseness
            vdmscrr_list = [vdmscrr_list[i] for i in range(len(vdmscrr_list)) 
                            if (i + 1) % 4 != 0] 
            vdmscrr_str = '_'.join(vdmscrr_list)
            newname = f'{biounit}_{vdmscrr_str}.pdb.gz'
            newpath = os.path.join(nr_dir, newname)
            
            assert not os.path.exists(newpath)
            shutil.copy(os.path.join(cgvdmbbclusdir, pdb), newpath) 
            break 

def cluster_vdgs_of_same_AA_comp(_vdgs, seq_sim_thresh, reordered_AAs, 
      symmetry_classes, vdglib_dir, align_cg_weight, num_flanking, 
      atomgroup_dict, print_flankbb, logfile):
   # vdG subsets that have identical vdm AA compositions may be redudant.
   # First, get all AA permutations for identical AAs. For example, if the vdms are 
   # [Ala1, Ala2, Glu], then we need to sample [Ala1, Ala2, Glu] and [Ala2, Ala1, 
   # Glu]. The easiest way to get the permutations is to label the index of the 
   # ordered AAs
   
   # Include all AA permutations (for identical AAs) of the vdgs.
   all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords, \
      all_AA_perm_flankingseqs, all_AA_perm_flankingCAs, \
      all_AA_perm_pdbpaths, all_AA_perm_vdm_scrr = \
      clust.get_vdg_AA_permutations(reordered_AAs, _vdgs)

   # Add in permutations of symmetric CG atoms.
   all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords, \
      all_AA_cg_perm_flankingseqs, all_AA_cg_perm_flankingCAs, \
      all_AA_cg_perm_pdbpaths, all_AA_cg_perm_vdm_scrr_cg_perm = \
      clust.get_vdg_AA_and_cg_perms(all_AA_perm_cg_coords, all_AA_perm_vdm_bbcoords, 
         all_AA_perm_flankingseqs, all_AA_perm_flankingCAs, all_AA_perm_pdbpaths, 
         all_AA_perm_vdm_scrr, symmetry_classes)
   
   all_AA_cg_perm_cg_and_vdmbb_coords = clust.combine_cg_and_vdmbb_coords(
      all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords)
   all_AA_cg_perm_flat_flankingseqs = clust.flatten_flanking_seqs(
      all_AA_cg_perm_flankingseqs)
   all_AA_cg_perm_flat_flankCAs = clust.flatten_flanking_CAs(
      all_AA_cg_perm_flankingCAs)
   # Count # of cg vdm bb atoms
   num_cg_atoms = len(all_AA_cg_perm_cg_coords[0])
   num_vdmbb_res = len(all_AA_cg_perm_vdm_bbcoords[0]) 
   num_vdmbb_atoms = num_vdmbb_res * 3 
   num_cgvdmbb_atoms = num_cg_atoms + num_vdmbb_atoms

   # First, cluster on rmsd of the cg and vdmbb atoms.
   # `cluster_assignments` dict: each key is a clus number, and each value
   # is a list of the indices belong to that cluster. Indexing is the relative position
   # in all_AA_cg_perm_cg_coords, all_AA_cg_perm_vdm_bbcoords, etc.).
   cgvdmbb_rmsd_cut = normalize_rmsd(num_cgvdmbb_atoms, 'cgvdmbb') 
   cgvdmbb_data_to_clus = zip([all_AA_cg_perm_cg_and_vdmbb_coords], ['cgvdmbb'])
   cgvdmbb_clus_assignments, cgvdmbb_clus_centroids = clust.get_hierarchical_clusters(
      cgvdmbb_data_to_clus, cgvdmbb_rmsd_cut, rmsd_max=1)

   # Output the cgvdmbb clusters
     
   # Define dir to output "temp" files (temp bc clusters need to be reassigned
   # later to handle degenerate binding sites for different AA and CG 
   # permutations of the same PDB).

   temp_cgvdmbb_clusdir = os.path.join(vdglib_dir, 'clusters', 'temp', 'cgvdmbb', 
      str(len(reordered_AAs)), '_'.join(reordered_AAs)) 
   utils.handle_existing_files(temp_cgvdmbb_clusdir)
   cgvdmbb_weights = clust.get_weights(num_cg_atoms, num_vdmbb_atoms, 
                                             align_cg_weight)
   
   first_pdb_out = None # update with pdb name of first pdb that's output for ref
   first_pdb_cg_vdmbb_coords = False 

   first_pdb_out, first_pdb_cg_vdmbb_coords, failed_pdbs = clust.write_out_clusters(
      temp_cgvdmbb_clusdir, cgvdmbb_clus_assignments, cgvdmbb_clus_centroids,
      all_AA_cg_perm_cg_coords, all_AA_cg_perm_pdbpaths, 
      all_AA_cg_perm_vdm_scrr_cg_perm, all_AA_cg_perm_cg_and_vdmbb_coords, 
      all_AA_cg_perm_flat_flankCAs, num_flanking, first_pdb_out, 
      first_pdb_cg_vdmbb_coords, cgvdmbb_weights, atomgroup_dict, print_flankbb)

   if len(failed_pdbs) > 0:
      with open(logfile, 'a') as file:
         for fail in failed_pdbs:
            file.write(f'\tFailed to parse {fail}.\n')

   # Clean up the cgvdmbb output dir by merging degenerate vdGs (diff AA and 
   # CG perms of the same PDB), deleting degenerate pdbs/clusters, and 
   # reassigning clus nums by size.
   tempdir = os.path.join(vdglib_dir, 'clusters', 'temp', 'cgvdmbb')
   cleandir = os.path.join(vdglib_dir, 'clusters', 'cgvdmbb')
   reassigned_cgvdmbb_clus = clust.rewrite_temp_clusters(tempdir, cleandir, 
                                    len(reordered_AAs), '_'.join(reordered_AAs))

   # Next, further subdivide each cluster on rmsd of the bb stretches flanking the vdms 
   # + their sequences. 
   for cgvdmbb_clusnum, indices_of_elements_in_cgvdmbb_cluster in \
      cgvdmbb_clus_assignments.items():
      # Select the data belonging to this cgvdmbb cluster.
      (cgvdmbb_clus_cg_coords, cgvdmbb_clus_vdmbb_coords, 
         cgvdmbb_clus_cgvdmbb_coords, cgvdmbb_clus_flat_flankingseqs, 
         cgvdmbb_clus_flat_flankingCAs, cgvdmbb_clus_pdbpaths, 
         cgvdmbb_clus_vdm_scrr_cg_perm) = clust.elements_in_clusters(
         indices_of_elements_in_cgvdmbb_cluster, all_AA_cg_perm_cg_coords, 
         all_AA_cg_perm_vdm_bbcoords, all_AA_cg_perm_cg_and_vdmbb_coords,
         all_AA_cg_perm_flat_flankingseqs, all_AA_cg_perm_flat_flankCAs, 
         all_AA_cg_perm_pdbpaths, 
         all_AA_cg_perm_vdm_scrr_cg_perm)   

      flankseq_and_bb_to_clus = zip([cgvdmbb_clus_flat_flankingseqs,
                        cgvdmbb_clus_flat_flankingCAs], ['flankseq', 'flankbb'])
      num_flank_bb_coords = len(cgvdmbb_clus_flat_flankingCAs[0])
      flankbb_rmsd_cut = normalize_rmsd(num_flank_bb_coords, 'flankbb')
      # Combine rmsd and seq into a single metric and threshold
      flankseq_and_bb_thresh = flankbb_rmsd_cut + (1 - seq_sim_thresh)
      flankingseq_and_bb_cluster_assignments, flankingseq_and_bb_clus_centroids \
         = clust.get_hierarchical_clusters(flankseq_and_bb_to_clus, flankseq_and_bb_thresh, 
                                           rmsd_max=1.5)
            
      # Output the flankseq+bb clusters
      flankseq_and_bb_clusdir_for_this_cgvdmbb_clus = os.path.join(vdglib_dir, 'clusters', 
         'flankseq_and_bb', str(len(reordered_AAs)), '_'.join(reordered_AAs), 
         f'flankseqclus_{cgvdmbb_clusnum}') 
      utils.handle_existing_files(flankseq_and_bb_clusdir_for_this_cgvdmbb_clus) 

      first_pdb_out, first_pdb_cg_vdmbb_coords, failed_pdbs = clust.write_out_clusters(
         flankseq_and_bb_clusdir_for_this_cgvdmbb_clus, 
         flankingseq_and_bb_cluster_assignments, flankingseq_and_bb_clus_centroids,
         cgvdmbb_clus_cg_coords, cgvdmbb_clus_pdbpaths, cgvdmbb_clus_vdm_scrr_cg_perm, 
         cgvdmbb_clus_cgvdmbb_coords, cgvdmbb_clus_flat_flankingCAs, num_flanking, 
         first_pdb_out, first_pdb_cg_vdmbb_coords, None, atomgroup_dict, print_flankbb)

      if len(failed_pdbs) > 0:
         with open(logfile, 'a') as file:
            for fail in failed_pdbs:
               file.write(f'\tFailed to parse {fail}.\n')

def delete_clusterdirs(vdglib_dir, logfile, size_subset):
   with open(logfile, 'a') as file:
      all_clusters = os.path.join(vdglib_dir, 'clusters')
      # delete temp folder
      tempdir = os.path.join(all_clusters, 'temp', 'cgvdmbb', str(size_subset))
      if not os.path.exists(tempdir):
         file.write(f"\t{tempdir} does not exist.\n")
      else:
         shutil.rmtree(tempdir)
      # delete cgvdmbb 
      cgvdmbbdir = os.path.join(all_clusters, 'cgvdmbb', str(size_subset))
      if not os.path.exists(cgvdmbbdir):
         file.write(f"\t{cgvdmbbdir} does not exist.\n")
      else:
         shutil.rmtree(cgvdmbbdir)

def normalize_rmsd(num_atoms, atoms):
   if atoms == 'flankbb':
      max_threshold = 1.5
      min_threshold = 0.5
   elif atoms == 'cgvdmbb':
      max_threshold = 1
      min_threshold = 0.5

   # Smoothly scale rmsd threshold by the number of atoms.
   min_atoms = 8
   max_atoms = 15
   
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
