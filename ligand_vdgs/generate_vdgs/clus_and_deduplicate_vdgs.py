'''
Removes redundant vdGs from the vdG library.
A vdG is considered redundant if it meets both criteria after clustering:

1) Binding pose equivalence
   vdGs are clustered by RMSD of [CG + vdM backbone atoms] using a size-normalized cutoff.
   Members belong to the same cluster when their distance ≤ cutoff.

2) Local environment similarity
   Within each pose cluster, vdGs are subclustered on a composite distance:
       `distance = flank-CA RMSD + (sequence dissimilarity)/2`
   with threshold:
       `threshold = size-normalized flank-CA RMSD cutoff + (1 − seq_sim_thresh)/2`, 
   where sequence dissimilarity = 1 − (percent identity).  

---
Workflow Overview

1. **Scan & stream**  
   Parse each PDB, extract vdGs, and write AA-combo buckets to disk.

2. **Primary clustering (pose)**  
   Within each AA bucket, cluster on RMSD of [CG + vdM backbone] using a size-scaled cutoff.

3. **Secondary subclustering (local environment)**  
   For each pose cluster, subcluster on  
       `distance = flank-CA RMSD + (sequence dissimilarity)/2`  
   using the threshold defined above.  The flank window is controlled by `--flank`.

4. **De-duplication & output**  
   Merge AA/CG-permutation duplicates, retain cluster centroids, and output resulting PDBs.
   Optionally keep intermediate clustered PDBs.

5. **Cleanup**  
   Remove temporary stream and clustering directories.
'''

import os
import sys
import argparse
import time
import tracemalloc
import shutil
import numpy as np
from numba import njit
import prody as pr
import multiprocessing as mp
import json, gzip, hashlib
import traceback
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import align_and_cluster as clust
import utils

def parse_args():
   parser = argparse.ArgumentParser()
   parser.add_argument('-c', "--cg", type=str, 
                       help="SMARTS pattern of the fragment/CG/FG. Informational only.")
   parser.add_argument('-v', "--vdglib-dir", type=str, required=True,
                       help="Directory for the vdms of this CG.")
   parser.add_argument('-k', "--keep-clustered-pdbs", action='store_true', 
                       help="Keep the `clusters/flankseq_and_bb/` dir. to keep track of "
                       "which clusters each nr_pdb came from.")
   parser.add_argument('-w', "--align-cg-weight", type=float, default=0.99, 
                       help="Weight assigned to CG atoms when superposing output vdGs "
                       "(not for clustering). Example: 0.5 assigns half of the weight " 
                       "to CG atoms and half to vdM backbone atoms. Defaults to 0.99.")
   parser.add_argument('-q', "--seq", default=0.40, type=float, 
                       help="Sequence similarity threshold (fraction between 0 and 1) "
                       "used to set the subclustering threshold. Higher similarity "
                       "indicates redundancy. Defaults to 0.40.")
   parser.add_argument('-f', "--flank", default=2, type=int, 
                       help="Number of residues on each side (+/-) of each vdM to " 
                       "include for computing sequence and backbone similarity. "
                       "Defaults to 2.")
   parser.add_argument('-s', "--symmetry-classes", nargs='+', type=int,
                       help="Integers defining CG symmetry. Ex: carboxylate CC(=O)[O-] "
                       "would be '0 1 2 2' b/c of its equivalent oxygens. If this arg " 
                       "is not provided, the CG is assumed to be asymmetric.")
   parser.add_argument('-n', "--size-subset", type=int, required=True,  
                       help="The number of residues in the vdG subset. The purpose "
                       "of this arg is for parallelization (running with different "
                       "-n values concurrently).")
   parser.add_argument('-l', "--logfile", type=str, default='log', help="Path to log.")
   parser.add_argument('-x', "--print-flankbb", action='store_true', 
                       help="Include flanking bb residues when writing out PDB.")
   parser.add_argument('-m', "--max-num-vdgs-to-clus", default=2500, type=int, 
                       help="Cap on # samples per AA-composition bucket to cluster.")
   parser.add_argument('--num-procs', type=int, default=10,
                       help="Number of AA composition buckets to run concurrently.")
   return parser.parse_args()

'''
Characterize all subsets of vdMs within each vdG. First, scan the parent dataset to 
identify and record all vdGs with the SMARTS pattern.  Each vdG record is immediately 
written ("spilled") to a temporary, per–AA–composition file on disk rather than stored 
in memory.  After the scan completes, those per–AA–bucket files are processed 
independently in parallel: each bucket is loaded, clustered by CG alignment and sequence 
similarity, and de-duplicated.

Contents of bucket JSONL (one vdG per line): 
{
  "cg_coords": [[x,y,z], ...],
  "bbcoords":  [ [[Nx,Ny,Nz],[CAx,CAy,CAz],[Cx,Cy,Cz]],  ... ],
  "flankseqs": [ ["-2AA","-1AA","vdm","+1AA","+2AA"], ... ],
  "flankCAs":  [ [[x,y,z],...], ... ],
  "pdbpath":   "vdg pdb path",
  "scrr":      [ [seg, chain, resnum, resname], ... ]
}

Each record contains several per-vdM lists that are kept strictly index-aligned: (AA 
identities, vdM backbone coords, flanking sequences, flanking CA coords, and vdM
identifiers (seg, chain, resnum, resname)). The i-th element across these lists refers 
to the same vdM residue. When we reorder vdMs (e.g., by AA name) or apply AA/CG 
permutations, we apply the same permutation to all lists, so downstream clustering can 
safely match elements by index.
'''

def _run_one_bucket_strict(args):
   # Run one AA-composition bucket per process; let any exceptions propagate back to 
   # parent to hard-fail and abort the entire run instead of being swallowed by the worker

   (bucket_path, seq_sim_thresh, _reordered_AAs, symmetry_classes,
     vdglib_dir, align_cg_weight, num_flanking, logfile, print_flankbb,
     max_num_to_clus) = args

   _vdgs = _load_bucket(bucket_path) 

   # Cluster vdG subsets with identical vdM AA compositions to identify redundancy. 
   # If there are > max_num_to_clus samples, select only the PDBs with highest PDB ID
   # diversity. Note that a single PDB can contribute multiple vdGs, so even after  
   # applying the max_num_to_clus cap (which limits PDBs, not vdGs), the resulting num. 
   # of vdGs may still exceed max_num_to_clus.
   if len(_vdgs) > max_num_to_clus:
      orig_num_vdgs = len(_vdgs)
      pdbIDs = [z[-2].split('/')[-1][:4] for z in _vdgs]
      diverse_pdbIDs = select_diverse_pdbIDs(pdbIDs, max_num_to_clus)
      _vdgs = [z for z in _vdgs if z[-2].split('/')[-1][:4] in diverse_pdbIDs]
      with open(logfile, 'a') as file:
         file.write(f'\tThere are {orig_num_vdgs} vdgs for the {tuple(_reordered_AAs)} '
            f'subset, so only {len(_vdgs)} vdgs with maximum PDB ID '
            f'diversity were selected.\n')
      # Note that len(_vdgs) may be > max_num_to_clus because >1 vdgs might 
      # belong to the same PDB.

   try:
      atomgroup_dict = {} # {pdb_path: prody.AtomGroup}

      # Cluster + write out for this AA bucket
      cluster_vdgs_of_same_AA_comp(
         _vdgs, seq_sim_thresh, _reordered_AAs, symmetry_classes, vdglib_dir,
         align_cg_weight, num_flanking, atomgroup_dict,
         print_flankbb=print_flankbb, logfile=logfile)

      # Copy centroids to nr output dir for this bucket
      copy_nr_to_outdir(vdglib_dir, os.path.join(vdglib_dir, 'nr_vdgs'), _reordered_AAs)
      return ('_'.join(_reordered_AAs), len(_vdgs))

   except Exception as e:
      aa_label = '_'.join(_reordered_AAs) if _reordered_AAs else 'UNKNOWN_AAs'
      tb = traceback.format_exc()
      raise RuntimeError(
         f"[WORKER ERROR] AA bucket: {aa_label} | count={len(_vdgs)}\n"
         f"Exception: {e}\nTraceback:\n{tb}")

def _aa_key_str(reordered_aas):
   return "_".join(reordered_aas)

def _aa_tmp_dir(vdglib_dir, size_subset):
   # Write per-AA-bucket records to disk to free memory
   return os.path.join(vdglib_dir, "stream_tmp", str(size_subset))

def _aa_tmp_path(vdglib_dir, size_subset, aa_key):
   h = hashlib.sha1(aa_key.encode()).hexdigest()[:16]
   fname = f"{aa_key[:80]}__{h}.jsonl.gz"
   return os.path.join(_aa_tmp_dir(vdglib_dir, size_subset), fname)

def _spill_record(vdglib_dir, size_subset, reordered_aas, record):
   # Append one vdG record to its AA bucket file (gzip JSONL)
   aa_key = _aa_key_str(reordered_aas)
   path = _aa_tmp_path(vdglib_dir, size_subset, aa_key)
   os.makedirs(os.path.dirname(path), exist_ok=True)
   with gzip.open(path, "at", encoding="utf-8") as f:
      f.write(json.dumps(record) + "\n")

def _load_bucket(path):
   # Load a whole AA bucket
   _vdgs = []
   with gzip.open(path, "rt", encoding="utf-8") as f:
      for line in f:
         rec = json.loads(line)
         _vdgs.append([
            np.asarray(rec["cg_coords"]),
            [np.asarray(r) for r in rec["bbcoords"]],
            rec["flankseqs"],
            [np.asarray(r) for r in rec["flankCAs"]],
            rec["pdbpath"],
            rec["scrr"],
         ])
   return _vdgs

def main():
   start_time = time.time()
   tracemalloc.start()
   args = parse_args()
   CG = args.cg
   seq_sim_thresh = args.seq
   num_flanking = args.flank
   symmetry_classes = args.symmetry_classes
   print_flankbb = args.print_flankbb
   keep_clustered_pdbs = args.keep_clustered_pdbs
   size_subset = args.size_subset
   max_num_to_clus = args.max_num_vdgs_to_clus
   logfile = args.logfile
   align_cg_weight = args.align_cg_weight
   if align_cg_weight < 0 or align_cg_weight > 1:
      raise ValueError('align_cg_weight must be a float between 0 and 1.')
   vdglib_dir = args.vdglib_dir
   vdg_pdbs_dir = os.path.join(vdglib_dir, 'vdg_pdbs')
   if not os.path.isdir(vdg_pdbs_dir):
      sys.exit(f"[ERROR] Missing directory: {vdg_pdbs_dir}")
   out_dir = os.path.join(vdglib_dir, 'nr_vdgs')
   
   # Iterate over the PDBs and CGs that were identified as containing the SMARTS group
   # Spill each AA bucket to disk instead of holding all in memory.
   vdg_pdbs_in_dir = sorted(os.listdir(vdg_pdbs_dir))
   for pdbname in vdg_pdbs_in_dir:
      pdbpath = os.path.join(vdg_pdbs_dir, pdbname)
      try:
         prody_obj = pr.parsePDB(pdbpath)
      except:
         with open(logfile, 'a') as file:
              file.write(f"Could not parse {pdbpath} \n")
         continue
      cg_coords = clust.get_cg_coords(prody_obj, pdbpath)
      if cg_coords is None:
         with open(logfile, 'a') as file:
            file.write(f'\tIncorrect number of occ >= 3 atoms in vdg_pdbs/{pdbname}.\n')
         continue
      # Define symmetry class if it's None so it's compatible with pdb output naming
      if symmetry_classes is None: 
         symmetry_classes = [i for i in range(len(cg_coords))] # global mutation is intentional
      # Characterize the vdM residues (bb coords, flanking residues, pdb paths, etc.)
      vdms_dict = clust.get_vdm_res_features(prody_obj, pdbpath, num_flanking)
      # Determine the vdM combinations, up to 4 residues
      vdm_resinds = list(vdms_dict.keys())
      vdg_subsets = clust.get_vdg_subsets(vdm_resinds, size_subset)
      # Iterate over subsets of resinds
      for vdg_subset in vdg_subsets:
         # Record these features in the same order as in vdg_subset. Then, sort all 
         # based on alphabetical order of the vdm AAs. If you find a bb-only vdm, 
         # assign the vdm resname as "bb".
         re_ordered_aas, re_ordered_bbcoords, re_ordered_flankingseqs, \
            re_ordered_CAs, re_ordered_scrr = clust.reorder_vdg_subset(
            vdg_subset, vdms_dict, prody_obj.select('occupancy > 2.9 and not element H'), 
            prody_obj) # exclude H b/c calculating dist to AA to determine bb or sc vdm
         
         # AA bucket stream: spill a single vdG record to its AA-bucket file to disk
         rec = {
            "cg_coords": np.asarray(cg_coords).tolist(),
            "bbcoords":  [np.asarray(x).tolist() for x in re_ordered_bbcoords],
            "flankseqs": re_ordered_flankingseqs,  # strings already fine
            "flankCAs":  [np.asarray(x).tolist() for x in re_ordered_CAs],
            "pdbpath":   str(pdbpath),
            "scrr":      [[str(s), str(ch), int(r), str(rn)] for (s, ch, r, rn) in 
                          re_ordered_scrr],
         }

         _spill_record(vdglib_dir, size_subset, re_ordered_aas, rec)

   # Evaluate the collection of vdGs of size_subset and determine redundancy. 
   # Each per-AA bucket is processed independently in parallel.
   stream_dir = _aa_tmp_dir(vdglib_dir, size_subset)
   if not os.path.isdir(stream_dir):
      with open(logfile, 'a') as f:
         f.write(f"\t[STREAM ERROR] No buckets found at {stream_dir}\n")
   else:
      jobs = []
      for fname in sorted(os.listdir(stream_dir)):
         if not fname.endswith(".jsonl.gz"):
            continue
         # recover AA key (the part before "__")
         aa_key = fname.split("__")[0]
         _reordered_AAs = aa_key.split("_")
         path = os.path.join(stream_dir, fname)

         jobs.append((
            path, seq_sim_thresh, list(_reordered_AAs), symmetry_classes,
            vdglib_dir, align_cg_weight, num_flanking, logfile,
            print_flankbb, max_num_to_clus))

      # Clean, isolated processes
      ctx = mp.get_context("spawn")
      try:
         with ctx.Pool(processes=int(args.num_procs), maxtasksperchild=1) as pool:
            # Do not swallow exceptions - any worker error should abort run
            for _ in pool.imap_unordered(_run_one_bucket_strict, jobs, chunksize=1):
               pass
      except Exception as e:
         err_text = f'\t[FATAL ERROR] Worker failed with exception:\n{e}\n'
         with open(logfile, 'a') as f:
            f.write(err_text)
         print(err_text, file=sys.stderr) # ensure it shows up in stdout/stderr
         sys.exit(1) # hard-fail the whole run if any worker fails

   # Remove per-bucket temp files that were used now that clustering is done
   try:
      stream_root = os.path.join(vdglib_dir, "stream_tmp")
      if os.path.isdir(stream_root):
         shutil.rmtree(stream_root)
   except Exception as _e:
      with open(logfile, 'a') as f:
         f.write(f'\t[STREAM WARNING]: failed to remove {stream_root}: {_e}\n')
   
   '''
   Clean up the entire output from clus_and_deduplicate_vdgs.py.
   At each clustering stage (cgvdmbb then flankseq_and_bb), the clustered PDBs are 
   output so that they can be loaded in for the next step. (The libraries 
   are usually way too large to store in memory.) At the end of the final stage of 
   clustering (the flankseq_and_bb level), the clustered output is organized by:
   CG/clusters/flankseq_and_bb/size_subset/vdg_AAs/cgvdmbb_clusnum/flankseq_and_bb_clusnum.
   Therefore, it's easy to trace the earlier cluster results, and that earlier output
   (clusters/temp/ and clusters/cgvdmbb) can be deleted. By default, the final stage of
   clustering (clusters/flankseq_and_bb) is deleted as well, because all of the clusters
   are quite memory-intensive. However, this final stage of clustering can be preserved 
   using the --keep-clustered-pdbs flag.
   '''
   delete_clusterdirs(vdglib_dir, logfile, size_subset, keep_clustered_pdbs)

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
   
   current, peak = tracemalloc.get_traced_memory()
   peak_mem = np.round(peak / (1024 * 1024 * 1024), 2)
   if peak_mem > 20:
      print(f"Peak memory usage at end of script for clustering subset size {size_subset}: "
            f"{peak_mem} GB")
   tracemalloc.stop()

def copy_nr_to_outdir(vdglib_dir, nr_dir, reordered_AAs):
   clustersdir = os.path.join(vdglib_dir, 'clusters', 'flankseq_and_bb', 
                              str(len(reordered_AAs)), '_'.join(reordered_AAs))
   nr_dir = os.path.join(nr_dir, str(len(reordered_AAs)), '_'.join(reordered_AAs))
   utils.handle_existing_files(nr_dir)
   os.makedirs(nr_dir, exist_ok=True)
   for cgvdmbbclus in os.listdir(clustersdir):
      cgvdmbbclusdir = os.path.join(clustersdir, cgvdmbbclus)
      for flankseq_and_bb_dir in os.listdir(cgvdmbbclusdir):
         flankseqandbbclusdir = os.path.join(cgvdmbbclusdir, flankseq_and_bb_dir)
         # make sure there's only 1 pdb in flankseqandbbclusdir that has 'centroid' in 
         # the name
         _centroid_pdbs = [pdb for pdb in os.listdir(flankseqandbbclusdir) 
                          if 'centroid' in pdb]
         if len(_centroid_pdbs) != 1:
            raise ValueError(f"Expected 1 centroid pdb in {flankseqandbbclusdir}, "
                             f"but found {len(_centroid_pdbs)}.")

         pdb = _centroid_pdbs[0]
         # copy the centroid to the output dir
         biounit = '_'.join(pdb.split('_')[1:6]).replace('.pdb.gz', '')
         vdmscrr_list = pdb.split('_')[6:-3]
         # remove resnames for conciseness
         vdmscrr_list = [vdmscrr_list[i] for i in range(len(vdmscrr_list)) 
                         if (i + 1) % 4 != 0] 
         vdmscrr_str = '_'.join(vdmscrr_list)
         newname = f'{biounit}_{vdmscrr_str}.pdb.gz'
         newpath = os.path.join(nr_dir, newname)
         assert not os.path.exists(newpath)
         # attempt to copy it to the newpath, but sge sometimes gives a 
         # communication send error when there's lag, so allow 100 attempts
         try:
            shutil.copy(os.path.join(flankseqandbbclusdir, pdb), newpath) 
         except: # retry
            for attempt in range(100):
               try:
                  shutil.copy(os.path.join(flankseqandbbclusdir, pdb), newpath) 
                  break
               except Exception as e:  
                  if attempt < 99:
                     time.sleep(10)  # wait before retrying
                  else:
                     print(f"[WARNING] clus_and_deduplicate failed to create "
                           f"{os.path.join(flankseqandbbclusdir, pdb)}: {e}")

def cluster_vdgs_of_same_AA_comp(_vdgs, seq_sim_thresh, reordered_AAs, 
   symmetry_classes, vdglib_dir, align_cg_weight, num_flanking,
   atomgroup_dict, print_flankbb, logfile):
   # Cluster vdG subsets with identical vdm AA compositions to identify redundancies. 
   # First, get all AA permutations for identical AAs. For example, if the vdms are 
   # [Ala1, Ala2, Glu], then we need to sample [Ala1, Ala2, Glu] and [Ala2, Ala1, 
   # Glu]. The easiest way to get the permutations is to label the index of the 
   # ordered AAs
   reordered_AAs_str = '_'.join(reordered_AAs)
   
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
   cgvdmbb_clus_assignments, cgvdmbb_clus_centroids = clust.get_leader_clusters(
      cgvdmbb_data_to_clus, cgvdmbb_rmsd_cut, reordered_AAs_str, len(reordered_AAs), 
      vdglib_dir, final_exact_medoid_pass=True, final_reassign_once=True)

   # Output the cgvdmbb clusters
     
   # Define dir to output "temp" files (temp bc clusters need to be reassigned
   # later to handle degenerate binding sites for different AA and CG 
   # permutations of the same PDB).
   temp_cgvdmbb_clusdir = os.path.join(vdglib_dir, 'clusters', 'temp', 'cgvdmbb', 
      str(len(reordered_AAs)), reordered_AAs_str) 
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
      first_pdb_cg_vdmbb_coords, cgvdmbb_weights, atomgroup_dict, print_flankbb, 
      symmetry_classes, reordered_AAs)

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
                                    len(reordered_AAs), reordered_AAs_str)

   # Next, further subdivide each cluster on rmsd of the bb stretches flanking the vdms 
   # + their sequences. 
   for cgvdmbb_clusnum, indices_of_elements_in_cgvdmbb_cluster in \
      reassigned_cgvdmbb_clus.items():
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
      # Then, do secondary clustering on a metric and threshold that combines rmsd and seq
      flankseq_and_bb_thresh = flankbb_rmsd_cut + (1 - seq_sim_thresh) / 2
      flankingseq_and_bb_cluster_assignments, flankingseq_and_bb_clus_centroids \
         = clust.get_leader_clusters(flankseq_and_bb_to_clus, flankseq_and_bb_thresh, 
         reordered_AAs_str, len(reordered_AAs), vdglib_dir, 
         final_exact_medoid_pass=True, final_reassign_once=True)
            
      # Output the flankseq+bb clusters
      flankseq_and_bb_clusdir_for_this_cgvdmbb_clus = os.path.join(vdglib_dir, 'clusters', 
         'flankseq_and_bb', str(len(reordered_AAs)), reordered_AAs_str, 
         f'cgvdmbb_{cgvdmbb_clusnum}') 
      utils.handle_existing_files(flankseq_and_bb_clusdir_for_this_cgvdmbb_clus) 

      first_pdb_out, first_pdb_cg_vdmbb_coords, failed_pdbs = clust.write_out_clusters(
         flankseq_and_bb_clusdir_for_this_cgvdmbb_clus, 
         flankingseq_and_bb_cluster_assignments, flankingseq_and_bb_clus_centroids,
         cgvdmbb_clus_cg_coords, cgvdmbb_clus_pdbpaths, cgvdmbb_clus_vdm_scrr_cg_perm, 
         cgvdmbb_clus_cgvdmbb_coords, cgvdmbb_clus_flat_flankingCAs, num_flanking, 
         first_pdb_out, first_pdb_cg_vdmbb_coords, cgvdmbb_weights, atomgroup_dict, 
         print_flankbb, symmetry_classes, reordered_AAs, clusterlabel='flankseq_and_bb')

      if len(failed_pdbs) > 0:
         with open(logfile, 'a') as file:
            for fail in failed_pdbs:
               file.write(f'\tFailed to parse {fail}.\n')

def delete_clusterdirs(vdglib_dir, logfile, size_subset, keep_clustered_pdbs):
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
      # delete flankseq_and_bb unless keep_clustered_pdbs is True
      if not keep_clustered_pdbs:
         flank_dir = os.path.join(all_clusters, 'flankseq_and_bb', str(size_subset))
         if not os.path.exists(flank_dir):
            file.write(f"\t{flank_dir} does not exist.\n")
         else:
            shutil.rmtree(flank_dir)

def normalize_rmsd(num_atoms, atoms):
   ''' Return a size-normalized RMSD threshold (Å) for the given atom set ('cgvdmbb' or 
   'flankbb'). The threshold scales linearly with the number of atoms between 8 and 15:
      flankbb: 0.5 Å → 1.5 Å
      cgvdmbb: 0.5 Å → 1.0 Å
   Below 8 atoms, use the minimum; above 15, use the maximum.'''

   if atoms == 'flankbb':
      max_threshold = 1.5
      min_threshold = 0.5
   elif atoms == 'cgvdmbb':
      max_threshold = 1
      min_threshold = 0.5

   # Smoothly scale rmsd threshold (Å) by the number of atoms.
   min_atoms = 8
   max_atoms = 15
   
   if num_atoms < min_atoms:
      return min_threshold  # below 8 atoms, use the minimum threshold (0.5 Å)
   if num_atoms > max_atoms:
      return max_threshold  # above 15 atoms, use the maximum threshold
   
   # Linear scaling for atoms 
   scaling_factor = (num_atoms - min_atoms) / (max_atoms - min_atoms)
   threshold = min_threshold + scaling_factor * (max_threshold - min_threshold)
   
   return threshold

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

if __name__ == "__main__":
   main()
