'''
Removes redundant vdGs from the vdG library. Redundant vdGs are defined as having the same: 
    1. CG and vdM AA identities/positions (backbones within a specified RMSD threshold)
    2. sequence similarity of residues flanking the vdMs (specified by seq. similarity 
       threshold)
    3. similar positions of the flanking residues (backbones within a specified RMSD 
       threshold)
'''

import os
import sys
import numpy as np
import prody as pr

# TODO: convert the following arguments to argparse-compatible 

script, CG = sys.argv

vdg_pdbs_dir = f'/wynton/group/degradolab/skt/docking/databases/vdg_lib/{CG}/vdg_pdbs' # TODO: remove hardcoding

rmsd_thresh = 0.75
seq_sim_thresh = 0.75
num_flanking = 5

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
                                             [  [CG coords], 
                                                [bb N-Ca-C coords of Ala1, Ala2, Cys3, Phe4], 
                                                [seq. +/- 5 of Ala1, Ala2, Cys3, Phe4], 
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
AA identity (such as having 2 Ala's in this example), all permutations must be sampled 
when checking for redundancy.
'''

def main():
   # Initialize vdm_combos dict to store the subsets of vdMs within a vdG
   vdm_combos = {}

   # Iterate over the PDBs and CGs that were identified as containing the 
   for pdbname in os.listdir(vdg_pdbs_dir):
      pdbpath = os.path.join(vdg_pdbs_dir, pdbname)
      prody_obj = pr.parsePDB(pdbpath)
      cg_coords = get_cg_coords(prody_obj)
      # Characterize the vdM residues (bb coords, flanking residues, pdb paths, etc.)
      vdms_dict = get_vdm_res_features(prody_obj, pdbpath)
      # Determine the vdM combinations, up to 4 residues

def get_vdm_res_features(prody_obj, pdbpath):
   # Identify the vdM residues (occ == 2). To be safe, select > 1.5 and < 2.5.
   vdm_residues = prody_obj.select('(occupancy) > 1.5 and (occupancy < 2.5)')
   vdm_resinds = set(vdm_residues.getResindices())
   # Record features of the vdm residues (bb coords, flanking residues, pdb paths, etc.)
   vdms_dict = {}
   print('-----------')
   print(pdbpath)
   for vdm_resind in vdm_resinds:
      vdm_obj = vdm_residues.select(f'resindex {vdm_resind}')
      print('--------- CENTRAL RES', vdm_obj.getResnums()[0])
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
      # LEFT OFF HERE. I guess now, need to go compare subsets of vdgs for redundancy

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

#def get_vdm_combos():
#   # Group by quadruples, triples, pairs, and singles.


if __name__ == "__main__":
    main()
