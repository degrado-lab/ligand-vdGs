import os
import sys
import time
import gzip
import pickle
import argparse
import glob

import numpy as np
import prody as pr
import pandas as pd

from io import StringIO
from copy import deepcopy
from functools import partial
from multiprocessing import Pool
from subprocess import check_output
from scipy.spatial.transform import Rotation
from itertools import permutations as permute

from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import Draw, AllChem

sys.path.append(os.path.dirname(__file__))

from atomtypes import *
from functions import utils

RDLogger.DisableLog('rdApp.*')

aa_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
               'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
               'SER', 'THR', 'TRP', 'TYR', 'VAL']

#def sort_and_screen(prody_pkls, validation_df, max_res, max_robs):
#    """Sort a list of pickled ProDy objects by resolution and R_obs.
#
#    Parameters
#    ----------
#    prody_pkls : list
#        List of pickled ProDy objects.
#    validation_df : pandas.DataFrame
#        Dataframe containing validation information on all protein 
#        structures in the PDB.
#    max_res : float
#        Maximum allowable resolution in the sorted list.
#    max_robs : float
#        Maximum allowable R_obs in the sorted list.
#
#    Return
#    ------
#    sort_pkls : list
#        List of pickled ProDy objects that satisfy the resolution and 
#        R_obs filters and are sorted in descending order by quality.
#    """
#    sort_df = validation_df.sort_values('pdb_acc')
#    pdb_acc = np.sort([pkl.split('/')[-1][:4] for pkl in prody_pkls])
#    sub_df = sort_df[(sort_df.resolution < max_res) & 
#                     (sort_df.R_obs < max_robs)]
#    mask = np.in1d(pdb_acc, np.array(sub_df.pdb_acc))
#    unique_acc, inv_idxs = np.unique(pdb_acc[mask], return_inverse=True)
#    sub_df = sub_df[sub_df.pdb_acc.isin(unique_acc)]
#    assert len(sub_df) == len(unique_acc)
#    score_all = (np.array(sub_df.R_obs) - 
#                 1. / np.array(sub_df.resolution))[inv_idxs]
#    # score_all = (np.array(sub_df.resolution) / max_res + 
#    #              np.array(sub_df.R_obs) / max_robs)[inv_idxs]
#    sort_idxs = np.argsort(score_all)[::-1]
#    return list(prody_pkls[mask][sort_idxs])


def assign_bond_orders_from_template(template, mol):
    """Assign bond orders to mol from template if one is substruct of another.

    Parameters
    ----------
    template : rdkit.Chem.rdchem.Mol
        The molecule object from which to extract bond orders.
        Can be a substructure or superstructure of mol.
    mol : rdkit.Chem.rdchem.Mol
        The molecule object to which to assign bond orders, assumed to 
        have all single bond orders.

    Returns
    -------
    assigned_mol : rdkit.Chem.rdchem.Mol
        The molecule object with reassigned bond orders.
    """
    if mol.HasSubstructMatch(template):
        assigned_mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        return assigned_mol
    assigned_mol = deepcopy(mol)
    [a.SetNoImplicit(True) for a in assigned_mol.GetAtoms()]
    mol_noH = Chem.RemoveHs(assigned_mol, sanitize=False)
    mol_Hmatch = \
        invert_substruct_match(assigned_mol.GetSubstructMatch(mol_noH), 
                               assigned_mol.GetNumAtoms()) + (-1,)
    template_noH = Chem.RemoveHs(template, sanitize=False)
    template_singleb = deepcopy(template_noH)
    for b in mol_noH.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
    for b in template_singleb.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
    if mol_noH.HasSubstructMatch(template_singleb):
        match = mol_noH.GetSubstructMatch(template_singleb)
    elif template_singleb.HasSubstructMatch(mol_noH):
        match = invert_substruct_match(
            template_singleb.GetSubstructMatch(mol_noH), 
            template_singleb.GetNumAtoms())
    else:
        raise ValueError('One input is not a substruct of the other.')
    for b in template_noH.GetBonds():
        a1 = b.GetBeginAtomIdx()
        a2 = b.GetEndAtomIdx()
        a1m = mol_Hmatch[match[a1]]
        a2m = mol_Hmatch[match[a2]]
        try:
            mol_b = assigned_mol.GetBondBetweenAtoms(a1m, a2m)
            mol_b.SetBondType(b.GetBondType())
        except:
            pass
    return assigned_mol


def invert_substruct_match(match, n_atoms):
    """Get the inverse mapping of a substructure match.

    Parameters
    ----------
    match : tuple
        Tuple of atom indices in a mol corresponding to zero-indexed 
        atoms of a substructure mol.
    n_atoms : int
        Number of atoms in mol.

    Returns
    -------
    inv_match : tuple
        Tuple of atom indices in the substructure corresponding to 
        zero-indexed atoms in the original mol.  Atoms in the mol 
        without matches in the substructure are assigned -1.
    """
    inv_match = []
    for i in range(n_atoms):
        if i in match:
            inv_match.append(match.index(i))
        else:
            inv_match.append(-1)
    inv_match = tuple(inv_match)
    return inv_match


def get_segs_chains_resnums(atomgroup, selection):
    """Get a set containing tuples of segments, chains, and resnums.

    Parameters
    ----------
    atomgroup : prody.atomic.atomgroup.AtomGroup
        ProDy AtomGroup for a protein structure.
    selection : str
        String specifying subset of atoms in the ProDy selection algebra.

    Returns
    -------
    segs_chains_resnums : set
        Set containing tuples of segments, chains, and resnums for 
        each residue that matches the selection string.
    """
    sel = atomgroup.select(selection)
    if sel is None:
        return set()
    return set(zip(sel.getSegnames(), sel.getChids(), sel.getResnums()))


def get_bond_elements(mol):
    """Get pairs of atomic numbers for all bonds in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.

    Returns
    -------
    bond_elem : list
        List of tuples containing two sorted integers representing the 
        atomic numbers of the two atoms in every bond in the molecule, 
        excluding those involving hydrogen or atoms with atomic number 0.
    """
    bond_elem = []
    for b in mol.GetBonds():
        n1 = b.GetBeginAtom().GetAtomicNum()
        n2 = b.GetEndAtom().GetAtomicNum()
        if n1 > 1 and n2 > 1: # ensure no bonds involving atomic number 0 or 1
            if n2 > n1:
                bond_elem.append((n1, n2))
            else:
                bond_elem.append((n2, n1))
    return bond_elem


def extract_coords(mol, patt):
    """Extract coordinates of atoms in a molecule that match a SMARTS pattern.

    Parameters
    ----------
    mol : rdkit.Chem.rdmol.Mol
        Mol object in which to search for matches to the SMARTS pattern.
    patt : rdkit.Chem.rdmol.Mol
        Mol object generated from a SMARTS pattern by RDKit.

    Returns
    -------
    coords : np.array [n_matches x n_atoms x 3]
        Atomic coordinates of atoms that match the SMARTS pattern.
    """
    n_atoms = patt.GetNumAtoms()
    if mol.HasSubstructMatch(patt):
        mol_coords = np.array(mol.GetConformer().GetPositions())
        matches = mol.GetSubstructMatches(patt)
        n_matches = len(matches)
        coords = np.zeros((n_matches, n_atoms, 3))
        for i, match in enumerate(matches):
             coords[i] = mol_coords[np.array(match)]
        return coords
    else:
        return np.zeros((0, n_atoms, 3))


def permute_coords(coords, ref_coords, perms):
    """Find the permutation of coords that aligns best with ref_coords.

    Parameters
    ----------
    coords : np.array [M x N]
        A set of N-dimensional coordinates.
    ref_coords : np.array [M x N]
        A set of N-dimensional reference coordinates.
    perms : [n_permutations x M]
        A set of allowable permutations to test for the best alignment 
        between coords and ref_coords.

    Returns
    -------
    perm_coords : np.array [M x N]
        Coords after the optimal permutation has been applied.
    min_perm : np.array [M]
        Indices of the optimal permutation.
    """
    mean_pos = np.mean(coords, axis=0)
    mean_ref = np.mean(ref_coords, axis=0)
    min_rmsd = np.inf
    min_perm = perms[0]
    assert len(coords) >= len(min_perm)
    for perm in perms:
        rot, rmsd = Rotation.align_vectors(
            coords[perm] - mean_pos, ref_coords - mean_ref)
        if rmsd < min_rmsd:
            min_rmsd = rmsd
            min_perm = perm
    return coords[min_perm], min_perm


def get_permutations(patt, sym_classes):
    """Get symmetry-equivalent permutations of a SMARTS pattern.

    Parameters
    ----------
    patt : rdkit.Chem.rdmol.Mol
        RDKit Mol object representing a SMARTS pattern.
    sym_classes : list
        List of integers denoting symmetry classes of the atoms. If the 
        length of this list is not n_atoms, RDKit symmetry classes are used.

    Returns
    -------
    perms : np.array [n_permutations x n_atoms]
        Symmetry-equivalent permutations of the atoms in the SMARTS pattern.
    """
    n_atoms = patt.GetNumAtoms()
    patt.UpdatePropertyCache()
    if len(sym_classes) == n_atoms:
        sym_classes = np.array(sym_classes)
    else:
        sym_classes = np.array(Chem.CanonicalRankAtoms(patt, breakTies=False)) 
    if len(sym_classes) == len(np.unique(sym_classes)):
        perms = np.arange(n_atoms).reshape(-1, n_atoms)
    else:
        dt = np.dtype([('', np.int64)] * n_atoms)
        perms_iter = permute(np.arange(n_atoms))
        n_perms = np.math.factorial(n_atoms)
        relevant_perms = []
        group_size = 10000000
        for i in range(n_perms // group_size + 1): # memory-efficient
            if i == n_perms // group_size:
                count = n_perms % group_size
            else:
                count = group_size
            perms = np.fromiter(perms_iter, dt, count).view(
                np.int64).reshape(-1, n_atoms)
            sym_perms = sym_classes[perms]
            if not len(relevant_perms):
                sym_perm_0 = sym_perms[0]
            relevant_perms.append(
                perms[np.all(sym_perms == sym_perm_0, axis=1)])
        perms = np.vstack(relevant_perms)
    return perms


def sdf_to_matching_mols(master_sdf, patt, debug):
    """Get all molecules from an SDF that may match a SMARTS pattern.

    Parameters
    ----------
    master_sdf : str
        Path to SDF file to search for matching molecules.
    patt : rdkit.Chem.rdmol.Mol
        Mol object generated from a SMARTS pattern by RDKit.

    Returns
    -------
    mols : list
        List of RDKit Mol objects that may match the SMARTS pattern.
    molnames : list
        List of three-letter accession codes for the matching ligands.
    """
    patt_bond_elem = get_bond_elements(patt)
    mols = []
    molnames = []
    for mol in Chem.SDMolSupplier(master_sdf, removeHs=False):
        if mol:
            molname = mol.GetProp("_Name")
            #if debug: 
            #    if molname != '0IV':
            #        continue
            # reject some molecules guaranteed not to match (to save time)
            mol_bond_elem = get_bond_elements(mol)
            if not all(pair in mol_bond_elem for pair in patt_bond_elem):
                continue
            mols.append(mol)
            molnames.append(molname)

    return mols, molnames



def has_substruct(mol, patt, mol_name, pdb_acc):
    ''' Substruct matching sometimes works for MolFromSmiles but not MolFromSmarts, 
        or vice versa, so try both.'''
    patt_from_smarts = patt
    patt_converted_to_smiles = Chem.MolFromSmiles(Chem.MolToSmiles(patt))
    patt_to_smiles_then_smarts = Chem.MolFromSmarts(Chem.MolToSmiles(patt))
    if mol.HasSubstructMatch(patt_from_smarts):
        return (True, patt_from_smarts)
    if mol.HasSubstructMatch(patt_converted_to_smiles):
        return (True, patt_converted_to_smiles)
    if mol.HasSubstructMatch(patt_to_smiles_then_smarts):
        print(f'Very odd case where substruct was only found when smiles converted to Mol converted to smiles again converted to Mol and then smarts, may need to investigate: {mol_name} in {pdb_acc}')
        return (True, patt_to_smiles_then_smarts)
    return (False, None)

def update_list_of_matches(mol, match_molname, segname, chain, resnum, match_mols, match_molnames,
                           segs_chains_resnums, new_patts, new_patt):
    match_mols.append(mol)
    match_molnames.append(match_molname)
    segs_chains_resnums.append((segname, chain, resnum))
    new_patts.append(new_patt)
    return match_mols, match_molnames, segs_chains_resnums, new_patts



def find_ref_coords(prody_pkl, patt, mols, molnames, validation_df, pdb_dataset_path, 
                    pdb_ligs_dict, outputfile, debug):
    """Find reference coordinates for a SMARTS pattern in a pickled ProDy obj.

    Parameters
    ----------
    prody_pkl : str
        Path to pickled ProDy object.
    patt : rdkit.Chem.rdmol.Mol
        Mol object generated from a SMARTS pattern by RDKit.
    mols : list
        List of RDKit Mol objects that may match the SMARTS pattern.
    molnames : list
        List of three-letter accession codes for the matching ligands.
    validation_df : pandas.DataFrame
        Dataframe containing validation information on all protein 
        structures in the PDB.

    Returns
    -------
    ref_coords : np.array [n_atoms x 3]
        Numpy array of reference atomic coordinates for the chemical group.
    """
    match_mols, _, _, new_patts = find_matches(prody_pkl, patt, mols, molnames, 
                                    validation_df, pdb_dataset_path, 
                                    pdb_ligs_dict, outputfile, debug, halt_early=True)
    if len(match_mols):
        return extract_coords(match_mols[0], new_patts[0])[0]
    else:
        return


def find_chain_clusters(prot_chains, pdb_acc):
    """Get mmseqs cluster numbers for a set of chains in a structure.
    
    Parameters
    ----------
    prot_chains : iterable
        One-letter chain IDs of all chains for which to find cluster numbers.
    pdb_acc : str
        Four-letter PDB accession code for the structure from which the 
        chains were pulled.
    pdb_cluster_file : str
        Path to file containing clusters for all protein chains in the PDB.

    Returns
    -------
    clustered_prot_chains : list
        List of one-letter IDs of chains in the structure that have been 
        successfully clustered.
    prot_chain_seq_clusters : list
        List of cluster IDs for each successfully clustered chain.
    """
    clustered_prot_chains = []
    prot_chain_seq_clusters = []
    fake_clusnum = 0
    # hacky: return everything, so that we do not filter on homology
    for prot_chain in prot_chains:
        clustered_prot_chains.append(prot_chain)
        prot_chain_seq_clusters.append(int(fake_clusnum))
        fake_clusnum += 1
        ## Ignore below to remove homology filtering
        #cmd = ['awk', 
        #       '/{}_{}/{{print NR}}'.format(pdb_acc.upper(), prot_chain), 
        #       pdb_cluster_file]
        #output = check_output(cmd, encoding=sys.stdout.encoding)
        #output = output.split('\n')[0]
        #if output:
        #    clustered_prot_chains.append(prot_chain)
        #    prot_chain_seq_clusters.append(int(output))
    return clustered_prot_chains, prot_chain_seq_clusters


def get_bonded_dict(mol):
    """For a molecule, get a dict pairing atoms to lists of bonded atom names.

    Parameters
    ----------
    mol : rdkit.Chem.rdmol.Mol
        RDKit Mol object for which to find the names of bonded atoms to 
        each atom.

    Returns
    -------
    bonded_dict : dict
        Dict that pairs names of atoms in the Mol object with lists of 
        names of atoms bonded to that atom.
    """
    names = [a.GetPDBResidueInfo().GetName().strip() for a in mol.GetAtoms()]
    bonded_dict = {name : [] for name in names}
    for bond in mol.GetBonds():
        name0 = names[bond.GetBeginAtomIdx()]
        name1 = names[bond.GetEndAtomIdx()]
        bonded_dict[name0].append(name1)
        bonded_dict[name1].append(name0)
    return bonded_dict


def pkl_to_df(prody_pkl, pdb_dataset_path, pdb_ligs_dict, patt, mols, molnames,
              probe_pkl_dir, ref_coords, perms, precedence, validation_df, 
              max_res, max_robs, require_lig_contacts, lig_contact_atypes,
              outfile_prefix, logdir, debug):
              
    """Given a pickled ProDy object and a SMARTS pattern in RDKit Mol form, 
       find all protein contacts with atoms that match the pattern and 
       return a dataframe with information on these contacts.

    Parameters
    ----------
    prody_pkl : str
        Path to pickled ProDy object.
    patt : rdkit.Chem.rdmol.Mol
        Mol object generated from a SMARTS pattern by RDKit.
    mols : list
        List of RDKit Mol objects that may match the SMARTS pattern.
    molnames : list
        List of three-letter accession codes for the matching ligands.
    pdb_cluster_file : str
        Path to file containing clusters of chains in the PDB.
    probe_pkl_dir : list
        Path to directory in which to find pickled Probe files for 
        ligands within the corresponding pickled ProDy objects.
    ref_coords : np.array [n_atoms x 3] 
        Numpy array of reference atomic coordinates for the chemical group.
    perms : np.array [n_permutations x n_atoms]
        Numpy array of symmetry-equivalent permutations of atoms within the 
        chemical group to test against the reference coordinates.
    precedence : list
        List of integers used to determine the precedence of atoms in the 
        chemical group, with higher scores denoting higher precedence.
    validation_df : pandas.DataFrame
        Dataframe containing validation information on all protein 
        structures in the PDB.
    max_res : float
        Maximum allowable resolution in the sorted list.
    max_robs : float
        Maximum allowable R_obs in the sorted list.
    require_lig_contacts : bool
        If True, require that all ligands included in the final pkl have at 
        least one contact with another ligand.
    lig_contact_atypes : list
        A list of atom types allowable for contacts between the chemical 
        group and another ligand, if require_lig_contacts is True. If this 
        argument is None, all atom types are allowed.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing information on contacts between chemical 
        groups matching the SMARTS pattern and protein chains in the 
        pickled ProDy object.
    """
    df_dict = {'probe_name' : [], 'lig_resnum' : [], 'lig_resname' : [], 
               'prot_chain' : [], 'prot_chain_cluster' : [], 
               'prot_resname' : [], 'score' : [], 'prec_score' : [], 
               'match_number' : [], 'generic_name' : [], 'name' : [], 
               'atom_type' : [], 'rosetta_atom_type' : [], 'bonded_to' : []}
 
    bio_name = prody_pkl.split('/')[-1][:-4]
    pdb_acc = bio_name[:4].lower()
    outputfile = f'{logdir}/{outfile_prefix}_{os.getpid()}.txt'  # for this specific thread


    ##validation_data = validation_df[
    ##    (validation_df.pdb_acc == pdb_acc) & 
    ##    (validation_df.resolution < max_res) & 
    ##    (validation_df.R_obs < max_robs)]
    #validation_data = validation_df[
    #    (validation_df.pdb_acc == pdb_acc)]
    #if len(validation_data):
    #    # Select only the first biounit as bug fix
    #    validation_data = validation_data.iloc[0]
    #    score = float(validation_data['R_obs']) - \
    #            1. / float(validation_data['resolution'])
    #    # score = float(validation_data['resolution']) / max_res + \
    #    #         float(validation_data['R_obs']) / max_robs
    #else:
    #    return
    score = 0 # hacky: set score to some number because we are no longer providing validation_df
    n_atoms = patt.GetNumAtoms()
    generic_names = ['atom' + str(i) for i in range(n_atoms)]
    # patt may be re-defined if it needed to be converted to a different form
    match_mols, match_molnames, segs_chains_resnums, new_patts = \
        find_matches(prody_pkl, patt, mols, molnames, validation_df, pdb_dataset_path, 
                     pdb_ligs_dict, outputfile, debug)


    if not len(match_mols):
        pdb_error = f'Failed to find SMARTS pattern in pdb {pdb_acc}'
        with open(outputfile, 'a') as outF: # write in output file
            outF.write('======================================= \n')
            outF.write(pdb_error + '\n')
        return
    # iterate over matching molecules to find contacts
    prev_name = ''
    for mol, resname, scr, patt in zip(match_mols, match_molnames, 
                                 segs_chains_resnums, new_patts):
        # get atom names and atom types of matching chemical groups
        if patt is None: 
            raise ValueError(f'patt is None for {resname} {scr} {pdb_acc}')
        all_coords = extract_coords(mol, patt)
        all_matches = mol.GetSubstructMatches(patt)
        all_match_names = \
            [[mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetName().strip() 
              for i in match] for match in all_matches]
        atypes, rtypes = assign_combs_types(mol, return_rosetta_types=True)
        all_match_atypes = \
            [[atypes[i] for i in match] for match in all_matches]
        all_match_rtypes = \
            [[rtypes[i] for i in match] for match in all_matches]
        bonded_dict = get_bonded_dict(mol)
        # unpickle probe file to determine ligand-protein contacts
        segment, chain, resnum = scr
        probe_name = bio_name + '_' + segment + '_' + chain
        probe_filename = os.path.join(probe_pkl_dir, probe_name+'.pkl')
        #probe_filename = probe_pkl_dir + pdb_acc[1:3] + '/' + \
        #                 probe_name + '.pkl'
        if not os.path.exists(probe_filename):
            # Sometimes, a probe file is missing for a certain chain, but there exists a probe
            # file for another chain within the same pdb; this is because it was determined to
            # be a redundant binding site in run_ligand_parsing.py. So only report the cases
            # where it's missing a probe file *and* there is no probe file for this pdb at all.
            if any(file.startswith(bio_name) for file in os.listdir(probe_pkl_dir)):
                continue
            else:  # report as error
                print(probe_filename, prody_pkl)
                raise ValueError()
            #    probe_error = f"File {probe_filename} does not exist!"
            #    with open(outputfile, 'a') as outF: # write in output file
            #        outF.write('======================================= \n')
            #        outF.write(probe_error + '\n')
        if probe_name != prev_name: # prevent unnecessary unpickling
            with open(probe_filename, 'rb') as inpkl:
                probe_df = pickle.load(inpkl)
        # extract contacts from probe dataframes
        prev_name = probe_name
        counter = 0
        for match_names, match_atypes, match_rtypes, coords \
                in zip(all_match_names, all_match_atypes, 
                       all_match_rtypes, all_coords):
            contact_df = probe_df[(probe_df.chain1 == chain) & 
                                  (probe_df.resnum1 == resnum) &
                                  (probe_df.name1.isin(match_names))]
            if require_lig_contacts:
                contact_lig_df = \
                    contact_df[~contact_df.resname2.isin(aa_resnames)]
                if not len(contact_lig_df):
                    continue # ensure CG contacts non-AA residues
                if lig_contact_atypes is not None:
                    disallowed = len(contact_lig_df[
                        ~contact_lig_df.atomtype2.isin(lig_contact_atypes)])
                    if disallowed:
                        continue # ensure only allowed atom types contact CG
            prec_score = 0
            if len(precedence) == len(match_names):
                for k, mname in enumerate(match_names):
                    prec_score -= list(contact_df.name1).count(mname) * \
                                  precedence[k]
            prot_chains = contact_df['orig_chain2'].unique()
            # TODO: clean up, this is only to get rid of deduplicate flag (clustering)
            clustered_prot_chains, prot_chain_seq_clusters = \
                find_chain_clusters(prot_chains, pdb_acc)
            #clustered_prot_chains = [] 
            #prot_chain_seq_clusters = []
            #if not len(clustered_prot_chains):
            #    continue
            coords, perm = permute_coords(coords, ref_coords, perms)
            perm_names = [match_names[i] for i in perm] 
            perm_rtypes = [match_rtypes[i] for i in perm]
            perm_atypes = [match_atypes[i] for i in perm]
            bonded = [bonded_dict[name] for name in perm_names]
            for prot_chain, cluster in zip(clustered_prot_chains, 
                                           prot_chain_seq_clusters):
                prot_resnames = contact_df[contact_df.orig_chain2 
                                           == prot_chain]['resname2'].unique()
                for prot_resname in prot_resnames:
                    df_dict['probe_name'].extend([probe_name] * n_atoms)
                    df_dict['lig_resnum'].extend([resnum] * n_atoms)
                    df_dict['lig_resname'].extend([resname] * n_atoms)
                    df_dict['prot_chain'].extend([prot_chain] * n_atoms)
                    df_dict['prot_chain_cluster'].extend([cluster] * n_atoms)
                    df_dict['prot_resname'].extend([prot_resname] * n_atoms)
                    df_dict['score'].extend([score] * n_atoms)
                    df_dict['prec_score'].extend([prec_score] * n_atoms)
                    df_dict['generic_name'].extend(generic_names)
                    df_dict['name'].extend(perm_names)
                    df_dict['atom_type'].extend(perm_atypes)
                    df_dict['rosetta_atom_type'].extend(perm_rtypes)
                    df_dict['bonded_to'].extend(bonded)
                    df_dict['match_number'].extend([counter] * n_atoms)
                    counter += 1
    return pd.DataFrame.from_dict(df_dict)


def parse_args():
    argp = argparse.ArgumentParser()
    argp.add_argument('--smarts', help="SMARTS pattern to search for.")
    argp.add_argument('--ligand-db', help="Path to ligand database generated in "
                      "the previous step.")
    argp.add_argument('-o', '--outdir', default="databases/vdms",
                      help="Directory at which to "
                      "output the final pkl and CSV files.")
    argp.add_argument('--master-sdf', help="Path to a file that "
                      "aggregates all the ligand sdf files in the dataset "
                      "you are mining. If your dataset is the entire PDB, "
                      "you can find this file at: "
                      "http://ligand-expo.rcsb.org/dictionaries/Components-pub.sdf", 
                      default='parse_PDB_FGs/Components-pub.sdf')
    argp.add_argument('-t', '--threads', type=int, default=8, 
                      help="Number of threads with which to run the "
                      "search in parallel.")
    argp.add_argument('--sym-classes', default=[], type=int, nargs='+', 
                      help="List of integers denoting symmetry classes of "
                      "atoms in the chemical group. By default, RDKit "
                      "determines these.")
    argp.add_argument('--precedence', default=[], type=int, nargs='+', 
                      help="List of integers used to determine the "
                      "precedence of atoms in the chemical group, with " 
                      "higher scores denoting higher precedence.")
    argp.add_argument('--require-lig-contacts', action='store_true', 
                      help="If True, require that all ligands included "
                      "in the final pkl have at least one contact with "
                      "another ligand.")
    argp.add_argument('--lig-contact-atypes', nargs='+', 
                      help="A space-separated list of allowed atom types "
                      "for the contacting ligands,  if --require-lig-contacts "
                      "is also set.")
    argp.add_argument('--debug', action='store_true',
                    help="This debug flag deactivates multiprocessing; don't "
                    "rely on t=1 to turn of MP, because even with t=1, "
                    "Exceptions will be processed in a way that does not "
                    "print stack trace, making debugging impossible. This "
                    "flag allows the script to run until an Exception is raised, "
                    "rather than passing over the Exception.")
    argp.add_argument('--num-debug-pdbs', default=50,
                      help="Number of pdbs to process when using the debug flag.")
    argp.add_argument('--csv', action='store_true', 
                      help="If True, output a CSV file in addition to "
                      "the pickled dataframe.")
    #argp.add_argument('--validation-pkl', help="Path to the validation "
    #                  "pickle file; this is generated by parse_PDB_ligands/")
    argp.add_argument("--pdb-dataset-path", default='databases/consolidated_BioLiP2',
                    help="Path to your dir. of PDBs to mine FGs; if mining the RCSB PDB, this "
                    "should be the path to the consolidated BioLiP2 dataset. Refer to "
                    "docs/database_generation.txt")
    argp.add_argument('--pdb-ligs-dict', default='resources/database_ligs_dict.pkl', 
                    help="Path to pickled dict mapping lig names to the PDBs containing them. "
                    "Defaults to a package-provided file; "
                    "refer to scripts/mapping_lig_to_pdb.py .")
    argp.add_argument('--overwrite', action='store_true',
                      help="Overwrite existing output files; otherwise, program will terminate and "
                      "require that you manually remove them, so as to prevent accidental data loss.")
    
    
    return argp.parse_args()


def main():

    args = parse_args()
    master_sdf = args.master_sdf 
    #pdb_cluster_file = args.homology_file
    # generated by ligand_database.py
    ligand_db_dir = args.ligand_db
    pdb_ligs_dict_path = args.pdb_ligs_dict
    pdb_dataset_path = args.pdb_dataset_path
    prody_pkl_dir =  os.path.join(ligand_db_dir, 'prody') 
    probe_pkl_dir =  os.path.join(ligand_db_dir, 'probe') 
    validation_pkl = os.path.join(ligand_db_dir, 'dataframe_validation.pkl')
    overwrite = args.overwrite
    debug = args.debug


    #with open(validation_pkl, 'rb') as inpkl:
    #    validation_df = pickle.load(inpkl)
    #validation_df = validation_df[
    #    (validation_df.resolution != 'NotAvailable') & 
    #    (validation_df.R_obs != 'NotAvailable')].astype({'resolution' : float, 
    #                                                     'R_obs' : float})
    validation_df = None

    prody_pkls = []
    prody_pklfiles = sorted(os.listdir(prody_pkl_dir))
    
    if debug:
        print('\n\t\t ===== Debug mode ===== \n')

        # contains benzene-containing ligands found in databases/smarts_matches/c1ccccc1/c1ccccc1_ligands.sdf,
        # like in 0IV in 1bb0
        if args.smarts == 'c1ccccc1':
            prody_pklfiles = prody_pklfiles[450:500] 


    patt = Chem.MolFromSmarts(args.smarts)
    n_atoms = patt.GetNumAtoms()
    perms = get_permutations(patt, args.sym_classes)
    precedence = args.precedence
    if len(precedence) and len(precedence) != n_atoms:
        print('Warning: incorrect length of precedence vector.')
    mols, molnames = sdf_to_matching_mols(master_sdf, patt, debug)
    # Filter the pdbs to search through, based on whether they contain the ligand(s) of interest
    pdb_ligs_dict = pickle.load(open(pdb_ligs_dict_path, 'rb'))
    for mol_name in molnames:
        pdbs_containing_molname = pdb_ligs_dict[mol_name]
        for pdb_containing_molname in pdbs_containing_molname:
            for f in prody_pklfiles:
                if pdb_containing_molname in f:
                    prody_pkls.append(os.path.join(prody_pkl_dir, f))
    assert len(prody_pkls) > 0

    # Define log file to write out failed pdbs and the errors/exceptions encountered
    filename_time = time.time()
    output_prefix = f'smarts_to_df_{filename_time}'
    logdir = os.path.join(args.outdir, 'logs')
    if not os.path.exists(logdir):
        os.makedirs(logdir) 
    
    merged_filename = os.path.join( # called "merged" bc will merge child pids
        logdir, f'pooled_log_{args.smarts}_{filename_time}.txt')
    


    ref_coords = None
    i = -1
    while ref_coords is None and i + 1 < len(prody_pkls):
        i += 1
        ref_coords = find_ref_coords(prody_pkls[i], patt, mols, molnames, 
                                     validation_df, pdb_dataset_path, pdb_ligs_dict, merged_filename, debug)
    if ref_coords is None:
        raise ValueError('SMARTS pattern cannot be found in the entire PDB.')



    func = partial(pkl_to_df, pdb_dataset_path=pdb_dataset_path, 
                   pdb_ligs_dict=pdb_ligs_dict, patt=patt, mols=mols, molnames=molnames, 
                   probe_pkl_dir=probe_pkl_dir, ref_coords=ref_coords, 
                   perms=perms, precedence=precedence,  
                   validation_df=validation_df, max_res=2.5, max_robs=0.3, 
                   require_lig_contacts=args.require_lig_contacts, 
                   lig_contact_atypes=args.lig_contact_atypes,
                   outfile_prefix=output_prefix, logdir=logdir, debug=debug)
    
    start_time = time.time()
    if debug: # Do not multithread
        dfs = [func(x) for x in prody_pkls[i:]]
    else: # Multiprocessing
        with Pool(args.threads) as pool:
            dfs = pool.map(func, prody_pkls[i:])

    dfs = [df for df in dfs if df is not None]

    if len(dfs):
        full_df = pd.concat(dfs).astype({'lig_resnum' : int, 
                                         'prot_chain_cluster' : int, 
                                         'match_number' : int, 
                                         'prec_score' : int})
        full_df.drop(columns=['prot_chain_cluster', 'match_number', 
                              'prec_score', 'score'], 
                     inplace=True)
        
        # Write out results (chains, resnums, etc. of FG (CG) instances)
        fg_dir = os.path.join(args.outdir, 'FG_instances')
        cg_pkl_file = os.path.join(fg_dir, f'{args.smarts}.pkl')
        cg_csv_file = os.path.join(fg_dir, f'{args.smarts}.csv')
        

        
        
        with open(cg_pkl_file, 'wb') as pkl: 
            pickle.dump(full_df, pkl, protocol=5)
        if args.csv:
            full_df.to_csv(cg_csv_file)
    else:
        raise ValueError('No matches found; please check inputs.')

    # Report stats
    num_prody_pklfiles = len(prody_pklfiles) # a separate pkl file for each biounit
    uniquepdbs = set([i[:4] for i in prody_pklfiles])
    num_matches = len(full_df) // n_atoms
    # Number of ligands that contain those matches might be a different number because if, for 
    # example, you have dichlorobenzene, it would return 2 matches.
    num_lig_matches = full_df.groupby(['probe_name', 'lig_resnum']).ngroups # ligand instances, not unique ligands
    num_probe_files = len(full_df.probe_name.unique())
    num_pdbs_with_ligs = len(prody_pkls)
    num_unique_ligs = len(set(full_df['lig_resname']))
    counts = f'Out of {num_prody_pklfiles} biounits in {len(uniquepdbs)} PDBs, there were {num_pdbs_with_ligs} PDBs that have ligands containing this smarts. ({num_probe_files} probe files, i.e. chains). {num_lig_matches} ligand instances were found in the whole database ({num_unique_ligs} unique ligands), and yielded {num_matches} successfully processed matches (in some cases, probe file does not exist, occupancy is missing, rdkit could not process the molecule bc internal (within lig) clashing leads rdkit to believe that there are more internal bonds than there are, etc.). Please note that this is still under development, as num_matches may be overcounted; please see TODO (debugging example)'
    print(counts)
    assert len(prody_pkls) > 0
    
    # Merge and clean up (remove) all output files, because there's a separate
    # output file for each thread. 

    with open(merged_filename, 'w') as merged_file:
        for filname in glob.glob(
            f'"{logdir}/{output_prefix}"*'): # asterisk has to be outside of quotations because of
                                             # brackets in smarts
            with open(filname, 'r') as inF:
                lines = inF.readlines()
                for line in lines:
                    merged_file.write(line)
            os.remove(filname) # clean up
        merged_file.write(counts) # report counts


    final_time = time.time() - start_time
    final_time_hrs = final_time/60/60
    print(f'Final time: {round(final_time_hrs, 1)} hrs using {args.threads} threads \n') 
    print('Please refer to this outfile to see how many PDBs passed/failed:')
    print('\t', merged_filename)

if __name__ == "__main__":
    # Print out entire command entered on command line
    print("\nFull command:", " ".join(sys.argv), '\n')
    main()
