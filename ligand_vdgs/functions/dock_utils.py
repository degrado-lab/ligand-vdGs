# dock_utils.py

from itertools import combinations, product
from functools import lru_cache
import numpy as np
import os
import re

def get_bsr_combinations(solved_struct, ligname, quiet=True):
    # Enumerate binding-site residue combinations for vdG matching.

    # 1) Find all binding-site residues once
    bindingsite_residues = get_bindingsite_residues(solved_struct, addl_residues=[],
        ligname=ligname, dist_from_lig=4.5, CA_only=False, quiet=quiet,)

    # 2) Precompute AA identity + backbone coords per residue
    #    key: (seg, chain, resnum)
    #    value: (AA_name, np.array[[N],[CA],[C]] with shape (3, 3))
    bb_cache = {}
    for seg, chain, resnum in bindingsite_residues:
        if seg == "":
            sele = f"chain {chain} and resnum {resnum}"
        else:
            sele = f"segment {seg} and chain {chain} and resnum {resnum}"

        res_obj = solved_struct.select(sele)
        if res_obj is None:
            continue

        bb_coords = []
        for atom_name in ["N", "CA", "C"]:
            atom = res_obj.select(f"name {atom_name}")
            if atom is None or atom.numAtoms() == 0:
                raise ValueError(f"Missing atom {atom_name} in residue "
                    f"{seg}:{chain}:{resnum}")
            bb_coords.append(atom.getCoords()[0])

        AA = get_res_AA_identity(res_obj)
        bb_cache[(seg, chain, resnum)] = (AA, np.asarray(bb_coords, dtype=np.float32),)

    # 3) Enumerate subsets
    bsr_combos = get_vdg_subsets(bindingsite_residues)

    # 4) Build all combinations of AA identities with 'bb' wildcards
    all_bsr_combos = []
    for bsr_combo in bsr_combos:
        bsr_AA_identities = []
        input_bsr_bb_coords = []

        for bsr in bsr_combo:
            if bsr not in bb_cache:
                # This should be rare; skip broken residues
                continue
            AA, bb_coords = bb_cache[bsr]
            bsr_AA_identities.append(AA)
            input_bsr_bb_coords.append(bb_coords)

        if not bsr_AA_identities:
            continue

        # Convert GLY -> bb
        bsr_AA_identities_conv_bb = [AA if AA != "GLY" else "bb" for AA in bsr_AA_identities]

        # Each AA can be itself or 'bb' (except literal 'bb', which stays 'bb')
        options = [(x,) if x == "bb" else (x, "bb") for x in bsr_AA_identities_conv_bb]
        combo_variants = list(product(*options))
        for c in combo_variants:
            bsr_aas_coords = (c, bsr_combo, bsr_AA_identities, input_bsr_bb_coords.copy())
            if bsr_aas_coords not in all_bsr_combos:
                all_bsr_combos.append(bsr_aas_coords)

    return all_bsr_combos

def get_bindingsite_residues(prody_obj, addl_residues, ligname, dist_from_lig=8, 
    CA_only=True, quiet=True):
    res = []
    # Use CA_only=True when you're doing blind docking and don't want to use sc info
    # Use CA_only=False when you want to use sc positions to define interactions
    if CA_only:
        _selection = "name CA"
    else:
        _selection = "protein"
    CAs = prody_obj.select(
        f"{_selection} within {dist_from_lig} of resname {ligname} and not element CA"
        )  # exclude calcium
    for ca in CAs:
        res_tup = (ca.getSegname(), ca.getChid(), ca.getResnum())
        if res_tup not in res:
            res.append(res_tup)
    if addl_residues:
        res += addl_residues
    for r in res:
        if not isinstance(r, tuple) or len(r) != 3:
            print(f"Invalid residue tuple: {r}")
    # Only print the pymol selection once (caller controls this via quiet)
    if not quiet:
        print("\nBinding site residues for pymol selection:\n")
        print("select bindingsite, " + " or ".join(
            f"(seg {seg} and chain {chain} and resi {resnum})" 
            for seg, chain, resnum in res) + "\n")
    return res

def get_vdg_subsets(input_list):
    # Initialize an empty list to store all subsets
    all_subsets = []
    # Loop through subset sizes 1& 2, generate combos, then add to list
    for r in range(1, 3):
        subsets = combinations(input_list, r)
        all_subsets.extend(subsets)
    return all_subsets

def select_residue(prody_obj, seg, chain, resnum):
    if seg == '':
        sele = f'chain {chain} and resnum {resnum}'
    else:
        sele = f'segment {seg} and chain {chain} and resnum {resnum}'
    res_obj = prody_obj.select(sele)
    return res_obj


def get_res_AA_identity(res_obj):
    assert len(set(res_obj.getResnames())) == 1
    AA = res_obj.getResnames()[0]
    return AA

def extract_elements(smiles: str):
    # Return a list of elements in the order they appear in the SMILES string.
    # Match:
    # - Two-letter uppercase elements (Cl, Br, Si, Na, Li)
    # - Single uppercase element (C, N, O, S, etc.)
    # - Single lowercase aromatic atom (c, n, o, s, p)
    pattern = r"(Cl|Br|Si|Na|Li|[A-Z]|[cnosp])"
    return re.findall(pattern, smiles)

def get_query_cg_coords(sub, sub_smiles):
    coords_list = []
    Mol_elements_list = []
    # Get atom names, etc.
    conf = sub.GetConformer()  # Get the 3D conformer to get coords
    for atom in sub.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())  # returns an RDKit Point3D object
        xyz = (pos.x, pos.y, pos.z)
        coords_list.append(xyz)
        Mol_elements_list.append(atom.GetSymbol())
    
    # Check that the atom orders are actually correct by checking the element names
    smiles_elements = extract_elements(sub_smiles)
    smiles_elements = [e.capitalize() for e in smiles_elements] # capital to match Mol elements
    if smiles_elements != Mol_elements_list:
        # Raise error
        raise ValueError(f"Element order mismatch between SMILES and RDKit Mol:\n"
                         f"SMILES elements: {smiles_elements}\n"
                         f"Mol elements: {Mol_elements_list}")
    return coords_list

def name_outdir(pdbfile, outdir, make_pdb_subfolder):
    # Name the outdir for this pdb query.
    # `make_pdb_subfolder` indicates whether to make a subfolder for each pdb within
    # the outdir (e.g. for multiple predictions of the same PDB).
    pdbname = os.path.basename(pdbfile)
    for ext in (".pdb.gz", ".pdb", ".cif.gz", ".cif"):
        if pdbname.endswith(ext):
            pdbname = pdbname[:-len(ext)]
            break
    pdb_id = pdbname[:4]
    if make_pdb_subfolder:
        output_dir = os.path.join(outdir, pdb_id, pdbname)
    else:
        output_dir = os.path.join(outdir, pdbname)
    return output_dir

@lru_cache(maxsize=8192)
def _map_aa_identities_to_vdg_resinds_cached(vdg_AAs_tuple, target_list_tuple, wildcard='bb'):
    vdg_AAs = list(vdg_AAs_tuple)
    target_list = list(target_list_tuple)
    AA_indices = {}

    # Handle regular AAs (non-wildcards, i.e. not `bb`)
    for AA in set(target_list):
        if AA != wildcard:
            AA_indices[AA] = [i for i, x in enumerate(vdg_AAs) if x == AA]
    
    # For wildcard, it can match any AA in the vdg (exclude lig)
    AA_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
                'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                'THR', 'TRP', 'TYR', 'VAL']
    if wildcard in target_list:
        AA_indices[wildcard] = [i for i, x in enumerate(vdg_AAs) if x in AA_codes]
    
    # Generate all permutations of indices
    result = []
    def backtrack(current_indices, target_position):
        # If we've matched all target elements, add the current permutation to result
        if target_position == len(target_list):
            result.append(current_indices.copy())
            return

        # Get the next element to match
        current_element = target_list[target_position]
        
        # Iterate over each possible index for this element
        for idx in AA_indices.get(current_element, []):
            # Make sure we don't use the same index twice
            if idx not in current_indices:
                current_indices.append(idx)
                backtrack(current_indices, target_position + 1)
                current_indices.pop()  # Backtrack
    backtrack([], 0)
    return result

def map_aa_identities_to_vdg_resinds(vdg_AAs, target_list, wildcard='bb'):
    '''
    Given a list of amino acids in the prody obj (`vdg_AAs`) and a list of AAs to 
    select for (`target_list`), find all permutations of indices (i.e., resindices) in the 
    vdg (`vdg_AAs`) that match the target list elements in order. The amino acid `bb` 
    is treated as a wildcard.
    '''
    return _map_aa_identities_to_vdg_resinds_cached(tuple(vdg_AAs), tuple(target_list), wildcard)
