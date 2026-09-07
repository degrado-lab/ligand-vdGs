# dock_utils.py

from itertools import combinations, product
from functools import lru_cache
import numpy as np

from ligand_vdgs.functions import vdg_struct_utils as struct_utils
from ligand_vdgs.functions.vdg_struct_utils import get_res_AA_identity
from ligand_vdgs.functions.utils import mol_from_fragment

# Largest vdG subset size the library is built for; caps BSR combo enumeration.
#
#   - AA-slot symmetry is minimized over *inside* the clustering distance
#     (vdg_fp_utils.build_perm_group), so a repeated-label bucket costs prod(k!)
#     permutations per pair-distance -- 2x at sizes 1-2, but 6x for a bb_bb_bb
#     bucket at size 3, on top of an already O(N^2) clustering.
#   - The CG_VDM_CONTACT_CUTOFF guard in clus_and_deduplicate_vdgs checks each
#     vdM slot but discards the *whole* environment when any one fails, so a
#     non-contacting slot takes its genuinely-contacting partners with it. At
#     size 2 that costs 0.133% of vdGs (measured over 1,501 size-2 nr vdGs),
#     which is why the simple version was kept. The loss grows with subset size:
#     more slots means more chances that at least one fails, and more good slots
#     discarded each time it does. Above 2, drop only the offending slot -- which
#     means re-deriving the subset's aa_bucket and de-duplicating the remainder
#     against the smaller-subset vdGs that generation already emits.
# Raising this means fixing both, not editing this line.
MAX_SUBSET_SIZE = 2

def get_bsr_combinations(solved_struct, ligname, quiet=True, pdbfile=""):
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

        # (seg, chain, resnum) is not unique when insertion codes are in use: the
        # selection pulls every icode variant at once, and get_bb_coords would take
        # N from one residue and CA from another. Skip rather than merge.
        if len(set(res_obj.getIcodes())) > 1:
            print(
                f"[WARNING] ({pdbfile}) Residue {seg}:{chain}:{resnum} has multiple "
                f"insertion codes; skipping (insertion codes are not part of the BSR key).",
                flush=True,)
            continue

        # Same backbone policy as the library side (best altloc, finite and
        # non-collinear N/CA/C); returns None if any of that fails.
        bb_coords = struct_utils.get_bb_coords(res_obj)
        if bb_coords is None:
            print(
                f"[WARNING] ({pdbfile}) Residue {seg}:{chain}:{resnum} has no usable "
                f"N/CA/C backbone; skipping residue.", flush=True,)
            continue

        AA = get_res_AA_identity(res_obj)
        if AA is None:
            continue
        bb_arr = np.asarray(bb_coords, dtype=np.float32)
        bb_arr.setflags(write=False)  # shared by every combo variant below
        bb_cache[(seg, chain, resnum)] = (AA, bb_arr)

    # 3) Enumerate subsets
    bsr_combos = get_vdg_subsets(bindingsite_residues)

    # 4) Build all combinations of sidechain / backbone labels
    seen_bsr_keys = set()  # (c, bsr_combo): both are tuples of hashable primitives
    all_bsr_combos = []
    for bsr_combo in bsr_combos:
        bsr_AA_identities = []
        input_bsr_bb_coords = []

        for bsr in bsr_combo:
            if bsr not in bb_cache:
                print(
                    f"[WARNING] BSR residue {bsr} not found in backbone cache; skipping.",
                    flush=True,)
                continue
            AA, bb_coords = bb_cache[bsr]
            bsr_AA_identities.append(AA)
            input_bsr_bb_coords.append(bb_coords)

        if len(bsr_AA_identities) != len(bsr_combo):
            continue  # one or more residues missing backbone atoms; skip whole combo

        # Each residue is tried both as its own sidechain identity and as the
        # single backbone label; whether a backbone geometry is physically
        # hostable is decided on the read path, not here (see
        # hit_finder_core.backbone_slots_can_host).
        options = [struct_utils.query_slot_labels(AA) for AA in bsr_AA_identities]

        combo_variants = list(product(*options))
        for c in combo_variants:
            key = (c, bsr_combo)
            if key not in seen_bsr_keys:
                seen_bsr_keys.add(key)
                # Coord arrays are the read-only ones from bb_cache, shared
                # across variants; consumers copy on conversion.
                bsr_aas_coords = (c, bsr_combo, bsr_AA_identities, tuple(input_bsr_bb_coords))
                all_bsr_combos.append(bsr_aas_coords)

    return all_bsr_combos

def get_bindingsite_residues(prody_obj, addl_residues, ligname, dist_from_lig=None,
    CA_only=True, quiet=True):
    res = []
    # Use CA_only=True when you're doing blind docking and want only C-alpha positions.
    # Use CA_only=False to select all protein atoms (not just CA) near the ligand.
    if dist_from_lig is None:
        dist_from_lig = 8 if CA_only else 4.5
    if CA_only:
        _selection = "name CA"
    else:
        _selection = "protein"
    atoms = prody_obj.select(
        f"({_selection} and not element Ca CA and not resname CA) "
        f"within {dist_from_lig} of resname {ligname}")  
          # exclude calcium (name CA can match Ca2+ ions, not just alpha carbons);
          # element casing ("Ca" vs "CA") isn't consistent across PDB writers
    seen = set()  # membership set alongside the list, which keeps atom order
    for atom in atoms:
        res_tup = (atom.getSegname(), atom.getChid(), atom.getResnum())
        if res_tup not in seen:
            seen.add(res_tup)
            res.append(res_tup)
    if addl_residues:
        res += addl_residues
    for r in res:
        if not isinstance(r, tuple) or len(r) != 3:
            print(f"[WARNING] Invalid residue tuple: {r}")
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
    for r in range(1, MAX_SUBSET_SIZE + 1):
        subsets = combinations(input_list, r)
        all_subsets.extend(subsets)
    return all_subsets

@lru_cache(maxsize=1024)
def cg_element_symbols(cg_smarts):
    """Per-atom element symbols for a CG pattern, in the pattern's own atom order.

    Read off the parsed mol rather than tokenized out of the string. The CG
    pattern is a SMARTS (that is what generation parsed it as), and
    utils.extract_elements is a SMILES tokenizer by its own docstring: it yields
    nothing at all for '[#6][#7]', and three symbols for the two atoms of
    'C[N,O]'. RDKit gets both right.

    ``None`` marks an atom whose element the pattern does not pin down -- an OR
    or negation query, which RDKit reports as atomic number 0. Callers must treat
    those slots as matching any element; only the query fragment itself knows
    what is really there.

    Depends only on the CG definition, so it is cached rather than recomputed per
    permutation.
    """
    mol = mol_from_fragment(cg_smarts)
    if mol is None:
        raise ValueError(f"Could not parse CG pattern as SMARTS: {cg_smarts!r}")
    return tuple(a.GetSymbol() if a.GetAtomicNum() else None for a in mol.GetAtoms())

def get_query_cg_coords(sub, cg_smarts):
    coords_list = []
    mol_elements = []
    conf = sub.GetConformer()  # Get the 3D conformer to get coords
    for atom in sub.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())  # returns an RDKit Point3D object
        coords_list.append((pos.x, pos.y, pos.z))
        mol_elements.append(atom.GetSymbol())

    # The guard stays per-Mol: atom order is established separately for every
    # permutation (by Frags, or by reorder_sub_to_target_smiles picking matches[0]),
    # so validating one representative and trusting the rest would not be equivalent.
    cg_elements = cg_element_symbols(cg_smarts)
    if len(mol_elements) != len(cg_elements) or any(
            expected is not None and expected != got
            for expected, got in zip(cg_elements, mol_elements)):
        raise ValueError(f"Element order mismatch between CG pattern and RDKit Mol:\n"
                         f"CG pattern elements: {list(cg_elements)} "
                         f"(None = element not pinned by the pattern)\n"
                         f"Mol elements: {mol_elements}")
    return coords_list
