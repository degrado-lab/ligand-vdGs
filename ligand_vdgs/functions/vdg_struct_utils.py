# vdg_struct_utils.py
#
# General ProDy structural extraction utilities shared across the vdG pipeline
# (library generation, hit-finding, analysis tools).

import numpy as np

_NAN3 = np.array([np.nan, np.nan, np.nan], dtype=np.float32)


# Occupancy protocol: CG slot i = 3.00 + .01*i (100 slots); vdM = 2, other ligand = 1.
# Keep in sync with the unimportable vdG-miner copy; readers rely on CG order/distinctness.
CG_OCC_BASE = 3.0
CG_OCC_STEP = 0.01
CG_OCC_CAPACITY = 100  # slots per band: 3.00 .. 3.99
VDM_OCC = 2.0
NONCG_LIGAND_OCC = 1.0
# Half-step margins tolerate PDB round trips without admitting adjacent bands.
_CG_OCC_MIN = CG_OCC_BASE - CG_OCC_STEP / 2.0
_CG_OCC_MAX = CG_OCC_BASE + (CG_OCC_CAPACITY - 0.5) * CG_OCC_STEP


def cg_slot_occupancy(slot_index):
    """Encode a CG slot occupancy; reject indices outside the 100-slot band."""
    if not 0 <= slot_index < CG_OCC_CAPACITY:
        raise ValueError(
            f'CG slot {slot_index} is outside the occupancy band '
            f'(0..{CG_OCC_CAPACITY - 1}); a CG this large cannot be encoded.')
    return CG_OCC_BASE + slot_index * CG_OCC_STEP


def select_cg_atoms(prody_obj):
    """The CG atoms of a vdG structure, unordered. None if there are none."""
    cg = prody_obj.select(f'occupancy > {_CG_OCC_MIN} and occupancy < {_CG_OCC_MAX}')
    return cg if cg is not None and len(cg) > 0 else None


def sort_cg_atoms_by_slot(cg):
    """Sort CG atoms by slot; duplicate two-decimal occupancies make order unknown."""
    atoms = sorted(cg, key=lambda a: a.getOccupancy())
    if len({round(float(a.getOccupancy()), 2) for a in atoms}) != len(atoms):
        return None
    return atoms


def get_res_iden(vdm_obj):
    seg, chain, resnum, resname = (list(set(vals)) for vals in (
        vdm_obj.getSegnames(), vdm_obj.getChids(),
        vdm_obj.getResnums(), vdm_obj.getResnames()))
    if max(map(len, (seg, chain, resnum, resname))) != 1:
        print(f'[WARNING] get_res_iden: expected single-residue object, got '
              f'segs={seg}, chains={chain}, resnums={resnum}, resnames={resname}. Skipping.')
        return None
    return seg[0], chain[0], resnum[0], resname[0]


def found_chain_break(flanking_seq_dict, chain_break_ind, label=None):
    # A break invalidates that flank and every farther position in its direction.
    if label is None:
        label = FLANK_CHAIN_BREAK
    if chain_break_ind == 0:
        return flanking_seq_dict  # should never be reached
    beyond = ((lambda n: n >= chain_break_ind) if chain_break_ind > 0
              else (lambda n: n <= chain_break_ind))
    for ind in [n for n in flanking_seq_dict if beyond(n)]:
        flanking_seq_dict[ind] = [label, _NAN3.astype(np.float64)]
    return flanking_seq_dict


def _pick_best_altloc(atom_sel):
    """One atom from a multi-altloc selection: highest occupancy, then altloc 'A'."""
    if len(atom_sel) == 1:
        return atom_sel[0]
    atoms = list(atom_sel)
    max_occ = max(a.getOccupancy() for a in atoms)
    top = [a for a in atoms if a.getOccupancy() == max_occ]
    return next((a for a in top if a.getAltloc() == 'A'), top[0])


def is_valid_backbone_coords(coords):
    """Whether coords are a finite, non-collinear N/CA/C triplet."""
    try:
        coords = np.asarray(coords, dtype=np.float32)
        if coords.shape != (3, 3) or not np.isfinite(coords).all():
            return False
        return np.linalg.matrix_rank(coords - coords.mean(axis=0, keepdims=True)) >= 2
    except (TypeError, ValueError, np.linalg.LinAlgError):
        return False


def build_flank_lookup_index(prody_obj):
    """Index protein residues with one resname and one CA; omit ambiguous residues."""
    try:
        resindices = prody_obj.getResindices()
        resnames = prody_obj.getResnames()
        names = prody_obj.getNames()
        coords = prody_obj.getCoords()
        protein = prody_obj.getFlags('protein')
    except Exception:
        return None
    if any(a is None for a in (resindices, resnames, names, coords, protein)):
        return None

    is_ca = names == 'CA'
    # One sort avoids scanning the full atomgroup once per residue.
    order = np.argsort(resindices, kind='stable')
    unique_resindices = np.unique(resindices)
    starts = np.searchsorted(resindices[order], unique_resindices, side='left')
    ends = np.searchsorted(resindices[order], unique_resindices, side='right')

    index = {}
    for resindex, start, end in zip(unique_resindices, starts, ends):
        rows = order[start:end]
        residue_resnames = set(resnames[rows])
        ca_rows = rows[is_ca[rows]]
        if not protein[rows].all() or len(residue_resnames) != 1 or len(ca_rows) != 1:
            continue
        index[int(resindex)] = (residue_resnames.pop(),
                                np.asarray(coords[ca_rows[0]], dtype=np.float32))
    return index


def get_AA_and_CA_coords(prody_obj, current_resindex, flank_index=None):
    """Return residue name and CA, or the unreadable-flank marker and NaNs."""
    unreadable = (FLANK_MISSING, _NAN3.copy())
    if flank_index is not None:
        hit = flank_index.get(int(current_resindex))
        if hit is not None:
            return hit[0], hit[1].copy()  # copy: callers get a fresh array

    # Backticks: ProDy needs them to parse a negative resindex.
    sel_str = (f'resindex `{current_resindex}`' if current_resindex < 0
               else f'resindex {current_resindex}')
    curr_resindex_obj = prody_obj.select(sel_str)
    if curr_resindex_obj is None or curr_resindex_obj.protein is None:
        return unreadable

    curr_res_AA = list(set(curr_resindex_obj.getResnames()))
    if len(curr_res_AA) != 1:
        print(f'[WARNING] get_AA_and_CA_coords: \nResindex {current_resindex} in '
              f'{prody_obj.getTitle()} contains >1 AA: {curr_res_AA}. Marking the flank unreadable.')
        return unreadable

    CA_sel = curr_resindex_obj.select('name CA')
    if CA_sel is None or len(CA_sel) == 0:
        return unreadable
    try:
        coords = np.asarray(_pick_best_altloc(CA_sel).getCoords(), dtype=np.float32)
    except Exception:
        print(f'[WARNING] Could not resolve CA for resindex {current_resindex} in '
              f'{prody_obj.getTitle()}. Marking the flank unreadable.')
        return unreadable
    return curr_res_AA[0], coords


def get_bb_coords(obj):
    bb_coords = []
    for atom_name in ('N', 'CA', 'C'):
        try:
            atom_obj = obj.select(f'name {atom_name}')
            if atom_obj is None or len(atom_obj) == 0:
                return None
            bb_coords.append(np.asarray(_pick_best_altloc(atom_obj).getCoords(),
                                        dtype=np.float32))
        except Exception:
            return None
    return bb_coords if is_valid_backbone_coords(bb_coords) else None


def get_bb_o_coords(obj):
    """Return a residue's carbonyl O, or float32 NaNs when absent."""
    try:
        atom_obj = obj.select('name O')
        if atom_obj is None or len(atom_obj) == 0:
            return np.full(3, np.nan, dtype=np.float32)
        return np.asarray(_pick_best_altloc(atom_obj).getCoords(),
                          dtype=np.float32).reshape(3)
    except Exception:
        return np.full(3, np.nan, dtype=np.float32)


def get_cg_atoms(prody_obj, pdbpath):
    """Return ordered CG data and its single-residue identity, or None if invalid."""
    cg = select_cg_atoms(prody_obj)
    if cg is None:
        return None
    sorted_atoms = sort_cg_atoms_by_slot(cg)
    if sorted_atoms is None:
        print(f'[WARNING] get_cg_atoms: CG atoms of {pdbpath} do not carry distinct '
              'slot occupancies; their order is undefined. Skipping.')
        return None

    coords = [a.getCoords() for a in sorted_atoms]
    names = [a.getName() for a in sorted_atoms]
    elements = [a.getElement() for a in sorted_atoms]
    residues = {(a.getSegname(), a.getChid(), int(a.getResnum()), a.getResname())
                for a in sorted_atoms}
    if len(residues) != 1:
        print(f'[WARNING] get_cg_atoms: CG of {pdbpath} spans {len(residues)} residues '
              f'{sorted(residues)}; a CG must lie within one ligand residue. Skipping.')
        return None
    seg, chain, resnum, resname = residues.pop()
    return (np.asarray(coords, dtype=np.float32),
            names, elements, seg, chain, resnum, resname)


def get_res_AA_identity(res_obj):
    resnames = list(set(res_obj.getResnames()))
    if len(resnames) != 1:
        print(f'[WARNING] get_res_AA_identity: expected single resname, got {resnames}. Skipping.')
        return None
    return resnames[0]


# Canonical heavy-atom names (mirrors vdg_miner/constants.py without hydrogens).
# Atom names, rather than residue labels or ProDy flags, identify modified chemistry.
CANONICAL_HEAVY_ATOMS = {
    'ALA': {'C', 'CA', 'CB', 'N', 'O'},
    'ARG': {'C', 'CA', 'CB', 'CD', 'CG', 'CZ', 'N', 'NE', 'NH1', 'NH2', 'O'},
    'ASN': {'C', 'CA', 'CB', 'CG', 'N', 'ND2', 'O', 'OD1'},
    'ASP': {'C', 'CA', 'CB', 'CG', 'N', 'O', 'OD1', 'OD2'},
    'CYS': {'C', 'CA', 'CB', 'N', 'O', 'SG'},
    'GLN': {'C', 'CA', 'CB', 'CD', 'CG', 'N', 'NE2', 'O', 'OE1'},
    'GLU': {'C', 'CA', 'CB', 'CD', 'CG', 'N', 'O', 'OE1', 'OE2'},
    'GLY': {'C', 'CA', 'N', 'O'},
    'HIS': {'C', 'CA', 'CB', 'CD2', 'CE1', 'CG', 'N', 'ND1', 'NE2', 'O'},
    'ILE': {'C', 'CA', 'CB', 'CD1', 'CG1', 'CG2', 'N', 'O'},
    'LEU': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CG', 'N', 'O'},
    'LYS': {'C', 'CA', 'CB', 'CD', 'CE', 'CG', 'N', 'NZ', 'O'},
    'MET': {'C', 'CA', 'CB', 'CE', 'CG', 'N', 'O', 'SD'},
    'PHE': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N', 'O'},
    'PRO': {'C', 'CA', 'CB', 'CD', 'CG', 'N', 'O'},
    'SER': {'C', 'CA', 'CB', 'N', 'O', 'OG'},
    'THR': {'C', 'CA', 'CB', 'CG2', 'N', 'O', 'OG1'},
    'TRP': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE2', 'CE3', 'CG', 'CH2', 'CZ2',
            'CZ3', 'N', 'NE1', 'O'},
    'TYR': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N', 'O',
            'OH'},
    'VAL': {'C', 'CA', 'CB', 'CG1', 'CG2', 'N', 'O'},
}

# Terminal carboxylate oxygens are canonical but outside backbone/sidechain.
_TERMINAL_HEAVY_ATOMS = {'OXT', 'OT1', 'OT2'}

BACKBONE_HEAVY_ATOMS = {'N', 'CA', 'C', 'O'}

# Se-Met and selenocysteine are renamed MET/CYS while retaining their SE atom.
_ALTERNATE_HEAVY_ATOMS = {'MET': {'SE'}, 'CYS': {'SE'}}

# Label for slots whose non-canonical atoms contact the CG; true resname is stored separately.
NONCANONICAL_AA_LABEL = 'X'

# One backbone label for all residues; hostability is checked per vdG against query virtual CB/Pro N.
BB_LABEL = 'bb'
BB_LABELS = frozenset((BB_LABEL,))


def bb_label_for(resname):
    """Return the generated backbone-contact label."""
    return BB_LABEL


def query_slot_labels(resname):
    """Return the identity and backbone labels allowed for a query residue."""
    return (resname, BB_LABEL)


# Flank markers describe unreadability, not chemistry, so they differ from 'X'.

FLANK_MISSING = '-'
FLANK_CHAIN_BREAK = '!'
FLANK_UNCOMPARABLE = frozenset((FLANK_MISSING, FLANK_CHAIN_BREAK))

# nr_slot_flag: bits 0-1 encode the nearest canonical moiety; bit 2 marks modified residues.
SLOT_SC          = 0      # sidechain is the closest canonical moiety
SLOT_NO_SC       = 1      # no sidechain heavy atoms present
SLOT_BB_CLOSER   = 2      # sidechain present, but a backbone atom is closer
SLOT_REASON_MASK = 0b011
SLOT_MODIFIED    = 0b100  # OR'd in when the residue carries non-canonical atoms


def slot_reason(flag):
    """The SLOT_SC / SLOT_NO_SC / SLOT_BB_CLOSER part of a slot flag."""
    return int(flag) & SLOT_REASON_MASK


def slot_is_modified(flag):
    """True if the slot's residue carries heavy atoms its resname does not have."""
    return bool(int(flag) & SLOT_MODIFIED)


def is_hydrogen(name, element):
    """Identify H/D from the element, falling back to the PDB atom name."""
    el = (element or '').strip().upper()
    if el:
        return el in ('H', 'D')
    return (name or '').strip().lstrip('0123456789')[:1].upper() in ('H', 'D')


def split_residue_heavy_atoms(res_obj, resname):
    """Split heavy atom coordinates into backbone, sidechain and non-canonical groups."""
    canonical = CANONICAL_HEAVY_ATOMS.get(resname)
    if canonical is None:
        return None, None, None
    allowed = (canonical | _TERMINAL_HEAVY_ATOMS
               | _ALTERNATE_HEAVY_ATOMS.get(resname, frozenset()))

    names = res_obj.getNames()
    coords = np.asarray(res_obj.getCoords(), dtype=np.float32).reshape(-1, 3)
    heavy = np.array([not is_hydrogen(n, e)
                      for n, e in zip(names, res_obj.getElements())],
                     dtype=bool).reshape(-1)
    known = heavy & np.array([n in allowed for n in names], dtype=bool).reshape(-1)
    is_bb = np.array([n in BACKBONE_HEAVY_ATOMS for n in names],
                     dtype=bool).reshape(-1)
    is_term = np.array([n in _TERMINAL_HEAVY_ATOMS for n in names],
                       dtype=bool).reshape(-1)
    return coords[known & is_bb], coords[known & ~is_bb & ~is_term], coords[heavy & ~known]

def virtual_cb(bb):
    """Ideal CB (ProteinMPNN constants) from (..., 3, 3) N/CA/C backbone coords."""
    bb = np.asarray(bb, np.float64)
    b, c = bb[..., 1, :] - bb[..., 0, :], bb[..., 2, :] - bb[..., 1, :]
    return (bb[..., 1, :] - 0.58273431 * np.cross(b, c) + 0.56802827 * b - 0.54067466 * c).astype(np.float32)

def vcb_slot_blockers(bsr_bb, slot_labels, slot_resnames):
    """Side-chain-free slot blockers for bb-labelled slots: the virtual CB (none for Gly) and the
    Pro N. bsr_bb: (n_slots, 3, 3) N/CA/C. Returns [(blocker_atoms, pro_n or None)] or None."""
    return [(virtual_cb(bb)[None] if rn != "GLY" else np.empty((0, 3), np.float32),
             np.asarray(bb[0], np.float32) if rn == "PRO" else None)
            for bb, label, rn in zip(bsr_bb, slot_labels, map(str, slot_resnames))
            if label in BB_LABELS] or None

def receptor_bb_vcb_coords(struct):
    """(P, 3) N/CA/C/O of every protein residue plus the virtual CB of every non-Gly residue with
    N/CA/C; the side-chain-free receptor for hit_finder_core `min_bb_dist`."""
    bb = struct.select("protein and name N CA C O")
    if bb is None:
        return np.empty((0, 3), np.float32)
    return np.concatenate([bb.getCoords().astype(np.float32), virtual_cb(np.reshape(
        [[a.getCoords() for a in atoms] for atoms in ([r.getAtom(n) for n in ("N", "CA", "C")]
                                                      for r in bb.getHierView().iterResidues()
                                                      if r.getResname() != "GLY") if None not in atoms],
        (-1, 3, 3)))])

def bsr_contact_atoms(bsr_bb, slot_resnames):
    """(P, 3) side-chain-free contact atoms of a BSR combo: N/CA/C plus the virtual CB of non-Gly slots."""
    bsr_bb = np.asarray(bsr_bb, np.float32)
    return np.concatenate((bsr_bb.reshape(-1, 3), virtual_cb(bsr_bb[[str(r) != "GLY" for r in slot_resnames]])))
