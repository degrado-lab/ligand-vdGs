# vdg_struct_utils.py
#
# General ProDy structural extraction utilities shared across the vdG pipeline
# (library generation, hit-finding, analysis tools).

import numpy as np

_NAN3 = np.array([np.nan, np.nan, np.nan], dtype=np.float32)


# ---------------------------------------------------------------------------
# Occupancy protocol
# ---------------------------------------------------------------------------
#
# The B-factor-free occupancy column encodes each atom's role in a vdG PDB, and
# for CG atoms also its slot index (which fixes their order):
#
#   CG atom i          3.00 + 0.01 * i   (band closed at 100 slots; >= 4.0 reserved)
#   vdM residue        2.0
#   other ligand atom  1.0
#
# Readers use only the *order* of the CG occupancies and their distinctness,
# never the value. The vdG-miner submodule is not
# importable and carries its own copy of these numbers -- change both together.
CG_OCC_BASE = 3.0
CG_OCC_STEP = 0.01
CG_OCC_CAPACITY = 100  # slots per band: 3.00 .. 3.99
VDM_OCC = 2.0
NONCG_LIGAND_OCC = 1.0
# Half a step outside the end slots, so float error in a written-then-parsed
# occupancy cannot move slot 0 or the last slot out of the band.
_CG_OCC_MIN = CG_OCC_BASE - CG_OCC_STEP / 2.0
_CG_OCC_MAX = CG_OCC_BASE + (CG_OCC_CAPACITY - 0.5) * CG_OCC_STEP


def cg_slot_occupancy(slot_index):
    """Occupancy encoding CG slot ``slot_index``. Raises past the top of the band."""
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
    """CG atoms in slot order, or None if their occupancies do not encode one.

    A shared slot is not a tie to break arbitrarily -- it means the encoding was
    lost (duplicated slot upstream, or a PDB round-trip that rounded two slots
    together). Rounding matches the PDB occupancy column's two decimals.
    """
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
    # Overwrite the flanking residue at chain_break_ind and every position beyond
    # it (in that direction) with `label`: past a break these are not the vdM's
    # sequence neighbours. chain_break_ind is +/-1..+/-(num_flanking+1); the
    # outermost value overwrites nothing. 0 is never passed.
    # Mutates and returns flanking_seq_dict.
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
    """Precompute the residues ``get_AA_and_CA_coords`` can answer without a select().

    Each ProDy selection costs ~250-300 us of fixed overhead, and the flank walk
    makes three per flank residue (~8 ms per environment at ``--flank 2``).

    Only residues with an unambiguous answer are indexed: all atoms flagged
    protein, exactly one resname, exactly one atom named CA. Everything else --
    non-protein, mixed resnames, no CA, altloc duplicates, unknown resindices
    (incl. the negative ones the flank walk produces at chain edges) -- is absent,
    so the caller falls through to the selection path with its warnings and altloc
    handling. Returns None if the atomgroup cannot be indexed, which is also a 
    fall-through.
    """
    try:
        resindices = prody_obj.getResindices()
        resnames = prody_obj.getResnames()
        names = prody_obj.getNames()
        coords = prody_obj.getCoords()
        # Resname-based, and verified to agree with `sel.protein is None` on
        # MSE/SEP/TPO/UNK/HOH/nucleic as well as the standard 20.
        protein = prody_obj.getFlags('protein')
    except Exception:
        return None
    if any(a is None for a in (resindices, resnames, names, coords, protein)):
        return None

    is_ca = names == 'CA'
    # Group by resindex through one sort rather than a per-residue mask over the
    # whole atomgroup, which would be O(atoms x residues).
    order = np.argsort(resindices, kind='stable')
    unique_resindices = np.unique(resindices)
    starts = np.searchsorted(resindices[order], unique_resindices, side='left')
    ends = np.searchsorted(resindices[order], unique_resindices, side='right')

    index = {}
    for resindex, start, end in zip(unique_resindices, starts, ends):
        rows = order[start:end]
        residue_resnames = set(resnames[rows])
        ca_rows = rows[is_ca[rows]]
        # len(ca_rows) != 1: none, or altloc copies to resolve by occupancy.
        if not protein[rows].all() or len(residue_resnames) != 1 or len(ca_rows) != 1:
            continue
        index[int(resindex)] = (residue_resnames.pop(),
                                np.asarray(coords[ca_rows[0]], dtype=np.float32))
    return index


def get_AA_and_CA_coords(prody_obj, current_resindex, flank_index=None):
    '''
    Return (AA, CA_coords) for a residue index.
    If the residue is missing, non-protein, ambiguous, or has no usable CA,
    returns AA=FLANK_MISSING and CA_coords = [nan, nan, nan]. That marker is
    deliberately not NONCANONICAL_AA_LABEL ('X').

    ``flank_index`` is an optional build_flank_lookup_index() result; a miss in
    it falls through to the selection path.
    '''
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
    """Backbone carbonyl O of one residue, float32 (3,), NaN-filled when absent.

    Not a Stage-1 atom: it is stored beside N/CA/C (``nr_vdm_o_coords``) so
    contact statistics can see the acceptor, and never enters an RMSD."""
    try:
        atom_obj = obj.select('name O')
        if atom_obj is None or len(atom_obj) == 0:
            return np.full(3, np.nan, dtype=np.float32)
        return np.asarray(_pick_best_altloc(atom_obj).getCoords(),
                          dtype=np.float32).reshape(3)
    except Exception:
        return np.full(3, np.nan, dtype=np.float32)


def get_cg_atoms(prody_obj, pdbpath):
    """Ordered CG atoms of one vdG.

    Returns ``(coords, names, elements, seg, chain, resnum, resname)``. The last
    four are *scalars*: a CG is a substructure of a single ligand residue, so a
    CG whose atoms disagree is an upstream bug (a SMARTS match straddling two
    residues, or a mis-set slot occupancy) and is skipped, not averaged away.
    """
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


# ---------------------------------------------------------------------------
# Residue-slot classification
# ---------------------------------------------------------------------------
#
# Canonical heavy-atom names of the 20 standard residues (mirrors
# vdg_miner/constants.py `protein_atoms` with hydrogens dropped.
#
# The point of this table is that a resname cannot be trusted on its own: a GFP
# chromophore is deposited as a residue named GLY carrying a fused imidazolinone
# ring, so ProDy calls it a GLY with a non-empty `sidechain`.
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

# C-terminal carboxylate oxygen (and depositors' alternate names): canonical, but
# neither backbone nor sidechain, matching ProDy's own split.
_TERMINAL_HEAVY_ATOMS = {'OXT', 'OT1', 'OT2'}

BACKBONE_HEAVY_ATOMS = {'N', 'CA', 'C', 'O'}

# Alternate heavy atoms a canonical resname legitimately carries. Se-Met is
# renamed MSE -> MET (by prepwizard, and by s01 via
# _prep_filters.MODIFIED_RESIDUE_RENAMES) while keeping SE in place of SD; calling
# it non-canonical would discard real MET observations (~51% of every `X` slot in
# the pre-change library). Selenocysteine is renamed SEC -> CYS the same way.
_ALTERNATE_HEAVY_ATOMS = {'MET': {'SE'}, 'CYS': {'SE'}}

# Label for a residue slot whose non-canonical atoms contact the CG. Not a
# resname (the true residue is in nr_scrr_resname). Unreachable by hit finding,
# which is the point.
NONCANONICAL_AA_LABEL = 'X'

# Label for a slot whose *backbone* contacts the CG -- one label covering every
# residue, GLY and PRO included. Their backbones really aren't substitutable, but
# that is a per-vdG property, not a per-bucket one, so the read path tests it
# directly (hit_finder_core.backbone_slot_blockers / backbone_slots_can_host)
# against the query's own sidechain. See docs/aa_bucket_treatment.md for the
# measurements. BB_LABELS stays a set for importers and for a possible future
# backbone label that is not a query-side property. 
BB_LABEL = 'bb'
BB_LABELS = frozenset((BB_LABEL,))


def bb_label_for(resname):
    """Which backbone label a generated slot on `resname` gets."""
    return BB_LABEL


def query_slot_labels(resname):
    """Labels a *query* residue named `resname` may be matched under.

    Counterpart of `bb_label_for` on the read path: every residue may be matched
    as its own sidechain identity or as a backbone contact. Whether a backbone
    geometry is physically hostable is decided per-vdG in
    hit_finder_core.backbone_slots_can_host, not here.

    Glycine's 'GLY' option matches nothing in practice (the write path cannot
    produce it), but is kept so BSR enumeration stays a uniform 2^n; the cost is
    one missing-file lookup per glycine per combo.
    """
    return (resname, BB_LABEL)


# Flanking-sequence markers, deliberately not 'X': a flanking position says
# nothing about chemistry. FLANK_MISSING = the residue is absent, non-protein,
# ambiguous, or has no usable CA (it may be perfectly ordinary, just unreadable).

FLANK_MISSING = '-'
FLANK_CHAIN_BREAK = '!'
FLANK_UNCOMPARABLE = frozenset((FLANK_MISSING, FLANK_CHAIN_BREAK))

# Per-slot provenance, written to the npz as `nr_slot_flag` (int8):
#
#   bits 0-1  which moiety of the *canonical* residue is nearest the CG
#   bit 2     whether the residue carries non-canonical heavy atoms at all
#
# Kept disjoint from the bucket label so nothing is recorded twice: the label
# already says which moiety contacts whenever that is not canonical ('X') or not
# a sidechain (BB_LABELS), leaving the flag free to describe the canonical
# residue underneath -- on an 'X' slot it says what the real sidechain was doing.
#
# SLOT_NO_SC is chemistry on a glycine and missing density anywhere else; since
# the backbone label no longer distinguishes them, check `nr_scrr_resname`
# against 'GLY' rather than pooling. The two are lopsided (253,992 vs 227 slots
# in one measured library), so do not design around the rare one.
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
    """Hydrogen/deuterium test that tolerates a blank element column.

    Falls back to the PDB atom-name convention (optional leading digit, then the
    element letter), common in older or hand-edited files.
    """
    el = (element or '').strip().upper()
    if el:
        return el in ('H', 'D')
    return (name or '').strip().lstrip('0123456789')[:1].upper() in ('H', 'D')


def split_residue_heavy_atoms(res_obj, resname):
    """Split a residue's heavy atoms into (backbone, sidechain, non-canonical) coords.

    Decided by atom *name* against CANONICAL_HEAVY_ATOMS, not ProDy's
    backbone/sidechain flags, which call a GFP chromophore's fused ring a
    sidechain. One table drives both
    questions, so an atom cannot be canonical for one and not the other. A
    canonical name is accepted whatever element the file declares (the element
    column is blank or wrong often enough that gating on it would relabel
    ordinary residues as non-canonical).

    Returns (None, None, None) when `resname` is outside the 20.
    Hydrogens are dropped here, so a stray H is never reported as non-canonical.
    """
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
