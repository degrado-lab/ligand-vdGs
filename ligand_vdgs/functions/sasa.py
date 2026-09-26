"""Buried-SASA contact gate: which residues bury a chemical group's surface.
A residue's `buried_area` is CG surface it occludes alone; `shared_area` is its
1/k split of surface occluded jointly with k-1 partners (see `ResidueContact`).
`contact_area = buried_area + shared_area` is the production criterion; one
leave-one-out pass over a Fibonacci-sphere point set computes all of it at once."""
# Heavy atoms only (H is neither surface nor occluder), radii fitted to Probe's own
# output so the recall comparison is meaningful. A point occluded by exactly one
# residue is that residue's exclusive `buried_area`; a point occluded by k>=1
# *reportable* residues jointly (and not also covered by an excluded/ligand atom,
# already inaccessible in the ligand-alone reference state) is split 1/k into each's
# `shared_area` (`shadowed_area` is the same points unsplit, for diagnostics). This
# matters because a contact whose whole patch is shared has `buried_area == 0` --
# invisible to a threshold on exclusive area alone. Summed over partners, buried_area
# + shared_area exactly partitions the surface accessible on the ligand alone and
# buried in the complex (asserted in tests/test_sasa_gate.py against two independent
# `accessible_area` runs).
import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree

PROBE_RADIUS = 1.4

# Frozen with theta (see module docstring): 512 points -> 0.19 A^2 quantum on
# an oxygen (4*pi*2.8^2/512), well under any usable theta.
N_SPHERE_POINTS = 512

# Membership threshold, A^2, against `contact_area`. THE one number: change the
# calibration only here.
#
# Pre-registered, not a judgement call: on the full mirror (2.67M CG sites,
# 65,593 structures) every contact-area band diverged from all three nulls, so
# "theta at the top of the flat region below the lowest divergent band" is 0.
# `contact_area(c) > theta` is strict, so 0.0 still requires >=1 attributed
# surface point, not "admit everything". A 0.25 A^2 "one point" floor was
# rejected: point area is element-dependent (4*pi*(r+probe)^2/n_points) and
# that floor fails for P/Br/I, letting single-point phosphate contacts
# through -- see `n_points` on `ResidueContact` for the radius-independent
# version. Calibrated water-free; re-sweep when waters return as occluders.
MIN_CONTACT_AREA = 0.0

# Water-bridge legs are polar-atom distance gates (N/O...O), not SASA.
WATER_BRIDGE_DIST = 3.5
POLAR_ELEMENTS = frozenset(('N', 'O'))
WATER_RESNAMES = frozenset(('HOH', 'DOD', 'WAT', 'H2O', 'TIP', 'TIP3', 'SOL'))

# Heavy-atom radii fitted to Probe's own output (scratch/probe_vs_distance_10k/report.txt,
# "used_in_sweep"), so the recall-vs-Probe sweep compares criteria, not radii.
FITTED_RADII = {
    'C': 1.750, 'N': 1.548, 'O': 1.400, 'S': 1.790, 'P': 1.800,
    'F': 1.300, 'CL': 1.770, 'BR': 1.850, 'I': 2.100}

# One pooled 'OTHER' bucket from the fit (1.569 A, +0.26 A residual on OTHER-S --
# an envelope, not a radius). Metals/Se are occluders only today. Bondi 1964 / CRC vdW.
OTHER_RADII = {
    'SE': 1.90, 'B': 1.92, 'SI': 2.10, 'AS': 1.85,
    'LI': 1.82, 'NA': 2.27, 'K': 2.75, 'MG': 1.73, 'CA': 2.31,
    'MN': 2.05, 'FE': 2.05, 'CO': 2.00, 'NI': 1.63, 'CU': 1.40, 'ZN': 1.39,
    'CD': 1.58, 'HG': 1.55, 'PT': 1.75, 'AU': 1.66, 'AG': 1.72, 'MO': 2.10,
    'W': 2.10, 'V': 2.05, 'CR': 2.05, 'PB': 2.02, 'SR': 2.49, 'BA': 2.68}

DEFAULT_RADIUS = 1.80

_TWO_LETTER_ELEMENTS = frozenset(OTHER_RADII) | {'CL', 'BR'}

# Atom names in these residues follow the protein/nucleic convention: the element is
# the leading character, so ``CA`` is C-alpha, not calcium, and ``CD`` is C-delta.
# MSE's selenium is the one two-letter element among them.
_STANDARD_RESNAMES = frozenset((
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEC', 'PYL', 'HID', 'HIE', 'HIP', 'CYX', 'CYM', 'ASH', 'GLH', 'LYN',
    'A', 'C', 'G', 'U', 'T', 'I', 'DA', 'DC', 'DG', 'DT', 'DI', 'DU',
    'HOH', 'DOD', 'WAT'))

def element_of(name, element='', resname=None):
    """Element symbol for one atom, falling back from the (often-blank) element field."""
    # Fallback order: (1) if `resname` is standard, atom names follow the protein/
    # nucleic convention -- element is the leading char (`CA` = C-alpha, not calcium)
    # -- so pass `resname` whenever known. (2) Else use the PDB column convention on
    # the raw 4-char `name`: a non-blank column 13 means a two-letter element, which
    # ProDy's stripped `getNames()` loses.
    elem = (element or '').strip().upper()
    if elem:
        return elem
    raw = name or ''
    if resname is not None and resname.strip().upper() in _STANDARD_RESNAMES:
        stripped = raw.strip().upper()
        if stripped == 'SE':
            return 'SE'
        return stripped[:1] or 'X'
    if len(raw) == 4:
        cand = raw[:2].strip().upper() if raw[0] != ' ' else raw[1:2].strip().upper()
    else:
        cand = raw.strip().upper()[:2]
    cand = ''.join(c for c in cand if c.isalpha())
    if len(cand) == 2 and cand not in _TWO_LETTER_ELEMENTS:
        cand = cand[0]
    return cand or 'X'

def elements_for(names, elements=None, resnames=None):
    """Element symbols for parallel name/element/resname arrays."""
    n = len(names)
    if elements is None:
        elements = [''] * n
    if resnames is None:
        resnames = [None] * n
    return [element_of(nm, el, rn) for nm, el, rn in zip(names, elements, resnames)]

def radius_of(elem):
    """vdW radius for an element symbol; `DEFAULT_RADIUS` for anything unrecognised."""
    elem = elem.upper()
    return FITTED_RADII.get(elem) or OTHER_RADII.get(elem, DEFAULT_RADIUS)

def radii_for(elements):
    """Radius array for a sequence of element symbols."""
    return np.array([radius_of(e) for e in elements], dtype=np.float64)

def _fibonacci_sphere(n_points):
    """`n_points` roughly equal-area points on the unit sphere. Deterministic."""
    i = np.arange(n_points, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * i / n_points
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = i * np.pi * (1.0 + 5.0 ** 0.5)
    return np.column_stack((r * np.cos(phi), r * np.sin(phi), z))

class ResidueContact(object):
    """Per-residue contact strength against one CG copy (see module docstring for
    buried/shared/contact_area). `shadowed_area` is the unsplit shared area (double-
    counts by construction; a recall sweep stratifies "how shared" on it) and
    `partner_shadowed_area` is the part of it shared with another *reportable*
    residue rather than the ligand's own atoms."""
    # n_co_occluders: other residues sharing occluded points with this one.
    # n_atom_pairs: heavy-atom pairs (CG atom, residue atom) inside the prefilter.
    # min_heavy_dist: closest CG-atom-to-residue-atom heavy distance, A.
    # n_points: CG surface points credited to this residue, the radius-independent
    # twin of contact_area -- a (CG atom, point index) pair, so the same point index
    # on two CG atoms counts twice.
    __slots__ = ('buried_area', 'shared_area', 'shadowed_area',
                 'partner_shadowed_area', 'n_co_occluders', 'n_atom_pairs',
                 'min_heavy_dist', 'n_points')

    def __init__(self, buried_area, shadowed_area, n_co_occluders,
                 n_atom_pairs, min_heavy_dist, partner_shadowed_area=0.0,
                 shared_area=0.0, n_points=0):
        self.buried_area = float(buried_area)
        self.shared_area = float(shared_area)
        self.shadowed_area = float(shadowed_area)
        self.partner_shadowed_area = float(partner_shadowed_area)
        self.n_co_occluders = int(n_co_occluders)
        self.n_atom_pairs = int(n_atom_pairs)
        self.min_heavy_dist = float(min_heavy_dist)
        self.n_points = int(n_points)

    def __repr__(self):
        return ('ResidueContact(buried_area={:.3f}, shared_area={:.3f}, '
                'shadowed_area={:.3f}, n_co_occluders={}, n_atom_pairs={}, '
                'n_points={}, min_heavy_dist={:.3f})'.format(
                    self.buried_area, self.shared_area, self.shadowed_area,
                    self.n_co_occluders, self.n_atom_pairs, self.n_points,
                    self.min_heavy_dist))

def contact_area(contact):
    """The quantity membership is thresholded on: exclusive + shared burial."""
    # Exclusive `buried_area` alone has a hard recall ceiling against Probe: a residue
    # whose whole occluded patch is also covered by another residue is credited zero
    # by leave-one-out at any threshold. Adding the 1/k shared share is what recovers
    # those contacts (see bench/theta_sweep.log).
    return contact.buried_area + contact.shared_area

def candidate_reach(radii, probe_radius=PROBE_RADIUS):
    """Largest centre-to-centre distance at which any atom can bury CG surface.
    A probe of radius ``p`` rolling on atom i's surface reaches r_i + 2p from
    i's centre, so an atom j beyond r_i + r_j + 2p buries no area of i."""
    r_max = float(radii.max()) if len(radii) else DEFAULT_RADIUS
    return r_max + r_max + 2.0 * probe_radius

def buried_area_by_residue(coords, radii, resindices, cg_atom_idxs,
                           exclude_resindices=(), n_points=N_SPHERE_POINTS,
                           probe_radius=PROBE_RADIUS, tree=None,
                           return_free_area=False):
    """Leave-one-out buried CG surface area per residue, in one pass. Returns
    resindex -> ResidueContact for every candidate residue (zero-burial rows
    included); with `return_free_area`, a `(dict, cg_free_sasa)` pair instead,
    where `cg_free_sasa` is how much CG surface was exposed to be buried (vs.
    the credited areas summing to how much partners actually took)."""
    # coords/radii: heavy-atom-only (n_atoms,[3]) arrays (radii via `radii_for`).
    # resindices: residue grouping, e.g. ProDy's getResindices(). cg_atom_idxs:
    # indices of the CG's atoms; everything else occludes them. exclude_resindices:
    # residues that still occlude but are never reported as partners -- the CG's own
    # residue and the rest of a multi-residue ligand. tree: prebuilt cKDTree over
    # `coords`, shared across a structure's CG copies.
    coords = np.asarray(coords, dtype=np.float64)
    radii = np.asarray(radii, dtype=np.float64)
    resindices = np.asarray(resindices)
    cg_atom_idxs = np.asarray(list(cg_atom_idxs), dtype=np.int64)
    if not len(cg_atom_idxs) or not len(coords):
        return ({}, 0.0) if return_free_area else {}
    if tree is None:
        tree = cKDTree(coords)
    excluded = set(int(r) for r in exclude_resindices)

    _unit = _fibonacci_sphere(n_points)
    _reach = candidate_reach(radii, probe_radius)
    _cg_free = 0.0

    sole, shared, weighted, partner_shared = (defaultdict(float) for _ in range(4))
    n_pts, n_pairs, co_occ = defaultdict(int), defaultdict(int), defaultdict(set)
    min_dist = defaultdict(lambda: np.inf)

    for cg_i in cg_atom_idxs:
        r_i = radii[cg_i]
        sphere_r = r_i + probe_radius
        point_area = 4.0 * np.pi * sphere_r * sphere_r / n_points

        nbr = np.asarray(tree.query_ball_point(coords[cg_i], _reach), dtype=np.int64)
        nbr = nbr[nbr != cg_i]
        d = np.linalg.norm(coords[nbr] - coords[cg_i], axis=1)
        # Exact prefilter: only atoms whose probe-expanded sphere intersects this
        # CG atom's can occlude any of its points, and only those count as atom pairs.
        keep = d <= (r_i + radii[nbr] + 2.0 * probe_radius)
        nbr, d = nbr[keep], d[keep]

        nbr_res = resindices[nbr]
        for rid, dist in zip(nbr_res, d):
            rid = int(rid)
            if rid in excluded:
                continue
            n_pairs[rid] += 1
            min_dist[rid] = min(min_dist[rid], float(dist))

        _points = coords[cg_i] + sphere_r * _unit
        # (n_points, n_nbr): is this point inside neighbour j's probe-expanded sphere?
        _dists = np.linalg.norm(_points[:, None, :] - coords[nbr][None, :, :], axis=2)
        occluded = _dists < (radii[nbr] + probe_radius)[None, :]
        if not occluded.any():
            # Also covers "no candidate atoms at all": empty `nbr` makes `occluded`
            # zero-width, so `.any()` is False here too -- one exit for both.
            _cg_free += point_area * n_points
            continue

        # Collapse atoms to residues: a point is occluded by residue R when any of
        # R's atoms occludes it. Residue count per point drives the attribution.
        res_ids = np.unique(nbr_res)
        per_res = np.array([occluded[:, nbr_res == rid].any(axis=1) for rid in res_ids])
        only_one = per_res.sum(axis=0) == 1
        reportable = np.array([int(r) not in excluded for r in res_ids])
        n_occ_partner = (per_res[reportable].sum(axis=0) if reportable.any()
                         else np.zeros(n_points, dtype=np.int64))
        # A point the ligand's own atoms already cover is inaccessible in the
        # ligand-alone reference state, so no partner buries it and leave-one-out
        # gives it to nobody -- it leaves the shared pool entirely.
        ligand_occluded = (per_res[~reportable].any(axis=0) if (~reportable).any()
                           else np.zeros(n_points, dtype=bool))
        _cg_free += point_area * int(np.count_nonzero(~ligand_occluded))

        for k, rid in enumerate(res_ids):
            rid = int(rid)
            if rid in excluded or not per_res[k].any():
                continue
            mask = per_res[k]
            n_sole = int(np.count_nonzero(mask & only_one))
            n_shared = int(np.count_nonzero(mask) - n_sole)
            if n_sole:
                sole[rid] += n_sole * point_area
                n_pts[rid] += n_sole
            if n_shared:
                shared_pts = mask & ~only_one
                shared[rid] += n_shared * point_area
                # 1/k of each jointly occluded point, over the k PARTNER residues
                # occluding it, excluding ligand-covered points -- summed over
                # partners this keeps buried_area + shared_area coherent with
                # buried_area's own reference state.
                credit_pts = shared_pts & ~ligand_occluded & (n_occ_partner >= 2)
                if credit_pts.any():
                    weighted[rid] += point_area * float((1.0 / n_occ_partner[credit_pts]).sum())
                    # Exactly the points that contributed area, counted rather than
                    # area-weighted, so n_points > 0 iff contact_area > 0.
                    n_pts[rid] += int(np.count_nonzero(credit_pts))
                with_partner = shared_pts & (n_occ_partner >= 2)
                if with_partner.any():
                    partner_shared[rid] += int(np.count_nonzero(with_partner)) * point_area
                _partner_res = res_ids[per_res[:, shared_pts].any(axis=1)]
                co_occ[rid].update(int(r) for r in _partner_res if int(r) != rid)

    out = {rid: ResidueContact(sole[rid], shared[rid], len(co_occ[rid]), n_pairs[rid],
                               min_dist[rid], partner_shared[rid], weighted[rid], n_pts[rid])
           for rid in n_pairs}
    return (out, _cg_free) if return_free_area else out
