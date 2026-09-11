"""Buried-SASA contact gate: which residues bury a chemical group's surface.

A residue R contacts the CG when it buries more than ``theta`` A^2 of the CG
atoms' SASA (all-atom SASA minus SASA with R removed). Replaces Probe/distance
cutoffs: zero for a residue hidden behind another, non-zero for a Met SD 5.5 A
off a ring face.

Heavy atoms only (H is neither surface nor occluder), with radii fitted to
Probe's own output so the recall comparison is meaningful.

**Leave-one-out in one pass.** A point becomes accessible when R is deleted
iff every atom occluding it belongs to R, so one pass recording occluders per
point yields every residue's buried area at once (`buried_area(R)` = area
occluded by >=1 atom of R and none outside R) -- exactly equal to the N+1-run
version, checked in `tests/test_sasa_gate.py`. Points occluded by >=2
residues belong to nobody under leave-one-out; `shared_area` splits each such
point 1/k among its k occluding **partner** residues (excluding points any
ligand atom also covers, since those are already inaccessible in the
ligand-alone reference state and so credited to nobody); `shadowed_area` is
the unsplit total. This matters because a contact with a fully-shared patch
has `buried_area == 0`, invisible to any threshold on exclusive area alone.
Summed over partners, `buried_area + shared_area` exactly partitions the
surface accessible on the ligand alone and buried in the complex:

    sum over partners of (buried_area + shared_area)
        == area accessible with only the CG's residue present
           minus area accessible with everything present

(asserted in `tests/test_sasa_gate.py` against two `accessible_area` calls).

Point set: a deterministic Fibonacci sphere of `N_SPHERE_POINTS` points,
quantising area to `4*pi*(r+probe)^2 / N_SPHERE_POINTS`. theta is only
meaningful at a fixed point count -- **frozen together**.
"""
import numpy as np
from scipy.spatial import cKDTree

PROBE_RADIUS = 1.4

# Frozen with theta (see module docstring). 512 points -> 0.19 A^2 quantum on
# an oxygen (4*pi*2.8^2/512), well under any usable theta.
N_SPHERE_POINTS = 512

# Membership threshold, A^2, against `contact_area` below. THE one number: change
# the calibration only here.
#
# CALIBRATED (pre-registered rule, not a judgement call). On the full mirror (job
# 4914928: 2.67M CG sites, 65,593 structures) every contact-area band diverged
# under all three nulls, so "theta at the top of the flat region below the lowest
# divergent band" yields 0. Other instruments agree: angular specificity is
# isotropic and partner order is flat throughout.
#
# 0 is the RECOVERABLE direction: everything above theta keeps buried_area/
# shared_area for stricter filtering later; below theta is never mined and needs
# a rebuild to recover. `contact_area(c) > theta` is strict, so 0.0 still requires
# >=1 attributed surface point, not "admit everything in the prefilter".
#
# Rejected: a 0.25 A^2 "one surface point" floor -- element-dependent since point
# area is 4*pi*(r+probe)^2/n_points; holds for C/N/O/F/Cl/S but FAILS for P/Br/I
# (0.2513/0.2592/0.3007), letting single-point phosphate contacts through. See the
# n_points field for the radius-independent version.
#
# Calibrated water-free; re-sweep when waters return as occluders.
MIN_CONTACT_AREA = 0.0

# Water-bridge legs are polar-atom distance gates (N/O...O), not SASA. Waters are
# absent from the current parent database, so only tests exercise this today.
WATER_BRIDGE_DIST = 3.5
POLAR_ELEMENTS = frozenset(('N', 'O'))
WATER_RESNAMES = frozenset(('HOH', 'DOD', 'WAT', 'H2O', 'TIP', 'TIP3', 'SOL'))

# Heavy-atom radii fitted to Probe's output (scratch/probe_vs_distance_10k/report.txt,
# "used_in_sweep"), so the recall-vs-Probe sweep compares criteria, not radii.
FITTED_RADII = {
    'C': 1.750, 'N': 1.548, 'O': 1.400, 'S': 1.790, 'P': 1.800,
    'F': 1.300, 'CL': 1.770, 'BR': 1.850, 'I': 2.100,
}

# One pooled 'OTHER' bucket from the fit (1.569 A, +0.26 A residual on OTHER-S --
# an envelope, not a radius). Metals/Se are occluders only today; re-examine once
# the metal/second-ligand rule makes them partners. Bondi 1964 / CRC vdW.
OTHER_RADII = {
    'SE': 1.90, 'B': 1.92, 'SI': 2.10, 'AS': 1.85,
    'LI': 1.82, 'NA': 2.27, 'K': 2.75, 'MG': 1.73, 'CA': 2.31,
    'MN': 2.05, 'FE': 2.05, 'CO': 2.00, 'NI': 1.63, 'CU': 1.40, 'ZN': 1.39,
    'CD': 1.58, 'HG': 1.55, 'PT': 1.75, 'AU': 1.66, 'AG': 1.72, 'MO': 2.10,
    'W': 2.10, 'V': 2.05, 'CR': 2.05, 'PB': 2.02, 'SR': 2.49, 'BA': 2.68,
}

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
    'HOH', 'DOD', 'WAT',
))


def element_of(name, element='', resname=None):
    """Element symbol for one atom, falling back from the (often-blank) element field.

    Fallback order: (1) if `resname` is a standard residue, atom names follow the
    protein/nucleic convention -- element is the leading char (`CA` = C-alpha, not
    calcium) -- so pass `resname` whenever known. (2) Else use the PDB column
    convention on the raw 4-char `name` field: a non-blank column 13 means a
    two-letter element, which ProDy's stripped `getNames()` loses.
    """
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
    return [element_of(nm, el, rn)
            for nm, el, rn in zip(names, elements, resnames)]


def radius_of(elem):
    """vdW radius for an element symbol; `DEFAULT_RADIUS` for anything unrecognised."""
    elem = elem.upper()
    if elem in FITTED_RADII:
        return FITTED_RADII[elem]
    return OTHER_RADII.get(elem, DEFAULT_RADIUS)


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
    """Per-residue contact strength against one CG copy.

    buried_area   : A^2 of CG SASA this residue alone occludes (leave-one-out delta,
                    exclusive and exact).
    shared_area   : this residue's 1/k share of points it occludes *jointly* with k-1
                    other partners (points any ligand atom also covers are excluded --
                    already inaccessible in the ligand-alone reference state, same as
                    for `buried_area`). Needed because `buried_area` alone is blind to
                    a contact whose whole patch is shared (3,311/64,299 Probe contacts
                    in the 2,000-structure sweep had `buried_area == 0` with a real
                    shared partner). `buried_area + shared_area` (`contact_area`) is
                    the production criterion (DR-3).
    shadowed_area : full unsplit area of those jointly occluded points -- not a
                    credit, double-counts by construction, but the denominator a
                    recall sweep stratifies "how shared" on.
    partner_shadowed_area :
                    the part of `shadowed_area` shared with another *reportable*
                    residue rather than with the ligand's own atoms (self-occlusion
                    swamps the signal; a second residue is the non-additivity that
                    matters).
    n_co_occluders: number of other residues (incl. excluded) it shares occluded
                    points with.
    n_atom_pairs  : heavy-atom pairs (CG atom, residue atom) inside the prefilter.
    min_heavy_dist: closest CG-atom-to-residue-atom heavy distance, A.
    """

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
        # Number of CG surface points credited to this residue -- the count whose
        # area-weighted form is `contact_area`. Unit is a (CG atom, point index)
        # pair, so the same point index on two CG atoms counts twice.
        #
        # Radius-independent twin of contact_area: a "one point" floor expressed
        # as area is element-dependent (fails for P/Br/I), as a count it isn't.
        # Distinct from `n_atom_pairs` (prefilter atom pairs) and `n_co_occluders`
        # (partner residue count).
        self.n_points = int(n_points)

    def __repr__(self):
        return ('ResidueContact(buried_area={:.3f}, shared_area={:.3f}, '
                'shadowed_area={:.3f}, n_co_occluders={}, n_atom_pairs={}, '
                'n_points={}, min_heavy_dist={:.3f})'.format(
                    self.buried_area, self.shared_area, self.shadowed_area,
                    self.n_co_occluders, self.n_atom_pairs, self.n_points,
                    self.min_heavy_dist))


def contact_area(contact):
    """The quantity membership is thresholded on: exclusive + shared burial.

    Criterion B of the theta sweep; the single definition, called by
    `structure_contacts`, the sweep, and the pi acceptance test, so calibration
    and build can't drift apart.

    Exclusive `buried_area` alone has a hard recall ceiling of 0.9485 against Probe
    (`bench/theta_sweep.log`): a residue whose whole occluded patch is also covered
    by another residue is credited zero by leave-one-out at any threshold (all
    3,311 unreachable pairs were this kind). Adding the 1/k shared share lifts
    recall to 1.000 through theta 7.5. See decision_records.md DR-3.
    """
    return contact.buried_area + contact.shared_area


def candidate_reach(radii, probe_radius=PROBE_RADIUS):
    """Largest centre-to-centre distance at which any atom can bury CG surface.

    A probe of radius ``p`` rolling on atom i's surface reaches r_i + 2p from i's
    centre, so an atom j beyond r_i + r_j + 2p buries no area of i by construction.
    With the fitted radii the worst case is ~6.5 A, which is the candidate prefilter
    radius the pi-recall test was run against.
    """
    r_max = float(radii.max()) if len(radii) else DEFAULT_RADIUS
    return r_max + r_max + 2.0 * probe_radius


def buried_area_by_residue(coords, radii, resindices, cg_atom_idxs,
                           exclude_resindices=(), n_points=N_SPHERE_POINTS,
                           probe_radius=PROBE_RADIUS, tree=None,
                           return_free_area=False):
    """Leave-one-out buried CG surface area per residue, in one pass.

    Parameters
    ----------
    coords, radii : (n_atoms,[3]) float arrays
        Heavy-atom-only (see module docstring) coordinates and vdW radii (`radii_for`).
    resindices : (n_atoms,) int array
        Residue grouping, e.g. ProDy's ``getResindices()``.
    cg_atom_idxs : sequence of int
        Indices of the chemical group's atoms; everything else in the structure
        (including the ligand's own non-CG atoms) occludes them.
    exclude_resindices : container of int
        Residues that still occlude but are never reported as partners -- the CG's
        own residue and the rest of a multi-residue ligand.
    tree : cKDTree, optional
        Prebuilt tree over `coords`, shared across a structure's CG copies.
    return_free_area : bool
        Also return the CG's SASA in the ligand-alone reference state (see Returns).

    Returns
    -------
    dict : resindex -> ResidueContact for every candidate residue, zero-burial rows
    included (`min_heavy_dist` on a zero row says whether the prefilter is wide
    enough). With `return_free_area`, a `(dict, cg_free_sasa)` pair instead:
    `cg_free_sasa` is how much CG surface was EXPOSED to be buried, vs. the credited
    areas summing to how much partners actually TOOK -- a small sum can mean either
    a small CG or an exposed one, and only the ratio distinguishes them.
    """
    coords = np.asarray(coords, dtype=np.float64)
    radii = np.asarray(radii, dtype=np.float64)
    resindices = np.asarray(resindices)
    cg_atom_idxs = np.asarray(list(cg_atom_idxs), dtype=np.int64)
    if not len(cg_atom_idxs) or not len(coords):
        return ({}, 0.0) if return_free_area else {}
    if tree is None:
        tree = cKDTree(coords)
    excluded = set(int(r) for r in exclude_resindices)

    unit = _fibonacci_sphere(n_points)
    reach = candidate_reach(radii, probe_radius)

    # SASA in the ligand-alone reference state (points no EXCLUDED/ligand atom
    # covers) -- the denominator `buried_area` is implicitly measured against, not
    # the buriable total. Free here since `ligand_occluded` is already computed.
    cg_free = 0.0

    # Accumulators keyed by residue index.
    sole = {}        # points occluded by this residue and nothing outside it
    shared = {}      # points it occludes jointly with any other occluder
    weighted = {}    # the same points, each split 1/k over its k occluding residues
    partner_shared = {}   # the subset of those shared with a reportable residue
    n_pts = {}       # credited POINT COUNT: the radius-independent twin of the area
    co_occ = {}      # set of residues it shares occluded points with
    n_pairs = {}
    min_dist = {}

    for cg_i in cg_atom_idxs:
        r_i = radii[cg_i]
        sphere_r = r_i + probe_radius
        point_area = 4.0 * np.pi * sphere_r * sphere_r / n_points

        nbr = np.asarray(tree.query_ball_point(coords[cg_i], reach),
                         dtype=np.int64)
        nbr = nbr[nbr != cg_i]
        if not len(nbr):
            cg_free += point_area * n_points
            continue
        d = np.linalg.norm(coords[nbr] - coords[cg_i], axis=1)
        # Exact prefilter: only atoms whose probe-expanded sphere intersects this
        # CG atom's can occlude any of its points, and only those are counted as
        # atom pairs.
        keep = d <= (r_i + radii[nbr] + 2.0 * probe_radius)
        nbr, d = nbr[keep], d[keep]
        if not len(nbr):
            cg_free += point_area * n_points
            continue

        # Contact-strength bookkeeping over the prefilter set, independent of
        # whether any surface point survives.
        nbr_res = resindices[nbr]
        for rid, dist in zip(nbr_res, d):
            rid = int(rid)
            if rid in excluded:
                continue
            n_pairs[rid] = n_pairs.get(rid, 0) + 1
            if dist < min_dist.get(rid, np.inf):
                min_dist[rid] = float(dist)

        points = coords[cg_i] + sphere_r * unit
        # (n_points, n_nbr): is this point inside neighbour j's probe-expanded sphere?
        dists = np.linalg.norm(points[:, None, :] - coords[nbr][None, :, :],
                               axis=2)
        occluded = dists < (radii[nbr] + probe_radius)[None, :]
        if not occluded.any():
            cg_free += point_area * n_points
            continue

        # Collapse atoms to residues: a point is occluded by residue R when any of
        # R's atoms occludes it. Residue count per point drives the attribution.
        res_ids = np.unique(nbr_res)
        per_res = np.empty((len(res_ids), n_points), dtype=bool)
        for k, rid in enumerate(res_ids):
            per_res[k] = occluded[:, nbr_res == rid].any(axis=1)
        n_occ_res = per_res.sum(axis=0)
        only_one = n_occ_res == 1
        reportable = np.array([int(r) not in excluded for r in res_ids])
        n_occ_partner = (per_res[reportable].sum(axis=0) if reportable.any()
                         else np.zeros(n_points, dtype=np.int64))
        # A point the ligand's own atoms already cover is inaccessible in the
        # ligand-alone reference state, so no partner buries it and leave-one-out
        # gives it to nobody -- it leaves the shared pool entirely.
        ligand_occluded = (per_res[~reportable].any(axis=0) if (~reportable).any()
                           else np.zeros(n_points, dtype=bool))
        cg_free += point_area * int(np.count_nonzero(~ligand_occluded))

        for k, rid in enumerate(res_ids):
            rid = int(rid)
            if rid in excluded:
                continue
            mask = per_res[k]
            if not mask.any():
                continue
            n_sole = int(np.count_nonzero(mask & only_one))
            n_shared = int(np.count_nonzero(mask) - n_sole)
            if n_sole:
                sole[rid] = sole.get(rid, 0.0) + n_sole * point_area
                n_pts[rid] = n_pts.get(rid, 0) + n_sole
            if n_shared:
                shared_pts = mask & ~only_one
                shared[rid] = shared.get(rid, 0.0) + n_shared * point_area
                # 1/k of each jointly occluded point, over the k PARTNER residues
                # occluding it, excluding points any ligand atom covers -- summed
                # over partners this is what makes buried_area + shared_area
                # coherent with buried_area's own reference state.
                credit_pts = shared_pts & ~ligand_occluded & (n_occ_partner >= 2)
                if credit_pts.any():
                    weighted[rid] = (weighted.get(rid, 0.0) + point_area *
                                     float((1.0 /
                                            n_occ_partner[credit_pts]).sum()))
                    # Exactly the points that contributed area, counted rather than
                    # area-weighted, so n_points > 0 iff contact_area > 0.
                    n_pts[rid] = (n_pts.get(rid, 0) +
                                  int(np.count_nonzero(credit_pts)))
                n_with_partner = int(np.count_nonzero(shared_pts &
                                                      (n_occ_partner >= 2)))
                if n_with_partner:
                    partner_shared[rid] = (partner_shared.get(rid, 0.0) +
                                           n_with_partner * point_area)
                partners = co_occ.setdefault(rid, set())
                for k2, rid2 in enumerate(res_ids):
                    if int(rid2) != rid and (per_res[k2] & shared_pts).any():
                        partners.add(int(rid2))

    out = {rid: ResidueContact(sole.get(rid, 0.0), shared.get(rid, 0.0),
                               len(co_occ.get(rid, ())), n_pairs[rid],
                               min_dist[rid], partner_shared.get(rid, 0.0),
                               weighted.get(rid, 0.0), n_pts.get(rid, 0))
           for rid in n_pairs}
    return (out, cg_free) if return_free_area else out
