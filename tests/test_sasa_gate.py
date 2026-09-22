import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ligand_vdgs.functions import sasa

def accessible_area(coords, radii, cg_atom_idxs, alive=None,
                    n_points=sasa.N_SPHERE_POINTS, probe_radius=sasa.PROBE_RADIUS):
    coords = np.asarray(coords, dtype=np.float64)
    radii = np.asarray(radii, dtype=np.float64)
    cg_atom_idxs = np.asarray(list(cg_atom_idxs), dtype=np.int64)
    if alive is None:
        alive = np.ones(len(coords), dtype=bool)
    unit = sasa._fibonacci_sphere(n_points)
    area = 0.0
    for cg_i in cg_atom_idxs:
        if not alive[cg_i]:
            continue
        sphere_r = radii[cg_i] + probe_radius
        points = coords[cg_i] + sphere_r * unit
        others = np.flatnonzero(alive)
        others = others[others != cg_i]
        free = np.ones(n_points, dtype=bool)
        for j in others:
            free &= np.linalg.norm(points - coords[j], axis=1) >= \
                (radii[j] + probe_radius)
        area += (4.0 * np.pi * sphere_r * sphere_r / n_points) * np.count_nonzero(free)
    return area

def buried_area_by_residue_bruteforce(coords, radii, resindices, cg_atom_idxs,
                                      exclude_resindices=(),
                                      n_points=sasa.N_SPHERE_POINTS,
                                      probe_radius=sasa.PROBE_RADIUS):
    coords = np.asarray(coords, dtype=np.float64)
    radii = np.asarray(radii, dtype=np.float64)
    resindices = np.asarray(resindices)
    excluded = {int(r) for r in exclude_resindices}

    def total_area(alive):
        return accessible_area(coords, radii, cg_atom_idxs, alive,
                               n_points, probe_radius)

    alive_all = np.ones(len(coords), dtype=bool)
    base = total_area(alive_all)
    out = {}
    for rid in np.unique(resindices):
        rid = int(rid)
        if rid in excluded:
            continue
        alive = alive_all.copy()
        alive[resindices == rid] = False
        out[rid] = total_area(alive) - base
    return out

def _loo_by_bruteforce(coords, radii, resindices, cg_idxs, exclude=(),
                       n_points=sasa.N_SPHERE_POINTS):
    return buried_area_by_residue_bruteforce(
        coords, radii, resindices, cg_idxs, exclude, n_points=n_points)

class PointSetGeometry(unittest.TestCase):
    def test_isolated_atom_area_matches_the_closed_form(self):
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.750])
        area = accessible_area(coords, radii, [0])
        self.assertAlmostEqual(area, 4 * np.pi * (1.750 + 1.4) ** 2, places=9)

    def test_two_sphere_overlap_matches_the_spherical_cap_formula(self):
        r1, r2, d = 1.750, 1.400, 3.4
        coords = np.array([[0.0, 0.0, 0.0], [d, 0.0, 0.0]])
        radii = np.array([r1, r2])
        big_r1, big_r2 = r1 + 1.4, r2 + 1.4
        h = big_r1 - (d * d + big_r1 ** 2 - big_r2 ** 2) / (2 * d)
        analytic_free = 4 * np.pi * big_r1 ** 2 - 2 * np.pi * big_r1 * h
        measured = accessible_area(coords, radii, [0], n_points=4096)
        self.assertLess(abs(measured - analytic_free) / analytic_free, 0.02,
                        msg='measured {:.3f} vs analytic {:.3f}'.format(
                            measured, analytic_free))

    def test_the_point_set_is_deterministic(self):
        a = sasa._fibonacci_sphere(512)
        b = sasa._fibonacci_sphere(512)
        self.assertTrue(np.array_equal(a, b))
        self.assertTrue(np.allclose(np.linalg.norm(a, axis=1), 1.0))

class _Scene(object):
    def __init__(self):
        self.xyz, self.elems, self.res = [], [], []

    def add(self, xyz, elem, resindex):
        self.xyz.append(xyz)
        self.elems.append(elem)
        self.res.append(resindex)
        return len(self.xyz) - 1

    def arrays(self):
        return (np.array(self.xyz, dtype=np.float64),
                sasa.radii_for(self.elems),
                np.array(self.res, dtype=np.int64))

class AttributionEqualsLeaveOneOut(unittest.TestCase):

    def _shared_patch_scene(self):
        s = _Scene()
        cg_a = s.add([0.0, 0.0, 0.0], 'C', 0)
        cg_b = s.add([1.5, 0.0, 0.0], 'O', 0)
        s.add([0.0, -3.4, 0.0], 'C', 0)
        s.add([0.0, 3.3, 0.0], 'O', 1)
        s.add([2.6, 3.2, 0.0], 'N', 2)
        s.add([0.0, 6.1, 0.0], 'C', 3)
        s.add([2.6, -3.5, 0.0], 'N', 4)
        s.add([0.0, 0.0, 5.9], 'C', 5)
        return s, [cg_a, cg_b]

    def test_buried_area_equals_two_run_leave_one_out(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        want = _loo_by_bruteforce(coords, radii, res, cg_idxs, exclude=(0,))
        for rid, expected in want.items():
            measured = got[rid].buried_area if rid in got else 0.0
            self.assertAlmostEqual(
                measured, expected, places=9,
                msg='residue {}: attribution {:.6f} vs leave-one-out {:.6f}'.format(
                    rid, measured, expected))

    def test_the_scene_actually_exercises_co_occlusion(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        self.assertGreater(got[1].shadowed_area, 0.0)
        self.assertGreater(got[2].shadowed_area, 0.0)
        self.assertIn(2, self._co_occluder_ids(coords, radii, res, cg_idxs, 1))
        self.assertGreater(got[4].shadowed_area, 0.0,
                           msg='residue 4 must share a patch with the ligand\'s '
                               'own non-CG atom')
        self.assertGreaterEqual(got[4].n_co_occluders, 1,
                                msg='the ligand occludes, so it counts as a '
                                    'co-occluder even though it is never a partner')

    @staticmethod
    def _co_occluder_ids(coords, radii, res, cg_idxs, rid):
        alive = np.ones(len(coords), dtype=bool)
        base = accessible_area(coords, radii, cg_idxs, alive)
        alive_r = alive.copy()
        alive_r[res == rid] = False
        solo = accessible_area(coords, radii, cg_idxs, alive_r)
        out = []
        for other in np.unique(res):
            if other == rid:
                continue
            alive_o = alive.copy()
            alive_o[res == other] = False
            other_solo = accessible_area(coords, radii, cg_idxs, alive_o)
            alive_both = alive_r.copy()
            alive_both[res == other] = False
            both = accessible_area(coords, radii, cg_idxs, alive_both)
            if both - base > (solo - base) + (other_solo - base) + 1e-9:
                out.append(int(other))
        return out

    def test_a_residue_behind_another_buries_nothing(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        self.assertIn(3, got, 'the shadowed residue must still be a candidate')
        self.assertEqual(got[3].buried_area, 0.0)
        self.assertGreater(got[3].min_heavy_dist, 4.5,
                           'and it is the kind of residue a 4.5 A gate also drops')

    def test_partners_partition_the_surface_the_ligand_does_not_already_bury(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        ligand_only = res == 0
        buriable = (accessible_area(coords, radii, cg_idxs, ligand_only)
                    - accessible_area(coords, radii, cg_idxs))
        credited = sum(sasa.contact_area(c) for c in got.values())
        self.assertGreater(buriable, 0.0, 'the scene must have buriable surface')
        self.assertAlmostEqual(credited, buriable, places=9)

    def test_n_points_is_nonzero_exactly_when_area_is(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        saw_shared_only = False
        for rid, c in got.items():
            if sasa.contact_area(c) > 0.0:
                self.assertGreaterEqual(c.n_points, 1,
                                        'residue {} has area but no points'.format(rid))
            else:
                self.assertEqual(c.n_points, 0,
                                 'residue {} has points but no area'.format(rid))
            if c.buried_area == 0.0 and c.shared_area > 0.0:
                saw_shared_only = True
                self.assertGreaterEqual(c.n_points, 1)
        self.assertTrue(saw_shared_only,
                        'scene must contain a fully shared patch or the test is '
                        'not discriminating')

    def test_n_points_union_equals_the_partner_occluded_surface(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        unit = sasa._fibonacci_sphere(sasa.N_SPHERE_POINTS)
        expected = 0
        for cg_i in cg_idxs:
            sphere_r = radii[cg_i] + sasa.PROBE_RADIUS
            pts = coords[cg_i] + sphere_r * unit
            by_partner = np.zeros(len(unit), dtype=bool)
            by_ligand = np.zeros(len(unit), dtype=bool)
            for j in range(len(coords)):
                if j == cg_i:
                    continue
                occ = (np.linalg.norm(pts - coords[j], axis=1)
                       < radii[j] + sasa.PROBE_RADIUS)
                if int(res[j]) == 0:
                    by_ligand |= occ
                else:
                    by_partner |= occ
            expected += int(np.count_nonzero(by_partner & ~by_ligand))
        total = sum(c.n_points for c in got.values())
        self.assertGreater(expected, 0, 'scene must bury something')
        self.assertGreaterEqual(
            total, expected,
            'the sum over partners must be at least the union; a smaller sum means '
            'points were dropped from the attribution')
        n_shared = sum(1 for c in got.values() if c.shared_area > 0.0)
        if n_shared:
            self.assertGreater(total, expected,
                               'a scene with shared patches must sum above its union')

    def test_n_points_does_not_scale_with_the_radius_table(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        base = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                           exclude_resindices={0})
        k = 1.35
        scaled = sasa.buried_area_by_residue(
            coords * k, radii * k, res, cg_idxs, exclude_resindices={0},
            probe_radius=sasa.PROBE_RADIUS * k)
        self.assertEqual({r: c.n_points for r, c in base.items()},
                         {r: c.n_points for r, c in scaled.items()},
                         'point counts must be invariant under a uniform rescale')
        grew = [sasa.contact_area(scaled[r]) / sasa.contact_area(base[r])
                for r in base
                if sasa.contact_area(base[r]) > 0]
        self.assertTrue(grew, 'scene must bury something')
        for ratio in grew:
            self.assertAlmostEqual(
                ratio, k * k, places=6,
                msg='areas must scale as (r + probe)^2 while counts do not -- that '
                    'difference is the whole argument for an n_points floor')

    def test_cg_free_sasa_is_the_ligand_alone_reference_state(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got, free = sasa.buried_area_by_residue(
            coords, radii, res, cg_idxs, exclude_resindices={0},
            return_free_area=True)
        ligand_only = res == 0
        self.assertAlmostEqual(
            free, accessible_area(coords, radii, cg_idxs, ligand_only),
            places=9)
        credited = sum(sasa.contact_area(c) for c in got.values())
        self.assertAlmostEqual(
            free - credited, accessible_area(coords, radii, cg_idxs),
            places=9)
        self.assertGreater(free, credited,
                           'this scene must leave some CG surface unburied')

    def test_cg_free_sasa_counts_a_wholly_exposed_cg_atom(self):
        coords = np.array([[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]])
        radii = np.array([1.75, 1.75])
        res = np.array([0, 1])
        _got, free = sasa.buried_area_by_residue(
            coords, radii, res, [0], exclude_resindices={0},
            return_free_area=True)
        sphere_r = 1.75 + sasa.PROBE_RADIUS
        self.assertAlmostEqual(free, 4.0 * np.pi * sphere_r ** 2, places=6)

    def test_area_the_ligand_already_buries_is_credited_to_no_partner(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        with_occluder = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                                   exclude_resindices={0})
        keep = np.array([i for i in range(len(coords))
                         if res[i] != 0 or i in cg_idxs])
        without = sasa.buried_area_by_residue(
            coords[keep], radii[keep], res[keep],
            [list(keep).index(i) for i in cg_idxs], exclude_resindices={0})
        self.assertGreater(sasa.contact_area(without[4]),
                           sasa.contact_area(with_occluder[4]),
                           'the ligand atom must be taking area away from residue 4, '
                           'not sharing it')

    def test_shared_area_rescues_a_fully_shadowed_contact(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        self.assertEqual(got[3].buried_area, 0.0)
        self.assertGreater(got[3].shared_area, 0.0,
                           'residue 3 is behind residues 1 and 2, so no theta on '
                           'exclusive area can admit it')

    def test_shadowed_area_is_not_the_credit(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        for rid, c in got.items():
            if c.shadowed_area > 0.0:
                self.assertLess(c.shared_area, c.shadowed_area,
                                'residue {}'.format(rid))

    def test_equality_holds_at_a_different_point_count(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        for n in (97, 512):
            got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                              exclude_resindices={0},
                                              n_points=n)
            want = _loo_by_bruteforce(coords, radii, res, cg_idxs, exclude=(0,),
                                      n_points=n)
            for rid, expected in want.items():
                measured = got[rid].buried_area if rid in got else 0.0
                self.assertAlmostEqual(measured, expected, places=9,
                                       msg='n_points={} residue={}'.format(n, rid))

    def test_excluded_residues_occlude_but_are_never_reported(self):
        s, cg_idxs = self._shared_patch_scene()
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, cg_idxs,
                                          exclude_resindices={0})
        self.assertNotIn(0, got)
        keep = np.array([i for i in range(len(coords)) if res[i] != 0 or i in cg_idxs])
        got2 = sasa.buried_area_by_residue(coords[keep], radii[keep], res[keep],
                                           [list(keep).index(i) for i in cg_idxs],
                                           exclude_resindices={0})
        self.assertGreater(got2[4].buried_area, got[4].buried_area)

class DiscriminatingContactCases(unittest.TestCase):
    def test_calpha_h_to_o_at_3_6_angstrom_is_kept(self):
        s = _Scene()
        cg = s.add([0.0, 0.0, 0.0], 'O', 0)
        s.add([3.6, 0.0, 0.0], 'C', 1)
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, [cg],
                                          exclude_resindices={0})
        self.assertGreater(got[1].buried_area, 1.0)
        self.assertAlmostEqual(got[1].min_heavy_dist, 3.6, places=6)
        self.assertEqual(got[1].n_atom_pairs, 1)

    def test_contact_through_a_non_cg_ligand_atom_only_is_dropped(self):
        s = _Scene()
        cg = s.add([0.0, 0.0, 0.0], 'C', 0)
        s.add([9.0, 0.0, 0.0], 'C', 0)
        s.add([13.4, 0.0, 0.0], 'N', 1)
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, [cg],
                                          exclude_resindices={0})
        self.assertNotIn(1, got, 'not even a candidate against the CG atom')

    def test_met_sd_5_5_angstrom_off_a_ring_is_kept(self):
        s = _Scene()
        ring = [s.add([1.4 * np.cos(t), 1.4 * np.sin(t), 0.0], 'C', 0)
                for t in np.linspace(0, 2 * np.pi, 6, endpoint=False)]
        s.add([0.0, 0.0, 5.5], 'S', 1)
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, ring,
                                          exclude_resindices={0})
        self.assertIn(1, got)
        self.assertGreater(got[1].min_heavy_dist, 4.5,
                           'the case only matters if a 4.5 A gate would drop it')
        self.assertGreater(sasa.contact_area(got[1]), sasa.MIN_CONTACT_AREA,
                           'Met SD at 5.5 A off a ring face must clear theta, not '
                           'merely bury something')

class PrefilterReach(unittest.TestCase):
    def _pair(self, d, elem='C'):
        s = _Scene()
        cg = s.add([0.0, 0.0, 0.0], 'C', 0)
        s.add([d, 0.0, 0.0], elem, 1)
        coords, radii, res = s.arrays()
        return sasa.buried_area_by_residue(coords, radii, res, [cg],
                                           exclude_resindices={0})

    def test_beyond_the_reach_nothing_is_buried(self):
        reach = 1.750 + 1.750 + 2 * 1.4
        got = self._pair(reach + 0.01)
        self.assertEqual(got, {}, 'an atom past the reach is not even a candidate')

    def test_just_inside_the_reach_area_is_buried(self):
        reach = 1.750 + 1.750 + 2 * 1.4
        got = self._pair(reach - 0.35, 'C')
        self.assertIn(1, got)
        self.assertGreater(got[1].buried_area, 0.0)

    def test_iodine_needs_the_wider_reach(self):
        d = 2.100 + 2.100 + 2 * 1.4 - 0.35
        s = _Scene()
        cg = s.add([0.0, 0.0, 0.0], 'I', 0)
        s.add([d, 0.0, 0.0], 'I', 1)
        coords, radii, res = s.arrays()
        got = sasa.buried_area_by_residue(coords, radii, res, [cg],
                                          exclude_resindices={0})
        self.assertGreater(got[1].buried_area, 0.0)
        self.assertGreater(d, 6.5, 'this case only bites beyond the nominal 6.5 A')

class ElementInference(unittest.TestCase):

    def test_c_alpha_is_carbon_not_calcium(self):
        self.assertEqual(sasa.element_of('CA', '', 'ALA'), 'C')
        self.assertEqual(sasa.element_of('CD1', '', 'LEU'), 'C')
        self.assertEqual(sasa.element_of('CA', '', 'CA'), 'CA')
        self.assertEqual(sasa.radius_of(sasa.element_of('CA', '', 'ALA')), 1.750)

    def test_the_element_column_wins_when_present(self):
        self.assertEqual(sasa.element_of('CA', 'C', 'CA'), 'C')

    def test_raw_four_character_names_follow_the_column_convention(self):
        self.assertEqual(sasa.element_of(' CA ', ''), 'C')
        self.assertEqual(sasa.element_of('CA  ', ''), 'CA')

    def test_selenomethionine_keeps_its_selenium(self):
        self.assertEqual(sasa.element_of('SE', '', 'MSE'), 'SE')
        self.assertEqual(sasa.radius_of('SE'), 1.90)

    def test_an_unknown_element_gets_the_documented_default(self):
        self.assertEqual(sasa.radius_of('XX'), sasa.DEFAULT_RADIUS)

if __name__ == '__main__':
    unittest.main()
