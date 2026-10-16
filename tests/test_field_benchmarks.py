"""Analytic tests deliberately independent of the solver's energy ledger."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from pulsar.field_benchmarks import (
    ExactDipoleBoundary, analytic_maxwell_residual, analytic_spherical_flux,
    deutsch_fields, deutsch_luminosity, deutsch_spherical_flux, finite_conductor_reference,
    gauge_sponge_check, point_dipole_luminosity, retarded_dipole,
    wave_packet_absorber,
)


class AnalyticRetardedDipoleTests(unittest.TestCase):
    def test_aligned_limit_is_static_magnetic_dipole_without_electric_field(self):
        xyz = np.array([[.4, 1.2], [.3, -.5], [1.1, .6]])
        r = np.linalg.norm(xyz, axis=0)
        n = xyz/r
        m = np.array([0., 0., 1.])[:, None]
        _, e, b = retarded_dipole(xyz, 1.23, inclination_deg=0.)
        assert_allclose(e, 0., atol=0.)
        assert_allclose(b, (3*n*np.sum(n*m, axis=0)-m)/r**3, rtol=2e-15)

    def test_temporal_gauge_potential_matches_electric_field(self):
        xyz = np.array([[.4, 1.2], [.3, -.5], [1.1, .6]])
        eps = 1e-5
        am, _, _ = retarded_dipole(xyz, .73-eps)
        ap, _, _ = retarded_dipole(xyz, .73+eps)
        _, e, _ = retarded_dipole(xyz, .73)
        assert_allclose(-(ap-am)/(2*eps), e, rtol=2e-9, atol=2e-10)

    def test_poynting_flux_matches_closed_form_at_all_radii_and_phases(self):
        for inclination in (0., 30., 65., 90.):
            for radius, phase in ((.5, 0.), (1., .73), (2.5, 2.1), (10., 5.)):
                self.assertAlmostEqual(analytic_spherical_flux(radius, phase, inclination),
                                       point_dipole_luminosity(inclination), places=12)

    def test_continuum_maxwell_residuals_converge_with_independent_fields(self):
        lo, hi = analytic_maxwell_residual(32), analytic_maxwell_residual(64)
        for key in ('faraday_relative_L2', 'ampere_relative_L2', 'curlA_relative_L2',
                    'divB_relative_L2', 'divE_relative_L2'):
            self.assertLess(hi[key], lo[key]/3.)
        self.assertLess(hi['faraday_relative_L2'], .02)

    def test_origin_rejected_in_exact_solution(self):
        with self.assertRaises(ValueError):
            retarded_dipole(np.zeros((3, 1)), 0.)


class BoundaryExperimentsTests(unittest.TestCase):
    def test_finite_conductor_matches_surface_conditions_and_published_luminosity(self):
        for R in (.1, .35, .6):
            ref = finite_conductor_reference(star_radius=R)
            self.assertLess(ref['radial_B_surface_relative_error'], 2e-15)
            self.assertLess(ref['tangential_E_surface_relative_error'], 3e-15)
            for row in ref['fluxes']:
                self.assertAlmostEqual(row['flux'], ref['finite_luminosity'], places=11)

    def test_finite_conductor_satisfies_independent_local_maxwell_derivatives(self):
        xyz = np.array([[.4, 1.2, -.7], [.3, -.5, 1.1], [1.1, .6, -.4]])
        eps, t = 1e-5, .73
        dE, dB = [], []
        for axis in range(3):
            delta = np.zeros_like(xyz)
            delta[axis] = eps
            ep, bp = deutsch_fields(xyz+delta, t)
            em, bm = deutsch_fields(xyz-delta, t)
            dE.append((ep-em)/(2*eps))
            dB.append((bp-bm)/(2*eps))
        local_curl = lambda d: np.array([d[1][2]-d[2][1], d[2][0]-d[0][2], d[0][1]-d[1][0]])
        ep, bp = deutsch_fields(xyz, t+eps)
        em, bm = deutsch_fields(xyz, t-eps)
        assert_allclose(local_curl(dE), -(bp-bm)/(2*eps), rtol=1e-7, atol=1e-8)
        assert_allclose(local_curl(dB), (ep-em)/(2*eps), rtol=1e-7, atol=1e-8)
        assert_allclose(sum(dE[j][j] for j in range(3)), 0., atol=1e-8)
        assert_allclose(sum(dB[j][j] for j in range(3)), 0., atol=1e-8)

    def test_pure_gauge_potential_sponge_creates_unphysical_magnetic_energy(self):
        ordinary = gauge_sponge_check(gauge_offset=1.)
        doubled = gauge_sponge_check(gauge_offset=2.)
        self.assertGreater(ordinary['added_magnetic_energy_per_transverse_area'], .01)
        self.assertEqual(doubled['added_magnetic_energy_per_transverse_area'],
                         4*ordinary['added_magnetic_energy_per_transverse_area'])
        self.assertEqual(ordinary['matched_fields_added_energy'], 0.)

    def test_outgoing_packet_exposes_potential_sponge_reflection(self):
        old = wave_packet_absorber(n=64)
        fixed = wave_packet_absorber(n=64, kind='matched_fields')
        self.assertGreater(old['peak_reflected_energy_fraction'], 1e-3)
        self.assertLess(fixed['peak_reflected_energy_fraction'], 1e-24)
        # Stronger matched damping should reduce periodic wrap transmission.
        strong = wave_packet_absorber(n=64, kind='matched_fields', strength=12.)
        self.assertLess(strong['peak_periodic_wrap_energy_fraction'],
                        fixed['peak_periodic_wrap_energy_fraction']/100)

    def test_exact_driver_sets_correct_stage_time_and_preserves_curl_constraint(self):
        model = ExactDipoleBoundary(n=32)
        model.step(.02)
        a, e, _ = model.exact(.02)
        assert_allclose(model.a[:, model.guard], a[:, model.guard], atol=0.)
        assert_allclose(model.e[:, model.guard], e[:, model.guard], atol=0.)
        self.assertEqual(model.time, .02)
        self.assertLess(model.benchmark_diagnostics()['divB_relative'], 1e-15)


if __name__ == '__main__':
    unittest.main()
