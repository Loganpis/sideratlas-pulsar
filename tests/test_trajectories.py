"""Independent geometric checks for the reduced slot-gap model.

These verify its mathematics, not the physical adequacy of vacuum geometry
for Vela. In particular the curvature reference integrates a trajectory in
the explicitly time-dependent laboratory field, independently of the
production rotating-coordinate curvature calculation.
"""
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp

from pulsar.trajectories import (
    Geometry, curvature, dipole_field, polar_cap_rim, rotate_z,
    trace_bundle, velocity, sky_profiles,
)


class FieldAndTrajectoryTests(unittest.TestCase):
    def setUp(self):
        self.chi = np.deg2rad(75.)
        self.points = np.array([[.11, .03, .12], [.5, -.25, .3], [.7, .4, .1]])

    def test_retarded_field_rotation_covariance_and_aligned_limit(self):
        angle, time = .73, -.21
        assert_allclose(
            dipole_field(rotate_z(self.points, angle), self.chi, time + angle),
            rotate_z(dipole_field(self.points, self.chi, time), angle),
            rtol=2e-14, atol=2e-13,
        )
        r = np.linalg.norm(self.points, axis=-1)
        n = self.points / r[:, None]
        static = (3 * n * n[:, 2, None] - np.array([0., 0., 1.])) / r[:, None]**3
        assert_allclose(dipole_field(self.points, 0., time), static, rtol=2e-15)

    def test_velocity_unit_speed_ideal_force_and_work_projection(self):
        q = np.concatenate((self.points, -self.points))
        pole = np.r_[np.ones(3), -np.ones(3)]
        v, relative, projection = velocity(q, self.chi, pole)
        b = dipole_field(q, self.chi)
        beta = np.column_stack((-q[:, 1], q[:, 0], np.zeros(6)))
        electric = -np.cross(beta, b)
        assert_allclose(np.linalg.norm(v, axis=-1), 1., atol=5e-15)
        assert_allclose(electric + np.cross(v, b), 0., atol=2e-13)
        assert_allclose(relative + beta, v, atol=1e-15)
        bhat = b / np.linalg.norm(b, axis=-1)[:, None]
        assert_allclose(np.sum(v * bhat, axis=-1), pole * projection, atol=5e-15)
        assert_allclose(v[:3], -v[3:], atol=2e-15)

    def test_curvature_against_integrated_inertial_trajectory(self):
        # Separate laboratory-frame ODE: B depends explicitly on emission time.
        # Merely differentiating the frozen magnetic line would fail this test.
        def laboratory_rhs(time, position):
            b = dipole_field(position, self.chi, time)
            b = b / np.linalg.norm(b)
            beta = np.array([-position[1], position[0], 0.])
            dot = np.dot(beta, b)
            f = -dot + np.sqrt(dot**2 + 1 - np.dot(beta, beta))
            return beta + f * b

        reference = 1 / curvature(self.points, self.chi)
        errors = []
        for h in (.004, .002, .001):
            accelerations = []
            for q in self.points:
                plus = solve_ivp(laboratory_rhs, (0., h), q,
                                 rtol=1e-12, atol=1e-14).y[:, -1]
                minus = solve_ivp(laboratory_rhs, (0., -h), q,
                                  rtol=1e-12, atol=1e-14).y[:, -1]
                acceleration = (laboratory_rhs(h, plus) - laboratory_rhs(-h, minus)) / (2*h)
                accelerations.append(np.linalg.norm(acceleration))
            errors.append(np.max(np.abs(np.asarray(accelerations) / reference - 1)))
        self.assertLess(errors[-1], 5e-7)
        self.assertGreater(errors[0] / errors[1], 3.5)
        self.assertGreater(errors[1] / errors[2], 3.5)

    def test_arrival_phase_from_explicit_observer_geometry(self):
        v = velocity(self.points, self.chi)[0]
        emission_time = -np.arctan2(v[:, 1], v[:, 0])
        emission_position = rotate_z(self.points, emission_time)
        observer = np.column_stack((np.sqrt(1 - v[:, 2]**2), np.zeros(3), v[:, 2]))
        assert_allclose(rotate_z(v, emission_time), observer, atol=2e-15)
        arrival_time = emission_time - np.sum(observer * emission_position, axis=-1)
        expected_phase = arrival_time / (2*np.pi) % 1
        actual_phase = (-np.arctan2(v[:, 1], v[:, 0])
                        - np.sum(v * self.points, axis=-1)) / (2*np.pi) % 1
        assert_allclose(actual_phase, expected_phase, atol=2e-15)

    def test_aligned_cap_matches_analytic_last_open_dipole(self):
        g = Geometry(inclination_deg=0.)
        az = np.array([.3, 2., 4.])
        exact = np.arcsin(np.sqrt(g.seed_radius))
        # r = r_equator sin^2(theta); the last open line has r_equator=1.
        rim = polar_cap_rim(g, az, trace_step=.005)
        assert_allclose(rim, exact, atol=2e-6, rtol=0.)

    def test_oblique_cap_step_convergence(self):
        g = Geometry()
        az = np.array([.3, 2., 4.])
        coarse = polar_cap_rim(g, az, trace_step=.012)
        fine = polar_cap_rim(g, az, trace_step=.006)
        assert_allclose(coarse, fine, atol=1e-5, rtol=0.)


class BundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g = Geometry()
        cls.bundle = trace_bundle(cls.g, n_azimuth=16, n_rings=1, dt=.006)

    def test_finite_samples_masks_pole_parity_and_seed_flux(self):
        b = self.bundle
        active = b['dt'] > 0
        self.assertTrue(np.all(np.diff(active.astype(int), axis=0) <= 0))
        self.assertTrue(np.all(np.isfinite(b['position'])))
        self.assertTrue(np.all(np.isfinite(b['rho'][active])))
        self.assertTrue(np.all(b['rho'][active] > 0))
        self.assertTrue(np.all(b['dt'] >= 0))
        self.assertTrue(np.all(b['dt'] <= b['step'] * (1 + 1e-12)))
        assert_allclose(np.linalg.norm(b['direction'][active], axis=-1), 1., atol=1e-14)
        assert_allclose(b['direction'][~active], 0., atol=0.)
        half = b['seed_weight'].size // 2
        assert_allclose(b['seed_weight'][:half], b['seed_weight'][half:], atol=1e-15)
        self.assertAlmostEqual(float(np.sum(b['seed_weight'])), 1., places=14)
        self.assertTrue(np.all(b['seed_weight'] > 0))
        assert_allclose(b['position'][:, :half], -b['position'][:, half:], atol=2e-12)
        assert_allclose(b['direction'][:, :half], -b['direction'][:, half:], atol=2e-12)
        assert_allclose(b['rho'][:, :half], b['rho'][:, half:], rtol=2e-10)

    def test_emission_respects_injection_and_escape_boundaries(self):
        b = self.bundle
        position = b['position'][b['dt'] > 0]
        radius = np.linalg.norm(position, axis=-1)
        cylindrical = np.linalg.norm(position[:, :2], axis=-1)
        self.assertGreaterEqual(float(np.min(radius)), self.g.seed_radius * .99 - 1e-10)
        self.assertLessEqual(float(np.max(radius)), self.g.radial_cutoff + 1e-10)
        self.assertLessEqual(float(np.max(cylindrical)), self.g.cylindrical_cutoff + 1e-10)

    def test_aligned_seed_weights_match_integrated_magnetic_flux(self):
        g = Geometry(inclination_deg=0.)
        nphi, nring = 8, 3
        b = trace_bundle(g, n_azimuth=nphi, n_rings=nring, dt=.008)
        rim = b['cap_rim'][0]
        edges = rim * np.linspace(g.annulus_inner, g.annulus_outer, nring + 1)
        # Exact axisymmetric flux through each ring is proportional to
        # sin²(theta_outer)-sin²(theta_inner), independent of field-line tracing.
        ring_flux = np.diff(np.sin(edges)**2)
        expected = np.repeat(ring_flux / ring_flux.sum() / (2*nphi), nphi)
        assert_allclose(b['seed_weight'][:nphi*nring], expected, rtol=5e-6, atol=0.)


class UnsmearedProfileTests(unittest.TestCase):
    def test_zero_phase_blur_preserves_single_bin_and_periodic_deposition(self):
        g = Geometry(observer_deg=90.,phase_sigma=0.)
        bundle = dict(direction=np.array([[[1.,0.,0.]]]),
                      seed_weight=np.ones(1),dt=np.ones((1,1)),
                      phase=np.array([[.25]]))
        _,profile = sky_profiles(bundle,np.ones((1,1,1)),g,nphase=16)
        expected = np.zeros((1,16)); expected[0,4] = 16
        assert_allclose(profile,expected,atol=0)
        # A photon halfway across the final-bin boundary splits equally
        # between the final bin and phase zero, without a Gaussian halo.
        bundle['phase'][0,0] = 15.5/16
        _,profile = sky_profiles(bundle,np.ones((1,1,1)),g,nphase=16)
        expected[:] = 0; expected[0,[0,15]] = 8
        assert_allclose(profile,expected,atol=0)


if __name__ == '__main__':
    unittest.main()
