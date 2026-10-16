"""Independent analytic checks of the resistive Maxwell prototype.

Passing these tests verifies selected numerical operators, not a resolved
pulsar magnetosphere, pair cascade, or a Vela emission calculation.
"""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from pulsar.plasma import (
    Magnetosphere, PlasmaConfig, curl, damp_parallel, divergence,
    drift_velocity, energy, filter4, invariants,
)


def vacuum_wave(n=32, mode=2, courant=.25):
    config = PlasmaConfig(n=n, plasma=False, conductor=False,
                          sponge_strength=0., filter_strength=0., courant=courant)
    model = Magnetosphere(config)
    k = 2*np.pi*mode/(2*config.half_extent)
    omega = np.sin(k*model.h)/model.h
    phase = k*model.xyz[0]
    # B_z = D_x A_y = cos(kx); E_y = -dt A_y = cos(kx).
    # A periodic central derivative has symbol i sin(kh)/h, not ik.
    model.a.fill(0.)
    model.e.fill(0.)
    model.a[1] = np.sin(phase)/omega
    model.e[1] = np.cos(phase)
    model.initial_energy = energy(model.e, model.b, model.h)
    return model, k, omega


class DifferentialOperatorTests(unittest.TestCase):
    def test_divergence_of_discrete_curl_is_zero(self):
        rng = np.random.default_rng(1049)
        potential = rng.normal(size=(3, 16, 16, 16))
        magnetic = curl(potential, .2)
        residual = divergence(magnetic, .2)
        self.assertLess(np.max(np.abs(residual)), 2e-13)

    def test_periodic_derivative_converges_quadratically_to_continuum(self):
        errors = []
        for n in (16, 32, 64):
            h = 8/n
            x = (np.arange(n)+.5)*h-4
            potential = np.zeros((3, n, 2, 2))
            k = 2*np.pi/8
            potential[1] = np.sin(k*x)[:, None, None]/k
            actual = curl(potential, h)[2, :, 0, 0]
            expected = np.cos(k*x)
            errors.append(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
        self.assertGreater(errors[0]/errors[1], 3.9)
        self.assertGreater(errors[1]/errors[2], 3.9)

    def test_filter_damps_fourier_modes_without_creating_divergence(self):
        n, strength = 16, .08
        i = np.arange(n)
        mode = np.cos(2*np.pi*(2*i[:, None, None]+3*i[None, :, None]+i[None, None, :])/n)
        potential = np.stack((mode, 2*mode, -mode))
        expected_factor = np.prod(1-strength*np.sin(np.pi*np.array([2,3,1])/n)**4)
        filtered = filter4(potential, strength)
        assert_allclose(filtered, potential*expected_factor, atol=8e-15)
        self.assertLess(np.max(np.abs(divergence(curl(filtered, .25), .25))), 2e-14)


class ClosureTests(unittest.TestCase):
    def test_parallel_decay_and_joule_energy_have_exact_analytic_limits(self):
        e = np.array([[2., -4., 0.], [3., 2., 0.], [5., -6., 0.]])
        b = np.array([[0., 0., 0.], [0., 0., 0.], [1., -2., 0.]])
        rate = np.array([0., 7., 3.])
        dt = .3
        actual = damp_parallel(e, b, rate, dt)
        expected = e.copy()
        expected[2, :2] *= np.exp(-rate[:2]*dt)
        assert_allclose(actual, expected, atol=9e-16)
        loss = .5*np.sum(e*e-actual*actual)
        analytic_loss = .5*np.sum(e[2,:2]**2*(-np.expm1(-2*rate[:2]*dt)))
        self.assertAlmostEqual(loss, analytic_loss, delta=7e-15)
        assert_allclose(damp_parallel(damp_parallel(e,b,rate,dt/2),b,rate,dt/2), actual, atol=9e-16)

    def test_uniform_parallel_field_in_stiff_limit_does_not_reverse(self):
        e = np.array([[1.], [2.], [3.]])
        b = np.array([[0.], [0.], [1.]])
        assert_allclose(damp_parallel(e,b,10000.,.1), [[1.],[2.],[0.]], atol=0.)

    def test_invariant_root_and_drift_remain_physical(self):
        # Magnetic dominance, electric dominance, null, parallel, B=0,
        # and fully zero fields all have independent simple limits.
        e = np.array([[0.,0.,0.,0.,2.,0.], [1.,3.,2.,0.,0.,0.], [0.,0.,0.,4.,0.,0.]])
        b = np.array([[0.,0.,0.,0.,0.,0.], [0.,0.,0.,0.,0.,0.], [2.,1.,2.,3.,0.,0.]])
        actual = drift_velocity(e,b)
        expected = np.zeros_like(e)
        expected[0,:3] = [.5, 1/3, 1.]
        assert_allclose(actual, expected, atol=2e-16)
        _,_,_,e0sq = invariants(e,b)
        assert_allclose(e0sq, [0.,8.,0.,16.,4.,0.], atol=0.)
        rng = np.random.default_rng(911)
        random_e = rng.normal(size=(3,10000))*10**rng.uniform(-6,6,size=10000)
        random_b = rng.normal(size=(3,10000))*10**rng.uniform(-6,6,size=10000)
        speed = np.linalg.norm(drift_velocity(random_e,random_b),axis=0)
        self.assertTrue(np.all(np.isfinite(speed)))
        self.assertLessEqual(float(speed.max()), 1.+4e-15)


class MaxwellEvolutionTests(unittest.TestCase):
    def test_vacuum_mode_reports_no_plasma_current_or_joule_power(self):
        model, _, _ = vacuum_wave(n=16)
        model.e.fill(0.)
        model.e[2] = model.b[2]
        result = model.diagnostics()
        self.assertGreater(result['parallel_E_fraction'], .99)
        self.assertEqual(result['current_norm_active_domain'], 0.)
        self.assertEqual(result['ohmic_power_active_domain'], 0.)

    def test_plane_wave_matches_discrete_dispersion_and_energy(self):
        model, k, omega = vacuum_wave()
        initial = model.initial_energy
        for _ in range(40):
            model.step()
        expected = np.cos(k*model.xyz[0]-omega*model.time)
        relative_error = np.linalg.norm(model.e[1]-expected)/np.linalg.norm(expected)
        self.assertLess(relative_error, 3e-6)
        assert_allclose(model.b[2],model.e[1],atol=3e-15)
        self.assertLess(abs(energy(model.e,model.b,model.h)/initial-1), 5e-7)
        # The continuum frequency would make a noticeably wrong prediction
        # here: this checks the central derivative's actual dispersion law.
        continuum = np.cos(k*model.xyz[0]-k*model.time)
        self.assertGreater(np.linalg.norm(model.e[1]-continuum)/np.linalg.norm(continuum), .05)
        self.assertEqual(model.ledger['ohmic_loss'],0.)
        self.assertEqual(model.ledger['conductor_exchange'],0.)
        self.assertEqual(model.ledger['sponge_exchange'],0.)
        self.assertEqual(model.ledger['filter_loss'],0.)

    def test_energy_ledger_accounts_for_all_applied_interventions(self):
        model = Magnetosphere(PlasmaConfig(n=16))
        for _ in range(5):
            model.step()
        result = model.diagnostics()
        self.assertLess(abs(result['energy_ledger_closure']), 3e-15)
        self.assertGreaterEqual(result['ledger']['ohmic_loss'], -1e-12)
        self.assertGreaterEqual(result['ledger']['filter_loss'], -1e-12)
        self.assertGreaterEqual(result['ledger']['inner_projection_loss'], -1e-12)
        self.assertLess(result['divB_relative'], 1e-15)
        self.assertLessEqual(result['max_drift_speed'], 1.+1e-14)
        # Accounting closure alone is not a physical accuracy assertion;
        # unlike this test, the wave test independently constrains evolution.

    def test_zero_vacuum_state_has_finite_diagnostics(self):
        model = Magnetosphere(PlasmaConfig(n=16,plasma=False,conductor=False,
                                          sponge_strength=0.,filter_strength=0.))
        model.a.fill(0.)
        model.e.fill(0.)
        model.initial_energy = 0.
        model.step()
        result = model.diagnostics()
        self.assertEqual(result['energy'],0.)
        self.assertEqual(result['divB_relative'],0.)
        self.assertEqual(result['max_drift_speed'],0.)


if __name__ == '__main__':
    unittest.main()
