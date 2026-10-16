"""Independent physical checks for laboratory-frame curvature transport."""
import unittest

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import kv

from pulsar.radiation import (
    C_CM_S, E_STATC, GEV_ERG, HBAR_ERG_S, MEC2_ERG,
    KERNEL_POWER_INTEGRAL, KERNEL_PHOTON_INTEGRAL,
    critical_energy_gev, curvature_kernel, curvature_photon_yield,
    curvature_band_energy_fraction, curvature_power_erg_s, transport_curvature,
)


class CurvatureRadiationTests(unittest.TestCase):
    def test_kernel_against_independent_bessel_quadrature(self):
        for x in [1e-3, .03, .3, 1., 4., 15., 40.]:
            reference = x*quad(lambda t: kv(5/3, t), x, np.inf,
                               epsabs=1e-28, epsrel=1e-10)[0]
            self.assertAlmostEqual(float(curvature_kernel(x))/reference, 1., delta=2e-5)
        self.assertEqual(float(curvature_kernel(0)), 0)
        self.assertEqual(float(curvature_kernel(np.inf)), 0)

    def test_total_power_and_photon_kernel_normalizations(self):
        power = quad(lambda z: float(curvature_kernel(np.exp(z)))*np.exp(z),
                     -35, np.log(80), epsabs=2e-7, limit=250)[0]
        photons = quad(lambda z: float(curvature_kernel(np.exp(z))),
                       -90, np.log(80), epsabs=1e-6, limit=250)[0]
        self.assertAlmostEqual(power/KERNEL_POWER_INTEGRAL, 1, delta=3e-6)
        self.assertAlmostEqual(photons/KERNEL_PHOTON_INTEGRAL, 1, delta=3e-6)

    def test_energy_units_and_photon_numbers(self):
        gamma, rho = 2e7, 3e8
        ec = critical_energy_gev(gamma, rho)
        self.assertAlmostEqual(ec*GEV_ERG,
                               3*HBAR_ERG_S*C_CM_S*gamma**3/(2*rho), delta=1e-18)
        radiated, lo, hi = 2.0, .1, 3.
        direct = quad(lambda energy: radiated/(ec*GEV_ERG)*
                      float(curvature_kernel(energy/ec))/KERNEL_POWER_INTEGRAL/energy,
                      lo, hi, epsabs=1e-7, limit=100)[0]
        tabulated = curvature_photon_yield(radiated, ec, [[lo, hi]])[0]
        self.assertAlmostEqual(tabulated/direct, 1, delta=5e-6)
        full_photons = curvature_photon_yield(radiated, ec, [[0, np.inf]])[0]
        expected = radiated/(ec*GEV_ERG)*KERNEL_PHOTON_INTEGRAL/KERNEL_POWER_INTEGRAL
        self.assertAlmostEqual(full_photons/expected, 1, delta=1e-12)
        self.assertAlmostEqual(curvature_band_energy_fraction(ec, [[0, np.inf]])[0], 1)

    def test_pure_loss_analytic_and_grid_independence(self):
        gamma0, rho, length = 8e7, 2e8, 3e8
        b = 2*E_STATC**2/(3*MEC2_ERG*rho**2)
        for n in [1, 8, 64]:
            result = transport_curvature(np.full(n, length/n), rho, 0., gamma0)
            s = np.arange(n+1)*length/n
            exact = (gamma0**-3+3*b*s)**(-1/3)
            np.testing.assert_allclose(result.gamma, exact, rtol=5e-14)
            self.assertTrue(np.all(result.gamma > 1))
            self.assertAlmostEqual(result.radiated_erg.sum(),
                                   MEC2_ERG*(gamma0-exact[-1]), delta=5e-12)

    def test_equilibrium_stiff_flow_and_work_accounting(self):
        rho, field = 2e8, 500.
        equilibrium = (3*field*rho**2/(2*E_STATC))**.25
        # Strong cooling over one coarse step must remain positive and not
        # overshoot the radiation-reaction equilibrium from either side.
        initial = equilibrium*np.array([.02, 1., 30.])
        result = transport_curvature([2e10], rho, field, initial, [[.1, 1], [1, 100]], 3)
        np.testing.assert_allclose(result.gamma[-1], equilibrium, rtol=2e-13)
        np.testing.assert_allclose(result.gamma[:, 1], equilibrium, rtol=2e-13)
        self.assertTrue(np.all(result.photon_yield >= 0))
        np.testing.assert_allclose(result.energy_balance_error_erg, 0, atol=2e-12)
        expected = curvature_power_erg_s(equilibrium, rho)*2e10/C_CM_S
        self.assertAlmostEqual(result.radiated_erg[0, 1]/expected, 1, delta=5e-14)

    def test_mixed_flow_against_independent_ode_solver(self):
        rho, field, length = 4e8, 50., 8e7
        a, b = E_STATC*field/MEC2_ERG, 2*E_STATC**2/(3*MEC2_ERG*rho**2)
        for initial in [1e4, 2e7, 1e8]:
            reference = solve_ivp(lambda s, g: a-b*g**4, (0, length), [initial],
                                  rtol=2e-12, atol=1e-5, method="DOP853").y[0, -1]
            result = transport_curvature([length], rho, field, initial)
            self.assertAlmostEqual(result.gamma[-1]/reference, 1, delta=2e-11)

    def test_padding_straight_motion_and_broadcast(self):
        ds = np.array([[1e5, 1e5], [0., 1e5], [0., 0.]])
        rho = np.array([[1e8, np.inf], [np.inf, np.inf], [np.inf, np.inf]])
        result = transport_curvature(ds, rho, 5., [1e6, 2e6], [[.1, 1.]])
        self.assertEqual(result.gamma.shape, (4, 2))
        np.testing.assert_array_equal(result.gamma[1:, 0], result.gamma[1, 0])
        np.testing.assert_array_equal(result.radiated_erg[:, 1], 0)
        np.testing.assert_array_equal(result.photon_yield[ds == 0], 0)
        np.testing.assert_array_equal(result.photon_yield[:, 1], 0)
        expected = 2e6+E_STATC*5*2e5/MEC2_ERG
        self.assertAlmostEqual(result.gamma[-1, 1], expected)

    def test_spectral_segment_refinement(self):
        # Rapid pure cooling: test spectrum convergence separately from
        # the exact bolometric cooling solution. Comparison has no detector
        # response and refers to emitted photons only.
        totals = []
        for n in [64, 128, 512]:
            result = transport_curvature(np.full(n, 1e8/n), 3e8, 0., 6e7,
                                         [[.1, .3], [.3, 1], [1, 3], [3, 100]], 3)
            totals.append(result.photon_yield.sum(axis=0))
        coarse = np.max(np.abs(totals[0]/totals[2]-1))
        fine = np.max(np.abs(totals[1]/totals[2]-1))
        self.assertLess(fine, coarse)
        self.assertLess(fine, 5e-5)

    def test_invalid_input_rejected(self):
        with self.assertRaises(ValueError):
            transport_curvature([1], 1e8, -1, 1e6)
        with self.assertRaises(ValueError):
            transport_curvature([1], 0, 1, 1e6)
        with self.assertRaises(ValueError):
            curvature_kernel(-1)


if __name__ == "__main__":
    unittest.main()
