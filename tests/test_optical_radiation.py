"""Physical limits and conservation tests, not observational validation of Vela."""
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.special import gamma as gamma_function

from pulsar.optical_radiation import (
    B_QUANTUM_T, C_M_S, E_CHARGE_C, EPSILON_0_F_M, M_E_KG,
    bin_lienard_radiation, bin_received_energy, incoherent_sum,
    lienard_angular_power_received_w_sr, lienard_angular_power_source_w_sr,
    lienard_total_power_w, retarded_arrival_times_s,
    synchrotron_critical_frequency_hz, synchrotron_f, synchrotron_power_w,
    synchrotron_spectrum_w_hz, synchrotron_wavelength_band_fraction,
)


def sphere_integral(beta, acceleration, n_mu=192, n_phi=64):
    """Independent sphere quadrature, with the polar axis along beta for tests."""
    mu, weight = np.polynomial.legendre.leggauss(n_mu)
    phi = 2*np.pi*np.arange(n_phi)/n_phi
    sin_theta = np.sqrt(1-mu*mu)[:, None]
    direction = np.stack(np.broadcast_arrays(sin_theta*np.cos(phi),
                                             sin_theta*np.sin(phi), mu[:, None]), axis=-1)
    power = lienard_angular_power_source_w_sr(beta, acceleration, direction)
    return np.sum(power*weight[:, None])*2*np.pi/n_phi


class LienardWiechertTests(unittest.TestCase):
    def test_nonrelativistic_dipole_pattern_integrates_to_larmor(self):
        acceleration = np.array([3., -4., 2.])
        expected = E_CHARGE_C**2*np.dot(acceleration, acceleration)/(6*np.pi*EPSILON_0_F_M*C_M_S)
        self.assertAlmostEqual(sphere_integral([0, 0, 0], acceleration)/expected, 1., delta=2e-14)
        self.assertAlmostEqual(float(lienard_total_power_w([0, 0, 0], acceleration))/expected,
                               1., delta=2e-15)
        # No radiation in the acceleration direction, maximum perpendicular.
        direction = acceleration/np.linalg.norm(acceleration)
        self.assertAlmostEqual(float(lienard_angular_power_source_w_sr([0, 0, 0], acceleration, direction)),
                               0., delta=expected*1e-28)

    def test_relativistic_perpendicular_and_parallel_lienard_limits(self):
        gamma = 8.
        beta = np.array([0., 0., np.sqrt(1-gamma**-2)])
        factor = E_CHARGE_C**2/(6*np.pi*EPSILON_0_F_M*C_M_S)
        for acceleration, expected in (([2., 0, 0], factor*4*gamma**4),
                                       ([0, 0, 2.], factor*4*gamma**6)):
            self.assertAlmostEqual(float(lienard_total_power_w(beta, acceleration))/expected,
                                   1., delta=5e-14)
            integrated = sphere_integral(beta, acceleration, n_mu=384)
            self.assertAlmostEqual(integrated/expected, 1., delta=2e-10)

    def test_angular_convergence_against_lienard_total(self):
        beta, acceleration = [0, 0, .91], [1.4, -.8, .3]
        coarse = sphere_integral(beta, acceleration, n_mu=64)
        fine = sphere_integral(beta, acceleration, n_mu=128)
        expected = float(lienard_total_power_w(beta, acceleration))
        self.assertAlmostEqual(coarse/fine, 1., delta=2e-10)
        self.assertAlmostEqual(fine/expected, 1., delta=2e-12)

    def test_source_and_observer_time_power_differ_by_arrival_jacobian(self):
        beta, acceleration, direction = [.3, 0, .8], [2., 1., .3], [0, 0, 1.]
        source = float(lienard_angular_power_source_w_sr(beta, acceleration, direction))
        received = float(lienard_angular_power_received_w_sr(beta, acceleration, direction))
        self.assertAlmostEqual(received/source, 5., delta=2e-14)
        self.assertEqual(float(lienard_total_power_w(beta, [0, 0, 0])), 0.)
        self.assertEqual(float(lienard_angular_power_source_w_sr(beta, [0, 0, 0], direction)), 0.)
        self.assertEqual(float(lienard_angular_power_source_w_sr(beta, acceleration, direction,
                                                              E_CHARGE_C)), source)

    def test_invalid_kinematics_and_directions_are_rejected(self):
        for beta in ([0, 0, 1.], [0, 0, np.nan], [0, 2, 0]):
            with self.assertRaises(ValueError):
                lienard_total_power_w(beta, [1., 0, 0])
        for direction in ([0, 0, 2.], [0, 0, 0], [0, np.nan, 1.]):
            with self.assertRaises(ValueError):
                lienard_angular_power_source_w_sr([0, 0, .8], [1., 0, 0], direction)


class ArrivalTimeTests(unittest.TestCase):
    def test_uniform_motion_has_the_analytic_retarded_time(self):
        time = np.linspace(0., 2., 31)
        position = np.zeros((len(time), 3))
        position[:, 2] = .8*C_M_S*time+3*C_M_S
        assert_allclose(retarded_arrival_times_s(time, position, [0, 0, 1]), .2*time-3., atol=9e-16)

    def test_weighted_energy_and_arrival_bin_power_conserve_total(self):
        result = bin_received_energy([0., .1, .3, .9, 1.], [2., 3., 5., 7., 11.],
                                     [0., .2, .5, 1.], particle_number_weights=4.)
        assert_allclose(result.energy_j_sr, [20., 20., 72.], atol=0.)
        self.assertEqual(float(result.energy_j_sr.sum()), 112.)
        self.assertEqual(float(np.sum(result.power_w_sr*np.diff(result.bin_edges_s))), 112.)
        self.assertFalse(result.energy_j_sr.flags.writeable)
        with self.assertRaisesRegex(ValueError, 'outside'):
            bin_received_energy([0., 2.], [1., 1.], [0., 1.])
        # An event with zero energy cannot discard any radiation.
        zero_outside = bin_received_energy([0., 2.], [1., 0.], [0., 1.])
        self.assertEqual(float(zero_outside.energy_j_sr.sum()), 1.)

    def test_lw_deposition_uses_source_time_once(self):
        # Consistent circular motion viewed along the instantaneous tangent
        # at t=0.  Evaluate source quadrature independently from the binner.
        intervals = np.full(1000, 1e-4)
        time = (np.arange(len(intervals))+.5)*intervals[0]
        speed, omega = .8, 2.
        phase = omega*time
        radius = speed*C_M_S/omega
        positions = radius*np.column_stack((np.cos(phase), np.sin(phase), np.zeros_like(phase)))
        beta = speed*np.column_stack((-np.sin(phase), np.cos(phase), np.zeros_like(phase)))
        acceleration = speed*omega*np.column_stack((-np.cos(phase), -np.sin(phase), np.zeros_like(phase)))
        direction = [0, 1., 0]
        arrival = retarded_arrival_times_s(time, positions, direction)
        edges = np.linspace(arrival.min()-1e-10, arrival.max()+1e-10, 51)
        result = bin_lienard_radiation(time, positions, beta, acceleration, direction,
                                       intervals, edges, particle_number_weights=17.)
        source_energy = np.sum(lienard_angular_power_source_w_sr(beta, acceleration, direction)*intervals)*17
        observer_energy = np.sum(lienard_angular_power_received_w_sr(beta, acceleration, direction)
                                 * (1-beta[:, 1])*intervals)*17
        self.assertAlmostEqual(float(result.energy_j_sr.sum())/source_energy, 1., delta=2e-15)
        self.assertAlmostEqual(observer_energy/source_energy, 1., delta=2e-15)
        self.assertAlmostEqual(float(np.sum(result.power_w_sr*np.diff(edges)))/source_energy,
                               1., delta=2e-15)

    def test_macroparticle_split_is_incoherent_and_linear(self):
        signal = np.array([[2., 4., 8.]])
        reference = incoherent_sum(signal, [90.])
        split = incoherent_sum(np.repeat(signal, 9, axis=0), np.full(9, 10.))
        assert_allclose(reference, split, atol=0.)
        assert_allclose(incoherent_sum(signal.T, [90.], axis=1), reference, atol=0.)
        self.assertEqual(float(bin_received_energy([.5]*9, [2.]*9, [0., 1.], [10.]*9)
                               .energy_j_sr.sum()), 180.)
        with self.assertRaises(ValueError):
            incoherent_sum(signal, [-1.])


class SynchrotronOpticalTests(unittest.TestCase):
    def test_fixed_pitch_power_matches_thomson_expression(self):
        gamma, field, pitch = 300., .004, .7
        sigma_t = 8*np.pi/3*(E_CHARGE_C**2/(4*np.pi*EPSILON_0_F_M*M_E_KG*C_M_S**2))**2
        # U_B = epsilon_0*c^2*B^2/2.  Fixed pitch has coefficient 2,
        # rather than the isotropic-pitch-averaged coefficient 4/3.
        energy_density = EPSILON_0_F_M*C_M_S**2*field**2/2
        expected = 2*sigma_t*C_M_S*energy_density*(gamma*gamma-1)*np.sin(pitch)**2
        self.assertAlmostEqual(float(synchrotron_power_w(gamma, field, pitch))/expected,
                               1., delta=2e-15)
        critical = 1.5*gamma*gamma*E_CHARGE_C*field*np.sin(pitch)/(2*np.pi*M_E_KG)
        self.assertAlmostEqual(float(synchrotron_critical_frequency_hz(gamma, field, pitch))/critical,
                               1., delta=2e-15)

    def test_spectrum_frequency_integral_recovers_fixed_pitch_larmor_power(self):
        gamma, field, pitch = 400., .04, .7
        critical = float(synchrotron_critical_frequency_hz(gamma, field, pitch))
        # Independent log-frequency quadrature covers essentially all power.
        x, weight = np.polynomial.legendre.leggauss(160)
        lower, upper = np.log(1e-12), np.log(200.)
        frequency = critical*np.exp((upper+lower)/2+(upper-lower)*x/2)
        spectrum = synchrotron_spectrum_w_hz(frequency, gamma, field, pitch)
        integral = np.sum(spectrum*frequency*weight)*(upper-lower)/2
        # Derive the bolometric power independently from Lorentz curvature.
        beta = np.sqrt(1-gamma**-2)
        acceleration = E_CHARGE_C*field*beta*np.sin(pitch)/(gamma*M_E_KG)
        expected = E_CHARGE_C**2*gamma**4*acceleration**2/(6*np.pi*EPSILON_0_F_M*C_M_S)
        self.assertAlmostEqual(float(integral)/expected, 1., delta=2e-8)

    def test_small_x_kernel_and_broad_spectrum(self):
        value = 1e-14
        asymptote = 2**(2/3)*gamma_function(2/3)*value**(1/3)
        self.assertAlmostEqual(float(synchrotron_f(value))/asymptote, 1., delta=5e-9)
        self.assertAlmostEqual(float(synchrotron_f(1.)), .651422815355364, delta=2e-12)
        self.assertEqual(float(synchrotron_f(0.)), 0.)
        self.assertEqual(float(synchrotron_f(900.)), 0.)

    def test_optical_band_fraction_matches_direct_frequency_integral(self):
        gamma, field, pitch = 800., .04, 1.
        low, high = C_M_S/700e-9, C_M_S/400e-9
        power = float(synchrotron_power_w(gamma, field, pitch))
        direct = quad(lambda frequency: float(synchrotron_spectrum_w_hz([frequency], gamma, field, pitch)[0]),
                      low, high, epsabs=1e-100, epsrel=2e-9)[0]/power
        fraction = float(synchrotron_wavelength_band_fraction(gamma, field, pitch))
        self.assertAlmostEqual(fraction/direct, 1., delta=5e-11)
        self.assertGreater(fraction, 0.)
        self.assertLess(fraction, 1.)
        fine = float(synchrotron_wavelength_band_fraction(gamma, field, pitch, nodes=128))
        self.assertAlmostEqual(fraction/fine, 1., delta=5e-12)

    def test_zero_emission_and_broadcasting(self):
        frequency = [0., 4e14, 5e14]
        spectrum = synchrotron_spectrum_w_hz(frequency, [300., 400.], [0., .01], [1., 0.])
        self.assertEqual(spectrum.shape, (2, 3))
        assert_allclose(spectrum, 0., atol=0.)
        assert_allclose(synchrotron_wavelength_band_fraction([300., 400.], [0., .01], [1., np.pi]),
                        0., atol=0.)

    def test_unsupported_emission_regimes_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'continuum'):
            synchrotron_spectrum_w_hz([5e14], 100., .01, 1e-4)
        with self.assertRaisesRegex(ValueError, 'Quantum'):
            synchrotron_power_w(100., B_QUANTUM_T, np.pi/2)
        with self.assertRaises(ValueError):
            synchrotron_spectrum_w_hz([-1.], 100., .01, 1.)
        with self.assertRaises(ValueError):
            synchrotron_wavelength_band_fraction(100., .01, 1., 700e-9, 400e-9)


if __name__ == '__main__':
    unittest.main()
