"""Independent limits, energy moments and convergence for radiation kernels."""
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.special import zeta

from pulsar.nonthermal import (
    B_QUANTUM_GAUSS, K_B_ERG_K, SIGMA_T_CM2, PhotonField,
    inverse_compton_rates, inverse_compton_spectrum_per_gev,
    isotropized_stellar_photon_field, stellar_dilution,
    synchrotron_critical_energy_gev, synchrotron_photon_rate,
    synchrotron_power_erg_s, thermal_photon_field,
)
from pulsar.radiation import C_CM_S, E_STATC, GEV_ERG, HBAR_ERG_S, M_E_G, MEC2_ERG


class SynchrotronTests(unittest.TestCase):
    def test_fixed_pitch_power_from_relativistic_larmor_formula(self):
        g, b, pitch = 1e4, 120., .2
        expected = 2*E_STATC**4*b*b*(g*g-1)*np.sin(pitch)**2/(3*M_E_G**2*C_CM_S**3)
        self.assertAlmostEqual(float(synchrotron_power_erg_s(g,b,pitch))/expected, 1., delta=2e-15)
        ec = 3*HBAR_ERG_S*E_STATC*b*g*g*np.sin(pitch)/(2*M_E_G*C_CM_S)/GEV_ERG
        self.assertAlmostEqual(float(synchrotron_critical_energy_gev(g,b,pitch))/ec, 1., delta=2e-15)

    def test_isotropic_pitch_average_has_four_thirds_coefficient(self):
        mu, weight = np.polynomial.legendre.leggauss(64)
        g, b = 1e4, 50.
        average = np.sum(weight*synchrotron_power_erg_s(g,b,np.arccos(mu)))/2
        expected = 4/3*SIGMA_T_CM2*C_CM_S*(b*b/(8*np.pi))*(g*g-1)
        self.assertAlmostEqual(average/expected, 1., delta=2e-14)

    def test_zero_pitch_and_zero_field_emit_nothing(self):
        assert_allclose(synchrotron_power_erg_s(1e4,10.,[0.,np.pi]), 0., atol=0)
        assert_allclose(synchrotron_photon_rate(1e4,10.,[0.,np.pi],[[0,np.inf]]), 0., atol=0)
        self.assertEqual(float(synchrotron_power_erg_s(1e4,0.,.5)), 0.)

    def test_band_photons_reproduce_independent_bolometric_power(self):
        g, b, pitch = 1e4, 100., .7
        ec = float(synchrotron_critical_energy_gev(g,b,pitch))
        # Numerically sum emitted photon energies, not the production power
        # primitive; broad log bands resolve the complete synchrotron shape.
        edges = ec*np.geomspace(1e-9, 100., 5001)
        rates = synchrotron_photon_rate(g,b,pitch,np.column_stack((edges[:-1],edges[1:])))
        power = np.sum(rates*np.sqrt(edges[:-1]*edges[1:])*GEV_ERG)
        expected = 2*E_STATC**4*b*b*(g*g-1)*np.sin(pitch)**2/(3*M_E_G**2*C_CM_S**3)
        self.assertAlmostEqual(power/expected, 1., delta=5e-6)

    def test_unsupported_quantum_and_cyclotron_regimes_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Quantum'):
            synchrotron_power_erg_s(1e6,B_QUANTUM_GAUSS/1e5,np.pi/2)
        with self.assertRaisesRegex(ValueError, 'continuum'):
            synchrotron_photon_rate(1e3,100.,1e-5,[[.1,1.]])


class SeedPhotonTests(unittest.TestCase):
    def test_planck_number_and_energy_density(self):
        temperature, dilution = 8e5, .0007
        field = thermal_photon_field(temperature,dilution,nodes=128)
        expected_n = dilution*2*zeta(3)/np.pi**2*(K_B_ERG_K*temperature/(HBAR_ERG_S*C_CM_S))**3
        expected_u = dilution*np.pi**2/15*(K_B_ERG_K*temperature)**4/(HBAR_ERG_S*C_CM_S)**3
        self.assertAlmostEqual(field.number_density_cm3/expected_n,1.,delta=2e-9)
        self.assertAlmostEqual(field.energy_density_erg_cm3/expected_u,1.,delta=2e-9)

    def test_stellar_solid_angle_and_far_field_luminosity_limit(self):
        radius = 1.2e6
        self.assertEqual(float(stellar_dilution(radius,radius)),.5)
        self.assertAlmostEqual(float(stellar_dilution(radius*1e8,radius))/(.25e-16),1.,delta=1e-15)
        field = isotropized_stellar_photon_field(1e6,radius,1e9)
        self.assertIn('ISOTROPIZED',field.description)
        self.assertGreater(field.energy_density_erg_cm3,0.)
        with self.assertRaises(ValueError):
            stellar_dilution(.9*radius,radius)


class InverseComptonTests(unittest.TestCase):
    def test_thomson_total_rate_and_net_electron_loss(self):
        g, energy, density = 1000., 1e-12*MEC2_ERG, 1e7
        field = PhotonField([energy],[density])
        result = inverse_compton_rates(g,field,[[0,np.inf]],q_nodes=128)
        expected_rate = SIGMA_T_CM2*C_CM_S*density
        expected_loss = 4/3*SIGMA_T_CM2*C_CM_S*energy*density*(g*g-1)
        self.assertAlmostEqual(float(result.photon_rate_s)/expected_rate,1.,delta=3e-6)
        self.assertAlmostEqual(float(result.electron_loss_erg_s)/expected_loss,1.,delta=3e-6)
        self.assertAlmostEqual(float(result.band_photon_rate_s[0])/float(result.photon_rate_s),1.,delta=1e-14)
        assert_allclose(result.emitted_power_erg_s,
                        result.electron_loss_erg_s+result.intercepted_seed_power_erg_s,rtol=2e-15)
        self.assertGreater(float(result.intercepted_seed_power_erg_s),0.)

    def test_kn_power_against_independent_energy_quadrature(self):
        # Direct integration in x=E_gamma/E_e, independent of the production
        # log-q quadrature and its Jacobian. Includes the incident energy term.
        g, eps, density = 1e7, 2e-6*MEC2_ERG, 1e5
        recoil = 4*g*eps/MEC2_ERG
        xmin = eps/(g*MEC2_ERG+eps)
        xmax = recoil/(1+recoil)
        pref = 3*SIGMA_T_CM2*C_CM_S*density/(4*g*g*eps)
        def integrand(x):
            q = x/(recoil*(1-x))
            kernel = 2*q*np.log(q)+(1+2*q)*(1-q)+.5*(recoil*q)**2*(1-q)/(1+recoil*q)
            return pref*kernel*(x*g*MEC2_ERG-eps)*g*MEC2_ERG
        reference = quad(integrand,xmin,xmax,epsabs=1e-35,epsrel=1e-10)[0]
        result = inverse_compton_rates(g,PhotonField([eps],[density]),q_nodes=128)
        self.assertAlmostEqual(float(result.electron_loss_erg_s)/reference,1.,delta=2e-9)
        thomson = 4/3*SIGMA_T_CM2*C_CM_S*eps*density*g*g
        self.assertLess(float(result.electron_loss_erg_s)/thomson,.01)

    def test_band_additivity_spectrum_and_kinematic_cutoff(self):
        g = 2e5
        field = PhotonField([1e-5*MEC2_ERG],[1e6])
        emax = g*MEC2_ERG/GEV_ERG
        edges = emax*np.array([0,.01,.1,.5,1.])
        rates = inverse_compton_rates(g,field,np.column_stack((edges[:-1],edges[1:])),q_nodes=128)
        self.assertAlmostEqual(rates.band_photon_rate_s.sum()/float(rates.photon_rate_s),1.,delta=2e-9)
        direct = quad(lambda loge: float(inverse_compton_spectrum_per_gev([np.exp(loge)],g,field)[0])*np.exp(loge),
                      np.log(edges[1]),np.log(edges[2]),epsabs=1e-25,epsrel=1e-8)[0]
        self.assertAlmostEqual(direct/rates.band_photon_rate_s[1],1.,delta=2e-8)
        assert_allclose(inverse_compton_spectrum_per_gev([emax,emax*2],g,field),0.,atol=0)

    def test_thermal_quadrature_convergence_and_density_scaling(self):
        g = np.array([1e3,1e7])
        bands = [[.1,1.],[1.,100.],[100.,10000.]]
        coarse = inverse_compton_rates(g,thermal_photon_field(7e5,1e-5,96),bands,q_nodes=96)
        fine = inverse_compton_rates(g,thermal_photon_field(7e5,1e-5,192),bands,q_nodes=192)
        assert_allclose(coarse.electron_loss_erg_s,fine.electron_loss_erg_s,rtol=3e-6)
        # Band boundaries vary across seed energies; this is a stricter test
        # than comparing a single bolometric normalization.
        assert_allclose(coarse.band_photon_rate_s,fine.band_photon_rate_s,rtol=.006,atol=1e-100)
        doubled = inverse_compton_rates(g,thermal_photon_field(7e5,2e-5,192),bands,q_nodes=192)
        assert_allclose(doubled.electron_loss_erg_s,2*fine.electron_loss_erg_s,rtol=1e-14)
        assert_allclose(doubled.band_photon_rate_s,2*fine.band_photon_rate_s,rtol=1e-14)

    def test_empty_field_and_invalid_inputs(self):
        result = inverse_compton_rates(1e4,PhotonField([1e-12],[0.]),[[.1,1.]])
        assert_allclose(result.band_photon_rate_s,0.,atol=0)
        self.assertEqual(float(result.electron_loss_erg_s),0.)
        with self.assertRaises(ValueError):
            inverse_compton_rates(1.,PhotonField([1e-12],[1.]))
        with self.assertRaises(ValueError):
            PhotonField([0.],[1.])


if __name__ == '__main__':
    unittest.main()
