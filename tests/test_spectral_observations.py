from pathlib import Path
import unittest
import numpy as np
from pulsar.spectral_observations import load_vela_spectrum,catalog_reference_sed,compare_integrated_spectrum,MEV_ERG

DATA=Path(__file__).resolve().parents[1]/'data/spectra/3PC_Catalog+SEDs_20230803.fits'

class SpectralChecks(unittest.TestCase):
    def test_real_catalog_and_units(self):
        d=load_vela_spectrum(DATA)
        self.assertEqual(d['source_name'],'4FGL J0835.3-4510')
        self.assertEqual(d['is_upper_limit'].tolist(),[False]*7+[True])
        self.assertAlmostEqual(d['period_s'],.08937108859,places=10)
        self.assertGreater(d['period_derivative'],1e-13)
        self.assertAlmostEqual(d['nufnu_erg_cm2_s'][3]/MEV_ERG,.001600264,places=9)
        p=d['plec_bfree']; e=p['Pivot_Energy']
        self.assertAlmostEqual(float(catalog_reference_sed(e/1000,d))/(e**2*p['PLEC_Flux_Density']*MEV_ERG),1,places=12)

    def test_one_normalization_recovered_without_fitting_upper_limit(self):
        d=load_vela_spectrum(DATA);m=d['flux_ph_cm2_s']*27
        before=compare_integrated_spectrum(m,d)
        self.assertAlmostEqual(before['normalization'],1/27,places=9)
        self.assertLess(before['descriptive_chi_squared'],1e-12)
        m[-1]=1e9
        after=compare_integrated_spectrum(m,d)
        self.assertAlmostEqual(after['normalization'],before['normalization'],places=8)
        self.assertTrue(after['upper_limit_exceeded'])

    def test_high_energy_excess_cannot_be_erased_by_per_band_scaling(self):
        d=load_vela_spectrum(DATA);m=d['flux_ph_cm2_s'].copy();m[6]*=2
        s=compare_integrated_spectrum(m,d)
        self.assertGreater(s['descriptive_chi_squared'],50)
        self.assertEqual(s['descriptive_dof'],6)

    def test_absolute_flux_mode_and_dark_model(self):
        d=load_vela_spectrum(DATA)
        s=compare_integrated_spectrum(d['flux_ph_cm2_s']*2,d,fit_normalization=False)
        self.assertEqual(s['normalization'],1)
        self.assertEqual(s['descriptive_dof'],7)
        self.assertGreater(s['descriptive_chi_squared'],100)
        s=compare_integrated_spectrum(np.zeros(8),d)
        self.assertTrue(np.isfinite(s['descriptive_chi_squared']))
        self.assertGreater(s['descriptive_chi_squared'],100)
