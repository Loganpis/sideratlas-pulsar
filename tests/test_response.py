"""Instrument response tests independent of Vela profile fitting."""
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma

from pulsar.radiation import curvature_photon_yield
from pulsar.response import (load_lat_response, incidence_scenarios, CurvatureResponse,
                             _cdf_asymmetric)


ROOT = Path(__file__).resolve().parents[1]
BANDS = ((.1, .3), (.3, 1.), (1., 3.), (3., 100.))


class ResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.response = load_lat_response()
        mu, w = incidence_scenarios(4)['uniform_cosine']
        cls.fold = CurvatureResponse(cls.response, BANDS, mu, w,
                                     intervals=128, critical_energy_nodes=2401)

    def test_calibration_integrity(self):
        provenance = json.loads((ROOT/'data/response/provenance.json').read_text())
        for info in provenance['files']:
            self.assertEqual(hashlib.sha256((ROOT/'data/response'/info['file']).read_bytes()).hexdigest(), info['sha256'])

    def test_area_units_sum_and_nonflatness(self):
        response = self.response
        from astropy.io import fits
        with fits.open(ROOT/'data/response/aeff_P8R3_SOURCE_V3_FB.fits') as f:
            t = f['EFFECTIVE AREA_FRONT'].data[0]
            e = np.sqrt(float(t['ENERG_LO'][30])*float(t['ENERG_HI'][30]))/1000
            mu = (float(t['CTHETA_LO'][20])+float(t['CTHETA_HI'][20]))/2
            exact = 10000 * float(t['EFFAREA'][20, 30])
        self.assertAlmostEqual(float(response.effective_area_cm2(e, mu, 'FRONT')), exact, places=8)
        e = np.array([.1, 1., 10.])
        total = response.effective_area_cm2(e, .8)
        np.testing.assert_allclose(total, response.effective_area_cm2(e, .8, 'FRONT')+response.effective_area_cm2(e, .8, 'BACK'))
        self.assertGreater(total[1], 2*total[0])

    def test_asymmetric_cdf_independent_integration(self):
        s, k, b, p = 1.8, .63, -.17, 1.55
        def density(x):
            u = x - b
            z = k*u/s if u >= 0 else -u/(k*s)
            return p/s/gamma(1/p)*k/(1+k*k)*np.exp(-z**p)
        for x in (-5., -1., b, .4, 3.):
            reference = quad(density, -np.inf, x, epsabs=1e-11)[0]
            self.assertAlmostEqual(float(_cdf_asymmetric(x, s, k, b, p)), reference, places=8)
        self.assertAlmostEqual(quad(density, -np.inf, np.inf, epsabs=1e-11)[0], 1., places=8)

    def test_dispersion_partition_and_energy_migration(self):
        e = np.array([.1, .3, 1., 3., 30.])
        bands = ((0., .1), (.1, .3), (.3, 1.), (1., 3.), (3., 100.), (100., np.inf))
        for part in ('FRONT', 'BACK'):
            p = self.response.band_probability(e, .7, bands, part)
            total = self.response.band_probability(e, .7, ((0., np.inf),), part)[:, 0]
            np.testing.assert_allclose(p.sum(axis=-1), total, atol=2e-15)
            self.assertTrue(np.all(total > .999))
            self.assertTrue(0 < p[1, 1] < 1 and 0 < p[1, 2] < 1)

    def test_distribution_interpolation_not_parameter_interpolation(self):
        # At fixed scaled energy, a midpoint must equal the mixture of the
        # neighbouring *distributions*, even when shape parameters differ.
        part = 'FRONT'
        le, mu, pars, c = self.response.dispersion[part]
        j, i = 3, 10
        log_e, m = (le[i]+le[i+1])/2, mu[j]
        e = 10**log_e
        l = log_e+3
        s = c[0]*l*l+c[1]*m*m+c[2]*l+c[3]*m+c[4]*l*m+c[5]
        xlo, xhi = -.5, 1.4
        actual = self.response.band_probability(e, m, ((e*(1+s*xlo), e*(1+s*xhi)),), part)[0]
        def one(p):
            return (p[0]*(_cdf_asymmetric(xhi,p[1],p[2],p[3],p[7])-_cdf_asymmetric(xlo,p[1],p[2],p[3],p[7]))
                    +(1-p[0])*(_cdf_asymmetric(xhi,p[5],p[6],p[4],p[8])-_cdf_asymmetric(xlo,p[5],p[6],p[4],p[8])))
        self.assertAlmostEqual(float(actual), .5*(one(pars[j,i])+one(pars[j,i+1])), places=12)

    def test_flat_response_limit(self):
        class FlatResponse:
            energy_limits_gev = (.01, 1000.)
            def band_area_cm2(self, e, bands, mu, w, energy_dispersion):
                return 1234. * ((e[:, None] >= bands[:, 0]) & (e[:, None] < bands[:, 1]))
        f = CurvatureResponse(FlatResponse(), BANDS, [.8], [1], energy_dispersion=False,
                              intervals=128, critical_energy_nodes=500)
        ec = np.array([.03, .7, 2.8, 15.])
        emitted = np.array([.4, 2., .8, .12])
        expected = 1234*curvature_photon_yield(emitted, ec, BANDS)
        np.testing.assert_allclose(f.direct_yield(emitted, ec), expected, rtol=5e-4, atol=1e-12)

    def test_fold_lookup_and_refinement(self):
        ec = np.array([.012, .113, .779, 1.92, 9.03, 55.1, 213.])
        direct = self.fold.direct_yield(np.ones_like(ec), ec)
        lookup = self.fold.photon_yield(np.ones_like(ec), ec)
        # Faint exponential tails require relative allowance; they contribute
        # negligibly to flux. Test meaningful response above 1e-8 of maximum.
        significant = direct > direct.max()*1e-8
        np.testing.assert_allclose(lookup[significant], direct[significant], rtol=.002)
        refined = CurvatureResponse(self.response, BANDS, self.fold.cos_theta, self.fold.weights,
                                     intervals=256, critical_energy_nodes=500)
        np.testing.assert_allclose(refined.direct_yield(np.ones_like(ec), ec), direct, rtol=.001, atol=1e-9)

    def test_incidence_weights_and_domains(self):
        for mu, w in incidence_scenarios().values():
            self.assertAlmostEqual(w.sum(), 1.)
            self.assertTrue(np.all((mu >= .2) & (mu <= 1)))
        for e, mu in ((.001, .8), (4000., .8), (1., .1)):
            with self.assertRaises(ValueError):
                self.response.effective_area_cm2(e, mu)
        np.testing.assert_equal(self.fold.photon_yield([0., 1.], [1., 0.]), np.zeros((2,4)))
        with self.assertRaises(ValueError):
            self.fold.photon_yield(1., -1.)
        with self.assertRaises(ValueError):
            self.response.band_area_cm2([1.], BANDS, [.8], [-1.])


if __name__ == '__main__':
    unittest.main()
