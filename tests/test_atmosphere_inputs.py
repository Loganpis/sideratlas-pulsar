"""Atmosphere coverage must fail closed; redshift laws retain spectral units."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest

import numpy as np
from astropy.io import fits

from pulsar.atmosphere_inputs import (
    INTENSITY_UNIT, audit_nsmaxg_fits, photon_energy_kev,
    schwarzschild_surface_factors, transform_bolometric_intensity,
    transform_specific_intensity_frequency,
    transform_specific_intensity_wavelength,
    validate_angledependent_atmosphere_table,
)


class AtmosphereInputTests(unittest.TestCase):
    def test_visible_photon_energy_units_and_order(self):
        energies = photon_energy_kev([400, 700])
        np.testing.assert_allclose(energies, [0.00309960496083, 0.00177120283476],
                                   rtol=1e-11)
        self.assertGreater(energies[0], energies[1])
        for value in (0, -1, np.inf, np.nan):
            with self.assertRaises(ValueError):
                photon_energy_kev(value)

    def test_surface_redshift_and_proper_gravity(self):
        surface = schwarzschild_surface_factors()
        self.assertAlmostEqual(surface.compactness_2gm_rc2, 0.3180415, delta=3e-6)
        self.assertAlmostEqual(surface.frequency_ratio_g, 0.8258078, delta=3e-6)
        self.assertAlmostEqual(surface.apparent_radius_km, 15.741, delta=.003)
        self.assertAlmostEqual(surface.redshift_z+1, 1/surface.frequency_ratio_g)
        self.assertGreater(surface.proper_surface_gravity_m_s2, 1e12)
        self.assertLess(surface.proper_surface_gravity_m_s2, 2e12)
        # A factor calculated for an assumed star is not an atmosphere map.
        self.assertNotIn("specific_intensity", surface.as_dict())
        for mass, radius in ((-1, 13), (1.4, 0), (1.4, 3), (np.nan, 13)):
            with self.assertRaises(ValueError):
                schwarzschild_surface_factors(mass, radius)
        # 5 km lies outside the horizon but inside the photon sphere for
        # this mass: R/g is not the visible silhouette in that regime.
        with self.assertRaisesRegex(ValueError, "photon sphere"):
            schwarzschild_surface_factors(1.4, 5.)

    def test_frequency_and_wavelength_invariance_agree(self):
        g = .81
        wavelength_em_nm = np.array([400., 550., 700.])
        intensity_nu_em = np.array([1., 2., .5])
        # I_lambda per nm = I_nu * |d nu / d lambda_nm|.
        jacobian_em = 299792458./(wavelength_em_nm*1e-9)**2*1e-9
        wavelength_obs_nm = wavelength_em_nm/g
        jacobian_obs = 299792458./(wavelength_obs_nm*1e-9)**2*1e-9
        from_frequency = transform_specific_intensity_frequency(intensity_nu_em, g)*jacobian_obs
        from_wavelength = transform_specific_intensity_wavelength(intensity_nu_em*jacobian_em, g)
        np.testing.assert_allclose(from_frequency, from_wavelength, rtol=2e-15)
        nu_em = 299792458./(wavelength_em_nm*1e-9)
        nu_obs = g*nu_em
        np.testing.assert_allclose(
            transform_specific_intensity_frequency(intensity_nu_em, g)/nu_obs**3,
            intensity_nu_em/nu_em**3, rtol=2e-15)

    def test_integrating_shifted_spectrum_gives_bolometric_g_four(self):
        # Arbitrary smooth numerical fixture, not a production atmosphere.
        wavelength = np.geomspace(1, 1e5, 12001)
        emitted = np.exp(-((np.log(wavelength)-6)/1.5)**2)
        g = .73
        emitted_integral = np.trapz(emitted, wavelength)
        observed_integral = np.trapz(
            transform_specific_intensity_wavelength(emitted, g), wavelength/g)
        self.assertAlmostEqual(observed_integral/emitted_integral, g**4, delta=2e-16)
        self.assertAlmostEqual(observed_integral,
                               transform_bolometric_intensity(emitted_integral, g),
                               delta=1e-11)
        for function in (transform_specific_intensity_frequency,
                         transform_specific_intensity_wavelength,
                         transform_bolometric_intensity):
            with self.assertRaises(ValueError):
                function(1, 0)
            with self.assertRaises(ValueError):
                function(-1, .8)

    def test_realistic_variable_length_nsmaxg_fits_audit(self):
        # Different-length PA/PE arrays reproduce the current official format.
        names = np.empty(2, dtype=object)
        names[:] = ["HB1200ThB00", "HB1230ThB00g1438"]
        energies = np.empty(2, dtype=object)
        energies[:] = [np.array([.047863, .1, 1, 10]),
                       np.array([.047863, .2, 10])]
        temperatures = np.empty(2, dtype=object)
        temperatures[:] = [np.array([5.5, 6.0]), np.array([5.5, 6.0, 6.8])]
        gravities = np.empty(2, dtype=object)
        gravities[:] = [np.array([13.6, 14.38]), np.array([14.38])]
        info = fits.BinTableHDU.from_columns([
            fits.Column(name="MODNUM", format="J", array=[1200, 1230]),
            fits.Column(name="MODNAME", format="PA()", array=names),
            fits.Column(name="TEMPS", format="PE()", array=temperatures),
            fits.Column(name="GRAVS", format="PE()", array=gravities),
            fits.Column(name="ENERGIES", format="PE()", array=energies),
        ], name="MODELINFO")
        spectra = []
        for name, energy in zip(names, energies):
            flux = np.empty(1, dtype=object)
            flux[0] = np.ones(len(energy))
            spectra.append(fits.BinTableHDU.from_columns([
                fits.Column(name="FLUX", format="PE()", array=flux)], name=name))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"synthetic-format-fixture.fits"
            fits.HDUList([fits.PrimaryHDU(), info, *spectra]).writeto(path)
            digest = sha256(path.read_bytes()).hexdigest()
            audit = audit_nsmaxg_fits(path, expected_sha256=digest)
            self.assertFalse(audit["any_model_covers_requested_wavelengths"])
            self.assertFalse(audit["local_angledependent_specific_intensity_available"])
            self.assertTrue(audit["expected_hash_verified"])
            self.assertEqual([row["energy_samples"] for row in audit["models"]], [4, 3])
            self.assertEqual(audit["models"][1]["name"], names[1])
            with self.assertRaisesRegex(ValueError, "SHA256"):
                audit_nsmaxg_fits(path, expected_sha256="0"*64)

    def _fixture(self, path):
        path.write_bytes(b"Clearly synthetic validator test fixture, not atmosphere model data")
        return {
            "wavelength_nm": np.array([350., 400., 550., 700., 750.]),
            "cos_emission_angle": np.array([0., .4, 1.]),
            "specific_intensity_w_m2_sr_nm": np.ones((3, 5)),
            "metadata": {
                "source_url": "https://example.invalid/synthetic-validator-fixture",
                "citation": "Synthetic unit-test fixture; no astrophysical claim",
                "source_sha256": sha256(path.read_bytes()).hexdigest(),
                "source_kind": "synthetic_test_fixture",
                "intensity_unit": INTENSITY_UNIT, "polarization": "total",
                "temperature_k": 1e6, "magnetic_field_t": 2e8,
                "magnetic_inclination_deg": 0.,
                "surface_gravity_m_s2": 1.3e12,
                "azimuthally_symmetric": True,
                "model_assumptions": ["Artificial values for structural validation only"],
            },
        }

    def test_synthetic_fixture_cannot_become_a_production_atmosphere(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"fixture.dat"
            table = self._fixture(path)
            with self.assertRaisesRegex(ValueError, "production"):
                validate_angledependent_atmosphere_table(table, source_path=path)
            report = validate_angledependent_atmosphere_table(
                table, source_path=path, allow_test_fixture=True)
            self.assertTrue(report["structurally_valid"])
            self.assertFalse(report["structurally_eligible_for_review"])
            self.assertFalse(report["scientific_validation"])
            self.assertNotIn("production_ready", report)
            self.assertFalse(report["extrapolation"])
            path.write_bytes(b"altered bytes")
            with self.assertRaisesRegex(ValueError, "SHA256"):
                validate_angledependent_atmosphere_table(
                    table, source_path=path, allow_test_fixture=True)

    def test_reject_flux_units_missing_provenance_and_optical_extrapolation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"fixture.dat"
            original = self._fixture(path)
            cases = []
            table = deepcopy(original)
            table["metadata"]["intensity_unit"] = "W m-2 nm-1"
            cases.append(table)
            table = deepcopy(original)
            del table["metadata"]["citation"]
            cases.append(table)
            table = deepcopy(original)
            table["wavelength_nm"] = np.array([10., 20., 30., 40., 50.])
            cases.append(table)
            table = deepcopy(original)
            table["cos_emission_angle"] = np.array([.1, .4, 1.])
            cases.append(table)
            table = deepcopy(original)
            table["specific_intensity_w_m2_sr_nm"][0, 0] = np.nan
            cases.append(table)
            table = deepcopy(original)
            table["metadata"]["magnetic_inclination_deg"] = 45
            cases.append(table)
            for table in cases:
                with self.subTest(table=table):
                    with self.assertRaises(ValueError):
                        validate_angledependent_atmosphere_table(
                            table, source_path=path, allow_test_fixture=True)

    def test_provenance_claim_and_hash_do_not_validate_model_science(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"fixture.dat"
            table = self._fixture(path)
            # A source-kind claim cannot be authenticated by an array/schema
            # checker. Passing it admits review, never a scientific claim.
            table["metadata"]["source_kind"] = "published_model"
            report = validate_angledependent_atmosphere_table(table, source_path=path)
            self.assertTrue(report["structurally_eligible_for_review"])
            self.assertFalse(report["scientific_validation"])
            self.assertNotIn("production_ready", report)

    def test_inclined_field_needs_periodic_azimuth_resolved_intensity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"fixture.dat"
            table = self._fixture(path)
            table["metadata"]["azimuthally_symmetric"] = False
            table["metadata"]["magnetic_inclination_deg"] = 45
            table["azimuth_deg"] = np.array([0., 180., 360.])
            table["specific_intensity_w_m2_sr_nm"] = np.ones((3, 3, 5))
            report = validate_angledependent_atmosphere_table(
                table, source_path=path, allow_test_fixture=True)
            self.assertEqual(report["shape"], [3, 3, 5])
            table["specific_intensity_w_m2_sr_nm"][-1] *= 2
            with self.assertRaisesRegex(ValueError, "periodic"):
                validate_angledependent_atmosphere_table(
                    table, source_path=path, allow_test_fixture=True)

    def test_observer_band_requires_redshifted_local_wavelength_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"fixture.dat"
            table = self._fixture(path)
            g = schwarzschild_surface_factors().frequency_ratio_g
            # A local 400..700 nm table does not cover the bluer local
            # photons observed at 400 nm after gravitational redshift.
            table["wavelength_nm"] = np.array([400., 450., 550., 650., 700.])
            with self.assertRaisesRegex(ValueError, "local wavelengths"):
                validate_angledependent_atmosphere_table(
                    table, source_path=path, allow_test_fixture=True,
                    frequency_ratio_g=g)
            table["wavelength_nm"] = np.array([300., 400., 550., 700., 750.])
            report = validate_angledependent_atmosphere_table(
                table, source_path=path, allow_test_fixture=True,
                frequency_ratio_g=g)
            np.testing.assert_allclose(report["required_local_wavelength_range_nm"],
                                       [330.3229169433525, 578.0651046508668],
                                       rtol=2e-15)
            self.assertEqual(report["observer_wavelength_range_nm"], [400., 700.])
            self.assertFalse(report["scientific_validation"])
            for invalid_g in (0., -1., np.nan, [g, g]):
                with self.assertRaises(ValueError):
                    validate_angledependent_atmosphere_table(
                        table, source_path=path, allow_test_fixture=True,
                        frequency_ratio_g=invalid_g)


if __name__ == "__main__":
    unittest.main()
