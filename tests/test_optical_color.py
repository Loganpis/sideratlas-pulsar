from pathlib import Path
import tempfile
import unittest

import numpy as np

from pulsar.optical_color import (
    BAND_NM, CIE1931_SHA256, C_M_PER_S, LUMINOUS_EFFICACY_LM_PER_W,
    XYZ_TO_LINEAR_SRGB, cie1931_cmf, frequency_density_to_per_nm,
    load_cie1931, spectrum_to_srgb, spectrum_to_xyz, xyz_to_srgb,
)


class OpticalColorChecks(unittest.TestCase):
    def test_official_table_and_hash_are_required(self):
        table = load_cie1931()
        self.assertEqual(table.shape, (471, 4))
        self.assertFalse(table.flags.writeable)
        np.testing.assert_array_equal(table[:, 0], np.arange(360, 831))
        np.testing.assert_allclose(cie1931_cmf(555), [.5120501, 1.0, .005749999],
                                   rtol=1e-7)
        self.assertEqual(len(CIE1931_SHA256), 64)
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / "altered.csv"
            altered.write_text("400,1,1,1\n700,1,1,1\n")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                load_cie1931(altered)

    def test_narrow_line_limit_recovers_cmf_without_arbitrary_color(self):
        # A normalized symmetric triangle wholly within one CIE linear
        # interval integrates to its centre CMF exactly.
        centre, halfwidth = 555.25, .001
        wavelength = [400, centre - halfwidth, centre, centre + halfwidth, 700]
        intensity = [0, 0, 1 / halfwidth, 0, 0]
        xyz = spectrum_to_xyz(wavelength, intensity)
        np.testing.assert_allclose(xyz,
                                   LUMINOUS_EFFICACY_LM_PER_W * cie1931_cmf(centre),
                                   rtol=3e-11, atol=1e-10)

    def test_red_green_blue_lines_produce_expected_channel_order(self):
        values = []
        for centre in (630.25, 530.25, 450.25):
            wavelength = [400, centre - .1, centre, centre + .1, 700]
            result = spectrum_to_srgb(wavelength, [0, 0, 10, 0, 0], exposure=1e-4)
            values.append(result.linear_srgb)
        self.assertEqual(int(np.argmax(values[0])), 0)
        self.assertEqual(int(np.argmax(values[1])), 1)
        self.assertEqual(int(np.argmax(values[2])), 2)

    def test_nonuniform_wavelengths_and_density_jacobian(self):
        wave = np.array([390., 401.25, 530.7, 699.9, 710.])
        intensity = 1e-5 * (2 + .003 * (wave - 400))
        nm = spectrum_to_xyz(wave, intensity)
        per_m = spectrum_to_xyz(wave * 1e-9, intensity * 1e9, wavelength_unit="m")
        per_angstrom = spectrum_to_xyz(wave * 10, intensity / 10,
                                       wavelength_unit="angstrom")
        np.testing.assert_allclose(nm, per_m, rtol=3e-15)
        np.testing.assert_allclose(nm, per_angstrom, rtol=3e-15)

        # Independent Gaussian quadrature over 1 nm CIE segments. The spectrum
        # is linear, so a three-point rule integrates this quadratic exactly.
        nodes, weights = np.polynomial.legendre.leggauss(3)
        sample = np.arange(400, 700)[:, None] + (nodes + 1)[None, :] / 2
        spectrum = 1e-5 * (2 + .003 * (sample - 400))
        cmf = cie1931_cmf(sample)
        expected = 683 * np.sum(spectrum[..., None] * cmf *
                                 weights[None, :, None] / 2, axis=(0, 1))
        np.testing.assert_allclose(nm, expected, rtol=3e-15)

    def test_frequency_density_jacobian(self):
        wave = np.array([400., 411.1, 505.4, 623.7, 700.])
        per_nm = np.array([1., 2., 3., 2., 1.]) * 1e-6
        per_hz = per_nm * (wave * 1e-9) ** 2 / (C_M_PER_S * 1e-9)
        converted = frequency_density_to_per_nm(wave, per_hz)
        np.testing.assert_allclose(converted, per_nm, rtol=2e-15)
        np.testing.assert_allclose(spectrum_to_xyz(wave, converted),
                                   spectrum_to_xyz(wave, per_nm), rtol=3e-15)

    def test_positive_homogeneity_and_zero_spectrum(self):
        wave = np.array([400., 450.5, 550.25, 700.])
        intensity = np.array([1., 2., .5, 3.]) * 1e-6
        xyz = spectrum_to_xyz(wave, intensity)
        np.testing.assert_allclose(spectrum_to_xyz(wave, 7.3 * intensity),
                                   7.3 * xyz, rtol=3e-14)
        np.testing.assert_array_equal(spectrum_to_xyz(wave, np.zeros(4)),
                                       np.zeros(3))

    def test_exact_band_excludes_outside_energy(self):
        wave = [360, 400, 500, 700, 830]
        base = spectrum_to_xyz(wave, [0, 1e-6, 1e-6, 1e-6, 0])
        bright_outside = spectrum_to_xyz(wave, [1e10, 1e-6, 1e-6, 1e-6, 1e10])
        np.testing.assert_array_equal(base, bright_outside)
        self.assertEqual(BAND_NM, (400., 700.))

    def test_fixed_exposure_preserves_linear_amplitude_and_frame_independence(self):
        wave = np.arange(400., 701.)
        spectrum = np.ones_like(wave) * 1e-6
        a = spectrum_to_srgb(wave, spectrum, exposure=.1)
        b = spectrum_to_srgb(wave, 2 * spectrum, exposure=.1)
        np.testing.assert_allclose(b.linear_srgb, 2 * a.linear_srgb, rtol=1e-15)
        combined = spectrum_to_srgb(wave, np.stack([spectrum, spectrum * 1e8]),
                                    exposure=.1)
        np.testing.assert_array_equal(combined.srgb[0], a.srgb)
        self.assertFalse(a.metadata()["automatic_normalization"])
        self.assertFalse(a.metadata()["white_balance_or_chromatic_adaptation"])
        self.assertEqual(a.exposure, .1)

    def test_d65_primaries_and_clipping_are_explicit(self):
        # XYZ of sRGB D65 white from its defined chromaticity.
        d65 = np.array([.3127 / .3290, 1, (1 - .3127 - .3290) / .3290])
        result = xyz_to_srgb(d65, exposure=1)
        np.testing.assert_allclose(result.linear_srgb, np.ones(3), atol=2e-15)
        out = xyz_to_srgb([0, 10, 0], exposure=1)
        self.assertTrue(np.any(out.clipped_low))
        self.assertTrue(np.any(out.clipped_high))
        np.testing.assert_array_equal(out.display_linear_srgb, [0, 1, 0])
        np.testing.assert_allclose(out.srgb, [0, 1, 0], atol=2e-16)
        np.testing.assert_array_equal(out.xyz, [0, 10, 0])
        np.testing.assert_array_equal(out.linear_srgb,
                                       XYZ_TO_LINEAR_SRGB @ [0, 10, 0])
        self.assertEqual(out.metadata()["clipped_low_channels"], 2)
        self.assertEqual(out.metadata()["clipped_high_channels"], 1)

    def test_reject_missing_support_instead_of_extrapolating(self):
        for wave in ([401, 700], [400, 699], [450, 550], [400.00001, 700]):
            with self.subTest(wave=wave):
                with self.assertRaisesRegex(ValueError, "full 400--700"):
                    spectrum_to_xyz(wave, [1, 1])
        for wavelength in (359.99, 830.01):
            with self.assertRaises(ValueError):
                cie1931_cmf(wavelength)

    def test_reject_invalid_physical_samples_and_exposure(self):
        invalid = [([400, 400, 700], [1, 1, 1]), ([700, 400], [1, 1]),
                   ([400, 700], [1, -1]), ([400, 700], [1, np.nan]),
                   ([400, 700], [1]), ([400, np.inf], [1, 1]),
                   ([[400, 700]], [1, 1]), ([400, 700], [1j, 1])]
        for wave, intensity in invalid:
            with self.subTest(wave=wave, intensity=intensity):
                with self.assertRaises(ValueError):
                    spectrum_to_xyz(wave, intensity)
        with self.assertRaises(ValueError):
            spectrum_to_xyz([400, 700], [1, 1], wavelength_unit="Hz")
        for exposure in (0, -1, np.inf, [1, 2], 1j):
            with self.subTest(exposure=exposure):
                with self.assertRaises(ValueError):
                    xyz_to_srgb([1, 1, 1], exposure=exposure)


if __name__ == "__main__":
    unittest.main()
