"""Band-limited physical spectra -> CIE XYZ -> fixed-exposure sRGB.

Inputs are observer-frame spectral radiance, not particle energy or a scalar
emissivity map. The 400--700 nm band is explicit and is not the full range of
human vision. This module supplies color integration, not a Vela emission,
atmosphere, eye-adaptation, telescope or detector model.

Both the supplied spectrum and official 1 nm CIE tables are interpolated
linearly in wavelength. Their product is integrated analytically on every
interval of their combined knot grid. No missing wavelength is extrapolated.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
from pathlib import Path

import numpy as np


BAND_NM = (400.0, 700.0)
LUMINOUS_EFFICACY_LM_PER_W = 683.0
C_M_PER_S = 299_792_458.0
CIE1931_SHA256 = "fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1"
CIE1931_DOI = "https://doi.org/10.25039/CIE.DS.xvudnb9b"
CIE1931_PATH = (Path(__file__).resolve().parents[2] /
                "data/color/CIE_xyz_1931_2deg.csv")

# D65-referenced sRGB primaries, W3C CSS Color 4 section 19 rational matrix.
# No chromatic adaptation or source-dependent white balance is applied.
XYZ_TO_LINEAR_SRGB = np.array([
    [12831 / 3959, -329 / 214, -1974 / 3959],
    [-851781 / 878810, 1648619 / 878810, 36519 / 878810],
    [705 / 12673, -2585 / 12673, 705 / 667],
], dtype=float)
XYZ_TO_LINEAR_SRGB.setflags(write=False)


def _real_array(value, name):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real")
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


@lru_cache(maxsize=4)
def load_cie1931(path=None):
    """Return readonly wavelength/xbar/ybar/zbar columns after hash checking.

    The unmodified table is CIE 2019, CC BY-SA 4.0; see data/color/PROVENANCE.md.
    A custom path must contain the same official bytes, not a fitted substitute.
    """
    path = CIE1931_PATH if path is None else Path(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != CIE1931_SHA256:
        raise ValueError("CIE table does not match the official SHA-256")
    table = np.loadtxt(path, delimiter=",")
    if (table.shape != (471, 4) or
            not np.array_equal(table[:, 0], np.arange(360, 831)) or
            not np.all(np.isfinite(table)) or np.any(table[:, 1:] < 0)):
        raise ValueError("Invalid CIE 1931 table")
    table.setflags(write=False)
    return table


def cie1931_cmf(wavelength_nm):
    """Linear interpolation of official CMFs; no out-of-table extrapolation."""
    wavelength = _real_array(wavelength_nm, "wavelength_nm")
    table = load_cie1931()
    if np.any(wavelength < table[0, 0]) or np.any(wavelength > table[-1, 0]):
        raise ValueError("CMF wavelengths must lie within 360--830 nm")
    return np.stack([np.interp(wavelength, table[:, 0], table[:, j])
                     for j in (1, 2, 3)], axis=-1)


def _validated_spectrum(wavelength, spectral_radiance, wavelength_unit):
    # The density is per selected wavelength unit, so its Jacobian must change
    # inversely with the wavelength axis. W m^-2 sr^-1 m^-1 / 1e9 = per nm.
    factors = {"nm": 1.0, "m": 1e9, "angstrom": 0.1}
    if not isinstance(wavelength_unit, str) or wavelength_unit not in factors:
        raise ValueError("wavelength_unit must be 'nm', 'm' or 'angstrom'")
    factor = factors[wavelength_unit]
    wave = _real_array(wavelength, "wavelength") * factor
    intensity = _real_array(spectral_radiance, "spectral_radiance") / factor
    if not np.all(np.isfinite(wave)) or not np.all(np.isfinite(intensity)):
        raise ValueError("Unit conversion produced nonfinite spectral samples")
    if wave.ndim != 1 or wave.size < 2 or np.any(wave <= 0):
        raise ValueError("wavelength must be a positive 1D array with >=2 samples")
    if np.any(np.diff(wave) <= 0):
        raise ValueError("wavelength samples must be strictly increasing")
    if intensity.ndim < 1 or intensity.shape[-1] != wave.size:
        raise ValueError("The last radiance axis must match wavelength")
    if np.any(intensity < 0):
        raise ValueError("Physical spectral radiance must be nonnegative")
    # Only conversion-roundoff-sized deviations at exact band endpoints are
    # snapped. This does not permit filling an actually unobserved band edge.
    wave = wave.copy()
    for edge in BAND_NM:
        wave[np.isclose(wave, edge, rtol=0, atol=1e-10)] = edge
    if np.any(np.diff(wave) <= 0):
        raise ValueError("Wavelength samples coincide at conversion precision")
    if wave[0] > BAND_NM[0] or wave[-1] < BAND_NM[1]:
        raise ValueError("Spectrum must cover the full 400--700 nm band")
    return wave, intensity


def spectrum_to_xyz(wavelength, spectral_radiance, *, wavelength_unit="nm"):
    """Integrate observer-frame radiance over exactly 400--700 nm.

    `spectral_radiance` has shape (..., N), units W m^-2 sr^-1 per wavelength
    unit. `wavelength_unit` specifies BOTH the wavelength axis and density
    denominator. With 'nm', the density is per nm; with 'm', it is per metre.

    X,Y,Z = 683 integral I_lambda(lambda) [xbar,ybar,zbar] d lambda

    Returned shape is (..., 3). Y is band-limited photopic luminance in cd/m²;
    X and Z use the corresponding photometric tristimulus scale. These are not
    exposure-normalized values. Piecewise-linear interpolation is an explicit
    sampling assumption: unresolved spectral lines cannot be recovered.
    """
    wave, intensity = _validated_spectrum(wavelength, spectral_radiance,
                                           wavelength_unit)
    table = load_cie1931()
    cie_knots = table[(table[:, 0] >= BAND_NM[0]) &
                      (table[:, 0] <= BAND_NM[1]), 0]
    knots = np.unique(np.concatenate((cie_knots,
                                     wave[(wave > BAND_NM[0]) &
                                          (wave < BAND_NM[1])], BAND_NM)))
    left = np.clip(np.searchsorted(wave, knots, side="right") - 1,
                   0, wave.size - 2)
    f = (knots - wave[left]) / (wave[left + 1] - wave[left])
    sampled = ((1 - f) * np.take(intensity, left, axis=-1) +
               f * np.take(intensity, left + 1, axis=-1))
    cmf = cie1931_cmf(knots)
    dx = np.diff(knots)[:, None]
    weights = np.zeros((knots.size, 3))
    # Exact integral of two linear functions, not trapezoidal integration of
    # their quadratic product: dx/6 [I0(2C0+C1)+I1(C0+2C1)].
    weights[:-1] += dx * (2 * cmf[:-1] + cmf[1:]) / 6
    weights[1:] += dx * (cmf[:-1] + 2 * cmf[1:]) / 6
    with np.errstate(over="raise", invalid="raise"):
        return LUMINOUS_EFFICACY_LM_PER_W * np.einsum("...n,nk->...k",
                                                       sampled, weights)


def frequency_density_to_per_nm(wavelength_nm, spectral_radiance_per_hz):
    """Convert paired I_nu(c/lambda) values to I_lambda per nm.

    I_lambda,nm = I_nu * c/lambda_m² * 1e-9. Wavelengths must be provided at
    the actual frequency sample locations; this function never invents samples
    or changes their order. The returned spectrum still needs full band support.
    """
    wave = _real_array(wavelength_nm, "wavelength_nm")
    intensity = _real_array(spectral_radiance_per_hz, "spectral_radiance_per_hz")
    if wave.ndim != 1 or wave.size < 1 or np.any(wave <= 0):
        raise ValueError("wavelength_nm must be a positive 1D sample array")
    if intensity.ndim < 1 or intensity.shape[-1] != wave.size:
        raise ValueError("The last radiance axis must match wavelength")
    if np.any(intensity < 0):
        raise ValueError("Physical spectral radiance must be nonnegative")
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        return intensity * C_M_PER_S * 1e-9 / (wave * 1e-9) ** 2


@dataclass(frozen=True)
class SpectralColor:
    xyz: np.ndarray
    linear_srgb: np.ndarray
    display_linear_srgb: np.ndarray
    srgb: np.ndarray
    clipped_low: np.ndarray
    clipped_high: np.ndarray
    exposure: float

    def metadata(self):
        return {
            "band_nm": list(BAND_NM),
            "observer": "CIE 1931 2-degree photopic standard observer",
            "cie_dataset_doi": CIE1931_DOI,
            "cie_dataset_sha256": CIE1931_SHA256,
            "xyz_scale": "683 lm/W; Y is band-limited luminance in cd/m^2",
            "display": "sRGB, D65-referenced primaries",
            "matrix_xyz_to_linear_srgb": XYZ_TO_LINEAR_SRGB.tolist(),
            "exposure": self.exposure,
            "exposure_units": "inverse photometric tristimulus scale",
            "clipping": "independent linear RGB channels clipped to [0,1]",
            "clipped_low_channels": int(np.count_nonzero(self.clipped_low)),
            "clipped_high_channels": int(np.count_nonzero(self.clipped_high)),
            "automatic_normalization": False,
            "white_balance_or_chromatic_adaptation": False,
            "tone_mapping": "none",
            "scope": "band-limited spectral color conversion; no Vela prediction",
        }


def xyz_to_srgb(xyz, *, exposure):
    """Apply ONE fixed scalar exposure, then sRGB channel clipping/encoding.

    The caller must retain the same exposure across all pixels and movie frames.
    It is a display setting, not a physically modeled eye/camera exposure.
    Negative RGB values identify out-of-gamut spectra and are retained in
    `linear_srgb` and `clipped_low`; clipping can change hue and saturation.
    """
    xyz = _real_array(xyz, "xyz")
    if xyz.ndim < 1 or xyz.shape[-1] != 3 or np.any(xyz < 0):
        raise ValueError("xyz must have a nonnegative final axis of length 3")
    exp = _real_array(exposure, "exposure")
    if exp.ndim != 0 or float(exp) <= 0:
        raise ValueError("exposure must be one finite positive scalar")
    exposure = float(exp)
    with np.errstate(over="raise", invalid="raise"):
        linear = exposure * np.einsum("...j,ij->...i", xyz,
                                      XYZ_TO_LINEAR_SRGB)
    clipped = np.clip(linear, 0, 1)
    srgb = np.where(clipped <= 0.0031308, 12.92 * clipped,
                     1.055 * clipped ** (1 / 2.4) - 0.055)
    return SpectralColor(xyz.copy(), linear, clipped, srgb,
                          linear < 0, linear > 1, exposure)


def spectrum_to_srgb(wavelength, spectral_radiance, *, exposure,
                     wavelength_unit="nm"):
    """Convenience wrapper preserving physical XYZ and clipping diagnostics."""
    xyz = spectrum_to_xyz(wavelength, spectral_radiance,
                          wavelength_unit=wavelength_unit)
    return xyz_to_srgb(xyz, exposure=exposure)
