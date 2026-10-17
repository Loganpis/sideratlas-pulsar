"""Auditable atmosphere inputs and surface redshift factors.

This module contains no atmosphere solver, blackbody substitute, angular
beaming prescription or optical extrapolation of an X-ray spectrum. NSMAXG
FITS inspection establishes what wavelengths and axes an actual file contains;
it does not turn integrated flux into local specific intensity.

An accepted external atmosphere table must provide total polarized-mode-summed
local specific intensity, not flux. Its source file is hashed independently.
Structural validation cannot authenticate an author's scientific claims: a
matching digest identifies bytes, not the physical adequacy of a model.
"""
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
import re

import numpy as np
from astropy.constants import G, M_sun, c, h
from astropy.io import fits


NSMAXG_URL = "https://heasarc.gsfc.nasa.gov/FTP/software/xspec/spectral/modelData/nsmaxg.fits"
NSMAXG_SHA256 = "848ec0c597e0a64627787e7783ba797cf43f52abc01fd406a017e1e2b89dffff"
NSMAXG_BYTES = 1301760
INTENSITY_UNIT = "W m-2 sr-1 nm-1"
VISIBLE_WAVELENGTH_NM = (400.0, 700.0)
_KEV_J = 1.602176634e-16


def _positive_finite(value, name):
    value = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(value)) or np.any(value <= 0):
        raise ValueError(f"{name} must be finite and positive")
    return value


def photon_energy_kev(wavelength_nm):
    """Vacuum photon energy hc/lambda, with explicit nm and keV units."""
    wavelength = _positive_finite(wavelength_nm, "Wavelength")
    return h.value*c.value/(wavelength*1e-9)/_KEV_J


@dataclass(frozen=True)
class SurfaceFactors:
    """Static spherical exterior factors; no rotating-star ray tracing implied."""

    mass_solar: float
    radius_km: float
    compactness_2gm_rc2: float
    frequency_ratio_g: float
    redshift_z: float
    apparent_radius_km: float
    proper_surface_gravity_m_s2: float

    def as_dict(self):
        return asdict(self)


def schwarzschild_surface_factors(mass_solar=1.4, radius_km=13.0):
    """Calculate g=sqrt(1-2GM/Rc²), R_inf=R/g and local gravity GM/(R²g).

    M and R are model choices, not measured Vela mass and radius. This static
    spherical approximation omits frame dragging, oblateness and surface
    Doppler shifts. It supplies no light-bending map or surface intensity.
    The apparent radius here is the limiting image impact parameter R/g
    for a surface strictly outside the photon sphere. Ultracompact stars
    at or within R=3GM/c² need a different silhouette treatment and are
    deliberately rejected by this API.
    """
    mass = _positive_finite(mass_solar, "Mass")
    radius = _positive_finite(radius_km, "Radius")
    if mass.ndim or radius.ndim:
        raise ValueError("Surface factors require scalar mass and radius")
    mass_kg = float(mass)*M_sun.value
    radius_m = float(radius)*1e3
    compactness = 2*G.value*mass_kg/(radius_m*c.value**2)
    if compactness >= 1:
        raise ValueError("Surface must lie outside the Schwarzschild radius")
    if compactness >= 2/3:
        raise ValueError("Apparent-radius API requires a surface outside the photon sphere (R > 3GM/c²)")
    g = float(np.sqrt(1-compactness))
    return SurfaceFactors(float(mass), float(radius), compactness, g, 1/g-1,
                          float(radius)/g, G.value*mass_kg/(radius_m**2*g))


def _shift_intensity(intensity, frequency_ratio_g, power):
    intensity = np.asarray(intensity, dtype=float)
    g = _positive_finite(frequency_ratio_g, "Frequency ratio")
    if np.any(~np.isfinite(intensity)) or np.any(intensity < 0):
        raise ValueError("Specific intensity must be finite and nonnegative")
    return intensity*g**power


def transform_specific_intensity_frequency(intensity_emitted, frequency_ratio_g):
    """I_nu,obs(nu_obs)=g³ I_nu,em(nu_obs/g), at corresponding frequencies.

    This changes the intensity value; the caller must also transform the
    spectral coordinate. It neither resamples a table nor extrapolates it.
    """
    return _shift_intensity(intensity_emitted, frequency_ratio_g, 3)


def transform_specific_intensity_wavelength(intensity_emitted, frequency_ratio_g):
    """I_lambda,obs(lambda_obs)=g⁵ I_lambda,em(g lambda_obs).

    The wavelength-density Jacobian gives g⁵, unlike the g³ frequency law.
    Input and output must use the same wavelength-density units.
    """
    return _shift_intensity(intensity_emitted, frequency_ratio_g, 5)


def transform_bolometric_intensity(intensity_emitted, frequency_ratio_g):
    """Frequency/wavelength-integrated specific intensity transforms as g⁴."""
    return _shift_intensity(intensity_emitted, frequency_ratio_g, 4)


def _file_sha256(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _axis(values, name, *, minimum=None, maximum=None):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or np.any(~np.isfinite(values)):
        raise ValueError(f"{name} must be a finite one-dimensional axis with at least two points")
    if np.any(np.diff(values) <= 0):
        raise ValueError(f"{name} must increase strictly")
    if minimum is not None and np.any(values < minimum):
        raise ValueError(f"{name} lies below its physical range")
    if maximum is not None and np.any(values > maximum):
        raise ValueError(f"{name} lies above its physical range")
    return values


def _text(value):
    # Astropy's variable-length PA column is a chararray in current files.
    if isinstance(value, (str, bytes)):
        return value.decode() if isinstance(value, bytes) else value
    return "".join(_text(item) for item in np.asarray(value).ravel())


def audit_nsmaxg_fits(path, *, expected_sha256=None,
                      wavelength_range_nm=VISIBLE_WAVELENGTH_NM):
    """Inspect actual variable-length MODELINFO energy arrays and spectral axes.

    Energy units follow the documented NSMAXG keV format. A file containing
    optical energies is not automatically a local atmosphere-intensity table.
    No network access or download occurs here.
    """
    band = _axis(wavelength_range_nm, "Requested wavelength", minimum=np.finfo(float).tiny)
    if len(band) != 2:
        raise ValueError("Requested wavelength range needs two endpoints")
    needed_energy = np.sort(photon_energy_kev(band))
    path = Path(path)
    digest = _file_sha256(path)
    if expected_sha256 is not None and digest != expected_sha256.lower():
        raise ValueError("NSMAXG source SHA256 does not match expected bytes")
    models = []
    with fits.open(path, memmap=False) as hdus:
        info = hdus["MODELINFO"]
        required = {"MODNUM", "MODNAME", "TEMPS", "GRAVS", "ENERGIES"}
        if not required.issubset(set(info.columns.names)):
            raise ValueError("NSMAXG MODELINFO is missing required columns")
        for row in info.data:
            energies = _axis(row["ENERGIES"], "NSMAXG energy", minimum=np.finfo(float).tiny)
            name = _text(row["MODNAME"]).strip()
            if name not in hdus:
                raise ValueError(f"NSMAXG spectral HDU missing: {name}")
            columns = list(hdus[name].columns.names)
            if "FLUX" not in columns:
                raise ValueError(f"NSMAXG FLUX column missing: {name}")
            covered = bool(energies[0] <= needed_energy[0]
                           and energies[-1] >= needed_energy[-1])
            models.append({
                "model_number": int(row["MODNUM"]), "name": name,
                "minimum_energy_kev": float(energies[0]),
                "maximum_energy_kev": float(energies[-1]),
                "energy_samples": len(energies), "spectral_columns": columns,
                "requested_wavelengths_covered": covered,
            })
    if not models:
        raise ValueError("NSMAXG MODELINFO contains no models")
    return {
        "source_sha256": digest, "source_bytes": path.stat().st_size,
        "expected_hash_verified": expected_sha256 is not None,
        "energy_unit": "keV", "requested_wavelength_range_nm": band.tolist(),
        "requested_energy_range_kev": needed_energy.tolist(),
        "any_model_covers_requested_wavelengths": any(
            model["requested_wavelengths_covered"] for model in models),
        "local_angledependent_specific_intensity_available": False,
        "quantity": "NSMAXG spectral flux, not local specific intensity",
        "models": models,
    }


def validate_angledependent_atmosphere_table(
        table, *, source_path, allow_test_fixture=False,
        wavelength_range_nm=VISIBLE_WAVELENGTH_NM, frequency_ratio_g=1.0):
    """Fail closed on missing optical coverage, angles, units or provenance.

    Required arrays are ``wavelength_nm``, ``cos_emission_angle`` and
    ``specific_intensity_w_m2_sr_nm``. Axisymmetric, normal-field tables use
    shape (n_mu,n_lambda). Inclined fields require ``azimuth_deg`` spanning
    0..360 and shape (n_phi,n_mu,n_lambda), including equal periodic endplanes.
    All emission angles mu=0..1 must be covered; no extrapolation is allowed.
    ``wavelength_range_nm`` is the requested observer band. For a supplied
    photon frequency ratio g=nu_obs/nu_em, the required local wavelengths
    are g times that observer band, not the observer wavelengths themselves.

    Metadata must state source/citation/hash, model type, surface T/B/gravity,
    magnetic inclination and the intensity units and summed polarization.
    ``synthetic_test_fixture`` is accepted only with the explicit test flag;
    it is never eligible as a scientific input. Structural eligibility only
    admits a file for independent review; it is never scientific validation
    or a claim of a Vela fit.
    """
    if not isinstance(table, dict) or not isinstance(table.get("metadata"), dict):
        raise ValueError("Atmosphere table requires arrays and metadata")
    meta = table["metadata"]
    required = {"source_url", "citation", "source_sha256", "source_kind",
                "intensity_unit", "polarization", "temperature_k",
                "magnetic_field_t", "magnetic_inclination_deg",
                "surface_gravity_m_s2", "azimuthally_symmetric", "model_assumptions"}
    if not required.issubset(meta):
        raise ValueError("Atmosphere metadata missing: "+", ".join(sorted(required-set(meta))))
    for name in ("source_url", "citation"):
        if not isinstance(meta[name], str) or not meta[name].strip():
            raise ValueError(f"Atmosphere metadata needs {name}")
    if not meta["source_url"].startswith("https://"):
        raise ValueError("Atmosphere source URL must identify an HTTPS source")
    if not isinstance(meta["source_sha256"], str) or not re.fullmatch(
            r"[0-9a-f]{64}", meta["source_sha256"]):
        raise ValueError("Atmosphere source SHA256 must be a lowercase 64-digit digest")
    kind = meta["source_kind"]
    fixture = kind == "synthetic_test_fixture"
    if kind not in {"published_model", "author_provided_model", "synthetic_test_fixture"}:
        raise ValueError("Unsupported atmosphere source kind")
    if fixture and not allow_test_fixture:
        raise ValueError("Synthetic atmosphere fixtures cannot be used for production")
    if _file_sha256(source_path) != meta["source_sha256"]:
        raise ValueError("Atmosphere source SHA256 does not match the supplied file")
    if meta["intensity_unit"] != INTENSITY_UNIT or meta["polarization"] != "total":
        raise ValueError("Provide total specific intensity in W m-2 sr-1 nm-1, not flux or a single mode")
    for name in ("temperature_k", "magnetic_field_t", "surface_gravity_m_s2"):
        if _positive_finite(meta[name], name).ndim:
            raise ValueError(f"Atmosphere {name} must be scalar for this local table")
    inclination = float(meta["magnetic_inclination_deg"])
    if not np.isfinite(inclination) or not 0 <= inclination <= 180:
        raise ValueError("Magnetic inclination must be within 0..180 degrees")
    if type(meta["azimuthally_symmetric"]) is not bool:
        raise ValueError("azimuthally_symmetric must be an explicit Boolean")
    if (not isinstance(meta["model_assumptions"], list)
            or not meta["model_assumptions"]
            or not all(isinstance(v, str) and v.strip() for v in meta["model_assumptions"])):
        raise ValueError("Atmosphere model assumptions must be explicitly listed")
    wave = _axis(table["wavelength_nm"], "Wavelength", minimum=np.finfo(float).tiny)
    mu = _axis(table["cos_emission_angle"], "Emission cosine", minimum=0, maximum=1)
    requested = _axis(wavelength_range_nm, "Requested wavelength", minimum=np.finfo(float).tiny)
    g = _positive_finite(frequency_ratio_g, "Frequency ratio")
    if g.ndim:
        raise ValueError("Atmosphere wavelength coverage requires a scalar frequency ratio")
    local_required = requested*float(g)
    if len(requested) != 2 or wave[0] > local_required[0] or wave[-1] < local_required[-1]:
        raise ValueError("Atmosphere table does not cover the local wavelengths required by the observer band and frequency ratio; extrapolation forbidden")
    if mu[0] != 0 or mu[-1] != 1:
        raise ValueError("Atmosphere emission-angle grid must cover mu=0..1 without extrapolation")
    intensity = np.asarray(table["specific_intensity_w_m2_sr_nm"], dtype=float)
    if (np.any(~np.isfinite(intensity)) or np.any(intensity < 0)
            or not np.any(intensity > 0)):
        raise ValueError("Atmosphere intensities must be finite, nonnegative and nonzero")
    if meta["azimuthally_symmetric"]:
        if inclination not in (0, 180):
            raise ValueError("Inclined magnetic fields need an azimuth-resolved intensity table")
        expected_shape = (len(mu), len(wave))
    else:
        phi = _axis(table.get("azimuth_deg"), "Azimuth", minimum=0, maximum=360)
        if phi[0] != 0 or phi[-1] != 360:
            raise ValueError("Azimuth must cover the periodic 0..360 degree interval")
        expected_shape = (len(phi), len(mu), len(wave))
        if intensity.shape == expected_shape and not np.allclose(
                intensity[0], intensity[-1], rtol=1e-10, atol=0):
            raise ValueError("Atmosphere azimuth endplanes must be periodic")
    if intensity.shape != expected_shape:
        raise ValueError(f"Atmosphere intensity shape must be {expected_shape}")
    return {
        "structurally_valid": True, "structurally_eligible_for_review": not fixture,
        "scientific_validation": False,
        "source_hash_verified": True, "source_kind": kind,
        "requested_wavelengths_covered": True,
        "observer_wavelength_range_nm": requested.tolist(),
        "required_local_wavelength_range_nm": local_required.tolist(),
        "frequency_ratio_g": float(g),
        "wavelength_range_nm": [float(wave[0]), float(wave[-1])],
        "shape": list(intensity.shape), "extrapolation": False,
        "scientific_adequacy": "Requires independent model and Vela validation",
    }
