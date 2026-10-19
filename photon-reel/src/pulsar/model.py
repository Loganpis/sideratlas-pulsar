"""Oblique split-monopole fields and optically thin, prescribed sheet emission.

Units: c = Omega = R_LC = 1; B_LC is an arbitrary field normalization.
The discontinuous field is valid on either face of the current sheet. Its
internal dissipative structure and emissivity are NOT solved here.
"""
from dataclasses import dataclass, asdict
import numpy as np
from numpy.polynomial.legendre import leggauss

TAU = 2 * np.pi
C_KM_S = 299792.458


@dataclass(frozen=True)
class Parameters:
    period_s: float = 0.089
    obliquity_deg: float = 60.0
    observer_deg: float = 63.6
    r_min: float = 1.5
    r_max: float = 5.0
    lorentz_factor: float = 10.0
    sheet_width: float = 0.05
    radial_index: float = 3.0
    spectral_index: float = 1.0

    def __post_init__(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError("All parameters must be finite")
        if self.period_s <= 0 or self.r_min < 1 or self.r_max <= self.r_min:
            raise ValueError("Require P>0 and 1<=r_min<r_max for this wind model")
        if self.lorentz_factor ** 2 <= 1 + self.r_max ** 2:
            raise ValueError("Gamma must exceed the maximum force-free drift Gamma")
        if not 0 < self.sheet_width < 1:
            raise ValueError("Sheet width must lie in (0,1)")
        if not 0 <= self.obliquity_deg <= 180 or not 0 <= self.observer_deg <= 180:
            raise ValueError("Angles must lie in [0,180] degrees")

    @property
    def chi(self):
        return np.deg2rad(self.obliquity_deg)

    @property
    def zeta(self):
        return np.deg2rad(self.observer_deg)

    @property
    def r_lc_km(self):
        return C_KM_S * self.period_s / TAU


def observer_direction(zeta):
    return np.array([np.sin(zeta), 0., np.cos(zeta)])


def sheet_function(xyz, t, chi):
    xyz = np.asarray(xyz, dtype=float)
    r = np.linalg.norm(xyz, axis=-1)
    if np.any(r == 0):
        raise ValueError("The wind solution excludes the origin")
    er = xyz / r[..., None]
    phase = t - r
    return (er[..., 2] * np.cos(chi) + np.sin(chi) *
            (er[..., 0] * np.cos(phase) + er[..., 1] * np.sin(phase)))


def fields(xyz, t, chi, polarity=None):
    """Return E and B in Cartesian components; E and B share Gaussian units.

    Optional polarity is for taking the limiting field on a sheet face.
    At u=0 the returned positive face is a convention, not a resolved sheet.
    """
    xyz = np.asarray(xyz, dtype=float)
    r = np.linalg.norm(xyz, axis=-1)
    if np.any(r == 0):
        raise ValueError("The wind solution excludes the origin")
    x, y, z = np.moveaxis(xyz, -1, 0)
    sign = np.where(sheet_function(xyz, t, chi) >= 0, 1., -1.) if polarity is None else polarity
    b = np.stack((x / r**3 + y / r**2,
                  y / r**3 - x / r**2, z / r**3), axis=-1)
    e = np.stack((-z*x / r**3, -z*y / r**3, (x*x+y*y) / r**3), axis=-1)
    return e * np.asarray(sign)[..., None], b * np.asarray(sign)[..., None]


def charge_current_4pi(xyz, t, chi):
    """4*pi*rho and 4*pi*J, excluding distributional sheet sources."""
    xyz = np.asarray(xyz, dtype=float)
    r = np.linalg.norm(xyz, axis=-1)
    sign = np.where(sheet_function(xyz, t, chi) >= 0, 1., -1.)
    rho = -2 * sign * xyz[..., 2] / r**3
    return rho, rho[..., None] * xyz / r[..., None]


def emitting_velocity(xyz, gamma):
    """Prescribed outward streaming on the two limiting force-free faces.

    v = E cross B / B^2 + parallel component, with |v|^2 = 1-Gamma^-2.
    This is an emission prescription: sheet acceleration is not calculated.
    """
    xyz = np.asarray(xyz, dtype=float)
    r = np.linalg.norm(xyz, axis=-1)
    er = xyz / r[..., None]
    cylindrical_r = np.hypot(xyz[..., 0], xyz[..., 1])
    q = cylindrical_r
    discriminant = 1 - (1 + q*q) / gamma**2
    if np.any(discriminant <= 0):
        raise ValueError("Chosen Gamma cannot support the local field drift")
    a = np.sqrt(discriminant)
    vr = (q*q + a) / (1 + q*q)
    vp = q * (1 - a) / (1 + q*q)
    ephi = np.zeros_like(xyz)
    np.divide(-xyz[..., 1], cylindrical_r, out=ephi[..., 0], where=cylindrical_r != 0)
    np.divide(xyz[..., 0], cylindrical_r, out=ephi[..., 1], where=cylindrical_r != 0)
    return vr[..., None] * er + vp[..., None] * ephi


def analytic_separation(chi, zeta):
    """Shorter separation (cycles), only for pencil-beam sheet crossings."""
    denominator = np.sin(chi) * np.sin(zeta)
    if abs(denominator) < 1e-14:
        return None
    a = abs(np.cos(chi) * np.cos(zeta) / denominator)
    return float(np.arccos(a) / np.pi) if a < 1 else None


def integration_grid(p, n_r=16, n_mu=64, n_az=128):
    """Full solid-angle quadrature about the observer, not a clipped beam cone."""
    if min(n_r, n_mu, n_az) < 2:
        raise ValueError("Quadrature dimensions must be at least two")
    xr, wr = leggauss(n_r)
    r = p.r_min + (xr + 1) * (p.r_max-p.r_min)/2
    wr *= (p.r_max-p.r_min)/2
    mu, wm = leggauss(n_mu)
    az = TAU * (np.arange(n_az) + 0.5) / n_az
    n = observer_direction(p.zeta)
    e1 = np.array([np.cos(p.zeta), 0., -np.sin(p.zeta)])
    e2 = np.array([0., 1., 0.])
    er = (mu[:, None, None] * n + np.sqrt(1-mu**2)[:, None, None] *
          (np.cos(az)[None, :, None]*e1 + np.sin(az)[None, :, None]*e2))
    xyz = r[:, None, None, None] * er[None, ...]
    weights = wr[:, None, None] * r[:, None, None]**2 * wm[None, :, None] * np.full((1, 1, n_az), TAU/n_az)
    return xyz.reshape(-1, 3), weights.ravel()


def pulse_profile(p, phases, resolution=(16, 64, 128), delays=True):
    """Integrate j_nu at emission time t_obs+n.r over the lab volume.

    Assumed comoving j'_nu' proportional to r^-q exp[-u^2/(2w^2)] nu'^-alpha.
    j_nu/nu^2 invariance gives D^(2+alpha), NOT the moving-blob D^(3+alpha).
    No attenuation, stellar occultation, GR, spectral cutoff or absolute flux.
    """
    xyz, weights = integration_grid(p, *resolution)
    n = observer_direction(p.zeta)
    r = np.linalg.norm(xyz, axis=-1)
    v = emitting_velocity(xyz, p.lorentz_factor)
    doppler = 1/(p.lorentz_factor * (1 - v @ n))
    weights = weights * r**(-p.radial_index) * doppler**(2+p.spectral_index)
    er = xyz/r[:, None]
    offset = r - (xyz @ n if delays else 0.)
    # u(t_obs) = constant + cos-coefficient cos(t_obs) + sin-coefficient sin(t_obs).
    c = er[:, 2] * np.cos(p.chi)
    a = np.sin(p.chi) * (er[:, 0]*np.cos(offset) - er[:, 1]*np.sin(offset))
    b = np.sin(p.chi) * (er[:, 0]*np.sin(offset) + er[:, 1]*np.cos(offset))
    result = np.empty(len(phases))
    for i, phase in enumerate(phases):
        t = TAU*phase
        u = c + a*np.cos(t) + b*np.sin(t)
        result[i] = np.sum(weights * np.exp(-0.5*(u/p.sheet_width)**2))
    return result


def normalize_profile(profile):
    mean = np.mean(profile)
    if not np.isfinite(mean) or mean <= 0:
        raise ValueError("Profile has no finite positive mean")
    return np.asarray(profile) / mean


def peak_separation(phases, flux):
    """Two strongest local maxima; returned in phase order."""
    indices = np.flatnonzero((flux > np.roll(flux, 1)) & (flux >= np.roll(flux, -1)))
    if len(indices) < 2:
        return None, []
    indices = indices[np.argsort(flux[indices])[-2:]]
    peaks = sorted(float(phases[i] % 1) for i in indices)
    delta = peaks[1]-peaks[0]
    return float(min(delta, 1-delta)), peaks
