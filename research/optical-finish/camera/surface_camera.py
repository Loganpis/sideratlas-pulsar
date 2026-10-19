"""Schwarzschild surface rays and a spectral camera for a spherical star.

The atmosphere is an external specific-intensity function, never an invented
glowing sphere. This module calculates outward direct rays, local SR rotation,
redshifts and relative travel times. It does not solve the atmosphere, stellar
structure, frame dragging, oblateness or magnetospheric radiative transfer.

References: Poutanen & Beloborodov (2006), doi:10.1111/j.1365-2966.2006.11088.x;
Poutanen (2020), doi:10.1051/0004-6361/202037471. Integrals are evaluated directly,
without the approximate cos(alpha) bending formula.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.special import roots_legendre

C_M_S = 299792458.0
G_SI = 6.67430e-11
M_SUN_KG = 1.988409870698051e30


def _positive(value, label, allow_zero=False):
    if not np.isfinite(value) or value < 0 or (not allow_zero and value == 0):
        raise ValueError(f'{label} must be finite and positive')
    return float(value)


def unit_vector(vector):
    a = np.asarray(vector, float)
    if a.shape != (3,) or not np.isfinite(a).all() or np.linalg.norm(a) == 0:
        raise ValueError('Expected one nonzero finite 3-vector')
    return a/np.linalg.norm(a)


@dataclass(frozen=True)
class Star:
    mass_solar: float = 1.4
    radius_m: float = 13000.

    def __post_init__(self):
        _positive(self.mass_solar, 'Mass', allow_zero=True)
        _positive(self.radius_m, 'Radius')
        if self.radius_m <= 1.5*self.schwarzschild_radius_m:
            raise ValueError('This outward-ray camera requires a surface above the photon sphere')

    @property
    def schwarzschild_radius_m(self):
        return 2*G_SI*self.mass_solar*M_SUN_KG/C_M_S**2

    @property
    def compactness(self):
        return self.schwarzschild_radius_m/self.radius_m

    @property
    def frequency_ratio_g(self):
        return np.sqrt(1-self.compactness)

    @property
    def apparent_radius_m(self):
        return self.radius_m/self.frequency_ratio_g


@lru_cache(maxsize=8)
def _quadrature(order):
    if not isinstance(order, int) or order < 16:
        raise ValueError('Quadrature order must be an integer >=16')
    z, w = roots_legendre(order)
    return (z+1)/2, w/2


def outward_geodesics(impact_m, star=Star(), *, order=192):
    """Return psi, cos(alpha), and travel delay relative to the radial ray.

    x=R/r and x=1-y² remove the integrable tangential-ray endpoint. With
    k=b/R, A=1-k²x²(1-u*x). The rationalized delay integrand avoids subtracting
    nearly equal values at infinity. No numerical geodesic time divergence is
    introduced by choosing a large but finite observer radius.
    """
    b = np.asarray(impact_m, float)
    if np.any(~np.isfinite(b)) or np.any(b < 0) or np.any(b > star.apparent_radius_m*(1+1e-13)):
        raise ValueError('Impact parameters must hit the outward-emitting surface')
    b = np.minimum(b, star.apparent_radius_m)
    k = (b/star.radius_m)[..., None]
    mu2 = np.maximum(0., 1-(b/star.apparent_radius_m)**2)[..., None]
    u = star.compactness
    y, weight = _quadrature(order)
    x = 1-y*y
    h = (1-u)*(1+x)-u*x*x
    root_a = np.sqrt(mu2+k*k*y*y*h)
    psi = np.sum(weight*2*y*k/root_a, axis=-1)
    delay = star.radius_m/C_M_S*np.sum(weight*2*y*k*k/(root_a*(1+root_a)), axis=-1)
    return psi, np.sqrt(mu2[..., 0]), delay


def camera_basis(observer_direction, up_hint=(0, 0, 1)):
    n = unit_vector(observer_direction)
    hint = unit_vector(up_hint)
    ey = hint-np.dot(hint, n)*n
    if np.linalg.norm(ey) < 1e-10:
        hint = np.eye(3)[np.argmin(np.abs(n))]
        ey = hint-np.dot(hint, n)*n
    ey /= np.linalg.norm(ey)
    ex = np.cross(ey, n)
    return ex, ey, n


@dataclass
class SurfaceRayMap:
    star: Star
    hit: np.ndarray
    normal_static: np.ndarray
    photon_static: np.ndarray
    cos_alpha: np.ndarray
    relative_delay_s: np.ndarray
    arrival_offset_s: np.ndarray
    observer_direction: np.ndarray
    reference_radius_m: float


def surface_raymap(x_m, y_m, *, star=Star(), observer_direction=(0, 0, 1),
                   up_hint=(0, 1, 0), reference_radius_m=None, order=192):
    """Map physical image-plane coordinates to surface rays at infinity.

    Arrival offsets subtract a common observer-distance/Shapiro constant:
    offset= -R/c -rs/c*log((R-rs)/reference_radius) + relative_delay.
    The reference sets a common time origin, not a measurable absolute delay.
    Other emitting layers must use the SAME propagation convention before
    their phases can be compared. The flat-space limit is -n.r/c.
    """
    x, y = np.broadcast_arrays(np.asarray(x_m, float), np.asarray(y_m, float))
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite image-plane coordinates')
    reference = star.radius_m if reference_radius_m is None else _positive(reference_radius_m, 'Time reference radius')
    ex, ey, n = camera_basis(observer_direction, up_hint)
    b = np.hypot(x, y)
    hit = b <= star.apparent_radius_m
    normal = np.zeros(x.shape+(3,)); photon = normal.copy()
    mu = np.zeros_like(x); relative = mu.copy(); offset = mu.copy()
    if hit.any():
        psi, mu_hit, delay = outward_geodesics(b[hit], star, order=order)
        bhat = (x[hit, None]*ex+y[hit, None]*ey)/np.maximum(b[hit, None], np.finfo(float).tiny)
        alpha = np.arccos(mu_hit)
        normal[hit] = np.cos(psi)[:, None]*n+np.sin(psi)[:, None]*bhat
        photon[hit] = np.cos(psi-alpha)[:, None]*n+np.sin(psi-alpha)[:, None]*bhat
        mu[hit], relative[hit] = mu_hit, delay
        offset[hit] = (-star.radius_m/C_M_S
                      -star.schwarzschild_radius_m/C_M_S*np.log((star.radius_m-star.schwarzschild_radius_m)/reference)
                      +delay)
    return SurfaceRayMap(star, hit, normal, photon, mu, relative, offset, n, reference)


def rotate_vectors(vectors, axis, angle):
    """Active Rodrigues rotation, broadcasting one angle per vector."""
    v = np.asarray(vectors, float)
    axis = unit_vector(axis)
    a = np.asarray(angle, float)[..., None]
    return v*np.cos(a)+np.cross(axis, v)*np.sin(a)+np.sum(v*axis, axis=-1, keepdims=True)*axis*(1-np.cos(a))


def aberrate_to_comoving(photon_static, beta):
    """Static orthonormal tetrad -> co-moving SR tetrad, with Doppler factor."""
    k, beta = np.broadcast_arrays(np.asarray(photon_static, float), np.asarray(beta, float))
    if k.shape[-1:] != (3,) or np.any(~np.isfinite(k)) or np.any(~np.isfinite(beta)):
        raise ValueError('Invalid photon/velocity vectors')
    if not np.allclose(np.sum(k*k, axis=-1), 1, rtol=0, atol=1e-10):
        raise ValueError('Photon directions must be unit vectors')
    beta2 = np.sum(beta*beta, axis=-1)
    if np.any(beta2 >= 1):
        raise ValueError('Surface velocities must be subluminal')
    gamma = 1/np.sqrt(1-beta2)
    dot = np.sum(beta*k, axis=-1)
    doppler = 1/(gamma*(1-dot))
    # (gamma-1)/beta² = gamma²/(gamma+1), including the zero-velocity limit.
    factor = gamma*gamma/(gamma+1)*dot-gamma
    comoving = (k+factor[..., None]*beta)*doppler[..., None]
    return comoving, doppler


def spectral_surface_image(raymap, observer_time_s, wavelength_observed_nm,
                           atmosphere_intensity, *, omega_rad_s=0., spin_axis=(0, 0, 1), phase0=0.):
    """Evaluate physical surface I_lambda on existing image-plane samples.

    Callback keyword contract: wavelength_emitted_nm [Nhits,Nlambda],
    cos_emission_angle [Nhits], normal_body/direction_body [Nhits,3], and
    emission_time_s [Nhits]. Return I_lambda in W m^-2 sr^-1 nm^-1 in the
    surface co-moving frame. This does not validate the callback's physics.

    Rotation uses the Schwarzschild-plus-Doppler spherical-star approximation;
    beta=(Omega cross r)/(c*g_gravity). Frame dragging and oblateness are absent.
    """
    wavelength = np.asarray(wavelength_observed_nm, float)
    if (wavelength.ndim != 1 or wavelength.size < 2 or np.any(~np.isfinite(wavelength))
            or np.any(wavelength <= 0) or np.any(np.diff(wavelength) <= 0)):
        raise ValueError('Wavelength axis must be positive, finite and increasing')
    if not np.isfinite(observer_time_s) or not np.isfinite(omega_rad_s) or not np.isfinite(phase0):
        raise ValueError('Invalid observation time or spin')
    axis = unit_vector(spin_axis)
    output = np.zeros(raymap.hit.shape+(len(wavelength),))
    if not raymap.hit.any():
        return output
    normal = raymap.normal_static[raymap.hit]
    photon = raymap.photon_static[raymap.hit]
    beta = np.cross(omega_rad_s*axis, raymap.star.radius_m*normal)/(C_M_S*raymap.star.frequency_ratio_g)
    k_em, doppler = aberrate_to_comoving(photon, beta)
    cos_emission = np.sum(normal*k_em, axis=-1)
    if np.any(cos_emission < -1e-10) or np.any(cos_emission > 1+1e-10):
        raise ValueError('Invalid aberrated outward emission direction')
    time_emission = observer_time_s-raymap.arrival_offset_s[raymap.hit]
    phase = phase0+omega_rad_s*time_emission
    g = raymap.star.frequency_ratio_g*doppler
    radiance = np.asarray(atmosphere_intensity(
        wavelength_emitted_nm=g[:, None]*wavelength,
        cos_emission_angle=np.clip(cos_emission, 0, 1),
        normal_body=rotate_vectors(normal, axis, -phase),
        direction_body=rotate_vectors(k_em, axis, -phase),
        emission_time_s=time_emission), float)
    if radiance.shape != (int(raymap.hit.sum()), len(wavelength)):
        raise ValueError('Atmosphere returned the wrong radiance shape')
    if np.any(~np.isfinite(radiance)) or np.any(radiance < 0):
        raise ValueError('Atmosphere returned negative/nonfinite radiance')
    output[raymap.hit] = g[:, None]**5*radiance
    return output


def integrate_surface_flux(spectral_image, sample_area_m2, distance_m):
    """Far-observer F_lambda from I_lambda and physical image-plane areas."""
    distance = _positive(distance_m, 'Observer distance')
    image = np.asarray(spectral_image, float)
    area = np.broadcast_to(np.asarray(sample_area_m2, float), image.shape[:-1])
    if np.any(~np.isfinite(image)) or np.any(image < 0) or np.any(~np.isfinite(area)) or np.any(area < 0):
        raise ValueError('Invalid radiance image or sample areas')
    return np.sum((image*area[..., None]).reshape(-1, image.shape[-1]), axis=0)/distance**2
