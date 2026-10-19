"""Laboratory-frame curvature radiation and reduced particle-energy transport.

The trajectories and accelerating work are inputs, not a solved plasma model.
Following the curvature-loss energy equation used by Barnard et al. (2022,
https://arxiv.org/abs/2111.03405), the independent variable is lab arclength:

    d gamma / ds = e E_work / (m_e c^2) - 2 e^2 gamma^4 / (3 m_e c^2 rho^2).

Here ``E_work`` is the electric work per unit positive charge and lab
arclength, E dot v / |v| with the charge sign included. It equals E_parallel
only for motion parallel to that electric field. Callers must project their
field onto their *actual particle velocity*. Curvature must likewise come
from the lab trajectory, not merely from a magnetic field line.

All emission is calculated in the lab frame. Applying another bulk Doppler
factor to these spectra would double-transform the emission. The equations
assume ultrarelativistic particles; this module is not a PIC calculation.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.special import gamma as gamma_function, kv


# Gaussian CGS; SI-to-CGS conversions of CODATA constants.
C_CM_S = 29979245800.0
E_STATC = 4.803204712570263e-10
M_E_G = 9.1093837139e-28
HBAR_ERG_S = 1.0545718176461565e-27
MEC2_ERG = M_E_G * C_CM_S**2
GEV_ERG = 1.602176634e-3
KERNEL_POWER_INTEGRAL = 8 * np.pi / (9 * np.sqrt(3))
KERNEL_PHOTON_INTEGRAL = 5 * np.pi / 3
_SMALL_X_COEFFICIENT = 2**(2/3) * gamma_function(2/3)
_SMALL_X_LINEAR = np.pi / np.sqrt(3)


@lru_cache(maxsize=1)
def _spectral_table():
    """Log-grid integration retaining the exponentially faint high-E tail.

    Integrate from high to low energy, avoiding subtraction of two nearly
    equal total integrals in the high-energy tail. 24,001 points give better
    than 1e-4 relative kernel accuracy through x=40 (tested independently).
    """
    logx = np.linspace(np.log(1e-10), np.log(700.0), 24001)
    x = np.exp(logx)

    def tail(values):
        return cumulative_simpson(values[::-1], x=-logx[::-1], initial=0)[::-1]

    f = x * tail(kv(5/3, x) * x)
    # scipy's Bessel implementation eventually underflows. These photons
    # carry negligible power, but keep their logarithmic tail continuous
    # using the large-x expansion of the *integrated* Bessel kernel.
    far_tail = x > 100
    xf = x[far_tail]
    f[far_tail] = (np.sqrt(np.pi*xf/2)*np.exp(-xf)
                   * (1+55/(72*xf)-10151/(10368*xf**2)))
    photon_tail = tail(f)
    power_tail = tail(f*x)
    return logx, np.log(f), photon_tail, power_tail


def curvature_kernel(x):
    """F(x)=x integral_x^infinity K_(5/3)(u) du, for x >= 0.

    The two-term small-x expansion is used below the integration table.
    Above x=700 the exponentially negligible kernel is set to zero.
    """
    x = np.asarray(x, dtype=float)
    if np.any(np.isnan(x)) or np.any(x < 0):
        raise ValueError("Curvature-kernel energies must be nonnegative")
    lx, lf, _, _ = _spectral_table()
    result = np.zeros_like(x)
    small = (x > 0) & (x < np.exp(lx[0]))
    result[small] = (_SMALL_X_COEFFICIENT*np.cbrt(x[small])
                     - _SMALL_X_LINEAR*x[small])
    table = (x >= np.exp(lx[0])) & (x <= np.exp(lx[-1]))
    result[table] = np.exp(np.interp(np.log(x[table]), lx, lf))
    return result


def critical_energy_gev(gamma, rho_cm):
    """Curvature critical photon energy 3 hbar c gamma^3 / (2 rho)."""
    gamma, rho_cm = np.broadcast_arrays(np.asarray(gamma, dtype=float),
                                        np.asarray(rho_cm, dtype=float))
    return 1.5 * HBAR_ERG_S * C_CM_S * gamma**3 / rho_cm / GEV_ERG


def curvature_power_erg_s(gamma, rho_cm):
    """Total curvature power 2 e^2 c gamma^4 / (3 rho^2)."""
    gamma, rho_cm = np.broadcast_arrays(np.asarray(gamma, dtype=float),
                                        np.asarray(rho_cm, dtype=float))
    return 2 * E_STATC**2 * C_CM_S * gamma**4 / (3*rho_cm**2)


def _tail_integral(x, photon):
    x = np.asarray(x, dtype=float)
    lx, _, ntail, ptail = _spectral_table()
    total = KERNEL_PHOTON_INTEGRAL if photon else KERNEL_POWER_INTEGRAL
    result = np.zeros_like(x)
    small = x < np.exp(lx[0])
    xs = x[small]
    if photon:
        lower = 3*_SMALL_X_COEFFICIENT*np.cbrt(xs) - _SMALL_X_LINEAR*xs
    else:
        lower = .75*_SMALL_X_COEFFICIENT*xs**(4/3) - .5*_SMALL_X_LINEAR*xs**2
    result[small] = total-lower
    table = (x >= np.exp(lx[0])) & (x < np.exp(lx[-1]))
    # Interpolate log-tail values to preserve the exponential tail. The
    # final (zero) tail value is represented below floating-point range.
    tail = ntail if photon else ptail
    logtail = np.log(np.maximum(tail, np.finfo(float).tiny))
    result[table] = np.exp(np.interp(np.log(x[table]), lx, logtail))
    return result


def _band_integrals(critical_gev, bands_gev, photon):
    bands = np.asarray(bands_gev, dtype=float)
    if bands.ndim != 2 or bands.shape[-1] != 2:
        raise ValueError("Energy bands must have shape (number_of_bands,2)")
    if (np.any(np.isnan(bands)) or np.any(bands[:, 0] < 0)
            or np.any(bands[:, 1] <= bands[:, 0])):
        raise ValueError("Energy bands require 0 <= minimum < maximum")
    ec = np.asarray(critical_gev, dtype=float)
    if np.any(~np.isfinite(ec)) or np.any(ec < 0):
        raise ValueError("Critical energy must be finite and nonnegative")
    safe_ec = np.where(ec > 0, ec, 1)[..., None]
    lo, hi = bands[:, 0]/safe_ec, bands[:, 1]/safe_ec
    integrals = _tail_integral(lo, photon)-_tail_integral(hi, photon)
    # Avoid cancellation when both band edges are in the very-low-x tail.
    small = hi < 1e-8
    low, high = lo[small], hi[small]
    if photon:
        integrals[small] = (3*_SMALL_X_COEFFICIENT*(np.cbrt(high)-np.cbrt(low))
                            - _SMALL_X_LINEAR*(high-low))
    else:
        integrals[small] = (.75*_SMALL_X_COEFFICIENT*(high**(4/3)-low**(4/3))
                            - .5*_SMALL_X_LINEAR*(high**2-low**2))
    return np.where(ec[..., None] > 0, np.maximum(integrals, 0), 0)


def curvature_photon_yield(radiated_erg, critical_gev, bands_gev):
    """Photon numbers in energy bands for a bolometric emitted energy.

    dN/dE = (E_radiated / E_c) F(E/E_c) / (I E), where
    I=integral F(x) dx. The GeV-to-erg conversion is applied once.
    These are emitted photon numbers, not detector weighted counts.
    """
    energy, ec = np.broadcast_arrays(radiated_erg, critical_gev)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0):
        raise ValueError("Radiated energy must be finite and nonnegative")
    integral = _band_integrals(ec, bands_gev, photon=True)
    denominator = np.where(ec > 0, ec*GEV_ERG, 1)
    return energy[..., None]/denominator[..., None]*integral/KERNEL_POWER_INTEGRAL


def curvature_band_energy_fraction(critical_gev, bands_gev):
    """Fraction of total radiated energy in each photon-energy band."""
    return _band_integrals(critical_gev, bands_gev, photon=False)/KERNEL_POWER_INTEGRAL


def _advance_constant(gamma, acceleration, loss, ds):
    """Exact positive flow of gamma'=A-B gamma^4 for frozen A,B >= 0.

    With y=gamma/(A/B)^(1/4), the separable integral is
    [atan(y)+atanh(y)]/2 below equilibrium, with acoth(y) above it.
    Substitution y=tanh(z) or coth(z) gives a monotone scalar inversion.
    Bisection is stable even when the loss relaxation is much shorter
    than one supplied segment; no explicit high-gamma instability occurs.
    """
    gamma, acceleration, loss, ds = np.broadcast_arrays(gamma, acceleration, loss, ds)
    out = gamma.copy()
    acceleration_only = (loss == 0) & (ds > 0)
    out[acceleration_only] = (gamma[acceleration_only]
                             + acceleration[acceleration_only]*ds[acceleration_only])
    cooling = (acceleration == 0) & (loss > 0) & (ds > 0)
    out[cooling] = gamma[cooling]/np.cbrt(1+3*loss[cooling]*gamma[cooling]**3*ds[cooling])
    mixed = (acceleration > 0) & (loss > 0) & (ds > 0)
    if np.any(mixed):
        a, b, h, g = acceleration[mixed], loss[mixed], ds[mixed], gamma[mixed]
        equilibrium = (a/b)**.25
        y = g/equilibrium
        close = np.abs(y-1) < 1e-13
        below = y < 1
        argument = np.where(below, y, 1/y)
        z0 = np.arctanh(np.minimum(argument, 1-1e-15))
        target = z0+np.arctan(y)+2*a*h/equilibrium
        lo = np.maximum(z0, target-np.pi/2)
        hi = target.copy()
        for _ in range(48):
            mid = (lo+hi)/2
            t = np.tanh(mid)
            ym = np.where(below, t, 1/np.maximum(t, np.finfo(float).tiny))
            left = mid+np.arctan(ym) < target
            lo = np.where(left, mid, lo)
            hi = np.where(left, hi, mid)
        t = np.tanh((lo+hi)/2)
        yn = np.where(below, t, 1/np.maximum(t, np.finfo(float).tiny))
        out[mixed] = np.where(close, g, equilibrium*yn)
    return out


@dataclass(frozen=True)
class CurvatureTransport:
    """Per-particle lab-frame outputs; axis zero indexes trajectory segments."""
    gamma: np.ndarray
    gamma_mid: np.ndarray
    work_erg: np.ndarray
    radiated_erg: np.ndarray
    critical_energy_gev: np.ndarray
    photon_yield: np.ndarray | None

    @property
    def energy_balance_error_erg(self):
        return (self.work_erg-self.radiated_erg
                - MEC2_ERG*np.diff(self.gamma, axis=0))


def transport_curvature(ds_cm, rho_cm, e_parallel_statvolt_cm, gamma0,
                        energy_bands_gev=None, spectral_nodes=1):
    """Transport particles along supplied curved lab trajectories.

    ``e_parallel_statvolt_cm`` is retained as a convenient input name, but
    MUST contain the charge-aligned work projection ``E_work`` described
    at module level. It must be nonnegative in this accelerating model.
    Geometry and work are frozen on each segment; their variation requires
    a trajectory/segment convergence test even though the frozen-coefficient
    energy propagation is exact. Gamma0 must be >=1 and modeled trajectories
    must stay ultrarelativistic. Radiation, work, and particle-energy change
    are conservatively balanced on every segment.

    Segment axis is first; remaining axes index particles. One-dimensional
    ds/rho/E arrays are interpreted as segment arrays and broadcast over
    seeds. Scalars are constant on every segment. ds=0 is an inactive padded
    segment and emits nothing. rho=infinity gives straight-line zero losses.

    ``spectral_nodes`` is 1 (midpoint) or 3 (Gauss-Legendre). Spectrum nodes
    are weighted by local curvature power, then normalized to the exactly
    integrated segment energy loss. Refinement is needed for spectra during
    rapid cooling even if the integrated energy loss is already exact.
    """
    if spectral_nodes not in (1, 3):
        raise ValueError("spectral_nodes must be 1 or 3")
    arrays = [np.asarray(a, dtype=float) for a in
              (ds_cm, rho_cm, e_parallel_statvolt_cm)]
    g0 = np.asarray(gamma0, dtype=float)
    ndim = max(1+g0.ndim, *(a.ndim for a in arrays), 1)
    arrays = [a.reshape(a.shape+(1,)*(ndim-1)) if a.ndim == 1 and ndim > 1
              else a for a in arrays]
    ds, rho, work_field, _ = np.broadcast_arrays(*arrays, g0[None, ...])
    if ds.ndim == 0:
        ds, rho, work_field = (a[None] for a in (ds, rho, work_field))
    if (np.any(~np.isfinite(ds)) or np.any(ds < 0)
            or np.any(np.isnan(rho)) or np.any(rho <= 0)
            or np.any(~np.isfinite(work_field)) or np.any(work_field < 0)
            or np.any(~np.isfinite(g0)) or np.any(g0 < 1)):
        raise ValueError("Require ds>=0, rho>0, E_work>=0 and gamma0>=1")
    acceleration = E_STATC*work_field/MEC2_ERG
    loss = 2*E_STATC**2/(3*MEC2_ERG*rho**2)
    shape = ds.shape
    gamma = np.empty((shape[0]+1,)+shape[1:])
    gamma[0] = g0
    gamma_mid = np.empty(shape)
    work_erg = acceleration*ds*MEC2_ERG
    radiated = np.empty(shape)
    ec = np.empty(shape)
    photons = None if energy_bands_gev is None else np.empty(
        shape+(len(energy_bands_gev),))
    if spectral_nodes == 3:
        node_locations = .5*(1+np.array([-np.sqrt(3/5), 0., np.sqrt(3/5)]))
        node_weights = np.array([5/18, 4/9, 5/18])
    else:
        node_locations, node_weights = [.5], [1.]
    for i in range(shape[0]):
        a, b, h, g = acceleration[i], loss[i], ds[i], gamma[i]
        gamma[i+1] = _advance_constant(g, a, b, h)
        gm = _advance_constant(g, a, b, h/2)
        gamma_mid[i] = gm
        # Exact work minus exact change of particle energy. Roundoff can
        # make a vanishing loss slightly negative; clip only at zero.
        radiated[i] = np.maximum(0., work_erg[i]-MEC2_ERG*(gamma[i+1]-g))
        ec[i] = critical_energy_gev(gm, rho[i])
        if photons is not None:
            if spectral_nodes == 1:
                photons[i] = curvature_photon_yield(radiated[i], ec[i], energy_bands_gev)
            else:
                node_g = np.stack([gm if fraction == .5 else
                                   _advance_constant(g, a, b, h*fraction)
                                   for fraction in node_locations])
                node_ec = critical_energy_gev(node_g, rho[i])
                power = curvature_power_erg_s(node_g, rho[i])
                weighted_power = power*np.reshape(node_weights, (3,)+(1,)*g.ndim)
                denominator = weighted_power.sum(axis=0)
                fractions = np.divide(weighted_power, denominator,
                                      out=np.zeros_like(weighted_power), where=denominator > 0)
                photons[i] = curvature_photon_yield(
                    radiated[i]*fractions, node_ec, energy_bands_gev).sum(axis=0)
    if np.any(gamma < 1) or not np.all(np.isfinite(gamma)):
        raise ValueError("Trajectory leaves the relativistic domain gamma>=1")
    return CurvatureTransport(gamma, gamma_mid, work_erg, radiated, ec, photons)
