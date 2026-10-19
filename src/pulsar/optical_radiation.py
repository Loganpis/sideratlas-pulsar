"""Verified radiation building blocks, not a Vela optical-camera model.

All quantities in this module use SI units.  The Lienard--Wiechert functions
calculate *bolometric* classical radiation of supplied physical trajectories
in flat spacetime.  The synchrotron functions are a separate, angle-integrated
ultrarelativistic continuum approximation in a local magnetic frame.  Multiplying
the two does NOT produce a directional optical spectrum.

Inputs must already describe physical particles.  A reduced-scale PIC momentum,
magnetic field, radius or macroparticle charge is not automatically a physical
Vela value.  No pair cascade, pitch evolution, coherent emission, atmosphere,
general-relativistic propagation, extinction or detector response is supplied.

References and the source/observer-time distinction are documented in
research/OBSERVER-RADIATION-REBUILD.md.  Relevant primary references include
Watt & O'Shea (2024), doi:10.1107/S1600577524000663, and Lai et al. (2025),
doi:10.1093/mnras/staf1295.  Incoherent optical radiation sums physical-particle
intensities with number weights, not squared macroparticle charges.
"""
from dataclasses import dataclass
from functools import lru_cache
from numbers import Integral

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as gamma_function, kv


C_M_S = 299_792_458.0
E_CHARGE_C = 1.602_176_634e-19
EPSILON_0_F_M = 8.854_187_8128e-12
M_E_KG = 9.109_383_7015e-31
HBAR_J_S = 1.054_571_817e-34
B_QUANTUM_T = M_E_KG**2 * C_M_S**2 / (E_CHARGE_C * HBAR_J_S)
SYNCHROTRON_F_INTEGRAL = 8 * np.pi / (9 * np.sqrt(3))


def _vector(value, name):
    value = np.asarray(value, dtype=float)
    if value.ndim == 0 or value.shape[-1] != 3 or np.any(~np.isfinite(value)):
        raise ValueError(f'{name} must have a final dimension of three finite components')
    return value


def _kinematics(beta, beta_dot_s):
    beta, acceleration = np.broadcast_arrays(_vector(beta, 'beta'),
                                             _vector(beta_dot_s, 'beta_dot_s'))
    beta2 = np.sum(beta*beta, axis=-1)
    if np.any(beta2 >= 1):
        raise ValueError('Require |beta|<1; a velocity rounded to c is not usable kinematics')
    return beta, acceleration, 1-beta2


def _direction(observer_direction):
    direction = _vector(observer_direction, 'observer_direction')
    norm = np.linalg.norm(direction, axis=-1)
    if not np.allclose(norm, 1., rtol=0, atol=1e-12):
        raise ValueError('observer_direction must be a unit vector, pointing toward the observer')
    # Unit roundoff must not spoil the cancellation-resistant Doppler factor.
    return direction / norm[..., None]


def _charge(charge_c):
    charge = np.asarray(charge_c, dtype=float)
    if np.any(~np.isfinite(charge)):
        raise ValueError('Physical charge in coulombs must be finite')
    return charge


def _doppler_denominator(beta, direction, inverse_gamma2):
    # Exactly 1-n.beta for unit n, while avoiding subtracting nearly equal
    # dot products.  beta itself must retain its physical subluminal speed.
    return .5*(np.sum((direction-beta)**2, axis=-1) + inverse_gamma2)


def lienard_total_power_w(beta, beta_dot_s, charge_c=-E_CHARGE_C):
    """Exact classical instantaneous emitted power, in W per physical particle.

    beta=v/c and beta_dot_s=d(beta)/dt in inverse seconds are laboratory-frame
    kinematics at the same source event.  No Lorentz factor or acceleration is
    inferred from track geometry, interpolated fields or a presumed pitch.
    """
    beta, acceleration, inv_g2 = _kinematics(beta, beta_dot_s)
    a2 = np.sum(acceleration*acceleration, axis=-1)
    parallel = np.sum(beta*acceleration, axis=-1)
    # Stable form of gamma^6*(a^2-|beta cross a|^2).
    return (_charge(charge_c)**2/(6*np.pi*EPSILON_0_F_M*C_M_S)
            * (a2/inv_g2**2 + parallel**2/inv_g2**3))


def lienard_angular_power_source_w_sr(beta, beta_dot_s, observer_direction,
                                     charge_c=-E_CHARGE_C):
    """Bolometric dE/(dt_source dOmega), W/sr per physical particle.

    This is the fifth-power-denominator Lienard--Wiechert expression, not an
    optical spectrum.  Time-integrate it in *source time* and then deposit the
    resulting energy at retarded arrival times.  Broadcasting follows NumPy's
    rules over the leading dimensions of the vector inputs.
    """
    beta, acceleration, inv_g2 = _kinematics(beta, beta_dot_s)
    direction = _direction(observer_direction)
    numerator = np.cross(direction, np.cross(direction-beta, acceleration))
    denominator = _doppler_denominator(beta, direction, inv_g2)
    return (_charge(charge_c)**2/(16*np.pi**2*EPSILON_0_F_M*C_M_S)
            * np.sum(numerator*numerator, axis=-1)/denominator**5)


def lienard_angular_power_received_w_sr(beta, beta_dot_s, observer_direction,
                                       charge_c=-E_CHARGE_C):
    """Bolometric dE/(dt_arrival dOmega), W/sr, for a fixed far observer.

    The sixth power includes dt_arrival/dt_source=1-n.beta.  Integrating this
    function against source-time intervals would count the Jacobian twice.
    Geometric dilution by the observer distance squared is not included.
    """
    checked_beta, _, inv_g2 = _kinematics(beta, beta_dot_s)
    direction = _direction(observer_direction)
    return (lienard_angular_power_source_w_sr(beta, beta_dot_s, direction, charge_c)
            / _doppler_denominator(checked_beta, direction, inv_g2))


def retarded_arrival_times_s(source_time_s, positions_m, observer_direction):
    """Far-observer arrival times t-n.r/c, subtracting the constant distance/c.

    The direction is fixed from the source toward the observer.  This is a
    flat-spacetime far-field construction, not finite-distance or GR ray tracing.
    For arrays shaped (particles,times,3), times may be a matching one-dimensional
    array; particle directions or other broadcasting axes must be explicit.
    """
    time = np.asarray(source_time_s, dtype=float)
    if np.any(~np.isfinite(time)):
        raise ValueError('Source times must be finite seconds')
    position = _vector(positions_m, 'positions_m')
    direction = _direction(observer_direction)
    return time-np.sum(position*direction, axis=-1)/C_M_S


def incoherent_sum(single_particle_values, particle_number_weights, axis=0):
    """Sum nonnegative physical-particle intensities with linear number weights.

    Weights are the number of statistically incoherent physical particles per
    sample, NOT charge in coulombs.  A length-N vector weights the selected
    axis; all remaining dimensions (time, frequency, pixels) are preserved.
    Identical trajectories with random optical phases contribute intensities,
    not phase-coherent electromagnetic amplitudes.
    """
    values = np.asarray(single_particle_values, dtype=float)
    weights = np.asarray(particle_number_weights, dtype=float)
    if values.ndim == 0:
        raise ValueError('Supply a particle axis, even for one physical particle')
    if not isinstance(axis, Integral) or isinstance(axis, bool):
        raise ValueError('axis must identify the particle dimension')
    axis = int(axis)
    if not -values.ndim <= axis < values.ndim:
        raise ValueError('axis must identify the particle dimension')
    axis %= values.ndim
    if (weights.ndim != 1 or len(weights) != values.shape[axis]
            or np.any(~np.isfinite(values)) or np.any(values < 0)
            or np.any(~np.isfinite(weights)) or np.any(weights < 0)):
        raise ValueError('Require nonnegative finite values and one number weight per particle')
    shape = [1]*values.ndim
    shape[axis] = len(weights)
    return np.sum(values*weights.reshape(shape), axis=axis)


@dataclass(frozen=True)
class ReceivedEnergyBins:
    """Bolometric energy per steradian and its arrival-bin mean power.

    power_w_sr is not a spectrum, pixel brightness, flux at Earth, or optical
    camera output.  Its arrival-time integral equals energy_j_sr.sum().
    """
    bin_edges_s: np.ndarray
    energy_j_sr: np.ndarray
    power_w_sr: np.ndarray


def bin_received_energy(arrival_time_s, emitted_energy_j_sr, bin_edges_s,
                        particle_number_weights=1.):
    """Conservatively histogram event energies into finite arrival-time bins.

    The sample energies have already been integrated over source time.  Number
    weights must broadcast to the event array (use weights[:,None] for an array
    shaped particles,times).  Nonzero energy outside the bins raises an error
    instead of disappearing silently.  The last edge is included, following
    numpy.histogram.  This bins point events; it is not a finite-exposure
    reconstruction of unresolved emission between events.
    """
    arrival, energy, weights = np.broadcast_arrays(np.asarray(arrival_time_s, float),
                                                  np.asarray(emitted_energy_j_sr, float),
                                                  np.asarray(particle_number_weights, float))
    edges = np.asarray(bin_edges_s, dtype=float)
    if (edges.ndim != 1 or len(edges) < 2 or np.any(~np.isfinite(edges))
            or np.any(np.diff(edges) <= 0)):
        raise ValueError('Arrival bin edges must be finite, one-dimensional and strictly increasing')
    if (np.any(~np.isfinite(arrival)) or np.any(~np.isfinite(energy)) or np.any(energy < 0)
            or np.any(~np.isfinite(weights)) or np.any(weights < 0)):
        raise ValueError('Require finite arrival times, nonnegative energies and number weights')
    weighted = energy*weights
    if np.any(~np.isfinite(weighted)):
        raise ValueError('Weighted emitted energies must be finite')
    outside = (arrival < edges[0]) | (arrival > edges[-1])
    if np.any(outside & (weighted != 0)):
        raise ValueError('Nonzero emitted energy lies outside the supplied arrival bins')
    received = np.histogram(arrival.ravel(), edges, weights=weighted.ravel())[0]
    edges = edges.copy()
    power = received/np.diff(edges)
    for value in (edges, received, power):
        value.flags.writeable = False
    return ReceivedEnergyBins(edges, received, power)


def bin_lienard_radiation(source_time_s, positions_m, beta, beta_dot_s,
                          observer_direction, source_intervals_s, bin_edges_s,
                          particle_number_weights=1., charge_c=-E_CHARGE_C):
    """Bin sampled bolometric LW radiation without duplicating the Jacobian.

    source_intervals_s are caller-supplied physical quadrature weights in
    source time.  This function does not infer intervals from an irregular
    trajectory or assume that accelerations and interpolated paths are mutually
    consistent.  Radiation and particle-sampling convergence remain caller work.
    """
    intervals = np.asarray(source_intervals_s, dtype=float)
    if np.any(~np.isfinite(intervals)) or np.any(intervals < 0):
        raise ValueError('Source quadrature intervals must be finite and nonnegative seconds')
    power = lienard_angular_power_source_w_sr(beta, beta_dot_s, observer_direction, charge_c)
    arrival = retarded_arrival_times_s(source_time_s, positions_m, observer_direction)
    return bin_received_energy(arrival, power*intervals, bin_edges_s, particle_number_weights)


def _synchrotron_inputs(gamma, magnetic_field_t, pitch_rad):
    g, b, pitch = np.broadcast_arrays(np.asarray(gamma, float),
                                      np.asarray(magnetic_field_t, float),
                                      np.asarray(pitch_rad, float))
    if (np.any(~np.isfinite(g)) or np.any(g < 10) or np.any(~np.isfinite(b))
            or np.any(b < 0) or np.any(~np.isfinite(pitch))
            or np.any(pitch < 0) or np.any(pitch > np.pi)):
        raise ValueError('Require finite gamma>=10, B>=0 tesla and 0<=pitch<=pi radians')
    sin_pitch = np.where((pitch == 0) | (pitch == np.pi), 0., np.sin(pitch))
    transverse = b*sin_pitch
    if np.any((transverse > 0) & (g*sin_pitch < 10)):
        raise ValueError('Synchrotron continuum requires gamma*sin(pitch)>=10; not cyclotron radiation')
    if np.any(g*transverse/B_QUANTUM_T > .1):
        raise ValueError('Quantum parameter chi_e>0.1; classical synchrotron is inadequate')
    return g, transverse


def synchrotron_power_w(gamma, magnetic_field_t, pitch_rad):
    """Angle-integrated fixed-pitch Larmor power in the local magnetic frame.

    P=e^4 B^2 (gamma^2-1) sin^2(pitch)/(6 pi epsilon_0 m_e^2 c).
    No isotropic-pitch average or bulk Doppler factor is assumed.
    """
    g, transverse = _synchrotron_inputs(gamma, magnetic_field_t, pitch_rad)
    return E_CHARGE_C**4*transverse**2*(g*g-1)/(6*np.pi*EPSILON_0_F_M*M_E_KG**2*C_M_S)


def synchrotron_critical_frequency_hz(gamma, magnetic_field_t, pitch_rad):
    """Ultrarelativistic nu_c=3 e B gamma^2 sin(pitch)/(4 pi m_e), in Hz."""
    g, transverse = _synchrotron_inputs(gamma, magnetic_field_t, pitch_rad)
    return 3*E_CHARGE_C*transverse*g*g/(4*np.pi*M_E_KG)


def synchrotron_f(x):
    """Broad kernel F(x)=x integral_x^infinity K_(5/3)(u) du.

    A logarithmic-variable adaptive numerical quadrature resolves x>=1e-10.
    Below that threshold the two-term small-x asymptote avoids Bessel overflow;
    at x>700 the exponentially negligible tail is returned as zero.  The latter
    convention cannot represent arbitrary-precision far-tail relative errors.
    """
    x = np.asarray(x, dtype=float)
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError('Synchrotron kernel requires finite x>=0')
    result = np.zeros_like(x)
    for index in np.ndindex(x.shape):
        value = float(x[index])
        if value == 0 or value > 700:
            continue
        if value < 1e-10:
            result[index] = (2**(2/3)*gamma_function(2/3)*value**(1/3)
                             - np.pi/np.sqrt(3)*value)
            continue
        # For large x the relevant integration interval is only O(1) wide;
        # exponential suppression bounds the omitted tail at x+80.
        lower, upper = np.log(value), np.log(max(80., value+80.))
        integral = quad(lambda logu: np.exp(logu)*kv(5/3, np.exp(logu)),
                        lower, upper, epsabs=1e-300, epsrel=2e-10, limit=150)[0]
        result[index] = value*integral
    return result


def synchrotron_spectrum_w_hz(frequency_hz, gamma, magnetic_field_t, pitch_rad):
    """ANGLE-INTEGRATED continuum P_nu, W/Hz per physical electron/positron.

    The one-dimensional frequency grid becomes the final output axis.  The
    ultrarelativistic F(nu/nu_c) shape is normalized to the exact fixed-pitch
    Larmor power proportional to gamma^2-1, as in this project's existing
    nonthermal kernel.  This explicitly retains a beta^2 normalization factor;
    it does not make the continuum shape exact for finite gamma.  Neither the
    spectral angular pattern nor photons reaching a camera are calculated.
    """
    frequency = np.asarray(frequency_hz, dtype=float)
    if frequency.ndim != 1 or np.any(~np.isfinite(frequency)) or np.any(frequency < 0):
        raise ValueError('Supply a one-dimensional finite frequency grid >=0 Hz')
    critical = synchrotron_critical_frequency_hz(gamma, magnetic_field_t, pitch_rad)
    power = synchrotron_power_w(gamma, magnetic_field_t, pitch_rad)
    shape = critical.shape + frequency.shape
    ratio = np.divide(frequency, critical[..., None], out=np.zeros(shape),
                      where=critical[..., None] > 0)
    normalization = np.divide(power, critical*SYNCHROTRON_F_INTEGRAL,
                              out=np.zeros_like(power), where=critical > 0)
    return normalization[..., None]*synchrotron_f(ratio)


@lru_cache(maxsize=8)
def _band_quadrature(nodes):
    if isinstance(nodes, bool) or not isinstance(nodes, Integral) or nodes < 16:
        raise ValueError('Band quadrature requires an integer number of nodes >=16')
    return np.polynomial.legendre.leggauss(nodes)


def synchrotron_wavelength_band_fraction(gamma, magnetic_field_t, pitch_rad,
                                        wavelength_min_m=400e-9,
                                        wavelength_max_m=700e-9, nodes=64):
    """ANGLE-INTEGRATED energy fraction; defaults to the 400--700 nm band.

    This is NOT the fraction observed along an arbitrary direction.  Spectral
    and angular beaming cannot generally be separated.  Wavelengths are vacuum
    wavelengths in the same local magnetic frame as the particle inputs.
    """
    if (not np.isscalar(wavelength_min_m) or not np.isscalar(wavelength_max_m)
            or not np.isfinite(wavelength_min_m) or not np.isfinite(wavelength_max_m)
            or not 0 < wavelength_min_m < wavelength_max_m):
        raise ValueError('Require finite 0<wavelength_min_m<wavelength_max_m')
    critical = synchrotron_critical_frequency_hz(gamma, magnetic_field_t, pitch_rad)
    x, weights = _band_quadrature(nodes)
    low, high = np.log(C_M_S/wavelength_max_m), np.log(C_M_S/wavelength_min_m)
    frequency = np.exp((high+low)/2+(high-low)*x/2)
    ratio = np.divide(frequency, critical[..., None],
                      out=np.zeros(critical.shape+frequency.shape),
                      where=critical[..., None] > 0)
    # integral F(nu/nu_c) dnu / [nu_c integral F(x) dx]
    fraction = np.sum(weights*ratio*synchrotron_f(ratio), axis=-1)*(high-low)/2
    fraction /= SYNCHROTRON_F_INTEGRAL
    if np.any((fraction < 0) | (fraction > 1+1e-10)):
        raise ArithmeticError('Synchrotron band quadrature failed to produce a physical fraction')
    return np.minimum(fraction, 1.)
