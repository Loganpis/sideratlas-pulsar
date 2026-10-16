"""Physical single-particle synchrotron and isotropic inverse-Compton rates.

These are independent kernels, not a chosen Vela emission mixture. All input
fields, particle energies/pitches and seed-photon densities must be supplied.
The IC kernel is the ultrarelativistic Jones/Blumenthal-Gould expression with
Klein-Nishina recoil (also Lefa et al. 2012, equations 2-5). It angle-averages
an ISOTROPIC seed field; replacing outward stellar photons by this field is a
declared approximation and cannot predict directional Vela IC emission.

No bulk Doppler factor, particle density normalization, pair cascade, detector
response, attenuation or self-consistent energy/pitch evolution is supplied.
CR and SR must not be summed when a resolved trajectory's curvature already
includes its gyromotion. A synchro-curvature treatment is then required.
"""
from dataclasses import dataclass
from functools import lru_cache
from numbers import Integral

import numpy as np

from .radiation import (C_CM_S, E_STATC, GEV_ERG, HBAR_ERG_S, M_E_G,
                        MEC2_ERG, curvature_photon_yield)

K_B_ERG_K = 1.380649e-16
SIGMA_T_CM2 = 8 * np.pi / 3 * (E_STATC**2 / MEC2_ERG)**2
B_QUANTUM_GAUSS = M_E_G**2 * C_CM_S**3 / (E_STATC * HBAR_ERG_S)


@lru_cache(maxsize=12)
def _legendre(nodes):
    if isinstance(nodes, bool) or not isinstance(nodes, Integral) or nodes < 16:
        raise ValueError('Quadrature requires an integer number of nodes >=16')
    return np.polynomial.legendre.leggauss(nodes)


def _gamma(gamma):
    gamma = np.asarray(gamma, dtype=float)
    if np.any(~np.isfinite(gamma)) or np.any(gamma < 10):
        raise ValueError('These ultrarelativistic kernels require finite gamma >=10')
    return gamma


def _synch_inputs(gamma, magnetic_field_gauss, pitch_rad):
    g, b, pitch = np.broadcast_arrays(_gamma(gamma), np.asarray(magnetic_field_gauss, float),
                                      np.asarray(pitch_rad, float))
    if (np.any(~np.isfinite(b)) or np.any(b < 0) or np.any(~np.isfinite(pitch))
            or np.any(pitch < 0) or np.any(pitch > np.pi)):
        raise ValueError('Require B>=0 and 0<=pitch<=pi in radians')
    sin_pitch = np.where((pitch == 0) | (pitch == np.pi), 0., np.sin(pitch))
    transverse = b * sin_pitch
    if np.any(g * transverse / B_QUANTUM_GAUSS > .1):
        raise ValueError('Quantum parameter chi_e>0.1: classical synchrotron is inadequate')
    # The broad synchrotron continuum is not a low-harmonic cyclotron solver.
    if np.any((transverse > 0) & (g * sin_pitch < 10)):
        raise ValueError('Require gamma*sin(pitch)>=10 for nonzero synchrotron continuum')
    return g, transverse


def synchrotron_power_erg_s(gamma, magnetic_field_gauss, pitch_rad):
    """Fixed-pitch classical power 2 sigma_T c U_B (gamma²-1) sin²(pitch).

    The factor 4/3 applies only after isotropic pitch averaging. The classical
    regime needs chi_e=gamma B sin(pitch)/B_Q <<1; values above .1 are rejected,
    not corrected. This cutoff is a guard, not an accuracy guarantee near .1.
    """
    g, transverse = _synch_inputs(gamma, magnetic_field_gauss, pitch_rad)
    return 2 * SIGMA_T_CM2 * C_CM_S * transverse**2 / (8*np.pi) * (g*g - 1)


def synchrotron_critical_energy_gev(gamma, magnetic_field_gauss, pitch_rad):
    """Ultrarelativistic critical energy 3 hbar e B gamma² sin(pitch)/(2 m c)."""
    g, transverse = _synch_inputs(gamma, magnetic_field_gauss, pitch_rad)
    return 1.5 * HBAR_ERG_S * E_STATC * transverse * g*g / (M_E_G*C_CM_S) / GEV_ERG


def synchrotron_photon_rate(gamma, magnetic_field_gauss, pitch_rad, bands_gev):
    """Angle-integrated photons/s/particle in supplied bands, final band axis.

    Synchrotron and curvature emission share F(x)=x integral K_(5/3), but have
    distinct powers and critical energies. Reusing the normalized kernel does
    not assume a curvature-radiation origin for these photons.
    """
    power = synchrotron_power_erg_s(gamma, magnetic_field_gauss, pitch_rad)
    ec = synchrotron_critical_energy_gev(gamma, magnetic_field_gauss, pitch_rad)
    return curvature_photon_yield(power, ec, bands_gev)


@dataclass(frozen=True)
class PhotonField:
    """Isotropic seed photons represented by positive energy quadrature.

    number_weights_cm3 contains n(epsilon_i)*d(epsilon_i), NOT a density per
    energy. A single element can represent a monochromatic field. Values are
    absolute local photon densities. Caller must justify angular isotropy.
    """
    energy_erg: np.ndarray
    number_weights_cm3: np.ndarray
    description: str = 'Caller-supplied isotropic photon field'

    def __post_init__(self):
        energy = np.array(self.energy_erg, dtype=float, copy=True)
        weight = np.array(self.number_weights_cm3, dtype=float, copy=True)
        if (energy.ndim != 1 or len(energy) == 0 or energy.shape != weight.shape
                or np.any(~np.isfinite(energy)) or np.any(energy <= 0)
                or np.any(~np.isfinite(weight)) or np.any(weight < 0)):
            raise ValueError('Positive finite seed energies and nonnegative density weights required')
        energy.flags.writeable = False
        weight.flags.writeable = False
        object.__setattr__(self, 'energy_erg', energy)
        object.__setattr__(self, 'number_weights_cm3', weight)

    @property
    def number_density_cm3(self):
        return float(np.sum(self.number_weights_cm3))

    @property
    def energy_density_erg_cm3(self):
        return float(self.energy_erg @ self.number_weights_cm3)


def thermal_photon_density_per_erg(energy_erg, temperature_k, dilution):
    """W epsilon²/[pi² (hbar c)³ (exp(epsilon/kT)-1)] in cm^-3 erg^-1."""
    energy = np.asarray(energy_erg, dtype=float)
    if (not np.isfinite(temperature_k) or temperature_k <= 0
            or not np.isfinite(dilution) or not 0 <= dilution <= 1
            or np.any(~np.isfinite(energy)) or np.any(energy < 0)):
        raise ValueError('Require T>0, 0<=dilution<=1 and finite energies >=0')
    x = energy / (K_B_ERG_K * temperature_k)
    occupation = np.divide(np.exp(-x), -np.expm1(-x), out=np.zeros_like(x), where=x > 0)
    return dilution * energy**2 * occupation / (np.pi**2 * (HBAR_ERG_S*C_CM_S)**3)


def thermal_photon_field(temperature_k, dilution, nodes=128):
    """Planck seed quadrature over 1e-7 <= epsilon/kT <=80.

    Temperature and dilution are explicit inputs, with no Vela defaults.
    Log-energy Gaussian quadrature preserves the low-energy tail; increasing
    nodes provides a convergence test. Omitted number/energy tails are tiny.
    """
    x, weight = _legendre(nodes)
    lower, upper = np.log(1e-7), np.log(80.)
    logenergy = (upper + lower)/2 + (upper - lower)/2*x
    energy = K_B_ERG_K*temperature_k*np.exp(logenergy)
    density = thermal_photon_density_per_erg(energy, temperature_k, dilution)
    return PhotonField(energy, density*energy*weight*(upper-lower)/2,
                       f'Isotropic diluted Planck field; T={temperature_k:g} K, W={dilution:g}')


def stellar_dilution(radius_cm, stellar_radius_cm):
    """Solid-angle fraction W=(1-sqrt(1-(R/r)^2))/2 for r>=R.

    This gives the correct flat-space energy density of a uniform spherical
    blackbody, but its photons travel outward rather than isotropically.
    General relativity, surface temperature variation and atmosphere beaming
    are absent. The algebra avoids cancellation when r greatly exceeds R.
    """
    r, star = np.broadcast_arrays(np.asarray(radius_cm, float), np.asarray(stellar_radius_cm, float))
    if (np.any(~np.isfinite(r)) or np.any(~np.isfinite(star))
            or np.any(star <= 0) or np.any(r < star)):
        raise ValueError('Require finite r>=stellar_radius>0')
    ratio2 = (star/r)**2
    return ratio2 / (2*(1 + np.sqrt(1-ratio2)))


def isotropized_stellar_photon_field(temperature_k, stellar_radius_cm, radius_cm, nodes=128):
    """Preserve stellar photon density while explicitly discarding directionality."""
    dilution = float(stellar_dilution(radius_cm, stellar_radius_cm))
    field = thermal_photon_field(temperature_k, dilution, nodes)
    return PhotonField(field.energy_erg, field.number_weights_cm3,
                       field.description + '; ISOTROPIZED stellar photons, not an anisotropic IC solution')


def _check_ic(gamma, field):
    g = _gamma(gamma)
    if not isinstance(field, PhotonField):
        raise TypeError('Supply a PhotonField with absolute local seed density')
    if np.max(field.energy_erg) > .01 * float(np.min(g)) * MEC2_ERG:
        raise ValueError('Seed photons must be soft relative to the electron energy')
    return g


def _kn_kernel(q, recoil):
    """Positive Jones/BG ultrarelativistic kernel on its kinematic support."""
    return np.maximum(0., 2*q*np.log(q) + (1+2*q)*(1-q)
                      + .5*(recoil*q)**2*(1-q)/(1+recoil*q))


def inverse_compton_spectrum_per_gev(energy_gev, gamma, field):
    """Angle-integrated dN/(dt dE_GeV), per particle, with full KN term.

    Output has gamma.shape + energy_gev.shape (energy must be one-dimensional).
    BG is ultrarelativistic and omits the negligible down-scattering branch;
    it is not the exact arbitrary-energy Compton redistribution function.
    """
    g = _check_ic(gamma, field)
    energy = np.asarray(energy_gev, float)
    if energy.ndim != 1 or np.any(~np.isfinite(energy)) or np.any(energy < 0):
        raise ValueError('Supply a one-dimensional finite photon-energy grid >=0')
    out = np.zeros(g.shape + energy.shape)
    eps = field.energy_erg
    output = out.reshape(-1, len(energy))
    for i, electron in enumerate(g.ravel()):
        recoil = 4*electron*eps/MEC2_ERG
        emitted = energy[:, None]*GEV_ERG
        remain = electron*MEC2_ERG-emitted
        q = np.divide(emitted, recoil*remain,
                      out=np.zeros((len(energy),len(eps))), where=remain > 0)
        valid = (q >= 1/(4*electron**2)) & (q <= 1) & (remain > 0)
        qsafe = np.where(valid, q, 1.)
        kernel = np.where(valid, _kn_kernel(qsafe, recoil), 0.)
        output[i] = (3*SIGMA_T_CM2*C_CM_S/(4*electron**2)
                     * (kernel @ (field.number_weights_cm3/eps)) * GEV_ERG)
    return out


@dataclass(frozen=True)
class InverseComptonRates:
    photon_rate_s: np.ndarray
    band_photon_rate_s: np.ndarray
    emitted_power_erg_s: np.ndarray
    intercepted_seed_power_erg_s: np.ndarray
    electron_loss_erg_s: np.ndarray


def inverse_compton_rates(gamma, field, bands_gev=(), q_nodes=96):
    """Total and band photon rates and NET electron losses in erg/s/particle.

    Integrates log(q) with exact kinematic band boundaries, avoiding energy
    histograms. Electron losses subtract incident seed energy from emitted
    photon energy. In extreme KN scattering this is an ensemble-average loss;
    continuous cooling alone does not reproduce a stochastic electron cascade.
    Arrays follow gamma.shape; band_photon_rate_s appends the band axis.
    """
    g = _check_ic(gamma, field)
    bands = np.asarray(bands_gev, float).reshape(-1, 2)
    if (np.any(np.isnan(bands)) or np.any(bands[:, 0] < 0)
            or np.any(bands[:, 1] <= bands[:, 0])):
        raise ValueError('Bands require 0<=minimum<maximum in GeV')
    node, weight = _legendre(q_nodes)
    total, power, incident = (np.zeros(g.size) for _ in range(3))
    band_rates = np.zeros((g.size, len(bands)))
    eps = field.energy_erg
    for i, electron in enumerate(g.ravel()):
        recoil = 4*electron*eps/MEC2_ERG
        qmin = 1/(4*electron**2)
        prefactor = 3*SIGMA_T_CM2*C_CM_S/(4*electron**2)*field.number_weights_cm3/eps

        def integrate(lower, upper):
            lower, upper = np.broadcast_arrays(lower, upper, subok=False)
            loglo, loghi = np.log(lower), np.log(upper)
            width = np.maximum(0., loghi-loglo)
            q = np.exp(loglo[..., None] + width[..., None]*(node+1)/2)
            emitted = electron*MEC2_ERG*recoil[..., None]*q/(1+recoil[..., None]*q)
            jacobian = electron*MEC2_ERG*recoil[..., None]*q/(1+recoil[..., None]*q)**2
            probability = (prefactor[..., None]*_kn_kernel(q, recoil[..., None])
                           *jacobian*weight*width[..., None]/2)
            return probability, emitted

        probability, emitted = integrate(np.full(len(eps), qmin), np.ones(len(eps)))
        total[i] = probability.sum()
        power[i] = np.sum(probability*emitted)
        incident[i] = np.sum(probability*eps[:, None])
        for j, (elo, ehi) in enumerate(bands):
            def q_at_energy(energy):
                if energy*GEV_ERG >= electron*MEC2_ERG:
                    return np.ones_like(eps)
                value = energy*GEV_ERG/(recoil*(electron*MEC2_ERG-energy*GEV_ERG))
                return np.clip(value, qmin, 1.)
            lower, upper = q_at_energy(elo), q_at_energy(ehi)
            band_rates[i, j] = integrate(lower, upper)[0].sum()
    return InverseComptonRates(total.reshape(g.shape), band_rates.reshape(g.shape+(len(bands),)),
                               power.reshape(g.shape), incident.reshape(g.shape),
                               (power-incident).reshape(g.shape))


def inverse_compton_loss_erg_s(gamma, field, q_nodes=96):
    """Net ensemble-average electron energy transferred to scattered photons."""
    return inverse_compton_rates(gamma, field, q_nodes=q_nodes).electron_loss_erg_s
