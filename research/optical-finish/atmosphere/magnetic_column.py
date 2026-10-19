"""Research-only, explicitly approximate normal-field ionized-H transfer column.

This is *not* a validated Vela atmosphere.  Unity Gaunt factors, transverse
linear X/O polarizations, an ideal-gas EOS, coherent scattering, no vacuum
polarization, refraction, partial ionization or conduction are model choices.
The emergent intensity is solved by transfer; it is not replaced or rescaled
by a blackbody. Planck B_nu/2 is the LTE per-mode thermal source.

Free-free electron/proton response: van Adelsberg & Lai 2006, eqs. 13, 24--28
(correcting the older unphysical electron-only low-frequency response).
Scattering: their eqs. 29--31, azimuth averaged in the chosen transverse basis.
LTE energy balance and hydrostatic balance: their eqs. 7--9.
Source: https://arxiv.org/abs/astro-ph/0607168

Use only as a bounded numerical experiment until independent validation,
physical opacity/EOS replacement, and resolution studies have passed.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import c, h, k, m_e, m_p, sigma, electron_volt
from scipy.optimize import least_squares
from scipy.sparse import bmat, csc_matrix, diags
from scipy.sparse.linalg import spsolve

C = c * 100.0
H = h * 1e7
KB = k * 1e7
MP = m_p * 1e3
SIGMA = sigma * 1e3
EV_ERG = electron_volt * 1e7
KEV_ERG = 1e3 * EV_ERG
KAPPA_T = 0.397726
ME_MP = m_e / m_p


@dataclass
class Parameters:
    teff_k: float = 0.7e6 / 0.8258072924
    gravity_cgs: float = 1.3313e14
    magnetic_field_g: float = 3e12
    depth_points: int = 40
    energy_points: int = 80
    angles_per_hemisphere: int = 4
    column_min: float = 1e-9
    column_max: float = 1e6
    energy_min_kev: float = 0.0001
    energy_max_kev: float = 10.0
    gaunt_factor: float = 1.0
    max_nfev: int = 100


def planck_nu(nu: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    x = H * np.asarray(nu) / (KB * np.asarray(temperature))
    return 2 * H * np.asarray(nu) ** 3 / C**2 / np.expm1(np.minimum(x, 700.0))


def cyclic_powers(mu: np.ndarray) -> np.ndarray:
    """Mode, angle, cyclic component(-,+,0); transverse linear X/O basis."""
    out = np.empty((2, len(mu), 3))
    out[0, :, :2] = 0.5
    out[0, :, 2] = 0.0
    out[1, :, 0] = out[1, :, 1] = mu**2 / 2
    out[1, :, 2] = 1 - mu**2
    return out.reshape((-1, 3))


def response_opacities(energy_kev, rho, temperature, magnetic_field_g, gaunt=1.0):
    """Absorption per cyclic component, plus electron+proton scattering factors.

    All opacities are cm²/g. Arrays end in cyclic(-,+,0).  Damping formula
    uses the cited published rounded numerical coefficients. Scattering
    damping follows the same oscillator approximation, with proton collision
    damping reduced by m_e/m_p; it matters only near a resonance.
    """
    if magnetic_field_g <= 0:
        raise ValueError('This unity-magnetic-Gaunt experiment requires B > 0; it has no validated B=0 limit')
    energy = np.asarray(energy_kev, dtype=float)
    rho = np.asarray(rho, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    omega = energy * KEV_ERG / H * (2 * np.pi)
    epe = 0.02871 * np.sqrt(rho)
    ve = (epe / energy) ** 2
    we = 11.577 * (magnetic_field_g / 1e12) / energy
    wi = we * ME_MP
    stim = -np.expm1(-energy * KEV_ERG / KB / temperature)
    gamma = 9.2e-5 * rho / np.sqrt(temperature / 1e6) / energy**2 * stim * gaunt
    gre = 9.5e-6 * energy
    gri = 5.2e-9 * energy
    common = omega / C / rho * ve * gamma
    absorption = np.empty((*np.shape(common), 3))
    scatter = np.empty_like(absorption)
    for j, alpha in enumerate((-1.0, 1.0)):
        gam_pm = gamma + (1 + alpha * we) * gri + (1 - alpha * wi) * gre
        denominator = ((1 + alpha * we) * (1 - alpha * wi)) ** 2 + gam_pm**2
        absorption[..., j] = common / denominator
        scatter[..., j] = (
            1 / ((1 + alpha * we) ** 2 + (gamma + gre) ** 2)
            + ME_MP**2 / ((1 - alpha * wi) ** 2 + (ME_MP * gamma + gri) ** 2)
        )
    absorption[..., 2] = common
    scatter[..., 2] = 1 / (1 + (gamma + gre) ** 2) + ME_MP**2 / (
        1 + (ME_MP * gamma + gri) ** 2
    )
    return absorption, scatter


class Column:
    def __init__(self, p: Parameters):
        self.p = p
        self.mass = np.geomspace(p.column_min, p.column_max, p.depth_points)
        # End points are real boundary nodes. Interior equations are finite
        # volume, using geometric interface placement and exact mass widths.
        self.faces = np.sqrt(self.mass[:-1] * self.mass[1:])
        self.volume = np.diff(np.r_[0.0, self.faces, self.mass[-1]])
        self.energy = np.geomspace(p.energy_min_kev, p.energy_max_kev, p.energy_points)
        self.nu = self.energy * KEV_ERG / H
        edges = np.r_[self.energy[0], np.sqrt(self.energy[:-1] * self.energy[1:]), self.energy[-1]]
        self.dnu = np.diff(edges) * KEV_ERG / H
        xi, wi = leggauss(p.angles_per_hemisphere)
        self.mu, self.w = (xi + 1) / 2, wi / 2
        self.q = cyclic_powers(self.mu)
        self.channel_mu = np.tile(self.mu, 2)
        self.channel_w = np.tile(self.w, 2)
        self.history = []
        self.calls = 0
        self.started = time.monotonic()

    def transfer(self, temperature, keep_intensity=False, boundary_temperature=None,
                 scattering_multiplier=1.0, bottom_intensity_factor=1.0):
        p, mass = self.p, self.mass
        nd, nc = len(mass), len(self.channel_mu)
        rho = MP * p.gravity_cgs * mass / (2 * KB * temperature)
        thermal = planck_nu(self.nu[:, None], temperature[None, :]) / 2
        bottom_b = planck_nu(self.nu, temperature[-1] if boundary_temperature is None else boundary_temperature) / 2 * bottom_intensity_factor
        abs_cyclic, scat_cyclic = response_opacities(
            self.energy[:, None], rho[None, :], temperature[None, :], p.magnetic_field_g, p.gaunt_factor
        )
        absorption = np.einsum('kda,ca->kdc', abs_cyclic, self.q)
        kernel = scattering_multiplier * 1.5 * KAPPA_T * np.einsum('kda,ca,ba,b->kdcb', scat_cyclic, self.q, self.q, self.channel_w)
        scattering_out = np.sum(kernel, axis=-1)
        total = absorption + scattering_out
        all_u = np.zeros((len(self.nu), nd, nc))
        heating = np.zeros(nd)
        emitted = np.zeros(nd)
        flux = np.zeros(nd + 1)
        surface_i = np.empty((len(self.nu), nc))
        for ki in range(len(self.nu)):
            # Exact constant-source cell integrals, including cell-averaged J.
            # Unlike a differential discretization in extreme optical depths,
            # this response contains only bounded exponentials.  No spectrum
            # is normalized after the transfer solve.
            tau = total[ki] * self.volume[:, None] / self.channel_mu[None, :]
            trans = np.exp(-np.minimum(tau, 750.0))
            absorb = -np.expm1(-np.minimum(tau, 750.0))
            phi = np.divide(absorb, tau, out=np.ones_like(tau), where=tau > 1e-300)
            local = 1 - phi
            tiny = tau < 1e-4
            local[tiny] = tau[tiny] / 2 - tau[tiny]**2 / 6 + tau[tiny]**3 / 24
            response = []
            boundary = np.empty((nc, nd))
            for channel in range(nc):
                lam = np.diag(local[:, channel])
                for j in range(nd):
                    attenuation = 1.0
                    for ell in range(j - 1, -1, -1):
                        lam[j, ell] = 0.5 * phi[j, channel] * attenuation * absorb[ell, channel]
                        attenuation *= trans[ell, channel]
                    attenuation = 1.0
                    for ell in range(j + 1, nd):
                        lam[j, ell] = 0.5 * phi[j, channel] * attenuation * absorb[ell, channel]
                        attenuation *= trans[ell, channel]
                    boundary[channel, j] = 0.5 * phi[j, channel] * attenuation * bottom_b[ki]
                response.append(lam)
            blocks = [[None] * nc for _ in range(nc)]
            rhs = np.empty((nc, nd))
            for a in range(nc):
                rhs[a] = response[a] @ (absorption[ki, :, a] / total[ki, :, a] * thermal[ki]) + boundary[a]
                for b in range(nc):
                    operator = -response[a] * (kernel[ki, :, a, b] / total[ki, :, a])[None, :]
                    if a == b:
                        operator = operator + np.eye(nd)
                    blocks[a][b] = csc_matrix(operator)
            solved = spsolve(bmat(blocks, format='csc'), rhs.ravel()).reshape(nc, nd).T
            all_u[ki] = solved
            source = (absorption[ki] * thermal[ki, :, None] + np.einsum('dab,db->da', kernel[ki], solved)) / total[ki]
            incoming = np.zeros((nd + 1, nc))
            outgoing = np.zeros_like(incoming)
            outgoing[-1] = bottom_b[ki]
            for j in range(nd):
                incoming[j + 1] = incoming[j] * trans[j] + source[j] * absorb[j]
            for j in range(nd - 1, -1, -1):
                outgoing[j] = outgoing[j + 1] * trans[j] + source[j] * absorb[j]
            surface_i[ki] = outgoing[0]
            flux += 2 * np.pi * np.sum(self.channel_w[None, :] * self.channel_mu[None, :] * (outgoing - incoming), axis=-1) * self.dnu[ki]
            aw = absorption[ki] * self.channel_w[None, :]
            heating += np.sum(aw * (solved - thermal[ki, :, None]), axis=-1) * self.dnu[ki]
            emitted += np.sum(aw * thermal[ki, :, None], axis=-1) * self.dnu[ki]
        residual = heating / np.maximum(emitted, 1e-100)
        output = dict(rho=rho, flux=flux, heating=heating, emitted=emitted,
                      residual=residual, absorption=absorption, total=total)
        if keep_intensity:
            output['u'] = all_u
            output['surface_i_nu_cgs'] = surface_i
        return output

    def residual(self, log_temperature):
        self.calls += 1
        t = np.exp(log_temperature[:-1])
        out = self.transfer(t, boundary_temperature=np.exp(log_temperature[-1]))
        target = SIGMA * self.p.teff_k**4
        residual = np.r_[out['residual'], np.log(max(out['flux'][0], 1e-100) / target)]
        if self.calls % (self.p.depth_points + 2) == 1:
            state = dict(calls=self.calls, seconds=time.monotonic() - self.started,
                         max_local_relative=float(np.max(np.abs(residual[:-1]))),
                         surface_flux_ratio=float(out['flux'][0] / target),
                         min_t=float(t.min()), max_t=float(t.max()))
            self.history.append(state)
            print(json.dumps(state), flush=True)
        return residual

    def solve(self, initial_temperature=None):
        if initial_temperature is None:
            initial_temperature = self.p.teff_k * (
                0.75 * (0.05 + 0.004 * self.mass)
            )**0.25
        initial_temperature = np.r_[initial_temperature, initial_temperature[-1] * 1.05]
        result = least_squares(self.residual, np.log(initial_temperature),
                               bounds=(np.log(self.p.teff_k * 0.03), np.log(self.p.teff_k * 100)),
                               ftol=2e-9, xtol=2e-9, gtol=2e-9,
                               diff_step=1e-4, max_nfev=self.p.max_nfev)
        return result, self.transfer(np.exp(result.x[:-1]), True, boundary_temperature=np.exp(result.x[-1]))


def save_run(destination: Path, p: Parameters, column: Column, result, output):
    destination.mkdir(parents=True, exist_ok=True)
    t = np.exp(result.x[:-1])
    target = SIGMA * p.teff_k**4
    flux_errors = output['flux'] / target - 1
    local_error = float(np.max(np.abs(output['residual'])))
    flux_error = float(np.max(np.abs(flux_errors)))
    # Numerical equilibrium is explicitly separate from validated physical
    # microphysics / a scientifically accepted atmosphere output.
    numerical_equilibrium = bool(result.success and local_error < 1e-4 and flux_error < 0.01)
    report = dict(
        created_utc=datetime.now(timezone.utc).isoformat(), parameters=asdict(p),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime_s=time.monotonic() - column.started, residual_evaluations=column.calls,
        optimizer_success=bool(result.success), optimizer_message=result.message,
        numerical_equilibrium_pass=numerical_equilibrium,
        scientific_validation_pass=False, production_optical_surface_ready=False,
        max_local_energy_balance_relative=local_error,
        max_flux_error_over_all_depths=flux_error,
        surface_flux_ratio=float(output['flux'][0] / target),
        min_temperature_k=float(t.min()), max_temperature_k=float(t.max()),
        bottom_boundary_temperature_k=float(np.exp(result.x[-1])),
        min_density_g_cm3=float(output['rho'].min()), max_density_g_cm3=float(output['rho'].max()),
        assumptions=['fully ionized H', 'ideal gas hydrostatic equilibrium',
                     'field normal to surface', 'unity magnetic Gaunt factors',
                     'linear transverse X/O basis', 'coherent mode-changing scattering',
                     'no partial ionization, bound opacities, vacuum conversion, refraction or conduction',
                     'no independently validated magnetic optical atmosphere benchmark'],
        history=column.history,
    )
    (destination / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    np.savez_compressed(destination / 'column.npz',
        mass_column_g_cm2=column.mass, temperature_k=t, bottom_boundary_temperature_k=np.exp(result.x[-1]), density_g_cm3=output['rho'],
        flux_cgs=output['flux'], local_relative_heating=output['residual'],
        frequency_hz=column.nu, d_frequency_hz=column.dnu, mu=column.mu, mu_weights=column.w,
        surface_i_nu_cgs=output['surface_i_nu_cgs'].reshape((-1, 2, len(column.mu))),
        absorption_cgs=output['absorption'], total_opacity_cgs=output['total'],
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--depths', type=int, default=40)
    parser.add_argument('--energies', type=int, default=80)
    parser.add_argument('--angles', type=int, default=4)
    parser.add_argument('--max-nfev', type=int, default=100)
    parser.add_argument('--initial', type=Path)
    args = parser.parse_args()
    p = Parameters(depth_points=args.depths, energy_points=args.energies,
                   angles_per_hemisphere=args.angles, max_nfev=args.max_nfev)
    column = Column(p)
    initial = None
    if args.initial:
        previous = np.load(args.initial)
        initial = np.exp(np.interp(np.log(column.mass), np.log(previous['mass_column_g_cm2']), np.log(previous['temperature_k'])))
    result, output = column.solve(initial)
    print(json.dumps(save_run(args.output, p, column, result, output), indent=2), flush=True)


if __name__ == '__main__':
    main()
