"""Linear-source variant of the bounded two-mode magnetic transfer experiment.

The source varies linearly in each interval's optical depth. A diffusion
boundary avoids turning deep, opaque cells into artificial uniform radiators.
The nonlinear temperature iteration enforces the calculated bolometric flux
at every depth, with a separate local LTE energy-balance residual audit.
Physics approximations remain those of magnetic_column.py; this is not a
validated Vela surface atmosphere.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import json
import numpy as np
from scipy.sparse import bmat, csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.optimize import least_squares

from magnetic_column import (Column, Parameters, response_opacities, planck_nu,
                             KAPPA_T, MP, KB, SIGMA, H, KEV_ERG, save_run)


class LinearColumn(Column):
    def __init__(self, p, optical_nodes=30):
        super().__init__(p)
        # Explicitly resolve local optical wavelengths; no X-ray extrapolation.
        optical_e = 1.239841984 / np.linspace(0.300, 0.650, optical_nodes) / 1000
        self.energy = np.unique(np.r_[self.energy, optical_e])
        self.nu = self.energy * KEV_ERG / H
        edges = np.r_[self.energy[0], np.sqrt(self.energy[:-1] * self.energy[1:]), self.energy[-1]]
        self.dnu = np.diff(edges) * KEV_ERG / H
        self.checkpoint = None

    def transfer(self, temperature, keep_intensity=False, **unused):
        p, mass = self.p, self.mass
        nd, nc = len(mass), len(self.channel_mu)
        rho = MP * p.gravity_cgs * mass / (2 * KB * temperature)
        thermal = planck_nu(self.nu[:, None], temperature[None, :]) / 2
        abs_cyclic, scat_cyclic = response_opacities(
            self.energy[:, None], rho[None, :], temperature[None, :], p.magnetic_field_g, p.gaunt_factor
        )
        absorption = np.einsum('kda,ca->kdc', abs_cyclic, self.q)
        kernel = 1.5 * KAPPA_T * np.einsum('kda,ca,ba,b->kdcb', scat_cyclic, self.q, self.q, self.channel_w)
        total = absorption + np.sum(kernel, axis=-1)
        all_u = np.zeros((len(self.nu), nd, nc))
        heating = np.zeros(nd)
        emitted = np.zeros(nd)
        flux = np.zeros(nd)
        surface_i = np.empty((len(self.nu), nc))
        identity = np.eye(nd)
        for ki in range(len(self.nu)):
            tau = 0.5 * (total[ki, :-1] + total[ki, 1:]) * np.diff(mass)[:, None] / self.channel_mu[None, :]
            trans = np.exp(-np.minimum(tau, 750.0))
            one_minus = -np.expm1(-np.minimum(tau, 750.0))
            phi = np.divide(one_minus, tau, out=np.ones_like(tau), where=tau > 1e-300)
            right = 1 - phi
            tiny = tau < 1e-4
            right[tiny] = tau[tiny] / 2 - tau[tiny]**2 / 6 + tau[tiny]**3 / 24
            left = one_minus - right
            lambda_u, lambda_plus, lambda_minus = [], [], []
            for a in range(nc):
                minus = np.zeros((nd, nd))
                plus = np.zeros((nd, nd))
                for j in range(1, nd):
                    minus[j] = trans[j - 1, a] * minus[j - 1] + left[j - 1, a] * identity[j - 1] + right[j - 1, a] * identity[j]
                # Diffusive incoming bottom intensity, with source thermalized.
                plus[-1] = identity[-1] + (identity[-1] - identity[-2]) / tau[-1, a]
                for j in range(nd - 2, -1, -1):
                    plus[j] = trans[j, a] * plus[j + 1] + right[j, a] * identity[j] + left[j, a] * identity[j + 1]
                lambda_minus.append(minus)
                lambda_plus.append(plus)
                lambda_u.append(0.5 * (minus + plus))
            # Scattering couples through only the three cyclic moments.
            # Solve 3*Ndepth unknown moments, then recover every angular ray.
            # This is algebraically identical to the larger mode/angle system.
            epsilon_b = absorption[ki] / total[ki] * thermal[ki, :, None]
            coefficient = (1.5 * KAPPA_T * scat_cyclic[ki, :, None, :]
                           * self.q[None, :, :] / total[ki, :, :, None])
            epsilon_b[-1, :] = thermal[ki, -1]
            coefficient[-1, :, :] = 0
            moment_matrix = np.eye(3 * nd)
            moment_rhs = np.zeros((3, nd))
            weighted_q = self.q * self.channel_w[:, None]
            for a in range(nc):
                direct = lambda_u[a] @ epsilon_b[:, a]
                for beta in range(3):
                    moment_rhs[beta] += weighted_q[a, beta] * direct
                    for alpha in range(3):
                        moment_matrix[beta*nd:(beta+1)*nd, alpha*nd:(alpha+1)*nd] -= (
                            weighted_q[a, beta] * lambda_u[a] * coefficient[:, a, alpha][None, :]
                        )
            from scipy.linalg import solve
            moments = solve(moment_matrix, moment_rhs.ravel(), check_finite=False).reshape(3, nd).T
            source = epsilon_b + np.einsum('dca,da->dc', coefficient, moments)
            solved = np.stack([lambda_u[a] @ source[:, a] for a in range(nc)], axis=-1)
            all_u[ki] = solved
            outgoing = np.stack([lambda_plus[a] @ source[:, a] for a in range(nc)], axis=-1)
            incoming = np.stack([lambda_minus[a] @ source[:, a] for a in range(nc)], axis=-1)
            surface_i[ki] = outgoing[0]
            flux += 2 * np.pi * np.sum(self.channel_w[None, :] * self.channel_mu[None, :] * (outgoing - incoming), axis=-1) * self.dnu[ki]
            aw = absorption[ki] * self.channel_w[None, :]
            heating += np.sum(aw * (solved - thermal[ki, :, None]), axis=-1) * self.dnu[ki]
            emitted += np.sum(aw * thermal[ki, :, None], axis=-1) * self.dnu[ki]
        output = dict(rho=rho, flux=flux, heating=heating, emitted=emitted,
                      residual=heating / np.maximum(emitted, 1e-100), absorption=absorption, total=total)
        if keep_intensity:
            output['u'] = all_u
            output['surface_i_nu_cgs'] = surface_i
        return output

    def residual(self, log_temperature):
        import time
        self.calls += 1
        t = np.exp(log_temperature)
        out = self.transfer(t)
        target = SIGMA * self.p.teff_k**4
        flux_residual = (out['flux'] - target) / target
        # The flux hardly constrains T in the optically thin upper layers.
        # There use the LTE integral energy equation; deep in the atmosphere
        # local J-B suffers cancellation, so use the constant-flux equation.
        residual = np.where(self.mass < 3.0, out['residual'], flux_residual)
        if self.calls % (self.p.depth_points + 1) == 1:
            state = dict(calls=self.calls, seconds=time.monotonic() - self.started,
                         max_flux_error=float(np.max(np.abs(flux_residual))),
                         max_local_relative=float(np.max(np.abs(out['residual'][:-1]))),
                         surface_flux_ratio=float(out['flux'][0] / target),
                         min_t=float(t.min()), max_t=float(t.max()))
            self.history.append(state)
            print(json.dumps(state), flush=True)
            if self.checkpoint:
                np.savez_compressed(self.checkpoint, mass_column_g_cm2=self.mass, temperature_k=t)
        return residual

    def solve(self, initial_temperature=None):
        if initial_temperature is None:
            initial_temperature = self.p.teff_k * (0.75 * (0.05 + 0.004 * self.mass))**0.25
        result = least_squares(self.residual, np.log(initial_temperature),
                               bounds=(np.log(self.p.teff_k * 0.03), np.log(self.p.teff_k * 100)),
                               ftol=2e-9, xtol=2e-9, gtol=2e-9,
                               diff_step=2e-5, max_nfev=self.p.max_nfev)
        return result, self.transfer(np.exp(result.x), True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--depths', type=int, default=24)
    parser.add_argument('--energies', type=int, default=48)
    parser.add_argument('--optical-nodes', type=int, default=30)
    parser.add_argument('--angles', type=int, default=3)
    parser.add_argument('--column-max', type=float, default=1e4)
    parser.add_argument('--max-nfev', type=int, default=80)
    parser.add_argument('--initial', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    p = Parameters(depth_points=args.depths, energy_points=args.energies,
                   angles_per_hemisphere=args.angles, max_nfev=args.max_nfev, column_max=args.column_max)
    column = LinearColumn(p, args.optical_nodes)
    column.checkpoint = args.output / 'checkpoint.npz'
    initial = None
    if args.initial:
        previous = np.load(args.initial)
        initial = np.exp(np.interp(np.log(column.mass), np.log(previous['mass_column_g_cm2']), np.log(previous['temperature_k'])))
    result, output = column.solve(initial)
    # Shared export helper's convention includes an explicit boundary
    # temperature; append it only for export, not to the actual unknowns.
    from types import SimpleNamespace
    exported_result = SimpleNamespace(x=np.r_[result.x, result.x[-1]], success=result.success, message=result.message)
    report = save_run(args.output, p, column, exported_result, output)
    report['solver'] = 'linear-source integral response, nonlinear constant-flux temperature solve'
    import hashlib
    report['solver_source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report['spectral_nodes_total'] = len(column.nu)
    report['optical_nodes_added'] = args.optical_nodes
    report['no_final_spectrum_renormalization'] = True
    report['refinement_pass'] = False
    (args.output / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
