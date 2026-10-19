"""Conservative linear-LTE-source / cell-mean-scattering magnetic column.

Research model only: retains fully-ionized ideal H, unity magnetic Gaunt
factors and transverse linear polarizations. No independent magnetic optical
microphysics validation is claimed. The discretization makes total face flux
jumps equal integrated cell heating, with no output renormalization.
"""
from __future__ import annotations
import argparse, hashlib, json, time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.optimize import least_squares
from linear_magnetic_column import LinearColumn
from magnetic_column import Parameters, MP, KB, KAPPA_T, SIGMA, planck_nu, response_opacities


class ConservativeColumn(LinearColumn):
    def __init__(self, p, optical_nodes=40):
        super().__init__(p, optical_nodes)
        self.cell_mass = .5 * (self.mass[:-1] + self.mass[1:])
        self.width = np.diff(self.mass)
        self.best_cost = np.inf

    def transfer(self, temperature, keep_intensity=False, **unused):
        nd = len(self.mass) - 1
        nc = len(self.channel_mu)
        tc = .5 * (temperature[:-1] + temperature[1:])
        rho = MP * self.p.gravity_cgs * self.cell_mass / (2 * KB * tc)
        bnode = planck_nu(self.nu[:, None], temperature[None, :]) / 2
        bmean = .5 * (bnode[:, :-1] + bnode[:, 1:])
        ac, sc = response_opacities(self.energy[:, None], rho[None, :], tc[None, :], self.p.magnetic_field_g, self.p.gaunt_factor)
        absorption = np.einsum('kda,ca->kdc', ac, self.q)
        kernel = 1.5 * KAPPA_T * np.einsum('kda,ca,ba,b->kdcb', sc, self.q, self.q, self.channel_w)
        total = absorption + kernel.sum(-1)
        flux = np.zeros(nd + 1)
        heating = np.zeros(nd)
        emitted = np.zeros(nd)
        all_u = np.empty((len(self.nu), nd, nc))
        surface_i = np.empty((len(self.nu), nc))
        for ki in range(len(self.nu)):
            tau = total[ki] * self.width[:, None] / self.channel_mu[None, :]
            tr = np.exp(-np.minimum(tau, 750))
            om = -np.expm1(-np.minimum(tau, 750))
            phi = om / tau
            right = 1 - phi
            tiny = tau < 1e-4
            right[tiny] = tau[tiny]/2 - tau[tiny]**2/6 + tau[tiny]**3/24
            left = om - right
            epsilon = absorption[ki] / total[ki]
            delta_s = epsilon * np.diff(bnode[ki])[:, None]
            thermal_mean = epsilon * bmean[ki, :, None]
            coeff = 1.5 * KAPPA_T * sc[ki, :, None, :] * self.q[None, :, :] / total[ki, :, :, None]
            bottom = bnode[ki, -1] + (bnode[ki, -1] - bnode[ki, -2]) / tau[-1]
            lam, direct, up_op, down_op, up_direct, down_direct = [], [], [], [], [], []
            for a in range(nc):
                down = np.zeros((nd + 1, nd))
                up = np.zeros_like(down)
                gd = np.zeros(nd + 1)
                gu = np.zeros(nd + 1)
                gu[-1] = bottom[a]
                for j in range(nd):
                    down[j+1] = tr[j,a] * down[j]
                    down[j+1,j] += om[j,a]
                    gd[j+1] = tr[j,a]*gd[j] + .5*(right[j,a]-left[j,a])*delta_s[j,a]
                for j in range(nd - 1, -1, -1):
                    up[j] = tr[j,a] * up[j+1]
                    up[j,j] += om[j,a]
                    gu[j] = tr[j,a]*gu[j+1] - .5*(right[j,a]-left[j,a])*delta_s[j,a]
                response = .5 * phi[:, a, None] * (down[:-1] + up[1:])
                response[np.arange(nd), np.arange(nd)] += right[:, a]
                lam.append(response)
                direct.append(.5 * phi[:,a] * (gd[:-1]+gu[1:]))
                up_op.append(up); down_op.append(down)
                up_direct.append(gu); down_direct.append(gd)
            mat = np.eye(3 * nd)
            rhs = np.zeros((3, nd))
            weighted_q = self.q * self.channel_w[:, None]
            for a in range(nc):
                initial = lam[a] @ thermal_mean[:,a] + direct[a]
                for beta in range(3):
                    rhs[beta] += weighted_q[a,beta] * initial
                    for alpha in range(3):
                        mat[beta*nd:(beta+1)*nd,alpha*nd:(alpha+1)*nd] -= weighted_q[a,beta] * lam[a] * coeff[:,a,alpha][None,:]
            moments = solve(mat, rhs.ravel(), check_finite=False).reshape(3, nd).T
            source_mean = thermal_mean + np.einsum('dca,da->dc', coeff, moments)
            u = np.stack([lam[a] @ source_mean[:,a] + direct[a] for a in range(nc)], axis=-1)
            outgoing = np.stack([up_op[a] @ source_mean[:,a] + up_direct[a] for a in range(nc)], axis=-1)
            incoming = np.stack([down_op[a] @ source_mean[:,a] + down_direct[a] for a in range(nc)], axis=-1)
            all_u[ki] = u
            surface_i[ki] = outgoing[0]
            flux += 2*np.pi*np.sum(self.channel_w[None,:]*self.channel_mu[None,:]*(outgoing-incoming),axis=-1)*self.dnu[ki]
            aw = absorption[ki]*self.channel_w[None,:]
            heating += np.sum(aw*(u-bmean[ki,:,None]),axis=-1)*self.dnu[ki]
            emitted += np.sum(aw*bmean[ki,:,None],axis=-1)*self.dnu[ki]
        out = dict(rho=rho, flux=flux, heating=heating, emitted=emitted,
                   residual=heating/np.maximum(emitted,1e-200), absorption=absorption,total=total)
        if keep_intensity:
            out.update(u=all_u,surface_i_nu_cgs=surface_i)
        return out

    def residual(self, log_t):
        self.calls += 1
        t = np.exp(log_t)
        o = self.transfer(t)
        target = SIGMA*self.p.teff_k**4
        fr = o['flux']/target-1
        # Flux equations and local equations are mathematically consistent
        # because the same exact cell average is used in transfer and heating.
        qr = o['residual']
        r = np.r_[fr,qr]
        cost = float(r@r)
        if self.calls % (len(t)+1)==1:
            state=dict(calls=self.calls,seconds=time.monotonic()-self.started,cost=cost,
                       max_flux_error=float(np.max(abs(fr))),max_local_relative=float(np.max(abs(qr))),
                       surface_flux_ratio=float(o['flux'][0]/target),min_t=float(t.min()),max_t=float(t.max()))
            self.history.append(state);print(json.dumps(state),flush=True)
        if cost<self.best_cost:
            self.best_cost=cost
            if self.checkpoint:
                np.savez_compressed(self.checkpoint,mass_column_g_cm2=self.mass,temperature_k=t)
        return r

    def solve(self, initial_temperature=None):
        if initial_temperature is None:
            initial_temperature=self.p.teff_k*(.75*(.05+.004*self.mass))**.25
        result=least_squares(self.residual,np.log(initial_temperature),
            bounds=(np.log(self.p.teff_k*.03),np.log(self.p.teff_k*100)),
            ftol=1e-10,xtol=1e-10,gtol=1e-10,diff_step=1e-4,
            max_nfev=self.p.max_nfev,x_scale='jac')
        return result,self.transfer(np.exp(result.x),True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--depths',type=int,default=32)
    p.add_argument('--energies',type=int,default=64)
    p.add_argument('--angles',type=int,default=4)
    p.add_argument('--optical-nodes',type=int,default=40)
    p.add_argument('--column-min',type=float,default=1e-9)
    p.add_argument('--column-max',type=float,default=1e4)
    p.add_argument('--max-nfev',type=int,default=100)
    p.add_argument('--initial',type=Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    params=Parameters(depth_points=a.depths,energy_points=a.energies,angles_per_hemisphere=a.angles,column_min=a.column_min,column_max=a.column_max,max_nfev=a.max_nfev)
    col=ConservativeColumn(params,a.optical_nodes);col.checkpoint=a.output/'checkpoint.npz'
    initial=None
    if a.initial:
        x=np.load(a.initial);initial=np.exp(np.interp(np.log(col.mass),np.log(x['mass_column_g_cm2']),np.log(x['temperature_k'])))
    result,o=col.solve(initial);t=np.exp(result.x);target=SIGMA*params.teff_k**4
    fr=o['flux']/target-1;qr=o['residual']
    ledger=(np.diff(o['flux'])-4*np.pi*o['heating']*col.width)/target
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),parameters=asdict(params),
        solver='conservative cell-mean scattering, linear LTE source, simultaneous flux/local balance',
        solver_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime_s=time.monotonic()-col.started,residual_evaluations=col.calls,
        optimizer_success=bool(result.success),optimizer_message=result.message,optical_nodes_added=a.optical_nodes,
        spectral_nodes_total=len(col.nu),max_flux_error=float(np.max(abs(fr))),max_local_relative_balance=float(np.max(abs(qr))),
        max_conservative_ledger_error_in_target_flux=float(np.max(abs(ledger))),
        numerical_equilibrium_pass=bool(result.success and np.max(abs(fr))<.01 and np.max(abs(qr))<1e-4),
        scientific_validation_pass=False,production_optical_surface_ready=False,refinement_pass=False,
        no_final_spectrum_renormalization=True,history=col.history)
    np.savez_compressed(a.output/'column.npz',mass_column_g_cm2=col.mass,cell_mass_g_cm2=col.cell_mass,
        temperature_k=t,density_g_cm3=o['rho'],flux_cgs=o['flux'],local_relative_heating=qr,
        frequency_hz=col.nu,d_frequency_hz=col.dnu,mu=col.mu,mu_weights=col.w,
        surface_i_nu_cgs=o['surface_i_nu_cgs'].reshape((-1,2,len(col.mu))),
        absorption_cgs=o['absorption'],total_opacity_cgs=o['total'])
    (a.output/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
