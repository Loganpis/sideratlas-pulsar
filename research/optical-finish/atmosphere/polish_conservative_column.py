"""Bounded numerical polish of an existing conservative atmosphere column.

Only the relative weighting of the same physical residuals changes. The
emergent spectrum is never rescaled. Stop when independently reported flux,
local balance and conservation thresholds have actually been reached.
"""
import argparse,json,time,hashlib
from pathlib import Path
from dataclasses import asdict
import numpy as np
from scipy.optimize import least_squares
from conservative_magnetic_column import ConservativeColumn
from magnetic_column import Parameters,SIGMA


class NumericalGateReached(Exception):
    pass


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--depths',type=int,default=64);p.add_argument('--max-nfev',type=int,default=20)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    params=Parameters(depth_points=a.depths,energy_points=64,angles_per_hemisphere=4,column_max=1e4,max_nfev=a.max_nfev)
    col=ConservativeColumn(params,40);x=np.load(a.checkpoint);initial=x['temperature_k']
    target=SIGMA*params.teff_k**4;calls=0;history=[];best=None;best_cost=np.inf;gate=False
    def residual(logt):
        nonlocal calls,best,best_cost,gate
        calls+=1;t=np.exp(logt);out=col.transfer(t)
        fr=out['flux']/target-1;qr=out['residual'];r=np.r_[fr,100*qr]
        cost=float(r@r)
        if cost<best_cost:
            best_cost=cost;best=t.copy()
            np.savez_compressed(a.output/'checkpoint.npz',mass_column_g_cm2=col.mass,temperature_k=t)
            state=dict(calls=calls,seconds=time.monotonic()-col.started,cost=cost,
                max_flux_error=float(np.max(abs(fr))),max_local_relative=float(np.max(abs(qr))))
            history.append(state);print(json.dumps(state),flush=True)
            if np.max(abs(fr))<.002 and np.max(abs(qr))<5e-5:
                gate=True;raise NumericalGateReached
        return r
    try:
        result=least_squares(residual,np.log(initial),
            bounds=(np.log(initial)-.25,np.log(initial)+.25),diff_step=1e-5,
            ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=a.max_nfev,x_scale='jac')
        stop=str(result.message)
    except NumericalGateReached:
        stop='Explicit numerical equation gates reached; not an optimizer gtol claim'
    t=best;out=col.transfer(t,True);fr=out['flux']/target-1;qr=out['residual']
    ledger=(np.diff(out['flux'])-4*np.pi*out['heating']*col.width)/target
    accepted=bool(np.max(abs(fr))<.01 and np.max(abs(qr))<1e-4 and np.max(abs(ledger))<1e-4 and np.all(out['surface_i_nu_cgs']>=0))
    rep=dict(parameters=asdict(params),solver='bounded conservative numerical polish, local residual weight100',
        solver_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_checkpoint_sha256=hashlib.sha256(a.checkpoint.read_bytes()).hexdigest(),
        optical_nodes_added=40,spectral_nodes_total=len(col.nu),runtime_s=time.monotonic()-col.started,
        residual_evaluations=calls,stop_reason=stop,max_flux_error=float(np.max(abs(fr))),
        max_local_relative_balance=float(np.max(abs(qr))),max_conservative_ledger_error_in_target_flux=float(np.max(abs(ledger))),
        temperature_monotonic=bool(np.all(np.diff(t)>=0)),max_negative_log_temperature_increment=float(max(0,-np.min(np.diff(np.log(t))))),
        numerical_equilibrium_pass=accepted,refinement_pass=False,scientific_validation_pass=False,
        production_optical_surface_ready=False,no_final_spectrum_renormalization=True,history=history)
    np.savez_compressed(a.output/'column.npz',mass_column_g_cm2=col.mass,cell_mass_g_cm2=col.cell_mass,
        temperature_k=t,density_g_cm3=out['rho'],flux_cgs=out['flux'],local_relative_heating=qr,
        frequency_hz=col.nu,d_frequency_hz=col.dnu,mu=col.mu,mu_weights=col.w,
        surface_i_nu_cgs=out['surface_i_nu_cgs'].reshape((-1,2,len(col.mu))),
        absorption_cgs=out['absorption'],total_opacity_cgs=out['total'])
    (a.output/'verification.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep,indent=2),flush=True)


if __name__=='__main__':main()
