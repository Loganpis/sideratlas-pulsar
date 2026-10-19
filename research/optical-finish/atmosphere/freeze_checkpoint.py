"""Freeze and independently reevaluate a conservative run's best checkpoint.

This does not claim optimizer termination. Physical validity and grid tests
remain separate. Useful for a costly ongoing refinement that already meets
the stated numerical equation gates.
"""
from pathlib import Path
import argparse,hashlib,json,shutil
from dataclasses import asdict
import numpy as np
from conservative_magnetic_column import ConservativeColumn
from magnetic_column import Parameters,SIGMA


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--depths',type=int,required=True);p.add_argument('--energies',type=int,default=64)
    p.add_argument('--angles',type=int,default=4);p.add_argument('--optical-nodes',type=int,default=40)
    p.add_argument('--column-min',type=float,default=1e-9);p.add_argument('--column-max',type=float,default=1e4)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(a.checkpoint,a.output/'input-checkpoint.npz')
    z=np.load(a.output/'input-checkpoint.npz');t=z['temperature_k']
    params=Parameters(depth_points=a.depths,energy_points=a.energies,angles_per_hemisphere=a.angles,column_min=a.column_min,column_max=a.column_max)
    col=ConservativeColumn(params,a.optical_nodes)
    np.testing.assert_allclose(z['mass_column_g_cm2'],col.mass,rtol=1e-13,atol=0)
    out=col.transfer(t,True);target=SIGMA*params.teff_k**4
    fr=out['flux']/target-1;qr=out['residual'];ledger=(np.diff(out['flux'])-4*np.pi*out['heating']*col.width)/target
    accepted=bool(np.max(abs(fr))<.01 and np.max(abs(qr))<1e-4 and np.max(abs(ledger))<1e-4 and np.all(out['surface_i_nu_cgs']>=0))
    rep=dict(parameters=asdict(params),solver='conservative checkpoint independently reevaluated',
        optical_nodes_added=a.optical_nodes,spectral_nodes_total=len(col.nu),
        optimizer_status='checkpoint only; no optimizer convergence claimed',
        max_flux_error=float(np.max(abs(fr))),max_local_relative_balance=float(np.max(abs(qr))),
        max_conservative_ledger_error_in_target_flux=float(np.max(abs(ledger))),
        temperature_monotonic=bool(np.all(np.diff(t)>=0)),max_negative_log_temperature_increment=float(max(0,-np.min(np.diff(np.log(t))))),
        numerical_equilibrium_pass=accepted,refinement_pass=False,scientific_validation_pass=False,
        production_optical_surface_ready=False,no_final_spectrum_renormalization=True,
        input_checkpoint_sha256=hashlib.sha256((a.output/'input-checkpoint.npz').read_bytes()).hexdigest())
    np.savez_compressed(a.output/'column.npz',mass_column_g_cm2=col.mass,cell_mass_g_cm2=col.cell_mass,
        temperature_k=t,density_g_cm3=out['rho'],flux_cgs=out['flux'],local_relative_heating=qr,
        frequency_hz=col.nu,d_frequency_hz=col.dnu,mu=col.mu,mu_weights=col.w,
        surface_i_nu_cgs=out['surface_i_nu_cgs'].reshape((-1,2,len(col.mu))),
        absorption_cgs=out['absorption'],total_opacity_cgs=out['total'])
    (a.output/'verification.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep,indent=2))


if __name__=='__main__':main()
