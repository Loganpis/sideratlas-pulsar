"""Explicit-angle SI export from a conservative approximate atmosphere run.

No spectrum is rescaled; no optical/RJ extrapolation is performed. All
physical-approximation and scientific-validation flags remain explicit.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
from scipy.constants import c
from conservative_magnetic_column import ConservativeColumn
from magnetic_column import Parameters, MP, KB, KAPPA_T, planck_nu, cyclic_powers, response_opacities


def emergent(column, temperature, u_mean, mu):
    tc=.5*(temperature[:-1]+temperature[1:])
    rho=MP*column.p.gravity_cgs*column.cell_mass/(2*KB*tc)
    bn=planck_nu(column.nu[:,None],temperature[None,:])/2
    ac,sc=response_opacities(column.energy[:,None],rho[None,:],tc[None,:],column.p.magnetic_field_g,column.p.gaunt_factor)
    q=cyclic_powers(mu)
    ab=np.einsum('kda,ca->kdc',ac,q)
    kernel=1.5*KAPPA_T*np.einsum('kda,ca,ba,b->kdcb',sc,q,column.q,column.channel_w)
    total=ab+kernel.sum(-1)
    eps=ab/total
    scatter=np.einsum('kdab,kdb->kda',kernel,u_mean)/total
    sm=eps*.5*(bn[:,:-1]+bn[:,1:])[...,None]+scatter
    ds=eps*np.diff(bn,axis=1)[...,None]
    positive=mu>0
    # mu=0 is the exact grazing limit of the outermost cell's source.
    answer=(eps[:,0]*bn[:,0,None]+scatter[:,0]).copy()
    channels=np.tile(positive,2)
    if np.any(channels):
        tau=total[:,:,channels]*column.width[None,:,None]/np.tile(mu,2)[None,None,channels]
        tr=np.exp(-np.minimum(tau,750));om=-np.expm1(-np.minimum(tau,750))
        phi=om/tau;right=1-phi
        tiny=tau<1e-4
        right[tiny]=tau[tiny]/2-tau[tiny]**2/6+tau[tiny]**3/24
        left=om-right
        outgoing=bn[:,-1,None]+np.diff(bn,axis=1)[:,-1,None]/tau[:,-1]
        for j in range(len(column.width)-1,-1,-1):
            outgoing=tr[:,j]*outgoing+om[:,j]*sm[:,j,channels]-.5*(right[:,j]-left[:,j])*ds[:,j,channels]
        answer[:,channels]=outgoing
    return answer.reshape((len(column.nu),2,len(mu))).sum(axis=1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    p.add_argument('--allow-unvalidated-physical-model',action='store_true')
    a=p.parse_args()
    report=json.loads((a.run/'verification.json').read_text())
    if not report['numerical_equilibrium_pass']:
        raise SystemExit('Refusing export: numerical equilibrium gate has not passed')
    if not a.allow_unvalidated_physical_model:
        raise SystemExit('The surface model is physically approximate; acknowledge --allow-unvalidated-physical-model for a disclosed research visualization')
    col=ConservativeColumn(Parameters(**report['parameters']),report['optical_nodes_added'])
    saved=np.load(a.run/'column.npz');t=saved['temperature_k'];out=col.transfer(t,True)
    mu=np.unique(np.r_[0,.0001,.001,.005,.01,.02,.04,.07,.1,np.linspace(.15,1,35)])
    i_nu=emergent(col,t,out['u'],mu)
    lam=c/col.nu;order=np.argsort(lam)
    ilam=i_nu*1e-3*c/lam[:,None]**2
    if not np.all(np.isfinite(ilam)) or np.any(ilam<0):
        raise ValueError('Nonfinite or negative emergent intensity')
    output=a.run/'surface-intensity-si.npz'
    np.savez_compressed(output,wavelength_m=lam[order],mu=mu,i_lambda_si=ilam[order],
        units='W m^-2 sr^-1 m^-1',mode_sum=True,numerical_equilibrium_pass=True,
        scientific_validation_pass=False,physical_model_approximate=True,
        temperature_map='uniform normal-field plane-parallel column repeated over surface',
        teff_k=col.p.teff_k,magnetic_field_g=col.p.magnetic_field_g,gravity_cgs=col.p.gravity_cgs)
    metadata=dict(file=output.name,sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        units='W m^-2 sr^-1 m^-1',axis_order=['wavelength_m','mu'],
        mode_sum='ordinary plus extraordinary',wavelength_min_m=float(lam.min()),wavelength_max_m=float(lam.max()),
        explicit_local_optical_interval_nm=[300,650],optical_nodes_added=report['optical_nodes_added'],
        mu_min=0,mu_max=1,mu_zero='analytic grazing limit of the outermost calculated cell',
        numerical_equilibrium_pass=True,scientific_validation_pass=False,physical_model_approximate=True,
        grid_refinement_pass=report.get('refinement_pass',False),
        no_final_spectrum_renormalization=True,
        assumptions=['fully ionized ideal-gas pure H','unity magnetic Gaunt factors',
                     'transverse linear X/O polarization','coherent mode-changing scattering',
                     'plane-parallel static column with normal homogeneous B',
                     'no partial ionization or atomic opacity','no refraction/dispersion or vacuum conversion',
                     'no conduction, irradiation or coupled return-current heating',
                     'independent magnetic optical spectrum benchmark absent; physical errors unbounded'])
    (a.run/'surface-intensity-si.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata,indent=2))


if __name__=='__main__':main()
