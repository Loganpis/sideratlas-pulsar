"""Independent boundary/energy diagnostics for the conservative column files.

Thermalization depth estimates are diagnostics, not proof of LTE accuracy or
an opacity error bound. Refining the bottom boundary remains required.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from scipy.constants import c,h,k


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);a=p.parse_args()
    rep=json.loads((a.run/'verification.json').read_text())
    x=np.load(a.run/'column.npz');m=x['mass_column_g_cm2'];dm=np.diff(m)
    t=x['temperature_k'];rho=x['density_g_cm3'];nu=x['frequency_hz'];dnu=x['d_frequency_hz']
    mu=x['mu'];w=x['mu_weights'];ab=x['absorption_cgs'];total=x['total_opacity_cgs']
    iv=x['surface_i_nu_cgs'];fnu=2*np.pi*np.sum(iv*np.tile(w*mu,2).reshape(1,2,len(mu)),axis=(1,2))
    weights=fnu*dnu;target=5.670374419e-5*rep['parameters']['teff_k']**4
    tau=np.sum(total*dm[None,:,None],axis=1)
    taueff=np.cumsum(np.sqrt(3*ab*total)*dm[None,:,None],axis=1)
    planck=lambda tt:2*h*1e7*nu**3/(c*100)**2/np.expm1(np.minimum(h*nu/(k*tt),700))/2
    bb=planck(t[-1]);deltab=bb-planck(t[-2])
    bottom_tau=total[:,-1]*dm[-1]/np.tile(mu,2)[None,:]
    diffusion_correction=np.abs(deltab[:,None]/bottom_tau/np.maximum(bb[:,None],1e-300))
    lam=c/nu
    optical=(lam>=300e-9)&(lam<=650e-9)
    support=weights>np.max(weights)*1e-7
    report=dict(input_sha256=hashlib.sha256((a.run/'column.npz').read_bytes()).hexdigest(),
        integrated_surface_flux_ratio=float(np.sum(weights)/target),
        optical_band_emitted_nm=[300,650],
        optical_bottom_total_tau_min=float(np.min(tau[optical])),
        optical_bottom_sqrt3_abs_total_tau_min=float(np.min(taueff[optical,-1])),
        supported_bolometric_bottom_sqrt3_abs_total_tau_min=float(np.min(taueff[support,-1])),
        optical_bottom_diffusion_intensity_correction_max=float(np.max(diffusion_correction[optical])),
        supported_bolometric_bottom_diffusion_intensity_correction_max=float(np.max(diffusion_correction[support])),
        zero_irradiation_boundary=True,
        lower_boundary_extension_test_pass=False,
        notes=['sqrt(3 kappa_abs kappa_tot) depth is a diagnostic in an anisotropic two-mode medium, not an exact independent solution',
               'Dense-cell pressure and opacity use the arithmetic mean of nodal temperatures',
               'Neutral and dispersive microphysics errors remain unbounded'])
    # Give the EOS audit consistent cell-centred arrays; do not misalign
    # nodal temperature with cell density.
    np.savez_compressed(a.run/'eos-cell-centres.npz',mass_column_g_cm2=.5*(m[:-1]+m[1:]),
        temperature_k=.5*(t[:-1]+t[1:]),density_g_cm3=rho)
    rows=[]
    for ki in np.flatnonzero(optical):
        mode_rows=[]
        for mode in range(2):
            channel=mode*len(mu)+int(np.argmin(abs(mu-.5)))
            j=int(np.searchsorted(taueff[ki,:,channel],1))
            j=min(j,len(rho)-1)
            ep_ev=28.71*np.sqrt(rho[j]);ph_ev=h*nu[ki]/1.602176634e-19
            mode_rows.append(dict(mode=mode,cell_index=j,density_g_cm3=float(rho[j]),
                temperature_k=float(.5*(t[j]+t[j+1])),plasma_energy_over_photon_energy=float(ep_ev/ph_ev)))
        rows.append(dict(wavelength_nm=float(lam[ki]*1e9),formation_diagnostics=mode_rows))
    report['optical_formation_rows']=rows
    (a.run/'boundary-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='optical_formation_rows'},indent=2))


if __name__=='__main__':main()
