#!/usr/bin/env python3
"""Real response / viewing sensitivity and independent radiation diagnostics.

This calculation explicitly uses the FROZEN LEGACY trajectories to isolate
instrument and smoothing effects. It does not claim that the new, not yet
converged numerical plasma field has been validated for particle transport.
"""
from pathlib import Path
import sys,os,json,hashlib
from dataclasses import replace,asdict
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pulsar.trajectories import Geometry,trace_bundle,sky_profiles
from pulsar.radiation import transport_curvature,MEC2_ERG,E_STATC,curvature_power_erg_s
from pulsar.response import load_lat_response,CurvatureResponse,incidence_scenarios
from pulsar.multiband import load_bands,score_profiles as _historical_score_profiles
from pulsar.injection import seed_injection_weights
from pulsar.nonthermal import (synchrotron_power_erg_s,isotropized_stellar_photon_field,
                              inverse_compton_rates)

OUT=ROOT/'outputs'/'emission-upgrade'
BANDS=((.1,.3),(.3,1.),(1.,3.),(3.,100.))


def score_profiles(*args,**kwargs):
    """Keep historical fitter mathematics, replace obsolete claim metadata."""
    scores=_historical_score_profiles(*args,**kwargs)
    for score in scores:
        score['role']='development'
        score['statistical_caveat']=(
            'Conditional shape-only diagnostic: independent nonnegative scale '
            'and constant per band; historical phase origin held fixed. '
            'All bands are development data. Response is explicit per '
            'scenario: flat, or calibrated LAT area and energy dispersion '
            'with assumed incidence. Source-specific exposure, spatial '
            'acceptance and catalog source weights are absent. No absolute '
            'flux fit, independent validation or calibrated significance.')
    return scores


def dump(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else float(x))+'\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    frozen=json.loads((ROOT/'outputs/particle-model/frozen-model.json').read_text())
    g=replace(Geometry(**frozen['geometry']),phase_sigma=0.)
    rate=frozen['acceleration_rate_cm_inverse']; phase_shift=frozen['phase_shift_cycles']
    plan={'legacy_geometry':asdict(g),'legacy_acceleration_cm_inverse':rate,
          'phase_shift_cycles':phase_shift,'phase_shift_refitted':False,
          'phase_blur_removed':True,'all_profiles_status':'development data; no new independent validation',
          'sky_bin':'0.5 degrees numerical Gaussian; compare 1 degree',
          'viewing_scenarios_deg':[62.,63.6,65.2],
          'viewing_scenarios_are_confidence_interval':False,
          'injection_gamma_scenarios':[1000.,100000.,1000000.],
          'primary_loading':'corotation_charge: |div E_corotation| |v_radial| dA at true seeds',
          'loading_limit':'Prescribed charge-loading hypothesis; pair supply and current closure are not solved',
          'loading_sensitivity':'Compare with legacy uniform magnetic-flux loading; no tuning to data',
          'instrument':'Official P8R3_SOURCE_V3 FRONT/BACK effective area and energy dispersion',
          'incidence_scenarios':'uniform cosine, near-axis, off-axis; none is measured exposure',
          'plasma_field_used_for_this_emission':'NO: frozen legacy field isolates response/smoothing effects'}
    dump(OUT/'protocol.json',plan)
    bundle=trace_bundle(g,n_azimuth=512,n_rings=7,dt=.0015)
    uniform_bundle=bundle.copy()
    bundle['seed_weight']=seed_injection_weights(bundle,g,mode='corotation_charge')
    loading_ratio=bundle['seed_weight']/uniform_bundle['seed_weight']
    work=rate*MEC2_ERG/E_STATC*bundle['projection']
    t=transport_curvature(bundle['dt']*g.rlc_cm,bundle['rho']*g.rlc_cm,work,1000.,BANDS)
    phase,flat=sky_profiles(bundle,t.photon_yield,g,nphase=1024,view_sigma_deg=.5)
    data=load_bands(ROOT/'data/J0835-4510_3PC_data.fits',rebin=8)
    for band in data:
        band['role']='development'
    response=load_lat_response(); profiles={}; scores={}
    folds={}
    for name,(mu,w) in incidence_scenarios().items():
        fold=CurvatureResponse(response,BANDS,mu,w); folds[name]=fold
        weighted=fold.photon_yield(t.radiated_erg,t.critical_energy_gev)
        _,profile=sky_profiles(bundle,weighted,g,nphase=1024,view_sigma_deg=.5)
        profiles[name]=profile
        scores[name]=score_profiles(phase,profile,data,phase_shift)
        print('Response scenario completed: '+name,flush=True)
        if name=='uniform_cosine': primary_yield=weighted
    _,view_wide=sky_profiles(bundle,primary_yield,g,nphase=1024,view_sigma_deg=1.)
    _,uniform_loading=sky_profiles(uniform_bundle,primary_yield,g,nphase=1024,view_sigma_deg=.5)
    loading_scores=score_profiles(phase,uniform_loading,data,phase_shift)
    viewing={}
    for angle in plan['viewing_scenarios_deg']:
        _,profile=sky_profiles(bundle,primary_yield,replace(g,observer_deg=angle),nphase=1024,view_sigma_deg=.5)
        viewing[str(angle)]=score_profiles(phase,profile,data,phase_shift)
    # All gamma variants use the same trajectory and detector assumptions.
    injection={}
    for initial in plan['injection_gamma_scenarios']:
        if initial==1000: other=t
        else: other=transport_curvature(bundle['dt']*g.rlc_cm,bundle['rho']*g.rlc_cm,work,initial)
        y=folds['uniform_cosine'].photon_yield(other.radiated_erg,other.critical_energy_gev)
        _,p=sky_profiles(bundle,y,g,nphase=1024,view_sigma_deg=.5)
        injection[str(initial)]=np.mean(np.abs(p-profiles['uniform_cosine']),axis=1)
        if other is not t: del other
    comparisons={'flat':score_profiles(phase,flat,data,phase_shift),**scores}
    plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':11})
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
    for j,(ax,d) in enumerate(zip(axes.ravel(),data)):
        norm=d['counts'].max()
        ax.errorbar(d['phase'],d['counts']/norm,d['error']/norm,fmt='.',color='#465b65',ms=4,label='Vela 3PC (development data)')
        ax.plot(d['phase'],comparisons['flat'][j]['prediction']/norm,color='#888888',ls='--',label='Flat response')
        ax.plot(d['phase'],scores['uniform_cosine'][j]['prediction']/norm,color='#aa351f',label='LAT area + energy migration')
        ax.plot(d['phase'],loading_scores[j]['prediction']/norm,color='#468f86',ls=':',label='LAT + legacy uniform-flux loading')
        alternatives=np.array([scores[k][j]['prediction']/norm for k in scores])
        ax.fill_between(d['phase'],alternatives.min(axis=0),alternatives.max(axis=0),color='#d89554',alpha=.3,label='Incidence scenarios')
        ax.set(title=d['band'],xlim=(0,1),ylim=(0,1.07),ylabel='Relative weighted counts')
        ax.spines[['top','right']].set_visible(False)
    axes[0,0].legend(fontsize=8); axes[1,0].set_xlabel('Rotation phase'); axes[1,1].set_xlabel('Rotation phase')
    fig.suptitle('Real detector response, without artificial phase blur',fontsize=18)
    fig.text(.5,.015,'Frozen legacy emission geometry · no phase refit · independent scale/constant per band · exposure and catalog weights remain unavailable',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.04,1,.95)); fig.savefig(OUT/'response-comparison.png',dpi=180); fig.savefig(OUT/'response-comparison.pdf'); plt.close(fig)
    # Evaluate independent physical channels at explicit illustrative inputs.
    # They are NOT mixed into the Vela pulse prediction with arbitrary weights.
    gamma=np.geomspace(1.1e4,1e8,100)
    radiation_inputs={'temperature_K':7e5,'stellar_radius_cm':1.2e6,'distance_cm':g.rlc_cm,
                      'magnetic_frame_B_gauss':5e4,'pitch_radians':1e-3,'curvature_radius_cm':g.rlc_cm,
                      'status':'illustrative conditions, not inferred Vela particle/seed distributions'}
    seeds=isotropized_stellar_photon_field(7e5,1.2e6,g.rlc_cm)
    ic=inverse_compton_rates(gamma,seeds,BANDS)
    sr=synchrotron_power_erg_s(gamma,5e4,1e-3)
    cr=curvature_power_erg_s(gamma,g.rlc_cm)
    fig,ax=plt.subplots(figsize=(9,5.5))
    ax.loglog(gamma,cr,label='Curvature; specified trajectory radius',color='#237b93')
    ax.loglog(gamma,sr,label='Synchrotron; specified magnetic-frame pitch',color='#a93723')
    ax.loglog(gamma,ic.electron_loss_erg_s,label='Inverse Compton; isotropized thermal field, KN',color='#836895')
    ax.set(xlabel=r'Particle Lorentz factor $\gamma$',ylabel='Single-particle energy loss [erg/s]',title='Calculated radiation channels — illustrative conditions')
    ax.legend(fontsize=10); ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.015,'B = 50 kG · pitch = 0.001 rad · T = 0.7 MK · R* = 12 km · r = Rc = 4246 km · not a fitted Vela spectrum',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.045,1,1)); fig.savefig(OUT/'radiation-channels.png',dpi=180); fig.savefig(OUT/'radiation-channels.pdf'); plt.close(fig)
    report={'protocol':plan,'scores':comparisons,'viewing_angle_sensitivity':viewing,
            'sky_bin_1_to_half_degree_L1_by_band':np.mean(np.abs(view_wide-profiles['uniform_cosine']),axis=1),
            'injection_gamma_L1_by_band':injection,'radiation_inputs':radiation_inputs,
            'uniform_flux_loading_scores':loading_scores,
            'loading_L1_by_band':np.mean(np.abs(uniform_loading-profiles['uniform_cosine']),axis=1),
            'seed_weight_relative_to_uniform_flux_range':[loading_ratio.min(),loading_ratio.max()],
            'new_channels_added_to_Vela_profile':False,
            'reason':'Particle pitch, pair distribution and anisotropic seed field are not solved; adding an arbitrary mixture would not be physical',
            'validated_for_Vela':False,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                             [ROOT/'run_emission_upgrade.py',*sorted((ROOT/'src/pulsar').glob('*.py')),
                              *sorted((ROOT/'data/response').glob('*.fits'))]}}
    dump(OUT/'validation.json',report)
    np.savez_compressed(OUT/'profiles.npz',phase=phase,flat=flat,uniform_flux_loading=uniform_loading,**profiles)
    np.savez_compressed(OUT/'radiation-channels.npz',gamma=gamma,curvature_power=cr,synchrotron_power=sr,IC_power=ic.electron_loss_erg_s)
    print('Completed detector, viewing, injection and radiation diagnostics',flush=True)


if __name__=='__main__': main()
