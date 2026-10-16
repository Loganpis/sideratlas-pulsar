#!/usr/bin/env python3
"""Calibrate a reduced curvature-emission hypothesis, then open reserved bands.

No synthetic pulse template or fitted coefficient per phase is used. This is
not a self-consistent PIC simulation or a claim of unique Vela reconstruction.
"""
from pathlib import Path
import os
import sys
import json
import hashlib
import gc
import platform
from dataclasses import asdict, replace
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
from pulsar.trajectories import Geometry,trace_bundle,sky_profiles,rotate_z,curvature,velocity
from pulsar.radiation import transport_curvature,MEC2_ERG,E_STATC
from pulsar.multiband import load_bands,fit_shared_phase,score_profiles,VALIDATION_CAVEAT,ERROR_CAVEAT
from pulsar.observations import bin_average

OUT=ROOT/'outputs'/'particle-model'
ENERGIES=((.1,.3),(.3,1.),(1.,3.),(3.,100.))


def native(obj):
    if isinstance(obj,np.ndarray): return obj.tolist()
    if isinstance(obj,np.generic): return obj.item()
    if isinstance(obj,Path): return str(obj)
    raise TypeError(type(obj).__name__)


def write_json(name,data):
    (OUT/name).write_text(json.dumps(data,default=native,indent=2,allow_nan=False)+'\n')


def calculate(bundle,g,rate,spectral_nodes=1):
    # rate=e |E_parallel|/(m_e c^2), cm^-1. Actual electric work has the
    # projection |v dot Bhat|; it is not e E_parallel c everywhere.
    work_field=rate*MEC2_ERG/E_STATC*bundle['projection']
    result=transport_curvature(bundle['dt']*g.rlc_cm,bundle['rho']*g.rlc_cm,
                               work_field,1000.,ENERGIES,spectral_nodes=spectral_nodes)
    phase,profiles=sky_profiles(bundle,result.photon_yield,g,nphase=1024)
    return phase,profiles,result


def shape_change(phase,a,b,data):
    changes=[]
    for first,second,d in zip(a,b,data):
        x=bin_average(phase,first,d['phase_lo'],d['phase_hi'],0)
        y=bin_average(phase,second,d['phase_lo'],d['phase_hi'],0)
        if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y)) and x.mean()>0 and y.mean()>0):
            raise ValueError('Cannot assess relative convergence of a dark or nonfinite profile')
        x/=x.mean(); y/=y.mean()
        changes.append(float(np.mean(np.abs(x-y))))
    return changes


def plot_comparison(data,scores,baseline):
    plt.rcParams.update({'font.family':'serif','mathtext.fontset':'cm','font.size':11})
    fig,axes=plt.subplots(4,2,figsize=(12,12),gridspec_kw={'width_ratios':[3,1]},sharex=True)
    for row,(d,s,old) in enumerate(zip(data,scores,baseline)):
        ax,res=axes[row]; factor=np.max(d['counts'])
        ax.errorbar(d['phase'],d['counts']/factor,d['error']/factor,fmt='.',color='#506170',ms=4,lw=.6,label='Fermi-LAT 3PC')
        ax.plot(d['phase'],old['prediction']/factor,color='#aaaaaa',lw=1.2,ls='--',label='Earlier wind benchmark')
        ax.plot(d['phase'],s['prediction']/factor,color='#b33d24',lw=1.8,label='Particle + curvature model')
        role='FITTED' if d['role']=='fit' else 'RESERVED FROM FITTING'
        ax.set_title(f"{d['band']}  ·  {role}",loc='left',fontsize=12)
        ax.set_ylabel('Relative weighted counts'); ax.set_ylim(bottom=0)
        res.axhline(0,color='#555555',lw=.7)
        res.plot(d['phase'],s['residual_sigma'],color='#b33d24',lw=1)
        res.set_title(r'$\chi^2/(N-2)$'+f" = {s['reduced_chi_squared']:.1f}",fontsize=11)
        res.set_ylabel('Residual / catalog error')
        for a in (ax,res):
            a.spines[['right','top']].set_visible(False); a.set_xlim(0,1)
    axes[0,0].legend(fontsize=9,loc='upper right')
    axes[-1,0].set_xlabel('Rotation phase'); axes[-1,1].set_xlabel('Rotation phase')
    fig.suptitle('@sideratlas  /  NO.002 — Vela emission-model test',fontsize=18,y=.992)
    fig.text(.5,.016,'Prescribed magnetic geometry · flat within-band response · normalization and constant fitted in every band',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.04,1,.97)); fig.savefig(OUT/'vela-energy-profiles.png',dpi=180)
    fig.savefig(OUT/'vela-energy-profiles.pdf'); plt.close(fig)


def plot_particles(bundle,g,transport):
    plt.rcParams.update({'font.family':'serif','mathtext.fontset':'cm','font.size':11})
    fig=plt.figure(figsize=(13,8)); ax=fig.add_subplot(121,projection='3d')
    b=bundle; colors=['#237b93','#a93723']; nseed=b['dt'].shape[1]
    picked=np.linspace(0,nseed-1,14,dtype=int)
    for j in picked:
        ok=b['dt'][:,j]>0
        # A Lagrangian path: accumulated time differs from instantaneous
        # rotating-field snapshots. Show the actual trajectory in the lab.
        times=np.cumsum(b['dt'][:,j])-.5*b['dt'][:,j]
        p=rotate_z(b['position'][:,j],times)[ok]
        ax.plot(*p.T,color=colors[int(b['pole'][j]<0)],alpha=.8,lw=1)
    ang=np.linspace(0,2*np.pi,200)
    for z in (-.7,.7): ax.plot(np.cos(ang),np.sin(ang),z,color='#cccccc',lw=.6)
    ax.scatter([0],[0],[0],color='#333333',s=20)
    ax.set(xlabel=r'$x/R_{\rm LC}$',ylabel=r'$y/R_{\rm LC}$',zlabel=r'$z/R_{\rm LC}$')
    ax.set_title('Calculated laboratory particle paths',pad=15)
    ax.set_box_aspect((1,1,1)); ax.view_init(23,-52)
    a=fig.add_subplot(222); c=fig.add_subplot(224)
    for j in picked:
        ok=b['dt'][:,j]>0; r=np.linalg.norm(b['position'][:,j],axis=-1)[ok]
        color=colors[int(b['pole'][j]<0)]
        a.plot(r,transport.gamma_mid[:,j][ok],color=color,alpha=.7,lw=.9)
        c.plot(r,transport.critical_energy_gev[:,j][ok],color=color,alpha=.7,lw=.9)
    a.set(xlabel=r'$r/R_{\rm LC}$',ylabel=r'Particle Lorentz factor $\gamma$',yscale='log')
    a.set_title('Acceleration and curvature-radiation losses')
    c.set(xlabel=r'$r/R_{\rm LC}$',ylabel=r'Critical photon energy $E_c$ [GeV]',yscale='log',ylim=(.01,100))
    c.set_title(r'$E_c=3\hbar c\gamma^3/(2R_c)$')
    for a in (a,c): a.spines[['right','top']].set_visible(False)
    fig.suptitle('Particle transport in a prescribed rotating dipole',fontsize=18)
    fig.text(.5,.015,'Two magnetic hemispheres · injection at 0.08 light-cylinder radii · no pair cascade or plasma feedback',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.04,1,.96)); fig.savefig(OUT/'particle-transport.png',dpi=180)
    fig.savefig(OUT/'particle-transport.pdf'); plt.close(fig)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    data=load_bands(ROOT/'data/J0835-4510_3PC_data.fits',rebin=8)
    plan={
        'model':'Retarded point-dipole slot-gap geometry + curvature energy transport',
        'training_band_indices':[0,2],'reserved_band_indices':[1,3],
        'inclination_grid_deg':[45.,55.,65.,75.,85.],
        'acceleration_grid_cm_inverse':[.01,.025,.06,.15,.4,1.],
        'fixed_geometry':asdict(Geometry()),'gamma_at_injection':1000.,
        'coarse_quadrature':{'n_azimuth':96,'n_rings':3,'dt':.006},
        'fine_quadrature':{'n_azimuth':256,'n_rings':5,'dt':.003},
        'refinement_quadrature':{'n_azimuth':512,'n_rings':7,'dt':.0015},
        'selection':'Grid calibration on training bands, then recheck three best candidates at fine quadrature; choose from training only',
        'final_phase':'Refit shared phase on training at refinement quadrature before opening reserved bands',
        'response':'Flat photon-number response within each energy band; independent nonnegative A,C per band',
        'convergence_target_rebinned_L1':.02,
        'heldout_used_for_model_redesign':False,
        'catalog_error_caveat':ERROR_CAVEAT,'validation_caveat':VALIDATION_CAVEAT,
    }
    resume='--resume' in sys.argv
    if resume and json.loads((OUT/'protocol.json').read_text()) != plan:
        raise ValueError('Saved protocol differs; cannot resume an incompatible calibration')
    write_json('protocol.json',plan)
    candidates=json.loads((OUT/'training-grid.json').read_text()) if resume else []
    for chi in ([] if resume else plan['inclination_grid_deg']):
        g=replace(Geometry(),inclination_deg=chi)
        bundle=trace_bundle(g,**plan['coarse_quadrature'])
        for rate in plan['acceleration_grid_cm_inverse']:
            phase,profiles,transport=calculate(bundle,g,rate)
            fit=fit_shared_phase(phase,profiles,data,n_shifts=128,evaluate_holdout=False)
            row={'inclination_deg':chi,'rate_cm_inverse':rate,'training_chi_squared':fit['training_chi_squared'],
                 'phase_shift_cycles':fit['phase_shift_cycles']}
            candidates.append(row)
            print(f"Training only: alpha={chi:g}, acceleration={rate:g}, score={row['training_chi_squared']:.1f}",flush=True)
            del transport
        del bundle; gc.collect()
    candidates.sort(key=lambda row:row['training_chi_squared'])
    write_json('training-grid.json',candidates)
    refined=json.loads((OUT/'frozen-model.json').read_text())['fine_candidates'] if resume else []
    for candidate in ([] if resume else candidates[:3]):
        g=replace(Geometry(),inclination_deg=candidate['inclination_deg'])
        bundle=trace_bundle(g,**plan['fine_quadrature'])
        phase,profiles,transport=calculate(bundle,g,candidate['rate_cm_inverse'])
        fit=fit_shared_phase(phase,profiles,data,n_shifts=256,evaluate_holdout=False)
        refined.append({**candidate,'training_chi_squared':fit['training_chi_squared'],
                        'phase_shift_cycles':fit['phase_shift_cycles']})
        print(f"Fine training check: {refined[-1]}",flush=True)
        np.savez_compressed(OUT/f"candidate-{g.inclination_deg:g}-{candidate['rate_cm_inverse']:g}.npz",phase=phase,profiles=profiles)
        del bundle,transport; gc.collect()
    best=min(refined,key=lambda row:row['training_chi_squared'])
    g=replace(Geometry(),inclination_deg=best['inclination_deg']); rate=best['rate_cm_inverse']
    bundle=trace_bundle(g,**plan['refinement_quadrature'])
    phase,profiles,transport=calculate(bundle,g,rate)
    final_fit=fit_shared_phase(phase,profiles,data,n_shifts=512,evaluate_holdout=False)
    shift=final_fit['phase_shift_cycles']
    selected={'geometry':asdict(g),'acceleration_rate_cm_inverse':rate,
              'parallel_electric_field_statvolt_cm':rate*MEC2_ERG/E_STATC,
              'gamma_at_injection':1000.,'phase_shift_cycles':shift,
              'fine_candidates':refined,'reserved_data_used_for_selection':False,
              'inclination_at_grid_boundary':g.inclination_deg in (plan['inclination_grid_deg'][0],plan['inclination_grid_deg'][-1]),
              'acceleration_at_grid_boundary':rate in (plan['acceleration_grid_cm_inverse'][0],plan['acceleration_grid_cm_inverse'][-1])}
    write_json('frozen-model.json',selected)
    print('Physical parameters and phase frozen; now evaluating reserved bands',flush=True)
    scores=score_profiles(phase,profiles,data,shift)
    # Fair baseline: the old >100 MeV power-law wind has an energy-independent
    # shape; fit its one phase from the same training bands, nuisance-only elsewhere.
    old=np.load(ROOT/'outputs/profile_32_128_256.npz')
    oldphase=old['phase'] if 'phase' in old else old['phases']
    oldcurve=old['normalized']
    oldprofiles=np.tile(oldcurve,(4,1))
    baseline_fit=fit_shared_phase(oldphase,oldprofiles,data,n_shifts=512,evaluate_holdout=False)
    baseline=score_profiles(oldphase,oldprofiles,data,baseline_fit['phase_shift_cycles'])
    middle=np.load(OUT/f'candidate-{g.inclination_deg:g}-{rate:g}.npz')['profiles']
    resolution=shape_change(phase,middle,profiles,data)
    # Sensitivity is reported, never used to redesign the frozen fit.
    _,narrow_view=sky_profiles(bundle,transport.photon_yield,g,nphase=1024,view_sigma_deg=.5)
    view_change=shape_change(phase,profiles,narrow_view,data)
    active=bundle['dt']>0
    energy_error=np.max(np.abs(transport.energy_balance_error_erg))
    diagnostics={
        'work_energy_balance_max_absolute_erg':float(energy_error),
        'work_energy_balance_relative_to_max_segment_work':float(energy_error/np.max(transport.work_erg)),
        'gamma_mid_percentiles_1_50_99':np.percentile(transport.gamma_mid[active],[1,50,99]),
        'critical_GeV_percentiles_1_50_99':np.percentile(transport.critical_energy_gev[active],[1,50,99]),
        'particle_curvature_RLC_percentiles_1_50_99':np.percentile(bundle['rho'][active],[1,50,99]),
        'lab_direction_unit_norm_max_error':float(np.max(np.abs(np.linalg.norm(bundle['direction'][active],axis=-1)-1))),
        'refinement_rebinned_profile_L1_by_band':resolution,
        'resolution_passed_2_percent':bool(max(resolution)<.02),
        'view_bin_1_to_half_degree_L1_by_band':view_change,
        'max_emission_cylindrical_radius':float(np.max(np.linalg.norm(bundle['position'][active][:,:2],axis=-1))),
        'max_emission_spherical_radius':float(np.max(np.linalg.norm(bundle['position'][active],axis=-1))),
        'trajectories':int(bundle['dt'].shape[1]),'active_segments':int(active.sum()),
        'emitted_photons_per_injected_particle_by_band':np.einsum('snb,n->b',transport.photon_yield,bundle['seed_weight']),
    }
    # A spectral-integration refinement is independent of observed residuals.
    # Reuse the selected fine bundle to keep this check modest in memory.
    plot_particles(bundle,g,transport)
    selected_seeds=np.linspace(0,bundle['dt'].shape[1]-1,96,dtype=int)
    np.savez_compressed(OUT/'particle-trajectories.npz',
        position_rotating=bundle['position'][:,selected_seeds],dt=bundle['dt'][:,selected_seeds],
        direction_reference=bundle['direction'][:,selected_seeds],rho_RLC=bundle['rho'][:,selected_seeds],
        gamma=transport.gamma[:,selected_seeds],critical_energy_GeV=transport.critical_energy_gev[:,selected_seeds],
        seed_weight=bundle['seed_weight'][selected_seeds],pole=bundle['pole'][selected_seeds])
    del transport,bundle; gc.collect()
    check=trace_bundle(g,**plan['fine_quadrature'])
    _,_,spec3=calculate(check,g,rate,spectral_nodes=3)
    _,three=sky_profiles(check,spec3.photon_yield,g,nphase=1024)
    diagnostics['spectral_midpoint_to_3node_rebinned_L1_by_band']=shape_change(phase,middle,three,data)
    diagnostics['spectral_quadrature_passed_2_percent']=bool(max(diagnostics['spectral_midpoint_to_3node_rebinned_L1_by_band'])<.02)
    diagnostics['viewing_bin_sensitivity_below_2_percent']=bool(max(view_change)<.02)
    diagnostics['all_three_refinement_checks_passed']=bool(
        diagnostics['resolution_passed_2_percent'] and diagnostics['spectral_quadrature_passed_2_percent']
        and diagnostics['viewing_bin_sensitivity_below_2_percent'])
    del check,spec3; gc.collect()
    files=[ROOT/'run_particle_model.py',ROOT/'data/J0835-4510_3PC_data.fits',
           *sorted((ROOT/'src/pulsar').glob('*.py'))]
    report={
        'model':selected,'diagnostics':diagnostics,'bands':scores,'baseline_bands':baseline,
        'physical_model_validated_for_Vela':False,
        'claim':'Vela-tested reduced particle and curvature-emission model; not a Vela reconstruction',
        'source_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
        'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},
        'statistical_caveat':VALIDATION_CAVEAT,
    }
    write_json('validation.json',report)
    np.savez_compressed(OUT/'energy-profiles.npz',phase=phase,profiles=profiles,energies_GeV=ENERGIES)
    for d,s in zip(data,scores):
        np.savetxt(OUT/f"comparison-band-{d['index']}.csv",np.column_stack((d['phase'],d['counts'],d['error'],s['prediction'],s['residual_sigma'])),
                   delimiter=',',header='phase,weighted_counts,catalog_error,model,residual_sigma',comments='')
    plot_comparison(data,scores,baseline)
    rows='\n'.join(f"| {d['band']} | {d['role']} | {s['reduced_chi_squared']:.2f} | {old['reduced_chi_squared']:.2f} |" for d,s,old in zip(data,scores,baseline))
    verdict=('Large residuals remain: the current model does not reproduce the observed pulse shapes within the catalog errors.'
             if max(s['reduced_chi_squared'] for s in scores)>10 else
             'The shape scores alone do not establish physical validation of this emission hypothesis.')
    better=[d['band'] for d,s,old in zip(data,scores,baseline) if s['chi_squared']<old['chi_squared']]
    worse=[d['band'] for d,s,old in zip(data,scores,baseline) if s['chi_squared']>=old['chi_squared']]
    boundary=(' The selected inclination is at the edge of the searched grid; this is not an established optimum or a physical inclination measurement.'
              if selected['inclination_at_grid_boundary'] else '')
    text=f'''# Vela particle/emission-model test

The model now calculates curved laboratory trajectories, electric work, radiation reaction, curvature spectra, and photon arrival phases. It is an original reduced emission hypothesis tested against Vela, not a validated reconstruction of its plasma.

## Frozen fit and reserved-data result

The training-only grid selected magnetic inclination **{g.inclination_deg:g}°**, acceleration scale **{rate:g} cm⁻¹**, and common phase offset **{shift:.6f} cycles**. Observer angle remains 63.6°. These choices are model dependent and are not measurements of those physical parameters.{boundary}

| Energy | Role | New descriptive χ²/(N−2) | Earlier wind χ²/(N−2) |
|---|---|---:|---:|
{rows}

These scores use 100 bins per cycle. A value near unity would be the rough scale of agreement under the supplied error model; these are descriptive scores, not calibrated significance levels. Each band has its own nonnegative scale and constant, including reserved bands. No per-band phase, peak ratio, or energy cutoff is adjusted. {verdict} The two reserved bands were excluded from optimization, though earlier work had inspected Vela's aggregate morphology.

Compared with the earlier wind, the scores improve in {', '.join(better) or 'no bands'} and worsen in {', '.join(worse) or 'no bands'}. More detailed particle physics has not produced uniformly better agreement. In the highest-energy band, the modeled second peak fades while Vela's observed second peak is strong. The current magnetic geometry and acceleration/injection prescription cannot be justified as Vela's actual configuration by these results.

## Numerical evidence

- {diagnostics['trajectories']:,} injected trajectories; {diagnostics['active_segments']:,} active transport segments
- Relative segment energy-balance error: {diagnostics['work_energy_balance_relative_to_max_segment_work']:.3g}
- Fine-to-refined rebinned profile L1 differences: {', '.join(f'{100*x:.2f}%' for x in resolution)}
- Halving the numerical viewing-bin width changes profiles by: {', '.join(f'{100*x:.2f}%' for x in view_change)}
- One-node to three-node spectral integration changes profiles by: {', '.join(f'{100*x:.2f}%' for x in diagnostics['spectral_midpoint_to_3node_rebinned_L1_by_band'])}

L1 differences are relative to the mean profile, evaluated over the observational bins. Each refinement target is 2% in every band. Trajectory sampling: **{'passed' if diagnostics['resolution_passed_2_percent'] else 'not passed'}**. Spectral integration: **{'passed' if diagnostics['spectral_quadrature_passed_2_percent'] else 'not passed'}**. Viewing-bin sensitivity: **{'passed' if diagnostics['viewing_bin_sensitivity_below_2_percent'] else 'not passed'}**. This tolerance is numerical and does not imply convergence to the much smaller catalog statistical errors. The narrower viewing-bin experiment does not alter the frozen prediction. Energy balance alone verifies bookkeeping, while analytic radiation and trajectory tests independently check the calculation.

## What the result supports

Use **“Vela-tested particle and curvature-radiation model”**. Do not use **“Vela reconstruction”**. The retarded vacuum magnetic geometry with corotation is prescribed; it does not calculate plasma currents, pair creation, screening, synchrotron/IC radiation, or magnetic feedback. A fitted constant can absorb residual off-pulse emission and background. The uniform response approximation also prevents spectral or absolute-luminosity validation.

The next physical upgrade requires a validated numerical force-free/dissipative or kinetic field solution and a justified injection model, then detector-aware comparison with independent photon data. Adding arbitrary brightness at observed phases would not resolve the missing physics.

See [model equations and sources](../../PARTICLE-MODEL.md), [frozen protocol](protocol.json), [machine-readable validation](validation.json), and [data audit](../../research/VELA-DATA-AUDIT.md).
'''
    (OUT/'RESULTS.md').write_text(text)
    print(json.dumps({'selected':selected,'diagnostics':diagnostics,'scores':[s['reduced_chi_squared'] for s in scores]},default=native,indent=2),flush=True)


if __name__=='__main__': main()
