#!/usr/bin/env python3
"""Connect supplied field -> traced cap -> outer transport -> spectra/pulses.

This runner is deliberately a diagnostic, with no Vela-validation release.
It can inspect unreleased fields without laundering their numerical errors.
Physical and numerical parameter choices are recorded before integration.
"""
from pathlib import Path
import argparse,sys,os,json,hashlib,time
from dataclasses import asdict
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
from scipy.optimize import minimize_scalar
from pulsar.snapshot_fields import RotatingSnapshot,SpinDownScale
from pulsar.field_injection import CapConfig,snapshot_injection
from pulsar.published_acceleration import Barnard2022Work,NearForceFreeFIDO
from pulsar.emission_projection import project_emission
from pulsar.spectral_observations import load_vela_spectrum,compare_integrated_spectrum,ENERGY_EDGES_GEV
from pulsar.multiband import load_bands,score_profiles
from pulsar.outer_transport import OuterTraceConfig,trace_outer,signed_transport_curvature
from pulsar.radiation import MEC2_ERG,E_STATC

PULSE_BANDS=[(.1,.3),(.3,1),(1,3),(3,100)]

def serial(x):
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    if isinstance(x,np.ndarray):return serial(x.tolist())
    if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
    if isinstance(x,np.generic):return x.item()
    return x

def dump(path,x):path.write_text(json.dumps(serial(x),indent=2,allow_nan=False)+'\n')

class ScreenedDrift:
    """Explicit FFE-geometry approximation; leaves B and E_perp unchanged."""
    def __init__(self,base):self.base=base
    def __call__(self,p,t):
        e,b=self.base(p,t);b2=np.sum(b*b,axis=-1)
        if np.any(b2<=0):raise ValueError('Undefined screened drift at magnetic null')
        return e-b*(np.sum(e*b,axis=-1)/b2)[...,None],b

def diagnostic_scores(phase,profiles,data):
    # All inspected bands participate in development; none is relabeled blind.
    for row in data:row['role']='development'
    def score(shift):return sum(s['chi_squared'] for s in score_profiles(phase,profiles,data,shift%1))
    grid=np.arange(256)/256;values=np.array([score(x) for x in grid])
    index=int(np.argmin(values));shift=grid[index];value=values[index]
    minima=np.flatnonzero((values<=np.roll(values,1))&(values<=np.roll(values,-1)))
    for i in minima[np.argsort(values[minima])[:4]]:
        fit=minimize_scalar(score,bounds=(grid[i]-1/256,grid[i]+1/256),method='bounded',options={'xatol':1e-9})
        if fit.fun<value:shift=fit.x%1;value=fit.fun
    scores=score_profiles(phase,profiles,data,shift)
    for s in scores:
        s['role']='development'
        s['statistical_caveat']='All bands inspected; shared phase fitted to all bands, separate nonnegative scale/constant in each. Flat detector response. Source-specific exposure and source weights absent. Descriptive shape diagnostic, not joint likelihood or validation.'
    return {'phase_shift_cycles':float(shift),'bands':scores,'independent_validation':False}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/outer-diagnostic')
    parser.add_argument('--inclination',type=float,required=True,help='Must match the supplied field; not fitted by this script')
    parser.add_argument('--observer',type=float,default=65.)
    parser.add_argument('--seed-radius',type=float,required=True,help='Explicit resolved finite-radius injection surface in RLC')
    parser.add_argument('--seed-inner-radius',type=float,required=True)
    parser.add_argument('--dt',type=float,default=.005)
    parser.add_argument('--azimuths',type=int,default=64)
    parser.add_argument('--rings',type=int,default=3)
    parser.add_argument('--view-width',type=float,default=2.)
    parser.add_argument('--gamma0',type=float,default=1000.)
    parser.add_argument('--loading',choices=['magnetic_flux','surface_area'],default='magnetic_flux')
    parser.add_argument('--charge-sign',type=int,choices=[-1,1],default=-1,help='Charge of every injected FIDO primary; no automatic local charge reversal')
    parser.add_argument('--interpolation',choices=['cubic','linear'],default='cubic')
    parser.add_argument('--magnetic-phase',type=float,help='Explicit magnetic-axis phase at snapshot time (radians); overrides source metadata')
    parser.add_argument('--prescription',choices=['barnard-scenario2','fido'],default='barnard-scenario2')
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    measured=load_vela_spectrum(ROOT/'data/spectra/3PC_Catalog+SEDs_20230803.fits')
    scale=SpinDownScale(measured['period_s'],measured['period_derivative'],1e45,args.inclination)
    base=RotatingSnapshot.from_npz(args.snapshot,amplitude_gauss=scale.dipole_scale_gauss,interpolation=args.interpolation)
    # Native COFFEE exporter records inclination, legacy prototype its config.
    with np.load(args.snapshot,allow_pickle=False) as raw:
        config=json.loads(str(raw['config'])) if 'config' in raw else base.provenance
    inclination=config.get('inclination_deg',base.provenance.get('inclination_deg'))
    if inclination is not None and not np.isclose(inclination,args.inclination):
        raise ValueError('Requested inclination does not match snapshot metadata')
    time0=base.reference_time
    if args.magnetic_phase is not None:magnetic_phase=args.magnetic_phase
    elif 'magnetic_phase_rad' in base.provenance:magnetic_phase=float(base.provenance['magnetic_phase_rad'])
    elif 'spinup_time' in config:
        # The legacy solver integrates its spin-up angular velocity; its
        # magnetic phase is t-T/2 after spin-up, not t.
        T=float(config['spinup_time']);u=min(max(time0/T,0),1)
        magnetic_phase=T*(u**3-.5*u**4)+max(0,time0-T)
    else:magnetic_phase=time0  # omega=1 native COFFEE pulsar convention
    cap=CapConfig(inclination_deg=args.inclination,seed_radius=args.seed_radius,
                  inner_radius=args.seed_inner_radius,field_time=time0,magnetic_phase=magnetic_phase,
                  step=min(.008,args.dt*2),azimuths=args.azimuths,rings=args.rings)
    geometry=OuterTraceConfig(dt=args.dt,r_min=args.seed_inner_radius,r_max=2.,max_time=12.,curvature_step=min(.0001,args.dt/10))
    model=Barnard2022Work()
    fields=(ScreenedDrift(base) if args.prescription=='barnard-scenario2' else
            NearForceFreeFIDO(base,derivative_step=float(np.min(base.h))/2))
    protocol={'scope':'Unreleased field/emission development diagnostic; NOT Vela reconstruction',
        'snapshot':base.provenance,'field_release_allowed':base.verified,'cap':asdict(cap),
        'transport':asdict(geometry),'physical_scale':asdict(scale),'b_lc_gauss':scale.dipole_scale_gauss,
        'prescription':args.prescription,'barnard_rates':asdict(model) if args.prescription=='barnard-scenario2' else None,
        'fido_sigma_over_omega':30 if args.prescription=='fido' else None,
        'source_spectrum':measured['provenance'],'observer_deg':args.observer,'view_width_deg':args.view_width,
        'seed_loading':args.loading+' at finite injection radius; no physical pair supply',
        'initial_gamma':args.gamma0,'charge_choice':args.charge_sign,
        'charge_choice_scope':'FIDO signed electric work only; Barnard Eq.2 prescribes acceleration per laboratory path length independently of charge',
        'barnard_work_convention':'d gamma / d l = R_acc - curvature losses; no extra velocity-field projection',
        'photon_visibility':'No stellar occultation or photon absorption; inward-ray energy fraction is explicitly reported. Enlarged numerical star is not treated as the physical opaque radius.',
        'departure_from_paper':'finite-radius angular cap fractions, specified loading, field solver/interpolation and gamma0 differ; no exact reproduction claim',
        'spectrum_statistic':'single free normalization; upper limits separately tested',
        'profile_statistic':'shared phase; per-band nuisance normalizations; flat response',
        'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*sorted((ROOT/'src/pulsar').glob('*.py'))]},
        'emission_release_allowed':False,'vela_validation':False}
    dump(out/'protocol.json',protocol)
    start=time.monotonic()
    try:
        injection=snapshot_injection(base,cap,loading=args.loading)
        print('Traced numerical-field caps; launching '+str(len(injection['position']))+' primaries',flush=True)
        bundle=trace_outer(fields,injection['position'],branch=injection['branch'],charge_sign=args.charge_sign,
                           config=geometry,start_time=time0,seed_weight=injection['seed_weight'])
        if not bundle['all_escaped']:
            raise RuntimeError('Particle horizon reached before all streams escaped; incomplete emission cannot be scored')
        terminations={str(status):int(np.count_nonzero(bundle['status']==status))
                      for status in np.unique(bundle['status'])}
        # r_min is an artificial injection-region cutoff, not a physical
        # stellar surface. Returning trajectories end with unmodelled
        # inward emission; do not discard them and renormalize survivors.
        if np.any(bundle['status']=='r_min'):
            dump(out/'transport-termination.json',{'counts':terminations,
                'returning_seed_weight':float(bundle['seed_weight'][bundle['status']=='r_min'].sum()),
                'emission_release_allowed':False})
            raise RuntimeError('Streams return to the artificial inner cutoff; remaining inward transport is undefined, so no emission fit is produced')
        ds=bundle['dt']*scale.rlc_cm;rho=bundle['rho']*scale.rlc_cm
        spectral=list(zip(ENERGY_EDGES_GEV[:-1],ENERGY_EDGES_GEV[1:]))
        bands=PULSE_BANDS+spectral
        if args.prescription=='barnard-scenario2':
            work=np.zeros_like(bundle['dt']);active=bundle['dt']>0
            p=bundle['position'][active];v=bundle['direction'][active];t=bundle['time'][active]
            b=base(p,t)[1];work[active]=model.work_field(p,v,b)
            transported=signed_transport_curvature(ds,rho,work,args.gamma0,bands,spectral_nodes=3)
        else:
            work=bundle['e_work']
            transported=signed_transport_curvature(ds,rho,work,args.gamma0,bands,spectral_nodes=3)
        projected=project_emission(bundle,transported.photon_yield,observer_deg=args.observer,
                                   view_width_deg=args.view_width,nphase=512)
        data=load_bands(ROOT/'data/J0835-4510_3PC_data.fits')
        profiles=diagnostic_scores(projected['phase'],projected['profiles'][:4],data)
        spectrum=compare_integrated_spectrum(projected['phase_averaged_yield_per_sr'][4:],measured)
        radius=np.linalg.norm(bundle['position'],axis=-1);cylinder=np.linalg.norm(bundle['position'][...,:2],axis=-1)
        energy=transported.radiated_erg*bundle['seed_weight'][None,:]
        total=float(energy.sum())
        outer=float(energy[radius>=1].sum()/total) if total else 0.
        outside_cylinder=float(energy[cylinder>=1].sum()/total) if total else 0.
        radial_direction=np.sum(bundle['position']*bundle['direction'],axis=-1)
        inward_energy=float(energy[radial_direction<0].sum()/total) if total else 0.
        active=bundle['dt']>0
        e,b=fields(bundle['position'][active],bundle['time'][active])
        magnetic_invariant=np.sum(b*b-e*e,axis=-1)
        # Conservative energy-related magnetic length, not a measured gyroradius:
        # pitch angles are not evolved by this model.
        bstrength=np.sqrt(np.maximum(magnetic_invariant,0))
        magnetization_length=np.divide(transported.gamma_mid[active]*MEC2_ERG,E_STATC*bstrength,
            out=np.full(bstrength.shape,np.inf),where=bstrength>0)
        active_energy=energy[active]
        quality={'fraction_energy_E_ge_B':float(active_energy[magnetic_invariant<=0].sum()/total) if total else 0.,
            'fraction_energy_magnetic_length_above_tenth_cell':float(active_energy[magnetization_length>.1*min(base.h)*scale.rlc_cm].sum()/total) if total else 0.,
            'fraction_energy_curvature_radius_below_three_cells':float(active_energy[bundle['rho'][active]<3*min(base.h)].sum()/total) if total else 0.,
            'magnetic_length_definition':'gamma*m_e*c^2/(e*sqrt(max(B^2-E^2,0))); conservative scale diagnostic, not the exact drift-frame gyroradius, no pitch-angle or gyration solution'}
        report={'protocol':protocol,'elapsed_seconds':time.monotonic()-start,'trace_status':bundle['status'],
                'segment_shape':list(bundle['dt'].shape),'selected_seeds':projected['selected_seeds'],
                'fraction_energy_spherical_r_ge_1':outer,'fraction_energy_outside_light_cylinder':outside_cylinder,
                'termination_counts':terminations,'fraction_energy_inward_photons':inward_energy,
                'max_energy_balance_error_erg':float(np.max(np.abs(transported.energy_balance_error_erg))),
                'gamma_range':[float(transported.gamma.min()),float(transported.gamma.max())],
                'guiding_centre_diagnostics':quality,
                'profiles':profiles,'spectrum':spectrum,'emission_release_allowed':False,'vela_validation':False}
        dump(out/'results.json',report)
        np.savez_compressed(out/'emission.npz',**{k:bundle[k] for k in ('position','direction','rho','dt','time','phase','seed_weight')},
                            gamma=transported.gamma,critical_energy_gev=transported.critical_energy_gev,
                            radiated_erg=transported.radiated_erg,e_work=work,photon_yield=transported.photon_yield,
                            phase_grid=projected['phase'],profiles=projected['profiles'],cap_rim=injection['cap_rim'])
        print(json.dumps({'output':str(out),'elapsed_seconds':report['elapsed_seconds'],
                          'outside_light_cylinder_energy_fraction':outside_cylinder,
                          'selected_seeds':projected['selected_seeds'],'release':False}),flush=True)
    except Exception as exc:
        dump(out/'failure.json',{'error':type(exc).__name__,'message':str(exc),'elapsed_seconds':time.monotonic()-start,
                                'protocol':protocol,'emission_release_allowed':False})
        raise

if __name__=='__main__':main()
