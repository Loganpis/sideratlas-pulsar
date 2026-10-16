#!/usr/bin/env python3
"""Measured emission quadrature on one declared, unreleased field snapshot."""
from pathlib import Path
import argparse,hashlib,json,sys,time
from dataclasses import asdict
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.optimize import minimize_scalar
from pulsar.snapshot_fields import RotatingSnapshot,SpinDownScale
from pulsar.field_injection import CapConfig,snapshot_injection
from pulsar.outer_transport import OuterTraceConfig,trace_outer,signed_transport_curvature
from pulsar.published_acceleration import Barnard2022Work
from pulsar.emission_quadrature import project_segment_emission,reconstruct_sky_endpoints,rebin_phase_averages
from pulsar.spectral_observations import load_vela_spectrum,compare_integrated_spectrum,ENERGY_EDGES_GEV
from pulsar.multiband import load_bands
from pulsar.observations import nonnegative_scale_and_offset


def serial(x):
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    if isinstance(x,np.ndarray):return serial(x.tolist())
    if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
    if isinstance(x,np.generic):return x.item()
    return x


def dump(path,x):path.write_text(json.dumps(serial(x),indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class ScreenedGeometry:
    def __init__(self,base):self.base=base
    def __call__(self,p,t):
        e,b=self.base(p,t);bb=np.sum(b*b,axis=-1)
        if np.any(bb<=0):raise ValueError('Magnetic null in screened geometry')
        return e-b*(np.sum(e*b,axis=-1)/bb)[...,None],b


def score_binned(edges,profiles,data):
    """Fit one shared phase and per-band nuisance A,C with exact overlaps."""
    def score(shift,details=False):
        rows=[]
        for band,profile in zip(data,profiles):
            m=rebin_phase_averages(edges,profile,band['phase_lo'],band['phase_hi'],shift)
            maximum=m.max();m=m/maximum if maximum>0 else m
            coef=nonnegative_scale_and_offset(m,band['counts'],band['error'])
            chi=float(np.sum(((coef[0]*m+coef[1]-band['counts'])/band['error'])**2))
            rows.append({'band':band['band'],'chi_squared':chi,'descriptive_dof':len(m)-2})
        return rows if details else sum(x['chi_squared'] for x in rows)
    grid=np.arange(256)/256;values=np.array([score(x) for x in grid]);best=int(np.argmin(values));phase=grid[best];value=values[best]
    minima=np.flatnonzero((values<=np.roll(values,1))&(values<=np.roll(values,-1)))
    for i in minima[np.argsort(values[minima])[:4]]:
        fit=minimize_scalar(score,bounds=(grid[i]-1/256,grid[i]+1/256),method='bounded',options={'xatol':1e-9})
        if fit.fun<value:phase=fit.x%1;value=fit.fun
    return {'phase_shift_cycles':float(phase),'bands':score(phase,True),
            'caveat':'Development only; common fitted phase, per-band nonnegative scale and constant; exact finite-bin overlaps, flat detector response; no independent validation or calibrated likelihood.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('snapshot',type=Path);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--azimuths',type=int,default=64);ap.add_argument('--rings',type=int,default=3)
    ap.add_argument('--dt',type=float,default=.004);ap.add_argument('--batch-size',type=int,default=256)
    ap.add_argument('--observer',type=float,default=65.);ap.add_argument('--widths',type=float,nargs='+',default=[2.,1.])
    ap.add_argument('--nphase',type=int,default=1000)
    ap.add_argument('--seed-radius',type=float,default=.7);ap.add_argument('--inner-radius',type=float,default=.68)
    ap.add_argument('--cap-boundary',choices=('light_cylinder','outgoing_sphere'),default='light_cylinder',
                    help='Finite-domain rim criterion; outgoing sphere requires every selected seed to reach that same outgoing boundary')
    ap.add_argument('--cache-root',type=Path,default=ROOT/'outputs/emission-quadrature/geometry-cache')
    args=ap.parse_args();out=args.output
    if (out/'results.json').exists():raise FileExistsError('Completed output exists; choose another directory')
    out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    measured=load_vela_spectrum(ROOT/'data/spectra/3PC_Catalog+SEDs_20230803.fits')
    with np.load(args.snapshot,allow_pickle=False) as raw:metadata=json.loads(str(raw['metadata']))
    alpha=float(metadata['inclination_deg']);t0=float(metadata['magnetic_phase_rad'])
    scale=SpinDownScale(measured['period_s'],measured['period_derivative'],1e45,alpha)
    base=RotatingSnapshot.from_npz(args.snapshot,amplitude_gauss=scale.dipole_scale_gauss)
    cap=CapConfig(inclination_deg=alpha,seed_radius=args.seed_radius,inner_radius=args.inner_radius,
                  field_time=base.reference_time,magnetic_phase=t0,step=.008,azimuths=args.azimuths,rings=args.rings,
                  rim_boundary=args.cap_boundary)
    trace=OuterTraceConfig(dt=args.dt,r_min=args.inner_radius,r_max=2.,max_time=12.,curvature_step=.0001)
    sources={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),*sorted((ROOT/'src/pulsar').glob('*.py'))]}
    keydata={'field_sha256':base.provenance['sha256'],'cap':asdict(cap),'transport':asdict(trace),
             'geometry':'screened E_perp/B rotating snapshot','scale':asdict(scale),'batch_size':args.batch_size,
             'geometry_source_hashes':{k:v for k,v in sources.items() if k in ['run_emission_quadrature.py',
                 'src/pulsar/outer_transport.py','src/pulsar/snapshot_fields.py','src/pulsar/field_injection.py','src/pulsar/trajectories.py']}}
    key=hashlib.sha256(json.dumps(keydata,sort_keys=True).encode()).hexdigest();cache=args.cache_root/key
    cache.mkdir(parents=True,exist_ok=True)
    if (cache/'manifest.json').exists():
        if json.loads((cache/'manifest.json').read_text())!=keydata:raise ValueError('Geometry cache key mismatch')
    else:dump(cache/'manifest.json',keydata)
    if (cache/'injection.npz').exists():
        with np.load(cache/'injection.npz',allow_pickle=False) as cached:injection={k:cached[k] for k in cached.files}
    else:
        inj=snapshot_injection(base,cap,loading='magnetic_flux')
        injection={k:inj[k] for k in ('position','branch','seed_weight','cap_rim','cap_azimuth')}
        np.savez_compressed(cache/'injection.npz',**injection)
    protocol={'scope':'Emission quadrature study on one unreleased field; no Vela reconstruction',
        'snapshot':base.provenance,'cap':asdict(cap),'transport':asdict(trace),'physical_scale':asdict(scale),
        'observer_deg':args.observer,'observer_widths_deg':args.widths,'phase_bins':args.nphase,
        'initial_gamma':1000.,'acceleration':asdict(Barnard2022Work()),
        'acceleration_convention':'Barnard Eq.2 per laboratory path; spherical radius switch; no velocity projection',
        'seed_loading':'magnetic flux at finite injection surface; normalized, no inferred primary supply',
        'cap_boundary_scope':('Strict light-cylinder crossing rim' if args.cap_boundary=='light_cylinder' else
            'Finite-domain boundary-open rim and seeds: cylinder crossing or verified outgoing escape at r=2.5 RLC; not a proof of opening at infinity; every selected annular seed separately reaches this outgoing boundary'),
        'projection':'conservative cell-time fractions through finite latitude/phase bins; linear midpoint reconstruction',
        'source_hashes':sources,'geometry_cache':str(cache.resolve()),'geometry_cache_key':key,
        'field_release_allowed':False,'emission_release_allowed':False,'vela_validation':False}
    dump(out/'protocol.json',protocol)
    bands=[(.1,.3),(.3,1.),(1.,3.),(3.,100.)]+list(zip(ENERGY_EDGES_GEV[:-1],ENERGY_EDGES_GEV[1:]))
    result={str(width):{'profiles':np.zeros((12,args.nphase)),'selected_seeds':0,'selected_segments':0,
                        'selected_cell_fraction_sum':0.} for width in args.widths}
    total_energy=0.;outside_energy=0.;inward_energy=0.;max_balance=0.;clipped=0;phase_jump=0.;statuses={}
    gamma_min=np.inf;gamma_max=0.;samples=0;fields=ScreenedGeometry(base);model=Barnard2022Work()
    try:
        count=len(injection['position'])
        for begin in range(0,count,args.batch_size):
            end=min(begin+args.batch_size,count);path=cache/f'trace-{begin:06d}-{end:06d}.npz'
            if path.exists():
                with np.load(path,allow_pickle=False) as saved:bundle={k:saved[k] for k in saved.files}
            else:
                b=trace_outer(fields,injection['position'][begin:end],branch=injection['branch'][begin:end],
                    charge_sign=-1,config=trace,start_time=base.reference_time,seed_weight=injection['seed_weight'][begin:end])
                bundle={k:b[k] for k in ('position','direction','rho','dt','time','phase','seed_weight','branch','status','end_position','end_time')}
                np.savez_compressed(path,**bundle)
            if not np.all(bundle['status']=='r_max'):raise RuntimeError('Not every stream reached outer boundary; no emission score')
            for label in bundle['status']:statuses[str(label)]=statuses.get(str(label),0)+1
            active=bundle['dt']>0;work=np.zeros_like(bundle['dt']);p=bundle['position'][active]
            work[active]=model.work_field(p,bundle['direction'][active],base(p,bundle['time'][active])[1])
            tr=signed_transport_curvature(bundle['dt']*scale.rlc_cm,bundle['rho']*scale.rlc_cm,work,1000.,bands,spectral_nodes=3)
            mu,phase,diag=reconstruct_sky_endpoints(bundle);clipped+=diag['mu_endpoints_limited'];phase_jump=max(phase_jump,diag['max_midpoint_phase_jump_cycles'])
            for width in args.widths:
                pr=project_segment_emission(bundle,tr.photon_yield,observer_deg=args.observer,view_width_deg=width,nphase=args.nphase,sky_endpoints=(mu,phase))
                for k in result[str(width)]:result[str(width)][k]+=pr[k]
            energy=tr.radiated_erg*bundle['seed_weight'][None];total_energy+=float(energy.sum())
            outside_energy+=float(energy[np.linalg.norm(bundle['position'][...,:2],axis=-1)>=1].sum())
            inward_energy+=float(energy[np.sum(bundle['position']*bundle['direction'],axis=-1)<0].sum())
            max_balance=max(max_balance,float(np.max(abs(tr.energy_balance_error_erg))))
            gamma_min=min(gamma_min,float(tr.gamma.min()));gamma_max=max(gamma_max,float(tr.gamma.max()));samples+=int(active.sum())
            print(json.dumps({'seeds_completed':end,'total':count,'elapsed_s':time.monotonic()-start}),flush=True)
        data=load_bands(ROOT/'data/J0835-4510_3PC_data.fits');edges=np.arange(args.nphase+1)/args.nphase
        arrays={'phase_edges':edges,'phase':(edges[:-1]+edges[1:])/2,'cap_rim':injection['cap_rim'],'seed_weight':injection['seed_weight']}
        summaries={}
        for width in args.widths:
            row=result[str(width)];profile=row.pop('profiles');means=profile.mean(axis=1)
            arrays[f'profiles_width_{width:g}']=profile
            summaries[str(width)]={**row,'phase_averaged_yield_per_sr':means,
                'spectrum':compare_integrated_spectrum(means[4:],measured),'pulse_shape':score_binned(edges,profile[:4],data)}
        report={'protocol':protocol,'elapsed_seconds':time.monotonic()-start,'seed_count':count,'active_segments':samples,
                'termination_counts':statuses,'gamma_range':[gamma_min,gamma_max],'energy_erg_per_primary':total_energy,
                'fraction_energy_outside_light_cylinder':outside_energy/total_energy,'fraction_energy_inward_photons':inward_energy/total_energy,
                'max_energy_balance_error_erg':max_balance,'mu_endpoints_limited':clipped,'max_midpoint_phase_jump_cycles':phase_jump,
                'width_scenarios':summaries,'emission_release_allowed':False,'vela_validation':False}
        np.savez_compressed(out/'profiles.npz',**arrays);dump(out/'results.json',report)
        print(json.dumps({'output':str(out),'elapsed_s':report['elapsed_seconds']}),flush=True)
    except Exception as exc:
        dump(out/'failure.json',{'error':type(exc).__name__,'message':str(exc),'protocol':protocol});raise


if __name__=='__main__':main()
