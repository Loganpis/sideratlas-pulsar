#!/usr/bin/env python3
"""Test global acceleration parameters on fixed, recorded particle geometry.

No phase-dependent weights or peak-specific gains are fitted. The directional
yield is integrated over the declared latitude bin and all arrival phases.
This development fit is conditional on its unreleased force-free geometry;
it cannot certify that geometry or reconstruct a particle supply.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time
from dataclasses import asdict
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from pulsar.emission_quadrature import reconstruct_sky_endpoints,project_segment_emission
from pulsar.outer_transport import signed_transport_curvature
from pulsar.snapshot_fields import SpinDownScale
from pulsar.radiation import MEC2_ERG,E_STATC
from pulsar.spectral_observations import load_vela_spectrum,compare_integrated_spectrum,ENERGY_EDGES_GEV
from pulsar.energetics import population_power_budget


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def serial(x):
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    if isinstance(x,np.ndarray):return serial(x.tolist())
    if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
    if isinstance(x,np.generic):return x.item()
    return x


def observer_cell_fraction(mu,observer,width):
    """Exact residence fraction for the same linear mu(t) reconstruction."""
    if (not np.isfinite([observer,width]).all() or not width>0
            or observer-width/2<0 or observer+width/2>180):
        raise ValueError('Finite sky bin must lie within 0 to 180 degrees')
    low,high=np.deg2rad([max(0,observer-width/2),min(180,observer+width/2)])
    minimum,maximum=np.cos(high),np.cos(low)
    a,b=mu[...,0],mu[...,1];delta=b-a
    start=np.minimum(a,b);end=np.maximum(a,b)
    fraction=np.divide(np.maximum(0,np.minimum(end,maximum)-np.maximum(start,minimum)),
        abs(delta),out=np.zeros_like(delta),where=delta!=0)
    fraction=np.where(delta==0,((a>=minimum)&(a<=maximum)).astype(float),fraction)
    return fraction,float(2*np.pi*(maximum-minimum))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('geometry_run',type=Path,help='Completed run_emission_quadrature output')
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--rate-pairs',type=float,nargs='+',default=[.04,.25,.08,.5,.16,1.,.32,2.])
    ap.add_argument('--observer',type=float,default=65.)
    ap.add_argument('--width',type=float,default=2.)
    ap.add_argument('--distance-pc',type=float,default=290.)
    args=ap.parse_args()
    if len(args.rate_pairs)%2 or any(not np.isfinite(x) or x<0 for x in args.rate_pairs):
        raise ValueError('Supply nonnegative finite inside/outside acceleration pairs')
    output=args.output
    if (output/'results.json').exists():raise FileExistsError('Choose a fresh scan output')
    output.mkdir(parents=True,exist_ok=True)
    source=args.geometry_run/'results.json';run=json.loads(source.read_text());protocol=run['protocol']
    cache=Path(protocol['geometry_cache']);manifest=cache/'manifest.json'
    if sha_bytes:=protocol.get('geometry_cache_key'):
        actual=hashlib.sha256(json.dumps(json.loads(manifest.read_text()),sort_keys=True).encode()).hexdigest()
        if actual!=sha_bytes:raise ValueError('Geometry manifest no longer matches recorded key')
    manifest_data=json.loads(manifest.read_text())
    if manifest_data['field_sha256']!=protocol['snapshot']['sha256']:
        raise ValueError('Field hash mismatch')
    with np.load(cache/'injection.npz') as saved:weights=saved['seed_weight']
    if not np.isclose(weights.sum(),1,rtol=1e-12,atol=1e-12):raise ValueError('Globally normalized seed population required')
    scale=SpinDownScale(**protocol['physical_scale'])
    measured=load_vela_spectrum(ROOT/'data/spectra/3PC_Catalog+SEDs_20230803.fits')
    bands=list(zip(ENERGY_EDGES_GEV[:-1],ENERGY_EDGES_GEV[1:]))
    chunks=[];total_weights=0.
    for path in sorted(cache.glob('trace-*.npz')):
        with np.load(path,allow_pickle=False) as raw:bundle={k:raw[k] for k in raw.files}
        if not np.all(bundle['status']=='r_max'):raise ValueError('Incomplete particle geometry cannot be fitted')
        mu,phase,diag=reconstruct_sky_endpoints(bundle)
        fraction,solid_angle=observer_cell_fraction(mu,args.observer,args.width)
        chunks.append((bundle,fraction,solid_angle,path))
        total_weights+=float(bundle['seed_weight'].sum())
    if not np.isclose(total_weights,1,rtol=1e-12,atol=1e-12):raise ValueError('Incomplete/duplicated geometry chunks')
    source_hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),
        ROOT/'src/pulsar/outer_transport.py',ROOT/'src/pulsar/radiation.py',
        ROOT/'src/pulsar/emission_quadrature.py',ROOT/'src/pulsar/spectral_observations.py',ROOT/'src/pulsar/energetics.py']}
    base=dict(scope='Conditional global acceleration response on unreleased fixed geometry; development only',
        source_run=str(source.resolve()),source_sha256=sha(source),geometry_protocol=protocol,
        geometry_files_sha256={str(p):sha(p) for *_,p in chunks},source_hashes=source_hashes,
        observer_deg=args.observer,view_width_deg=args.width,
        free_parameters='Two global nonnegative gains dgamma/dl in cm^-1 (spherical switch r=RLC); one spectrum-wide population normalization',
        distance_pc=args.distance_pc,distance_scope='Explicit Barnard-like adoption, not a newly inferred distance',
        projection='Exact finite latitude-bin residence fraction of linear sky endpoints; integrated over all phase',
        changes_geometry=False,emission_release_allowed=False,vela_validation=False)
    (output/'protocol.json').write_text(json.dumps(serial(base),indent=2,allow_nan=False)+'\n')
    rows=[];started=time.monotonic()
    for inside,outside in np.array(args.rate_pairs).reshape(-1,2):
        yields=np.zeros(8);energy=0.;work_total=0.;initial=0.;escape=0.;ledger=0.
        for bundle,fraction,solid_angle,path in chunks:
            active=bundle['dt']>0
            gain=np.where(np.linalg.norm(bundle['position'],axis=-1)<1,inside,outside)
            work=np.where(active,gain*MEC2_ERG/E_STATC,0.)
            tr=signed_transport_curvature(bundle['dt']*scale.rlc_cm,bundle['rho']*scale.rlc_cm,
                work,protocol['initial_gamma'],bands,spectral_nodes=3)
            weight=bundle['seed_weight'];angular_weight=fraction*weight[None]/solid_angle
            yields+=np.einsum('snb,sn->b',tr.photon_yield,angular_weight)
            energy+=float(np.sum(tr.radiated_erg*weight[None]))
            work_total+=float(np.sum(tr.work_erg*weight[None]))
            initial+=float(np.sum((tr.gamma[0]-1)*MEC2_ERG*weight))
            escape+=float(np.sum((tr.gamma[-1]-1)*MEC2_ERG*weight))
            ledger=max(ledger,float(np.max(abs(tr.energy_balance_error_erg))))
        fit=compare_integrated_spectrum(yields,measured)
        budget=population_power_budget(fit['normalization'],args.distance_pc,scale.luminosity_erg_s,energy,
            work_erg_per_primary=work_total,injected_kinetic_erg_per_primary=initial,
            escaping_kinetic_erg_per_primary=escape)
        row=dict(rate_inside_cm_inverse=float(inside),rate_outside_cm_inverse=float(outside),
            phase_averaged_yield_per_sr=yields,spectrum=fit,power_budget=budget,
            max_energy_ledger_error_erg=ledger,elapsed_seconds=time.monotonic()-started)
        rows.append(row)
        print(json.dumps({'rates':[float(inside),float(outside)],'spectrum_chi2':fit['descriptive_chi_squared'],
            'supply_fraction':budget['modeled_supply_fraction_of_spin_down'],'elapsed_s':time.monotonic()-started}),flush=True)
        (output/'results.json').write_text(json.dumps(serial({**base,'trials':rows}),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
