#!/usr/bin/env python3
"""Reuse cached trajectories and original 3-node gamma-ray emission transport.

No field or particle trajectory is recomputed. Save actual emission positions
and directions, a deterministic display quadrature and the unmodified pulse.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from pulsar.snapshot_fields import SpinDownScale
from pulsar.outer_transport import signed_transport_curvature
from pulsar.radiation import MEC2_ERG,E_STATC
from pulsar.emission_quadrature import project_segment_emission
from pulsar.photon_scene import source_frame,observer_event,systematic_sample

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,default=ROOT/'physics/diagnostics/n97-final-illustration-a128-r6')
    p.add_argument('--output',type=Path,default=ROOT/'scene')
    p.add_argument('--geometry-cache',type=Path,help='Original six trace chunks; no particle integration is performed')
    p.add_argument('--packets',type=int,default=12288)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    target=a.output/'photons.npz'
    if target.exists():raise FileExistsError('Choose a fresh scene or reuse the existing one')
    proto=json.loads((a.run/'protocol.json').read_text());scale=SpinDownScale(**proto['physical_scale'])
    cache=a.geometry_cache or Path(proto['geometry_cache']);gain=proto['acceleration'];bands=[(.1,.3),(.3,1),(1,3),(3,100)]
    if len(list(cache.glob('trace-*.npz')))!=6:raise FileNotFoundError('Provide the six original trajectory-cache chunks with --geometry-cache; the included scene can be rendered directly')
    rlist=[];nlist=[];ylist=[];idlist=[];hashes={};profile=np.zeros((4,1000));seeds=0;segments=0
    maximum_phase_error=0.;maximum_norm_error=0.;start=time.monotonic()
    for path in sorted(cache.glob('trace-*.npz')):
        with np.load(path) as d:b={k:d[k] for k in d.files}
        if not np.all(b['status']=='r_max'):raise ValueError('Incomplete original trajectories')
        hashes[path.name]=sha(path);active=b['dt']>0
        rate=np.where(np.linalg.norm(b['position'],axis=-1)<gain['transition_radius_rlc'],gain['rate_inside_cm_inverse'],gain['rate_outside_cm_inverse'])
        work=np.where(active,rate*MEC2_ERG/E_STATC,0.)
        tr=signed_transport_curvature(b['dt']*scale.rlc_cm,b['rho']*scale.rlc_cm,work,proto['initial_gamma'],bands,spectral_nodes=3)
        profile+=project_segment_emission(b,tr.photon_yield,observer_deg=proto['observer_deg'],view_width_deg=2.,nphase=1000)['profiles']
        y=tr.photon_yield*b['seed_weight'][None,:,None]
        used=active&(y.sum(axis=-1)>0);ids=np.nonzero(used)
        r,n=source_frame(b['position'][used],b['direction'][used],b['time'][used])
        _,_,_,phase=observer_event(r,n)
        residual=(phase-b['phase'][used]+.5)%1-.5
        maximum_phase_error=max(maximum_phase_error,float(abs(residual).max()))
        maximum_norm_error=max(maximum_norm_error,float(abs(np.linalg.norm(n,axis=-1)-1).max()))
        rlist.append(r);nlist.append(n);ylist.append(y[used]);idlist.append(np.stack((ids[0],ids[1]+seeds),axis=-1))
        seeds+=len(b['seed_weight']);segments+=int(active.sum())
        print(json.dumps(dict(seeds=seeds,active_segments=segments,elapsed_seconds=time.monotonic()-start)),flush=True)
    r=np.concatenate(rlist);n=np.concatenate(nlist);y=np.concatenate(ylist);source_ids=np.concatenate(idlist)
    with np.load(a.run/'profiles.npz') as d:
        original=d['profiles_width_2'][:4];edges=d['phase_edges']
    relative=float(np.max(abs(profile-original))/np.max(original))
    if relative>1e-10 or maximum_phase_error>1e-10 or maximum_norm_error>1e-10:
        raise ValueError('Reused calculation failed original profile/direction/arrival audit')
    # Sample the calculated photon number, stratified in polar direction and
    # arrival phase. Repetition of a cell changes its display weight, never its
    # direction or spectrum. Four channels retain the original energy bands.
    phase=observer_event(r,n)[-1]
    group=np.clip(((n[:,2]+1)*32).astype(int),0,63)*240+(phase*240).astype(int)
    order=np.lexsort((phase,group));weights=y.sum(axis=-1)
    selected,packet_weight=systematic_sample(weights[order],a.packets)
    selected=order[selected]
    fractions=y[selected]/weights[selected,None]
    # Audit the full sky at declared, coarser diagnostic bins, not a fitted curve.
    key=np.clip(((n[:,2]+1)*8).astype(int),0,15)*48+(phase*48).astype(int)
    full=np.bincount(key,weights=weights,minlength=16*48)
    sampled=np.bincount(key[selected],weights=np.full(a.packets,packet_weight),minlength=16*48)
    sky_l1=float(abs(full-sampled).sum()/full.sum())
    meta=dict(scope='Vela-informed calculated 0.1–100 GeV outgoing photon visualization; not a quantitative Vela reconstruction',
        source_run=str(a.run.resolve()),source_profiles_sha256=sha(a.run/'profiles.npz'),
        protocol_sha256=sha(a.run/'protocol.json'),geometry_sha256=hashes,
        reused_physics='Same force-free snapshot, cached trajectories, prescribed acceleration and signed curvature transport with original 3-node spectral integration',
        new_field_solves=0,new_particle_trajectories=0,active_segments=segments,seed_count=seeds,
        bands_gev=bands,period_s=scale.period_s,rlc_km=scale.rlc_cm/1e5,
        observer_deg=proto['observer_deg'],observer_width_deg=2.,
        source_extent_rlc=[.7,2.],packet_count=a.packets,packet_weight=float(packet_weight),
        sampling='Deterministic equal-photon-number systematic quadrature, ordered by 64 polar-direction and 240 arrival-phase strata; each packet retains an actual source-cell midpoint and exact lab direction; no angular convolution',
        total_photons_per_normalized_primary=float(weights.sum()),band_photons=y.sum(axis=0).tolist(),
        sky_diagnostic_bins=[16,48],sampled_sky_relative_L1=sky_l1,
        packet_spectra='Original four-band photon fractions attached to each selected cell; additive photon-number weights',
        phase_convention='Phi=(t_em-azimuth(n)-n.r)/(2pi) modulo 1; distant angular observer, original 2deg polar pixel',
        packet_geometry='Source midpoint quadrature; no interpolation into artificial cones; directions fixed after emission',
        quantitative_vela_validation=False)
    audit=dict(passed=True,original_profile_max_relative_difference=relative,
        all_source_midpoint_phase_max_error_cycles=maximum_phase_error,
        all_source_direction_max_norm_error=maximum_norm_error,
        sampled_sky_relative_L1=sky_l1,
        note='Receiver uses full conservative cell/latitude/phase quadrature; finite display packets are not used to estimate the pulse')
    np.savez_compressed(target,position0=r[selected],direction0=n[selected],band_fraction=fractions,
        packet_weight=np.full(a.packets,packet_weight),source_cell=source_ids[selected],
        phase_edges=edges,pulse_bands=original,metadata=json.dumps(meta))
    (a.output/'scene.json').write_text(json.dumps(meta,indent=2)+'\n')
    (a.output/'calculation-verification.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2));print(target)

if __name__=='__main__':main()
