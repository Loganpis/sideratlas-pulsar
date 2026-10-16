#!/usr/bin/env python3
"""Prepare a recorded all-direction emissivity scene, not a camera photograph.

The movie shows a rigidly rotating snapshot of outer-magnetosphere curvature
emission. Particle paths, gain/loss transport and the pulse are calculated;
the scene colors, exposure, voxel/pixel footprint and time scale are display
choices. This does not upgrade the source model's observational validation.
"""
from pathlib import Path
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.ndimage import rotate
from pulsar.snapshot_fields import SpinDownScale,rotate_z
from pulsar.outer_transport import signed_transport_curvature
from pulsar.radiation import MEC2_ERG,E_STATC,curvature_band_energy_fraction


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def deposit_voxels(grid,points,values,extent):
    """Conservative trilinear deposit; value has one additive channel axis."""
    n=grid.shape[0];u=(points+extent)*(n-1)/(2*extent)
    base=np.floor(u).astype(int);fraction=u-base
    if np.any(base<0) or np.any(base+1>=n):raise ValueError('Emissivity outside scene volume')
    flat=grid.reshape(-1,grid.shape[-1])
    for dx in (0,1):
        for dy in (0,1):
            for dz in (0,1):
                off=np.array([dx,dy,dz]);ij=base+off
                w=np.prod(np.where(off,fraction,1-fraction),axis=-1)
                ids=(ij[:,0]*n+ij[:,1])*n+ij[:,2]
                for channel in range(values.shape[-1]):
                    flat[:,channel]+=np.bincount(ids,weights=w*values[:,channel],minlength=n**3)


def continuous_display_volume(volume,azimuths,samples=9):
    """Positive finite-angle interpolation for displaying sampled streams.

    A triangular angular footprint of one source-azimuth cell on each side
    suppresses visible sampling spokes. It is applied to the display volume,
    not to trajectories, radiation spectra or arrival-phase calculations.
    Rotation interpolation is linear; one global conservation correction
    compensates voxel resampling roundoff, never a frame/peak adjustment.
    """
    if azimuths<4 or samples<3 or samples%2!=1:raise ValueError('Valid angular display quadrature required')
    original=float(volume[...,0].sum());output=np.zeros_like(volume)
    offsets=np.linspace(-1,1,samples);weight=1-abs(offsets);weight/=weight.sum()
    for offset,w in zip(offsets,weight):
        if w:
            output+=w*rotate(volume,float(offset*360/azimuths),axes=(0,1),reshape=False,
                             order=1,mode='constant',cval=0.,prefilter=False)
    retained=float(output[...,0].sum())
    if retained<=0:raise ValueError('Display angular interpolation erased the scene')
    correction=original/retained;output*=correction
    return output,dict(type='Positive triangular angular interpolation of display volume about spin axis',
        half_width_deg=360/azimuths,quadrature_samples=samples,global_energy_conservation_correction=correction,
        scope='Display-only finite angular footprint; no change to model trajectories, spectra or pulse profiles')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--grid',type=int,default=144)
    p.add_argument('--width',type=float,default=2.)
    p.add_argument('--angular-display-samples',type=int,default=9)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if (a.output/'scene.npz').exists():raise FileExistsError('Choose a fresh scene directory')
    report=json.loads((a.run/'results.json').read_text());proto=report['protocol']
    scale=SpinDownScale(**proto['physical_scale']);cache=Path(proto['geometry_cache'])
    n=a.grid;extent=2.08
    if n<24:raise ValueError('At least24display voxels required')
    volume=np.zeros((n,n,n,2),float);traces=[];hashes={};supplied_weight=0.;total=0.;max_error=0.
    gain=proto['acceleration'];seeds=report['seed_count'];shown_tracks=0
    for path in sorted(cache.glob('trace-*.npz')):
        hashes[str(path.resolve())]=sha(path)
        with np.load(path,allow_pickle=False) as saved:b={k:saved[k] for k in saved.files}
        if not np.all(b['status']=='r_max'):raise ValueError('Incomplete streams cannot form a final scene')
        active=b['dt']>0;rate=np.where(np.linalg.norm(b['position'],axis=-1)<gain['transition_radius_rlc'],
            gain['rate_inside_cm_inverse'],gain['rate_outside_cm_inverse'])
        work=np.where(active,rate*MEC2_ERG/E_STATC,0.)
        tr=signed_transport_curvature(b['dt']*scale.rlc_cm,b['rho']*scale.rlc_cm,work,
                                      proto['initial_gamma'],spectral_nodes=3)
        # This display channel uses each cell's midpoint characteristic energy;
        # scientific pulse/spectral fits retain their original 3-node quadrature.
        energy=tr.radiated_erg*curvature_band_energy_fraction(tr.critical_energy_gev,[(.1,100.)])[...,0]
        energy*=b['seed_weight'][None]
        valid=active&(energy>0)&np.isfinite(tr.critical_energy_gev)&(tr.critical_energy_gev>0)
        points=rotate_z(b['position'][valid],-b['time'][valid])
        values=np.stack((energy[valid],energy[valid]*np.log10(tr.critical_energy_gev[valid])),axis=-1)
        deposit_voxels(volume,points,values,extent)
        total+=float(energy.sum());supplied_weight+=float(b['seed_weight'].sum())
        max_error=max(max_error,float(np.max(abs(tr.energy_balance_error_erg))))
        # Paths are faint geometric guides, with no claim that they are B lines.
        for seed in range(b['dt'].shape[1]):
            global_seed=shown_tracks+seed
            if global_seed%max(1,seeds//64)==0:
                mask=active[:,seed];path_points=rotate_z(b['position'][mask,seed],-b['time'][mask,seed])
                if len(path_points)>1:
                    ids=np.unique(np.linspace(0,len(path_points)-1,min(150,len(path_points))).astype(int))
                    traces.append(path_points[ids])
        shown_tracks+=b['dt'].shape[1]
        print(json.dumps({'paths_processed':shown_tracks,'paths_total':seeds}),flush=True)
    if not np.isclose(supplied_weight,1,rtol=1e-12,atol=1e-12):raise ValueError('Incomplete or duplicated seed population')
    if not np.isclose(volume[...,0].sum(),total,rtol=2e-12,atol=1e-14):raise ValueError('Display deposition lost emitted energy')
    volume,angular_display=continuous_display_volume(volume,proto['cap']['azimuths'],a.angular_display_samples)
    nonzero=volume[...,0]>0;ids=np.argwhere(nonzero)
    points=ids*(2*extent/(n-1))-extent;power=volume[...,0][nonzero]
    logenergy=volume[...,1][nonzero]/power
    trackpoints=np.concatenate(traces) if traces else np.empty((0,3))
    offsets=np.r_[0,np.cumsum([len(v) for v in traces])]
    with np.load(a.run/'profiles.npz') as saved:
        edges=saved['phase_edges'];profile=saved[f'profiles_width_{a.width:g}'][:4].sum(axis=0)
    meta=dict(scope='User-approved Vela-informed visualization; explicitly approximate, not a quantitative reconstruction',
        source_run=str(a.run.resolve()),source_result_sha256=sha(a.run/'results.json'),
        source_profiles_sha256=sha(a.run/'profiles.npz'),geometry_sha256=hashes,
        physical_protocol=proto,quantity='0.1–100GeV curvature energy emitted in all directions per steady primary supply; optically thin volume display',
        display_band_gev=[.1,100.],display_band_quadrature='Cell midpoint critical energy applied to integrated cell energy; pulse retains original 3-node spectral quadrature',
        color_quantity='log10(cell-energy-weighted characteristic curvature photon energy in GeV), illustrative palette',
        projection='Orthographic line-of-sight integral of angle-integrated emitted energy; not a photograph or beamed camera intensity',
        pattern='Rigid rotation of one numerical snapshot, not time-dependent PIC plasma',
        stellar_marker='Enlarged position marker; no surface brightness or physical stellar radius reconstructed',
        modeled_radius_rlc=[proto['cap']['seed_radius'],proto['transport']['r_max']],
        voxel_grid=n,voxel_extent_rlc=extent,voxel_deposition='Conservative trilinear display grid',
        angular_display_interpolation=angular_display,
        display_energy_erg_per_primary=total,display_deposition_relative_error=float(abs(volume[...,0].sum()-total)/total),
        max_particle_energy_ledger_error_erg=max_error,
        seed_count=seeds,observer_deg=proto['observer_deg'],observer_width_deg=a.width,
        period_s=scale.period_s,rlc_km=scale.rlc_cm/1e5,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'src/pulsar/radiation.py',ROOT/'src/pulsar/outer_transport.py']},
        quantitative_vela_validation=False)
    np.savez_compressed(a.output/'scene.npz',points=points.astype('f4'),energy=power,
        log10_critical_energy_gev=logenergy.astype('f4'),track_points=trackpoints.astype('f4'),track_offsets=offsets,
        phase_edges=edges,pulse_photons_per_sr=profile,metadata=json.dumps(meta))
    (a.output/'scene.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(a.output/'scene.npz',flush=True)


if __name__=='__main__':main()
