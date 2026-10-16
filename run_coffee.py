#!/usr/bin/env python3
"""Reproducible local COFFEE runs and explicit collocated HDF5 conversion.

The imported upstream equations are force free. Every snapshot remains gated
against emission use until separate field-release checks pass.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time
import tomllib

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent
UPSTREAM_COMMIT = "9986522d08b4d887e5ef9a53e12cecae4e376d22"
BIN = ROOT / "third_party/coffee-portable/bin/coffee-3d"
RUNS = ROOT / "third_party/coffee-runs"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_config(n=32, end=0.1, alpha=75., extent=3., radius=.5,
                output_count=2, courant=.3, omega=1., ko=.05):
    """Node-centred odd N avoids placing the dipole singularity on a node."""
    h = 2 * extent / n
    steps = max(1, math.ceil(end / (courant * h)))
    interval = max(1, math.ceil(steps / output_count))
    steps = interval * output_count
    dt = end / steps
    pml = max(3, round(.75 / h))
    params = dict(dt=dt, max_steps=steps, data_interval=interval,
        snapshot_interval=steps, periodic_boundary=[False]*3,
        N=[n]*3, guard=[3]*3, lower=[-extent]*3, size=[2*extent]*3,
        nodes=[1., 1., 1.], shift_ghost=0, downsample=1,
        calc_current=True, clean_ep=True, check_egb=True,
        use_edotb_damping=False, pml=[pml]*3, pmllen=pml, sigpml=.5,
        divB_clean=True, ch2=1., tau=2., KOeps=ko, KO_geometry=False,
        problem=1, radius=radius, omega=omega, b0=1.,
        p1=math.sin(math.radians(alpha)), p2=0., p3=math.cos(math.radians(alpha)),
        q11=0., q12=0., q13=0., q22=0., q23=0.)
    return params


def write_config(path, config):
    # JSON values form a compatible TOML subset for these primitive fields.
    Path(path).write_text("\n".join(f"{k} = {json.dumps(v)}" for k,v in config.items())+"\n")


def convert(run_dir, source=None):
    run_dir = Path(run_dir).resolve()
    config = tomllib.loads((run_dir / "config.toml").read_text())
    manifest_path=run_dir/'run-result.json'
    manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    omega = float(config['omega'])
    if omega <= 0:
        raise ValueError('Canonical R_LC conversion requires positive omega')
    if any(config.get(k,0) != 0 for k in ('q11','q12','q13','q22','q23')):
        raise ValueError('Unit dipole normalization is implemented only for zero quadrupole')
    mu = config['b0'] * np.linalg.norm([config['p1'],config['p2'],config['p3']])
    scale = mu * omega**3
    paths = [Path(source)] if source else sorted((run_dir/'Data').glob('fld.*.h5'))
    output = []
    for path in paths:
        step = int(path.stem.split('.')[-1]) * config['data_interval']
        with h5py.File(path) as f:
            e = np.stack([np.array(f[k]).transpose(2,1,0) for k in ('Ex','Ey','Ez')]) / scale
            b = np.stack([np.array(f[k]).transpose(2,1,0) for k in ('Bx','By','Bz')]) / scale
        if not np.all(np.isfinite(e)) or not np.all(np.isfinite(b)):
            raise ValueError(f'Non-finite field in {path}')
        axes = [omega*(config['lower'][j]+np.arange(e.shape[j+1])*
                config['size'][j]/config['N'][j]*config['downsample']) for j in range(3)]
        meta = dict(solver='COFFEE CPU EZ', source_commit=UPSTREAM_COMMIT,
            source_url='https://github.com/fizban007/CoffeeGPU',
            coordinate_units='R_LC', time_units='Omega*t',
            field_units='mu/R_LC^3; common Gaussian c=1 E and B scale',
            staggering='All E/B components collocated at Stagger(111) nodes; no averaging',
            hdf5_storage_order='z,y,x; transposed to component,x,y,z',
            export_precision='Upstream fld HDF5 casts double evolution to float32',
            output_time_method='filename index * data_interval * dt; output precedes evolution',
            dipole_phase_zero='mu in +x,+z plane; positive rotation about +z',
            magnetic_phase_rad=omega*step*config['dt'],
            inclination_deg=float(np.degrees(np.arccos(config['p3']/np.linalg.norm(
                [config['p1'],config['p2'],config['p3']])))),
            star_radius=float(config['radius']*omega),
            normalization_factor=float(scale), raw_hdf5_sha256=digest(path),
            config=config, emission_release_allowed=False,
            caveat='Numerical force-free field. No physical E_parallel or radiation population is inferred.')
        patch = ROOT/'third_party/coffee-portable.patch'
        if patch.exists(): meta['portable_patch_sha256']=digest(patch)
        meta['boundary_model']=manifest.get('boundary_model','upstream')
        meta['binary_sha256']=manifest.get('binary_sha256')
        if 'boundary_patch_sha256' in manifest:
            meta['boundary_patch_sha256']=manifest['boundary_patch_sha256']
            meta['boundary_model_caveat']='Numerical boundary experiment, not the unmodified upstream boundary model'
        for key in ('repair_patch_sha256','repair_model','physical_glm_damping_time',
                    'physical_sponge_peak_rate','absorber_patch_sha256','absorber_model',
                    'source_file_sha256'):
            if key in manifest:
                meta[key]=manifest[key]
        target=run_dir/f'snapshot-{step:06d}.npz'
        np.savez_compressed(target,e=e,b=b,x=axes[0],y=axes[1],z=axes[2],
                            time=omega*step*config['dt'],metadata=json.dumps(meta))
        output.append(str(target))
    return output


def run(run_dir, config, binary=BIN, ranks=1, boundary_model='upstream'):
    run_dir=Path(run_dir).resolve()
    run_dir.mkdir(parents=True,exist_ok=True)
    if (run_dir/'Data').exists():
        raise FileExistsError(f'Existing data in {run_dir}; choose a new run directory')
    write_config(run_dir/'config.toml',config)
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OMPI_MCA_btl='self,tcp',OMPI_MCA_pml='ob1')
    command=[str(binary)]
    if ranks>1:
        command=[str(ROOT/'third_party/coffee-tools/env/bin/mpirun'),'-n',str(ranks)]+command
    started=time.perf_counter()
    with (run_dir/'run.log').open('w') as log:
        process=subprocess.run(command,cwd=run_dir,env=env,stdout=log,stderr=subprocess.STDOUT)
    elapsed=time.perf_counter()-started
    result=dict(wall_seconds=elapsed,returncode=process.returncode,
                command=command,config=config,source_commit=UPSTREAM_COMMIT,
                binary_sha256=digest(binary),emission_release_allowed=False)
    result['boundary_model']=boundary_model
    if boundary_model=='surface_anchored_experiment':
        result['boundary_patch_sha256']=digest(ROOT/'third_party/coffee-surface-anchored.patch')
    if boundary_model in ('curl_conductor_experiment','curl_ohmic_experiment','curl_radiative_experiment'):
        result['repair_patch_sha256']=digest(ROOT/'third_party/coffee-curl-conductor.patch')
        result['repair_model']='D4-curl finite C4 dipole initialization; E-only rotating conductor; uniform GLM; collocated background'
        result['physical_glm_damping_time']=config['tau']*config['dt']
        result['physical_sponge_peak_rate']=config['sigpml']/config['dt']
    if boundary_model=='curl_ohmic_experiment':
        result['absorber_patch_sha256']=digest(ROOT/'third_party/coffee-curl-ohmic.patch')
        result['repair_model']+='; E-only Ohmic outer absorber (no direct magnetic damping)'
        result['absorber_model']='Cubic spatial electric conductivity profile; exponential full-step E/P update; B changed only by Faraday curl, GLM gradient and KO dissipation'
        result['source_file_sha256']={f:digest(ROOT/'third_party/coffee-curl-ohmic'/f) for f in ('src/algorithms/field_solver_EZ.cpp','src/algorithms/damping_boundary.cpp','src/algorithms/regularized_dipole.h','src/main-pulsar3d.cpp','src/main-verification.cpp')}
    if boundary_model=='curl_radiative_experiment':
        result['absorber_patch_sha256']=digest(ROOT/'third_party/coffee-curl-radiative.patch')
        result['repair_model']+='; electric incoming-characteristic penalty layer (no direct magnetic damping)'
        result['absorber_model']='Tangential E relaxed toward -n cross (B-B_initial), normal E preserved; cubic face-weighted layer; full-step exponential splitting; B changed only by Faraday curl, uniform GLM gradient and KO dissipation'
        result['source_file_sha256']={f:digest(ROOT/'third_party/coffee-curl-radiative'/f) for f in ('src/algorithms/field_solver_EZ.cpp','src/algorithms/damping_boundary.cpp','src/algorithms/regularized_dipole.h','src/main-pulsar3d.cpp','src/main-verification.cpp')}
    (run_dir/'run-result.json').write_text(json.dumps(result,indent=2)+'\n')
    if process.returncode: raise RuntimeError(f'COFFEE failed; see {run_dir}/run.log')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('pulsar')
    q.add_argument('--n',type=int,default=33)
    q.add_argument('--end',type=float,default=.2)
    q.add_argument('--alpha',type=float,default=75.)
    q.add_argument('--extent',type=float,default=3.)
    q.add_argument('--radius',type=float,default=.5)
    q.add_argument('--output-count',type=int,default=2)
    q.add_argument('--ranks',type=int,default=1)
    q.add_argument('--binary',type=Path,default=BIN)
    q.add_argument('--boundary-model',choices=('upstream','surface_anchored_experiment'),default='upstream')
    q.add_argument('--run-dir',type=Path,required=True)
    q=sub.add_parser('convert')
    q.add_argument('run_dir',type=Path)
    q.add_argument('--source',type=Path)
    a=p.parse_args()
    if a.command=='convert':
        print(json.dumps(convert(a.run_dir,a.source),indent=2));return
    if a.n<9 or a.n%2==0:
        raise ValueError('Use odd N >=9 to avoid the stellar point singularity')
    config=make_config(a.n,a.end,a.alpha,a.extent,a.radius,a.output_count)
    result=run(a.run_dir,config,binary=a.binary.resolve(),ranks=a.ranks,boundary_model=a.boundary_model)
    result['canonical_snapshots']=convert(a.run_dir)
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
