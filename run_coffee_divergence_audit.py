#!/usr/bin/env python3
"""Independent D4/GLM/source diagnostics; does not modify simulation fields."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
AS = (0., -.4178904745, -1.192151694643, -1.697784692471, -1.514183444257)
BS = (.1496590219993, .3792103129999, .8229550293869, .6994504559488, .1530572479681)


def derivative(f, h, axis, order=4):
    coefficients = {2: ((-1, -.5), (1, .5)),
                    4: ((-2, 1/12), (-1, -2/3), (1, 2/3), (2, -1/12)),
                    6: ((-3, -1/60), (-2, 3/20), (-1, -3/4),
                        (1, 3/4), (2, -3/20), (3, 1/60))}[order]
    return sum(weight*np.roll(f, -offset, axis) for offset, weight in coefficients)/h


def divergence(b, h, order=4):
    return sum(derivative(b[j], h, j, order) for j in range(3))


def curl(b, h):
    return np.stack([derivative(b[2], h, 1)-derivative(b[1], h, 2),
                     derivative(b[0], h, 2)-derivative(b[2], h, 0),
                     derivative(b[1], h, 0)-derivative(b[0], h, 1)])


def regular_dipole_potential(xyz, moment, radius=.5):
    """C1 radial potential, exact dipole outside radius, finite at origin."""
    r = np.linalg.norm(xyz, axis=0)
    # A=(m cross r) f(r); f,f' match r^-3 at R and are regular at r=0.
    f = np.where(r < radius, 2.5/radius**3-1.5*r*r/radius**5,
                 1/np.maximum(r, 1e-30)**3)
    return np.moveaxis(np.cross(moment, np.moveaxis(xyz, 0, -1)), -1, 0)*f


def glm_mode(n=64, dt=.01, duration=2., tau_steps=2., physical_tau=None, mode=1):
    """Actual grid D4 and five-stage low-storage RK versus exact Fourier ODE.

    The physical damping time is tau_steps*dt unless explicitly overridden.
    The source-free periodic test has s(r)=1 and disables boundary/KO terms.
    """
    h=2*np.pi/n; x=np.arange(n)*h
    steps=int(round(duration/dt)); dt=duration/steps
    gamma=1/(tau_steps*dt if physical_tau is None else physical_tau)
    b=np.cos(mode*x); p=np.zeros(n); db=np.zeros(n); dp=np.zeros(n)
    for _ in range(steps):
        for a, c in zip(AS, BS):
            db=a*db-dt*derivative(p,h,0)
            dp=a*dp-dt*(derivative(b,h,0)+gamma*p)
            b=b+c*db; p=p+c*dp
    k=(8*np.sin(mode*h)-np.sin(2*mode*h))/(6*h)
    matrix=np.array([[0.,-1j*k],[-1j*k,-gamma]],complex)
    amplitude=expm(duration*matrix)@np.array([1.,0.])
    exact_b=np.real(amplitude[0]*np.exp(1j*mode*x))
    exact_p=np.real(amplitude[1]*np.exp(1j*mode*x))
    return dict(n=n,dt=dt,steps=steps,duration=duration,mode=mode,
                tau_physical=1/gamma,modified_wavenumber=k,
                max_abs_B_error=float(np.max(np.abs(b-exact_b))),
                max_abs_P_error=float(np.max(np.abs(p-exact_p))),
                divergence_amplitude_ratio=float(abs(amplitude[0])),
                measured_divergence_amplitude_ratio=float(np.linalg.norm(derivative(b,h,0))/np.linalg.norm(derivative(np.cos(mode*x),h,0))),
                characteristic_eigenvalues_real=np.linalg.eigvals(matrix).real.tolist(),
                characteristic_eigenvalues_imag=np.linalg.eigvals(matrix).imag.tolist())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot_audit(path):
    import h5py
    with np.load(path,allow_pickle=False) as d:
        b=np.asarray(d['b'],float);e=np.asarray(d['e'],float)
        axes=tuple(d[k] for k in ('x','y','z'));time=float(d['time'])
        metadata=json.loads(str(d['metadata']))
    h=axes[0][1]-axes[0][0];xyz=np.stack(np.meshgrid(*axes,indexing='ij'))
    r=np.linalg.norm(xyz,axis=0);bmag=np.linalg.norm(b,axis=0)
    div=divergence(b,h);rot=curl(b,h);score=np.linalg.norm(rot,axis=0)/np.maximum(bmag,1e-30)
    alpha=np.deg2rad(metadata['inclination_deg']);phase=metadata['magnetic_phase_rad']
    m=np.array([np.sin(alpha)*np.cos(phase),np.sin(alpha)*np.sin(phase),np.cos(alpha)])
    nn=xyz/r
    reference=(3*nn*np.einsum('i,ijkl->jkl',m,nn)-m[:,None,None,None])/r**3
    refdiv=divergence(reference,h)
    # Exported diagnostic divB is from before the last RK/source/KO operations.
    # Recompute from final B, and keep that old field separate below.
    data=path.parent/'Data/fld.00012.h5'
    with h5py.File(data) as d:
        p=np.asarray(d['P'],float).transpose(2,1,0)
        saved_div=np.asarray(d['divB'],float).transpose(2,1,0)
    if p.shape!=r.shape:raise ValueError('HDF5 P grid differs from collocated snapshot')
    s=.5*(1-np.tanh((r-3)/.5))
    gradp=np.stack([derivative(p,h,j) for j in range(3)])
    clean_rhs=-divergence(s*gradp,h)
    cfg=metadata['config'];dt=cfg['dt'];tau=cfg['tau'];gamma=1/(tau*dt)
    w=np.where(r<.75,.5*(1-np.tanh((r-.5)/(.5*h))),0.)
    w=np.where(r<.25,1.,w)
    delta=w*(reference-b);delta_div=divergence(delta,h)
    target_div=divergence(curl(regular_dipole_potential(xyz,m),h),h)
    valid=np.ones(r.shape,bool)
    for j in range(3):valid &= (xyz[j]>=axes[j][3])&(xyz[j]<=axes[j][-4])
    region=(r>.75)&(r<2)&valid
    denom=np.sum(div[region]**2)
    rows=[]
    for lo,hi in ((0,.25),(.25,.35),(.35,.5),(.5,.6),(.6,.75),(.75,1),(1,1.5),(1.5,2),(2,2.25),(2.25,2.6)):
        mask=(r>lo)&(r<hi)&valid
        if not mask.any():continue
        norm=lambda f:float(np.sqrt(np.mean(f[mask]**2)))
        rows.append(dict(rmin=lo,rmax=hi,cells=int(mask.sum()),divB_rms=norm(div),B_rms=norm(bmag),
            divB_over_B_rms=norm(div)/norm(bmag),reference_dipole_divB_rms=norm(refdiv),
            regular_curl_target_divB_rms=norm(target_div),P_rms=norm(p),
            scalar_damping_rms=gamma*norm(p),scalar_divergence_drive_rms=norm(div),
            scalar_damping_vs_drive_cosine=float(np.sum(p[mask]*div[mask])/(np.linalg.norm(p[mask])*np.linalg.norm(div[mask]))),
            cleaning_divB_rhs_rms=norm(clean_rhs),reapplied_source_divB_increment_rms=norm(delta_div),
            saved_stage_divB_rms=norm(saved_div)))
    outer=(r>1)&(r<2);threshold=float(np.quantile(score[outer],.9))
    sheet=outer&(score>=threshold);smooth=outer&~sheet
    def stat(mask):
        return dict(cells=int(mask.sum()),divB_over_B_rms=float(np.sqrt(np.sum(div[mask]**2)/np.sum(bmag[mask]**2))),
                    share_exterior_divB_square=float(np.sum(div[mask]**2)/denom))
    return dict(path=str(path.relative_to(ROOT)),snapshot_sha256=sha(path),hdf5_sha256=sha(data),
                time=time,h=h,dt=dt,tau_steps=tau,tau_physical=tau*dt,
                exterior=stat(region),inner_exterior=stat(region&(r<1)),
                sheet_proxy_threshold=threshold,sheet_proxy_top10percent=stat(sheet),
                remaining_outer=stat(smooth),shells=rows,
                derivative_comparison={str(order):float(np.sqrt(np.sum(divergence(b,h,order)[region]**2)/np.sum(bmag[region]**2))) for order in (2,4,6)},
                div_curl_E_max_interior=float(np.max(np.abs(divergence(curl(e,h),h)[region]))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/field-benchmarks/coffee-divergence-audit.json')
    parser.add_argument('--analytic-only',action='store_true')
    a=parser.parse_args()
    paths=[ROOT/'third_party/coffee-runs'/f'pulsar-a75-n{n}-3rot-surface-anchored'/f'snapshot-{steps:06d}.npz' for n,steps in ((65,684),(97,1020))]
    sources=['run_coffee_divergence_audit.py','third_party/coffee-surface-anchored/src/algorithms/field_solver_EZ.cpp',
             'third_party/coffee-surface-anchored/src/algorithms/finite_diff_simd.h','third_party/coffee-surface-anchored/src/sim_params.h']
    out=dict(audited_at_utc=datetime.now(timezone.utc).isoformat(),software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
             source_sha256={p:sha(ROOT/p) for p in sources},
             scope='Independent analytic and spatial diagnostics; no release thresholds, clipping or production modifications.',
             glm_modes=[],snapshots=[])
    for n,steps in ((65,684),(97,1020)):
        dt=6*np.pi/steps
        for physical_tau in (None,2*(6*np.pi/684),1.):
            result=glm_mode(n=n,dt=dt,duration=6*np.pi,physical_tau=physical_tau)
            result['damping_policy']='2 timesteps' if physical_tau is None else 'fixed physical time'
            out['glm_modes'].append(result)
    if not a.analytic_only:
        for path in paths:
            print('Auditing',path.name,path.parent.name,flush=True)
            out['snapshots'].append(snapshot_audit(path))
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print('Saved',a.output,flush=True)


if __name__=='__main__':main()
