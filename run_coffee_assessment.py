#!/usr/bin/env python3
"""Measure COFFEE pulsar resolution, periodicity and numerical constraints.

Two resolutions do not establish asymptotic convergence. This report keeps
the emission gate closed, including when some individual tolerances pass.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from run_coffee import ROOT


def load(path):
    with np.load(path,allow_pickle=False) as d:
        return dict(e=np.array(d['e'],float),b=np.array(d['b'],float),
                    axes=tuple(d[k] for k in ('x','y','z')),
                    time=float(d['time']),metadata=json.loads(str(d['metadata'])),path=str(path))


def sample(d,key,points):
    return RegularGridInterpolator(d['axes'],d[key].transpose(1,2,3,0),
                                   bounds_error=True)(points)


def norms(d):return np.sqrt(np.sum(d*d,axis=0))


def comparison(e0,b0,e1,b1,mask):
    def rel(a,b):return float(np.sqrt(np.sum((a[:,mask]-b[:,mask])**2)/np.sum(b[:,mask]**2)))
    a=b0[:,mask];b=b1[:,mask]
    cos=np.sum(a*b,axis=0)/np.maximum(norms(a)*norms(b),1e-30)
    angle=np.degrees(np.arccos(np.clip(cos,-1,1)))
    return dict(cells=int(np.sum(mask)),relative_l2_E=rel(e0,e1),relative_l2_B=rel(b0,b1),
        median_B_angle_deg=float(np.median(angle)),p90_B_angle_deg=float(np.quantile(angle,.9)))


def derivative(f,h,axis):
    return (np.roll(f,2,axis)-8*np.roll(f,1,axis)+8*np.roll(f,-1,axis)-np.roll(f,-2,axis))/(12*h)


def curls(b,h):
    return np.stack((derivative(b[2],h,1)-derivative(b[1],h,2),
                     derivative(b[0],h,2)-derivative(b[2],h,0),
                     derivative(b[1],h,0)-derivative(b[0],h,1)))


def quality(d):
    e,b=d['e'],d['b'];h=d['axes'][0][1]-d['axes'][0][0]
    xyz=np.stack(np.meshgrid(*d['axes'],indexing='ij'));r=norms(xyz)
    mask=(r>.75)&(r<2.)
    div=sum(derivative(b[j],h,j) for j in range(3))
    bmag=norms(b);emag=norms(e);eb=np.sum(e*b,axis=0)
    scaled_div=np.abs(div)/np.maximum(bmag,1e-30)
    cos=np.abs(eb)/np.maximum(emag*bmag,1e-30)
    return dict(shell=[.75,2.],relative_rms_divB_RLC=float(np.sqrt(np.mean(div[mask]**2)/np.mean(bmag[mask]**2))),
        p95_abs_divB_RLC_over_B=float(np.quantile(scaled_div[mask],.95)),
        h_times_relative_rms_divB=float(h*np.sqrt(np.mean(div[mask]**2)/np.mean(bmag[mask]**2))),
        p95_abs_EB_over_EBmag=float(np.quantile(cos[mask],.95)),
        max_abs_EB_over_EBmag=float(np.max(cos[mask])),
        E_exceeds_B_by_1e_6_fraction=float(np.mean(emag[mask]>bmag[mask]*(1+1e-6))))


def surface_integrals(d,radius):
    u,w=np.polynomial.legendre.leggauss(64);phi=np.arange(128)*2*np.pi/128
    u,phi=np.meshgrid(u,phi,indexing='ij')
    direction=np.stack((np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u),axis=-1)
    e=sample(d,'e',radius*direction);b=sample(d,'b',radius*direction)
    s=np.sum(np.cross(e,b)*direction,axis=-1)
    br=np.sum(b*direction,axis=-1)
    luminosity=radius**2*np.sum(s*w[:,None])*(2*np.pi/128)/(4*np.pi)
    flux=radius**2*np.sum(np.abs(br)*w[:,None])*(2*np.pi/128)/2
    return dict(radius=radius,luminosity_mu2_omega4_over_c3=float(luminosity),unsigned_flux_one_hemisphere=float(flux))


def rotation(points,angle):
    c,s=np.cos(angle),np.sin(angle)
    return np.stack((c*points[...,0]-s*points[...,1],s*points[...,0]+c*points[...,1],points[...,2]),axis=-1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coarse',type=Path,default=ROOT/'third_party/coffee-runs/pulsar-a75-n65-3rot-cpu-fixed')
    p.add_argument('--fine',type=Path,default=ROOT/'third_party/coffee-runs/pulsar-a75-n97-3rot-complete')
    p.add_argument('--output',type=Path,default=ROOT/'third_party/coffee-runs/assessment.json')
    a=p.parse_args()
    paths=[sorted(folder.glob('snapshot-*.npz')) for folder in (a.coarse,a.fine)]
    if any(len(files)!=13 for files in paths):
        raise ValueError('This assessment requires thirteen quarter-period snapshots from zero through three rotations')
    finals=[load(files[-1]) for files in paths]
    c,f=finals
    if any(abs(d['time']-6*np.pi)>1e-10 for d in finals):
        raise ValueError('Both runs must finish at exactly three rotations')
    for key in ('source_commit','boundary_model','inclination_deg','star_radius','coordinate_units','field_units'):
        if c['metadata'].get(key)!=f['metadata'].get(key):
            raise ValueError(f'Resolution comparison requires matching {key}')
    xyz=np.stack(np.meshgrid(*c['axes'],indexing='ij'));r=norms(xyz);points=xyz.transpose(1,2,3,0)
    fe=sample(f,'e',points).transpose(3,0,1,2);fb=sample(f,'b',points).transpose(3,0,1,2)
    masks={f'{lo:g}<r<{hi:g}':(r>lo)&(r<hi) for lo,hi in ((.75,1.),(1.,1.5),(1.5,2.),(.75,2.))}
    sheet_score=norms(curls(fb,c['axes'][0][1]-c['axes'][0][0]))/np.maximum(norms(fb),1e-30)
    outer=(r>1)&(r<2);threshold=float(np.quantile(sheet_score[outer],.9))
    masks['outer_top10percent_curlB_over_B']=outer&(sheet_score>=threshold)
    result=dict(emission_release_allowed=False,
        source_commit=c['metadata']['source_commit'],coarse_path=c['path'],fine_path=f['path'],
        final_snapshot_sha256={key:hashlib.sha256(Path(d['path']).read_bytes()).hexdigest() for key,d in zip(('coarse','fine'),finals)},
        boundary_model=c['metadata'].get('boundary_model','upstream boundary'),
        assessment_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        comparison_time=c['time'],resolution_comparison={k:comparison(c['e'],c['b'],fe,fb,m) for k,m in masks.items()},
        current_sheet_proxy='Top decile of |curlB|/|B| on fine fields interpolated to coarse nodes, 1<r<2; not a geometric sheet identification.',
        sheet_proxy_threshold=threshold,quality={},periodicity={},quarter_rotation_pattern={},surfaces={})
    rng=np.random.default_rng(20261015)
    direction=rng.normal(size=(12000,3));direction/=np.linalg.norm(direction,axis=-1)[:,None]
    positions=direction*rng.uniform(.75**3,2**3,len(direction))[:,None]**(1/3)
    for folder,d in zip((a.coarse,a.fine),finals):
        key=str(len(d['axes'][0]));snapshots=sorted(folder.glob('snapshot-*.npz'))
        result['quality'][key]=quality(d)
        result['surfaces'][key]=[surface_integrals(d,rr) for rr in (1.,1.5,2.)]
        previous=load(snapshots[8]);xyz1=np.stack(np.meshgrid(*d['axes'],indexing='ij'));r1=norms(xyz1)
        if abs(previous['time']-4*np.pi)>1e-10:raise ValueError('Missing second-rotation snapshot')
        result['periodicity'][key]=dict(from_time=previous['time'],to_time=d['time'],
            **comparison(previous['e'],previous['b'],d['e'],d['b'],(r1>.75)&(r1<2.)))
        previous_e=sample(previous,'e',positions);previous_b=sample(previous,'b',positions)
        rows=[]
        for path in snapshots[9:13]:
            next_d=load(path);angle=next_d['time']-previous['time']
            if abs(angle-(len(rows)+1)*np.pi/2)>1e-10:raise ValueError('Unexpected quarter-period snapshot time')
            moved=rotation(positions,angle)
            ne=sample(next_d,'e',moved);nb=sample(next_d,'b',moved)
            pe=rotation(previous_e,angle);pb=rotation(previous_b,angle)
            rows.append(dict(time=next_d['time'],relative_l2_E=float(np.linalg.norm(ne-pe)/np.linalg.norm(ne)),
                             relative_l2_B=float(np.linalg.norm(nb-pb)/np.linalg.norm(nb))))
        result['quarter_rotation_pattern'][key]=rows
    result['reference_spindown']=float(1+np.sin(np.deg2rad(c['metadata']['inclination_deg']))**2)
    result['limitations']=['Only two grids; no asymptotic convergence estimate.',
        'Numerically enlarged radius 0.5 R_LC; domain plus/minus 3 R_LC with cubic sponge starting about2.25.',
        'No outer-domain size sensitivity or oblique-wave sponge benchmark completed.',
        'A normal-incidence planar test is insufficient to certify pulsar sheet/outflow boundary behavior.',
        'No self-consistent dissipative E_parallel; force-free cleanup removes it.',
        'Comparisons use linear interpolation of fine snapshots; all quoted times are dimensionless Omega*t.']
    if result['boundary_model']=='curl_radiative_experiment':
        result['limitations']=[
            'Only two grids; no asymptotic convergence estimate; sheet structure remains resolution sensitive.',
            'Numerically enlarged radius 0.5 R_LC; finite cube plus/minus 3 R_LC.',
            'Finite C4 dipole potential initialized by the native D4 curl; a rotating electric conductor evolves rather than hard-resets the magnetic source.',
            'Source moment, inclination and phase must be measured separately; source diffusion is not a physical spin-down prediction.',
            'Outer layer penalizes the incoming electric characteristic without directly damping magnetic components; it is not a perfectly matched layer.',
            'A 45-degree pulse benchmark is separate from this report and does not establish domain-size independence of the pulsar wind or current sheet.',
            'No self-consistent dissipative E_parallel; force-free cleanup removes it.',
            'Comparisons use linear interpolation of fine snapshots; all quoted times are dimensionless Omega*t.']
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
