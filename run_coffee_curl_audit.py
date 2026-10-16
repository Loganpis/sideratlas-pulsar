#!/usr/bin/env python3
"""Measure discrete divergence and source moment drift for curl-conductor runs."""
import argparse
import json
from pathlib import Path
import numpy as np
from run_coffee_assessment import load,quality,derivative,norms,surface_integrals

def potential(points,moment,radius):
    r2=np.sum(points**2,axis=-1);u=r2/radius**2-1
    poly=(1+u*(-1.5+u*(1.875+u*(-2.1875+u*2.4609375))))/radius**3
    f=np.where(r2<radius**2,poly,np.maximum(r2,1e-30)**-1.5)
    return np.cross(moment,points)*f[...,None]

def source_basis(points,h,radius):
    basis=[]
    for moment in np.eye(3):
        ds=[]
        for axis in range(3):
            out=np.zeros_like(points)
            for offset,weight in ((-2,1),(-1,-8),(1,8),(2,-1)):
                moved=points.copy();moved[:,axis]+=offset*h
                out+=weight*potential(moved,moment,radius)/(12*h)
            ds.append(out)
        basis.append(np.stack((ds[1][:,2]-ds[2][:,1],ds[2][:,0]-ds[0][:,2],ds[0][:,1]-ds[1][:,0]),axis=-1).ravel())
    return np.stack(basis,axis=-1)

def audit(directory):
    rows=[];basis=None
    for path in sorted(directory.glob('snapshot-*.npz')):
        d=load(path);h=d['axes'][0][1]-d['axes'][0][0]
        xyz=np.stack(np.meshgrid(*d['axes'],indexing='ij'));r=norms(xyz)
        mask=r<.35
        if basis is None:basis=source_basis(xyz[:,mask].T,h,d['metadata']['star_radius'])
        values=d['b'][:,mask].T.ravel();moment=np.linalg.lstsq(basis,values,rcond=None)[0]
        magnitude=np.linalg.norm(moment);phase=np.arctan2(moment[1],moment[0])
        phase_error=np.angle(np.exp(1j*(phase-d['time'])))
        div=sum(derivative(d['b'][j],h,j) for j in range(3));bm=norms(d['b'])
        shells={}
        for lo,hi in ((0,.35),(.35,.75),(.75,1),(1,1.5),(1.5,2),(2,2.5)):
            m=(r>lo)&(r<hi)
            shells[f'{lo:g}<r<{hi:g}']=dict(div_rms=float(np.sqrt(np.mean(div[m]**2))),
                div_rms_over_B_RLC=float(np.sqrt(np.sum(div[m]**2)/np.sum(bm[m]**2))))
        rows.append(dict(snapshot=str(path),time=d['time'],source_fit=dict(
            radius_max=.35,moment=moment.tolist(),magnitude=float(magnitude),
            phase_error_deg=float(np.degrees(phase_error)),
            inclination_deg=float(np.degrees(np.arccos(moment[2]/magnitude))),
            relative_fit_residual=float(np.linalg.norm(basis@moment-values)/np.linalg.norm(values))),
            quality=quality(d),divergence_shells=shells,
            luminosity=[surface_integrals(d,rr) for rr in (1.,1.5,2.)]))
    result=dict(run_dir=str(directory),emission_release_allowed=False,rows=rows,
        source_fit_method='Least-squares fit to the actual finite C4 dipole D4-curl basis at native grid nodes r<0.35; fitted moment need not imply a globally dipolar field.')
    (directory/'curl-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);a=p.parse_args()
    result=audit(a.directory)
    print(json.dumps([dict(time=row['time'],source=row['source_fit'],quality=row['quality']) for row in result['rows']],indent=2))
