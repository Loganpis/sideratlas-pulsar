#!/usr/bin/env python3
"""Bounded analytical checks of the actual compiled COFFEE CPU field update.

Guide-field wave: E=(0,f,0), B=(1,0,f), f(x,t)=A sin[2pi(x-t)].
The exact solution has J=0, E.B=0, B^2-E^2=1 and propagates at c=1.
Pulse reflection uses the same solution with a Gaussian f and measures the
left-moving characteristic (Ey-Bz)/2 after the right-moving pulse exits.
"""
import argparse
import json
import math
from pathlib import Path
import h5py
import numpy as np
from run_coffee import ROOT, make_config, run

BIN=ROOT/'third_party/coffee-portable/bin/coffee-verification'


def read(path):
    with h5py.File(path) as f:
        return {k:np.array(f[k],float).transpose(2,1,0) for k in ('Ex','Ey','Ez','Bx','By','Bz')}


def execute(folder,params):
    if not (folder/'run-result.json').exists():
        run(folder,params,binary=BIN)
    return json.loads((folder/'run-result.json').read_text())


def wave(n,directory):
    # Small transverse directions still exercise the real 3-D MPI guard exchange.
    steps=math.ceil(n/.3)
    p=make_config(n,1.)
    p.update(N=[n,8,8],size=[1.,1.,1.],lower=[0.,0.,0.],
        periodic_boundary=[True]*3,pml=[0,0,0],pmllen=1,problem=0,
        dt=1/steps,max_steps=steps,data_interval=steps,KOeps=0.,divB_clean=True)
    folder=directory/f'guide-wave-n{n}'
    result=execute(folder,p)
    f=read(folder/'Data/fld.00001.h5')
    exact=.2*np.sin(2*np.pi*(np.arange(n)/n-1))[:,None,None]
    err=np.sqrt(np.mean((f['Ey']-exact)**2))/np.sqrt(np.mean(exact**2))
    return dict(n=n,relative_l2_Ey=float(err),relative_l2_Bz=float(
        np.sqrt(np.mean((f['Bz']-exact)**2))/np.sqrt(np.mean(exact**2))),
        max_E_dot_B=float(np.max(np.abs(f['Ex']*f['Bx']+f['Ey']*f['By']+f['Ez']*f['Bz']))),
        wall_seconds=result['wall_seconds'],run_dir=str(folder))


def absorber(n,directory,ko=.05):
    h=8/n;steps=math.ceil(7.5/(.3*h));interval=math.ceil(steps/30);steps=interval*30
    p=make_config(n,7.5)
    p.update(N=[n,8,8],size=[8.,2.,2.],lower=[-4.,-1.,-1.],
        periodic_boundary=[False,True,True],pml=[round(1/h),0,0],
        pmllen=round(1/h),problem=3,dt=7.5/steps,max_steps=steps,
        data_interval=interval,KOeps=ko,divB_clean=True)
    folder=directory/f'absorber-n{n}'
    result=execute(folder,p)
    x=-4+np.arange(n)*h;inside=np.abs(x)<2
    initial=read(folder/'Data/fld.00000.h5')
    incident=float(np.sum(((initial['Ey']+initial['Bz'])/2)**2))
    rows=[]
    for k,path in enumerate(sorted((folder/'Data').glob('fld.*.h5'))):
        f=read(path)
        left=(f['Ey']-f['Bz'])/2;right=(f['Ey']+f['Bz'])/2
        rows.append(dict(time=k*interval*p['dt'],
            left_energy_fraction=float(np.sum(left[inside]**2)/incident),
            right_energy_fraction=float(np.sum(right[inside]**2)/incident)))
    return dict(n=n,maximum_returned_left_energy_fraction=max(
        r['left_energy_fraction'] for r in rows if r['time']>=4.5),
        interior_final_wave_energy_fraction=rows[-1]['left_energy_fraction']+rows[-1]['right_energy_fraction'],
        trace=rows,wall_seconds=result['wall_seconds'],run_dir=str(folder))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,default=ROOT/'third_party/coffee-runs/verification')
    a=p.parse_args();a.directory.mkdir(parents=True,exist_ok=True)
    result=dict(description=__doc__,source_equations_changed=False,
                wave=[wave(n,a.directory) for n in (16,32,64)],
                absorber=[absorber(n,a.directory) for n in (96,192)])
    for left,right in zip(result['wave'][:-1],result['wave'][1:]):
        right['measured_order']=math.log(left['relative_l2_Ey']/right['relative_l2_Ey'],2)
    result['analytic_wave_pass']=all(r['relative_l2_Ey']<.01 for r in result['wave']) and all(
        r.get('measured_order',4)>3.5 for r in result['wave'])
    result['planar_absorber_pass']=all(r['maximum_returned_left_energy_fraction']<.001 for r in result['absorber'])
    result['limitations']=['Planar guide-field wave is not a pulsar solution.',
        'Normal-incidence pulse with constant guide field does not establish oblique/current-sheet absorption.',
        'No global field convergence or emission release follows from these checks.']
    (a.directory/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('description','absorber')},indent=2))
    print(json.dumps([ {k:v for k,v in r.items() if k!='trace'} for r in result['absorber']],indent=2))


if __name__=='__main__':main()
