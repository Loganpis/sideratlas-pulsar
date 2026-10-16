#!/usr/bin/env python3
"""Exact 45-degree, source-free Gaussian pulse through the native backend.

n=(1,1,0)/sqrt(2), E=(0,0,f), B=n+(n_y,-n_x,0)f and
f=.2 exp[-((n.r-t+1.5)/.35)^2]. Equal x/y spacings make the sampled initial
magnetic field D4-solenoidal exactly. Compare energy after t=4 in |x|,|y|<1
with an otherwise identical wider domain whose boundary is causally excluded
there through t=10. This tests a finite cube and its actual guard treatment.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from run_coffee import ROOT,make_config,run,digest
from run_coffee_verification import read

BIN=ROOT/'third_party/coffee-curl-radiative/bin/coffee-verification'

def d4(v,h,axis):
    out=np.zeros_like(v);s=[slice(None)]*3;s[axis]=slice(2,-2)
    v0=[]
    for off in (-2,-1,1,2):
        ss=[slice(None)]*3;ss[axis]=slice(2+off,(-2+off) or None)
        v0.append(v[tuple(ss)])
    out[tuple(s)]=(v0[0]-8*v0[1]+8*v0[2]-v0[3])/(12*h)
    return out

def check(folder,n,extent,rate,cleaning):
    end=10.;h=2*extent/n;count=24
    interval=math.ceil(end/(.3*h)/count);steps=interval*count;dt=end/steps
    p=make_config(n,end)
    p.update(N=[n,n,8],size=[2*extent,2*extent,2.],lower=[-extent,-extent,-1.],
        periodic_boundary=[False,False,True],pml=[round(1/h),round(1/h),0],
        pmllen=round(1/h),problem=5,dt=dt,max_steps=steps,data_interval=interval,
        KOeps=.05,divB_clean=cleaning,tau=.5/dt,sigpml=rate*dt)
    if not (folder/'run-result.json').exists():run(folder,p,binary=BIN,boundary_model='curl_radiative_experiment')
    x=-extent+np.arange(n)*h;inside=(abs(x[:,None])<1)&(abs(x[None,:])<1)
    divmask=(abs(x[:,None])<2)&(abs(x[None,:])<2)
    rows=[];c=1/np.sqrt(2.)
    for k,path in enumerate(sorted((folder/'Data').glob('fld.*.h5'))):
        f=read(path);energy=sum(f['E'+a]**2 for a in 'xyz')+(f['Bx']-c)**2+(f['By']-c)**2+f['Bz']**2
        div=d4(f['Bx'],h,0)+d4(f['By'],h,1)
        rows.append(dict(time=k*interval*dt,interior_perturbation_energy=float(np.sum(energy[inside])),
            interior_divB_rms=float(np.sqrt(np.mean(div[divmask]**2))),
            max_interior_divB=float(np.max(abs(div[divmask])))))
    incident=max(r['interior_perturbation_energy'] for r in rows if r['time']<=3)
    return dict(n=n,extent=extent,h=h,cleaning=cleaning,run_dir=str(folder),trace=rows,
        late_interior_energy_over_incident_peak=max(r['interior_perturbation_energy'] for r in rows if r['time']>=4)/incident,
        max_interior_divB_rms=max(r['interior_divB_rms'] for r in rows),
        late_window=[4.,10.],normalization='Maximum perturbation energy inside |x|,|y|<1 during 0<=t<=3')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,default=ROOT/'third_party/coffee-runs/verification-curl-oblique')
    p.add_argument('--n',type=int,default=96)
    p.add_argument('--rate',type=float,default=.5/.02755783029464731)
    p.add_argument('--no-cleaning',action='store_true')
    a=p.parse_args();a.directory.mkdir(parents=True,exist_ok=True)
    rows=[check(a.directory/f'cube-n{a.n}',a.n,4.,a.rate,not a.no_cleaning),
          check(a.directory/f'causal-reference-n{3*a.n}',3*a.n,12.,a.rate,not a.no_cleaning)]
    result=dict(description=__doc__,binary_sha256=digest(BIN),results=rows,
        emission_release_allowed=False,limitations=['Single 45-degree pulse with guide field; not a current-sheet or global wind certification',
            'Boundary reflection and corner effects are both included in the small-domain late perturbation energy'])
    (a.directory/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in row.items() if k!='trace'} for row in rows],indent=2))
