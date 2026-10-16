#!/usr/bin/env python3
"""Exact longitudinal Fourier-mode check of the compiled CPU GLM update.

B=(1+a(t) cos(kx),0,0), E=0, P=p(t) sin(kx), J=0.
a'=-k*p; p'=ch2*k*a-p/tau_physical.  The existing configuration tau
is measured in timesteps, hence tau=tau_physical/dt at every resolution.
"""
import json
import math
from pathlib import Path
import h5py
import numpy as np
from scipy.linalg import expm
from run_coffee import ROOT, make_config, run

BIN=ROOT/'third_party/coffee-curl-conductor/bin/coffee-verification'
DEST=ROOT/'third_party/coffee-runs/curl-conductor-glm-verification'

def check(n):
    end=.2; length=1.; k=2*np.pi; h=length/n; tau=.2
    steps=math.ceil(end/(.25*h)); dt=end/steps
    params=make_config(n,end)
    params.update(N=[n,8,8],lower=[0.,0.,0.],size=[1.,1.,1.],
        periodic_boundary=[True]*3,pml=[0]*3,problem=4,KOeps=0.,
        dt=dt,max_steps=steps,data_interval=steps,tau=tau/dt)
    folder=DEST/f'n{n}'
    if not (folder/'run-result.json').exists():
        run(folder,params,binary=BIN,boundary_model='curl_conductor_experiment')
    with h5py.File(folder/'Data/fld.00001.h5') as f:
        bx=np.array(f['Bx'],float).transpose(2,1,0)-1
        p=np.array(f['P'],float).transpose(2,1,0)
    def exact(kk):
        a,pamp=expm(np.array([[0.,-kk],[kk,-1/tau]])*end)@np.array([.2,0.])
        return a*np.cos(k*np.arange(n)*h)[:,None,None],pamp*np.sin(k*np.arange(n)*h)[:,None,None]
    continuous=exact(k)
    modified_k=(8*np.sin(k*h)-np.sin(2*k*h))/(6*h)
    discrete=exact(modified_k)
    return dict(n=n,dt=dt,physical_damping_time=tau,configured_tau=tau/dt,
        continuum_relative_l2=float(np.sqrt(np.mean((bx-continuous[0])**2+(p-continuous[1])**2))/.2),
        semidiscrete_relative_l2=float(np.sqrt(np.mean((bx-discrete[0])**2+(p-discrete[1])**2))/.2),
        modified_wavenumber=modified_k,run_dir=str(folder))

if __name__=='__main__':
    rows=[check(n) for n in (16,32,64)]
    for left,right in zip(rows[:-1],rows[1:]):
        right['continuum_order']=math.log(left['continuum_relative_l2']/right['continuum_relative_l2'],2)
    result=dict(test=__doc__,results=rows,
        pass_check=all(row['continuum_relative_l2']<.001 for row in rows)
        and all(row.get('continuum_order',4)>3.3 for row in rows),
        limitations=['Tests uniform-coefficient GLM only, not a pulsar or sponge boundary.'])
    DEST.mkdir(parents=True,exist_ok=True)
    (DEST/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
