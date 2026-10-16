#!/usr/bin/env python3
"""Compiled pulse reflection test of the separate curl/Ohmic CPU backend.

E=(0,f,0), B=(1,0,f) is an exact force-free wave at c=1.  The absorber
modifies E only; its physical conductivity profile is held fixed under grid
refinement.  Returned characteristic energy is measured inside |x|<2.
This normal-incidence test does not certify pulsar current-sheet absorption.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from run_coffee import ROOT, make_config, run, digest
from run_coffee_verification import read

BIN=ROOT/'third_party/coffee-curl-ohmic/bin/coffee-verification'

def absorber(n,directory,rate):
    h=8/n; steps=math.ceil(7.5/(.3*h)); interval=math.ceil(steps/40); steps=interval*40
    params=make_config(n,7.5)
    params.update(N=[n,8,8],size=[8.,2.,2.],lower=[-4.,-1.,-1.],
        periodic_boundary=[False,True,True],pml=[round(1/h),0,0],
        pmllen=round(1/h),problem=3,dt=7.5/steps,max_steps=steps,
        data_interval=interval,KOeps=.05,divB_clean=True,
        tau=.5/(7.5/steps),sigpml=rate*(7.5/steps))
    folder=directory/f'absorber-n{n}-rate{rate:g}'
    if not (folder/'run-result.json').exists():
        result=run(folder,params,binary=BIN,boundary_model='curl_ohmic_experiment')
    else: result=json.loads((folder/'run-result.json').read_text())
    x=-4+np.arange(n)*h;inside=np.abs(x)<2
    initial=read(folder/'Data/fld.00000.h5')
    incident=float(np.sum(((initial['Ey']+initial['Bz'])/2)**2))
    rows=[]
    for k,path in enumerate(sorted((folder/'Data').glob('fld.*.h5'))):
        f=read(path);left=(f['Ey']-f['Bz'])/2;right=(f['Ey']+f['Bz'])/2
        # The guide Bx stays uniform if the magnetic boundary is untouched.
        rows.append(dict(time=k*interval*params['dt'],
            left_energy_fraction=float(np.sum(left[inside]**2)/incident),
            right_energy_fraction=float(np.sum(right[inside]**2)/incident),
            max_guide_Bx_error=float(np.max(np.abs(f['Bx']-1)))))
    return dict(n=n,physical_peak_rate=rate,run_dir=str(folder),
        maximum_returned_left_energy_fraction=max(r['left_energy_fraction'] for r in rows if r['time']>=4.5),
        interior_final_wave_energy_fraction=rows[-1]['left_energy_fraction']+rows[-1]['right_energy_fraction'],
        max_guide_Bx_error=max(r['max_guide_Bx_error'] for r in rows),
        trace=rows,wall_seconds=result['wall_seconds'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,default=ROOT/'third_party/coffee-runs/verification-curl-ohmic')
    p.add_argument('--rates',type=float,nargs='+',default=[4.,10.,18.1449])
    p.add_argument('--resolutions',type=int,nargs='+',default=[96,192])
    args=p.parse_args();args.directory.mkdir(parents=True,exist_ok=True)
    result=dict(description=__doc__,binary_sha256=digest(BIN),
        results=[absorber(n,args.directory,rate) for rate in args.rates for n in args.resolutions],
        emission_release_allowed=False,
        limitations=['Normal incidence only; no independent oblique-wave or global-current-sheet certification',
                     'Electric conductivity absorber is not a perfectly matched layer'])
    (args.directory/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in row.items() if k!='trace'} for row in result['results']],indent=2))
