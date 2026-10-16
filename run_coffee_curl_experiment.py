#!/usr/bin/env python3
"""Bounded finite-source, curl-only conducting-star COFFEE experiment."""
import argparse
import json
from pathlib import Path
from run_coffee import ROOT,make_config,run,convert

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=65)
    p.add_argument('--end',type=float,default=6.283185307179586)
    p.add_argument('--output-count',type=int,default=4)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--extent',type=float,default=3.)
    p.add_argument('--tau-physical',type=float,default=.5)
    a=p.parse_args()
    config=make_config(a.n,a.end,75.,a.extent,.5,a.output_count)
    config['tau']=a.tau_physical/config['dt']
    # Match the old65 sponge's physical rate, rather than its damping per step.
    config['sigpml']=(.5/.02755783029464731)*config['dt']
    result=run(a.directory,config,
        binary=ROOT/'third_party/coffee-curl-conductor/bin/coffee-3d',
        boundary_model='curl_conductor_experiment')
    result['canonical_snapshots']=convert(a.directory)
    print(json.dumps(result,indent=2))
