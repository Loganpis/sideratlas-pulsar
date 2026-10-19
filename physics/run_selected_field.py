#!/usr/bin/env python3
"""Bounded finite-source, curl-only conducting-star / electric characteristic layer COFFEE experiment."""
import argparse
import json
from pathlib import Path
from run_coffee import ROOT,make_config,run,convert

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=97)
    p.add_argument('--end',type=float,default=18.84955592153876)
    p.add_argument('--output-count',type=int,default=12)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--extent',type=float,default=3.)
    p.add_argument('--tau-physical',type=float,default=.5)
    p.add_argument('--sponge-rate',type=float,default=.5/.02755783029464731)
    a=p.parse_args()
    config=make_config(a.n,a.end,75.,a.extent,.5,a.output_count)
    config['tau']=a.tau_physical/config['dt']
    config['divB_clean']=False
    # Match the old65 sponge's physical rate, rather than its damping per step.
    config['sigpml']=a.sponge_rate*config['dt']
    result=run(a.directory,config,
        binary=ROOT/'third_party/coffee-curl-radiative/bin/coffee-3d',
        boundary_model='curl_radiative_experiment')
    result['canonical_snapshots']=convert(a.directory)
    print(json.dumps(result,indent=2))
