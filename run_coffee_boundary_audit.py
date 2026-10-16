#!/usr/bin/env python3
"""Audit imposed stellar boundary and compare a deliberately changed boundary model."""
import argparse,json
from pathlib import Path
import numpy as np
from run_coffee_assessment import load,sample,comparison,norms,quality,surface_integrals

def audit(path):
 d=load(path); alpha=np.deg2rad(d['metadata']['inclination_deg']);phase=d['metadata']['magnetic_phase_rad']
 m=np.array([np.sin(alpha)*np.cos(phase),np.sin(alpha)*np.sin(phase),np.cos(alpha)])
 rows=[]
 for r in (.3,.4,.5,.6,.7,.9,1.1,1.5,2.):
  e=sample(d,'e',r*m);b=sample(d,'b',r*m)
  b=np.reshape(b,3);e=np.reshape(e,3)
  rows.append(dict(radius=r,B=b.tolist(),B_radial=float(b@m),
   B_radial_over_instantaneous_dipole=float((b@m)/(2/r**3)),
   B_radial_cosine=float((b@m)/np.linalg.norm(b)),E=e.tolist()))
 return dict(path=str(path),time=d['time'],metadata=d['metadata'],pole_ray=rows,
  quality=quality(d),surfaces=[surface_integrals(d,r) for r in (1.,1.5,2.)]),d

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--original',type=Path,required=True);p.add_argument('--variant',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 before,b=audit(a.original);after,f=audit(a.variant)
 if any(not np.array_equal(x,y) for x,y in zip(b['axes'],f['axes'])):raise ValueError('Same-grid comparison required')
 r=norms(np.stack(np.meshgrid(*b['axes'],indexing='ij')))
 out=dict(classification='Numerical stellar-boundary-model sensitivity; not a portability fix',
  original=before,variant=after,emission_release_allowed=False,
  unforced_exterior_comparison=comparison(b['e'],b['b'],f['e'],f['b'],(r>.75)&(r<2.)),
  criterion='Pole radial B should keep the prescribed dipole polarity within the conducting star.',
  pole_polarity_inside_pass=all(row['B_radial']>0 for row in after['pole_ray'] if row['radius']<=.5),
  limitations=['Boundary variant changes the numerical source model and requires independent refinement.',
   'Enlarged .5R_LC star and finite sponge remain. No physical emission release follows from polarity repair.'])
 a.output.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({k:v for k,v in out.items() if k not in ('original','variant')},indent=2))
 print('Original pole ray:',before['pole_ray']);print('Variant pole ray:',after['pole_ray'])

if __name__=='__main__':main()
