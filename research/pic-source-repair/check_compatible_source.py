#!/usr/bin/env python3
"""Independent algebra checks of the compatible rotating-source electric field.

These test the discrete source, not convergence of the plasma magnetosphere.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np


def phase(t, spinup=.8):
    s=np.clip(t/spinup,0,1)
    return spinup*(s**3-.5*s**4)+max(0,t-spinup)


def omega(t, spinup=.8):
    s=np.clip(t/spinup,0,1)
    return s*s*(3-2*s)


def potential(x,t,R=.5):
    angle=phase(t);chi=np.radians(75);mu=R**3/2
    m=mu*np.array([np.sin(chi)*np.cos(angle),np.sin(chi)*np.sin(angle),np.cos(chi)])
    r2=np.sum(x*x,axis=-1)
    f=np.empty_like(r2);outside=r2>=R*R
    f[outside]=r2[outside]**-1.5
    f[~outside]=(2.5-1.5*r2[~outside]/R**2)/R**3
    return np.cross(m,x)*f[...,None]


def phi(x,t):
    A=potential(x,t)
    return omega(t)*(-x[...,1]*A[...,0]+x[...,0]*A[...,1])


def compatible_e(x,c,t,h,dt):
    d=np.eye(3)[c]*h/2
    return -(potential(x,t+dt/2)[...,c]-potential(x,t-dt/2)[...,c])/dt-(phi(x+d,t)-phi(x-d,t))/h


def curl_vector_function(function,points,h):
    result=np.empty_like(points)
    for c in range(3):
        a,b=(c+2)%3,(c+1)%3
        db=np.eye(3)[b]*h/2;da=np.eye(3)[a]*h/2
        result[...,c]=(function(points+db)[...,a]-function(points-db)[...,a])/h-(function(points+da)[...,b]-function(points-da)[...,b])/h
    return result


def checks():
    rng=np.random.default_rng(1741)
    points=rng.uniform(-.25,.25,(2048,3))
    errors=[]
    for t in (0,.24,.8,1.9,6.65):
        h=5/48;dt=.35*h/np.sqrt(3)
        Ep=lambda x: np.stack([compatible_e(x,c,t,h,dt) for c in range(3)],axis=-1)
        before=curl_vector_function(lambda x:potential(x,t-dt/2),points,h)
        after=curl_vector_function(lambda x:potential(x,t+dt/2),points,h)
        update=before-dt*curl_vector_function(Ep,points,h)
        error=float(np.max(np.abs(update-after)))
        assert error<2e-14,error
        errors.append(error)
    # The compatible source tends to the ideal co-rotation electric field.
    # Independent high-accuracy derivatives reconstruct B at each sample.
    t=2.;ref_h=1e-5
    B=curl_vector_function(lambda x:potential(x,t),points,ref_h)
    v=omega(t)*np.cross(np.array([0.,0.,1.]),points)
    exact=-np.cross(v,B)
    relative=[]
    for h in (5/48,5/96,5/192):
        dt=.35*h/np.sqrt(3)
        E=np.stack([compatible_e(points,c,t,h,dt) for c in range(3)],axis=-1)
        relative.append(float(np.linalg.norm(E-exact)/np.linalg.norm(exact)))
    assert relative[1]<.27*relative[0] and relative[2]<.27*relative[1],relative
    return {'scope':'discrete source algebra; not plasma or optical validation',
            'native_one_step_reference_dt':.35*(5/48)/np.sqrt(3),
            'max_B_update_absolute_error_each_time':errors,
            'ideal_conductor_E_relative_L2_at_dx_5_over_48_96_192':relative,
            'second_order_ratios':[relative[0]/relative[1],relative[1]/relative[2]],
            'checks_passed':7}

if __name__=='__main__':
    result=checks()
    target=Path(__file__).resolve().parents[2]/'outputs/pic-source-repair-algebra.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
