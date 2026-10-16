"""Conservative finite-bin integration of trajectory emission.

Each transport cell already owns a band-integrated photon yield. Reconstruct
its sky latitude and unwrapped arrival phase linearly in emission time, then
integrate the fractional cell residence in each observer and phase bin. No
Gaussian smoothing, phase-dependent weights or band renormalization is used.
The reconstruction, orbit step, seed grid and observer width must converge
separately. This does not recover unresolved field or particle physics.
"""
import numpy as np


def rebin_phase_averages(edges,profiles,lower,upper,shift=0.):
    """Conservatively integrate periodic phase-bin averages over new bins.

    Values are piecewise constant within their declared edges, never treated
    as nodal samples. Positive shift translates the model to later phase.
    """
    edges=np.asarray(edges,float);p=np.asarray(profiles,float)
    lo,hi=np.broadcast_arrays(np.asarray(lower,float),np.asarray(upper,float))
    if (edges.ndim!=1 or len(edges)<2 or edges[0]!=0 or edges[-1]!=1 or np.any(np.diff(edges)<=0)
            or p.shape[-1]!=len(edges)-1 or not all(np.all(np.isfinite(a)) for a in (edges,p,lo,hi))
            or np.any(hi<=lo) or not np.isfinite(shift)):
        raise ValueError('Finite periodic phase averages and ordered target bins required')
    prefix=np.concatenate((np.zeros(p.shape[:-1]+(1,)),np.cumsum(p*np.diff(edges),axis=-1)),axis=-1)
    def integral(x):
        period=np.floor(x);r=x-period
        i=np.minimum(np.searchsorted(edges,r,side='right')-1,len(edges)-2)
        return prefix[...,-1,None]*period+prefix[...,i]+p[...,i]*(r-edges[i])
    return (integral(hi-shift)-integral(lo-shift))/(hi-lo)


def reconstruct_sky_endpoints(bundle):
    """Return cell endpoint mu/phase from active midpoint samples.

    Second-order finite differences reconstruct a slope inside each cell;
    endpoint mu is clipped to the physical [-1,1] interval, with the count
    reported. Phase is unwrapped independently on each continuous stream.
    A single-cell stream has constant sky coordinates. Cells retain their
    own integrated yield, even if neighboring reconstructions do not join.
    """
    dt=np.asarray(bundle['dt'],float);time=np.asarray(bundle['time'],float)
    phase=np.asarray(bundle['phase'],float);direction=np.asarray(bundle['direction'],float)
    if dt.ndim!=2 or time.shape!=dt.shape or phase.shape!=dt.shape or direction.shape!=dt.shape+(3,):
        raise ValueError('Expected segment/seed midpoint arrays')
    mu=np.zeros(dt.shape+(2,));phi=np.zeros_like(mu);clipped=0;largest_jump=0.;unresolved=0
    for seed in range(dt.shape[1]):
        ids=np.flatnonzero(dt[:,seed]>0)
        if not len(ids):continue
        if not np.array_equal(ids,np.arange(len(ids))):raise ValueError('Active stream cells must form a prefix')
        tm=time[ids,seed];m=direction[ids,seed,2];p=phase[ids,seed]
        if not np.all(np.isfinite([tm,m,p])) or np.any(np.diff(tm)<=0):
            raise ValueError('Finite increasing active times and nonaxial sky phases required')
        bad=np.zeros(len(ids),bool)
        if len(ids)>1:
            jump=(np.diff(p)+.5)%1-.5
            largest_jump=max(largest_jump,float(np.max(abs(jump))))
            # Near a spin-axis ray phase can change rapidly while the ray is
            # invisible to this observer. Mark adjacent cells undefined, and
            # reject them only if their latitude interval enters the sky bin.
            bad[:-1]|=abs(jump)>=.25;bad[1:]|=abs(jump)>=.25
            p=np.unwrap(p*2*np.pi)/(2*np.pi)
            order=2 if len(ids)>2 else 1
            dm=np.gradient(m,tm,edge_order=order);dp=np.gradient(p,tm,edge_order=order)
        else:dm=np.zeros(1);dp=np.zeros(1)
        h=dt[ids,seed]/2
        endpoints=m[:,None]+dm[:,None]*h[:,None]*np.array([-1.,1.])
        clipped+=int(np.count_nonzero(abs(endpoints)>1))
        mu[ids,seed]=np.clip(endpoints,-1,1)
        phi[ids,seed]=p[:,None]+dp[:,None]*h[:,None]*np.array([-1.,1.])
        bad|=abs(phi[ids,seed,1]-phi[ids,seed,0])>=.25
        phi[ids[bad],seed]=np.nan;unresolved+=int(bad.sum())
    return mu,phi,{'mu_endpoints_limited':clipped,'max_midpoint_phase_jump_cycles':largest_jump,
                   'undefined_phase_cells':unresolved}


def project_segment_emission(bundle,photon_yield,*,observer_deg=65.,view_width_deg=2.,
                             nphase=100,sky_endpoints=None):
    """Integrate top-hat sky and phase bins; return bin *averages* per sr.

    sky_endpoints may supply exact (mu,phase) cell endpoints for verification.
    Otherwise they are reconstructed from trajectory midpoints. Phase endpoint
    differences must already be unwrapped. Partial seed-weight sums are valid
    for batching: weights are preserved and never renormalized here.
    Emissivity per band is constant within each cell, so fractions of each
    integrated yield equal fractions of cell time. Refining transport cells
    tests this approximation independently of the exact bin intersection.
    """
    y=np.asarray(photon_yield,float);dt=np.asarray(bundle['dt'],float)
    w=np.asarray(bundle['seed_weight'],float)
    if dt.ndim!=2 or y.ndim!=3 or y.shape[:2]!=dt.shape or w.shape!=(dt.shape[1],):
        raise ValueError('Inconsistent emission arrays')
    if not all(np.all(np.isfinite(a)) for a in (y,dt,w)) or np.any(y<0) or np.any(dt<0) or np.any(w<0):
        raise ValueError('Finite nonnegative emission, durations and weights required')
    if (not np.isfinite(observer_deg) or not np.isfinite(view_width_deg) or view_width_deg<=0
            or observer_deg-view_width_deg/2<0 or observer_deg+view_width_deg/2>180):
        raise ValueError('Finite sky bin must lie within 0 to 180 degrees')
    if isinstance(nphase,bool) or not isinstance(nphase,int) or nphase<1:raise ValueError('Positive integer phase bins required')
    if sky_endpoints is None:mu,phase,diagnostic=reconstruct_sky_endpoints(bundle)
    else:
        mu,phase=(np.asarray(a,float) for a in sky_endpoints)
        diagnostic={'mu_endpoints_limited':0,'max_midpoint_phase_jump_cycles':None}
    if (mu.shape!=dt.shape+(2,) or phase.shape!=mu.shape or not np.all(np.isfinite(mu))
            or np.any(abs(mu)>1)):
        raise ValueError('Finite physical sky endpoints required for every cell')
    low=np.cos(np.deg2rad(observer_deg+view_width_deg/2));high=np.cos(np.deg2rad(observer_deg-view_width_deg/2))
    solid=2*np.pi*(high-low);dm=mu[...,1]-mu[...,0]
    moving=dm!=0;lower=np.zeros_like(dt);upper=np.ones_like(dt)
    x1=np.divide(low-mu[...,0],dm,out=np.zeros_like(dt),where=moving)
    x2=np.divide(high-mu[...,0],dm,out=np.ones_like(dt),where=moving)
    lower[moving]=np.maximum(0,np.minimum(x1,x2))[moving]
    upper[moving]=np.minimum(1,np.maximum(x1,x2))[moving]
    stationary_inside=(mu[...,0]>=low)&((mu[...,0]<high)|((high==1)&(mu[...,0]==1)))
    upper[~moving&~stationary_inside]=0
    fraction=np.maximum(0,upper-lower)*(dt>0)
    selected=fraction>0
    if np.any(~np.isfinite(phase[selected])):
        raise ValueError('Unresolved arrival-phase cell reaches observer bin; reduce orbit step')
    ids=np.nonzero(fraction>0);a=phase[...,0][ids];delta=(phase[...,1]-phase[...,0])[ids]
    begin=a+lower[ids]*delta;end=a+upper[ids]*delta
    pmin=np.minimum(begin,end)*nphase;pmax=np.maximum(begin,end)*nphase
    values=y[ids]*w[ids[1],None];frac=fraction[ids]
    histogram=np.zeros((nphase,y.shape[-1]));constant=abs(pmax-pmin)<1e-13
    if np.any(constant):np.add.at(histogram,np.floor(pmin[constant]).astype(np.int64)%nphase,values[constant]*frac[constant,None])
    variable=~constant
    if np.any(variable):
        lo=pmin[variable];hi=pmax[variable];v=values[variable];mass=frac[variable]
        first=np.floor(lo).astype(np.int64);last=np.ceil(hi).astype(np.int64)-1
        crossings=int(np.max(last-first))+1
        if crossings>4*nphase:raise ValueError('Unresolved multi-rotation cell; reduce orbit step')
        for offset in range(crossings):
            bins=first+offset;overlap=np.maximum(0,np.minimum(hi,bins+1)-np.maximum(lo,bins))
            amount=mass*overlap/(hi-lo);used=amount>0
            np.add.at(histogram,bins[used]%nphase,v[used]*amount[used,None])
    profile=histogram.T*nphase/solid
    return {'phase':(np.arange(nphase)+.5)/nphase,'phase_edges':np.arange(nphase+1)/nphase,
            'profiles':profile,'phase_averaged_yield_per_sr':profile.mean(axis=1),
            'solid_angle_sr':float(solid),'selected_segments':int(len(ids[0])),
            'selected_seeds':int(len(np.unique(ids[1]))),'selected_cell_fraction_sum':float(fraction.sum()),
            'seed_weight_sum':float(w.sum()),'observer_deg':float(observer_deg),
            'view_width_deg':float(view_width_deg),'physical_phase_blur':0.,
            'profile_convention':'finite-bin phase averages, piecewise-constant cell emissivity',
            'reconstruction':diagnostic,'absolute_population_known':False}
