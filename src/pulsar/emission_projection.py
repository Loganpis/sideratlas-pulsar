"""Project lab-frame photon yields without renormalizing energy bands.

The steady rotating-pattern approximation replicates each injected trajectory
around the spin axis. Its sky phase already includes propagation delay.
Counts are divided by the sky bin's solid angle, not by their per-band mean;
this retains the phase-averaged spectrum and the relative peak amplitudes.
"""
import numpy as np


def project_emission(bundle,photon_yield,*,observer_deg=65.,view_width_deg=1.,nphase=256):
    """Finite sky-bin average, CIC phase deposition, no extra phase smoothing.

    Output profiles are photons per injected primary per steradian as a
    function of phase. Their mean is the phase-averaged yield per steradian.
    Multiply by a physically justified primary injection rate/distance² for
    flux at Earth. Normalized seed weights alone do not provide that rate.
    Angular bin width and phase bins are numerical resolutions, not physical
    emission-cone widths, and must be refined with the trajectory quadrature.
    """
    y=np.asarray(photon_yield,float);phase=np.asarray(bundle['phase'],float)
    direction=np.asarray(bundle['direction'],float);duration=np.asarray(bundle['dt'],float)
    weight=np.asarray(bundle['seed_weight'],float)
    if (phase.ndim!=2 or direction.shape!=phase.shape+(3,) or duration.shape!=phase.shape or
        weight.shape!=(phase.shape[1],) or y.shape[:2]!=phase.shape or y.ndim!=3):
        raise ValueError('Inconsistent segment/seed/band arrays')
    if not all(np.isfinite(a).all() for a in (direction,duration,weight,y)) or np.any(duration<0) or np.any(y<0) or np.any(weight<0) or not np.isclose(weight.sum(),1,rtol=1e-12):
        raise ValueError('Finite nonnegative emission and normalized injection required')
    if not (np.isfinite(observer_deg) and np.isfinite(view_width_deg) and
            view_width_deg>0 and observer_deg-view_width_deg/2>=0 and observer_deg+view_width_deg/2<=180):
        raise ValueError('Sky bin must lie within 0–180 degrees')
    if isinstance(nphase,bool) or not isinstance(nphase,int) or nphase<8:raise ValueError('At least eight phase bins required')
    mu_low=np.cos(np.deg2rad(observer_deg+view_width_deg/2))
    mu_high=np.cos(np.deg2rad(observer_deg-view_width_deg/2))
    mu=direction[...,2]
    selected=(duration>0)&(mu>=mu_low)&(mu<mu_high)
    if mu_high==1:selected|=(duration>0)&(mu==1)
    if np.any(~np.isfinite(phase[selected])):
        raise ValueError('Selected axial/undefined sky phase needs a dedicated polar-observer treatment')
    if np.any(np.abs(np.linalg.norm(direction[duration>0],axis=-1)-1)>1e-7):
        raise ValueError('Active emission directions must be unit lab velocities')
    solid_angle=2*np.pi*(mu_high-mu_low)
    p=(np.where(selected,phase,0)%1)*nphase;left=np.floor(p).astype(int)%nphase;fraction=p-np.floor(p)
    profiles=[]
    weights=weight[None,:]*selected
    for band in range(y.shape[-1]):
        w=y[...,band]*weights*nphase/solid_angle
        value=np.bincount(left.ravel(),weights=(w*(1-fraction)).ravel(),minlength=nphase)
        value+=np.bincount(((left+1)%nphase).ravel(),weights=(w*fraction).ravel(),minlength=nphase)
        profiles.append(value)
    profiles=np.array(profiles)
    return {'phase':np.arange(nphase)/nphase,'profiles':profiles,
            'phase_averaged_yield_per_sr':profiles.mean(axis=1),'solid_angle_sr':float(solid_angle),
            'selected_segments':int(selected.sum()),
            'selected_seeds':int(np.count_nonzero(np.any(selected,axis=0))),
            'observer_deg':float(observer_deg),'view_width_deg':float(view_width_deg),
            'physical_phase_blur':0.,'absolute_population_known':False}
