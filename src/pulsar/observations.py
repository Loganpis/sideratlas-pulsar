"""Read genuine Fermi 3PC data and compare shape without claiming a spectral fit."""
from pathlib import Path
import numpy as np
from astropy.io import fits


def load_vela(path, rebin=4):
    with fits.open(path, checksum=True) as hdul:
        table=hdul['GAMMA_LC'].data
        lo=np.array(table['Ph_Min'],dtype=float)
        hi=np.array(table['Ph_Max'],dtype=float)
        y=np.array(table['GT100_WtCnt'],dtype=float)
        err=np.array(table['Unc_GT100_WtCnt'],dtype=float)
    use=(lo>=0)&(lo<1)&(hi<=1+1e-6)
    second=(lo>=1)&(lo<2)
    repeats=(int(second.sum())==int(use.sum()) and
             np.array_equal(y[second],y[use]) and np.array_equal(err[second],err[use]))
    if not repeats:
        raise ValueError("Unexpected Fermi display-cycle structure: inspect before fitting")
    lo,hi,y,err=lo[use],hi[use],y[use],err[use]
    if not np.all(err>0) or not np.all(np.isfinite(y)) or len(y)%rebin:
        raise ValueError("Invalid observations or rebin factor")
    if not np.allclose(hi[:-1],lo[1:],atol=1e-7):
        raise ValueError("Noncontiguous phase bins")
    original=len(y)
    lo=lo.reshape(-1,rebin)[:,0]; hi=hi.reshape(-1,rebin)[:,-1]
    y=y.reshape(-1,rebin).sum(axis=1)
    err=np.sqrt((err.reshape(-1,rebin)**2).sum(axis=1))
    return {'phase_lo':lo,'phase_hi':hi,'phase':(lo+hi)/2,'counts':y,'error':err,
            'band':'Fermi-LAT >100 MeV; GT100_WtCnt',
            'original_independent_bins':original,'fit_bins':len(y),
            'duplicate_display_cycle_removed':bool(repeats),'rebin_factor':rebin}


def bin_average(phases, profile, lo, hi, shift, samples=16):
    sample=lo[:,None]+(hi-lo)[:,None]*(np.arange(samples)+0.5)/samples
    return np.interp((sample-shift)%1,phases,profile,period=1.).mean(axis=1)


def nonnegative_scale_and_offset(model,y,err):
    """Exact weighted least squares with A>=0, C>=0 for A*model+C."""
    design=np.column_stack((model,np.ones_like(model)))
    unconstrained=np.linalg.lstsq(design/err[:,None],y/err,rcond=None)[0]
    w=1/err**2
    candidates=[np.array([max(0,np.sum(w*model*y)/np.sum(w*model*model)),0.]),
                np.array([0.,max(0,np.sum(w*y)/np.sum(w))])]
    if np.all(unconstrained>=0):
        candidates.append(unconstrained)
    return min(candidates,key=lambda coef:np.sum(((design@coef-y)/err)**2))


def compare_to_vela(phases,profile,data,n_shifts=4096):
    best=None
    for shift in np.arange(n_shifts)/n_shifts:
        m=bin_average(phases,profile,data['phase_lo'],data['phase_hi'],shift)
        coef=nonnegative_scale_and_offset(m,data['counts'],data['error'])
        predicted=coef[0]*m+coef[1]
        chi2=float(np.sum(((predicted-data['counts'])/data['error'])**2))
        if best is None or chi2<best['chi_squared']:
            best={'chi_squared':chi2,'phase_shift_cycles':float(shift),
                  'scale':float(coef[0]),'constant_component':float(coef[1]),
                  'prediction':predicted}
    best['degrees_of_freedom']=len(data['counts'])-3
    best['reduced_chi_squared']=best['chi_squared']/best['degrees_of_freedom']
    best['residual_sigma']=(data['counts']-best['prediction'])/data['error']
    best['statistical_caveat']='Shape diagnostic using supplied errors; no instrument/spectral forward model or model-error covariance'
    best['fitted_parameters']='Phase origin, nonnegative normalization and nonnegative constant only; physical parameters fixed'
    return best
