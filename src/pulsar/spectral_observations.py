"""Measured 3PC Vela spectrum, kept separate from the pulse-profile data.

The eight integrated photon-flux bins can constrain a model's relative
spectrum with ONE population normalization. nuFnu is retained for display;
it is model-dependent within each broad catalog energy bin. Supplied errors
are asymmetric. Upper limits are reported separately, never fitted as
detections or interpreted as an exact likelihood without the likelihood scan.
"""
from pathlib import Path
import hashlib
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize_scalar

CATALOG_URL = ('https://fermi.gsfc.nasa.gov/ssc/data/access/lat/'
               '3rd_PSR_catalog/3PC_Catalog%2BSEDs_20230803.fits')
ENERGY_EDGES_GEV = np.array([.05,.1,.3,1.,3.,10.,30.,100.,1000.])
MEV_ERG = 1.602176634e-6


def load_vela_spectrum(path):
    """Join PSRJ to the exact 4FGL counterpart; retain catalog units/errors."""
    path=Path(path)
    with fits.open(path, checksum=True) as h:
        psrs=h['PULSARS_BIGFILE'].data
        selected=np.flatnonzero(np.char.strip(psrs['PSRJ'])=='J0835-4510')
        if len(selected)!=1: raise ValueError('Require exactly one Vela timing row')
        p=psrs[selected[0]]
        name=str(p['FGLname']).strip()
        spectra=h['LAT_Point_Source_Catalog'].data
        selected=np.flatnonzero(np.char.strip(spectra['Source_Name'])==name)
        if len(selected)!=1: raise ValueError('Require exactly one matching Vela spectrum')
        row=spectra[selected[0]]
        flux=np.array(row['Flux_Band'],float)
        unc=np.array(row['Unc_Flux_Band'],float)
        significance=np.array(row['Sqrt_TS_Band'],float)
        nufnu=np.array(row['nuFnu_Band'],float)*MEV_ERG
        plec={k:float(row[k+'_bfr']) for k in ('PLEC_Flux_Density',
             'PLEC_IndexS','Pivot_Energy','PLEC_ExpfactorS','PLEC_Exp_Index')}
        period=float(p['P0']); pdot=-float(p['F1'])/float(p['F0'])**2
    if flux.shape!=(8,) or unc.shape!=(8,2) or nufnu.shape!=(8,):
        raise ValueError('Unexpected catalog energy bin format')
    upper=np.isnan(unc[:,0]) | (significance<2)
    lower_error=-unc[:,0]; upper_error=unc[:,1]
    if (np.any(~np.isfinite(flux)) or np.any(flux<=0) or
        np.any(~np.isfinite(upper_error)) or np.any(upper_error<=0) or
        np.any(~np.isfinite(lower_error[~upper])) or np.any(lower_error[~upper]<=0)):
        raise ValueError('Invalid catalog flux or asymmetric errors')
    return dict(pulsar='J0835-4510',source_name=name,
        energy_edges_gev=ENERGY_EDGES_GEV.copy(),
        energy_centers_gev=np.sqrt(ENERGY_EDGES_GEV[:-1]*ENERGY_EDGES_GEV[1:]),
        flux_ph_cm2_s=flux,error_lower_ph_cm2_s=lower_error,
        error_upper_ph_cm2_s=upper_error,is_upper_limit=upper,
        upper_limit_ph_cm2_s=np.where(upper,flux+2*upper_error,np.nan),
        nufnu_erg_cm2_s=nufnu,nufnu_error_lower=nufnu*lower_error/flux,
        nufnu_error_upper=nufnu*upper_error/flux,plec_bfree=plec,
        period_s=period,period_derivative=pdot,
        provenance={'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                    'url':CATALOG_URL,'data_role':'development',
                    'observable':'phase-averaged catalog source spectrum; not phase-resolved pulsed-only flux',
                    'errors':'catalog statistical errors; no full bin covariance or likelihood scans',
                    'upper_limit_convention':'flux + 2 upper errors, matching official 3PC plotter'})


def catalog_reference_sed(energy_gev,data):
    """Official fitted PLEC4 curve as a DATA reference, never a prediction.

    Catalog IndexS is negative; dN/dE uses IndexS+d/b as its power. Input
    energy GeV, catalog density per MeV, output E^2 dN/dE in erg/cm2/s.
    """
    e=np.asarray(energy_gev,float)*1000
    if np.any(~np.isfinite(e)) or np.any(e<=0):raise ValueError('Positive finite energies required')
    p=data['plec_bfree'];b=p['PLEC_Exp_Index'];d=p['PLEC_ExpfactorS']
    x=e/p['Pivot_Energy']
    dnde=p['PLEC_Flux_Density']*np.exp((p['PLEC_IndexS']+d/b)*np.log(x)+(d/b**2)*(1-x**b))
    return e**2*dnde*MEV_ERG


def compare_integrated_spectrum(model_flux,data,*,fit_normalization=True):
    """Descriptive asymmetric-error comparison of integrated flux bins.

    One nonnegative normalization for the entire spectrum; no per-band or
    peak-specific weights. When fit_normalization=False, model_flux must be
    physical photon flux at Earth. A zero model is valid and fails naturally.
    Detections alone set the normalization; upper-limit exceedances remain
    separate diagnostics. This is not a calibrated likelihood or p-value.
    """
    m=np.asarray(model_flux,float);y=np.asarray(data['flux_ph_cm2_s'],float)
    use=~np.asarray(data['is_upper_limit'],bool)
    if m.shape!=y.shape or np.any(~np.isfinite(m)) or np.any(m<0):
        raise ValueError('Finite nonnegative model in every catalog energy bin required')
    if np.count_nonzero(use)<2: raise ValueError('At least two spectral detections required')
    low=np.asarray(data['error_lower_ph_cm2_s'],float)[use]
    high=np.asarray(data['error_upper_ph_cm2_s'],float)[use]
    if np.any(~np.isfinite([low,high])) or np.any(low<=0) or np.any(high<=0):
        raise ValueError('Positive detection errors required')
    maximum=np.max(m)
    shape=m/maximum if maximum>0 else m.copy()
    def score(amplitude):
        delta=amplitude*shape[use]-y[use]
        error=np.where(delta>=0,high,low)
        return float(np.sum((delta/error)**2))
    if fit_normalization and maximum>0 and np.any(shape[use]>0):
        # Normalize flux units before the scalar minimization so absolute
        # optimizer tolerance cannot swallow astronomical ~1e-6 flux values.
        reference=float(np.max(y[use])); positive=shape[use]>0
        bound=float(np.max(y[use][positive]/shape[use][positive])/reference)
        result=minimize_scalar(lambda a:score(a*reference),bounds=(0,bound*1.01),
                               method='bounded',options={'xatol':1e-12})
        amplitude=float(result.x*reference)
        if score(0)<score(amplitude):amplitude=0.
        normalization=amplitude/maximum
    else:
        normalization=1. if not fit_normalization else 0.
        amplitude=normalization*maximum
    prediction=m*normalization
    residual=(y[use]-prediction[use])/np.where(prediction[use]>=y[use],high,low)
    limits=np.asarray(data['upper_limit_ph_cm2_s'],float)
    upper_ratio=np.full(m.shape,np.nan)
    upper_ratio[~use]=prediction[~use]/limits[~use]
    dof=int(use.sum())-int(fit_normalization and maximum>0 and np.any(shape[use]>0))
    return dict(normalization=float(normalization),prediction_ph_cm2_s=prediction,
        detection_indices=np.flatnonzero(use),residual_sigma=residual,
        descriptive_chi_squared=float(residual@residual),descriptive_dof=dof,
        normalization_fitted=bool(fit_normalization),upper_limit_ratio=upper_ratio,
        upper_limit_exceeded=bool(np.any(upper_ratio[~use]>1)),
        caveat='Development diagnostic; catalog bin errors without full covariance/likelihood; upper limits excluded from fitted statistic; source spectrum is not a phase-resolved pulsar-only separation')
