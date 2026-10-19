"""Compare two actually calculated emergent optical tables, with no scaling."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def read(path):
    z=np.load(path)
    return z['wavelength_m'],z['mu'],z['i_lambda_si']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('coarse',type=Path);p.add_argument('fine',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--minimum-nm',type=float,default=330)
    p.add_argument('--maximum-nm',type=float,default=580)
    a=p.parse_args()
    w1,m1,i1=read(a.coarse);w2,m2,i2=read(a.fine)
    wave=np.linspace(a.minimum_nm*1e-9,a.maximum_nm*1e-9,251)
    mu=np.unique(np.r_[0,.0001,.001,.005,.01,.02,.04,np.linspace(.05,1,96)])
    xx,yy=np.meshgrid(wave,mu,indexing='ij');points=np.c_[xx.ravel(),yy.ravel()]
    interp=lambda w,m,i:RegularGridInterpolator((w,m),i,bounds_error=True)(points).reshape(xx.shape)
    first=interp(w1,m1,i1);second=interp(w2,m2,i2)
    relative=second/np.maximum(first,1e-300)-1
    band1=np.trapz(first,wave,axis=0);band2=np.trapz(second,wave,axis=0)
    flux1=2*np.pi*np.trapz(band1*mu,mu);flux2=2*np.pi*np.trapz(band2*mu,mu)
    nongrazing=mu>=.02
    report=dict(coarse_file=str(a.coarse),fine_file=str(a.fine),
        coarse_sha256=hashlib.sha256(a.coarse.read_bytes()).hexdigest(),fine_sha256=hashlib.sha256(a.fine.read_bytes()).hexdigest(),
        emitted_band_nm=[a.minimum_nm,a.maximum_nm],
        comparison='absolute SI radiance; no exposure adjustment, fitted scale or normalization',
        max_pointwise_relative_change_non_grazing=float(np.max(abs(relative[:,nongrazing]))),
        median_pointwise_relative_change_non_grazing=float(np.median(abs(relative[:,nongrazing]))),
        max_band_radiance_relative_change_non_grazing=float(np.max(abs(band2[nongrazing]/band1[nongrazing]-1))),
        hemispheric_band_flux_coarse_w_m2=float(flux1),hemispheric_band_flux_fine_w_m2=float(flux2),
        hemispheric_band_flux_relative_change=float(flux2/flux1-1),
        one_percent_optical_refinement_pass=bool(np.max(abs(relative[:,nongrazing]))<.01 and abs(flux2/flux1-1)<.01),
        excluded_grazing_area_fraction=.02**2,
        scientific_validation_pass=False,
        angular_band_rows=[dict(mu=float(m),coarse=float(f),fine=float(g),relative_change=float(g/f-1)) for m,f,g in zip(mu,band1,band2)])
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='angular_band_rows'},indent=2))


if __name__=='__main__':main()
