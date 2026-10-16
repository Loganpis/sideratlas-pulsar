#!/usr/bin/env python3
"""Write the observed spectral constraint and its provenance; no model fitting."""
from pathlib import Path
import sys,json,os
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from pulsar.spectral_observations import load_vela_spectrum,catalog_reference_sed

def serial(obj):
    if isinstance(obj,np.ndarray):return [serial(x) for x in obj.tolist()]
    if isinstance(obj,(np.floating,float)) and not np.isfinite(obj):return None
    if isinstance(obj,np.generic):return obj.item()
    if isinstance(obj,dict):return {k:serial(v) for k,v in obj.items()}
    if isinstance(obj,list):return [serial(v) for v in obj]
    return obj

def main():
    path=ROOT/'data/spectra/3PC_Catalog+SEDs_20230803.fits'
    data=load_vela_spectrum(path)
    with fits.open(path,checksum=False) as h:
        checksum={x.name:{'CHECKSUM':int(x.verify_checksum()),'DATASUM':int(x.verify_datasum())} for x in h}
    data['provenance']['embedded_fits_verification']=checksum
    data['provenance']['download_verification']='Two separate HTTPS downloads matched SHA256 df6e610774d3de90213eff933662b30046709c85d37aa3486bc8f352bb5bab10; original upstream HDU3 CHECKSUM/DATASUM fail and are preserved, not repaired'
    output=ROOT/'outputs/spectral-reference';output.mkdir(parents=True,exist_ok=True)
    (output/'vela-spectrum.json').write_text(json.dumps(serial(data),indent=2,allow_nan=False)+'\n')
    (ROOT/'data/spectra/provenance.json').write_text(json.dumps(serial(data['provenance']),indent=2,allow_nan=False)+'\n')
    plt.rcParams.update({'font.family':'serif','mathtext.fontset':'cm','font.serif':['CMU Serif','DejaVu Serif']})
    fig,ax=plt.subplots(figsize=(8.5,5.4))
    x=data['energy_centers_gev'];y=data['nufnu_erg_cm2_s'];use=~data['is_upper_limit']
    edges=data['energy_edges_gev'];xerr=np.array([x-edges[:-1],edges[1:]-x])
    ax.errorbar(x[use],y[use],xerr=xerr[:,use],yerr=np.array([data['nufnu_error_lower'][use],data['nufnu_error_upper'][use]]),fmt='o',color='#244c70',label='3PC source spectrum')
    limit=y[~use]+2*data['nufnu_error_upper'][~use]
    ax.errorbar(x[~use],limit,xerr=xerr[:,~use],yerr=limit*.45,uplims=True,fmt='v',color='#244c70')
    e=np.geomspace(.05,1000,500)
    ax.plot(e,catalog_reference_sed(e,data),color='#ad7833',label='Catalog fitted curve (data reference)')
    ax.set(xscale='log',yscale='log',xlim=(.04,1200),ylim=(1e-14,1e-8),xlabel='Photon energy (GeV)',ylabel=r'$E^2\,dN/dE$ (erg cm$^{-2}$ s$^{-1}$)',title='Vela: an observed constraint for NO.002')
    ax.legend(frameon=False,loc='lower left',fontsize=9)
    fig.text(.12,.015,'Phase-averaged source flux • development data • no simulation curve shown',fontsize=9,color='#555555')
    fig.tight_layout(rect=(0,.035,1,1));fig.savefig(output/'vela-spectrum.png',dpi=180);fig.savefig(output/'vela-spectrum.pdf');plt.close(fig)
    (output/'README.md').write_text('''# Vela spectral reference

Seven detections from 0.05–100 GeV and a 100–1000 GeV upper limit were imported from the official Fermi 3PC catalog. The plotted curve is the catalog fit, **not a simulation result**. These are phase-averaged source measurements; they do not isolate an independent second-peak spectrum or an exact pulsed-only component.

The comparison module integrates model photon yields over the same eight energy bins. It allows one population normalization across the spectrum. Upper limits are reported separately and cannot be treated as Gaussian detections. Pulse shapes still use conditional per-band normalizations until a source-specific exposure/weighting calculation exists, so adding the two descriptive statistics is not an exact joint likelihood. All inspected measurements remain development data.

The distributed source catalog has invalid embedded CHECKSUM/DATASUM values in its spectral extension. Two independent downloads from the official URL produced the same SHA256; the original bytes and warning are preserved. Details, errors, units and original source are in `vela-spectrum.json` and `data/spectra/provenance.json`.
''')
    print(json.dumps({'output':str(output),'detections':int(use.sum()),'upper_limits':int((~use).sum()),'period_s':data['period_s']}))

if __name__=='__main__':main()
