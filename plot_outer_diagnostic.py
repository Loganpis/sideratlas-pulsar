#!/usr/bin/env python3
"""Inspect completed outer-transport diagnostics without implying validation."""
from pathlib import Path
import os,sys,json,argparse
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pulsar.multiband import load_bands
from pulsar.spectral_observations import load_vela_spectrum

def main():
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);args=parser.parse_args()
    out=args.output;report=json.loads((out/'results.json').read_text())
    data=load_bands(ROOT/'data/J0835-4510_3PC_data.fits')
    sed=load_vela_spectrum(ROOT/'data/spectra/3PC_Catalog+SEDs_20230803.fits')
    plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':11})
    fig,axes=plt.subplots(2,3,figsize=(14,8.3))
    for ax,obs,pred in zip(axes.ravel()[:4],data,report['profiles']['bands']):
        maximum=obs['counts'].max()
        ax.errorbar(obs['phase'],obs['counts']/maximum,yerr=obs['error']/maximum,fmt='.',ms=2.5,color='#647583',label='Vela 3PC')
        ax.plot(obs['phase'],np.array(pred['prediction'])/maximum,color='#a94c24',label='Diagnostic model')
        ax.set(xlim=(0,1),xlabel='Rotation phase',ylabel='Relative weighted counts',title=obs['band'])
    axes[0,0].legend(fontsize=9,frameon=False)
    ax=axes.ravel()[4];x=sed['energy_centers_gev'];y=sed['flux_ph_cm2_s'];use=~sed['is_upper_limit']
    ax.errorbar(x[use],y[use],yerr=np.array([sed['error_lower_ph_cm2_s'][use],sed['error_upper_ph_cm2_s'][use]]),fmt='o',color='#647583',label='3PC integrated bin flux')
    limit=sed['upper_limit_ph_cm2_s'][~use]
    ax.errorbar(x[~use],limit,yerr=.5*limit,uplims=True,fmt='v',color='#647583')
    ax.plot(x,report['spectrum']['prediction_ph_cm2_s'],'o-',color='#a94c24',label='One spectral normalization')
    ax.set(xscale='log',yscale='log',xlabel='Photon energy (GeV)',ylabel='Bin photon flux (cm$^{-2}$ s$^{-1}$)',title='Spectrum across eight energy bins')
    ax.legend(fontsize=8,frameon=False)
    ax.text(.98,.98,f'$\\chi^2$ = {report["spectrum"]["descriptive_chi_squared"]:,.0f} (6 descriptive dof)',
            transform=ax.transAxes,fontsize=9,va='top',ha='right')
    ax=axes.ravel()[5]
    with np.load(out/'emission.npz') as d:
        active=d['dt']>0;radius=np.linalg.norm(d['position'],axis=-1)
        emission=d['radiated_erg']*d['seed_weight'][None,:]
        ax.hist(radius[active],bins=np.linspace(0,2,50),weights=emission[active],color='#2c7790')
    ax.axvline(1,color='black',lw=1,ls=':')
    ax.set(xlabel='Spherical radius ($R_{LC}$)',ylabel='Energy per injected primary (erg)',title='Where curvature energy is emitted')
    for ax in axes.ravel():ax.spines[['top','right']].set_visible(False)
    grid=report['protocol']['snapshot'].get('config',{}).get('N',['?'])[0]
    fig.suptitle('Outer-magnetosphere calculation — development diagnostic\n'
                 f'{grid}³ field · {report["segment_shape"][1]} primaries · {report["selected_seeds"]} reach the observer bin',fontsize=17)
    fig.text(.5,.015,'Unreleased numerical field · prescribed primary loading and acceleration · flat pulse response · no independent Vela validation',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.035,1,.95));fig.savefig(out/'comparison.png',dpi=180);fig.savefig(out/'comparison.pdf');plt.close(fig)
    print(str(out/'comparison.png'))

if __name__=='__main__':main()
