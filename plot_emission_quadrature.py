#!/usr/bin/env python3
"""Render conservative pulse averages and catalog spectrum without a validation claim."""
from pathlib import Path
import argparse, json, os, sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pulsar.emission_quadrature import rebin_phase_averages
from pulsar.observations import nonnegative_scale_and_offset
from pulsar.multiband import load_bands
from pulsar.spectral_observations import load_vela_spectrum


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run',type=Path)
    ap.add_argument('--width',type=float,default=2.)
    args=ap.parse_args();out=args.run
    report=json.loads((out/'results.json').read_text())
    row=report['width_scenarios'][str(args.width)]
    data=load_bands(ROOT/'data/J0835-4510_3PC_data.fits')
    sed=load_vela_spectrum(ROOT/'data/spectra/3PC_Catalog+SEDs_20230803.fits')
    with np.load(out/'profiles.npz') as d:
        edges=d['phase_edges'];profiles=d[f'profiles_width_{args.width:g}']
    plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':11})
    fig,axes=plt.subplots(2,3,figsize=(14,8.1))
    for ax,obs,profile in zip(axes.ravel()[:4],data,profiles[:4]):
        m=rebin_phase_averages(edges,profile,obs['phase_lo'],obs['phase_hi'],
                               row['pulse_shape']['phase_shift_cycles'])
        m=m/m.max() if m.max()>0 else m
        a,c=nonnegative_scale_and_offset(m,obs['counts'],obs['error'])
        top=obs['counts'].max()
        ax.errorbar(obs['phase'],obs['counts']/top,yerr=obs['error']/top,
                    fmt='.',color='#526f84',ms=3,label='Vela / Fermi 3PC')
        ax.stairs((a*m+c)/top,np.r_[obs['phase_lo'],obs['phase_hi'][-1]],
                  color='#b4562c',lw=1.4,label='Calculated / finite-bin average')
        ax.set(xlim=(0,1),xlabel='Rotation phase',ylabel='Relative weighted counts',title=obs['band'])
    axes[0,0].legend(frameon=False,fontsize=8)
    ax=axes.ravel()[4];x=sed['energy_centers_gev'];use=~sed['is_upper_limit']
    ax.errorbar(x[use],sed['flux_ph_cm2_s'][use],
        yerr=[sed['error_lower_ph_cm2_s'][use],sed['error_upper_ph_cm2_s'][use]],
        fmt='o',ms=5,color='#526f84',label='Catalog integrated flux')
    lim=sed['upper_limit_ph_cm2_s'][~use]
    ax.errorbar(x[~use],lim,yerr=lim*.5,uplims=True,fmt='v',color='#526f84')
    ax.plot(x,row['spectrum']['prediction_ph_cm2_s'],'o-',color='#b4562c',ms=4,
            label='One population normalization')
    ax.set(xscale='log',yscale='log',ylim=(1e-13,3e-5),xlabel='Photon energy (GeV)',
           ylabel='Bin photon flux (cm$^{-2}$ s$^{-1}$)',title='Phase-averaged spectrum')
    ax.legend(frameon=False,fontsize=8,loc='lower left')
    ax.text(.97,.97,f'Descriptive $\\chi^2$ = {row["spectrum"]["descriptive_chi_squared"]:,.0f}',
            transform=ax.transAxes,ha='right',va='top',fontsize=9)
    ax=axes.ravel()[5];ax.axis('off')
    p=report['protocol'];cfg=p['snapshot'].get('config',{});n=cfg.get('N',['?'])[0]
    lines=['Calculation and limits', '', f'Field grid: {n}³',
        f'Injected paths: {report["seed_count"]:,}',
        f'Observer bin: {p["observer_deg"]:g}° ± {args.width/2:g}°',
        f'Path step: {p["transport"]["dt"]:g} $R_{{LC}}$',
        'Finite observer/phase crossing integration',
        'Published, prescribed particle acceleration', '',
        'Numerical field and population remain approximate',
        'All observations used for development',
        'No independent Vela validation']
    ax.text(0,.98,'\n'.join(lines),transform=ax.transAxes,va='top',linespacing=1.65)
    for ax in axes.ravel()[:5]:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('@sideratlas · NO.002 — Vela reconstruction development\nCalculated emission against measured gamma-ray data',fontsize=17)
    fig.text(.5,.018,'Diagnostic model · illustrative plot colors · flat pulse-response approximation · no final Reel release',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.035,1,.94))
    target=out/f'comparison-width-{args.width:g}'
    fig.savefig(target.with_suffix('.png'),dpi=180);fig.savefig(target.with_suffix('.pdf'));plt.close(fig)
    print(target.with_suffix('.png'))


if __name__=='__main__':main()
