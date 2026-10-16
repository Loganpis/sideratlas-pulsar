#!/usr/bin/env python3
"""Show what changing global accelerating gains can and cannot repair."""
from pathlib import Path
import argparse,json,os,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pulsar.spectral_observations import load_vela_spectrum


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output',type=Path);args=ap.parse_args()
    report=json.loads((args.output/'results.json').read_text())
    data=load_vela_spectrum(ROOT/'data/spectra/3PC_Catalog+SEDs_20230803.fits')
    plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':12})
    fig,(ax,table)=plt.subplots(1,2,figsize=(12,5.6),gridspec_kw={'width_ratios':[1.2,1]})
    x=data['energy_centers_gev'];use=~data['is_upper_limit']
    ax.errorbar(x[use],data['flux_ph_cm2_s'][use],
        yerr=[data['error_lower_ph_cm2_s'][use],data['error_upper_ph_cm2_s'][use]],
        fmt='o',color='black',label='Vela / Fermi 3PC',zorder=10)
    lim=data['upper_limit_ph_cm2_s'][~use]
    ax.errorbar(x[~use],lim,yerr=lim/2,uplims=True,fmt='v',color='black')
    colors=['#7d8f9c','#bb5628','#216c81','#75638e']
    cells=[]
    for i,row in enumerate(report['trials']):
        inside=row['rate_inside_cm_inverse'];outside=row['rate_outside_cm_inverse']
        ax.plot(x,row['spectrum']['prediction_ph_cm2_s'],'o-',ms=3,
            color=colors[i%len(colors)],label=f'$R_{{in}}={inside:g},\\ R_{{out}}={outside:g}$ cm$^{{-1}}$')
        cells.append([f'{inside:g} / {outside:g}',
            f'{row["spectrum"]["descriptive_chi_squared"]:,.0f}',
            f'{100*row["power_budget"]["modeled_supply_fraction_of_spin_down"]:.2f}%'])
    ax.set(xscale='log',yscale='log',ylim=(1e-13,3e-5),xlabel='Photon energy (GeV)',
        ylabel='Integrated bin photon flux (cm$^{-2}$ s$^{-1}$)')
    ax.legend(frameon=False,fontsize=9,loc='lower left');ax.spines[['top','right']].set_visible(False)
    table.axis('off')
    t=table.table(cellText=cells,colLabels=['Inside / outside gain','Spectrum $\\chi^2$','Modeled supply / $L_{sd}$'],
        loc='upper center',colWidths=[.38,.28,.4],cellLoc='center',bbox=[0,.57,1,.35])
    t.auto_set_font_size(False);t.set_fontsize(10)
    table.text(.02,.47,'Each spectrum uses one fitted population normalization\nSeven detections, six descriptive degrees of freedom\nUpper limit checked separately\n\nFixed particle paths and observer geometry\nTwo global gain values, no peak-specific weights\nDistance 290 pc explicitly adopted\n\nThe trials still disagree with the spectrum\nA sufficient energy budget does not establish a fit',
        transform=table.transAxes,va='top',linespacing=1.5,fontsize=11)
    fig.suptitle('Vela development test — response to global particle acceleration',fontsize=17)
    fig.text(.5,.02,'Conditional on an unreleased numerical field · no independent observational validation · colors identify trials',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(args.output/'comparison.png',dpi=180)
    fig.savefig(args.output/'comparison.pdf');plt.close(fig)
    print(args.output/'comparison.png')


if __name__=='__main__':main()
