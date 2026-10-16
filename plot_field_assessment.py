#!/usr/bin/env python3
"""Plot measured field differences; these are diagnostics, not visible light."""
from pathlib import Path
import argparse,json,os
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from run_coffee_assessment import load,sample


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,default=ROOT/'outputs/coffee-assessment/results.json')
    p.add_argument('--corrected',type=Path,required=True)
    p.add_argument('--output',type=Path,default=ROOT/'outputs/current-assessment')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    reports=[json.loads(path.read_text()) for path in (a.original,a.corrected)]
    plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':11})
    fig=plt.figure(figsize=(12,13),layout='constrained')
    grid=fig.add_gridspec(3,2,height_ratios=(1,1,.85))
    axes=[]
    line=np.linspace(-2,2,220);x,z=np.meshgrid(line,line)
    points=np.stack((x,np.zeros_like(x),z),axis=-1)
    for j,report in enumerate(reports):
        for k,key in enumerate(('coarse_path','fine_path')):
            d=load(Path(report[key]));ax=fig.add_subplot(grid[j,k]);axes.append(ax)
            b=sample(d,'b',points);mag=np.linalg.norm(b,axis=-1)
            image=ax.pcolormesh(x,z,np.log10(np.maximum(mag,1e-5)),
                cmap='magma',vmin=-2,vmax=2,shading='auto',rasterized=True)
            sub=slice(5,None,14);bx=b[sub,sub,0];bz=b[sub,sub,2]
            norm=np.hypot(bx,bz);outside=np.hypot(x[sub,sub],z[sub,sub])>.58
            ax.quiver(x[sub,sub][outside],z[sub,sub][outside],
                (bx/np.maximum(norm,1e-30))[outside],(bz/np.maximum(norm,1e-30))[outside],
                color='white',alpha=.65,scale=26,width=.0023)
            ax.add_patch(Circle((0,0),d['metadata']['star_radius'],facecolor='#b7bac2',edgecolor='white',zorder=4))
            for lc in (-1,1):ax.axvline(lc,color='cyan',lw=.8,ls='--',alpha=.65)
            n=len(d['axes'][0]);label='Original boundary' if j==0 else 'Revised boundary'
            ax.set(title=f'{label} · {n}³ cells',xlabel=r'$x/R_{\rm LC}$',ylabel=r'$z/R_{\rm LC}$',aspect='equal',xlim=(-2,2),ylim=(-2,2))
    fig.colorbar(image,ax=axes,shrink=.75,pad=.025,label=r'$\log_{10}[|B|/(\mu/R_{\rm LC}^3)]$')
    ax=fig.add_subplot(grid[2,0]);loc=np.arange(2);width=.3
    for j,(report,color,label) in enumerate(zip(reports,('#9d3c29','#207b8d'),('Original boundary','Revised boundary'))):
        row=report['resolution_comparison']['0.75<r<2']
        values=[100*row[f'relative_l2_{field}'] for field in ('E','B')]
        ax.bar(loc+(j-.5)*width,values,width,label=label,color=color)
        for xx,y in zip(loc+(j-.5)*width,values):ax.text(xx,y+1,f'{y:.1f}%',ha='center',fontsize=10)
    ax.set(xticks=loc,xticklabels=('Electric field','Magnetic field'),ylabel='Relative L2 grid difference (%)',title=r'Resolution sensitivity, $0.75<r/R_{\rm LC}<2$')
    ax.set_ylim(0,max(ax.get_ylim()[1]*1.12,10));ax.legend(frameon=False,fontsize=9)
    ax=fig.add_subplot(grid[2,1])
    for j,report in enumerate(reports):
        for k,(n,rows) in enumerate(report['surfaces'].items()):
            label=('Original' if j==0 else 'Revised')+f' {n}³'
            ax.plot([r['radius'] for r in rows],
                [r['luminosity_mu2_omega4_over_c3'] for r in rows],
                marker='o',ls='--' if j==0 else '-',color=('#aa5348','#d88d68','#207b8d','#162c51')[2*j+k],label=label)
    ax.axhline(reports[0]['reference_spindown'],color='gray',ls=':',label='Point-dipole reference')
    ax.set(xlabel=r'Sphere radius ($R_{\rm LC}$)',ylabel=r'$L/(\mu^2\Omega^4/c^3)$',title='Outgoing electromagnetic power')
    ax.legend(frameon=True,facecolor='white',framealpha=.95,fontsize=8)
    fig.suptitle('@sideratlas · NO.002 / Pulsar field verification',fontsize=21)
    fig.supxlabel('Three rotations · enlarged numerical star · arrows project B into the y=0 plane\nColors show field strength, not visible emission · no validated Vela reconstruction',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(a.output/f'field-assessment.{ext}',dpi=180)
    plt.close(fig)
    print(a.output/'field-assessment.png')


if __name__=='__main__':main()
