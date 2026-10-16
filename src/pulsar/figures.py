"""Static scientific figures; colors encode geometry, not visible radiation."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .model import TAU


def make_figures(out,p,phases,norm,fine_phases,fine_norm,no_delay,sensitivity,data,comparison,report):
    plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix',
                         'font.size':12,'axes.titlesize':16,'axes.labelsize':12,
                         'savefig.dpi':180,'figure.constrained_layout.use':True})
    cyan='#31b9d4'; gold='#edab4c'; red='#cc504f'; muted='#a5b2bd'
    with plt.rc_context({'figure.facecolor':'#080d13','axes.facecolor':'#080d13',
                         'text.color':'#edf2f5','axes.labelcolor':'#edf2f5',
                         'axes.edgecolor':'#6b7d8c','xtick.color':'#c7d1da','ytick.color':'#c7d1da'}):
        fig=plt.figure(figsize=(13.5,7.3)); gs=fig.add_gridspec(2,2,width_ratios=[1.12,1])
        ax=fig.add_subplot(gs[:,0],projection='3d')
        r=np.linspace(p.r_min,p.r_max,90)[:,None]
        eta=np.linspace(0,TAU,150)[None,:]
        psi=-r
        # Each radius cuts the retarded sheet in a great circle.
        x=r*(np.cos(eta)*np.cos(p.chi)*np.cos(psi)-np.sin(eta)*np.sin(psi))
        y=r*(np.cos(eta)*np.cos(p.chi)*np.sin(psi)+np.sin(eta)*np.cos(psi))
        z=np.broadcast_to(-r*np.cos(eta)*np.sin(p.chi),x.shape)
        ax.plot_surface(x,y,z,color=gold,alpha=.55,linewidth=0,rcount=70,ccount=100,shade=True)
        rr=np.linspace(p.r_min,p.r_max,140)
        for theta in np.deg2rad([25,55,90,125,155]):
            for phi0 in [0,TAU/3,2*TAU/3]:
                phi=phi0-(rr-p.r_min)
                ax.plot(rr*np.sin(theta)*np.cos(phi),rr*np.sin(theta)*np.sin(phi),
                        rr*np.cos(theta),color=cyan,lw=.65,alpha=.8)
        ax.plot([0,0],[0,0],[-5.4,5.4],color=muted,lw=.7,ls='--')
        ax.text(0,0,5.5,'spin axis',color=muted,fontsize=10)
        for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
            axis.pane.fill=False; axis.pane.set_edgecolor('#24303b')
            axis._axinfo['grid']['color']='#23313b'
        ax.set(xlabel=r'$x/R_{\rm LC}$',ylabel=r'$y/R_{\rm LC}$',zlabel=r'$z/R_{\rm LC}$',
               xlim=(-5,5),ylim=(-5,5),zlim=(-5,5))
        ax.set_box_aspect((1,1,1)); ax.view_init(elev=24,azim=-57)
        ax.set_title('Calculated outer-wind geometry',pad=10)
        ax.text2D(.01,.97,'Amber: current sheet   •   Cyan: magnetic field lines\nInner region excluded; explanatory colors; equal simulation time',
                  transform=ax.transAxes,fontsize=10,color=muted,va='top')
        a=fig.add_subplot(gs[0,1]); a.plot(fine_phases,fine_norm,color=gold,lw=2)
        a.set(title='Radiation reaching a fixed observer',ylabel='Relative flux (mean = 1)',xlim=(0,1))
        a.text(.5,.78,r'$\chi=60^\circ,\ \zeta=63.6^\circ,\ \Gamma=10$'+'\nPrescribed sheet emissivity',ha='center',va='top',transform=a.transAxes,fontsize=10,color=muted)
        a.grid(alpha=.15)
        b=fig.add_subplot(gs[1,1]); b.plot(phases,norm,color=gold,lw=2,label='With photon travel times')
        b.plot(phases,no_delay,color=cyan,lw=1.5,ls='--',label='Delays omitted: diagnostic only')
        b.set(xlabel='Arrival phase / rotation',ylabel='Relative flux',xlim=(0,1),title='Why arrival time matters')
        b.legend(frameon=False,fontsize=10); b.grid(alpha=.15)
        fig.suptitle('@sideratlas  /  NO.002     Field-and-pulse benchmark',fontsize=21)
        fig.savefig(out/'field-and-pulse.png'); fig.savefig(out/'field-and-pulse.pdf'); plt.close(fig)

    fig,(ax,res)=plt.subplots(2,1,figsize=(11,7.4),sharex=True,gridspec_kw={'height_ratios':[3,1]})
    scale=np.max(data['counts'])
    ax.errorbar(data['phase'],data['counts']/scale,yerr=data['error']/scale,fmt='.',ms=3,
                color='#283c53',elinewidth=.6,label='Fermi-LAT Vela, >100 MeV')
    ax.plot(data['phase'],comparison['prediction']/scale,color=red,lw=2,label='Fixed physical model; phase, scale and constant fitted')
    ax.set(ylabel='Weighted counts / observed maximum',ylim=(0,1.13),
           title='Observation test: the simple emission prescription is insufficient')
    ax.legend(frameon=False,fontsize=11,loc='upper right')
    ax.text(.98,.72,rf'$\chi^2/\mathrm{{dof}}={comparison["reduced_chi_squared"]:.1f}$'+'\n200 independent bins; 3 fitted parameters',
            ha='right',va='top',transform=ax.transAxes,fontsize=12)
    res.axhline(0,color='black',lw=.7); res.axhspan(-3,3,color='#eeeeee')
    res.plot(data['phase'],comparison['residual_sigma'],color=red,lw=1)
    res.set(xlabel='Observed rotational phase',ylabel=r'Residual / $\sigma$',xlim=(0,1))
    for a in (ax,res): a.grid(alpha=.15)
    fig.suptitle('NO.002  •  Numerical verification does not establish a Vela reconstruction',fontsize=17)
    fig.savefig(out/'vela-comparison.png'); fig.savefig(out/'vela-comparison.pdf'); plt.close(fig)

    fig,(a,b)=plt.subplots(1,2,figsize=(12,4.8))
    conv=report['quadrature_convergence'][1:]
    a.loglog([int(np.prod(v['resolution'])) for v in conv],[v['normalized_L1_change'] for v in conv],
             'o-',color='#187a96',label='Change from previous spatial grid')
    a.axhline(1e-4,color=red,ls='--',lw=1,label='Acceptance threshold')
    a.set(xlabel='Spatial quadrature points',ylabel='Mean absolute profile change',title='Spatial convergence')
    a.legend(frameon=False,fontsize=10); a.grid(alpha=.2,which='both')
    b.plot(phases,norm,color='#283c53',lw=2,label=r'Baseline: $\chi=60^\circ$, width 0.05')
    b.plot(phases,sensitivity['tilt_70_deg'],color='#187a96',label=r'Change tilt to $70^\circ$')
    b.plot(phases,sensitivity['sheet_width_0.10'],color=red,ls='--',label='Change width to 0.10')
    b.set(xlabel='Arrival phase / rotation',ylabel='Relative flux (mean = 1)',xlim=(0,1),title='Dependence on physical assumptions')
    b.legend(frameon=False,fontsize=10); b.grid(alpha=.2)
    fig.savefig(out/'convergence.png'); plt.close(fig)
