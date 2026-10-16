#!/usr/bin/env python3
"""Compare fixed-field quadrature runs using conservative phase-bin averages."""
from pathlib import Path
import argparse,json,sys,os
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
from pulsar.emission_quadrature import rebin_phase_averages
from pulsar.snapshot_fields import SpinDownScale
from pulsar.radiation import MEC2_ERG
from pulsar.energetics import population_power_budget
from pulsar.multiband import load_bands
from run_emission_quadrature import dump,sha


def rebinned(arrays,width,bins,shift=0):
    edges=np.arange(bins+1)/bins
    return rebin_phase_averages(arrays['phase_edges'],arrays[f'profiles_width_{width:g}'],edges[:-1],edges[1:],shift)


def comparison_invariants(protocol):
    """Only trajectory step and cap quadrature may differ in this study."""
    signature={k:v for k,v in protocol.items() if k not in
        ('snapshot','cap','transport','geometry_cache','geometry_cache_key','source_hashes')}
    signature['snapshot_sha256']=protocol['snapshot']['sha256']
    signature['cap']={k:v for k,v in protocol['cap'].items() if k not in ('azimuths','rings')}
    signature['transport']={k:v for k,v in protocol['transport'].items() if k!='dt'}
    return signature


def peaks(profile,bins,windows=((.05,.25),(.45,.75))):
    """Declared development windows, no smoothing or individual peak shifts."""
    phase=(np.arange(bins)+.5)/bins;answer=[]
    for p in profile[:4]:
        ids=[]
        for low,high in windows:
            available=np.flatnonzero((phase>=low)&(phase<high));ids.append(int(available[np.argmax(p[available])]))
        answer.append({'peak_bin_centres_cycles':phase[ids].tolist(),'half_bin_cycles':.5/bins,
                       'separation_cycles':float(phase[ids[1]]-phase[ids[0]]),
                       'peak2_over_peak1':float(p[ids[1]]/p[ids[0]]) if p[ids[0]]>0 else None,
                       'both_are_local_maxima':bool(all(p[i]>=p[(i-1)%bins] and p[i]>=p[(i+1)%bins] for i in ids))})
    return answer


def compare_runs(left,a,right,b,widths):
    """Compare unshifted finite-bin profiles; a fitted phase never enters."""
    changed=[k for k in ('azimuths','rings','dt','field_sha256') if left[k]!=right[k]]
    comparison={'reference':left['name'],'refined':right['name'],'changed_parameters':changed,'width_scenarios':{}}
    for width in widths:
        if str(width) not in right['width_scenarios']:continue
        answer={}
        for bins in (100,50):
            x=rebinned(a,width,bins);y=rebinned(b,width,bins)
            mx=x.mean(axis=1);my=y.mean(axis=1)
            xunit=np.divide(x,mx[:,None],out=np.zeros_like(x),where=mx[:,None]>0)
            yunit=np.divide(y,my[:,None],out=np.zeros_like(y),where=my[:,None]>0)
            answer[f'pulse_unit_mean_L1_{bins}_bins']=np.mean(abs(xunit[:4]-yunit[:4]),axis=1)
            # These windows are in the physical model phase before any fitted
            # alignment; they are the same for every member of this study.
            px=peaks(x,bins,((.15,.4),(.6,.85)));py=peaks(y,bins,((.15,.4),(.6,.85)))
            drift=np.asarray([p['peak_bin_centres_cycles'] for p in py])-np.asarray([p['peak_bin_centres_cycles'] for p in px])
            answer[f'unshifted_peak_bin_centre_change_{bins}_bins']=drift
            answer[f'unshifted_peak_separation_change_{bins}_bins']=[v['separation_cycles']-u['separation_cycles'] for u,v in zip(px,py)]
            answer[f'peak_bin_width_cycles_{bins}_bins']=1/bins
        answer['phase_mean_relative_change_all_bands']=(my-mx)/np.maximum(my,1e-300)
        sx=np.asarray(left['width_scenarios'][str(width)]['spectral_prediction_ph_cm2_s'])
        sy=np.asarray(right['width_scenarios'][str(width)]['spectral_prediction_ph_cm2_s'])
        answer['normalized_spectral_flux_relative_change']=(sy-sx)/np.maximum(sy,1e-300)
        pa=left['width_scenarios'][str(width)]['peaks_100_bins'];pb=right['width_scenarios'][str(width)]['peaks_100_bins']
        answer['peak_separation_change_cycles_100_bins']=[y['separation_cycles']-x['separation_cycles'] for x,y in zip(pa,pb)]
        comparison['width_scenarios'][str(width)]=answer
    return comparison


def observed_peaks(data):
    """Return observed 3PC peak bins in the same declared data-phase windows."""
    answer=[]
    for band in data:
        selected=[]
        for low,high in ((.05,.25),(.45,.75)):
            candidates=np.flatnonzero((band['phase']>=low)&(band['phase']<high))
            selected.append(int(candidates[np.argmax(band['counts'][candidates])]))
        answer.append({'band':band['band'],'peak_bin_centres_cycles':band['phase'][selected],
            'peak_bin_edges_cycles':np.stack((band['phase_lo'][selected],band['phase_hi'][selected]),axis=-1),
            'peak_weighted_counts':band['counts'][selected],
            'peak_count_errors':band['error'][selected],
            'separation_cycles':float(band['phase'][selected[1]]-band['phase'][selected[0]])})
    return answer


def write_markdown(path,assessment):
    """Readable, reproducible account of measured refinement differences."""
    rows={r['name']:r for r in assessment['runs']}
    lines=['# Fixed-field emission quadrature assessment','',
        'This study integrates calculated photons through finite observer and rotation-phase bins. '
        'It tests numerical sampling on one unreleased field; it does not validate a Vela reconstruction.','',
        'The primary observable is the observer-latitude average over 64°–66° with 100 phase bins. '
        'One-degree observer-width sensitivity and 50-bin comparisons are retained in the JSON report. '
        'No physical phase smoothing or separate alignment of numerical profiles is applied.','',
        '## Controlled refinement differences','',
        'Percentages are unit-mean profile L1 differences at the same physical phase. '
        'They measure differences between finite calculations, not a rigorous error bound against the continuum solution.','',
        '| Refinement | Changed parameter | 0.1–0.3 GeV | 0.3–1 GeV | 1–3 GeV | 3–100 GeV | Largest mean-flux change |',
        '|---|---|---:|---:|---:|---:|---:|']
    for c in assessment['controlled_comparisons']:
        parameter=c['changed_parameters'][0];a=rows[c['reference']];b=rows[c['refined']]
        if parameter in ('azimuths','rings') and b[parameter]!=2*a[parameter]:continue
        if parameter=='dt' and not np.isclose(b[parameter],a[parameter]/2):continue
        if '2.0' not in c['width_scenarios']:continue
        v=c['width_scenarios']['2.0'];dif=np.asarray(v['pulse_unit_mean_L1_100_bins'])*100
        flux=100*max(abs(np.asarray(v['phase_mean_relative_change_all_bands'][:4])))
        label=f"{a[parameter]:g} → {b[parameter]:g}"
        lines.append(f"| {label} | {parameter} | "+' | '.join(f'{x:.3f}%' for x in dif)+f' | {flux:.3f}% |')
    fine=max(assessment['runs'],key=lambda r:(r['azimuths'],r['rings']))
    v=fine['width_scenarios']['2.0'];budget=v['power_budget']
    lines+=['','## Finest azimuth calculation','',
        f"{fine['azimuths']} azimuths × {fine['rings']} annular rings × 2 poles = {fine['seed_count']:,} streams; path step {fine['dt']:g} light-cylinder radii.",
        '',f"The one-normalization spectrum has descriptive χ² = {v['spectrum_chi_squared']:,.3f} for six descriptive degrees of freedom. "
        'This remains a failed spectral reconstruction; the number is not a calibrated detector likelihood.','',
        '| Band | Measured peak-bin centres | Model centres after one fitted common phase |',
        '|---|---|---|']
    observed=assessment['development_observation']['peak_bins_100']
    for obs,model in zip(observed,v['peaks_100_bins']):
        op=', '.join(f'{x:.3f}' for x in obs['peak_bin_centres_cycles'])
        mp=', '.join(f'{x:.3f}' for x in model['peak_bin_centres_cycles'])
        lines.append(f"| {obs['band']} | {op} | {mp} |")
    lines+=['','Peak locations above are bin maxima in fixed, declared development windows. '
        'Each bin is 0.01 rotation wide. A weak or broad feature can switch between competing local maxima under refinement; '
        'the bin width is representation resolution, not a confidence interval or a bound on the underlying continuous peak. '
        'The JSON report separately preserves unshifted peak locations and their refinement changes.','',
        f"At the explicitly adopted distance of {assessment['distance_assumption_pc']:g} pc, the fitted population is "
        f"{budget['primaries_per_second']:.6g} primaries s⁻¹ across both poles. "
        f"Calculated curvature radiation uses {100*budget['radiation_fraction_of_spin_down']:.3f}% of measured spin-down power; "
        f"prescribed acceleration work plus injected kinetic energy uses {100*budget['modeled_supply_fraction_of_spin_down']:.3f}%. "
        'This necessary energy check does not determine the pair supply or complete the pulsar energy budget.','',
        'The magnetic field, enlarged numerical star, finite injection surface, angular-fraction annulus, flux loading, '
        'prescribed acceleration, and curvature-only emission remain explicit model assumptions. '
        'All four pulse bands and the spectral data are development observations. These results do not establish an independently validated Vela reconstruction for a final Reel.','',
        '## Provenance','',f"The finest run is `{fine['path']}`. Its result SHA-256 is `{fine['results_sha256']}` and profile SHA-256 is `{fine['profile_sha256']}`. "
        'The adjacent `assessment.json` records every run, field and source hash, observer-width sensitivity, spectral prediction, and power budget. '
        'Regenerate this report with `assess_emission_quadrature.py` and the same ordered run directories.','']
    path.write_text('\n'.join(lines))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('runs',type=Path,nargs='+')
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--distance-pc',type=float,default=290.)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    observation_path=ROOT/'data/J0835-4510_3PC_data.fits'
    observed=observed_peaks(load_bands(observation_path))
    runs=[];invariants=None
    for path in args.runs:
        r=json.loads((path/'results.json').read_text());p=r['protocol']
        signature=comparison_invariants(p)
        if invariants is None:invariants=signature
        elif signature!=invariants:
            raise ValueError('Fixed-field comparison may change only cap azimuths, rings and trajectory dt')
        with np.load(path/'profiles.npz') as d:arrays={k:d[k] for k in d.files}
        scale=SpinDownScale(**p['physical_scale']);work=0.
        for trace in sorted(Path(p['geometry_cache']).glob('trace-*.npz')):
            with np.load(trace) as g:
                radius=np.linalg.norm(g['position'],axis=-1)
                rate=np.where(radius<p['acceleration']['transition_radius_rlc'],p['acceleration']['rate_inside_cm_inverse'],p['acceleration']['rate_outside_cm_inverse'])
                work+=float(np.sum(rate*MEC2_ERG*g['dt']*scale.rlc_cm*g['seed_weight'][None]))
        row={'name':path.name,'path':str(path.resolve()),'seed_count':r['seed_count'],'azimuths':p['cap']['azimuths'],
             'rings':p['cap']['rings'],'dt':p['transport']['dt'],'field_sha256':p['snapshot']['sha256'],
             'profile_sha256':sha(path/'profiles.npz'),'results_sha256':sha(path/'results.json'),
             'source_hashes':p['source_hashes'],'width_scenarios':{}}
        for width in p['observer_widths_deg']:
            summary=r['width_scenarios'][str(width)];shift=summary['pulse_shape']['phase_shift_cycles']
            budget=population_power_budget(summary['spectrum']['normalization'],args.distance_pc,scale.luminosity_erg_s,
                r['energy_erg_per_primary'],work_erg_per_primary=work,injected_kinetic_erg_per_primary=(p['initial_gamma']-1)*MEC2_ERG)
            row['width_scenarios'][str(width)]={'selected_seeds':summary['selected_seeds'],
                'spectrum_chi_squared':summary['spectrum']['descriptive_chi_squared'],
                'spectral_prediction_ph_cm2_s':summary['spectrum']['prediction_ph_cm2_s'],
                'phase_averaged_yield_per_sr':summary['phase_averaged_yield_per_sr'],
                'fitted_common_phase_cycles':shift,
                'peaks_100_bins':peaks(rebinned(arrays,width,100,shift),100),
                'peaks_50_bins':peaks(rebinned(arrays,width,50,shift),50),
                'unshifted_peaks_100_bins':peaks(rebinned(arrays,width,100),100,((.15,.4),(.6,.85))),
                'unshifted_peaks_50_bins':peaks(rebinned(arrays,width,50),50,((.15,.4),(.6,.85))),
                'power_budget':budget}
            fitted=row['width_scenarios'][str(width)]['peaks_100_bins']
            row['width_scenarios'][str(width)]['fitted_minus_observed_peak_bin_centres_cycles']=np.asarray(
                [x['peak_bin_centres_cycles'] for x in fitted])-np.asarray([x['peak_bin_centres_cycles'] for x in observed])
        runs.append((row,arrays,r))
    comparisons=[]
    for (left,a,ra),(right,b,rb) in zip(runs[:-1],runs[1:]):
        comparisons.append(compare_runs(left,a,right,b,ra['protocol']['observer_widths_deg']))
    controlled=[]
    for i,(left,a,ra) in enumerate(runs):
        for right,b,rb in runs[i+1:]:
            comparison=compare_runs(left,a,right,b,ra['protocol']['observer_widths_deg'])
            if len(comparison['changed_parameters'])==1:controlled.append(comparison)
    width_changes=[]
    for row,arrays,report in runs:
        widths=report['protocol']['observer_widths_deg']
        for broad,narrow in zip(widths[:-1],widths[1:]):
            x=rebinned(arrays,broad,100);y=rebinned(arrays,narrow,100)
            width_changes.append({'run':row['name'],'from_width_deg':broad,'to_width_deg':narrow,
                'pulse_unit_mean_L1_100_bins':np.mean(abs(x[:4]/x[:4].mean(axis=1)[:,None]-y[:4]/y[:4].mean(axis=1)[:,None]),axis=1),
                'phase_mean_relative_change_all_bands':(y.mean(axis=1)-x.mean(axis=1))/np.maximum(y.mean(axis=1),1e-300)})
    assessment={'scope':'Fixed-field numerical sampling study; development data, no Vela validation',
        'fixed_comparison_settings':invariants,
        'assessment_source_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),
            *sorted((ROOT/'src/pulsar').glob('*.py'))]},
        'pulse_L1_definition':'mean absolute difference of separately unit-mean phase profiles; common physical phase, no independent alignment',
        'peak_definition':'Peak-bin centres after one shared fitted phase, in declared development windows [0.05,0.25) and [0.45,0.75); no smoothing or separate peak adjustment',
        'unshifted_peak_definition':'Peak-bin centres at physical model phase in fixed development windows [0.15,0.4) and [0.6,0.85); bin width reports representation resolution, not a confidence interval or a continuous-peak accuracy bound',
        'distance_assumption_pc':args.distance_pc,'distance_scope':'Explicit Barnard-like distance scenario, not an independently exact measured distance',
        'development_observation':{'path':str(observation_path.resolve()),'sha256':sha(observation_path),
            'peak_bins_100':observed,'scope':'All four pulse bands are development data; peak-bin maxima are descriptive and do not include a fitted peak-position uncertainty'},
        'runs':[x[0] for x in runs],'adjacent_comparisons':comparisons,'controlled_comparisons':controlled,'observer_width_changes':width_changes,
        'physical_phase_blur':0.,'field_release_allowed':False,'emission_release_allowed':False,'vela_validation':False}
    dump(args.output/'assessment.json',assessment)
    write_markdown(args.output/'RESULTS.md',assessment)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':11})
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    labels=['0.1–0.3 GeV','0.3–1 GeV','1–3 GeV','3–100 GeV'];phase=(np.arange(100)+.5)/100
    for row,arrays,report in runs:
        p=rebinned(arrays,report['protocol']['observer_widths_deg'][0],100)[:4];p/=p.mean(axis=1)[:,None]
        for i,ax in enumerate(axes.flat):
            ax.stairs(p[i],np.arange(101)/100,label=f"{row['azimuths']} az × {row['rings']} rings, dt={row['dt']:g}",linewidth=1)
            ax.set(xlabel='Physical model phase',ylabel='Unit-mean finite-bin profile',title=labels[i],xlim=(0,1))
    axes[0,0].legend(fontsize=7,ncol=2);fig.suptitle('Fixed-field quadrature sensitivity · no physical phase smoothing')
    fig.savefig(args.output/'profiles.png',dpi=160);fig.savefig(args.output/'profiles.pdf');plt.close(fig)
    print(json.dumps({'output':str(args.output),'runs':len(runs)}))


if __name__=='__main__':main()
