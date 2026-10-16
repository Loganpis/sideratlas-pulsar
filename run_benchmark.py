#!/usr/bin/env python3
"""Reproduce the field, pulse, data comparison, plots and verification report."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import sys
import time
from dataclasses import asdict,replace

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl-cache'))
import numpy as np
from pulsar.model import Parameters,TAU,pulse_profile,normalize_profile,peak_separation
from pulsar.validation import field_validation,timing_validation
from pulsar.observations import load_vela,compare_to_vela


def json_ready(obj):
    if isinstance(obj,np.ndarray): return obj.tolist()
    if isinstance(obj,np.generic): return obj.item()
    if isinstance(obj,dict): return {k:json_ready(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)): return [json_ready(v) for v in obj]
    return obj


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=ROOT/'outputs')
    args=ap.parse_args(); out=args.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    p=Parameters(); phases=np.arange(512)/512
    print('Checking fields and arrival times...',flush=True)
    field=field_validation(p); timing=timing_validation(p)
    profiles=[]; convergence=[]
    for res in ((12,32,64),(20,64,128),(24,96,192),(32,128,256)):
        t0=time.monotonic(); raw=pulse_profile(p,phases,res); norm=normalize_profile(raw)
        item={'resolution':res,'seconds':time.monotonic()-t0,'mean_unnormalized_flux':float(raw.mean())}
        if profiles:
            item['normalized_L1_change']=float(np.mean(abs(norm-profiles[-1]['norm'])))
            item['mean_flux_relative_change']=float(abs(raw.mean()/profiles[-1]['raw'].mean()-1))
        profiles.append({'raw':raw,'norm':norm}); convergence.append(item)
        np.savez_compressed(out/('profile_'+'_'.join(map(str,res))+'.npz'),
                            phases=phases,raw=raw,normalized=norm)
        print('Quadrature',res,'done',flush=True)
    norm=profiles[-1]['norm']; raw=profiles[-1]['raw']
    # Independent radial and phase refinements resolve different error sources.
    radial=normalize_profile(pulse_profile(p,phases,(48,128,256)))
    half_phases=(np.arange(512)+0.5)/512
    half_raw=pulse_profile(p,half_phases,(32,128,256))
    half=normalize_profile(half_raw)
    # A common normalization is required to interlace two half grids.
    fine_phases=np.arange(1024)/1024
    fine_raw=np.empty(1024); fine_raw[::2]=raw; fine_raw[1::2]=half_raw
    fine_norm=normalize_profile(fine_raw)
    interp_half=np.interp(half_phases,phases,norm,period=1.)
    convergence_extra={
        'radial_normalized_L1_change_32_to_48':float(np.mean(abs(radial-norm))),
        'phase_interpolation_L1_error_512_to_1024':float(np.mean(abs(half-interp_half))),
        'thresholds':{'angular_L1':1e-4,'angular_mean_relative':1e-4,'radial_L1':1e-4,'phase_L1':0.003},
    }
    converged=(convergence[-1]['normalized_L1_change']<1e-4 and
               convergence[-1]['mean_flux_relative_change']<1e-4 and
               convergence_extra['radial_normalized_L1_change_32_to_48']<1e-4 and
               convergence_extra['phase_interpolation_L1_error_512_to_1024']<0.003)
    separation,peaks=peak_separation(fine_phases,fine_norm)
    print('Comparing with independent Fermi data...',flush=True)
    data_path=ROOT/'data/J0835-4510_3PC_data.fits'
    data=load_vela(data_path)
    comparison=compare_to_vela(fine_phases,fine_norm,data)
    # Change one physical assumption at a time, at already converged angular resolution.
    sensitivity={}
    for name,params in [('tilt_70_deg',replace(p,obliquity_deg=70.)),
                        ('sheet_width_0.10',replace(p,sheet_width=0.10))]:
        curve=normalize_profile(pulse_profile(params,phases,(24,96,192)))
        sensitivity[name]=curve
    # Removing retardation makes the radial integrand much less smooth; it needs
    # its own radial convergence check instead of the main profile's settings.
    no_delay=normalize_profile(pulse_profile(p,phases,(80,96,192),delays=False))
    no_delay_refined=normalize_profile(pulse_profile(p,phases[::2],(128,96,192),delays=False))
    no_delay_error=float(np.mean(abs(normalize_profile(no_delay[::2])-no_delay_refined)))
    converged=converged and no_delay_error<1e-4
    convergence_extra['no_delay_diagnostic_radial_L1_change_80_to_128']=no_delay_error
    import matplotlib,astropy
    report={
        'model':'Analytic oblique split-monopole outer wind + prescribed optically thin sheet emissivity',
        'status':'Numerical benchmark only; not a validated reconstruction of Vela or a global PIC model',
        'parameters':asdict(p),'light_cylinder_radius_km':p.r_lc_km,
        'field_validation':field,'timing_validation':timing,'quadrature_convergence':convergence,
        'additional_convergence':convergence_extra,'quadrature_converged':bool(converged),
        'numerical_verification_passed':bool(field['passed'] and timing['passed'] and converged),
        'finite_beam_peak_separation_cycles':separation,'finite_beam_peak_phases':peaks,
        'data':{k:v for k,v in data.items() if not isinstance(v,np.ndarray)},
        'data_sha256':hashlib.sha256(data_path.read_bytes()).hexdigest(),
        'observation_comparison':{k:v for k,v in comparison.items() if not isinstance(v,np.ndarray)},
        'physical_model_validated_for_Vela':False,
        'software':{'python':platform.python_version(),'numpy':np.__version__,
                    'matplotlib':matplotlib.__version__,'astropy':astropy.__version__},
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in [ROOT/'run_benchmark.py',*sorted((ROOT/'src/pulsar').glob('*.py'))]},
        'sources':{
          'field':'https://arxiv.org/abs/astro-ph/9907051',
          'pulse_geometry':'https://arxiv.org/abs/1011.3431',
          'radiative_transfer':'https://academic.oup.com/mnras/article/434/3/2636/1045892',
          'viewing_angle':'https://arxiv.org/abs/0710.4168',
          'Fermi_data':'https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/',
        },
    }
    np.savetxt(out/'pulse-profile.csv',np.column_stack((fine_phases,fine_norm,fine_raw)),delimiter=',',
               header='arrival_phase_cycles,relative_flux_mean_one,arbitrary_flux',comments='')
    np.savetxt(out/'vela-comparison.csv',np.column_stack((data['phase_lo'],data['phase_hi'],data['counts'],data['error'],comparison['prediction'],comparison['residual_sigma'])),delimiter=',',
               header='phase_lo,phase_hi,weighted_counts,uncertainty,scaled_model,residual_sigma',comments='')
    (out/'validation.json').write_text(json.dumps(json_ready(report),indent=2,allow_nan=False)+'\n')
    (out/'parameters.json').write_text(json.dumps(asdict(p),indent=2)+'\n')
    from pulsar.figures import make_figures
    make_figures(out,p,phases,norm,fine_phases,fine_norm,no_delay,sensitivity,data,comparison,report)
    write_report(out,report)
    print(json.dumps({'numerical_verification_passed':report['numerical_verification_passed'],
                      'separation_cycles':separation,'reduced_chi_squared':comparison['reduced_chi_squared'],
                      'output':str(out)},indent=2),flush=True)
    if not report['numerical_verification_passed']:
        raise SystemExit('Verification failed: inspect outputs/validation.json')


def write_report(out,r):
    c=r['observation_comparison']; f=r['field_validation']; t=r['timing_validation']
    last=r['quadrature_convergence'][-1]
    text=f'''# NO.002 — Field and pulse benchmark

**Numerical verification: {"PASS" if r['numerical_verification_passed'] else "FAIL"}**  
**Physical reconstruction of Vela: NOT VALIDATED**

The calculation evaluates an analytic, oblique split-monopole wind and integrates a prescribed emitting current sheet over three-dimensional space at equal photon arrival time. It is a reduced outer-wind benchmark. It does not evolve a dipolar magnetosphere, solve particle acceleration, or reproduce NASA's kinetic simulation.

## Results

| Check | Result |
| --- | --- |
| Maxwell equations away from sheet | Pass; finite differences converge at order {min(f['observed_orders'].values()):.3f} or better |
| Force-free identities and field-compatible streaming | Pass; maximum relative identity error {max(f['identities'].values()):.2e} |
| Moving-sheet Maxwell boundary conditions | Pass; maximum residual {max(f['sheet_Maxwell_jumps'].values()):.2e} |
| Radial electromagnetic luminosity | {f['poynting_luminosities'][0]:.10f}; analytic value 2/3 in normalized units |
| Independent pencil-beam peak separation | Numerical {t['numerical_pencil_beam_separation']:.6f}; analytic {t['analytic_pencil_beam_separation']:.6f} cycles |
| Angular refinement, mean-normalized profile L1 difference | {last['normalized_L1_change']:.3g} |
| Radial refinement L1 difference | {r['additional_convergence']['radial_normalized_L1_change_32_to_48']:.3g} |
| 512-to-1024 phase interpolation L1 error | {r['additional_convergence']['phase_interpolation_L1_error_512_to_1024']:.3g} |
| Finite-beam pulse separation | {r['finite_beam_peak_separation_cycles']:.6f} cycles |
| Vela shape comparison | Reduced chi-squared {c['reduced_chi_squared']:.2f}, with {c['degrees_of_freedom']} degrees of freedom |

L1 differences are averages of absolute differences between profiles normalized to mean one. They describe discretization stability, not accuracy relative to a real pulsar. Maxwell differences are normalized locally by |B|/r. Sheet force balance is excluded.

## What the observation test means

The official Fermi >100 MeV profile contains two copies of one rotation for display. We verified they match exactly and used only the first. Its 800 independent bins were combined in groups of four; weighted counts were summed and supplied variances added. The result contains 200 bins. The original FITS file and its provenance are retained in `data/`.

Only the phase origin, nonnegative amplitude and nonnegative constant component were fitted. Obliquity, width, bulk Lorentz factor and viewing angle stayed fixed. The constant can include background and unmodeled unpulsed emission. The predicted curve was averaged over each observation bin. The fitted phase is {c['phase_shift_cycles']:.6f} cycles; it is not a prediction of the radio/gamma lag.

The residuals reject this simple prescribed-emission configuration as a quantitative model of Vela. The comparison is a shape diagnostic using catalog uncertainties, not an instrument-response or spectral fit. A more detailed particle/emission model is needed before claiming a Vela reconstruction. A physically sensible double pulse alone does not validate that claim.

## Chosen and inferred inputs

| Parameter | Value | Evidential status |
| --- | --- | --- |
| Period | 89 ms | Rounded Vela reference; no epoch-resolved timing prediction |
| Light-cylinder radius | {r['light_cylinder_radius_km']:.2f} km | Derived from this rounded period |
| Magnetic obliquity | 60 degrees | Chosen benchmark, not a Vela measurement |
| Observer angle | 63.6 degrees | Adopted from a model-dependent nebular-torus inference |
| Emission radius | 1.5–5 light-cylinder radii | Chosen outer-wind region |
| Bulk Lorentz factor | 10 | Chosen streaming prescription, not particle acceleration output |
| Sheet width | 0.05 in the dimensionless sheet function | Chosen emission width; not a fixed length in km |
| Radial emissivity | r^-3 | Prescribed |
| Comoving spectrum | j'_nu proportional to nu^-1 | Prescribed, without a spectral cutoff |
| Absolute field and luminosity | Unspecified | Only normalized quantities are calculated |

The field's discontinuous surface and the finite-width emission layer are different: the former is an analytic electromagnetic boundary, the latter a phenomenological radiating region. Dissipation and plasma supply within that layer are not solved. The selected bulk velocity is compatible with the limiting field outside the sheet, but is only prescribed within the emission layer. The radiation normalization is arbitrary and is not coupled to the electromagnetic energy budget. GR, stellar surface radiation, opacity, radio emission, polarization, and feedback of radiation on the flow are omitted.

## Deliverables

- `field-and-pulse.png` / `.pdf`: field geometry, pulse profile and timing comparison
- `vela-comparison.png` / `.pdf`: Fermi data, scaled model and residuals
- `convergence.png`: spatial refinement and physical sensitivity
- `pulse-profile.csv`: 1024 arrival-phase samples
- `vela-comparison.csv`: measured data, errors, model and residuals
- `validation.json`: complete machine-readable diagnostics
- `parameters.json`: inputs used

See [MODEL.md](../MODEL.md) for equations, source roles, and reproduction instructions. This benchmark establishes tested numerical components for the animation. Captions using it should identify the reduced analytic field and prescribed emission model.
'''
    (out/'RESULTS.md').write_text(text)


if __name__=='__main__':
    main()
