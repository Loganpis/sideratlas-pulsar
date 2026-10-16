#!/usr/bin/env python3
"""Reproducible independent field and absorbing-boundary experiments."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
import time

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.mpl-cache'))

from pulsar.field_benchmarks import (
    ExactDipoleBoundary, analytic_maxwell_residual, analytic_spherical_flux,
    conductor_gauge_seam, finite_conductor_reference, gauge_sponge_check, point_dipole_luminosity,
    wave_packet_absorber,
)


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def plots(out, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), constrained_layout=True)
    for name, symbol in (('e_relative_L2', 'E'), ('b_relative_L2', 'B')):
        axes[0].plot([x['n'] for x in result['source_evolution']],
                     [100*x['final'][name] for x in result['source_evolution']], 'o-', label=symbol)
    axes[0].set(xlabel='Grid N (N³ cells)', ylabel='Relative exterior field error (%)',
                title='Exact dipole boundary / vacuum evolution')
    axes[0].legend()
    for row in result['packets']:
        if row['n'] == 128 and row['filter_strength'] == 0:
            label = ('A/E sponge' if row['kind'] == 'potential' else 'Matched E/B')+f", strength {row['strength']:g}"
            axes[1].semilogy([x['time'] for x in row['history']],
                             [max(x['reflected_energy_fraction'], 1e-17) for x in row['history']], label=label)
            axes[2].semilogy([x['time'] for x in row['history']],
                             [max(x['total_energy_fraction'], 1e-17) for x in row['history']], label=label)
    axes[1].set(xlabel='Time (Ω⁻¹)', ylabel='Leftgoing energy / initial energy',
                title='1-D packet reflection (zero clipped at 10⁻¹⁷)', ylim=(1e-17, 1))
    axes[2].set(xlabel='Time (Ω⁻¹)', ylabel='Total energy / initial energy',
                title='Absorption and periodic leakage', ylim=(1e-9, 2))
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.15)
    fig.savefig(out/'field-benchmarks.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'outputs/field-benchmarks')
    parser.add_argument('--duration', type=float, default=2*np.pi)
    parser.add_argument('--resolutions', type=int, nargs='+', default=[48, 64])
    args = parser.parse_args()
    if not np.isfinite(args.duration) or args.duration <= 0:
        parser.error('duration must be finite and positive')
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    source_paths = ['src/pulsar/plasma.py', 'src/pulsar/field_benchmarks.py', 'run_field_benchmarks.py']
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths}
    protocol = dict(source='Uniformly rotating point magnetic dipole; exact A/E inner and outer Dirichlet data',
                    source_radius=.65, resolutions=args.resolutions, duration=args.duration,
                    source_evolution_shell=[1., 2.4], packet_resolutions=[48, 64, 128, 256],
                    packet_geometry='1-D normal incidence, baseline 8-wide periodic domain and cubic edge damping',
                    source_sha256=hashes,
                    generated_utc=datetime.now(timezone.utc).isoformat(),
                    runtime=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
                    gates=dict(source_field_relative_L2=.10, source_flux_relative_error=.10,
                               packet_reflected_energy_fraction=.001, packet_wrap_energy_fraction=.001),
                    gate_scope='Source gates apply at the specified final time, measured in 1<r<2.4; packet gates apply to maximum interior characteristic energies over their stated windows',
                    scope='Numerical vacuum/boundary verification only; no Vela or plasma-emission validation')
    dump(out/'protocol.json', protocol)
    result = dict(protocol=protocol,
                  analytic_maxwell=[analytic_maxwell_residual(n) for n in (32, 48, 64, 96)],
                  analytic_flux=[dict(radius=r, value=analytic_spherical_flux(r),
                                      expected=point_dipole_luminosity()) for r in (.5, 1., 2., 4., 10.)],
                  source_gauge_seam=[conductor_gauge_seam(n) for n in (48, 64, 96)],
                  gauge_sponge=[gauge_sponge_check(n) for n in (48, 64, 128)],
                  finite_conductor=finite_conductor_reference(),
                  source_evolution=[], packets=[])
    baseline_path = ROOT/'outputs/plasma-upgrade/field-validation.json'
    if baseline_path.exists():
        baseline = json.loads(baseline_path.read_text())
        mean = baseline['poynting']['Vacuum 64³']['last_rotation_mean'][1]
        expected = result['finite_conductor']['finite_luminosity']
        result['legacy_vacuum_comparison'] = dict(
            source_file=str(baseline_path.relative_to(ROOT)),
            source_file_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
            legacy_snapshot_generation_source_sha256=baseline.get('snapshot_provenance', {}).get('vacuum-n64.npz', {}).get('generation_source_sha256'),
            radius=1.5, last_rotation_mean_flux=mean,
            finite_conductor_expected=expected,
            relative_deviation=abs(mean/expected-1),
            claim='Legacy result compared with the physically appropriate ideal finite-conductor luminosity; sharp-surface convergence is not demonstrated')
    for n in args.resolutions:
        start = time.monotonic()
        model = ExactDipoleBoundary(n=n)
        history = [model.benchmark_diagnostics()]
        interval = max(1, round(.5/(.3*model.h)))
        while model.time < args.duration:
            model.step(min(.3*model.h, args.duration-model.time))
            if model.step_number % interval == 0 or model.time >= args.duration:
                diag = model.benchmark_diagnostics()
                history.append(diag)
                print(f"Exact source {n}³: t={model.time:.3f}, E error={diag['e_relative_L2']:.4f}, B error={diag['b_relative_L2']:.4f}", flush=True)
        diag = history[-1]
        accepted = max(diag['e_relative_L2'], diag['b_relative_L2']) < .1 and max(diag['flux_relative_error']) < .1
        result['source_evolution'].append(dict(n=n, elapsed_seconds=time.monotonic()-start,
                                              source_sha256=hashes, accepted=accepted,
                                              final=diag, history=history))
        dump(out/f'exact-source-n{n}.json', result['source_evolution'][-1])
        del model
    for n in protocol['packet_resolutions']:
        for kind, strength, filt in [('potential', 4., 0.), ('matched_fields', 4., 0.),
                                     ('matched_fields', 12., 0.), ('potential', 4., .08),
                                     ('matched_fields', 12., .08)]:
            value = wave_packet_absorber(n, kind, strength, filt)
            value['source_sha256'] = hashes
            value['accepted'] = value['peak_reflected_energy_fraction'] < .001 and value['peak_periodic_wrap_energy_fraction'] < .001
            result['packets'].append(value)
            print(f"Packet {n}: {kind}, σ={strength:g}, filter={filt:g}: reflection={value['peak_reflected_energy_fraction']:.3g}, wrap={value['peak_periodic_wrap_energy_fraction']:.3g}", flush=True)
    finest = max(result['source_evolution'], key=lambda row: row['n'])['final']
    result['source_propagation_passed'] = max(finest['e_relative_L2'], finest['b_relative_L2']) < .1 and max(finest['flux_relative_error']) < .1
    result['baseline_absorber_passed'] = all(x['peak_reflected_energy_fraction'] < .001 and x['peak_periodic_wrap_energy_fraction'] < .001 for x in result['packets'] if x['kind'] == 'potential' and x['filter_strength'] == 0)
    result['matched_1d_absorber_passed'] = all(x['peak_reflected_energy_fraction'] < .001 and x['peak_periodic_wrap_energy_fraction'] < .001 for x in result['packets'] if x['kind'] == 'matched_fields' and x['strength'] == 12 and x['filter_strength'] == 0)
    result['plasma_solver_revalidated'] = False
    result['vela_validation'] = False
    dump(out/'results.json', result)
    plots(out, result)
    print(json.dumps({k: v for k, v in result.items() if k.endswith('passed')}, indent=2))


if __name__ == '__main__':
    main()
