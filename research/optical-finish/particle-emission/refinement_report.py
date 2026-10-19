#!/usr/bin/env python3
"""Compare actual radiation quadrature runs without asserting convergence."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(base, fine):
    q0 = np.asarray(base['source_band_power_per_phase_w_sr'])
    q1 = np.asarray(fine['source_band_power_per_phase_w_sr'])
    if len(q1) % len(q0):
        raise ValueError('Fine source phase grid must be an integer refinement')
    averaged = q1.reshape(len(q0), -1).mean(axis=1)
    params0, params1 = base['model'], fine['model']
    changed = {k: [v, params1[k]] for k, v in params0.items() if params1[k] != v}
    e0, e1 = base['total_source_frame_band_energy_j_sr'], fine['total_source_frame_band_energy_j_sr']
    return dict(changed_parameters=changed,
        baseline_band_energy_j_sr=e0, refined_band_energy_j_sr=e1,
        relative_energy_change_vs_baseline=float((e1-e0)/e0),
        relative_energy_difference_normalized_to_refined=float(abs(e1-e0)/e1),
        normalized_source_phase_profile_L1=float(np.abs(averaged-q0).sum()/averaged.sum()),
        largest_source_phase_power_difference_over_refined_peak=float(np.max(np.abs(averaged-q0))/averaged.max()),
        phase_profile_comparison='Fine source-time phase bins averaged to baseline bins; no artificial phase shifting, fitting, peak alignment or pulse smoothing',
        one_percent_energy_and_profile_gate_passed=bool(abs(e1-e0)/e1 < .01 and np.abs(averaged-q0).sum()/averaged.sum() < .01),
        spatial_image_convergence_established=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--baseline', default='baseline-v2')
    p.add_argument('--variants', nargs='+', default=['phase-refined', 'azimuth-refined', 'cooling-refined', 'radial-refined', 'annulus-refined', 'combined-refined', 'production-spatial'])
    args = p.parse_args()
    base = json.loads((args.root/args.baseline/'physical-pair-report.json').read_text())
    result = dict(schema='sideratlas.physical-pair-refinement.v1', baseline=args.baseline,
        status='Finite quadrature sensitivity study, not a Vela-model or image validation',
        no_parameters_fit_to_optical_data=True, comparisons={})
    for name in args.variants:
        path = args.root/name/'physical-pair-report.json'
        if path.is_file():
            result['comparisons'][name] = compare(base, json.loads(path.read_text()))
        else:
            result['comparisons'][name] = dict(status='not_completed_or_missing; no convergence claim')
    result['limitations'] = [
        'One-at-a-time changes are not a full multidimensional convergence proof',
        'Source phase profile comparison does not validate final finite-exposure image rasterization',
        'Band-integrated scalar convergence can coexist with unresolved spatial/color structure',
        'The explicit population and field assumptions dominate scientific uncertainty even when quadrature errors are small']
    out = args.root/'refinement-report.json'
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: {a: b for a, b in v.items() if a in ('relative_energy_change_vs_baseline', 'normalized_source_phase_profile_L1', 'one_percent_energy_and_profile_gate_passed', 'status')} for k, v in result['comparisons'].items()}, indent=2))


if __name__ == '__main__':
    main()
