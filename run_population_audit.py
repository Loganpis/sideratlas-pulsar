#!/usr/bin/env python3
"""Check absolute population and power implied by a completed diagnostic.

Distance=290 pc is an explicit model adoption following Barnard et al. (2022),
not a new fitted/parallax estimate. A sensitivity bracket may be supplied.
This script never promotes an unreleased field or a failed spectral fit.
"""
from pathlib import Path
import argparse, hashlib, json, sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
from pulsar.energetics import population_power_budget
from pulsar.snapshot_fields import SpinDownScale
from pulsar.radiation import E_STATC, MEC2_ERG


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run', type=Path)
    ap.add_argument('--distance-pc', type=float, default=290.)
    ap.add_argument('--distance-bracket-pc', type=float, nargs=2, default=[270.,310.])
    args = ap.parse_args()
    path = args.run/'results.json'
    report = json.loads(path.read_text())
    scale = SpinDownScale(**report['protocol']['physical_scale'])
    if (args.run/'emission.npz').exists():
        with np.load(args.run/'emission.npz') as d:
            weights = d['seed_weight']
            if not np.isclose(weights.sum(), 1, rtol=1e-12, atol=1e-12):
                raise ValueError('Seed population must be globally normalized')
            radiation = float(np.sum(d['radiated_erg']*weights[None]))
            work = float(np.sum(E_STATC*d['e_work']*d['dt']*scale.rlc_cm*weights[None]))
            initial = float(np.sum((d['gamma'][0]-1)*MEC2_ERG*weights))
            final = float(np.sum((d['gamma'][-1]-1)*MEC2_ERG*weights))
        scenarios = {'point_sampled':report['spectrum']}
        extra = dict(work_erg_per_primary=work, injected_kinetic_erg_per_primary=initial,
                     escaping_kinetic_erg_per_primary=final)
    else:
        radiation = report['energy_erg_per_primary']
        scenarios = {key:row['spectrum'] for key,row in report['width_scenarios'].items()}
        extra = {}
        for key in ('work_erg_per_primary','injected_kinetic_erg_per_primary',
                    'escaping_kinetic_erg_per_primary'):
            if key in report:extra[key] = report[key]
    results = {}
    for key, spectrum in scenarios.items():
        results[key] = population_power_budget(spectrum['normalization'],args.distance_pc,
                                               scale.luminosity_erg_s,radiation,**extra)
        results[key]['spectrum_descriptive_chi_squared'] = spectrum['descriptive_chi_squared']
        results[key]['spectrum_descriptive_dof'] = spectrum['descriptive_dof']
        results[key]['distance_sensitivity'] = [population_power_budget(
            spectrum['normalization'],dist,scale.luminosity_erg_s,radiation,**extra)
            for dist in args.distance_bracket_pc]
    document = dict(source_run=str(path.resolve()),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        model_distance_reference='https://arxiv.org/html/2111.03405v1',
        distance_choice='290 pc adopted from Barnard model; bracket is a sensitivity test, not a confidence interval',
        normalization_scope='Total primary supply across both hemispheres for the declared seed loading; not pair multiplicity',
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in [Path(__file__),ROOT/'src/pulsar/energetics.py']},
        scenarios=results,emission_release_allowed=False,vela_validation=False)
    output = args.run/'population-audit.json'
    output.write_text(json.dumps(document,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'output':str(output),'scenarios':{
        key:{k:v for k,v in value.items() if k in ('primaries_per_second',
            'radiation_fraction_of_spin_down','modeled_supply_fraction_of_spin_down',
            'energy_ledger_residual_erg_per_primary')}
        for key,value in results.items()}},indent=2))


if __name__ == '__main__':main()
