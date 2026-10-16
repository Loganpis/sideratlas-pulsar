#!/usr/bin/env python3
"""Summarize completed artifacts without promoting diagnostics to validation."""
from pathlib import Path
import json,hashlib,datetime
ROOT=Path(__file__).resolve().parent


def read(relative):return json.loads((ROOT/relative).read_text())
def pct(x):return f'{100*x:.2f}%'


def main():
    original=read('outputs/coffee-assessment/results.json')
    revised=read('outputs/coffee-surface-assessment/results.json')
    verification=read('outputs/current-verification.json')
    sampling=read('outputs/current-assessment/sampling-comparison.json')
    # Refuse a stale test certificate after any source/test change.
    for relative,digest in verification['source_sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'Tested file changed: {relative}')
    names=['outer-coffee-n65-surface-anchored','outer-coffee-n65-refined',
           'outer-coffee-n97-surface-anchored']
    runs=[]
    for name in names:
        folder=ROOT/'outputs'/name
        path=folder/'results.json'
        if not path.exists():
            failure=folder/'failure.json'
            if not failure.exists():raise ValueError(f'Unfinished diagnostic: {name}')
            runs.append({'name':name,'completed':False,'failure':json.loads(failure.read_text())['message']})
            continue
        result=json.loads(path.read_text())
        runs.append({'name':name,'completed':True,'seeds':result['segment_shape'][1],
            'selected_seeds':result['selected_seeds'],
            'energy_outside_light_cylinder':result['fraction_energy_outside_light_cylinder'],
            'energy_inward':result['fraction_energy_inward_photons'],
            'max_energy_balance_error_erg':result['max_energy_balance_error_erg'],
            'spectrum_chi2':result['spectrum']['descriptive_chi_squared'],
            'spectrum_dof':result['spectrum']['descriptive_dof'],
            'pulse_chi2_per_descriptive_dof':[b['reduced_chi_squared'] for b in result['profiles']['bands']],
            'guiding_centre_diagnostics':result['guiding_centre_diagnostics']})
    status={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'completed_scope':'Boundary repair experiment, numerical verification, and outer-emission development diagnostics',
        'emission_release_allowed':False,'vela_reconstruction_validated':False,
        'final_reel_produced':False,'tests_passed':verification['tests'],
        'original_grid_difference':original['resolution_comparison']['0.75<r<2'],
        'revised_grid_difference':revised['resolution_comparison']['0.75<r<2'],
        'sampling':sampling,'runs':runs,
        'remaining_requirements':['Converged global field with controlled divergence, stellar-boundary, domain and numerical-star-size sensitivity',
            'Converged particle/observer quadrature and resolved trajectory curvature',
            'Justified particle supply and acceleration; resolved FIDO interface if that alternative is used',
            'Consistent observational response and successful spectral/pulse comparison; independent validation data']}
    out=ROOT/'outputs/current-assessment';out.mkdir(exist_ok=True)
    (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    lines=['# @sideratlas · NO.002 / Current scientific status','',
        '**The implementation and verification pass is complete. The Vela reconstruction and final Reel are not ready.** The revised boundary repairs reversed polar fields and allows the full field → particle → radiation calculation to run. The measured field, sampling and observational checks still do not justify a Vela claim. All outputs retain `emission_release_allowed=false`.','',
        f'**{verification["tests"]} tests pass.** The [test manifest](outputs/current-verification.json) records exact source hashes and dependency versions; [full log](outputs/current-tests.txt). These tests verify stated equations and numerical identities, not the global astrophysical model.','',
        '## What was fixed','',
        '- Native CPU initialization now assigns the field arrays correctly; superseded runs remain identified',
        '- A separate, fully recorded boundary experiment anchors B and E at the same numerical stellar surface, repairing the reversed polar field',
        '- Cap tracing recognizes closure on return to the launch sphere, while rejecting nulls, unresolved paths and sampled disconnected cap topology',
        '- Particle motion extends outside the light cylinder with laboratory trajectory curvature, signed work, energy losses and propagation delays',
        '- The Barnard scenario uses its published gain per path length without an extra velocity projection',
        '- The sharp FIDO screening boundary now rejects crossing derivative/integration stencils instead of generating artificial curvature radiation',
        '- Genuine 3PC spectral bins use one normalization across all energies; no manual second-peak brightening is inserted','',
        '## Measured field outcome','',
        'Both grid sizes start independently and run for three rotations. The boundary experiment changes the numerical source model; it is not a claim of exact reproduction of upstream COFFEE. The star is enlarged to 0.5 light-cylinder radii.','',
        '| Check | Original boundary | Revised boundary |','|---|---:|---:|']
    for field in ('E','B'):
        key=f'relative_l2_{field}'
        lines.append(f'| 65³ → 97³ {field} relative L2 difference, 0.75 < r < 2 | {pct(original["resolution_comparison"]["0.75<r<2"][key])} | {pct(revised["resolution_comparison"]["0.75<r<2"][key])} |')
    for n in ('65','97'):
        lines.append(f'| {n}³ B change between rotations 2 and 3 | {pct(original["periodicity"][n]["relative_l2_B"])} | {pct(revised["periodicity"][n]["relative_l2_B"])} |')
        lines.append(f'| {n}³ RMS divergence / RMS B, in R_LC⁻¹ | {original["quality"][n]["relative_rms_divB_RLC"]:.4f} | {revised["quality"][n]["relative_rms_divB_RLC"]:.4f} |')
    lines+=['','A nearly repeating solution can still have substantial spatial error. In the revised model, the physical divergence metric does not decrease under refinement (0.4733 → 0.4744 per R_LC), and the current-sheet proxy still differs by 43–45%. This needs a divergence-preserving source/boundary treatment; more sampling alone will not repair it. Two grids cannot establish asymptotic convergence. Source-size, outer-domain and oblique-boundary checks remain incomplete.','',
        '![Measured field assessment](outputs/current-assessment/field-assessment.png)','',
        'Full evidence: [original grid comparison](outputs/coffee-assessment/results.json), [revised comparison](outputs/coffee-surface-assessment/results.json), [boundary experiment](research/COFFEE-BOUNDARY.md), [cap/source audit](outputs/field-benchmarks/CAP-AUDIT.md).','',
        '## Completed particle and observation comparisons','',
        'The development scenarios use inclination 75°, observer 65°, finite-radius injection at 0.7 R_LC, prescribed magnetic-flux loading and Barnard scenario-2 rates. None predicts pair creation or an absolute population.','',
        '| Calculation | Injected / observer-selected seeds | Energy beyond cylindrical light cylinder | Spectral χ² (6 descriptive dof) |','|---|---:|---:|---:|']
    for run in runs:
        label=run['name'].replace('outer-coffee-','')
        if run['completed']:
            lines.append(f'| [{label}](outputs/{run["name"]}/results.json) | {run["seeds"]} / {run["selected_seeds"]} | {pct(run["energy_outside_light_cylinder"])} | {run["spectrum_chi2"]:.1f} |')
        else:lines.append(f'| {label} | Stopped: {run["failure"]} | — | — |')
    lines+=['','The spectral score compares seven detections with one fitted normalization; the eighth bin is an upper limit. It is descriptive, not a calibrated joint likelihood. Pulse fits use a shared phase and separate scale/background per band with a flat detector response. All inspected data are development data.','',
        'On the same 65³ field, changing from 128 to 384 seeds and halving the orbit step changes the four pulse profiles by **'+', '.join(pct(x) for x in sampling['relative_L1_profiles'])+'** in relative L1 at 100 phase bins, without re-aligning them. This combined refinement does not establish pulse convergence.','',
        'The independent [emission audit](outputs/current-assessment/emission-audit.md) reproduces the physical scale, energy ledger, phase projection and spectral normalization. The coarse diagnostic has too few high-energy photons and substantial unresolved curvature. No fit has been used to compensate for those field errors.','',
        '[View the higher-resolution data comparison](outputs/outer-coffee-n97-surface-anchored/comparison.png) and [refined coarse-grid comparison](outputs/outer-coffee-n65-refined/comparison.png). Photons are not rendered as visible-light colors here; these are numerical and gamma-ray diagnostics.','',
        '## What still prevents a final Vela Reel','',
        '1. Establish a reliable global field: control divergence and numerical stellar/sponge boundaries, then demonstrate resolution, domain and stellar-size sensitivity',
        '2. Establish converged observer sampling and resolved trajectory curvature before tuning emission physics',
        '3. Justify particle populations and acceleration, then compare spectra and pulses with consistent observation selections/response and an independent validation test','',
        'The optional FIDO path cannot yet cross its sharp interface. A physically stated transition and a sensitivity study are required; choosing an arbitrary smoothing width merely to finish would add an unsupported assumption. The adopted moment of inertia, force-free torque normalization, enlarged star, fixed geometry, injection annulus, initial energy and acceleration remain model assumptions. Synchrotron/Compton kernels alone do not solve a pair cascade.','',
        'The source spectral FITS has invalid embedded checksum cards. Independent official downloads had identical bytes; the original data and warning are preserved in [spectral provenance](data/spectra/provenance.json).','',
        '## Reproduce and continue','',
        'Use Python 3.12 with NumPy, SciPy ≥ 1.12, Matplotlib, Astropy and h5py. Exact tested versions are in the manifest. The retained local interpreter was `/opt/anaconda3/bin/python`.','',
        '```sh','PYTHONPATH=src /opt/anaconda3/bin/python -m unittest discover -s tests -v',
        '/opt/anaconda3/bin/python run_outer_diagnostic.py \\\n  third_party/coffee-runs/pulsar-a75-n97-3rot-surface-anchored/snapshot-001020.npz \\\n  --inclination 75 --seed-radius .7 --seed-inner-radius .68 \\\n  --azimuths 64 --rings 3 --dt .004 --output outputs/new-diagnostic',
        '```','',
        'The separate native builds, source patches, licenses and exact field-run commands are recorded in [COFFEE](research/COFFEE.md) and [the boundary experiment](research/COFFEE-BOUNDARY.md). The earlier black-hole release is separate.','']
    (ROOT/'CURRENT-STATUS.md').write_text('\n'.join(lines))
    print(ROOT/'CURRENT-STATUS.md')


if __name__=='__main__':main()
