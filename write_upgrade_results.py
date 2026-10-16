"""Summarize saved numerical evidence; never infer success from a rendered image."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent


def read(relative):
    return json.loads((ROOT/relative).read_text())


def main():
    field=read('outputs/plasma-upgrade/field-validation.json')
    rows='\n'.join(f"| {name.replace('_to_','³ → ')}³ | {100*v['e_relative_L2']:.2f}% | {100*v['b_relative_L2']:.2f}% |"
                   for name,v in field['field_comparisons'].items())
    path=ROOT/'outputs/plasma-upgrade/refinement-validation.json'
    if path.exists():
        finer=json.loads(path.read_text())
        v=finer['96_to_128_restarted']
        refinement=f"""The 96³ snapshot was interpolated onto 128³ and evolved for one further rotation. Compared with its initial 96³ state at the same rotation phase, the electric-field difference is **{100*v['e_relative_L2']:.2f}%** and the magnetic-field difference is **{100*v['b_relative_L2']:.2f}%**. This comparison combines spatial refinement with extra relaxation time; it cannot by itself establish convergence. [Machine-readable check](refinement-validation.json)"""
        matched_path=ROOT/'outputs/plasma-upgrade/matched-time-validation.json'
        if matched_path.exists():
            matched=json.loads(matched_path.read_text());m=matched['same_time_96_to_128'];p=matched['coarse_3_to_4_periods']
            refinement+=f"""\n\nA separate 96³ control evolved for the same extra rotation. Comparing both grids at $t=4P$ gives **{100*m['e_relative_L2']:.2f}% electric-field change** and **{100*m['b_relative_L2']:.2f}% magnetic-field change**. The 10% same-time field target {'is met in this check' if matched['same_time_field_change_below_target'] else 'is not met'}. The 96³ field alone changes by {100*p['e_relative_L2']:.2f}% (electric) and {100*p['b_relative_L2']:.2f}% (magnetic) between three and four periods, measuring residual coarse-grid evolution. Even an improving warm-start result does not replace an independent fine startup, timestep/boundary studies or resolved pair physics. [Matched-time evidence](matched-time-validation.json)"""
    else:
        refinement='The 128³ warm-start check has not completed. No result is inferred from a partial history.'
    flux=field['poynting_relative_change_64_to_96']
    field_text=f"""# Evolving plasma fields — measured results

**Status: prototype verified on analytic cases, global magnetosphere not converged. Not released for Vela photon transport.**

The 48³, 64³ and 96³ calculations evolve three nominal rotation periods with the same configuration. The table compares electric and magnetic fields at the same final rotation phase over $0.8<r/R_{{\\rm LC}}<2.4$, interpolating the coarser grid onto the finer grid. Relative L2 differences use the finer field norm.

| Grid comparison | Electric field | Magnetic field |
| --- | ---: | ---: |
{rows}

The declared target is 10% for each field. The 64³ → 96³ comparison fails. The final-rotation mean outward Poynting flux changes by **{100*flux[0]:.2f}%, {100*flux[1]:.2f}% and {100*flux[2]:.2f}%** at radii 1, 1.5 and 2 $R_{{\\rm LC}}$, also exceeding the 10% target. These thresholds are development gates, not astrophysical error bars.

## Finer restart

{refinement}

## Independent checks and limitations

The periodic vacuum-wave benchmark has relative field error $2.67\\times10^{{-6}}$ against its discrete analytic solution, with relative energy change $-4.26\\times10^{{-7}}$. The magnetic divergence constraint is preserved to floating-point precision in the global runs. These checks validate particular numerical operations, not the astrophysical field geometry.

The rotating-star vacuum comparator has mean luminosity {field['poynting']['Vacuum 64³']['last_rotation_mean'][1]:.4f} at $1.5R_{{\\rm LC}}$, versus the point-dipole reference {field['vacuum_point_dipole_reference']:.4f}: a **{100*field['vacuum_64_relative_deviation_at_r1p5']:.2f}%** discrepancy. Finite enlarged stellar size, forcing, filtering, resolution and boundary effects all differ from the ideal point-source reference; this is an additional warning, not a passed vacuum calibration.

The star spans only 2.8 grid cells in radius at 64³ and 4.2 at 96³. Conductivity is prescribed. Pair production, particle inertia, kinetic reconnection and general-relativistic surface physics are absent. The sponge has not passed a reflection test, and its signed energy exchange can include injection. An exactly closing energy ledger is bookkeeping; numerical energy changes remain reported separately.

The field inclination is 65°, whereas the separate, frozen legacy emission diagnostic retains 45°. They are different experiments, not successive stages of one Vela reconstruction.

![Calculated fields and resolution dependence](evolved-fields.png)

See [full assumptions and equations](../../REALISM-UPGRADE.md), [independent solver review](../../research/PLASMA-SOLVER-REVIEW.md), [protocol](protocol.json), and [numerical evidence](field-validation.json). Native snapshots and time histories are saved alongside these files.
"""
    (ROOT/'outputs/plasma-upgrade/RESULTS.md').write_text(field_text)
    emission=read('outputs/emission-upgrade/validation.json')
    rows='\n'.join(f"| {v['band']} | {emission['scores']['flat'][j]['reduced_chi_squared']:.2f} | {v['reduced_chi_squared']:.2f} | {emission['uniform_flux_loading_scores'][j]['reduced_chi_squared']:.2f} |"
                 for j,v in enumerate(emission['scores']['uniform_cosine']))
    loading=emission['loading_L1_by_band']; sky=emission['sky_bin_1_to_half_degree_L1_by_band']
    gamma=emission['injection_gamma_L1_by_band']['1000000.0']
    view_rows='\n'.join(f"| {angle}° | "+' | '.join(f"{v['reduced_chi_squared']:.2f}" for v in values)+' |'
                        for angle,values in emission['viewing_angle_sensitivity'].items())
    emission_text=f"""# Loading, detector response and radiation — measured results

**Status: implemented physics changes, but Vela's measured pulse shape remains badly reproduced. No validated Vela reconstruction.**

The diagnostic traces 7,168 seeds in the frozen legacy rotating-dipole geometry. It uses the previously selected 45° inclination, 63.6° viewing angle, acceleration rate 0.06 cm⁻¹ and phase shift 0.9770689775403691 cycles. **The new numerical plasma field is not used for these photons.** Geometry and phase were not refitted.

The upgraded profile uses charge-based loading $|\\rho_{{\\rm corotation}}|\\,|v_r|\\,dA$ with the full curl(B) correction at the injection surface, zero additional phase blur and a 0.5° numerical viewing bin. The detector fold uses official P8R3_SOURCE_V3 FRONT/BACK effective area and energy dispersion, averaged over an explicitly assumed incidence distribution.

## Conditional shape comparison

All four energy bands are **development data**. Two nonnegative nuisance parameters—scale and constant—are fitted separately per band. The values below are descriptive chi-squared divided by 98 for 100 bins, with the catalog's supplied weighted-count errors. They are not calibrated likelihood probabilities, absolute flux fits or fresh independent validation.

| Energy band | Charge loading, flat response | Charge loading, LAT response | Uniform-flux loading, LAT response |
| --- | ---: | ---: | ---: |
{rows}

The response correction alone does not fix the broad/incorrect peak structure. Charge-based loading improves some bands and worsens others relative to uniform-flux loading. In particular, the calculated high-energy second peak remains much too weak. These failures are retained in the report rather than hidden by phase-dependent brightness adjustments.

![Vela profiles and model scenarios](response-comparison.png)

## Sensitivity, not confidence intervals

The profile comparisons below use mean absolute differences between unit-mean phase profiles. Percentages are relative to that unit mean, not parameter uncertainties.

- Charge-based versus uniform-flux loading changes profiles by **{100*min(loading):.1f}–{100*max(loading):.1f}%**, depending on energy band. Individual normalized seed weights are multiplied by {emission['seed_weight_relative_to_uniform_flux_range'][0]:.3f}–{emission['seed_weight_relative_to_uniform_flux_range'][1]:.3f}
- Narrowing the numerical observer bin from 1° to 0.5° changes profiles by **{100*min(sky):.2f}–{100*max(sky):.2f}%**. No artificial Gaussian phase blur is applied
- Changing injection Lorentz factor from $10^3$ to $10^6$ changes profiles by **{100*min(gamma):.1f}–{100*max(gamma):.1f}%**. The injection energy is therefore not an innocuous choice
- Uniform-cosine, near-axis and off-axis detector scenarios are saved separately. None is measured Vela exposure, and their spread is not a confidence band

The viewing-angle scenarios retain the same geometry and phase, fitting only the per-band nuisance terms:

| Viewing angle | 0.1–0.3 GeV | 0.3–1 GeV | 1–3 GeV | 3–100 GeV |
| --- | ---: | ---: | ---: | ---: |
{view_rows}

## Additional radiation channels

Classical synchrotron and isotropic Klein–Nishina inverse-Compton kernels have independent analytic-limit, spectral-integral and convergence tests. Their illustrative comparison specifies $B=50$ kG, pitch $=0.001$ rad, $T=0.7$ MK, stellar radius 12 km, and distance/curvature radius 4246.5 km. **These are declared inputs, not inferred Vela populations or a fitted Vela spectrum.** No arbitrary component mixture was added to the pulse calculation. Stellar seed-photon directions were isotropized for this kernel demonstration, which is not adequate for a predictive anisotropic inverse-Compton model.

![Single-particle radiation powers for declared illustrative conditions](radiation-channels.png)

The next physical dependency is a converged plasma field and justified particle injection, pair supply, pitch evolution and seed-photon geometry. A complete observation comparison also needs source-specific exposure, spatial acceptance and catalog weighting. See [radiation details](../../research/RADIATION-UPGRADE.md), [detector provenance and limitations](../../research/DETECTOR-UPGRADE.md), [protocol](protocol.json), and [complete numerical results](validation.json).
"""
    (ROOT/'outputs/emission-upgrade/RESULTS.md').write_text(emission_text)
    print('Wrote evidence-based field and emission reports')


if __name__=='__main__': main()
