# Calculated physical pair-radiation component

This folder adds a **Vela-informed prescribed physical radiation model**. It
does not relabel the reduced local PIC particles as physical Vela particles.
The output is observer-specific optical photon power and CIE tristimulus
power, rather than species-colored tracks. It remains an assumed emission
model, not a fitted or validated reconstruction of Vela's optical image.

## Physical inputs and what is solved

The explicit reference inputs are P=0.08937108859 s, R=13 km, polar B=3×10¹² G,
inclination 75° and observer 65°. Existing pairs are assigned gamma 300 and
pitch 0.1 rad in a prescribed activation layer between 0.08 and 0.18 light-cylinder
radii. Their supply is a chosen static dipolar Goldreich–Julian cap flux,
multiplied by 6000 pairs/primary (12000 total leptons/primary), across cap
coordinates 0.8–0.9, both poles and all azimuths. Equal electron/positron
intensities are summed incoherently with physical number weights.

Harding et al. 2018 motivate cap pairs, the multiplicity, a low-altitude optical
synchrotron component and this viewing geometry. Their separate cascade pair
spectrum and resonant radio-absorption calculation are **not reproduced** by
our monoenergetic, prescribed activation layer. The initial spectrum,
activation distribution and pitch are assumptions, not measured Vela data.

After activation, guiding centres follow the physical static dipole in a local
magnetic frame. The momentum equations include adiabatic focusing and
gyro-averaged classical radiation reaction. With dimensionless parallel and
perpendicular momenta p_parallel and p_perp,

```
gamma² = 1 + p_parallel² + p_perp²
a = e⁴ B² / (6 pi epsilon0 m_e³ c³)
d p_parallel/dt = -(p_perp²/(2 p_parallel)) dlnB/dt
                  - a p_parallel p_perp²/gamma
d p_perp/dt = +(p_perp/2) dlnB/dt
              - a p_perp(1+p_perp²)/gamma
d E_radiated/dt = m_e c² a p_perp²
```

The focusing terms conserve energy; the drag terms reproduce the Liénard
power. A separate energy ledger checks their combined numerical solution.
Cooling continues to p_perp=0.001, with the remaining perpendicular-energy
bound recorded. It does not stop at the synchrotron continuum boundary.
Gyroperiod/cooling time, gyroperiod/field-variation time, gyro-radius/field
gradient length and quantum parameter chi are reported for every column.

Rigid rotation of these prescribed cells adds a local Lorentz transformation
of the photons. It is an emission-pattern assumption, not a global Maxwell
solution or a relativistic continuity solution for the pair supply. Native
PIC gamma values, charges and fields are never multiplied into this model.
No field feedback, pair cascade, resonant radio heating, inverse Compton,
curvature loss, absorption, extinction, GR or surface atmosphere is solved in
this component. Those omissions must remain in the bundled description.

## Why exact harmonics are used

Even a relativistic population can radiate optical photons in its first few
gyro harmonics when B is high. A smooth synchrotron continuum can therefore
give the wrong directional optical spectrum. The production calculation uses
the classical **Schott harmonic spectrum**, not a smooth curve or a Gaussian
beam. The two Bessel polarization terms and their delta-line resonance are
evaluated with physical gamma, B, pitch and observer angle.

With A=1-beta_parallel cos(theta) and nu_B=eB/(2pi gamma m_e),

```
nu_n = n nu_B/A
x_n = n beta_perp sin(theta)/A
dP_n/dOmega = [e² nu_B²/(2 epsilon0 c)] [n²/A³]
              * {[(cos(theta)-beta_parallel)/sin(theta)]² J_n(x_n)²
                 + beta_perp² J'_n(x_n)²}
```

This is emitted power per **source time**, integrated over the delta line.
It follows the resonance integration of Marcowith & Malzac 2003 equations 1, 2, 6
converted to SI. No linewidth is invented. A separate high-harmonic continuum
kernel is retained for a benchmark; it does not produce the movie data.

For a magnetic cell moving with local corotation velocity, photons are
aberrated into its magnetic frame. D=nu_lab/nu_local. Source-time spectral
power transforms as P_nu=D²/Gamma_bulk P'_nu'(nu/D); consequently a line's
**frequency-integrated** power transforms as D³/Gamma_bulk P'_n. Arrival-time
compression is not included in this power and must be applied through
retarded photon deposition, not an additional arbitrary beaming factor.

Very weak off-axis harmonic sums are skipped only using a generating-function
upper bound on Bessel functions. The discarded per-lepton band power is below
10⁻⁴⁵W/sr; its number-weighted energy bound is retained in the report. This is
a numerical error control, not a widened beam or a display brightness cutoff.

## Camera interface

Run from the repository root:

```
python research/optical-finish/particle-emission/run_component.py \
  --output outputs/optical-finish/particle-emission/baseline-v2 \
  --phases 128 --azimuths 64 --cooling 32 --radii 3 --wavelength-bins 64
```

A short smoke run can use
`--phases 8 --azimuths 16 --cooling 24 --radii 2`.

The higher spatial sampling run uses:

```
python research/optical-finish/particle-emission/run_component.py \
  --output outputs/optical-finish/particle-emission/production-spatial \
  --phases 256 --azimuths 128 --cooling 64 --radii 5 --annuli 2 \
  --wavelength-bins 64 --compact-camera
```

The compact schema retains each event's band power and exact-line XYZ,
aggregate per-phase wavelength-bin spectra, reference-cell geometry and event
indices. Geometry, time, number weights and pattern velocity reconstruct
deterministically with `load_camera_events.py`. Calculations use float64;
compact XYZ storage uses float32 and records the actual number-weighted XYZ
energy error bound. This reduces storage; it does not change the photon
kernel, broaden beams or smooth pulses. The compact archive omits individual
per-event spectra, so use its phase-summed spectra for the optical audit.

The NPZ contains:

| Array | Meaning |
|---|---|
| positions_m | Physical Cartesian source positions at source time |
| source_time_s, source_interval_s | Midpoint and full lab source-time quadrature interval |
| number_weights | Physical lepton count, not macroparticle charge; multiply once |
| power_bins_w_sr | Per-lepton line power summed into each wavelength bin, W/sr/bin |
| wavelength_edges_m | Actual spectral bin edges; these are not smooth spectral samples |
| xyz_lm_sr | Per-lepton CIE1931 XYZ luminous power/sr; 683 lm/W already applied once |
| pattern_beta | Rigidly rotating cell velocity/c; not individual gyro velocity |
| phase_index, phase_angles_rad | Prescribed source-pattern rotation phase |
| cell_index | Index into fixed reference-cell geometry |
| reference_cell_positions_m | Cell centres at reference rotation phase 0 |
| cell_footprint_vertices_reference_m | Eight physical injection-coordinate support corners |
| observer_direction | Fixed source-to-observer unit vector; direction is already used in spectra |
| period_s | One physical pattern rotation, seconds |

XYZ is evaluated from **each individual harmonic's actual wavelength** against
the official CIE table; it is not reconstructed from bin-averaged spectra.
The official CIE data SHA256 is retained in the JSON report. A camera must
multiply either XYZ or radiant power by `number_weights*source_interval_s`
once and deposit its energy at `tau=t-n.r/c`, modulo period. Use moving-pattern
interval endpoints for conservative finite-exposure accumulation. The included
simple midpoint arrival histogram is an energy-accounting diagnostic only.

The finite footprints describe the physical injection-coordinate quadrature.
They are not an optical point-spread function. Treating each support cell as
uniform for rasterization is an explicitly finite spatial approximation;
directional emission varies inside it. Refinement is required to assess it.
Rotating a precomputed directional spectrum for a different camera is invalid;
rerun its observer kernel for the new viewing direction.

## Verification and remaining uncertainty

`test_physical_synchrotron.py` checks exact line integration against independent
Liénard and direct time-domain Liénard–Wiechert power, the nonrelativistic
angular dipole limit, low-harmonic and high-harmonic limits, source-time boost
factors, exact-line XYZ normalization, spectral-bin energy conservation,
zero-pitch emission and the cooling ledger. These checks validate numerical
ingredients, not Vela's adopted injection or global field.

The run report records actual phase sampling, number of active cells, pulse
coverage, total optical energy and omitted harmonic-power bounds. A separate
refinement report must quantify phase/azimuth/cooling/radial sampling. The
Vela optical slope and absolute flux must be compared as predictions, with
no normalization tuning silently promoted into validation.

The initial refinement study finds only 0.4–3.1% changes in total band energy,
but 4.7–8.2% changes in source-phase profile L1 norm. Neither scalar stability
nor exact CIE integration establishes resolved image or spectral convergence.
Narrow spectral bins vary more than band-integrated power for the assumed
monoenergetic distribution. At 290 pc this chosen population's integrated
optical flux is roughly twenty times an approximate intrinsic Vela reference,
and its band-averaged spectral slope is too blue. The final per-run audit retains the
actual numbers; display exposure must not be mistaken for physical fitting.
No Vela optical validation is claimed.

The completed higher spatial grid contains 5,271,172 nonzero source events.
It predicts 9.617×10⁻¹⁷ W/m² over 400–700 nm at the chosen 290 pc distance,
about 19.55 times the rough intrinsic reference. Its binned spectral index is
alpha≈−0.544 versus the reported Vela value −0.04±0.04; the narrow-bin model
spectrum itself remains sensitive to quadrature. Compared with the preceding
256-phase / 128-azimuth / 64-cooling run, increasing the radial/annulus grid
changes total energy by only −0.109%, but changes the source-phase profile
L1 norm by 9.68%. The image/pulse is therefore a finite-resolution model,
not a converged reconstruction. Compact color storage contributes at most
2.17×10⁻⁸ relative error to number-weighted XYZ energy, far below this physical
and spatial uncertainty. Full results and hashes are in the production audit.

## Primary references

- [Harding et al. 2018](https://arxiv.org/abs/1811.11157), DOI 10.3847/2041-8213/aaf3b2
- [Harding & Kalapotharakos 2015](https://arxiv.org/abs/1508.06251), DOI 10.1088/0004-637X/811/1/63
- [Marcowith & Malzac 2003](https://arxiv.org/abs/astro-ph/0309107), DOI 10.1051/0004-6361:20031060
- [Lai et al. 2025](https://arxiv.org/abs/2509.04572), DOI 10.1093/mnras/staf1295
- [Mignani et al. 2007](https://arxiv.org/abs/0707.2036), DOI 10.1051/0004-6361:20077774
- [Official CIE1931 2° dataset](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer), DOI 10.25039/CIE.DS.xvudnb9b; CC BY-SA4.0
