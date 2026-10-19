# @sideratlas · NO.002 · THE PULSAR · v1

This release is an original **Vela-informed scientific visualization** of a reduced outer-magnetosphere emission model. Its scope was selected explicitly after the numerical and observation comparisons. The image is a projected map of calculated emission, with illustrative colors and slowed rotation. It is not a photograph, a unique reconstruction of Vela, or a quantitative fit to Vela's observed pulses and spectrum.

The Reel uses the full vertical canvas, with only `@sideratlas`, `THE PULSAR` and `NO.002` over the animation. Scientific explanation and credits accompany the video in `instagram-description.text` and this document.

## Physical model

The field is the 97³, three-rotation `curl-radiative / no-GLM` candidate documented in [`research/COFFEE-CURL-REPAIR.md`](../research/COFFEE-CURL-REPAIR.md). It evolves a finite force-free electromagnetic field around a rotating conductor, with magnetic inclination **75°**. The model extends to a cube of half-width **3 light-cylinder radii**; the numerical star radius is **0.5 light-cylinder radii**. The light cylinder is the cylindrical distance at which rigid corotation would reach the speed of light.

The initial magnetic field is a fourth-order discrete curl of a finite vector potential. The rotating conductor imposes the electric relation `E = −(Ω × r) × B`; it does not reset evolved magnetic components to an analytic dipole. The outer electric boundary damps an incoming radiation characteristic without directly damping magnetic components. Magnetic-gradient GLM correction is disabled in this selected candidate. These choices improve the native magnetic constraint, but do not establish global physical convergence. The finite conductor has measured moment decay and phase lag, the transition spans two grid cells, and the outer boundary is an approximation with nonzero oblique-wave reflection.

The period **0.08937108859019084 s** and its catalog derivative set the dimensional scale, using an assumed moment of inertia of `10^45 g cm²` and an ideal force-free oblique-point-dipole torque normalization. This gives a light-cylinder radius of approximately **4,264.2 km**. An enlarged finite numerical star need not share the exact torque of the ideal point-dipole reference.

The selected run injects **1,536 modeled streams** at **0.7 R_LC**, from **128 azimuths × 6 rings × 2 hemispheres** near a sampled finite-domain cap rim. The annulus fractions are **0.90–0.96** on that injection surface. This differs from injection at a physical neutron-star surface. The outward velocity branch is selected explicitly, and seed weights follow normalized magnetic flux; a self-consistent electron/positron supply is not predicted.

The v1 run explicitly opts into a **finite-domain boundary-open** criterion: a magnetic trace either crosses the cylindrical light cylinder or escapes outward through a sphere of radius **2.5 R_LC**. This does not prove that a line opens at infinity. Of the selected annular magnetic traces, **1,376** cross the light cylinder and **160** reach the outgoing sphere first. Those 160 retain their outgoing classification when the magnetic-trace step is reduced from **0.008 to 0.002 R_LC**; their radial direction cosines at escape span **0.6777–0.9651**. The same criterion is applied to the sampled rim and selected seeds, with no silently discarded spokes. The code's default stricter light-cylinder criterion remains available. These magnetic traces define injection; they are distinct from the subsequent particle trajectories.

Particle directions use the ultrarelativistic guiding-center prescription of Kalapotharakos, Harding & Kazanas (2014), equation 12. Laboratory trajectory curvature is calculated from the evolving position and direction; it is not replaced with magnetic-field-line curvature. This approximation does not resolve gyration, solve the complete Lorentz dynamics, or feed particle energy and pair currents back into the electromagnetic field.

The selected Barnard et al. acceleration scenario prescribes `dγ/ds` gains of **0.04 cm⁻¹** inside spherical radius `R_LC` and **0.25 cm⁻¹** outside. Curvature losses are integrated along the path,

\[
\frac{d\gamma}{ds}=R_{\rm acc}-\frac{2e^2\gamma^4}{3m_ec^2\rho_c^2},\qquad
E_c=\frac{3\hbar c\gamma^3}{2\rho_c}
\]

The spherical acceleration switch is distinct from the cylindrical light cylinder. The gain is an effective energy input per laboratory path length; it is not inferred from the ideal force-free field and does not have an additional velocity projection. Initial Lorentz factor is **1,000**, and modeled transport ends at **2 R_LC**. No phase-dependent gain or manually brightened observed peak is inserted.

## What the image represents

The scene uses **0.1–100 GeV curvature-radiation energy emitted in all directions**, deposited into a three-dimensional display grid and projected along the viewing direction. It shows where this model radiates. It is not the directional specific intensity that a resolved gamma-ray camera would see. The final Reel contains the emission map and three title labels; there are no separate path guides, pulse chart, or added beaming flashes.

The palette represents characteristic photon energy and is illustrative: gamma rays do not have visible-light colors. Its fixed color stops are 0.03, 0.15, 0.7 and 3 GeV. Exposure is fixed throughout the movie. Tone mapping, display interpolation and diagonal composition are presentation choices. The image is rolled **−28°** and intentionally cropped to fill the vertical frame from edge to edge. A positive triangular angular display footprint spans one 128-spoke azimuth cell (±2.8125°), sampled at nine offsets, to represent the finite angular sampling as a continuous map. It uses one global energy-conservation correction. This interpolation changes neither the calculated particle trajectories nor their spectra or pulse profiles, and does not resolve physical structures between the calculated streams. The trilinear display voxels are projected using their Fourier-slice kernel, with a subpixel raster footprint; these are rendering operations, not physical plasma diffusion. An enlarged central marker indicates the star's position; it does not reconstruct the physical stellar surface or its thermal brightness. The black background omits unrelated sky objects.

The field and emission pattern rotate rigidly from one recorded field snapshot. This is not a time-dependent particle-in-cell simulation of the pair plasma. One physical rotation is displayed over **8 seconds**, approximately **89.5 times slower** than the adopted Vela period. Four identical cycles make the **32-second** Reel. Any diagrammatic motion or projection should be interpreted with these conventions.

The independently calculated pulse uses relativistic beaming directions and propagation delays, with arrival phase `(t − n·r/c)/P mod 1`. The outer diagnostic adopts an observer angle **65°** from the spin axis and a finite **2°** observer-latitude interval, using **1,000 phase bins** and conservative trajectory-segment fractions through finite sky/phase cells. Particle integration uses a time step of **0.002 Ω⁻¹** and a laboratory curvature derivative step of **0.0001 R_LC**. These pulse diagnostics are preserved in the package but are not drawn over the Reel. The displayed emission map itself is angle-integrated, so its total apparent brightness is not an observer pulse measurement.

## Recorded v1 outcome

The selected final snapshot is `pulsar-a75-n97-3rot-curl-radiative-noglm/snapshot-001020.npz`, SHA-256 `fea546506dad0c969559c198a36bd7e7fddc0003b66d7b268d5971aa8616d215`. The emission protocol is recorded as `n97-final-illustration-a128-r6`. Final evidence is included under `diagnostics/final-field-audit/` and `diagnostics/n97-final-illustration-a128-r6/` in the v1 package.

| Numerical diagnostic | Measured value |
| --- | ---: |
| Native RMS divergence / RMS B, `0.75 < r/R_LC < 2` | `9.206 × 10⁻⁶ R_LC⁻¹` |
| 65³ → 97³ E / B relative L2 difference in that shell | 10.27% / 9.49% |
| E / B grid difference in the high-curl current-sheet proxy | 41.03% / 43.76% |
| 97³ E / B change from rotation 2 to rotation 3 | 4.32% / 3.97% |
| Fitted final source moment, initial moment = 1 | 0.97726 |
| Final source phase offset | −1.7608° |

The small native divergence is one numerical constraint check. The grid differences, source drift and unresolved sheet remain material limitations; they are not canceled by a visually smooth loop. The frozen final snapshot is rotated for playback, so the movie does not show the measured secular source drift. A subsequent 512-spoke refinement revealed a northern spoke with disconnected open intervals, stable when the tracing step was reduced from 0.008 to 0.001 R_LC. That attempted 6,144-stream calculation was rejected before emission; no failed spokes were removed or validated result claimed. V1 retains the explicit 128-spoke, single-rim injection approximation and its 1,536 individually outgoing streams. Display interpolation must not be interpreted as that rejected higher-resolution physical model.

All **1,536 particle trajectories** reach the declared outer transport surface at **2 R_LC**, through **1,418,970 active integration segments**. The maximum path-energy-ledger residual is **8.88 × 10⁻¹⁶ erg**. About **76.15%** of modeled curvature energy is emitted beyond the cylindrical light cylinder. About **0.889%** is emitted into inward-pointing photon directions, and no occultation correction is applied. These energy fractions refer to the computed all-direction curvature emission, not a measured gamma-ray flux.

The final development spectrum comparison still has a descriptive χ² of **32,743 for 6 descriptive degrees of freedom**, after one shared spectral normalization. The pulse comparisons are also substantially discrepant. These results are retained to make the approximation's limits reviewable; the visualization does not claim observational agreement.

## Numerical and observational limits

- Native discrete-divergence checks apply to the native grid derivative. Interpolated fields and the display grid do not preserve that identity exactly
- Two resolutions, a finite domain and an enlarged source do not establish asymptotic convergence, current-sheet resolution, or independence from the stellar size and boundaries
- Curvature radiation is one selected GeV emission prescription. The model does not supply pair creation, self-consistent plasma loading, synchrotron pitch-angle evolution, or the inverse-Compton component needed for Vela's complete broadband emission
- No physical stellar occultation or photon absorption is calculated from the enlarged numerical star
- The 3PC pulse and spectral files are development comparisons. The pulse shape and spectrum remain discrepant; no successful fit or independent predictive validation is claimed
- Weighted catalog pulse counts are not ordinary Poisson counts. Current comparisons use simplified response assumptions and are descriptive, rather than a full event-level likelihood with exposure and source weights
- The retained catalog spectral FITS has invalid embedded checksum cards. Independent downloads from the official source gave the same bytes; the warning and file hash are recorded in [`data/spectra/provenance.json`](../data/spectra/provenance.json)

The visualization-release scope is separate from quantitative validation. Existing scientific flags such as `vela_validation=false` and `emission_release_allowed=false` remain meaningful: a polished video does not turn this approximate model into a validated reconstruction. Exact final scene, source and render hashes accompany the release manifest. Earlier outputs from the magnetic-blending/sponge variants are historical diagnostics and must not be substituted for the selected field.

## Studies, data and software actually used

1. **Spitkovsky (2006)**, *Time-dependent force-free pulsar magnetospheres: axisymmetric and oblique rotators*. Provides the force-free magnetosphere context and ideal oblique-dipole spin-down scaling. [Primary paper](https://arxiv.org/abs/astro-ph/0603147) · [DOI](https://doi.org/10.1086/507518)
2. **Kalapotharakos, Harding & Kazanas (2014)**, *Gamma-Ray Emission in Dissipative Pulsar Magnetospheres: From Theory to Fermi Observations*. Provides the adopted ultrarelativistic particle-direction prescription and laboratory radiation/arrival-phase framework. This release does not reproduce its full dissipative field solutions. [Primary paper](https://arxiv.org/abs/1310.3545) · [DOI](https://doi.org/10.1088/0004-637X/793/2/97)
3. **Barnard, Venter, Harding, Kalapotharakos & Johnson (2022)**, *Probing the High-energy γ-ray Emission Mechanism in the Vela Pulsar via Phase-resolved Spectral and Energy-dependent Light Curve Modeling*. Provides the prescribed acceleration scenario and curvature-radiation energy equation. The present geometry, enlarged source and finite-radius injection differ from the paper's complete model; this is not a reproduction of its Vela fit. [Primary paper](https://arxiv.org/abs/2111.03405) · [DOI](https://doi.org/10.3847/1538-4357/ac2a3d)
4. **Smith et al. (2023)**, *The Third Fermi Large Area Telescope Catalog of Gamma-ray Pulsars*. Supplies Vela timing information and the pulse/spectral data used in development comparisons. The catalog is not a three-dimensional emissivity map. [Primary paper](https://arxiv.org/abs/2307.11132) · [Official 3PC data](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/) · [DOI](https://doi.org/10.3847/1538-4357/acee67)
5. **COFFEE / CoffeeGPU**, developed by **Yajie Yuan and Alex Chen**, with repository contributors. The field calculation adapts this code at commit `9986522d08b4d887e5ef9a53e12cecae4e376d22`; it includes substantial local CPU, source and boundary changes. Upstream code and modifications retain **GPL-2.0**, with the original license and complete patch chain preserved. [Upstream source](https://github.com/fizban007/CoffeeGPU/tree/9986522d08b4d887e5ef9a53e12cecae4e376d22) · [Developer attribution](https://yuanyajie.github.io/) · [Local provenance](../third_party/coffee-curl-radiative-provenance.json)

## Visual research reference

NASA Goddard's **[Simulations Create New Insights Into Pulsars](https://svs.gsfc.nasa.gov/13058/)** (2018) was a visual and research reference: visualization by Tom Bridgman, production by Scott Wiessinger, scientific consulting by Gabriele Brambilla. Its associated study is **Brambilla et al. (2018)**, *Electron–Positron Pair Flow and Current Composition in the Pulsar Magnetosphere* ([paper](https://arxiv.org/abs/1710.03536), [DOI](https://doi.org/10.3847/1538-4357/aab3e1)). No NASA frames, simulation data, particle population or music are incorporated into this original Reel. NASA's displayed pair-flow calculation is more self-consistent than the prescribed radiation population used here, and should not be described as the same model.

## Licensing and attribution

Project code is distributed under **GNU GPL version 2**; see the repository [`LICENSE`](../LICENSE) and [`NOTICE`](../NOTICE). Upstream native code, its dependencies and bundled fonts retain their original notices and applicable licenses. Latin Modern uses the GUST Font License; STIX Two Math uses SIL Open Font License 1.1. Observational data retain their original attribution and distribution terms; public access is not presented as a new project license grant.

The original rendered movie, still images and Instagram caption are **Copyright © 2026 Loganpis. All rights reserved**. The code license does not grant a separate license to those media or caption assets. No NASA media or music license is required for this original rendering because none is incorporated.

