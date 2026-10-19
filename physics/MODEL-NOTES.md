# Physical foundation of the calculated gamma-ray scene

These are the physical assumptions, numerical outcomes and citations for the selected force-free/particle/curvature-radiation calculation. Current rendering conventions are in [the animation notes](../docs/MODEL-NOTES.md). The saved photon scene and original field/particle data are reused without a new physical solve

## Physical model

The field is the 97³, three-rotation `curl-radiative / no-GLM` candidate documented in [`research/COFFEE-CURL-REPAIR.md`](research/COFFEE-CURL-REPAIR.md). It evolves a finite force-free electromagnetic field around a rotating conductor, with magnetic inclination **75°**. The model extends to a cube of half-width **3 light-cylinder radii**; the numerical star radius is **0.5 light-cylinder radii**. The light cylinder is the cylindrical distance at which rigid corotation would reach the speed of light.

The initial magnetic field is a fourth-order discrete curl of a finite vector potential. The rotating conductor imposes the electric relation `E = −(Ω × r) × B`; it does not reset evolved magnetic components to an analytic dipole. The outer electric boundary damps an incoming radiation characteristic without directly damping magnetic components. Magnetic-gradient GLM correction is disabled in this selected candidate. These choices improve the native magnetic constraint, but do not establish global physical convergence. The finite conductor has measured moment decay and phase lag, the transition spans two grid cells, and the outer boundary is an approximation with nonzero oblique-wave reflection.

The period **0.08937108859019084 s** and its catalog derivative set the dimensional scale, using an assumed moment of inertia of `10^45 g cm²` and an ideal force-free oblique-point-dipole torque normalization. This gives a light-cylinder radius of approximately **4,264.2 km**. An enlarged finite numerical star need not share the exact torque of the ideal point-dipole reference.

The selected run injects **1,536 modeled streams** at **0.7 R_LC**, from **128 azimuths × 6 rings × 2 hemispheres** near a sampled finite-domain cap rim. The annulus fractions are **0.90–0.96** on that injection surface. This differs from injection at a physical neutron-star surface. The outward velocity branch is selected explicitly, and seed weights follow normalized magnetic flux; a self-consistent electron/positron supply is not predicted.

The selected run explicitly opts into a **finite-domain boundary-open** criterion: a magnetic trace either crosses the cylindrical light cylinder or escapes outward through a sphere of radius **2.5 R_LC**. This does not prove that a line opens at infinity. Of the selected annular magnetic traces, **1,376** cross the light cylinder and **160** reach the outgoing sphere first. Those 160 retain their outgoing classification when the magnetic-trace step is reduced from **0.008 to 0.002 R_LC**; their radial direction cosines at escape span **0.6777–0.9651**. The same criterion is applied to the sampled rim and selected seeds, with no silently discarded spokes. The code's default stricter light-cylinder criterion remains available. These magnetic traces define injection; they are distinct from the subsequent particle trajectories.

Particle directions use the ultrarelativistic guiding-center prescription of Kalapotharakos, Harding & Kazanas (2014), equation 12. Laboratory trajectory curvature is calculated from the evolving position and direction; it is not replaced with magnetic-field-line curvature. This approximation does not resolve gyration, solve the complete Lorentz dynamics, or feed particle energy and pair currents back into the electromagnetic field.

The selected Barnard et al. acceleration scenario prescribes `dγ/ds` gains of **0.04 cm⁻¹** inside spherical radius `R_LC` and **0.25 cm⁻¹** outside. Curvature losses are integrated along the path,

\[
\frac{d\gamma}{ds}=R_{\rm acc}-\frac{2e^2\gamma^4}{3m_ec^2\rho_c^2},\qquad
E_c=\frac{3\hbar c\gamma^3}{2\rho_c}
\]

The spherical acceleration switch is distinct from the cylindrical light cylinder. The gain is an effective energy input per laboratory path length; it is not inferred from the ideal force-free field and does not have an additional velocity projection. Initial Lorentz factor is **1,000**, and modeled transport ends at **2 R_LC**. No phase-dependent gain or manually brightened observed peak is inserted.

## Selected calculation and numerical checks

The selected final snapshot is `pulsar-a75-n97-3rot-curl-radiative-noglm/snapshot-001020.npz` (included under `inputs/`), SHA-256 `fea546506dad0c969559c198a36bd7e7fddc0003b66d7b268d5971aa8616d215`. The emission protocol is recorded as `n97-final-illustration-a128-r6`. Final evidence is included under `diagnostics/final-field-audit/` and `diagnostics/n97-final-illustration-a128-r6/` in this folder.

| Numerical diagnostic | Measured value |
| --- | ---: |
| Native RMS divergence / RMS B, `0.75 < r/R_LC < 2` | `9.206 × 10⁻⁶ R_LC⁻¹` |
| 65³ → 97³ E / B relative L2 difference in that shell | 10.27% / 9.49% |
| E / B grid difference in the high-curl current-sheet proxy | 41.03% / 43.76% |
| 97³ E / B change from rotation 2 to rotation 3 | 4.32% / 3.97% |
| Fitted final source moment, initial moment = 1 | 0.97726 |
| Final source phase offset | −1.7608° |

The small native divergence is one numerical constraint check. The grid differences, source drift and unresolved sheet remain material limitations; they are not canceled by a visually smooth loop. The frozen final snapshot is rotated for playback, so the movie does not show the measured secular source drift. A subsequent 512-spoke refinement revealed a northern spoke with disconnected open intervals, stable when the tracing step was reduced from 0.008 to 0.001 R_LC. That attempted 6,144-stream calculation was rejected before emission; no failed spokes were removed or validated result claimed. The selected calculation retains the explicit 128-spoke, single-rim injection approximation and its 1,536 individually outgoing streams. Display interpolation must not be interpreted as that rejected higher-resolution physical model.

All **1,536 particle trajectories** reach the declared outer transport surface at **2 R_LC**, through **1,418,970 active integration segments**. The maximum path-energy-ledger residual is **8.88 × 10⁻¹⁶ erg**. About **76.15%** of modeled curvature energy is emitted beyond the cylindrical light cylinder. About **0.889%** is emitted into inward-pointing photon directions, and no occultation correction is applied. These energy fractions refer to the computed all-direction curvature emission, not a measured gamma-ray flux.

The final development spectrum comparison still has a descriptive χ² of **32,743 for 6 descriptive degrees of freedom**, after one shared spectral normalization. The pulse comparisons are also substantially discrepant. These results are retained to make the approximation's limits reviewable; the visualization does not claim observational agreement.

## Numerical and observational limits

- Native discrete-divergence checks apply to the native grid derivative. Interpolated fields and the display grid do not preserve that identity exactly
- Two resolutions, a finite domain and an enlarged source do not establish asymptotic convergence, current-sheet resolution, or independence from the stellar size and boundaries
- Curvature radiation is one selected GeV emission prescription. The model does not supply pair creation, self-consistent plasma loading, synchrotron pitch-angle evolution, or the inverse-Compton component needed for Vela's complete broadband emission
- No physical stellar occultation or photon absorption is calculated from the enlarged numerical star
- The 3PC pulse and spectral files are development comparisons. The pulse shape and spectrum remain discrepant; no successful fit or independent predictive validation is claimed
- Weighted catalog pulse counts are not ordinary Poisson counts. Current comparisons use simplified response assumptions and are descriptive, rather than a full event-level likelihood with exposure and source weights
- The retained catalog spectral FITS has invalid embedded checksum cards. Independent downloads from the official source gave the same bytes; the warning and file hash are recorded in [`data/spectra/provenance.json`](data/spectra/provenance.json)

The visualization-release scope is separate from quantitative validation. Existing scientific flags such as `vela_validation=false` and `emission_release_allowed=false` remain meaningful: a polished video does not turn this approximate model into a validated reconstruction. Exact final scene, source and render hashes accompany the release manifest. Earlier outputs from the magnetic-blending/sponge variants are historical diagnostics and must not be substituted for the selected field.

## Studies, data and software actually used

1. **Spitkovsky (2006)**, *Time-dependent force-free pulsar magnetospheres: axisymmetric and oblique rotators*. Provides the force-free magnetosphere context and ideal oblique-dipole spin-down scaling. [Primary paper](https://arxiv.org/abs/astro-ph/0603147) · [DOI](https://doi.org/10.1086/507518)
2. **Kalapotharakos, Harding & Kazanas (2014)**, *Gamma-Ray Emission in Dissipative Pulsar Magnetospheres: From Theory to Fermi Observations*. Provides the adopted ultrarelativistic particle-direction prescription and laboratory radiation/arrival-phase framework. This release does not reproduce its full dissipative field solutions. [Primary paper](https://arxiv.org/abs/1310.3545) · [DOI](https://doi.org/10.1088/0004-637X/793/2/97)
3. **Barnard, Venter, Harding, Kalapotharakos & Johnson (2022)**, *Probing the High-energy γ-ray Emission Mechanism in the Vela Pulsar via Phase-resolved Spectral and Energy-dependent Light Curve Modeling*. Provides the prescribed acceleration scenario and curvature-radiation energy equation. The present geometry, enlarged source and finite-radius injection differ from the paper's complete model; this is not a reproduction of its Vela fit. [Primary paper](https://arxiv.org/abs/2111.03405) · [DOI](https://doi.org/10.3847/1538-4357/ac2a3d)
4. **Smith et al. (2023)**, *The Third Fermi Large Area Telescope Catalog of Gamma-ray Pulsars*. Supplies Vela timing information and the pulse/spectral data used in development comparisons. The catalog is not a three-dimensional emissivity map. [Primary paper](https://arxiv.org/abs/2307.11132) · [Official 3PC data](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/) · [DOI](https://doi.org/10.3847/1538-4357/acee67)
5. **COFFEE / CoffeeGPU**, developed by **Yajie Yuan and Alex Chen**, with repository contributors. The field calculation adapts this code at commit `9986522d08b4d887e5ef9a53e12cecae4e376d22`; it includes substantial local CPU, source and boundary changes. Upstream code and modifications retain **GPL-2.0**, with the original license and complete patch chain preserved. [Upstream source](https://github.com/fizban007/CoffeeGPU/tree/9986522d08b4d887e5ef9a53e12cecae4e376d22) · [Developer attribution](https://yuanyajie.github.io/) · [Local provenance](third_party/coffee-curl-radiative-provenance.json)

## Visual research reference

NASA Goddard's **[Simulations Create New Insights Into Pulsars](https://svs.gsfc.nasa.gov/13058/)** (2018) was a visual and research reference: visualization by Tom Bridgman, production by Scott Wiessinger, scientific consulting by Gabriele Brambilla. Its associated study is **Brambilla et al. (2018)**, *Electron–Positron Pair Flow and Current Composition in the Pulsar Magnetosphere* ([paper](https://arxiv.org/abs/1710.03536), [DOI](https://doi.org/10.3847/1538-4357/aab3e1)). No NASA frames, simulation data, particle population or music are incorporated into this original Reel. NASA's displayed pair-flow calculation is more self-consistent than the prescribed radiation population used here, and should not be described as the same model.

