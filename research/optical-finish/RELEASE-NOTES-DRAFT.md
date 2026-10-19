> Historical draft from the earlier Vela-reconstruction scope, superseded by release/optical-v1/MODEL-NOTES.md and the final caption

# Optical release notes — working draft

This is a caption and citation inventory for editing, **not release approval**. The atmosphere is still being solved numerically. No accepted combined Reel or Vela reconstruction is asserted here. Recheck every status against the final run receipts before publishing

## Caption draft

The following block is 1912 characters including line breaks. Numerical mismatch statements refer specifically to the recorded `baseline-v2` photon model and should be recalculated if its parameters or sampling change

```text
NO.002 — THE PULSAR
@sideratlas
Original optical model — Vela-informed, not a Vela reconstruction
P = 89.37 ms · R = 13 km · B_p = 3×10¹² G
Magnetic inclination 75° · Observer 65°
Magnetic fields curve charged-particle motion
Accelerated pairs emit optical harmonics
Radiation cools their transverse motion
Rotation directs pulsed light toward the observer
ν_n = n eB/[2πγm_e(1−β∥cosθ)]
τ = t−n·r/c · Visible band: 400–700 nm
CALCULATED
Directional classical harmonic emission, adiabatic focusing and radiative cooling
Photon Doppler shifts, aberration and relative arrival times
CIE1931 spectral color with fixed display exposure and sRGB tone mapping
Surface light bending and shifts implemented; atmosphere still being solved
ASSUMED
Static dipole and rigidly rotating emission pattern, separate from the reduced PIC calculation
Pairs activated at 0.08–0.18 R_LC with γ=300 and pitch=0.1 rad
6000 pairs per primary; cap annulus 0.8–0.9; both poles
Optically thin plasma; no cascade or radio-heating calculation
No plasma absorption, extinction, inverse Compton or curvature loss
Surface: uniform normal-field, fully ionized H; unity magnetic Gaunt factors
No atomic opacity, vacuum mode conversion or dipolar surface library
LIMITS
Baseline approximate 550 nm flux ≈16× Vela at assumed 290 pc
α≈−0.56 versus −0.04±0.04 observed; F_ν∝ν^(−α)
Different spectral intervals/band conventions; no fitted normalization
No converged or validated combined optical reconstruction yet
SOURCES
Harding et al. 2018; Harding & Kalapotharakos 2015
Marcowith & Malzac 2003; van Adelsberg & Lai 2006
Poutanen & Beloborodov 2006; Poutanen 2020
Mignani et al. 2007; Ofengeim & Zyuzin 2018; CIE1931
Motion/context reference: NASA / Goddard 2018
https://svs.gsfc.nasa.gov/13058/
No NASA footage or code used
Full citations, equations, assumptions and numerical checks in the repository
https://github.com/loganpis/sideratlas-pulsar
```

## Source and assumption inventory

| Component | Actually calculated / implemented | Prescribed inputs and material limits | Primary basis |
|---|---|---|---|
| Optical plasma | Classical Schott/Bessel harmonic line power in the observer direction; physical number weights; guiding-centre adiabatic focusing and radiation-reaction cooling; local SR photon boost | Static physical dipole, rigidly rotating emission cells, gamma 300, pitch .1, uniform activation .08–.18 RLC, cap annulus .8–.9, 6000 pairs/primary, both poles; not a cascade, radio-absorption solve, or PIC-derived physical field | Harding et al. 2018; Harding & Kalapotharakos 2015 motivate geometry/population; Marcowith & Malzac 2003 supplies harmonic formulation; Lai et al. 2025 is an additional directional-radiation reference listed by the component |
| Photon camera | Conserved source-time energy deposited through rotating-pattern arrival intervals; fixed observer; shared physical pixel scale; optical XYZ integration | Flat propagation for magnetospheric emission, optically thin, no extinction, absorption or GR plasma transport; finite quadrature/pixels, CIC rasterization and .6-pixel display PSF; fixed gain, common RGB tone map and clipping of negative out-of-gamut channels; no naked-eye exposure claim | Component equations and camera implementation; official CIE1931 2° table |
| Surface | Direct outward Schwarzschild ray integrals, relative delays, local SR Doppler/aberration and I_lambda shift factor (gD)^5 implemented | M=1.4 Msun, R=13 km, spherical Schwarzschild-plus-Doppler surface; no frame dragging/oblateness; a uniform normal-B column cannot represent the full dipolar atmosphere; surface/plasma phase alignment is not a validated pulse prediction | Poutanen & Beloborodov 2006; Poutanen 2020 |
| Magnetic atmosphere | Local polarized transfer, hydrostatic ideal-H column, mode-changing coherent scattering and broad-band equilibrium iterations are being solved; conservative flux/heating audit is separate from optimizer exit | B=3e12 G, T_eff≈847655 K, g≈1.3313e14 cm/s²; fully ionized ideal H, transverse linear modes, unity magnetic Gaunt factors; omit bound opacity, dispersion/refraction, vacuum conversion, Compton redistribution, conduction/irradiation; provisional-column neutral-fraction audit found non-negligible atoms; no accepted independent magnetic optical spectrum or refinement result yet | van Adelsberg & Lai 2006 §3.1.2 for corrected electron–proton response; author Ioffe magnetic-H EOS tables for the independent neutral-fraction audit; Ofengeim & Zyuzin 2018 supplies a chosen Vela thermal prior, not an observed surface map |
| Reduced PIC research | Native coupled field/particle updates; compatible rotating-conductor E; source-core, charge-export, feedback and 48³/96³ controls | Separate numerical study; enlarged star .5 RLC, reduced kinetic scales, prescribed injection, startup and boundary limits; no conversion of its gamma/charge/field values into physical optical emitters; surface phase error remains about 2.4° at 96³ | Entity/Kokkos pinned receipts, source-repair report; Philippov & Spitkovsky 2014 and Belyaev 2015 for conducting-star modelling context |

The photon model and reduced PIC model must remain visibly distinct in the repository. The successful interior source repair does not validate the separate prescribed photon population or remove its observed flux/slope discrepancy

## Observational mismatch to retain

`outputs/optical-finish/particle-emission/baseline-v2/independent-optical-audit.json` records approximate unextincted 550 nm Fnu = 2.5105e−5 Jy at 290 pc, versus a rough dereddened reference 1.5277e−6 Jy, ratio 16.4336. This compares a monochromatic AB-like model quantity with a Johnson-V-derived estimate, so it is an order-of-magnitude check. Its model alpha = −.55649 over 400–700 nm contrasts with Mignani et al.'s −.04±.04 over 400–800 nm. Neither comparison is a data likelihood. Display exposure cannot be used to conceal these physical mismatches

## Full citation links

- Harding et al. 2018, https://arxiv.org/abs/1811.11157, DOI 10.3847/2041-8213/aaf3b2
- Harding & Kalapotharakos 2015, https://arxiv.org/abs/1508.06251, DOI 10.1088/0004-637X/811/1/63
- Marcowith & Malzac 2003, https://arxiv.org/abs/astro-ph/0309107, DOI 10.1051/0004-6361:20031060
- Lai et al. 2025, https://arxiv.org/abs/2509.04572, DOI 10.1093/mnras/staf1295 — retain the specific benchmark role rather than imply its model was reproduced
- van Adelsberg & Lai 2006, https://arxiv.org/abs/astro-ph/0607168 — approximate opacity response only; not their complete atmosphere solver
- Poutanen & Beloborodov 2006, DOI 10.1111/j.1365-2966.2006.11088.x
- Poutanen 2020, DOI 10.1051/0004-6361/202037471
- Mignani, Zharikov & Caraveo 2007, https://arxiv.org/abs/0707.2036, DOI 10.1051/0004-6361:20077774
- Ofengeim & Zyuzin 2018, https://www.ioffe.ru/astro/Stars/Paper/ofengeim_zyuzin18.pdf
- Ioffe author magnetic-H EOS tables/documentation, https://www.ioffe.ru/astro/NSG/Hmagnet/hmagtab.html — EOS/neutral fractions, not optical atomic opacities; retain recorded table hashes and redistribution review
- CIE1931 2° colour-matching functions, https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer, DOI 10.25039/CIE.DS.xvudnb9b; CC BY-SA4.0 with attribution retained
- Philippov & Spitkovsky 2014, DOI 10.1088/2041-8205/785/2/L33; Belyaev 2015, https://academic.oup.com/mnras/article/449/3/2759/1132376 — source boundary context, not validation of this implementation
- NASA/Goddard 2018, https://svs.gsfc.nasa.gov/13058/ — motion/context reference; no NASA footage or code used

McPHAC is a separate **nonmagnetic** native baseline, not the visible surface input and not an independent passed magnetic benchmark. Its ultraviolet-limited output is not extrapolated into the optical. Keep unlicensed vendor source out of publication pending review. NSMAXG and Rosseland-mean tables are likewise not optical intensity inputs
