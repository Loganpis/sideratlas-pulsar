# SiderAtlas NO.002 — model and limits

This release depicts a plausible young rotation-powered pulsar through an approximate optical radiation model. Vela supplies useful reference scales; its luminosity, optical slope and pulse profile are **not fitting targets or release acceptance criteria**. The model is not a unique description of all pulsars or an observational reconstruction of Vela. Numerical checks support the stated equations without removing the physical assumptions below.

## Physical inputs

| Quantity | Adopted value |
|---|---|
| Period | 0.08937108859 s |
| Mass / circumferential radius | 1.4 M☉ / 13 km |
| Polar surface field | 3 × 10¹² G |
| Redshifted effective temperature | T∞ = 0.70 MK; local T_eff ≈ 0.848 MK |
| Magnetic / observer inclination to rotation axis | 75° / 65° |
| Light-cylinder radius | R_LC = cP/(2π), approximately 4264 km |
| Observer optical band | 400–700 nm |
| Activated pairs | Initial γ = 300, pitch angle 0.1 rad, both poles |
| Activation region | 0.08–0.18 R_LC, dipolar cap annulus coordinates 0.8–0.9 |
| Supply | Chosen Goldreich–Julian cap supply × 6000 pairs per primary |

These describe an energetic young pulsar rather than a population average. Temperature, field, surface gravity and gravitational redshift are used consistently across atmosphere and stellar camera. Multiplicity, activation layer and initial pitch/energy are prescribed inputs; no pair-creation mechanism or global spin-down solution is inferred from them.

## Energy consistency

With the polar-field convention, a force-free spin-down reference is 1.90×10³⁰ W (vacuum-dipole reference 6.13×10²⁹ W). The chosen total lepton kinetic supply is 2.92×10²⁵ W and its calculated local synchrotron loss is 2.83×10²⁵ W. Thus the source respects its injection-energy budget and requires only a small part of the rotational-power scale. These are analytic consistency scales, not a global spin-down solution or a fit to Vela. See `research/optical-finish/PHYSICAL-CONSISTENCY.md` for conventions and the scalar receipt

## Light calculation

The optical component uses a physical dipole in the near zone, below approximately 0.2 R_LC. Guiding centres follow it with adiabatic focusing and classical gyro-averaged radiation reaction. An energy ledger tracks cooling losses. The integration continues into the low transverse-momentum regime, where a smooth synchrotron continuum would be inaccurate for optical photons.

Directional classical helical harmonics are evaluated with Bessel functions. With A = 1 − β∥ cos θ, line frequencies obey ν_n = n eB/(2πγm_e A). Physical lepton-number weights are summed incoherently. A local Lorentz transformation accounts for a rotating emission cell's photon Doppler shift and aberration. Field lines, species-colored tracks and arbitrary luminous cones are not substituted for optical photons.

A second component comes from a calculated two-mode LTE hydrogen atmosphere with hydrostatic balance, magnetic electron–proton free–free absorption and coherent mode-changing scattering. Its intensity is integrated through the lensed stellar image. It is neither a blackbody replacement nor an extrapolation of an X-ray-only atmosphere table.

The surface camera directly integrates outward Schwarzschild rays, relative delays, redshift and local rotational Doppler/aberration. Wavelength intensity transforms as Iλ,obs = (gD)⁵ Iλ,em at the corresponding emitted wavelength. Plasma photons use flat-space retardation τ = t − n·r/c. Source-time energy is conservatively deposited in arrival intervals; an extra observer-time beaming factor is not applied afterward.

## Separate PIC research

The repository also preserves a native coupled particle/Maxwell calculation with source, charge and feedback audits. Its enlarged star and reduced scales make it a separate numerical study. Its Lorentz factors, fields and charges are **not rescaled into physical optical photons**. The compatible conducting-source repair improves that study's interior drift but does not validate the optical injection population.

The released optical motion and emission are calculated in a prescribed physical dipole. They are not a self-consistent global PIC optical solution: global current closure, evolving field feedback and pair cascades are absent.

## Material approximations

- **Field and population:** static dipole in a local magnetic frame; rigidly rotating emission pattern. The near-zone restriction motivates this approximation, but finite-radius sweepback, current-sheet light and force-free/PIC distortions are absent. The supply is not derived from relativistic continuity or current closure
- **Particle physics:** monoenergetic activation and prescribed pitch, with no pair cascade, resonant radio heating or inverse Compton calculation. Static-dipole curvature power is negligible for these recorded cooling columns; this is not a claim about all particles throughout real pulsars
- **Propagation:** optically thin plasma, with no absorption, scattering, extinction or general-relativistic plasma-photon transport. Relativistic surface transport is separate; their approximate common phase convention is not calibrated pulsar timing
- **Atmosphere:** one plane-parallel homogeneous normal-field column used uniformly over the sphere, fully ionized ideal hydrogen, transverse linear modes and unity magnetic Gaunt factors. Atomic opacity/partial ionization, dispersion/refraction, vacuum mode conversion, Compton redistribution, conduction and return-current irradiation are omitted. Neutral-fraction audits show that complete ionization is substantive. Independent magnetic optical-spectrum validation is absent, so physical opacity errors are not bounded
- **Star:** spherical Schwarzschild plus local Doppler/aberration, without frame dragging or oblateness. The uniform normal-field column is not a latitude-dependent dipolar atmosphere library
- **Sampling:** 5,271,172 nonzero production events, 256 source phases, 128 azimuths, 64 cooling nodes, 5 activation radii and 2 cap-annulus nodes. Uniform emitting-cell rasterization is finite spatial quadrature, not resolved microscopic structure. Narrow harmonics and pulse/image details retain sampling sensitivity

## Numerical and display interpretation

The selected 64-depth atmosphere satisfies its numerical equation gates: maximum flux error approximately 0.095%, local relative balance residual 4.86 × 10⁻⁵ and conservative ledger error 2.47 × 10⁻⁷ of the target flux. Optical output changes by approximately 5% from 32 depths. Its temperature profile is not everywhere monotone. These support the approximate transfer problem, **not complete atmosphere validation**; stored scientific-validation and refinement flags remain false.

Transport checks cover energy accounting, small quantum parameter and guiding-centre scale separation. Existing sampling controls and historical Vela comparisons remain research records. Differences from Vela do not fail this generic release; no physical normalization or pulse is tuned to hide them. Further expensive fitting or refinement is outside the clarified visualization scope.

The 1080 × 1920 output covers a 108 × 192 km physical image plane at 100 m per native pixel, with a −10° roll. The stellar image is calculated at native scale. Smooth finite-cell plasma radiance uses a 135 × 240 grid at 800 m per pixel, then conservative bilinear resampling onto the shared native image. Plasma structure below that grid is not resolved. A normalized 0.6-native-pixel PSF is applied after adding both layers.

Cooling cells are grouped along their ordered path at the same source phase if their nodes and eight support corners remain within 200 m. The group sums already-calculated `XYZ * number_weights` exactly; no angular or spectral averaging is added. Group geometry is approximate, with actual displacement and timing bounds recorded. Frozen-motion pieces have an arc limit of two plasma-grid pixels, equivalent to 16 native pixels. Source/arrival overlap is analytic along each static ray for the chosen 2×2 plasma-pixel quadrature. Wholly offscreen volumes do not affect the image; their cropped-energy timing ledger uses 4×4×4 Gauss nodes rather than exact spatial timing. Covering arrival bins retain the full source-energy ledger.

Spectra use official CIE1931 2° color functions and one fixed enhanced exposure: every linear RGB value is divided by half the maximum positive surface RGB radiance before the common intensity tone map and sRGB encoding. This leaves the layers' linear light ratios unchanged before tone mapping. Negative out-of-gamut linear RGB components are clipped. Video pixels are transformed to tagged Rec.709. This is a spectral optical basis, not a naked-eye brightness/adaptation or particular telescope prediction.

One physical rotation is stretched to eight seconds and repeated four times: constructed slow motion, not 32 seconds of physical evolution. A shared physical image scale fills the portrait canvas and can crop edge emission; the star is not enlarged independently. Only the three title labels appear. The final export receipt records the delivered movie properties and repeated-loop checks.

## Sources and roles

- Spitkovsky 2006 — force-free spin-down energy scale: https://arxiv.org/abs/astro-ph/0603147; DOI 10.1086/507518
- Smith et al. 2023, third Fermi-LAT pulsar catalog — period/catalog reference: https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/
- Harding et al. 2018 — pair optical-emission and reference geometry: https://arxiv.org/abs/1811.11157; DOI 10.3847/2041-8213/aaf3b2
- Harding & Kalapotharakos 2015 — pair-population context; cascade/radio heating not reproduced: https://arxiv.org/abs/1508.06251; DOI 10.1088/0004-637X/811/1/63
- Marcowith & Malzac 2003 — classical helical harmonic formulation: https://arxiv.org/abs/astro-ph/0309107; DOI 10.1051/0004-6361:20031060
- Lai et al. 2025 — directional-radiation benchmark context, not a reproduced global model: https://arxiv.org/abs/2509.04572; DOI 10.1093/mnras/staf1295
- van Adelsberg & Lai 2006 — magnetic opacity/transfer basis, not their complete solver: https://arxiv.org/abs/astro-ph/0607168
- Ofengeim & Zyuzin 2018 — chosen mass/radius/thermal reference: https://www.ioffe.ru/astro/Stars/Paper/ofengeim_zyuzin18.pdf
- Poutanen & Beloborodov 2006 — relativistic surface radiation: DOI 10.1111/j.1365-2966.2006.11088.x
- Poutanen 2020 — Schwarzschild light bending: DOI 10.1051/0004-6361/202037471
- Mignani, Zharikov & Caraveo 2007 — historical optical comparison, not a fit target: https://arxiv.org/abs/0707.2036; DOI 10.1051/0004-6361:20077774
- Ioffe magnetic-H EOS tables — independent neutral-fraction audit, not optical atomic opacity: https://www.ioffe.ru/astro/NSG/Hmagnet/hmagtab.html
- CIE1931 2° observer data — spectral color: https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer; DOI 10.25039/CIE.DS.xvudnb9b; CC BY-SA 4.0 attribution retained
- Philippov & Spitkovsky 2014 — PIC/source-boundary context: DOI 10.1088/2041-8205/785/2/L33
- Belyaev 2015 — conducting-star/PIC context: https://academic.oup.com/mnras/article/449/3/2759/1132376
- NASA / Goddard 2018 — visual motion/context: https://svs.gsfc.nasa.gov/13058/; no NASA footage or code used

Unlicensed McPHAC vendor files and external Ioffe downloads are excluded from publication; source URLs, hashes and derived audit results remain where applicable. NSMAXG and Rosseland means are not local angle-dependent optical surface inputs.
