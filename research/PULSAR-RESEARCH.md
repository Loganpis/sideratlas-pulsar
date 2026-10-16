# NO.002 — A scientifically grounded pulsar

Research brief for @sideratlas · 12 October 2026

## Recommendation

Build an **obliquely rotating, plasma-filled pulsar magnetosphere**, with a swept-back magnetic structure and a rippling current sheet, and calculate a pulse profile for a specified distant observer. Use **Vela as the initial high-energy observational benchmark**. If visible-light color becomes the main objective, use Crab instead. If the star's surface is the subject, choose a thermal X-ray model such as PSR J0030+0451.

This is a design recommendation drawn from the studies below, not a claim that one universal pulsar model has been established. A useful image must identify whether it represents radiation reaching a camera, plasma density, or emitted power integrated over directions. Those quantities produce different pictures.

The two complementary anchors are **Brambilla et al. (2018)** for the particle flows in the NASA visualization supplied by the user, and **Cerutti, Figueiredo & Dubus (2025)** for synthetic radiation and observational comparisons. They are supported by force-free geometry, independent emission models, and surface-transfer benchmarks. The preferred route to maximum physical fidelity is rendering verified three-dimensional kinetic simulation outputs. A smaller original calculation can use a solved force-free field with an explicit emission prescription, but it must retain that qualification.

## Essential reading and its role

| Priority | Primary study | Contribution to our calculation | Limit to preserve |
| --- | --- | --- | --- |
| 1 | **Cerutti, Figueiredo & Dubus (2025)** — *Synthetic pulsar light curves from global kinetic simulations and comparison with the Fermi-LAT catalog* · [paper](https://arxiv.org/abs/2412.02307) · [full text](https://arxiv.org/html/2412.02307v2) · [journal](https://doi.org/10.1051/0004-6361/202451948) | Three-dimensional particle-in-cell (PIC) radiation, current-sheet geometry, energy-dependent pulses, and comparison with observations | Reduced scale separation: light cylinder at five stellar radii, simplified pair injection, no GR; many profile fits remain statistically poor |
| 2 | **Spitkovsky (2006)** — *Time-dependent force-free pulsar magnetospheres: axisymmetric and oblique rotators* · [paper](https://arxiv.org/abs/astro-ph/0603147) · [DOI](https://doi.org/10.1086/507518) | Plasma-filled, magnetically dominated structure, open/closed field regions, current sheet, and spin-down benchmark | Ideal force-free fields alone do not predict particle acceleration or brightness |
| 3 | **Cerutti, Philippov & Spitkovsky (2016)** — *Modelling high-energy pulsar light curves from first principles* · [paper](https://arxiv.org/abs/1511.01785) · [DOI](https://doi.org/10.1093/mnras/stw124) | Particle radiation, beaming, arrival-time sorting, and formation of sharp pulses | Injected plasma and compressed physical scales; distinguish angle-integrated emission from radiation reaching one observer |
| 4 | **Kalapotharakos et al. (2023)** — *The Gamma-Ray Pulsar Phenomenology in View of 3D Kinetic Global Magnetosphere Models* · [paper](https://arxiv.org/abs/2303.04054) | Independent current-sheet models and gamma-ray observational comparisons, consistent with curvature radiation | Plasma injection changes the result; synchrotron versus curvature dominance remains model dependent |
| 5 | **Bai & Spitkovsky (2010)** — *Modeling of Gamma-Ray Pulsar Light Curves with Force-Free Magnetic Field* · [paper](https://arxiv.org/abs/0910.5741) · [DOI](https://doi.org/10.1088/0004-637X/715/2/1282) | Practical link between electromagnetic geometry, emission zones, and observer pulse profiles | Emissivity is prescribed, not derived from a kinetic plasma solution |
| 6 | **Bogdanov et al. (2019)** — *Constraining the Neutron Star Mass-Radius Relation and Dense Matter Equation of State with NICER. II. Emission from Hot Spots on a Rapidly Rotating Neutron Star* · [paper](https://arxiv.org/abs/1912.05707) | Relativistic surface transfer and independent pulse-profile benchmarks | Applies to surface X-rays, not an entire gamma-ray magnetosphere |
| 7 | **Smith et al. (2023)** — *The Third Fermi Large Area Telescope Catalog of Gamma-Ray Pulsars* · [paper](https://arxiv.org/abs/2307.11132) · [data](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/) | Measured profiles, spectra, and timing models for testing predictions | These are observations of unresolved sources, not three-dimensional plasma maps |

## NASA reference supplied by the user

**[NASA SVS 13058 — Simulations Create New Insights Into Pulsars](https://svs.gsfc.nasa.gov/13058/)** is a strong reference for our educational presentation. It visualizes simulated macroparticle trajectories: blue denotes electrons, red positrons, and lighter colors higher particle energies. These are explanatory colors, not a visible-light spectrum. A particle-flow rendering and an observer's radiation image are different deliverables.

Its linked primary study is **Brambilla, Kalapotharakos, Timokhin, Harding & Kazanas (2018), *Electron–Positron Pair Flow and Current Composition in the Pulsar Magnetosphere***. [Paper](https://arxiv.org/abs/1710.03536) · [DOI](https://doi.org/10.3847/1538-4357/aab3e1). It follows electron–positron plasma injected near the stellar surface and examines how its currents build the magnetosphere. Different injection prescriptions can give similar electromagnetic structures yet different particle densities. Therefore, matching a field silhouette is insufficient to reproduce the depicted particle flow.

Numerical limits include $R_\star=0.36R_{\rm LC}$, a prescribed pair supply, compressed particle energies, and enhanced radiation reaction. The paper cautions against simply rescaling its particle-energy distribution to real pulsars. Reproducing its appearance would inherit these assumptions.

NASA provides useful separate references:

| Component | Official reference | Use for NO.002 |
| --- | --- | --- |
| Vacuum field | [SVS 4637](https://svs.gsfc.nasa.gov/4637/) | Comparison case illustrating the missing plasma feedback |
| Plasma-modified field | [SVS 4638](https://svs.gsfc.nasa.gov/4638/) | Field structure including currents carried by the particles |
| Bulk trajectories | [SVS 4644](https://svs.gsfc.nasa.gov/4644/) | Motion of the lower-energy particle population |
| Combined particle populations | [SVS 4648](https://svs.gsfc.nasa.gov/4648/) | Camera tours, polar views, and distinctions between stationary and rotating viewpoints |

The combined-flow page offers 4K movies and rendered frame sequences, but explicitly says SVS does not store the underlying model dataset. Those frames cannot substitute for three-dimensional particle and field data. The generic model should not be relabeled as a uniquely determined reconstruction of Vela.

**Recommended application:** a labeled scientific visualization of the rotating plasma environment, accompanied by a separately calculated pulse profile for a fixed observer. Use NASA's treatment to explain particle motion and later emission studies to calculate the signal. Preserve the distinction between particle energy, particle density, and radiated intensity; do not use one as an undocumented proxy for another.

## Further studies that change implementation decisions

**Magnetospheric structure and plasma**

- **Goldreich & Julian (1969), *Pulsar Electrodynamics*.** [Original paper](https://articles.adsabs.harvard.edu/pdf/1969ApJ...157..869G). Establishes why the rotating magnetized star is surrounded by charges and introduces the corotation charge-density scale. Its aligned idealization is a foundation, not an emission map.
- **Bogovalov (1999), *On the physics of cold MHD winds from oblique rotators*.** [Paper](https://arxiv.org/abs/astro-ph/9907051). Supplies an analytic inclined split-monopole wind and corrugated current-sheet geometry. Use outside the inner dipolar magnetosphere, with explicit radial validity and an independent emissivity prescription.
- **Li, Spitkovsky & Tchekhovskoy (2012), *Resistive Solutions for Pulsar Magnetospheres*.** [Paper](https://arxiv.org/abs/1107.0979). Conductivity bridges vacuum and force-free regimes and permits non-ideal electric fields. Conductivity is an added physical assumption that must be specified.
- **Philippov & Spitkovsky (2018), *Ab-Initio Pulsar Magnetosphere: Particle Acceleration in Oblique Rotators and High-energy Emission Modeling*.** [Paper](https://arxiv.org/abs/1707.04323). Includes pair creation, ion extraction, frame-dragging corrections, and photon propagation. Important for linking near-surface plasma supply to outer emission.
- **Soudais, Cerutti & Contopoulos (2024), *Scaling up global kinetic models of pulsar magnetospheres using a hybrid force-free-PIC numerical approach*.** [Paper](https://arxiv.org/abs/2406.14512) · [full text](https://arxiv.org/html/2406.14512v1). Addresses unrealistic scale compression. Its realistic-scale example is an aligned, axisymmetric weak millisecond pulsar; it cannot be used unchanged as an oblique rotating lighthouse.

**Radiation and optical appearance**

- **Bai & Spitkovsky (2010), *Uncertainties of Modeling Gamma-Ray Pulsar Light Curves with Vacuum Dipole Magnetic Field*.** [Paper](https://arxiv.org/abs/0910.5740). Identifies a frame inconsistency in older aberration prescriptions. A laboratory-frame photon direction must not be boosted a second time as though it were a corotating-frame direction.
- **Philippov, Timokhin & Spitkovsky (2020), *Origin of Pulsar Radio Emission*.** [Paper](https://arxiv.org/abs/2001.02236). Kinetic calculations demonstrate a candidate coherent-emission process from intermittent, nonuniform pair discharges. This is not a complete agreed model of every pulsar's radio beam.
- **Cerutti, Mortier & Philippov (2016), *Polarized synchrotron emission from the equatorial current sheet in gamma-ray pulsars*.** [Paper](https://arxiv.org/abs/1609.00021). Phase-dependent polarization offers an additional diagnostic of the emission geometry. Their Crab-like viewing and magnetic angles are model choices, not universal measured parameters.
- **Sollerman et al. (2019), *The optical and near-infrared spectrum of the Crab pulsar with X-shooter*.** [Paper](https://arxiv.org/abs/1906.09878). A flux-calibrated 300–2400 nm spectrum supports a spectrum-derived visible color after declaring extinction and exposure. The reported continuum is approximately $F_\nu\propto\nu^{0.16\pm0.07}$; it is spatially unresolved and phase averaged.
- **H.E.S.S. Collaboration (2023), *Discovery of a Radiation Component from the Vela Pulsar Reaching 20 Teraelectronvolts*.** [Paper](https://arxiv.org/abs/2310.06181) · [journal](https://doi.org/10.1038/s41550-023-02052-3). Establishes an additional very-high-energy component that challenges a single-process explanation of Vela. Begin validation with a declared Fermi energy band; reproducing a GeV pulse profile would not establish that the model also explains the TeV radiation.
- **Mignani, Zharikov & Caraveo (2007), *The Optical Spectrum of the Vela Pulsar*.** [Paper](https://arxiv.org/abs/0707.2036). Optical spectroscopy supports a magnetospheric continuum. A thermal glowing stellar surface would not explain that observed optical component.
- **Ng & Romani (2008), *Fitting Pulsar Wind Tori. II. Error Analysis and Applications*.** [Paper](https://arxiv.org/abs/0710.4168). Nebular geometry constrains the spin-axis viewing angle: approximately 61–63° for Crab, 63.6° for Vela in their fits. This is not the magnetic obliquity. Modeling/systematic uncertainties must accompany any adopted value.

**Surface gravity and rapid rotation**

- **Beloborodov (2002), *Gravitational Bending of Light Near Compact Objects*.** [Paper](https://arxiv.org/abs/astro-ph/0201117). Useful analytic light-bending check; it is an approximation, not an exact replacement for geodesic integration.
- **Poutanen & Beloborodov (2006), *Pulse profiles of millisecond pulsars and their Fourier amplitudes*.** [Paper](https://arxiv.org/abs/astro-ph/0608663). Independent checks for Doppler effects, angular emission laws, light bending, and time delays.
- **Morsink et al. (2007), *The Oblate Schwarzschild Approximation for Light Curves of Rapidly Rotating Neutron Stars*.** [Paper](https://arxiv.org/abs/astro-ph/0703123). Includes the rotating star's flattened surface and changed surface normals with Schwarzschild ray propagation and rotational Doppler effects. Visually squashing a sphere is insufficient.
- **Nättilä & Pihajoki (2018), *Radiation from rapidly rotating oblate neutron stars*.** [Paper](https://arxiv.org/abs/1709.07292). A more complete treatment of the rotating neutron-star exterior, including frame dragging and quadrupole effects. The Kerr–Newman black-hole metric from NO.001 is not an interchangeable neutron-star metric.
- **Kini et al. (2026), *A NICER View of PSR J0030+0451: Updated Constraints from Six Years of NICER Observations*.** [Manuscript](https://arxiv.org/abs/2602.23743). A current surface-model reference: $M=1.43^{+0.20}_{-0.17}M_\odot$ and $R_{\rm eq}=12.68^{+1.31}_{-1.04}$ km, with 68% credible intervals and temperature-gradient evidence. These are inferred, model-dependent quantities. The accessed arXiv record is marked submitted to ApJ.

**Recent preprints to evaluate, rather than silently treat as settled physics**

- **Pétri (31 July 2026), *Polar cap plasma loading and the morphology of pulsar gamma-ray light curves*.** [Preprint](https://arxiv.org/abs/2607.29126). An analytically specified current sheet with asymmetric emissivity can reproduce varied pulse shapes. Promising for a tractable original model; the emissivity and loading are prescribed. The accessed record says submitted to A&A.
- **Nättilä & Salmi (22 August 2026), *Pair Discharges and Radio Emission from Pulsar Magnetospheres*.** [Preprint](https://arxiv.org/abs/2608.21920). One-dimensional QED/PIC discharge calculations address radio microphysics. They do not supply a global three-dimensional magnetosphere or a universal radio beam.

## What the animation should actually calculate

Use consistent units, distinguish the magnetic obliquity $\chi$ from the observer angle $\zeta$, and keep physical time separate from screen playback time. The following is an implementation specification synthesized from the cited work, not a simulation already performed.

**Rotation and scale.** Set $\Omega=2\pi/P$ and $R_{\rm LC}=c/\Omega$. The light cylinder is a cylindrical corotation limit around the spin axis, not a material shell. Plasma beyond it cannot rigidly corotate. Representative periods of 33.65 ms (Crab), 89 ms (Vela), and 4.87 ms (J0030) give light-cylinder radii of roughly 1,606, 4,246, and 232 km. These radii are calculated from the rounded periods, not separate measurements; they do not depend on an assumed stellar radius. A 12 km stellar radius is a separate illustrative choice for comparing these scales. [Crab timing](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/3PC_HTML/J0534%2B2200.html) · [J0030 timing](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/3PC_HTML/J0030%2B0451.html) · [Vela reference](https://arxiv.org/abs/2310.06181)

At young-pulsar wind scales the star will occupy very few pixels. Use a separately labeled close-up, or disclose any enlargement. Do not set the star to a large fraction of the light cylinder merely to make the composition attractive.

**Electromagnetic field.** In Gaussian units, a force-free solver should evolve Maxwell's equations and enforce $\mathbf E\cdot\mathbf B=0$, $\rho_e\mathbf E+\mathbf J\times\mathbf B/c=0$, and magnetic domination away from dissipative sheets. The spin-down test is approximately $L=\mu^2\Omega^4(1+\sin^2\chi)/c^3$. Define the magnetic-moment convention explicitly: for a centered dipole with polar field $B_p$, $\mu=B_pR_\star^3/2$. Do not confuse polar and equatorial field strengths. [Spitkovsky 2006](https://arxiv.org/abs/astro-ph/0603147)

**Outer current sheet.** A useful ultrarelativistic split-monopole limit is

$$\cos\theta\cos\chi+\sin\theta\sin\chi\cos[\phi-\Omega(t-r/c)]=0.$$

This provides a traveling corrugated sheet, not its thickness, particle distribution, heating, or brightness. It is an outer-wind approximation and a validation case; prescribing it throughout the star's surroundings would not reproduce the inner closed magnetic region. [Bogovalov 1999](https://arxiv.org/abs/astro-ph/9907051)

**Photon arrival time.** For straight propagation toward a distant observer in direction $\mathbf n$, remove the common distance term and use

$$t_{\rm obs}=t_{\rm em}-\mathbf n\cdot\mathbf r/c.$$

Accumulate photons by arrival phase, with directions and radiation calculated in a consistent frame. Relativistic emission is concentrated around particle velocity, so the entire current sheet will not appear equally luminous to every observer. Surface photons require the appropriate GR bending and delay treatment. [Cerutti et al. 2016](https://arxiv.org/abs/1511.01785) · [Bogdanov et al. 2019](https://arxiv.org/abs/1912.05707)

**Emission.** Select a band before assigning color or brightness. For the kinetic route, use particle momentum and electromagnetic fields to compute the selected radiation process and its spectrum. For a reduced model, declare the assumed particle-energy distribution, pitch angles, sheet width, spatial emissivity, and emission radius range. Merely painting a solved magnetic field does not calculate its radiation. Treat competing curvature and synchrotron prescriptions as model alternatives rather than blend them invisibly. [Kalapotharakos et al. 2023](https://arxiv.org/abs/2303.04054)

**Surface view, if included.** Specify mass, radius, rotation, atmosphere, and temperature map. Trace surface photons, include occultation and the angle-dependent atmospheric intensity, and verify Doppler/redshift transformations. Start with controlled spherical or oblate-Schwarzschild benchmarks; use a neutron-star metric when the chosen rotation and required accuracy warrant it. [Morsink et al. 2007](https://arxiv.org/abs/astro-ph/0703123) · [Nättilä & Pihajoki 2018](https://arxiv.org/abs/1709.07292)

## Available observations and software

**Fermi 3PC:** the [official portal](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/) provides FITS/ASCII phase profiles, uncertainties, energy bands, spectra, and ephemerides. Use the July 2025 corrected profile archive. Read FITS extensions by name, since radio profiles and fitted profiles are not available for every object. Compare peak separation, widths, bridge emission, energy dependence, and radio/gamma phase offsets. These rotational profiles are different from the Fermi repository's multi-day variability curves. [Profile format](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/README_3PC_ProfileData.txt)

**Entity:** [public code](https://github.com/entity-toolkit/entity), BSD-3-Clause, with CPU/GPU-portable PIC infrastructure. Its [special-relativistic method paper](https://arxiv.org/abs/2511.17710) and [general-relativistic method paper](https://arxiv.org/abs/2511.17701) were posted in November 2025. A serious candidate for future simulation work; the existence of this framework does not establish that a selected pulsar setup runs locally or has reproduced the reference figures.

**Aperture4:** [public code](https://github.com/fizban007/Aperture4) and [setup documentation](https://www.pulsarmagnetosphere.com/Aperture4_docs/1-setup/). CPU/GPU PIC framework with a GPL-3.0 license. Another candidate requiring an explicit stellar boundary, plasma injection, radiation model, and numerical convergence study.

**Zeltron:** [official download page](https://ipag.osug.fr/~ceruttbe/Zeltron/download.html) and [repository](https://github.com/bcerutti/Zeltron-code). The download page describes public Cartesian 2D/3D v1.0 and GPLv3-or-later licensing. The exact spherical pulsar setup and full particle/field dumps underlying the 2025 study were not located as a verified public download in this search. Do not claim the public release reproduces that paper unchanged.

**Surface checks:** [X-PSI](https://xpsi-group.github.io/xpsi/) and the [X-PSI / Illinois–Maryland / Alberta comparison repository](https://github.com/DevarshiChoudhury/Neutron-star-relativistic-ray-tracing-code-comparisons) supply a much firmer verification route for a surface-focused calculation than eyeballing a glowing sphere.

No codes were installed or simulations run during this research. Availability was checked from the linked primary repositories and documentation. No complete reference three-dimensional PIC dataset has yet been secured.

## Target choice

| Direction | Why choose it | What we must not imply |
| --- | --- | --- |
| **Vela-like gamma-ray pulsar — recommended for NO.002** | Strong observed pulses and a direct connection to the modern current-sheet models | A unique reconstructed magnetic obliquity, stellar mass/radius, or visible-color image |
| **Crab optical pulsar** | Actual visible/near-IR spectrum, pulse information, and polarization comparisons | The measured integrated spectrum resolves the local color of every plasma element |
| **PSR J0030+0451 surface X-rays** | Surface pulse modeling, inferred hot regions, GR tests, and current NICER constraints | Hotspots are a directly photographed map, or are necessarily two identical antipodal circles |

The displayed period must refer to a named timing model/epoch. Use rounded values for educational labels and the selected ephemeris for any quantitative phase comparison. Surface mass/radius, magnetic tilt, field topology, pair supply, and observer geometry have different evidential status; record each as measured, inferred, or chosen.

## Production and validation plan

1. **Lock the observable.** For the recommended path, use a high-energy emission model. Choose either a directional radiation image or a labeled volume rendering of plasma/emitted power, accompanied by a fixed-observer pulse trace.
2. **Establish reproducible physical input.** Prefer a documented 3D PIC dataset or a reproducible research-code run. If that proves impractical, solve a force-free field and declare the emission prescription. A wholly analytic current sheet is a lower-fidelity wind model and must carry that name.
3. **Test before styling.** Check electromagnetic divergence and energy accounting, force-free limits, boundary effects, spatial/temporal refinement, and the spin-down benchmark. Verify beaming and arrival-time signs using analytic cases. If surface emission is included, reproduce published pulse-profile benchmarks with their exact inputs.
4. **Compare with observations.** Fit or compare a declared energy band and fixed observer angle. Keep observational uncertainty and model residuals visible in the research output. An attractive double peak alone is not a quantitative validation.
5. **Make the presentation faithful.** Use false color for radio/X-ray/gamma-ray intensity. Any field lines are explicitly diagnostic overlays. A moving educational camera must not silently change the fixed observer used for the light curve. The light cylinder is an optional diagram overlay, not luminous matter.
6. **Handle time and looping honestly.** Display a slowdown factor. For a periodically prescribed model, repeat an integer number of rotations. A noisy PIC flow is not exactly periodic; a phase-averaged cycle or edited repeat must be labeled, rather than claimed to be an exact dynamical loop.
7. **Reuse the established visual identity.** Retain @sideratlas, NO.002, THE PULSAR, black background, scientific fonts, readable Instagram margins, and a companion caption containing every study actually used. Keep the full bibliography and assumption table with the code.

The next useful milestone is a small validated field-and-pulse calculation with a fixed observer. That establishes whether the physical model supports the intended image before spending time on a polished 32-second export.
