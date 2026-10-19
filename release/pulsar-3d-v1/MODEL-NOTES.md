# SiderAtlas NO.002 — 3D pulsar edition

This edition explains a pulsar's three-dimensional structure: a prominent central neutron star, an oblique rotating magnetic field and paired radiation-direction envelopes. It is a scientifically motivated conceptual view of a plausible young rotation-powered pulsar. Vela supplies reference scales, not an observational fitting target. The earlier optical observer calculation remains available as separate quantitative research.

## Retained physical inputs

| Quantity | Adopted value |
|---|---|
| Rotation period | 0.08937108859 s |
| Mass / circumferential radius | 1.4 M☉ / 13 km |
| Polar surface field | 3 × 10¹² G |
| Redshifted effective temperature | T∞ = 0.70 MK; local T_eff ≈ 0.848 MK |
| Magnetic inclination to rotation axis | 75° |
| Camera inclination to rotation axis | 65° |
| Light-cylinder radius | R_LC = cP/(2π) ≈ 4264 km |

The field is strong and the star rotates rapidly, but these inputs are representative of an energetic young pulsar rather than a population average. Their existing energy, radiation-transport and scale audits are preserved. This rendering revision adds no particle evolution, atmosphere solve, Vela fitting or expensive convergence study.

## What the visible layers mean

**Stellar surface — calculated approximate radiation.** The star reuses the existing normalized optical surface image from the two-mode LTE hydrogen atmosphere. The surface camera integrates outward Schwarzschild light rays, gravitational redshift and local special-relativistic Doppler/aberration. The image is uniform and axisymmetric; no invented rotating hot spots are added. Its peak linear RGB channel is normalized to 4.5 before adding diagram layers and applying a fixed per-channel display shoulder, x/(1+x), followed by sRGB encoding. The shoulder compresses highly exposed colors toward white; the displayed color is therefore tone-mapped rather than a direct spectrum. This enhanced exposure is a display choice, not a naked-eye brightness prediction. The homogeneous normal-field, fully ionized hydrogen atmosphere uses simplified magnetic free–free and scattering opacities; atomic opacity, dispersion, vacuum conversion and return-current irradiation remain omitted. Stored scientific-validation flags remain false. See the [optical model notes](../optical-v1/MODEL-NOTES.md) for the complete calculation and limitations.

**Magnetic curves — calculated geometry, illustrative visibility.** In a frame aligned with the magnetic unit vector m, the prescribed flat-space dipole is

`B = (B_p/2)(R/r)^3 [3(m·r_hat)r_hat − m]`

Its integral curves satisfy `r = L sin²θ`, with θ measured from the magnetic axis. Representative closed families have L/R = 2.0, 3.2, 4.7 and 5.6. A few polar segments indicate the local direction of flux that a global rotating magnetosphere could carry outward. A segment leaving the close view is not proof that a static dipole line is open at infinity. The chosen curve count and luminous widths do not encode magnetic flux density or particle density. Magnetic field lines are mathematical guides and are not intrinsically glowing wires.

The star retains its 13 km geometric radius relative to these curves; it is not enlarged independently. The close view shows the near zone: the largest closed guide reaches 5.6R ≈ 0.0171 R_LC, and local polar guides extend to 11R ≈ 0.0335 R_LC. The polar guides fade between 7R and 11R. Rigid rotation is therefore a useful geometric approximation here. It is not a force-free or PIC solution, and omits current-driven distortion, sweepback, current sheets and wind structure. The flat-space dipole also omits general-relativistic field distortion near the compact star. Field and beam guides use Euclidean projection; only the stellar radiation uses the relativistic surface camera. Field curves fade in between 1.3R and 1.7R. Their near-surface fade avoids suggesting that the Euclidean guides were ray-traced through the atmosphere's lensed image.

**Polar envelopes — illustrative radiation directions.** Two smooth blue-white envelopes co-rotate with opposite magnetic poles. Their angular scale is motivated by the near-zone dipole tangent at an assumed reference radius `r_em = 6R = 78 km` from the center (65 km above the surface). Using the approximate open-cap relation `sin²θ = r_em/R_LC`, the tangent angle to the magnetic axis is

`ρ = θ + atan(tanθ/2) = atan2(3 sinθ cosθ, 3 cos²θ − 1) ≈ 11.7°`

This sets a nominal angular width scale, not a hard beam edge. The volume renderer uses a transverse Gaussian, `exp[−d_perp²/(2σ²)]`, with `σ/R = 0.055 + tanρ max(|s|/R − 1, 0)`, where s is axial distance from the center. A narrow core has 0.3 times this width. A smooth axial window starts near the star and fades between 10R and 14R; additional axial brightness decay is a display choice. The tangent angle therefore sets an approximate Gaussian σ scale, not a 1/e half-width or a measured pulse width.

The open-cap relation is a geometric estimate, not an oblique force-free polar-cap calculation. The selected reference radius, smooth profile, exposure, color and finite length are illustrative. The envelopes show how a rotating directional radiation pattern can produce a lighthouse effect; they are not individual photon trajectories, luminous material cones, calculated particle paths, or a prediction of radio/optical/gamma-ray emission at 6R. Real emission sites and beam patterns depend on band and plasma physics. The visual analogy does not imply that optical radiation in our retained particle model is emitted in these envelopes.

**Quantitative plasma radiation — preserved separately.** No magnetospheric photon-event archive is displayed in this 3D edition. Existing guiding-centre paths, cooling, directional helical harmonics, Lorentz transformations and photon arrival-time calculations remain unchanged in the research package. The local native PIC experiments are also separate. Their numerical particles are not rescaled into physical optical photons, and their enlarged numerical star is not used for the geometry of this movie. Guide light and surface radiation do not share an absolute radiance calibration.

## Camera and motion

The orthographic camera looks toward the center with a −32° roll, an inclination of 65° and a smooth ±12° azimuth wobble. A geometric scale of 105 native pixels per stellar radius at 1080 × 1920 makes the star prominent while magnetic curves and beam envelopes span the portrait image. The true geometric radius and lensed optical apparent radius are distinct: `g = sqrt[1−2GM/(Rc²)] ≈ 0.8258`, so `R_apparent = R/g ≈ 1.211R`. Opaque stellar masking uses that apparent outline; the surrounding field and beam guides use Euclidean depth. This hybrid is a teaching composition, not a common general-relativistic image of all layers. A small screen-space stellar halo is a presentation effect. The guide brightness is independent of the stellar exposure; the combined layers do not have an absolute radiometric calibration.

The camera motion reveals depth; it is not an observer-timing experiment. An 8-second constructed loop represents one physical 89.37 ms rotation and is repeated four times for a 32-second Reel. The source rotates, rather than evolving through 32 physical seconds. Only `@sideratlas`, `NO.002` and `THE PULSAR` are overlaid.

## Sources and roles

- Spitkovsky 2006, *Time-dependent force-free pulsar magnetospheres: axisymmetric and oblique rotators* — physical context and preserved spin-down scale; no force-free solution is recreated by the guide renderer: https://arxiv.org/abs/astro-ph/0603147, DOI 10.1086/507518
- Brambilla et al. 2018, *Electron–Positron Pair Flow and Current Composition in the Pulsar Magnetosphere* — PIC context for NASA's reference visualization, not the source of our illustrated curves or beams: DOI 10.3847/1538-4357/aab3e1
- Smith et al. 2023, third Fermi-LAT pulsar catalog — period/reference context: https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/
- Harding & Kalapotharakos 2015, *Synchrotron Self-Compton Emission from the Crab and Other Pulsars* — preserved pair-emission research context: https://arxiv.org/abs/1508.06251, DOI 10.1088/0004-637X/811/1/63
- Harding et al. 2018, *Multi-TeV Emission From the Vela Pulsar* — preserved pair optical-emission context; its fitted broadband model is not reproduced: https://arxiv.org/abs/1811.11157, DOI 10.3847/2041-8213/aaf3b2
- Marcowith & Malzac 2003 — preserved classical helical-harmonic emission: https://arxiv.org/abs/astro-ph/0309107, DOI 10.1051/0004-6361:20031060
- van Adelsberg & Lai 2006 — magnetic atmosphere transfer basis, not their complete solver: https://arxiv.org/abs/astro-ph/0607168
- Ofengeim & Zyuzin 2018 — retained stellar thermal/mass/radius reference: https://www.ioffe.ru/astro/Stars/Paper/ofengeim_zyuzin18.pdf
- Poutanen & Beloborodov 2006 — relativistic stellar radiation: DOI 10.1111/j.1365-2966.2006.11088.x
- Poutanen 2020 — Schwarzschild surface light bending: DOI 10.1051/0004-6361/202037471
- NASA / Goddard 2018, *Simulations Create New Insights Into Pulsars* — visual/context reference. Its blue/red tracks are PIC macroparticle diagnostics, not visible optical photons: https://svs.gsfc.nasa.gov/13058/
- NASA / Goddard 2008, *NASA's Fermi Telescope Discovers First Gamma-Ray-Only Pulsar* — conceptual lighthouse reference: https://www.nasa.gov/universe/nasas-fermi-telescope-discovers-first-gamma-ray-only-pulsar/

No NASA footage, simulation dataset or code is used. All prior optical source roles, spectral-color attribution and research checks remain in the linked optical model notes. The 3D guide palette is illustrative and is not derived from CIE spectral integration.
