# Optical observer rendering: requirements and claim boundaries

Research checked 2026-10-17 against primary papers and NASA's own descriptions. This document is a design/audit, not evidence that these calculations have been implemented.

## The reference is not an optical camera image

NASA SVS [13058](https://svs.gsfc.nasa.gov/13058/) and the underlying [4648 all-particle visualization](https://svs.gsfc.nasa.gov/4648/) display simulation macroparticle tracks. Blue/red distinguish electrons/positrons, with whitening indicating higher Lorentz factors. Those colors are particle diagnostics, not visible photons. NASA's supplied paths through the simulation also include camera tours and a transition from inertial to corotating views. They cannot establish what a fixed optical observer would receive.

NASA links these animations to Brambilla et al. (2018), [Electron–Positron Pair Flow and Current Composition in the Pulsar Magnetosphere](https://arxiv.org/abs/1710.03536), DOI [10.3847/1538-4357/aab3e1](https://doi.org/10.3847/1538-4357/aab3e1). NASA's [dataset page/4648 description](https://svs.gsfc.nasa.gov/4648/) says the site does not store the underlying datasets. Its absence there is not proof that no archive exists elsewhere.

## What the NASA-linked calculation actually assumes

Brambilla et al., section 2, uses C-3PA, a relativistic electromagnetic PIC code with particle charge/current feeding Maxwell's equations and the updated fields pushing both species. Important numerical scales are:

- Effective conducting star radius 0.36 light-cylinder radii, with conductor below 0.28 and a transition layer
- Cartesian cube side 9.6 light-cylinder radii, spatial step 0.02, temporal step 0.003 light-cylinder light-crossing times
- Surface pair supply controlled by a magnetization condition and prescribed global injection rate; pairs start at rest
- Voltage Lorentz-factor scale about 500, unlike about a billion for a young physical pulsar
- Enhanced radiation reaction to cool perpendicular momentum; incompletely resolved inner skin depth and acknowledged numerical heating

These choices resolve a global kinetic problem at reduced physical scale. They are not a physical Vela atmosphere or an absolute optical radiation prediction. They must be disclosed in any attempted reproduction. They cannot be replaced by an inexpensive coarse run while claiming the NASA calculation has been matched quantitatively.

The companion paper, Kalapotharakos et al. (2018), [3D Kinetic Pulsar Magnetosphere Models: Connecting to Gamma-Ray Observations](https://arxiv.org/abs/1710.03170), section 5, explains why simply multiplying PIC energies is invalid. Real radiative cooling and gyro scales cannot simultaneously be represented by the reduced run. Its physical gamma-ray treatment recalculates particle energies in the model fields in an asymptotic drift approximation and assumes curvature emission. This procedure does not supply physical optical synchrotron emission or a pair pitch-angle distribution.

## Vela optical data require different physics from the old GeV map

Mignani, Zharikov & Caraveo (2007), [The Optical Spectrum of the Vela Pulsar](https://arxiv.org/abs/0707.2036), DOI [10.1051/0004-6361:20077774](https://doi.org/10.1051/0004-6361:20077774), provide a useful optical check: the extinction-corrected 4000–8000 Å continuum is approximately flat in frequency flux, with their convention `F_nu ∝ nu^(-alpha)` and `alpha = -0.04 ± 0.04`. The measured visual magnitude is about 23.6. Their result identifies this optical component as magnetospheric. This is a check on the integrated spectrum; it does not determine the spatial image uniquely.

Harding, Kalapotharakos, Barnard & Venter (2018), [Multi-TeV Emission From the Vela Pulsar](https://arxiv.org/abs/1811.11157), DOI [10.3847/2041-8213/aaf3b2](https://doi.org/10.3847/2041-8213/aaf3b2), explicitly models the IR/optical component as synchrotron from cascade pairs whose pitch angles are maintained by resonant radio absorption. This is a model rather than a unique established emission mechanism, but it identifies essential missing inputs: pair energy distribution, multiplicity, pitch/gyro phase, radio intensity/direction, and radiation evolution. A zero-pitch, primary-only GeV curvature calculation filtered to optical wavelengths does not reconstruct that population.

Pavlov et al. (2001), [The X-ray Spectrum of the Vela Pulsar Resolved with Chandra](https://arxiv.org/abs/astro-ph/0103171), DOI [10.1086/320342](https://doi.org/10.1086/320342), models the thermal component with a magnetized hydrogen atmosphere. Its thermal fit must not be interpreted as a measured spatial surface texture. A blackbody with the same temperature is not an implemented magnetized-atmosphere transfer calculation. This surface component and nonthermal optical emission are separate physical components.

## Minimal observer equations

Use physical positions, momenta and coordinate times from a consistent physical model, rather than post hoc colors assigned to track energy. For a fixed far observer with source-to-observer unit vector `n`, source time `t`, and position `r(t)`, subtract the irrelevant constant distance delay:

`tau = t - n · r(t) / c`

`d tau / dt = 1 - n · beta`

For a finite-distance camera use the exact distance `tau = t + |x_cam-r(t)|/c` and the event-dependent direction. Interpolate emission into **arrival-time** bins. Do not rotate a frozen simultaneous-volume emissivity and call the result a retarded image.

Classical Liénard–Wiechert radiation in flat spacetime gives, in SI units:

`E_rad = [q / (4 pi epsilon_0 c R)] [n × ((n-beta) × beta_dot)] / (1-n·beta)^3`

The flux per observer-time and solid angle therefore has the sixth power in the denominator. The **emitted energy per source-time** and solid angle has the fifth power because of the arrival-time Jacobian. Specifically, with `A=n×((n-beta)×beta_dot)`:

`dE / (dt_source dOmega) = q^2 |A|^2 / [16 pi^2 epsilon_0 c (1-n·beta)^5]`

Deposit this quantity times the source-time interval in the appropriate arrival bin; divide by the arrival-bin duration at the end. Avoid applying the Jacobian twice.

For a spectral camera the far-field positive-angular-frequency spectrum from one trajectory can be written:

`d^2E/(domega dOmega) = q^2/(16 pi^3 epsilon_0 c) |integral [A/(1-n·beta)^2] exp{i omega[t-n·r(t)/c]} dt|^2`

This is a frequency-domain trajectory integral, not an instantaneous optical power multiplier. It assumes classical point-charge radiation in flat spacetime. Window endpoints need inertial continuation or a validated equivalent treatment; abruptly creating/deleting tracks produces spurious radiation. The full frequency integral also requires sufficient trajectory accuracy and oscillatory quadrature.

An implementable primary reference for a trajectory-to-detector solver is Watt & O'Shea (2024), [A differentiable simulation package for performing inference of synchrotron-radiation-based diagnostics](https://journals.iucr.org/s/issues/2024/02/00/tv5053/), DOI [10.1107/S1600577524000663](https://doi.org/10.1107/S1600577524000663). Their SYRIPY implementation separates particle tracking, retarded radiation and detector propagation, handles artificial endpoint radiation, and uses specialized oscillatory quadrature. Their accelerator benchmarks do not validate a pulsar model.

For incoherent optical radiation, sum **single-physical-particle intensities** with multiplicity weights: `I_total = sum w_i I_i`. Do not use a macroparticle charge `w_i e` in a squared radiation amplitude and thereby produce an artificial `w_i^2` signal. Coherent radio emission is a different calculation and cannot be inferred from this incoherent optical treatment.

## Synchrotron approximation: what can and cannot be implemented quickly

A local ultrarelativistic synchrotron calculation is practical only after a physical pair energy/pitch distribution is supplied. It uses

`omega_c = 3 gamma^3 c / (2 rho)`

`F(x) = x integral_x^infinity K_(5/3)(u) du`

where `rho` is the actual trajectory curvature radius, including gyro motion if it is the radiation source. The observer-angle spectrum requires the two polarization terms containing `K_(2/3)^2(xi)` and `K_(1/3)^2(xi)`, with `xi=(omega rho/3c)(gamma^-2+theta^2)^(3/2)`. It is not generally separable into an angle-integrated spectrum times a bolometric beaming factor.

Lai et al. (2025), [Spectropolarimetry of synchrotron radiation from relativistic electrons with anisotropic pitch-angle and various energy distributions](https://doi.org/10.1093/mnras/staf1295), equations 2–9, give the angular kernels and show why an isotropic momentum assumption can fail. In their helical-orbit formula the per-pulse spectral energy is multiplied by the gyro frequency to obtain time-averaged power. Do not reuse the per-pulse formula at every trajectory step: that double counts the formation interval. A tangent-following beaming Gaussian of arbitrary width is not this calculation.

The existing `nonthermal.py` is an angle-integrated synchrotron building block, as documented in `research/RADIATION-UPGRADE.md`. It is neither a physical particle population nor an observer solver. Its applicability checks, particularly `gamma*sin(pitch)` and quantum parameter, must remain active. In a frame with strong `E×B` drift, transform into the local magnetic/drift frame before applying a magnetic-gyration kernel; transform photons and intensity back consistently. Do not identify laboratory `B` and laboratory pitch with the magnetic-frame values unconditionally.

Do not add a separate synchrotron loss to curvature radiation calculated from fully resolved total acceleration: that counts the same radiation twice. A guiding-centre plus unresolved gyro treatment can split them only when the separation and energy/pitch evolution are defined consistently.

## The scale constraint in a full-frame Reel

For period about 0.08937 s, `R_LC=cP/(2pi) ≈ 4264 km`. An assumed 13 km star has `R_star/R_LC≈0.00305`. A global view spanning four light-cylinder radii across 1080 pixels gives a stellar radius about **0.82 pixel**, before the modest general-relativistic increase of its apparent radius. A conspicuous 50–150 pixel star beside a globally visible magnetosphere would be a scale illustration.

The correct way to fill the Reel is a clearly specified camera field of view and crop of the simulated domain, retaining one common physical scale. Rotating the camera, changing exposure, showing slow motion, adding a finite point-spread function, or mapping photons into RGB are presentation assumptions that must be recorded. They do not license making particles into luminous blue/red wires.

For a direct optical view, field lines are not rendered as emitters. Gamma rays cannot become visible purple cones. Actual photon colors are species independent. Particle trails may be a separately labeled diagnostic product, but must not be composited into the optical-camera product while claiming they are what the eye sees.

Near the stellar surface a physical atmosphere image requires gravitational redshift, ray bending, occultation and angle-dependent atmosphere intensity. A flat-spacetime outer PIC volume does not include those effects. If a surface atmosphere and outer-volume image use different numerical radii or units, their composite is not a single physical camera image without an explicit matching construction.

## Required validation before claiming a Vela optical reconstruction

1. Lorentz pusher: analytic uniform-field helices, constant electric acceleration, time-step convergence and kinetic/field energy accounting
2. PIC: charge conservation, Gauss residuals, plasma/skin/gyro resolution, boundary reflections, grid/time/particle-number and domain convergence
3. Radiation: constant velocity produces no radiation including endpoints; uniform circular motion recovers the analytic synchrotron spectrum and angular pattern; angular/energy integration recovers the Liénard power
4. Retardation: demonstrate that source-time deposition and direct observer-time flux agree and conserve the emitted energy
5. Macroparticles: splitting one particle into equal-weight copies leaves incoherent flux unchanged
6. Camera: convergence of temporal exposure, frequency quadrature, angular integration, pixels and particle sampling; do not broaden physical beams to hide undersampling
7. Observation: extinction-corrected optical continuum and normalization, optical pulse structure/phases, and atmosphere/X-ray component compared separately; record which data tuned the model versus independently assessed it
8. Error budget: scaled-PIC versus physical-Vela extrapolation must be bounded or explicitly left unvalidated, not hidden behind successful unit tests

## Caption wording boundary

Use the requested NO.002 format, but only include statements supported by the implemented result. A genuine completed physical-camera product could list under CALCULATED the Lorentz trajectories, electromagnetic feedback, physical radiation toward the camera, optical wavelength integration and arrival delays. Each item requires evidence, not just presence of an equation in the repository.

Under ASSUMED record injection/pair supply, composition, momentum distribution, field strength, mass/radius/temperature, viewing geometry, magnetic inclination, outer boundary, atmosphere model/grid, reduced-scale extrapolation, exposure/tone mapping, and any constructed periodic loop. State the actual simulated time; an eight-second visual cycle or repeated clip is not eight seconds of physical pulsar evolution.

Until these gates are met, defensible wording is: **“Reduced-scale kinetic pulsar simulation with a synthetic optical-radiation calculation; Vela's period/geometry used as reference, optical spectrum and pulse profile not validated.”** If even the observer solver is not implemented, remove the synthetic-optical claim. Do not write “NASA-matched Vela reconstruction,” “true optical appearance,” or “Vela atmosphere simulated” for a tracer, prescribed-cone, blackbody or scaled-PIC diagnostic render.
