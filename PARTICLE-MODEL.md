# NO.002: a Vela-tested particle and emission hypothesis

This stage replaces a prescribed power-law emissivity with particle-energy transport and a calculated curvature spectrum. It does **not** solve the plasma magnetosphere. The existing split-monopole benchmark remains intact; its Maxwell/force-free checks do not validate this different magnetic geometry.

## Geometry and actual particle paths

We use the retarded magnetic **point dipole in the laboratory frame**, with rotating moment m and retarded time t−r/c:

\[
\mathbf B={3\hat n(\hat n\cdot\mathbf m)-\mathbf m\over r^3}
+{3\hat n(\hat n\cdot\dot{\mathbf m})-\dot{\mathbf m}\over cr^2}
+{\hat n(\hat n\cdot\ddot{\mathbf m})-\ddot{\mathbf m}\over c^2r}.
\]

The transverse electric prescription is corotation, not the vacuum dipole's electric field. In units c=Ω=R_LC=1:

\[
\boldsymbol\beta=\boldsymbol\beta_{\rm rot}+f\hat b,\qquad
f=-\boldsymbol\beta_{\rm rot}\cdot\hat b
\pm\sqrt{(\boldsymbol\beta_{\rm rot}\cdot\hat b)^2+1-\beta_{\rm rot}^2},
\qquad \boldsymbol\beta_{\rm rot}=\hat z\times\mathbf q.
\]

Choose the outward root separately at each magnetic pole, retaining that root along the path. Integrate dq/dτ=f b̂ in rotating coordinates. These conventions follow [Bai & Spitkovsky (2010), vacuum-dipole analysis](https://arxiv.org/abs/0910.5740). Their analysis also documents sensitivity to polar-cap geometry and the vacuum approximation near the light cylinder. The implementation traces the cap numerically rather than substituting a circular rim.

The laboratory trajectory is x(τ)=R_z(τ)q(τ). Its curvature is calculated from

\[
{R_{\rm LC}\over R_c}=\left|(f\hat b\cdot\nabla)\boldsymbol\beta+
\hat z\times\boldsymbol\beta\right|.
\]

This includes rotation; using magnetic-field-line curvature alone would give a different answer. [Bai & Spitkovsky's force-free companion, Appendix C](https://arxiv.org/abs/0910.5741), supplies the corresponding trajectory-curvature expression and distinguishes particle-path from field-line weighting. The present code does **not** import that paper's force-free field solution.

Seeds occupy a thin annulus of geometrically open tubes at **0.08 R_LC**, weighted by unsigned magnetic flux. This surface is about 340 km from the centre for the rounded Vela period; it is **not the neutron-star surface**. The nominal annulus spans 0.95–0.995 times the traced rim angle. Both hemispheres have equal flux-normalized injection; no phase-dependent source brightness is fitted. Particles begin at γ=1000 and stop at cylindrical radius 0.95 R_LC or spherical radius 1.2 R_LC. Terminal steps are clipped to the escape surface. Upstream acceleration and pair injection are not calculated.

Geometrically “open” here means crossing the light cylinder: vacuum dipole lines can close farther out. The rim finder assumes a single boundary representable as polar angle versus magnetic azimuth. Its bracket and actual cylinder crossings are checked, with independent step-convergence tests. It is not a general solver for arbitrary cap topology.

## Acceleration, energy losses, and photon spectrum

For each laboratory arclength segment:

\[
{d\gamma\over ds}={eE_{\rm work}\over m_ec^2}
-{2e^2\gamma^4\over3m_ec^2R_c^2},\qquad
P_c={2e^2c\gamma^4\over3R_c^2},\qquad
E_c={3\hbar c\gamma^3\over2R_c}.
\]

The fitted acceleration scale is A_parallel=e|E_parallel|/(m_e c²), in cm⁻¹. Because the particle velocity includes drift, electric work uses

\[
E_{\rm work}=|E_\parallel|\,|\boldsymbol\beta\cdot\hat b|.
\]

Its sign is chosen to accelerate the outward population at each pole. This is prescribed charge-aligned work, not a computed gap potential or screening solution. Direction feedback from E_parallel is neglected. The guide directions use |v|=c; γ≥1000 makes the initial speed approximation below 5.1×10⁻⁷ relative. Strong local fields still require the guiding-centre approximation; this project has not solved particle gyro-motion or pair distributions.

[Barnard et al. (2022)](https://arxiv.org/abs/2111.03405) motivates this curvature-energy transport and energy-dependent pulse analysis, but uses numerical force-free trajectories rather than the geometry here. Agreement with its equations is not replication of its Vela model. Acceleration is held constant along B in the current model family, and fitted on a declared finite grid; no posterior uncertainties or unique parameter inference are claimed.

Transport uses the exact positive solution for frozen acceleration/curvature on each segment. Work minus particle-energy change gives the emitted energy. Spectral weights use

\[
F(x)=x\int_x^\infty K_{5/3}(u)\,du,\qquad
{dP\over dE}={P_c\over E_c}\,{F(E/E_c)\over I},\qquad
I={8\pi\over9\sqrt3}.
\]

Photon counts integrate (dP/dE)/E over each energy interval. All these quantities are laboratory-frame quantities; applying another Doppler factor would double-transform them. Tests compare the numerical kernel with independent Bessel quadrature, verify total power and photon integrals, compare energy transport to analytic cooling/equilibrium and a separate ODE solver, and refine the spectral integration during rapid cooling.

## Photons reaching the observer

For a photon direction β at rotating-coordinate position q, its sky colatitude is arccos β_z and its arrival phase is

\[
\phi_{\rm obs}={-\operatorname{atan2}(\beta_y,\beta_x)-\mathbf q\cdot\boldsymbol\beta\over2\pi}\pmod1.
\]

Both aberration and propagation delay enter this mapping. Steady injection times per-particle emitted photon yield supplies residence-time weighting automatically; field-line length is not substituted for flight time. Directly constructed inertial trajectories and observer geometry independently test curvature and the phase sign.

The observer angle is fixed at 63.6° as in the existing project, based on the nebular geometry discussion in [Ng & Romani (2008)](https://arxiv.org/abs/0710.4168). It is a model input with unpropagated uncertainty, not an exact known viewing direction. A Gaussian sky bin of 1° is numerical angular averaging, **not** the single-particle relativistic beam. The output reports sensitivity to halving this bin. A fixed phase width σ=0.002 cycles represents unresolved directional/population structure; it is an explicit smoothing assumption, not a measured instrumental blur. Quantitative light curves are binned at 0.01 cycle for comparison.

## Calibration and test protocol

- Fit the **0.1–0.3 GeV and 1–3 GeV** profiles
- Reserve **0.3–1 GeV and 3–100 GeV** from optimization
- Fit one shared phase origin and select inclination/acceleration from the stated grid
- Freeze physical parameters and phase before scoring reserved data
- Fit only a nonnegative scale and constant in each reserved band

The software writes `protocol.json` before parameter search and `frozen-model.json` before scoring reserved data. Three leading training candidates are rechecked on a finer trajectory grid; the selected model is refined once more. No brightness coefficient per observed phase, catalog empirical template, or independent phase shift per band is used. A data-leakage test changes reserved counts and confirms the selected phase is unaffected.

The four profiles contain disjoint photon-energy ranges, with the repeated display cycle removed. [Smith et al. (2023), 3PC](https://arxiv.org/abs/2307.11132) describes the weighted light curves and their error prescription. Summing supplied variances retains the small-count correction from each original bin, so rebinned errors are conservative rather than exact event-level errors. See the [data audit](research/VELA-DATA-AUDIT.md).

**This is conditional pulse-shape testing.** The archive used here lacks detector response and exposure information required for quantitative spectral/absolute-flux fitting. Within-band response is assumed flat; separate normalizations erase cross-band flux constraints. A constant can represent both background and unmodeled steady emission. The χ²/(N−2) numbers are descriptive diagnostics, not calibrated probabilities after model selection. Earlier inspection of Vela's aggregate morphology also means this is not a fully blind discovery test.

## What remains necessary for a reconstruction claim

This calculation does not solve Maxwell evolution coupled to plasma, pair creation, field screening, reconnection, synchrotron radiation, inverse Compton scattering, or global spin-down energy accounting. It also does not predict absolute injection rates, luminosity, optical colours, radio/X-ray emission, or the TeV component. Gamma rays have no visible true colour; any animation palette would be illustrative.

A convincing Vela model needs physically validated magnetospheric fields and injection, detector-aware spectral comparisons, quantified degeneracies and parameter uncertainty, and genuinely independent observations. Even a successful pulse-shape match would not uniquely recover a three-dimensional plasma from one viewing direction. Until those tests pass, use **“Vela-tested reduced particle and curvature-emission model”**, not an unqualified **“Vela reconstruction.”**
