# Generic rotation-powered pulsar: bounded consistency audit

The adopted parameters are mutually plausible for a young rotation-powered
pulsar. No energy-budget or scalar-parameter problem requires a new particle
run. The product should be described as a **calculated, explicitly approximate
optical visualization of a pulsar**. Vela supplies physical reference scales;
its luminosity, spectrum and pulse profile are not acceptance targets.

This audit reads the existing `production-spatial` report and archive. It
performs scalar calculations and a retained-event occultation check only. It
does not rerun the model, fit observations or establish convergence. Exact
inputs, values, source hashes and scope are recorded in
`outputs/optical-finish/physical-consistency.json`.

## Compatible adopted scales

The model uses period 89.371 ms, mass 1.4 solar masses, radius 13 km, **polar**
field 3 × 10¹² G, magnetic inclination 75° and observer inclination 65°. The
light-cylinder radius is 4,264 km. The classical near-zone optical population
is activated between 0.08 and 0.18 light-cylinder radii, below the light
cylinder; the retained archive's maximum corotation speed is 0.191 c.

The surface camera and atmosphere agree on compactness and redshift:
`2GM/(Rc²) = 0.3180423`, `g_grav = 0.8258073`. The atmosphere's local effective
temperature is 847,655 K, giving 700,000 K at infinity. Its input gravity
1.3313 × 10¹⁴ cm s⁻² agrees with `GM/(R² g_grav)` to 3.9 parts per million.
This agreement is bookkeeping consistency, not validation of the approximate
atmospheric opacities. The spherical-star approximation is reasonable at this
spin: the Newtonian centrifugal/gravitational acceleration ratio at the
equator is only 5.84 × 10⁻⁵; surface rotational speed is 0.00369 c in the local
static frame.

The chosen pair multiplicity is 6,000 **pairs** per primary, hence 12,000
leptons per primary. Initial gamma 300 and pitch angle 0.1 rad are explicit
injection assumptions, not measured population parameters. The existing
transport report has maximum quantum parameter chi = 9.29 × 10⁻⁵,
gyro-radius/field-gradient-length = 1.38 × 10⁻¹⁰, and gyroperiod/cooling-time =
8.50 × 10⁻⁵. These support the classical, gyro-averaged guiding-centre treatment
for this particular prescribed population. They do not establish a global
magnetospheric solution. The low-harmonic radiation kernel is retained rather
than replacing it with an invalid smooth synchrotron continuum.

## Rotational and particle energy scales

For the **polar** field convention, the Gaussian dipole moment is
`mu_G = B_p R³ / 2`. Converting the force-free reference luminosity of
[Spitkovsky (2006)](https://arxiv.org/abs/astro-ph/0603147) to SI gives

```
K = pi B_p² R⁶ Omega⁴ / (mu_0 c³)
L_force-free = K (1 + sin² alpha) = 1.9035 × 10³⁰ W
L_vacuum = K (2/3) sin² alpha = 6.1250 × 10²⁹ W
```

Using an equatorial field in this polar-field formula would cause a factor-four
mistake. These luminosities are analytic reference scales for the adopted star;
the optical model does not calculate its global electromagnetic spin-down.

The existing selected-cap supply is 1.19154 × 10³⁶ leptons s⁻¹. Multiplying once
by `(gamma - 1) m_e c²` gives local injected kinetic power
2.9168 × 10²⁵ W: 1.53 × 10⁻⁵ of the force-free scale, or 4.76 × 10⁻⁵ of the
vacuum scale. Its computed local synchrotron loss is 2.8285 × 10²⁵ W, or
96.97% of injected kinetic power. It stays below the supplied kinetic energy.
The total rest-energy supply is another 9.76 × 10²² W. Even including rest
energy and modest corotation-frame corrections leaves a wide margin below
the rotational-power scale.

These are local-source versus reference-scale comparisons, not an exact
global conservation ledger. The injection cap flux is normalized with the
chosen static-dipole Goldreich–Julian density
`rho_GJ = -2 epsilon_0 Omega dot B`, not a solved current-closure condition.
The initial gamma, pitch and activation layer represent unspecified earlier
acceleration/heating. Their physical motivation is discussed with
[Harding et al. (2018)](https://arxiv.org/abs/1811.11157); its full cascade and
radio absorption calculation is not reproduced.

The effective-temperature surface luminosity at infinity is 4.2398 × 10²⁵ W.
It is small relative to spin-down and can be residual stellar cooling; it need
not be deducted from a rotation-powered particle accelerator. No stellar age
or measured thermal luminosity is being inferred.

## Material limits to retain in the release

- The optical field is a prescribed rotating near-zone dipole. Its particles
  follow calculated guiding-centre transport with cooling and focusing; its
  fields, pair creation, acceleration and currents are not evolved together.
  The separate reduced PIC calculation must not be advertised as the source
  of these physical optical photons
- Flat outer-photon propagation is combined with a Schwarzschild surface
  camera and local SR aberration/Doppler shifts. At the innermost retained
  optical source `r_s/r = 0.01145`, a useful weak-field scale, **not** a rigorous
  lensing or phase error bound. No optical/thermal/radio phase alignment is
  predicted
- All 5,271,172 retained event centres were checked for far-side straight rays
  intersecting the apparent stellar radius 15.742 km. None do; the blocked
  centre-source Y energy fraction is zero. This does not bound occultation of
  every finite support corner or strongly bent ray
- The uniform normal-field hydrogen atmosphere has passed its declared
  numerical equilibrium check, with documented physical opacity/EOS omissions
  and finite depth sensitivity. It is an approximate transfer calculation,
  not a validated surface-atmosphere library. A solved equilibrium is not
  evidence that these omissions are negligible
- Finite spatial supports, finite phase sampling, fixed exposure, a camera PSF,
  crop, tone mapping and slowed playback are declared approximations or
  presentation choices. The result is not a naked-eye or calibrated telescope
  image, a Vela reconstruction, or a reproduction of NASA's full PIC model

These limitations are important disclosures, but do not require fitting Vela
or spending more computation before producing the requested generic-pulsar
Reel.

## Limited camera review

No material units or double-counting issue was found in the reviewed camera
paths. Per-lepton XYZ power is multiplied by physical lepton number and source
time exactly once. The source's factor 683 is not reapplied. SI surface
specific intensity per metre is divided by 10⁹ to obtain intensity per nm.
The surface uses `I_lambda,obs = g⁵ I_lambda,em` at emitted wavelength
`g lambda_obs`. Outer arrival-time compression is integrated through source
and arrival windows; it is not an extra arbitrary Doppler multiplier.
This is a bounded code review, not independent end-to-end validation.
