# Independent review of the resistive Maxwell prototype

Review date: 2026-10-13

The new field solver evolves the electromagnetic field with a prescribed
macroscopic plasma current. This is a substantive advance over a fixed vacuum
dipole. It is **not** a kinetic particle-in-cell calculation, a pair-cascade
calculation, or an independently validated Vela magnetosphere.

## Primary-source correspondence

[Kalapotharakos, Harding & Kazanas (2014), Equation 9](https://arxiv.org/pdf/1310.3545)
specifies a charge-drift current plus a conducting current parallel to the
magnetic field. Its invariant electric-field term makes the drift subluminal
even in electrically dominated regions. The paper's Gaussian-unit Maxwell
equation contains a factor of `4π` multiplying current. In this code's
rationalized units, the damping parameter is therefore
`κ = 4π σ/Ω`, not `σ/Ω`. Default damping values 80 and 5 correspond to
`σ/Ω ≈ 6.37` and `0.398`. Agreement in closure form does not reproduce the
paper's conductivity map, boundaries, mesh, or production solution. Conductivity
still represents unresolved plasma processes; the paper does not establish a
unique microphysical conductivity for Vela.

## Analytic and discrete benchmarks

`tests/test_plasma.py` independently checks:

- A periodic electromagnetic plane wave, with conductor, plasma, sponge, and
  filter disabled. Its expected numerical frequency is
  `ω_h = sin(kh)/h`, the central-derivative dispersion relation. For `N=32`,
  `k=π/2`, `h=0.25`, `Δt=0.0625`, and 40 steps, relative field error is
  `2.6714×10⁻⁶`; relative electromagnetic-energy change is `−4.2551×10⁻⁷`.
  Comparing to the continuum frequency instead would incorrectly attribute
  spatial dispersion to time-integration error.
- Second-order convergence of the magnetic-field derivative to a smooth
  continuum solution under mesh refinement.
- `div(curl(A))=0` for arbitrary periodic grid data, to floating-point accuracy.
  This is a useful structural constraint; it does not show the field geometry
  is astrophysically correct.
- Exact exponential damping of a field component parallel to a fixed magnetic
  field, including zero and stiff conductivity, composition of half-steps,
  and the independent analytic electric-energy loss.
- Drift speeds for magnetic-dominated, electric-dominated, null, parallel,
  magnetic-zero, fully zero, and random fields. Speeds stay at or below `c`.
- The filter's known Fourier-mode attenuation and preservation of the discrete
  divergence constraint.
- Complete accounting of the energy changes made by the separate numerical
  interventions, and finite diagnostics in the zero-field limit.
- Absence of closure current and Joule power in vacuum mode.

The energy ledger telescopes changes measured before and after each operation.
Its closure is an accounting check, **not** an independent conservation law
test. The periodic wave's independently predicted energy is a separate test.
The ledger's `maxwell_drift_error` measures time-integration error for an ideal
periodic drift update, because drift current does no work. A small total ledger
residual cannot excuse a large value of that term.

## Unresolved numerical and physical assumptions

1. **Macroscopic current closure.** The conductivity is prescribed, not derived
   from evolving distribution functions. No pair creation, charge starvation,
   thermal pressure, particle inertia, or kinetic reconnection is solved.
   Evolving a conducting field is not the same as solving particle feedback.
2. **Enlarged and poorly resolved star.** With `R*=0.35 R_LC`, the default star
   is roughly 124 times Vela's `12 km / 4246.5 km` radius ratio. At `N=64` in a
   domain of width `8 R_LC`, the stellar radius spans only 2.8 grid cells. The
   surface, polar caps, open flux, and current sheet cannot be assumed resolved.
3. **Boundary forcing.** A smooth mask replaces the stellar vector potential
   and relaxes the electric field toward corotation. This injects/removes
   energy and changes the field outside a mathematical conductor surface.
   Conductor work is recorded; its convergence remains to be established.
4. **Absorbing layer.** The outer domain is periodic with a sponge. Damping
   `A` also contributes `grad(damping) × A` to `B`; it is gauge-dependent and
   can add magnetic energy. The signed `sponge_exchange` term is appropriate.
   This layer is not a perfectly matched boundary, and a dedicated outgoing
   wave-packet reflection test is still required. Initial dipole fields are
   also not exactly periodic at the outer faces.
5. **Numerical regularization.** The fourth-difference filter and an optional
   inner-region `|E|≤|B|` projection remove energy separately from resistive
   heating. Their effect on current-sheet thickness, accelerating fields,
   open flux, and pulse predictions needs convergence tests.
6. **Finite conductivity profile.** The transition is spherical at `r=R_LC`,
   while the light cylinder is cylindrical. This is an explicit model choice,
   not a computed boundary between open and closed field lines. It does not
   enforce exact force-free conditions throughout the inner region.
7. **Rotating steady state and radiation.** A short stable run does not establish
   a periodic magnetosphere. Multiple rotations, mesh/boundary/conductivity
   variations, particle coupling, resolved trajectory curvature, emission
   energy balance, and comparison to Vela data are separate requirements.

The supported claim is **“verified numerical resistive-Maxwell prototype.”**
The field benchmarks do not replace the emission model's failed observational
validation and must not be presented as evidence of a Vela reconstruction.
