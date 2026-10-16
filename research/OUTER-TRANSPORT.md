# Transport through supplied fields beyond the light cylinder

`src/pulsar/outer_transport.py` supplies a prescribed-velocity, ultrarelativistic
guiding-centre calculation. The old corotation construction and its 0.95 R_LC
cutoff are unchanged. This module requires an explicit field provider and does
not create a plasma solution, current closure, pair distribution, or absolute
particle injection rate.

## Published prescription and its limits

Kalapotharakos, Harding & Kazanas (2014), [Eqs. 6, 12 and 13](https://arxiv.org/html/1310.3545v2),
specify

\[
\boldsymbol\beta={\mathbf E\times\mathbf B\over B^2+E_0^2}
 + f\widehat{\mathbf B},\qquad
f=s\sqrt{1-|\boldsymbol\beta_d|^2},\quad s=\pm1,
\]

where `E0>=0` is the invariant electric field. Their outgoing-particle
prescription fixes the speed at c. Here the field-parallel branch and charge
sign are explicit independent inputs. This allows examination of either
stream and preserves negative electric work. The energy equation is

\[
{d\gamma\over ds}={e E_{\rm work}\over m_ec^2}
 -{2e^2\gamma^4\over3m_ec^2\rho_c^2},\qquad
 E_{\rm work}=(q/e)\mathbf E\cdot\boldsymbol\beta.
\]

The paper's `f E_parallel` factor is thus retained through the actual velocity
projection. Radiation uses laboratory trajectory curvature and the existing
laboratory curvature spectrum; an additional Doppler transformation would
double-count the frame transformation.

This is not an integration of the full Lorentz force or resolved gyration.
Energy does not feed back into the prescribed unit velocity. A magnetic null
raises because Bhat is undefined. A resolved near-null current sheet may still
invalidate the guiding-centre assumption even though the algebra returns a
finite speed. Neither magnetization nor pair support is guaranteed by these
equations. Fields must be continuous and adequately resolved over the orbit
and derivative stencils; a discontinuous sheet face cannot be treated as
resolved dissipation.

The optional postprocessed `NearForceFreeFIDO` field has a deliberately sharp
screening interface at cylindrical radius 1. A direct derivative across that
interface produced artificial curvature growing as the derivative step was
reduced. The field now supplies a stencil validator: RK stages and curvature
stencils touching or crossing that interface raise an error. It cannot be used
for complete crossing trajectories until a resolved transition prescription
and its sensitivity are supplied. No transition thickness is invented by the
code. This restriction does not affect the smooth-geometry Barnard diagnostic.

The alternative `Barnard2022Work` follows
[Barnard et al., equation 2](https://arxiv.org/html/2111.03405v1): the gain per
laboratory path length is directly 0.04 cm⁻¹ inside spherical radius 1 and
0.25 cm⁻¹ outside for the selected scenario. These are prescribed effective
gain rates, not a global electric field or a predicted population. They do
not receive an additional field-parallel velocity projection. An earlier
implementation applied that extra factor; it has been corrected before any
successful end-to-end result. This convention is distinct from the signed
work of a supplied physical electric field in the equation above.

The Aristotelian/radiation-dominated direction
`[E×B ± (B0 B + E0 E)]/(B²+E0²)` is a different velocity assumption
([Gruzinov 2013, Eqs. 1–2](https://arxiv.org/html/1303.4094v1)). It is **not** an
automatic fallback here. Its charge-associated choice would impose a
radiation-dominated force-balance interpretation; adopting it would require
separate validity and energy treatment rather than silently changing Eq.12.
This alternative still needs a plasma-production model to determine the
electron and positron populations.

## Interface and units

`fields(position,time) -> (E,B)` accepts Cartesian positions in R_LC and time
in Omega^-1, with array time broadcast over the position batch. E is returned
in statvolt/cm and B in gauss. Shared field normalization is essential; there
is no hidden scaling. `RotatingSnapshot` and `NearForceFreeFIDO` implement this
interface. A source's diagnostic/release status is unchanged by transport.

`select_branch(fields,seeds,time,flow=+1/-1)` chooses the most outward/inward
injection root and rejects an unavailable radial direction. In a wind outside
the light cylinder both roots can move outward. The chosen branch remains
fixed along its trajectory; charge is never flipped to manufacture positive
work. This does not enforce outflow at every later point, nor claim that a
chosen electron/positron stream is supplied by the magnetosphere.

`trace_outer(fields,seeds,branch,charge_sign,config=OuterTraceConfig(...),
start_time=0.,seed_weight=...)` integrates fourth-order Runge–Kutta laboratory
trajectories. Returned arrays use `(segment,seed[,xyz])`. `position`, `time`,
`direction`, `rho`, `e_work`, and phases are midpoint quantities. `dt` is the
segment residence time; multiply it by R_LC in cm to obtain ds. Multiply `rho`
by the same R_LC for radiation. Inactive padded segments have zero dt/work,
zero direction, infinite rho and NaN phases. Supplied seed weights are
preserved, including their normalization.

Spherical `r_min`/`r_max` exits are clipped by event bisection. The provider
must cover the RK and curvature stencils, including a dt-wide buffer beyond
the termination surfaces. `status` distinguishes `r_min`, `r_max`, and
`horizon`; a time horizon is not physical escape. Curvature is computed from
the material derivative `partial_t beta + beta·grad beta`; its finite
difference step must be checked separately from the trajectory step. Only
floating-point derivatives smaller than `32*epsilon/h` are treated as zero.

The low-level `all_escaped` flag means all paths reached one of the two radial
boundaries. The outer diagnostic additionally rejects any `r_min` return:
that radius is an artificial injection-region cutoff and the remaining
inward transport is undefined. It records the returning seed weight instead
of dropping those seeds and renormalizing the survivors. The projection does
not yet include photon absorption or stellar occultation; successful runs
report the fraction of emitted energy directed inward. The enlarged
numerical star is not automatically treated as Vela's physical opaque radius.

`signed_transport_curvature(ds_cm,rho_cm,e_work_statvolt_cm,gamma0,...)`
returns the existing `CurvatureTransport` energy and photon-yield structure.
It freezes work and curvature on each segment. Positive work uses the
existing exact solver; negative work uses a monotone inversion of the
analytic integral for `gamma'=-a-b gamma^4`, evaluated with a reciprocal tail
to avoid catastrophic cancellation. Integrated emitted energy follows the
same work/particle-energy ledger. Midpoint or three-node power-weighted
spectral quadrature is available. Rapid spectral evolution still requires
segment refinement.

`gamma_floor=10` is a declared validity guard for the c-speed approximation,
not physical thermalization. Crossing the selected floor raises an error.
No rest energy is added and no gamma is silently clamped. A caller that needs
slower particles must supply another transport approximation or explicitly
truncate the modeled stream at its validity boundary.

## Arrival observables

For a fixed observer direction n, `arrival_phase` returns
`(t-n·r)/(2*pi) mod 1`; visibility remains a separate angular selection.
`fixed_phase` in a bundle uses n equal to each emitted photon's local
direction. For a stationary pattern rotating at Omega, `sky_phase` maps a
trajectory sample into the observer-azimuth-zero sky map:

`(t-atan2(v_y,v_x)-v·r)/(2*pi) mod 1`, with `mu=v_z`.

The explicit t term is required for laboratory-integrated trajectories.
The extra azimuth shift requires rotating-pattern symmetry and must not be
applied to arbitrary transient fields. Axial rays have undefined azimuth;
`sky_phase_valid=False` and NaN phase make that limit explicit. Filter active
and valid samples before binning phases.

## Independent checks and interpretation

Tests cover uniform crossed/parallel/electric-dominated fields, branch and
charge signs, outgoing/incoming streams, escape-time clipping, explicit
time-dependent curvature, independently integrated laboratory trajectories,
observer geometry, and signed energy flow against adaptive ODE integration.
The critical split-monopole regression verifies that its outgoing Eq.12
trajectory is radial, its laboratory curvature is zero, and curvature
emission vanishes despite spiral magnetic field lines.

Passing these tests validates implementation of the stated approximation.
It does not establish convergence of a particular numerical magnetosphere,
current-sheet resolution, observed pulse agreement, or a kinetic pair model.
Any full prediction still needs explicit field provenance, physical scaling,
injection, convergence of geometry/energy/spectral bins, and independent
observational comparison.
