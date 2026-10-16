# Field source and absorbing-boundary benchmarks

The added experiments isolate two errors that a closed energy ledger cannot
detect. The baseline sponge creates magnetic fields from a pure-gauge vector
potential, and it reflects approximately 12% of the energy of a specified
normally incident wave packet. The reflected fraction persists under mesh
refinement. A separate exact-source vacuum experiment shows that the Maxwell
propagator converges once its source and outer data are specified consistently.
The legacy production solver has **not** been repaired or revalidated by these
experiments. Its source and sponge code are unchanged.

Implementation: `src/pulsar/field_benchmarks.py`; reproducible driver:
`run_field_benchmarks.py`; results and the complete time histories:
`outputs/field-benchmarks/results.json`. All numbers below come from those runs.

## Two distinct analytic references

The **point source** is an independently evaluated retarded magnetic dipole,
with unit moment, inclination 65°, and `c=Ω=R_LC=1`. If `n=r_vec/r`, dots mean time
derivatives of the rotating moment, and every moment is evaluated at `t-r`,

```
A = m×n/r² + m_dot×n/r
E = -m_dot×n/r² - m_ddot×n/r
B = [3n(n·m)-m]/r³ + [3n(n·m_dot)-m_dot]/r²
    + [n(n·m_ddot)-m_ddot]/r
L_point = (8π/3) sin²χ
```

`B` is evaluated directly from this formula when checking the continuum
reference; it is not defined as the numerical curl on both sides of an
identity. Independent centered time derivatives and the baseline spatial
derivatives check both vacuum Maxwell equations and both divergences. Over
`1<r<2.4`, Faraday/Ampère relative residuals fall from 4.63%/26.16% at N=32
to 0.56%/3.13% at N=96, consistent with second-order spatial convergence.
The relatively large Ampère error comes from differentiating the steep
near-zone magnetic field. The aligned limit, `E=-∂A/∂t`, and the integrated
Poynting luminosity at several phases and radii are independently tested.

The **finite conducting sphere** is the appropriate physical reference for
the intended baseline conductor. `deutsch_fields` implements the outgoing
magnetic dipole and induced electric quadrupole, with no net electric charge.
It obeys the prescribed radial magnetic field and tangential ideal-conductor
electric field at `r=R*`. The finite-radius luminosity ratio is

```
L_finite/L_point = 1/(1+a²) + 3a⁴/[5(a⁶-3a⁴+36)], a=R*/R_LC.
```

These modes and the flux formula correspond to [Pétri (2015), equations 24
and 46a](https://arxiv.org/pdf/1503.05307). The implementation is checked
against local vacuum Maxwell derivatives, surface conditions, and flux
integration at five radii. Surface relative residuals at R*=0.35 are below
`7×10⁻¹⁶`; integrated flux agrees with the closed expression to `10⁻¹¹` or
better. A freely chosen stationary electric monopole changes the charge but
not this radial luminosity. This distinction matters because the baseline
does not explicitly enforce a neutral conductor.

At R*=0.35 and χ=65°, the finite-sphere reference is **6.13205**, versus
**6.88129** for the point source. The legacy N=64 vacuum mean flux at r=1.5
was **3.98682**: its deficit relative to the appropriate finite-size reference
is **34.98%**, not the previously quoted 42.06% point-source discrepancy.
Finite radius therefore explains part, but not all, of that discrepancy.
The baseline's diffuse mask and poorly resolved surface do not exactly impose
the mathematical Deutsch boundary, so this is a target for source refinement,
not an assertion that every remaining percent has one known cause.

## Controlled exterior evolution

`ExactDipoleBoundary` retains the baseline curl operators and RK4 integration.
It starts from analytic point-dipole data, prescribes analytic A/E inside
r=0.65 and in the outer three grid layers at every RK substage, and evolves
the intervening vacuum. The bounded continuation inside r=0.2 is excluded
from the validated exterior. Plasma, conductor relaxation, spinup, filtering,
and the sponge are disabled in this controlled experiment.

| Grid | E error, exterior L2 | B error, exterior L2 | Maximum flux error |
|---|---:|---:|---:|
| 48³ | 8.40% | 8.44% | 2.66% |
| 64³ | 3.17% | 3.80% | 0.96% |

Values are at t=2π; the field shell is `1<r<2.4` and flux radii are
1.25, 1.75, and 2.25. Sampling and interpolating the analytic fields itself
introduces up to 0.48%/0.28% flux error at N=48/64. Discrete div(B) is near
machine precision because the numerical B remains a compatible curl.
The acceptance gate is 10% for final E/B error and every final flux error.
Both resolutions pass that gate; the complete early transient is retained
(N=48 reaches 17.15% E error). This is one-rotation propagation/maintenance
verification with externally supplied exact boundaries, not the generation
of radiation from a self-consistent plasma or a proof of long-term periodicity.

## A reproducible sponge failure and a bounded repair

The packet test uses the baseline 8-wide periodic domain, cubic layer starting
at |x|=2.8, RK4, central differences, and optional fourth-difference filter.
Initially `A_y=exp[-(x+1)²/(2×0.45²)]`, `B_z=D_x A_y`, and `E_y=B_z`.
Thus the initial wave has exactly zero leftgoing discrete characteristic.
Interior leftgoing energy is `∫[(E_y-B_z)/2]² dx`, divided by initial total
energy. Its maximum is measured after t=3 in |x|<2.5. Periodic leakage is
measured separately from rightgoing energy returning to the left interior
after t=6. The run ends at t=10.

| N | A/E sponge, strength 4: reflected energy | Matched E/B, strength 12: reflected energy | Matched E/B: periodic leakage |
|---|---:|---:|---:|
| 48 | 12.4% | 0 | 1.75×10⁻⁴ |
| 64 | 12.2% | 0 | 1.84×10⁻⁶ |
| 128 | 12.0% | 0 | 7.83×10⁻⁷ |
| 256 | 12.0% | 0 | 6.07×10⁻⁷ |

These rows disable the filter to isolate the absorber. Enabling the baseline
filter strength 0.08 leaves its reflected fraction close to 12%; filtering
does not fix this error. The new normal-incidence prototype evolves the
physical E/B pair and damps both equally. It preserves `E=B` exactly in this
one-dimensional problem, so its opposite characteristic remains zero to
floating-point precision. The plot clips those zeros at 10⁻¹⁷ for display.

Merely matching the two fields with the original strength 4 still leaks
0.81% of the packet energy through both periodic edge layers at N=256.
The analytic continuum transmission is `exp(-1.2×strength)`: 0.823% for
strength 4 and `5.57×10⁻⁷` for strength 12. The higher-resolution measured
transmission approaches those predictions. The gate is 0.1% for both peak
reflected energy and periodic leakage. The baseline fails; the strength-12
matched one-dimensional prototype passes at all tested resolutions.

There is also an exact gauge counterexample. A constant A_y has E=B=0.
Applying the baseline sponge for Δt=0.1 creates magnetic energy 0.10841 per
unit transverse area at N=64. The created field scales linearly with the
arbitrary gauge offset and its energy quadratically. This follows from
`curl(d A)=d curl(A)+grad(d)×A`, and cannot be explained by physical absorption.

## Source seam and safe production work

The baseline evolves `∂A/∂t=-E`, a zero-scalar-potential gauge, then resets A
through a finite-width mask to the rotating static potential while imposing
corotation E. The unmasked difference `∂A_static/∂t+E_corotation` is a gradient
in the continuum; physical rotating B and ideal-conductor E are compatible.
Masking that gradient can introduce curl. The diagnostic corresponding to
its first-order magnetic correction is 2.12%, 1.74%, and 1.92% of local seam
B per timestep at N=48, 64, and 96. This is not the measured total production
step error and does not uniquely explain the luminosity deficit. The actual
star radius spans only 2.1, 2.8, and 4.2 cells, respectively.

The demonstrated E/B absorber cannot be pasted into the current 3-D A solver.
In three dimensions, multiplying B by a spatial damping function generally
creates `div(B)=grad(d)·B`; reconstructing A after such a step also changes the
field globally. The normal-incidence characteristic identity is insufficient
for oblique incidence, spherical near-zone fields, or corner intersections.

A safe next production implementation should therefore:

1. Choose a consistent physical-field boundary formulation—compatible
   constrained transport or an explicit divergence-control formulation with
   a carefully specified absorbing layer. Preserve the baseline for comparison.
2. Implement a multidimensional matched layer or characteristic open boundary
   and its required auxiliary variables. The established PML construction is
   [Bérenger (1994)](https://doi.org/10.1006/jcph.1994.1159); this prototype is
   only a one-dimensional impedance-matching demonstration, not that algorithm.
3. Repeat gauge, normal/oblique packet, low-frequency, static dipole, and
   corner tests; measure active-domain divergence, reflection, and energy
   separately. A larger domain alone does not demonstrate absorption.
4. Impose stellar radial B and tangential conductor E with a compatible gauge
   or field update. Validate against the finite-sphere Deutsch reference,
   including source-radius and mesh refinement, before restoring plasma.
5. Re-run the plasma resolution, boundary-location, conductivity, and
   periodicity studies. The previous failed plasma gates remain failed until
   these new production runs pass; no emission or Vela gate is unlocked here.

## Reproduction and provenance

```
PYTHONPATH=src .venv/bin/python -m unittest tests.test_field_benchmarks -v
.venv/bin/python run_field_benchmarks.py
```

Ten tests cover the analytic limits, temporal gauge, both Maxwell equations,
both divergences, conductor surface conditions, published luminosity, source
substage timing, gauge failure, reflected characteristics, and absorption.
The driver records source SHA-256 values per run, the protocol and gates,
Python/NumPy/SciPy versions, UTC generation time, and full scalar histories.
The reanalysed legacy result records its input JSON hash and retains a null
legacy generation-source hash where that information is unknown. New analysis
hashes are never attributed to legacy snapshot generation. No large 3-D
snapshots are needed for these bounded-memory verification experiments.
