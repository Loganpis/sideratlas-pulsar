# Conservative atmosphere transfer: actual implementation and acceptance

The initial nodal hybrid iteration was rejected because its abrupt switch
between local and flux equations admitted a discontinuous transported flux.
The later nodal simultaneous solver was useful diagnostically but did not
enforce a common interval energy equation when opacities changed.

`conservative_magnetic_column.py` instead uses these finite-volume choices:

- Mass/temperature nodes bound each cell; pressure and density are evaluated
  at the cell centre using the ideal fully ionized H EOS
- Absorption, cyclic scattering response and polarization powers are held
  constant within a cell
- The LTE thermal source `epsilon B_nu(T)/2` varies linearly between the two
  cell faces; coherent scattering uses the **exact cell-averaged** radiation
- For `t=kappa Delta_m/mu`, transmission is `exp(-t)` and the exact mean
  transmission is `phi=(1-exp(-t))/t`, with cancellation-safe small-t series
- The mean source is `epsilon (B_left+B_right)/4 + scattering(J_mean)`;
  its known thermal slope is `epsilon (B_right-B_left)/2`
- Upward and downward formal integrals include that slope analytically
- The scattering system closes through three cyclic radiation moments,
  reducing the matrix to `3*Ncell` without discarding the mode exchange

For each frequency, angle and mode the formal solution is the exact solution
of this **specified finite-cell approximation**. Summing gives

`F_right - F_left = 4 pi Delta_m sum_i integral dnu dmu kappa_abs,i (J_mean,i - B_mean/2)`

Coherent scattering cancels in this sum. This identity is checked on an
arbitrary temperature profile, independently of the equilibrium iteration.
At typical test values the relative ledger discrepancy is below `7e-11`;
in a deep converged model cancellation of large LTE intensities gives a
residual around `5e-7` of the target flux.

The nonlinear solve retains **both** flux errors at every face and local
relative energy-balance errors in every cell. No emergent spectrum is
renormalized. The 32-node run converged with optimizer `gtol` termination in
348 s: maximum face-flux error `3.19e-13`, local relative balance `3.36e-12`.
The discretization can admit an alternating nodal temperature mode, so these
equation residuals alone do not establish a resolved atmosphere. A separate
monotone-temperature variant and a 64-node refinement evaluate this issue.

`export_conservative_surface.py` evaluates the solved scattering source at
explicit outgoing emission cosines, including `mu=1` and the analytic
grazing limit. The exported two-mode sum is `I_lambda` in SI units:

`I_lambda[W m^-2 sr^-1 m^-1] = 1e-3 I_nu[erg s^-1 cm^-2 sr^-1 Hz^-1] c/lambda_m^2`

It retains the actual broad-spectrum solution and explicit local 300–650 nm
nodes. The user-facing camera's redshift/Doppler interpolation must stay
inside those local optical nodes. No X-ray extrapolation or RJ patch is used.

`compare_surface_intensities.py` compares absolute radiance tables under
refinement: wavelength/angle point errors, band-integrated intensity by
angle, and hemispheric band flux. It neither fits an exposure nor rescales
the spectra. The default local 330–580 nm interval covers the 400–700 nm
observer band for the chosen compactness, with a small rotational margin.
The declared one-percent gate excludes only `mu<.02`, covering `.0004` of
projected emitting area; grazing values remain in the export and the
hemispheric integral.

The first 32-node bottom-boundary diagnostic gives optical total depth above
411 and `integral sqrt(3 kappa_abs kappa_tot) dm` above 712. The optical
diffusion-boundary correction is below `.00068`; these are **diagnostics**,
not a substitute for extending the bottom boundary and comparing the output.

Physical validity remains independent of numerical balance. The atmosphere
uses unity prescribed magnetic Gaunt factors, fully ionized ideal H, linear
transverse modes, coherent scattering and a normal homogeneous field. Atomic
opacity, dispersive/refraction effects, vacuum conversion, conduction and
irradiation are omitted. No independent magnetic optical spectrum benchmark
has passed. Along the first converged-equation column, public EOS tables give
atomic fractions up to 62% at an intermediate `T≈202 kK` layer, while the
ordinary optical mode forms near 90 kK below the tables' temperature domain.
The extraordinary optical mode forms deeper at roughly 1.1–1.9 MK, with
small tabulated neutral fractions, but propagates where plasma energy greatly
exceeds optical photon energy. These are material unbounded approximations.
Numerical convergence does not make the output a validated Vela atmosphere.

## Selected final 64-node approximate input

The selected file is
`outputs/optical-finish/atmosphere/polished-conservative-64/surface-intensity-si.npz`,
SHA256 `0915abb9b0e9f4cedfd328e1a2e5e905a42b7c56ac004bd02245b688552d6f37`.
It has 64 mass/temperature nodes, 63 material cells, four Gauss emission
cosines per hemisphere, and 104 frequencies (64 broad-grid points plus 40
explicit local optical points). The export adds explicit outgoing angles
using the calculated incoming radiation field, rather than extrapolating
the four quadrature intensities. It sums both modes and retains SI units.

The first 64-node iteration was frozen and passed to a bounded polish of the
same conservative equations. Weighting local residuals by 100 changes the
numerical conditioning only. The polish ended after 131 residual evaluations
and 57.50 seconds when the explicit acceptance criteria were satisfied:

- Maximum face-flux error: `9.4601933e-4`
- Maximum relative local energy-balance error: `4.8590916e-5`
- Maximum energy-ledger discrepancy: `2.4716138e-7` of target flux
- Independently integrated emergent flux: `1.000051655` times target

The radiance was never normalized after calculation. The thermal structure
is **not smooth/monotone**: the largest downward nodal log-temperature step
is 0.25885. A separately restricted monotone sensitivity run was stopped
unfinished; it supplies no accepted alternative surface table.

The direct absolute-radiance comparison against the converged 32-node case
over local 330–580 nm, without changing exposure, gives:

- Hemispheric optical flux: **+4.6802%**
- Median non-grazing pointwise intensity change: **4.6988%**
- Worst non-grazing pointwise change: **5.8824%**
- Worst angle's band-integrated radiance change: **5.3995%**

Thus `one_percent_optical_refinement_pass=false`. The comparison excludes
only `mu<.02` from pointwise maxima; those rays remain in the hemispheric
integral, and their projected area fraction is `.0004`.

The final bottom diagnostics give optical total depth above 338 and the
thermalization diagnostic above 585, with optical diffusion-boundary
correction below `.000648`. A fully relaxed bottom-extension test is not
claimed. Source-microphysics, partial-ionization and refraction errors remain
unbounded. All final records keep `scientific_validation_pass=false`,
`refinement_pass=false`, and `physical_model_approximate=true`.

Reproduce the bounded polish in a **new output folder**, preserving the
published table bytes:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/optical-finish/atmosphere/polish_conservative_column.py outputs/optical-finish/atmosphere/checkpoint-refinement-64/input-checkpoint.npz --depths 64 --max-nfev 20 --output /tmp/reproduced-atmosphere-64
python research/optical-finish/atmosphere/export_conservative_surface.py /tmp/reproduced-atmosphere-64 --allow-unvalidated-physical-model
```

Archive and copy only the files in `FINAL-SURFACE-HANDOFF.json`. The local
vendor checkout and external EOS tables are deliberately excluded; no
redistribution license was established for them in this work.
