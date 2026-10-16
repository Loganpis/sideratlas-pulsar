# Absolute energetics and acceleration response

The acceleration gains used by [Barnard et al. (2022)](https://arxiv.org/html/2111.03405v1) prescribe energy gain per laboratory path length. Their example values cannot be assumed to reconstruct Vela when the numerical magnetosphere, injection surface and particle loading differ. We retain those departures explicitly and test the gains as model parameters.

`run_acceleration_scan.py` reuses recorded particle geometry. It transports the energy again for each pair of **global** inside/outside gains, includes curvature losses, and integrates the directional photon yield through the declared observer-latitude bin. There are no phase-dependent gain parameters or manual peak weights. A single population normalization is fitted to the seven detected catalog spectral bins; the eighth upper-limit bin is checked separately. Every trial remains a development calculation on an unreleased field.

The phase-integrated sky calculation is independently checked against summing the full conservative phase projection, including forward/reversed sky-bin crossings and stationary rays. It does not assume that stored bin averages are nodal samples.

## Absolute particle supply

For a phase-averaged yield per injected primary per steradian, the photon flux at Earth is

\[
F_b=\frac{\dot N}{d^2}Y_b
\]

Both magnetic hemispheres are already included in the globally normalized seed population. An extra factor of \(4\pi\) would double-count solid angle. The fitted spectral normalization therefore gives \(\dot N=d^2 A\). The all-direction radiated power is \(L_\gamma=\dot N\sum_i w_i E_{\gamma,i}\), and the modeled supplied power includes prescribed electric work and injected kinetic energy. A positive final kinetic-energy flux is retained in the ledger.

The audit adopts **290 pc**, following the Barnard model, and reports 270/310 pc sensitivity scenarios as requested by the audit code. That bracket is not a measured confidence interval. The adopted moment of inertia remains \(10^{45}\,\mathrm{g\,cm^2}\). A budget below the spin-down power is a necessary consistency check; it does not establish pair supply, current closure, a successful spectrum or a complete broadband model.

## Field amplitude

The current geometry runs use the ideal force-free oblique-point-dipole spin-down normalization. An enlarged finite numerical star need not obey that torque exactly. The separate `field_scale_from_source_power` function provides the dimensionally consistent alternative

\[
B_{\rm scale}=\left[\frac{L_{\rm sd}}{cR_{\rm LC}^{2}L_{\rm code}}\right]^{1/2}
\]

It may be used only after measuring stable source power outside the conductor, before artificial absorption, with source-size, resolution and boundary sensitivity established. It has **not** been silently applied to the historical emission runs. Analytic tests show that it recovers the point-dipole limit and the expected square-root scaling.

## Reproduce

```sh
/opt/anaconda3/bin/python run_acceleration_scan.py \
  outputs/emission-quadrature/n65-a64-r3-dt002-v2 \
  --output outputs/acceleration-response/new-scan
/opt/anaconda3/bin/python plot_acceleration_scan.py outputs/acceleration-response/new-scan
/opt/anaconda3/bin/python run_population_audit.py outputs/outer-coffee-n97-surface-anchored
```

The initial recorded scan is in `outputs/acceleration-response/n65-initial-global-scan/`. Its best of four tested gain pairs improves the descriptive spectral statistic from about 25,580 to 5,428, still a substantial mismatch. These four trials are a sensitivity scan, not an optimization or a successful Vela fit. Further tuning on a field that fails numerical checks would not certify the reconstruction.
