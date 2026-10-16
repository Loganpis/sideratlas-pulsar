# NO.002 — Outer-wind field and pulse calculation

This milestone is a verified **reduced outer-wind benchmark**. The field is an analytic solution evaluated numerically; it is not the result of evolving a dipolar force-free solver or a PIC plasma. The photon calculation integrates an explicitly prescribed emitter. No NASA particle data or simulation code has been reused.

The calculation and its independent checks are reproducible. The observed Vela pulse shape is a separate test, and the baseline configuration fails that test. This distinction must remain in any caption or animation based on these outputs.

## Units and domain

Use Gaussian electromagnetic units with c = Omega = R_LC = 1, and an arbitrary magnetic amplitude B_LC = 1. Time in the code is Omega times physical time; an arrival phase of one is 2 pi in code time. Distances are in light-cylinder radii. A rounded period of 89 ms gives R_LC = 4246.497 km. It only supplies a scale, not an epoch-resolved Vela ephemeris.

The emission domain is 1.5 <= r/R_LC <= 5. There is no resolved stellar surface, closed dipolar region, pair cascade, or accelerating gap. Calling the wind approximation exact refers to its idealized field equations, not to the physical magnetosphere of Vela.

## Electromagnetic solution

In spherical coordinates about the spin axis, define magnetic obliquity chi and

\[
u=\cos\theta\cos\chi+\sin\theta\sin\chi\cos[\phi-t+r],\qquad S=\operatorname{sgn}(u).
\]

The current sheet is u = 0. In the outgoing, ultrarelativistic split-monopole limit,

\[
B_r=S/r^2,\quad B_\theta=0,\quad B_\phi=-S\sin\theta/r,
\qquad E_\theta=B_\phi,\quad E_r=E_\phi=0.
\]

Away from the sheet,

\[
4\pi\rho_e=-2S\cos\theta/r^2,\qquad4\pi\mathbf J=(4\pi\rho_e)\mathbf e_r,
\]

\[
\mathbf E\cdot\mathbf B=0,\quad B^2-E^2=r^{-4}>0,
\quad\rho_e\mathbf E+\mathbf J\times\mathbf B=0.
\]

The integrated radial electromagnetic luminosity is 2/3 in these normalized units, independently of radius and magnetic tilt for fixed monopolar flux. The net magnetic flux through a sphere is zero. The dipolar spin-down factor (1 + sin^2 chi) is **not** the appropriate benchmark for this fixed-flux split-monopole solution.

We evaluate Maxwell residuals with centered Cartesian differences at 512 reproducible points safely away from u = 0. Three step sizes test second-order convergence. At the moving discontinuity, we independently test the magnetic-normal and Faraday jump conditions. These checks do not solve sheet pressure balance or reconnection.

Field source: [Bogovalov (1999)](https://arxiv.org/abs/astro-ph/9907051), in the outgoing force-free limit. The test equations above follow directly by differentiating the stated fields.

## Emission prescription

The bulk streaming velocity is selected from the velocities compatible with the limiting ideal field:

\[
\mathbf v=\frac{\mathbf E\times\mathbf B}{B^2}
+s_\parallel\sqrt{1-\Gamma^{-2}-\frac{E^2}{B^2}}\,\widehat{\mathbf B},
\]

where the parallel sign selects outward motion. We check both |v|^2 = 1-Gamma^-2 and E + v cross B = 0. Gamma = 10 is a **chosen bulk Lorentz factor**, not the Lorentz factor of individually accelerated gamma-ray-emitting particles. It must exceed sqrt(1 + r_max^2) for this prescription to exist throughout the domain. Within the radiating sheet layer, this velocity is prescribed by continuation of the exterior solution.

Assume isotropic comoving spectral emissivity

\[
j'_{\nu'}=j_0\,r^{-q}\exp[-u^2/(2w^2)]\,{\nu'}^{-\alpha},
\qquad q=3,\quad w=0.05,\quad\alpha=1.
\]

The width w is in the dimensionless sheet function, not km or a constant physical thickness. The Gaussian layer is an emission prescription; smoothing it does not replace the discontinuous field by a self-consistent dissipative solution. No synchrotron or curvature spectrum, cutoff, or particle distribution has been calculated. Its normalization is arbitrary, and the radiated energy is not subtracted from the field.

For observer direction n and Doppler factor D = [Gamma(1-v dot n)]^-1, emissivity invariance j_nu/nu^2 gives j_nu = D^(2+alpha) j'_nu apart from the factored frequency dependence. At fixed arrival time,

\[
F_\nu(t_{\rm obs})\propto\int j_\nu(\mathbf r,t_{\rm em})\,d^3r,
\qquad t_{\rm em}=t_{\rm obs}+\mathbf n\cdot\mathbf r.
\]

The common source-to-observer distance is removed. The integration uses the laboratory volume at emission events selected by arrival time; adding a separate moving-blob boosting factor would double count a transformation. The spectral power law makes the normalized shape independent of energy in this simplified model. Comparing it with Fermi counts tests morphology only; it does not test a predicted spectrum or detector response.

Transfer and geometry references: [Pétri (2011)](https://arxiv.org/abs/1011.3431) and [Pétri (2013)](https://academic.oup.com/mnras/article/434/3/2636/1045892). These motivate the retarded volume integration; our simple isotropic Gaussian emissivity is a declared choice, not a reproduction of their complete emission models.

## Pulse and quadrature checks

For an infinitely narrow radial beam, the two sheet-crossing phases obey

\[
\Delta=\pi^{-1}\arccos|\cot\chi\cot\zeta|.
\]

This applies only when |cot chi cot zeta| < 1; otherwise no two distinct sheet crossings exist. The code independently evaluates retarded emission events along the sightline and measures the resulting peak separation. It also checks periodicity, an aligned rotator, a polar observer, and a deliberate missing-delay negative control.

The full finite-beam pulse uses Gaussian quadrature in radius and observer-centered cos(angle), plus periodic azimuthal quadrature over **the entire sphere**. Four progressively finer grids test shape and mean-flux stability. Additional radial and phase refinements test separate errors. The output contains 1024 arrival-phase samples, with one cycle represented on [0,1); the duplicate endpoint is omitted.

## Observation comparison and provenance

The data come from [Fermi 3PC](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/), using the [July 2025 corrected archive](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/3PC_ProfileData_20250718.tgz). Read extension GAMMA_LC by name and use GT100_WtCnt with Unc_GT100_WtCnt. The original file, archive URL and SHA-256 hashes are saved under data/.

The file repeats one rotational cycle for plotting. Only the first cycle is used; its 800 bins are rebinned to 200 by summing counts and variances. A nonnegative scale, a nonnegative constant and the relative phase origin are fitted by weighted least squares with a 4096-point phase search. Bin-integrated predictions are approximated with 16 subpoints per observed bin. No physical parameters are fitted, and no absolute radio/gamma phase lag is claimed.

The adopted observer angle of 63.6 degrees comes from [Ng & Romani (2008)](https://arxiv.org/abs/0710.4168), a model-dependent Vela nebular-torus inference. The chosen 60-degree magnetic obliquity is a benchmark value. The observer angle and magnetic obliquity are different quantities.

The chi-squared result demonstrates the inadequacy of this **specific configuration** for the measured profile. It does not disprove current-sheet emission, the split-monopole field as an outer approximation, or other parameter/emissivity choices. The next physical upgrade is to constrain the emission distribution with a plasma model and test its predicted pulse shape and energy dependence.

## Reproduction

Tested with Python 3.12. Create a virtual environment and install requirements.txt, then run:

```sh
python run_benchmark.py
PYTHONPATH=src python -m unittest discover -s tests -v
```

The run takes tens of seconds on the development machine. It writes figures, CSV data, parameters and all diagnostics to outputs/ and exits unsuccessfully if the numerical verification criteria fail. The observation mismatch is reported separately; it must never be hidden by the numerical pass status.

NASA's [particle visualization](https://svs.gsfc.nasa.gov/13058/) and [Brambilla et al. (2018)](https://arxiv.org/abs/1710.03536) remain references for a future particle-flow model. Their particle trajectories are not calculated by this milestone. See research/PULSAR-RESEARCH.md for the full development bibliography.
