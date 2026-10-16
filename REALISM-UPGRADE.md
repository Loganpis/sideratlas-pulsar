# NO.002 — revising the pulsar's physical assumptions

> Historical phase report (October 13–14). See [CURRENT-STATUS.md](CURRENT-STATUS.md) for subsequent source-boundary audits, COFFEE calculations, outer transport and the current release decision. The test count and field backend described below belong to this earlier phase.

2026-10-13–14 · @sideratlas

The project now has an evolving resistive electromagnetic-field solver, charge-based primary loading, calibrated Fermi response folding, and tested synchrotron and inverse-Compton radiation calculations. These are concrete physics upgrades. **They do not yet establish a realistic Vela reconstruction.** The global field remains resolution-sensitive; a predictive particle population and its coupling to that field remain unfinished.

Read the measured [field results](outputs/plasma-upgrade/RESULTS.md) and [emission results](outputs/emission-upgrade/RESULTS.md). All previously inspected Vela profiles are now development data, including the bands reserved in the earlier experiment. Their discrepancies informed this redesign, so another comparison with them is not independent validation.

**Verification: 74 tests passed**, covering analytic fields, trajectories and timing, radiation energy integrals, charge-loading identities, detector response, data integrity and numerical limits. [Recorded test evidence](outputs/upgrade-verification.json) and [full test log](outputs/upgrade-tests.txt) identify the tested source files. Passing these tests does not imply the global magnetosphere or Vela pulse model is correct.

## What changed

| Earlier simplification | Implemented replacement | What remains assumed |
| --- | --- | --- |
| Prescribed rotating vacuum-dipole geometry | New backend evolves Maxwell's equations with charge drift and parallel conduction, extending beyond the light cylinder | Macroscopic conductivity, enlarged star, boundary treatment and numerical regularization; not a kinetic pair calculation |
| Uniform loading per magnetic-flux tube | Alternative injection weights from the charge required by the assumed corotation electric field and the radial particle flux | Net charge does not determine electron/positron populations, pair multiplicity or actual extraction from the surface |
| Curvature radiation only | Independent classical synchrotron and isotropic inverse-Compton kernels with Klein–Nishina recoil, energy accounting and regime checks | No arbitrary mixture added to the Vela profile; pitch evolution and anisotropic seed radiation remain unsolved |
| Flat detector response within each energy band | Official LAT FRONT/BACK effective area and energy migration | Incidence scenarios replace unavailable source-specific observing history; not exact catalog weighted-count predictions |
| Artificial phase smoothing | Removed from the new diagnostic; angular-bin refinement is measured separately | Finite angular sampling still needs convergence |
| One assumed injection energy and viewing direction | Explicit injection-energy, loading, incidence and viewing-angle sensitivity calculations | Scenario spread is not a measured confidence interval |

The two calculations are deliberately distinct. **The new field backend is undergoing verification. The emission comparison still uses the frozen legacy trajectories**, allowing loading, response and smoothing changes to be isolated without claiming transport through an unvalidated field. The old inclination, acceleration prescription and phase alignment are retained in that comparison; no new parameter fit was performed.

## Field equations and interpretation

In rationalized units with $c=\Omega=R_{\rm LC}=1$,

$$\partial_t\mathbf A=-\mathbf E,\qquad \mathbf B=\nabla\times\mathbf A,\qquad
\partial_t\mathbf E=\nabla\times\mathbf B-\mathbf J,$$

$$\mathbf J=(\nabla\cdot\mathbf E)\frac{\mathbf E\times\mathbf B}{B^2+E_0^2}
             +\kappa\mathbf E_\parallel.$$

The invariant $E_0$ keeps the drift speed at or below light speed. Compatible discrete derivatives preserve $\nabla\cdot\mathbf B=0$. The current closure follows [Kalapotharakos, Harding & Kazanas (2014), equation 9](https://arxiv.org/pdf/1310.3545). Our damping parameter is **$\kappa=4\pi\sigma/\Omega$** for that paper's Gaussian-unit conductivity. The exploratory values 80 inside and 5 outside do not reproduce its production Vela calculation.

The code separately records conduction losses, the inner magnetic-dominance projection, numerical filtering, stellar forcing, and signed sponge energy exchange. A closing ledger verifies bookkeeping, not physical correctness. Independent wave and derivative tests address numerical correctness; grid studies address the global solution. The outer sponge still needs an outgoing-wave reflection test.

The numerical stellar radius is $0.35R_{\rm LC}$, while a 12 km Vela star would be about $0.00283R_{\rm LC}$. This solver therefore cannot predict the surface, polar-cap microphysics or general-relativistic effects. See the [independent field review](research/PLASMA-SOLVER-REVIEW.md).

## Charge-based loading

For the *assumed* corotation field $\mathbf E=-(\boldsymbol\Omega\times\mathbf r)\times\mathbf B$, the implemented charge has

$$4\pi\rho=-2\boldsymbol\Omega\cdot\mathbf B
 +(\boldsymbol\Omega\times\mathbf r)\cdot(\nabla\times\mathbf B),$$

in normalized $c=1$ units. The seed weights are proportional to $|\rho|\,|v_r|\,dA$, calculated on the actual injection surface with the curl term retained. This identity is checked against an independent electric-field divergence. The leading Goldreich–Julian option and the old uniform-flux option remain available for sensitivity tests. [Bai & Spitkovsky (2010), equation 3](https://arxiv.org/pdf/0910.5741)

Normalized loading does not supply an absolute current or luminosity. The legacy field and accelerating electric work remain prescribed; this charge calculation does not make them a self-consistent plasma solution.

## Radiation and observations

The [radiation note](research/RADIATION-UPGRADE.md) describes formulae, supported regimes and analytic checks. The [detector note](research/DETECTOR-UPGRADE.md) records the version-pinned calibration, units, interpolation and missing exposure information. The official response computes area-weighted photon yield; independent scale/background fits in each energy band make the Vela comparison a conditional **profile-shape diagnostic**, not an absolute spectral likelihood.

## Before an updated Reel can represent this model

The field must converge with resolution, boundaries, stellar size, time step and run duration. Particle transport can then use those fields with an explicit, justified injection and pair prescription, consistent radiation energy loss and light-travel timing. Multi-band profiles and spectra must be compared using consistent observation selections and response. An independent observational test requires data or observables that did not guide model selection.

The current release is a reproducible physics-development package with explicit failure criteria. It contains no newly validated Vela animation.

## Reproduce

From this project directory, using the existing environment:

```sh
.venv/bin/python run_plasma_upgrade.py
.venv/bin/python refine_plasma.py
.venv/bin/python check_refinement.py
.venv/bin/python run_emission_upgrade.py
.venv/bin/python write_upgrade_results.py
PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -v
```

The field runner checks configuration, final time and available generation fingerprints before reusing snapshots. The initial legacy snapshots predate source-fingerprint recording; their generating source hash is explicitly unknown, while their file hashes identify the analyzed data. Current analysis code is recorded separately. `refine_plasma.py` interpolates the 96³ snapshot onto 128³ and evolves one more rotation; it saves resumable checkpoints. `check_refinement.py` continues the 96³ field for the same extra rotation, separating coarse time evolution from a comparison at equal final times. The fine run still lacks an independent high-resolution startup. The emission runner retains the earlier frozen calibration, removes phase blur, changes loading, folds the response, and saves all scenarios and source hashes. No network access is needed after dependencies are installed.
