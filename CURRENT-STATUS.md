> Historical research status retained from October 17. Current 3D delivery: [CURRENT-3D-STATUS.md](CURRENT-3D-STATUS.md). The completed optical research edition is recorded separately in CURRENT-OPTICAL-STATUS.md.

# v1 visualization release scope — 16 October 2026

The v1 package now presents the user-approved **Vela-informed visualization**: a 32-second full-frame Reel with three title labels. Scientific explanation and limitations accompany the video in [the model notes](release/MODEL-NOTES.md) and [Instagram description](release/instagram-description.text). The image is an illustrative projection of calculated outer emission; it does not claim a quantitative fit to Vela

The measurements below are retained as a **historical diagnostic record from 15 October 2026**. They describe earlier field and emission variants and must not be substituted for the selected v1 scene or its final provenance. In particular, historical statements that no Reel exists refer to the state of the project on that date

---

# @sideratlas · NO.002 / Current scientific status

**The implementation and verification pass is complete. The Vela reconstruction and final Reel are not ready.** The revised boundary repairs reversed polar fields and allows the full field → particle → radiation calculation to run. The measured field, sampling and observational checks still do not justify a Vela claim. All outputs retain `emission_release_allowed=false`.

**121 tests pass.** The [test manifest](outputs/current-verification.json) records exact source hashes and dependency versions; [full log](outputs/current-tests.txt). These tests verify stated equations and numerical identities, not the global astrophysical model.

## What was fixed

- Native CPU initialization now assigns the field arrays correctly; superseded runs remain identified
- A separate, fully recorded boundary experiment anchors B and E at the same numerical stellar surface, repairing the reversed polar field
- Cap tracing recognizes closure on return to the launch sphere, while rejecting nulls, unresolved paths and sampled disconnected cap topology
- Particle motion extends outside the light cylinder with laboratory trajectory curvature, signed work, energy losses and propagation delays
- The Barnard scenario uses its published gain per path length without an extra velocity projection
- The sharp FIDO screening boundary now rejects crossing derivative/integration stencils instead of generating artificial curvature radiation
- Genuine 3PC spectral bins use one normalization across all energies; no manual second-peak brightening is inserted

## Measured field outcome

Both grid sizes start independently and run for three rotations. The boundary experiment changes the numerical source model; it is not a claim of exact reproduction of upstream COFFEE. The star is enlarged to 0.5 light-cylinder radii.

| Check | Original boundary | Revised boundary |
|---|---:|---:|
| 65³ → 97³ E relative L2 difference, 0.75 < r < 2 | 104.81% | 16.24% |
| 65³ → 97³ B relative L2 difference, 0.75 < r < 2 | 104.77% | 13.37% |
| 65³ B change between rotations 2 and 3 | 48.28% | 0.54% |
| 65³ RMS divergence / RMS B, in R_LC⁻¹ | 1.1696 | 0.4733 |
| 97³ B change between rotations 2 and 3 | 4.67% | 1.05% |
| 97³ RMS divergence / RMS B, in R_LC⁻¹ | 0.0381 | 0.4744 |

A nearly repeating solution can still have substantial spatial error. In the revised model, the physical divergence metric does not decrease under refinement (0.4733 → 0.4744 per R_LC), and the current-sheet proxy still differs by 43–45%. This needs a divergence-preserving source/boundary treatment; more sampling alone will not repair it. Two grids cannot establish asymptotic convergence. Source-size, outer-domain and oblique-boundary checks remain incomplete.

![Measured field assessment](outputs/current-assessment/field-assessment.png)

Full evidence: [original grid comparison](outputs/coffee-assessment/results.json), [revised comparison](outputs/coffee-surface-assessment/results.json), [boundary experiment](research/COFFEE-BOUNDARY.md), [cap/source audit](outputs/field-benchmarks/CAP-AUDIT.md).

## Completed particle and observation comparisons

The development scenarios use inclination 75°, observer 65°, finite-radius injection at 0.7 R_LC, prescribed magnetic-flux loading and Barnard scenario-2 rates. None predicts pair creation or an absolute population.

| Calculation | Injected / observer-selected seeds | Energy beyond cylindrical light cylinder | Spectral χ² (6 descriptive dof) |
|---|---:|---:|---:|
| [n65-surface-anchored](outputs/outer-coffee-n65-surface-anchored/results.json) | 128 / 22 | 79.53% | 23064.4 |
| [n65-refined](outputs/outer-coffee-n65-refined/results.json) | 384 / 66 | 79.76% | 25730.7 |
| [n97-surface-anchored](outputs/outer-coffee-n97-surface-anchored/results.json) | 384 / 71 | 80.18% | 11114.5 |

The spectral score compares seven detections with one fitted normalization; the eighth bin is an upper limit. It is descriptive, not a calibrated joint likelihood. Pulse fits use a shared phase and separate scale/background per band with a flat detector response. All inspected data are development data.

On the same 65³ field, changing from 128 to 384 seeds and halving the orbit step changes the four pulse profiles by **57.40%, 64.27%, 77.38%, 84.29%** in relative L1 at 100 phase bins, without re-aligning them. This combined refinement does not establish pulse convergence.

The independent [emission audit](outputs/current-assessment/emission-audit.md) reproduces the physical scale, energy ledger, phase projection and spectral normalization. The coarse diagnostic has too few high-energy photons and substantial unresolved curvature. No fit has been used to compensate for those field errors.

[View the higher-resolution data comparison](outputs/outer-coffee-n97-surface-anchored/comparison.png) and [refined coarse-grid comparison](outputs/outer-coffee-n65-refined/comparison.png). Photons are not rendered as visible-light colors here; these are numerical and gamma-ray diagnostics.

## What still prevents a final Vela Reel

1. Establish a reliable global field: control divergence and numerical stellar/sponge boundaries, then demonstrate resolution, domain and stellar-size sensitivity
2. Establish converged observer sampling and resolved trajectory curvature before tuning emission physics
3. Justify particle populations and acceleration, then compare spectra and pulses with consistent observation selections/response and an independent validation test

The optional FIDO path cannot yet cross its sharp interface. A physically stated transition and a sensitivity study are required; choosing an arbitrary smoothing width merely to finish would add an unsupported assumption. The adopted moment of inertia, force-free torque normalization, enlarged star, fixed geometry, injection annulus, initial energy and acceleration remain model assumptions. Synchrotron/Compton kernels alone do not solve a pair cascade.

The source spectral FITS has invalid embedded checksum cards. Independent official downloads had identical bytes; the original data and warning are preserved in [spectral provenance](data/spectra/provenance.json).

## Reproduce and continue

Use Python 3.12 with NumPy, SciPy ≥ 1.12, Matplotlib, Astropy and h5py. Exact tested versions are in the manifest. The retained local interpreter was `/opt/anaconda3/bin/python`.

```sh
PYTHONPATH=src /opt/anaconda3/bin/python -m unittest discover -s tests -v
/opt/anaconda3/bin/python run_outer_diagnostic.py \
  third_party/coffee-runs/pulsar-a75-n97-3rot-surface-anchored/snapshot-001020.npz \
  --inclination 75 --seed-radius .7 --seed-inner-radius .68 \
  --azimuths 64 --rings 3 --dt .004 --output outputs/new-diagnostic
```

The separate native builds, source patches, licenses and exact field-run commands are recorded in [COFFEE](research/COFFEE.md) and [the boundary experiment](research/COFFEE-BOUNDARY.md). The earlier black-hole release is separate.
