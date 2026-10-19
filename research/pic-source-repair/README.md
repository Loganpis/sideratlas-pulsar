# Compatible rotating conductor experiment

This is an isolated repair of the reduced local 3D PIC source. It is **not** an optical Vela prediction. Upstream Entity/Kokkos and `research/pic-local` are unchanged. The authored PGen, source, configuration and new binary hashes are recorded per run

## Source equation

The finite regularized seed vector potential is rotated by its prescribed inclination and integrated spin phase. Let `A_rot(t)` be that analytic potential evaluated on the electric Yee edges, and `Phi(t) = (Omega(t) × r) · A_rot(t)` on vertices. In the maintained core the source electric field at the integer electric-field time is

`E_edge(t) = -[A_rot(t+dt/2) - A_rot(t-dt/2)]/dt - grad_h Phi(t)`

There is no magnetic-field reset, fit correction, or poststep magnetic impulse. Entity's two Faraday half-updates use the same saved E(t). Their sum advances B from t-dt/2 to t+dt/2. The compatible curl of the vertex gradient vanishes algebraically, so the native Faraday update advances a maintained core exactly through the intended discrete rotating-potential family

In the continuum, a rigidly rotated vector potential satisfies `dA_rot/dt = v_rot × B - grad(v_rot · A_rot)`. Thus the above expression tends to the ideal conducting-star electric field `E = -v_rot × B`. The independent algebra check measures second-order convergence of E toward this continuous condition; it does not validate the kinetic magnetosphere

The compatible source is blended with the original interpolated corotation E in a compact region **inside the numerical star**. The weight is 1 inside `anchor_core_fraction * Rstar`, a cubic smoothstep to 0 at Rstar, and zero outside Rstar. Only electric edges inside the star are imposed. All magnetic components are always evolved by native Faraday; all external electric fields retain the native Ampere/current updates. Particle current feedback is enabled. This is a prescribed maintained stellar magnetization, not a simulation of the neutron-star interior current generation. The imposed E represents an interior source current/charge that balances the native Ampere update; the exported J channels still contain only the deposited plasma current. The existing conductor electric-energy change is recorded, while a complete boundary/source energy balance remains unimplemented

For a face with radius less than `anchor_core_fraction * Rstar - dx/2`, every electric edge of its native curl belongs to the maintained core. The audit fits only these guaranteed faces. Surface/exterior fits are reported separately and may remain affected by staircasing, the transition layer, plasma loading and outer boundary behavior. An exact core is not evidence that those physical/numerical errors have converged

## Independent algebra check

Run `python research/pic-source-repair/check_compatible_source.py`. It independently evaluates the edge and vertex formulas. Five checks at startup and rotating times reproduce B updates to better than 1e-14. The conductor E error decreases by approximately four on each doubling of resolution on interior points. This checks sign, units, gradient/curl compatibility and the temporal staggering

## Native run

From the repository root:

```
python research/pic-source-repair/run_source_repair.py --build --cmake /absolute/path/to/cmake --jobs 2 --run-dir outputs/my-source-repair
```

The default route performs current clean-source/pin scans and builds a separate executable under `outputs/pic-source-repair-build`. No old build executable is overwritten. To avoid a slow repeated filesystem worktree scan, `--upstream-receipt outputs/local-pic-startup-v3/verification.json` can instead document a prior direct successful pinned build and check current HEAD pins. This optional route explicitly records that upstream worktree bytes were not freshly rescanned. The external PGen is still freshly built and hashed. CMake can also query source Git metadata during configuration. `--kokkos-package /path/to/existing/KokkosConfig-directory` imports an already built dependency instead of rebuilding it; the runner records its static-library hashes. `--build-directory outputs/a-separate-build` selects a distinct compiler workspace

```
python research/pic-source-repair/audit_source.py outputs/my-source-repair/native --output outputs/my-source-repair/source-audit.json
python research/pic-local/audit_native.py outputs/my-source-repair/native --output outputs/my-source-repair/independent-audit.json
python research/pic-local/audit_charge.py outputs/my-source-repair/native --require-match --output outputs/my-source-repair/independent-charge-audit.json
```

The source audit explicitly compares at B time `t-dt/2`. Newborn momentum, source supply, enlarged star, boundary damping, reduced gyro/plasma scales and all other limitations remain those of the original prototype. No radiation reaction, cascade, physical pair-energy/pitch mapping, surface atmosphere, optical transfer or steady-state/convergence claim is added

The initial E is initialized using the compatible conductor too. Hence this initial E is not bit-identical to the older prototype, although the initial B and particle supply remain the same

## Relation to literature

The prescribed rotating magnetized conducting-star boundary is a standard global pulsar-model assumption, e.g. Philippov & Spitkovsky (2014), https://doi.org/10.1088/2041-8205/785/2/L33, and Belyaev (2015), https://academic.oup.com/mnras/article/449/3/2759/1132376. The specific compact compatible edge-potential implementation here is original; these references are not asserted to validate it or to use this exact discretization. Native Yee curl compatibility is what preserves discrete magnetic divergence; measured source and kinetic validation are separate requirements

## Bounded native results, 18 October 2026

The main successful run is `outputs/pic-source-repair-imported-one-turn-v2`. Its fresh PGen/Entity compile imports the previously built exact-pinned Kokkos static libraries; their byte hashes and the historical direct-build source receipt are in `verification.json`. Earlier interrupted/failed build attempts are preserved and are not counted as successful runs

The aggregate scalar report is `outputs/pic-source-repair-report.json`. Raw field arrays, particle snapshots and event logs remain local. Keep only source/configuration, verification receipts, ledger, metadata, small audit JSON and their hashes in a compact publication package

| Measurement | Original 48³ | Compatible E, 48³, core .8R | Compatible E, 96³, core .8R |
|---|---:|---:|---:|
| Steps | 317 | 317 | 634 |
| Final position/E time | 6.67260545 | 6.67260545 | 6.67260545 |
| Stellar radius in cells | 4.8 | 4.8 | 9.6 |
| Max original 48³ interior-fit magnitude drift | 7.94% | 0.434% | Not a same-mask comparison |
| Max original 48³ interior-fit phase error | 11.04° | 0.222° | Not a same-mask comparison |
| Guaranteed core magnitude drift | Not maintained | <3e−14 fractional | <5e−14 fractional |
| Final surface radial-fit amplitude ratio | 1.0274 | 0.9204 | 0.9707 |
| Final surface radial-fit phase error | −15.92° | −2.23° | −2.39° |
| Final surface radial-fit L2 residual | 8.82% | 39.19% | 19.87% |

The target of <1% interior moment drift and <1° interior phase error is met in the same 48³ interior mask. The guaranteed core result is expected from the imposed source equation and checks the numerical boundary implementation. It is not an inference of Vela magnetization

The 96³ run takes 166.6s on this local CPU setup. Its 256 pairs/step, half timestep and 1/8 macro charge factor preserve the normalized injected species number flux `q0*pairs/dt` exactly. The final position/E times match to 3.4e−14. It has 180,108 live particles, gamma_max≈3.09. The initial dt-proportional batch, injection shell location `R+.7dx`, sampling, staircasing and kinetic resolution also change; this is one controlled refinement, not a complete convergence study

All 17 main-run and 5 refined-run charge reconstructions match the separately exported CIC vertex charge, and all particle inventories match birth/death events. Final `rms(divB)/rms(B)` is about 2.26e−14/RLC in the 48³ diagnostic shell and 6.90e−14/RLC in the 96³ shell. Their symmetric Gauss residuals are .00125 and .00155 in resolution-dependent declared masks. These discrete constraints are not a physical accuracy certificate and do not establish a converged charge solution

The current-deposition control preserves initial fields/particles exactly. At the final matched step, disabling deposited plasma current changes exterior-shell E and B L2 norms by 55.2% and 47.1%, and leaves deposited current zero. Only the output snapshot interval is additionally reduced for storage. Plasma feedback remains active in the main run

The `.95R` core-radius sensitivity test worsens the surface dipole fit residual to 70.5%, despite maintaining the interior exactly. It is not promoted as an accuracy fix. The `.8R` no-current control has surface radial amplitude ratio .9169 and phase error −.288°, compared with .9204 and −2.23° when plasma current is coupled

**Remaining limitation:** refining improves the surface amplitude deficit and dipole residual, but its phase error does not converge in these two runs. The source interior is repaired; the physical stellar surface boundary, outer boundary effects, relaxation, kinetic population and optical scaling remain unvalidated. No optical Reel has been generated by this experiment
