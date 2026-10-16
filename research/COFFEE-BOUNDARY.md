# Controlled stellar-boundary experiment

The boundary variant repairs the magnetic-pole reversal and permits the guarded open-cap tracer to work, but it still fails numerical release. The same-boundary 65³–97³ comparison differs by 16.24% in E and 13.37% in B over 0.75<r<2 R_LC, with about 43–45% differences in the current-sheet proxy. The normalized divergence does not decrease under refinement. These are diagnostic fields, not a validated Vela solution.

The original COFFEE CPU boundary anchors B with a transition centered at `radius - 4*dx`, anchors E at `radius`, and overwrites both below `radius/2`. At the same numerical radius `0.5 R_LC`, the B transition center is only `0.130769 R_LC` on the 65³ grid and `0.252577 R_LC` on the 97³ grid. The coarse transition lies inside the hard-reset radius `0.25 R_LC`. This is a materially different source condition under grid refinement.

The original corrected 65³ run reverses the prescribed polar magnetic field inside the stated stellar surface after three rotations. Its disagreement with 97³ is about 105% in the exterior field, so this pair does not support an emission release.

`third_party/coffee-surface-anchored` is a separate experimental source copy. **This is a numerical boundary-model change, not merely an ARM portability repair.** The single difference from `coffee-portable/src` is in `boundary_pulsar`: `d1=4*dx` becomes `d1=0`. B and E then use transition centers at the same stated stellar radius. The stellar transition widths, inner overwrite, field equations, current closure, time integration, divergence cleaner, KO dissipation and outer sponge are otherwise unchanged. `boundary_alfven` is not changed.

The change is recorded in `third_party/coffee-surface-anchored.patch`; its parent portable patch, source, and binary hashes are saved in `third_party/coffee-surface-anchored-provenance.json`. The original source, binary and outputs remain available. The GPL-2.0 license is preserved. Every experimental NPZ carries `boundary_model=surface_anchored_experiment`, the distinct boundary patch hash, and `emission_release_allowed=false`.

The bounded experiment runs exactly the same 75-degree dipole, enlarged `0.5 R_LC` star, domain ±3 R_LC, KOeps .05 and three rotations. The first-rotation check at the magnetic pole gives Br=30.982 at r=.4 and Br=16.111 at r=.5, versus instantaneous dipole values 31.25 and16, with radial direction cosine >.9999. This repairs the observed interior polarity reversal at that checkpoint. The exterior remains a numerical solution to be checked: the divergence RMS divided by field RMS is still .494/R_LC in .75<r<2 at one rotation, so polarity alone does not establish usable accuracy.

Build with `./run_coffee_build_boundary.sh`, then run the separate executable without replacing the original:

```sh
/opt/anaconda3/bin/python run_coffee.py pulsar --n 65 \
  --end 18.84955592153876 --alpha 75 --extent 3 --radius .5 \
  --output-count 12 \
  --binary third_party/coffee-surface-anchored/bin/coffee-3d \
  --boundary-model surface_anchored_experiment \
  --run-dir third_party/coffee-runs/pulsar-a75-n65-3rot-surface-anchored
```

`run_coffee_boundary_audit.py` compares pole radial fields and unforced exterior r>.75 against the original same-resolution run. A same-boundary 65³/97³ comparison uses `run_coffee_assessment.py`. The numerical radius enlargement, finite domain and sponge limitations still apply. An additional retained boundary caveat is that the upstream solver initializes its dipolar sponge background before the pulsar setup changes evolving E/B to collocated111 staggering; the sponge background therefore retains the former component-coordinate offsets. The experiment deliberately does not combine a second boundary change with the source test.

## Completed coarse-grid result

The 65³ surface-anchored run completed all three rotations in 140.85 s. Its final canonical file is `third_party/coffee-runs/pulsar-a75-n65-3rot-surface-anchored/snapshot-000684.npz`. The boundary audit is `third_party/coffee-runs/surface-boundary-audit-n65.json`. At the final magnetic pole Br is 30.9819 at r=.4 and 16.1110 at r=.5, preserving the imposed polarity and approximately the expected dipole strength. The force-free constraints remain near the float32 output tolerance outside the star.

The exterior divergence metric improves from 1.1696/R_LC in the original coarse model to 0.4733/R_LC, which remains large. The normalized luminosities are 0.9066, 0.8747, 0.8675 at r=1,1.5,2, well below the ideal point-dipole reference 1.9330. Differences between rotations 2 and 3 are 0.7307% in E and 0.5357% in B, so the field is fairly periodic without being accurate enough for release. A stable numerical solution can still be inaccurate. The unforced exterior differs by order unity from the original coarse solution; the variant is not an interchangeable continuation of that model.

The independent field-line audit also improves: the guarded tracer reaches the light cylinder with outward crossing for all 128 selected open-cap seeds (both poles, 32 azimuths; seed radius 0.7, inner boundary 0.68). Its resolved cap rim spans 32.04°–86.57° in this coarse numerical model. See `outputs/field-benchmarks/surface-anchored-cap-status.json`. This verifies usable field-line topology for those seed samples; it does not remove the large divergence or resolution uncertainty.

## Completed two-grid comparison

Both runs completed three full rotations at exactly the same output phases. The 65³ run took 140.85 s; the 97³ run took 476.19 s. The final files are:

- `third_party/coffee-runs/pulsar-a75-n65-3rot-surface-anchored/snapshot-000684.npz`
- `third_party/coffee-runs/pulsar-a75-n97-3rot-surface-anchored/snapshot-001020.npz`

The complete machine-readable comparison, including final snapshot hashes and completeness checks, is `outputs/coffee-surface-assessment/results.json`. A fine-grid comparison of the original and changed source models is retained separately in `third_party/coffee-runs/surface-boundary-audit-n97.json`.

| Measured quantity | 65³ | 97³ |
|---|---:|---:|
| RMS divB / RMS B, in 1/R_LC | 0.47326 | 0.47439 |
| Same-phase E drift, rotation 2 to 3 | 0.7307% | 1.3778% |
| Same-phase B drift, rotation 2 to 3 | 0.5357% | 1.0453% |
| Poynting luminosity at r=1 | 0.90656 | 1.00446 |
| Poynting luminosity at r=1.5 | 0.87469 | 0.98429 |
| Poynting luminosity at r=2 | 0.86751 | 0.97345 |

Luminosity units are mu² Omega⁴/c³; the ideal point-dipole reference at 75° is 1.93301. Both the enlarged numerical star and boundary effects limit this comparison. The fine-grid polar Br at r=.5 is 16.01075 versus the imposed dipole value 16, confirming that the targeted polarity/source repair persists on the fine grid.

The exterior resolution difference is E 16.2381%, B 13.3716%; in the current-sheet proxy it is E 45.0097%, B 43.2119%. Quarter-rotation rigid-pattern tests reach E 9.8310% and B 7.5413% error on the fine grid. Small full-rotation drift therefore does not make the frozen rotating-pattern approximation exact.

The spatial-divergence error is essentially unchanged in physical 1/R_LC units. Multiplying it by the smaller cell width produces a smaller number, but that alone is not evidence of convergence. The controlled experiment is an actionable repair of the stellar source condition and a usable diagnostic backend; it is not a globally converged magnetosphere or an emission-release result. No further production runs were started for this bounded pass.
