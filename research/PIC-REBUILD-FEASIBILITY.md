# Rebuilding the pulsar bundle with real particle-in-cell calculations

Research date: 2026-10-17

## Finding

The existing rendered magnetosphere is not a self-consistent particle-in-cell (PIC) solution. Replacing its graphics with more numerous moving dots would not change that. The original `src/pulsar/plasma.py` says explicitly that it is a resistive Maxwell prototype; it evolves fields using an assumed current closure, not particle charge/current deposited back onto the mesh. The later COFFEE field solution plus test-particle emission is also not PIC. Neither supplies a calculated Vela atmosphere or a visible-light image of Vela.

There is a practical open-source route to real PIC on this computer: **Entity with its Serial CPU backend**. Its existing PIC kernels can be built and checked locally. However, its public aligned, axisymmetric magnetosphere problem is not an oblique three-dimensional Vela production model. A new problem setup, physical scaling, resolution study and synthetic optical-emission calculation are required before claiming the requested result.

No approved production PIC dataset for Vela has been located in the NASA references. This is a missing scientific input, not a license to substitute an illustration and label it a reconstruction.

## What NASA actually provides

- [NASA SVS 13058](https://svs.gsfc.nasa.gov/13058/) describes two-way coupling of particles and electromagnetic fields, particle injection from the star and simulations on Discover and Pleiades supercomputers. Its electron/positron track colors encode species and energy. They are not camera-visible emitted colors.
- The underlying [SVS 4648 all-particle-flow page](https://svs.gsfc.nasa.gov/4648/) offers rendered movies, rendered frame sets and color scales. Its dataset entry identifies the **Brambilla Pulsar Model, ID 1013**, and explicitly states that SVS does **not** store the datasets themselves or further dataset details. Thus the frame download is not a raw particle/field download.
- The linked study is Brambilla et al., *Electron–Positron Pair Flow and Current Composition in the Pulsar Magnetosphere*, ApJ 858, 81 (2018), [DOI 10.3847/1538-4357/aab3e1](https://doi.org/10.3847/1538-4357/aab3e1), [arXiv:1710.03536](https://arxiv.org/abs/1710.03536). It investigates surface pair injection and particle/current flow; the public animation is not an observer-specific optical image of Vela.
- The [2008 NASA Fermi discovery illustration](https://www.nasa.gov/universe/nasas-fermi-telescope-discovers-first-gamma-ray-only-pulsar/) is a conceptual magnetic-field-and-lighthouse view. Matching its appearance would not validate an optical radiation model.

The requested visible-light camera and the NASA particle diagnostic are distinct output quantities. An optical camera needs directional spectral emission and radiative transfer; particle paths alone do not define image brightness.

## Source-code inspection

The following sources were read from pinned local checkouts. They were downloaded for research; they are not automatically part of the public release.

| Code | Source and pin | License | Relevant capability | Limitation for this task |
|---|---|---|---|---|
| Entity | [entity-toolkit/entity](https://github.com/entity-toolkit/entity), `62ab34fa0a0d1d6965918abb16d1d0c4b4740cf3` | BSD-3-Clause | Real relativistic PIC, charge/current deposition, field evolution, Boris and GCA pushers; Serial/OpenMP/GPU portability through Kokkos; 1D/2D/3D Cartesian streaming problems | `pgens/magnetosphere/pgen.hpp` explicitly supports only 2D spherical/QSpherical coordinates and an aligned dipole/monopole. Its supplied atmosphere is a numerical pair-injection boundary, not Vela's measured radiating atmosphere |
| Aperture4 | [fizban007/Aperture4](https://github.com/fizban007/Aperture4), `59a975f2202202b9442fbcbce66709541539f989` | GPL-3.0 | PIC with particle, field and radiation modules; core can build on CPU and CUDA/HIP GPUs | Its public `problems/pulsar/src/main.cu` uses `Config<2>` spherical geometry and CUDA classes/device lambdas. The pulsar target is CUDA-gated and commented out of default problem selection. It is not a ready 3D CPU pulsar executable |
| Tristan-MP v2 | [PrincetonUniversity/tristan-mp-v2](https://github.com/PrincetonUniversity/tristan-mp-v2), public v2.9 inspected through primary README/LICENSE | BSD-3-Clause | Established relativistic PIC with MPI and particle/radiation facilities; Fortran | Requires MPI Fortran toolchain and task-specific pulsar user setup; public README does not establish a ready reproducible oblique Vela atmosphere/radiation run |
| Zeltron | [official code page](https://ipag.osug.fr/~ceruttbe/Zeltron/index.html), [download and license](https://ipag.osug.fr/~ceruttbe/Zeltron/download.html), [public source](https://github.com/bcerutti/Zeltron-code) | GPL-3.0-or-later per official download page | Explicit relativistic electromagnetic PIC, radiation reaction; public 2D/3D Cartesian release; official page lists OS X support | The public download is not the authors' ready pulsar production setup. Magnetospheric boundary/injection/scale choices and optical transport still need implementation and validation |

A code being capable of PIC is necessary but is not evidence that a particular pulsar run is accurate. Licensing must be retained with any redistributed source. NASA's name and source credits do not imply NASA endorsement or that our implementation reproduces its results.

## Local benchmark route

Entity was selected for a minimal local foundation check because it can run on the available Apple CPU without a CUDA device and without installing a global MPI stack. The local dependency is Kokkos 5.2.1, commit `89f0b12e43583d02bd19739810ceadcedf1ccd56`, under its included Apache-2.0-with-LLVM-exception license.

Configuration uses AppleClang 21, C++20, double precision, Serial backend, no MPI, no ADIOS2 output. The `streaming` problem and upstream pusher/deposition/Maxwell tests are built. Passing those checks is a toolchain and algorithm regression result, not a validation of pulsar production physics.

The local 3D smoke input is `research/pic-benchmark/periodic-pairs-3d.toml`, adapted from Entity's public Weibel setup: 16³ periodic cells, four electron/positron species, four reference particles per cell, counterstreaming drift four-velocities ±0.3, dimensionless temperature 10⁻⁴, duration 8 code units. This deliberately small run does not resolve a pulsar atmosphere, skin-depth hierarchy, current sheet or Vela emission. No result from it may be used as the final Reel.

Measured result: the unmodified Entity source compiled, all four upstream kernel checks passed, and the coupled 3D run completed successfully. There were **16,384 macroparticles (4,096 per species), 55 timesteps, Δt = 0.1443375673**, and a final simulated time of **7.9385662014** code units (the requested runtime is rounded down to a whole number of steps). The recorded rerun took about **0.394 seconds** on the Serial CPU backend. Particle counts remained fixed in the periodic box. The run log confirms active field solver, current deposition and Boris-pusher stages, with radiation drag and emission policies absent.

The four upstream tests are `KERNELS::pusher`, `KERNELS::deposit`, `KERNELS::faraday_mink` and `KERNELS::ampere_mink`. The pusher test compares two particles' three-dimensional momentum components and Lorentz factors with analytic magnetic gyromotion using the upstream double-precision tolerance of 10⁻³. The other checks exercise their upstream operator/deposition assertions. All four returned exit code zero, both when run directly and through CTest. These checks do not establish collective instability growth, global energy conservation, convergence or astrophysical validity for the 3D smoke run.

Build and run records are under `outputs/entity-pic-*`; measured status and hashes are in `outputs/entity-pic-benchmark/verification.json`. `research/pic-benchmark/run_local_benchmark.py` checks the exact source pins and unmodified tracked sources, reruns the upstream checks and smoke input, and writes a machine-readable record. It neither downloads sources nor claims pulsar validation.

### Reproduction configuration

Use the pinned Entity and Kokkos source checkouts named above. Configure with CMake and the following options, pointing the dependency source directory at the local pinned Kokkos checkout:

```text
-Doutput=OFF -Dmpi=OFF -Dpgen=streaming -Dprecision=double
-DTESTS=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF
-DKokkos_ARCH_ARMV8=ON -DKokkos_VERSION=5.2.1
-DFETCHCONTENT_SOURCE_DIR_KOKKOS=<absolute path to pinned Kokkos>
-DCONNECTION_CHECKED=TRUE -DFETCHCONTENT_FULLY_DISCONNECTED=OFF
-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_C_COMPILER=/usr/bin/clang
```

`Kokkos_VERSION` is supplied explicitly because this Entity commit otherwise tries to query a hard-coded `_deps/kokkos-src` directory when `FETCHCONTENT_SOURCE_DIR_KOKKOS` points elsewhere. No code changes were needed. The dependency source override prevents any network fetch during this configuration.

Build targets `entity.xc`, `test-kernels-pusher.xc`, `test-kernels-deposit.xc`, `test-kernels-faraday_mink.xc` and `test-kernels-ampere_mink.xc`. Then run `research/pic-benchmark/run_local_benchmark.py` with Python 3. The retained CMake configuration log and build log are `outputs/entity-pic-configure.log` and `outputs/entity-pic-build.log`; the CTest log is `outputs/entity-pic-tests.log`.

## Requirements before a scientifically defensible production release

1. Obtain a verified 3D particle/field dataset with particle IDs, charges, weights, momenta, times, mesh geometry, boundaries and physical normalization, or generate one with a documented pulsar problem generator
2. Establish kinetic scales and show grid, particle-number, timestep, domain-size and runtime convergence for the target diagnostics; report unresolved ratios explicitly
3. Distinguish the injected numerical pair atmosphere from a physical stellar atmospheric radiation model
4. Specify Vela priors with provenance, including period/epoch, geometry, mass/radius and atmosphere assumptions; preserve uncertainties rather than treating fitted or conventional values as measurements
5. Compute observer-specific optical emission from particle distributions, pitch angles, radiation reaction and absorption, and combine a physically defined atmosphere spectrum with gravitational redshift/light bending as appropriate
6. Compare independent synthetic observables with Vela optical flux, pulse shape and polarization; use only relevant measured quantities as validation
7. Render the validated image into the full 1080×1920 frame without changing the solved star/magnetosphere scale ratio or inventing visible field lines, colored particle tracks, beams or thermal surface structure

Until these requirements are met, a complete scientifically rigorous Vela visible-light reconstruction is not available for publication. The caption must reflect the actual achieved level, not the desired one.
