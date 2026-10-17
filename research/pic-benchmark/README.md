# Reproduce the local PIC foundation check

This is a small, real coupled particle-in-cell (PIC) calculation and four upstream algorithm checks. It is **not a pulsar run, Vela dataset, visible-light image or finished Reel**

Requirements: Python 3, Git, CMake ≥3.22, a C++20 compiler and a native build tool. The recorded run used AppleClang 21 on macOS ARM64, Kokkos Serial CPU, double precision, no MPI and no ADIOS2. No GPU is required for this small check

From the repository root, inspect the plan without downloading or building:

```sh
python research/pic-benchmark/build_local_benchmark.py
```

Fetch the exact public commits and build only the benchmark executable and four upstream tests:

```sh
python research/pic-benchmark/build_local_benchmark.py --fetch --build --jobs 2 --cmake /absolute/path/to/cmake
```

If `cmake` is already on PATH, omit `--cmake`. In the original local workspace the bundled executable is `third_party/coffee-tools/env/bin/cmake`; that toolchain is not redistributed in the publication repository

Run the checks and 3D smoke test:

```sh
python research/pic-benchmark/run_local_benchmark.py
```

The runner records `outputs/entity-pic-benchmark/verification.json`, individual test logs and `run.log`. The builder keeps a uniquely named build log alongside those files. It does not overwrite source checkouts, change non-pinned checkouts, reset local changes or delete directories. Repeated builds may update the existing dedicated build directory. Missing source directories are fetched only with `--fetch`; test execution is a separate explicit command

## Pins and numerical scope

- Entity: [entity-toolkit/entity](https://github.com/entity-toolkit/entity), commit `62ab34fa0a0d1d6965918abb16d1d0c4b4740cf3`, BSD-3-Clause
- Kokkos 5.2.1: [kokkos/kokkos](https://github.com/kokkos/kokkos), commit `89f0b12e43583d02bd19739810ceadcedf1ccd56`, Apache-2.0 WITH LLVM-exception

The original upstream license files remain in their downloaded checkouts. The builder fetches exact commit hashes over HTTPS and never disables certificate verification. It skips optional Git LFS image downloads. CMake receives the pinned local Kokkos directory, so it does not need to download a moving branch or dependency tag

The 3D smoke input adapts Entity's supplied neutral counterstreaming-pair problem to a **16³ periodic box, four species and 16,384 total macroparticles**, with a requested runtime of 8 code units. It completes **55 timesteps**, with Δt = 0.1443375673, ending at 7.9385662014 code units. It uses actual current deposition, field evolution and Boris particle updates. Radiation drag and emission are absent

The four upstream checks cover analytic 3D magnetic gyromotion, current deposition, Faraday and Ampere operators. Their runtime checks use explicit throws or Entity's `raise::ErrorIf`/`raise::Error` exception path; they are not C `assert` checks disabled by `NDEBUG`. The pusher check uses the upstream 10⁻³ double-precision tolerance. All four passed in the recorded run

These results establish local execution of the PIC foundation. They do not measure collective-instability growth, global energy conservation, kinetic convergence, pulsar boundaries, pair cascades or Vela's atmosphere/radiation. See [the feasibility report](../PIC-REBUILD-FEASIBILITY.md) for the production requirements still outstanding
