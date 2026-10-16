# COFFEE CPU backend: provenance, verification and limits

This project now builds and runs the native CPU, three-dimensional EZ field solver from [COFFEE](https://github.com/fizban007/CoffeeGPU), pinned to commit `9986522d08b4d887e5ef9a53e12cecae4e376d22`. It is a force-free field backend. It does not supply a self-consistent dissipative electric field, particle population, pair cascade, or Vela reconstruction. Every canonical snapshot exports `emission_release_allowed=false`.

## Source and licensing

- `third_party/coffee-src` is the unmodified upstream checkout, including its GPL-2.0 `LICENSE` and dependency licenses.
- `third_party/coffee-portable` is the buildable copy; `third_party/coffee-portable.patch` contains every source change relative to the pinned upstream source, including the independent verification driver.
- `third_party/coffee-provenance.json` records the compiler, source file hashes, patch hash and executable hashes. Each run records its executable hash and exact configuration. `third_party/coffee-dependencies.json` records downloaded dependency versions, URLs and available SHA256 hashes.
- The imported and modified COFFEE code retains GPL-2.0. Keep the license and corresponding source/patch with redistributed COFFEE-derived binaries. The numerical outputs are not a relicensing of the solver.

The native ARM compatibility changes replace the unsupported x86 VectorClass API with a one-lane wrapper, omit an empty CPU spherical implementation, remove unused OpenMP runtime includes and the obsolete `stdc++fs` link, and supply dependency locations explicitly because the MPI wrapper mishandles spaces in its installation path. The wrapper uses upstream equations, difference stencils, time integration, force-free constraints and boundary formulas unchanged. It does not claim bitwise equivalence with an AVX/FMA build.

A separate upstream CPU correctness bug was found and repaired: `multi_array::assign_dev()` was empty, although `vector_field::initialize()` and `assign()` call it on CPU. It now assigns host memory. Diagnostic scalar arrays are explicitly initialized to zero. Runs made before this fix are marked `SUPERSEDED.md` and must not be used as verification or physical field inputs.

## Build and reproduce

The completed build uses AppleClang 21, C++17, native arm64, double-precision evolution, Release optimization, OpenMPI 5.0.11 and parallel HDF5 2.2.0. No CUDA, cloud resources or global package installation is used. Project-local dependency tools are in `third_party/coffee-tools/env`; all downloaded packages retain their license metadata.

If that local environment is missing, install it into the project (requires network access):

```sh
mkdir -p third_party/coffee-tools
curl -L --fail https://micro.mamba.pm/api/micromamba/osx-arm64/latest \
  -o third_party/coffee-tools/micromamba.tar.bz2
tar -xjf third_party/coffee-tools/micromamba.tar.bz2 \
  -C third_party/coffee-tools bin/micromamba
third_party/coffee-tools/bin/micromamba create -y \
  -r third_party/coffee-tools/mamba-root -p third_party/coffee-tools/env \
  -c conda-forge --override-channels cmake openmpi 'hdf5=*=mpi_openmpi*'
```

The exact installed package identities are saved in `third_party/coffee-dependencies.json`; the command above is a convenient fresh solve, not a version lock. The retained environment is the version-pinned build used for the reported runs.

```sh
./run_coffee_build.sh
python run_coffee_verification.py --directory third_party/coffee-runs/verification-new
python run_coffee.py pulsar --n 65 --end 18.84955592153876 \
  --alpha 75 --extent 3 --radius .5 --output-count 12 \
  --run-dir third_party/coffee-runs/pulsar-new-n65
python run_coffee.py pulsar --n 97 --end 18.84955592153876 \
  --alpha 75 --extent 3 --radius .5 --output-count 12 \
  --run-dir third_party/coffee-runs/pulsar-new-n97
python run_coffee_assessment.py \
  --coarse third_party/coffee-runs/pulsar-new-n65 \
  --fine third_party/coffee-runs/pulsar-new-n97
```

Python requires NumPy, SciPy and h5py. The actual runs used the already-installed `/Users/Logan/anaconda3/bin/python`; the project's `.venv` did not contain h5py. The build script uses the already installed macOS command-line SDK and deliberately fails if its project-local build dependencies are missing. Existing output directories are never silently overwritten.

## Exact numerical setup and import contract

The pulsar runs set `problem=1`, `omega=1`, `b0=1`, `(p1,p2,p3)=(sin75°,0,cos75°)`, zero quadrupole, radius `0.5`, and cubic domain `[-3,3)` in each direction. **The radius is an enlarged numerical star, not Vela's physical radius.** Odd grids avoid placing the singular dipole origin on a node. The step is at most `0.3 dx` and adjusted to put each output exactly a quarter rotation apart. Each run covers `Omega t=6pi`, three rotations, starting from a static dipole and zero electric field. A three-rotation runtime alone does not establish relaxation; same-phase and rotated-pattern comparisons are mandatory.

The CPU solver uses fourth-order central spatial differences, five-stage low-storage RK4, sixth-order Kreiss–Oliger dissipation (`KOeps=.05`), divergence cleaning, and force-free current with explicit `E.B` removal and `E<=B` clipping. `calc_current`, `clean_ep` and `check_egb` configuration toggles are not obeyed by this upstream CPU implementation: the current and combined cleaning run unconditionally. A vacuum test cannot be selected by turning these flags off.

The sponge width is approximately `0.75 R_LC` at each box face, starting at about `|x|=2.25`, with a cubic profile and `sigpml=.5`. It damps E toward zero and B toward the initial dipole once per full time step. It is not an outgoing-characteristic boundary or a demonstrated perfectly matched layer, despite the upstream `pml` parameter names. Its damping strength is per step, so the effective physical damping changes under resolution at fixed `sigpml`. Nonperiodic outer guard cells retain initial values. Domain-size and oblique-outflow boundary sensitivity remain unmeasured. The scalar divergence-cleaning gradient also has the upstream fixed radial taper near coordinate radius 3; it is not rescaled by `omega` in the code.

All physical E and B components in the 3D EZ setup are collocated at `Stagger(111)` nodes. Their coordinates are `lower + i*size/N`, without a half-cell shift. The generic library supports other staggerings; the importer does not assume those are collocated. The standard `fld.NNNNN.h5` outputs store `Ex,Ey,Ez,Bx,By,Bz` as float32 with NumPy shape `(Nz,Ny,Nx)` even when evolution uses double precision. With `downsample=1`, the 111-to-111 export does not average E/B values. The canonical NPZ transposes each component to `(Nx,Ny,Nz)` and stacks them as `(3,Nx,Ny,Nz)`.

The standard field HDF5 has no time dataset. Its exact output time is `file_index * data_interval * dt`, before that step's update. The main upstream loop performs one extra update after the last saved snapshot; this does not change the saved time. NPZ metadata records this derivation, the exact TOML configuration, source and patch hashes, HDF5 hash, inclination, `star_radius` in R_LC, and `magnetic_phase_rad=omega*t`. The phase-zero dipole lies in the +x,+z plane and rotates positively about +z.

For zero quadrupole, upstream `B=b0 [3 rhat(p.rhat)-p]/r^3`; there is **no** extra stellar-radius-cubed factor. The dipole moment is `mu_code=b0*|p|`. Conversion multiplies coordinates/time by `omega` and divides E/B by `mu_code*omega^3`. Thus the canonical fields have unit `mu/R_LC^3` normalization, directly compatible with an explicit physical B_LC multiplier. For the reported runs, `omega=b0=|p|=1`, so no numerical rescaling occurs. E and B use a common Gaussian `c=1` scale. No hidden parallel electric field is inserted.

## Completed elementary benchmarks

Results are in `third_party/coffee-runs/verification-cpu-init-fixed/verification.json`, each linked to its actual executable/configuration/log. The exact periodic force-free solution is `E=(0,f,0)`, `B=(1,0,f)`, `f=.2 sin[2pi(x-t)]`, evolved for one period. It has zero current, E.B=0 and B²−E²=1. Small periodic transverse dimensions still exercise the actual 3D field update and MPI guard exchange.

| x resolution | relative L2 wave error | measured order |
|---|---:|---:|
| 16 | 4.8939703e-3 | — |
| 32 | 3.1012046e-4 | 3.9801 |
| 64 | 1.9446696e-5 | 3.9952 |

For a right-moving Gaussian guide-field pulse crossing the cubic sponge at normal incidence, the maximum returned left-characteristic energy fraction was 0 at Nx96 and 2.00e-30 at Nx192 (numerical zero for this symmetry). Remaining interior wave energy at t=7.5 was 4.98e-10 and 1.75e-19 of the initial incident energy. This test passes the specified 0.1% reflected-energy bound, **but its exact planar symmetry is particularly favorable**: it does not certify oblique waves, a dipolar damping background, a current sheet, or a pulsar wind.

The wave verifies the compiled update and expected order for a smooth exact force-free solution. These elementary checks do not validate the global pulsar field or authorize an emission release.

## Pulsar convergence and release status

The completed original-boundary comparison is in `outputs/coffee-assessment/results.json`. It compares the 65³ run ending `n65-3rot-cpu-fixed` with the complete 97³ run ending `n97-3rot-complete`, at exactly three rotations. The earlier 97³ directory ending `cpu-fixed` was interrupted and is not the final comparison. The analysis linearly interpolates the fine grid to coarse nodes; measures shells from 0.75 to 2 R_LC and the top decile of |curl B|/|B| in the outer shell as a current-sheet proxy; measures divergence and force-free constraints; compares the second and third rotation; and rotates quarter-period snapshots against the same initial field pattern. Surface Poynting fluxes are compared with the ideal point-dipole reference `1+sin²75°`, with the enlarged star and finite box stated explicitly.

The original 65³ versus 97³ fields differ by 104.81% in E and 104.77% in B in the combined shell. At 65³ the polar B reverses near the nominal stellar surface and some interpolated magnetic lines approach artificial attracting nulls. The initial grid-node field agrees with the analytic dipole to float32 export precision; this is not an axis or export sign error. The smooth B forcing was centered at `Rstar-4dx`, while E forcing was centered at `Rstar`. At 65³ that B center is 0.1308 R_LC, inside the separate hard-reset radius 0.25; immediately outside the hard-reset region, smooth B anchoring is almost absent. These diagnoses are recorded in `outputs/field-benchmarks/native-source-audit.json` and `native-null-audit.json`.

A **separate numerical boundary-model experiment**, `third_party/coffee-surface-anchored`, changes only `d1=4dx` to `d1=0` in `boundary_pulsar`, placing the smooth B and E forcing centers at the same nominal stellar surface. This is not a pure portability change or an exact reproduction of upstream COFFEE. The original executable and runs are preserved. Distinct patch, binary and snapshot provenance identify the experiment. Its 65³ three-rotation run repairs polar polarity and reduces the two-to-three-rotation drift below 1%, but has substantial remaining divergence. The current status records the full comparison and release decision.

An additional outer-boundary limitation remains: the solver constructs its background dipole before the pulsar initialization changes the evolving E/B staggering to 111. Consequently the sponge's background components use slightly shifted coordinates. The exported/evolving fields are collocated; the sponge target mismatch is distinct and has not been silently fixed in the source-boundary experiment.

The release gate remains closed until field resolution, relaxation, boundary/domain sensitivity, and physical model limitations are assessed together. No radiation fit or Vela match can retroactively validate these fields.
