# Physical foundation of the gamma-ray animation

This folder contains the calculation behind the current 0.1–100 GeV photon
scene. The selected **97³ force-free snapshot, 1,536 particle streams,
prescribed acceleration, curvature emission and original observer profile**
are the same inputs used for the delivered animation. The animation's
ready-to-render photon scene is in [`../scene/`](../scene/).

- [Physical assumptions, numerical limits and source roles](MODEL-NOTES.md)
- [Selected field snapshot and run configuration](inputs/)
- [Selected emission protocol, results and pulse](diagnostics/n97-final-illustration-a128-r6/)
- [Field and injection audits](diagnostics/final-field-audit/)
- [Force-free source repair and checks](research/COFFEE-CURL-REPAIR.md)
- [Outer particle transport](research/OUTER-TRANSPORT.md)
- [Catalog data audit](research/VELA-DATA-AUDIT.md)

## Inspect the exact native source

The pinned upstream COFFEE source, dependency notices and applied patch
chain are included. Preparing its source needs no simulation or compilation:

```sh
python physics/build_native_release.py curl-radiative --prepare-only
```

The default creates `physics/native-work/curl-radiative/source`. To compile,
omit `--prepare-only`; a C++17 compiler, CMake, OpenMPI and parallel HDF5
are required. See the source provenance and numerical notes before changing
the field or its boundaries. A visually smooth animation does not establish
resolution convergence or agreement with Vela.

## Particle and emission calculation

The original runner and data are retained for inspection. Its Python modules
are shared with the renderer under `../src/pulsar/`. Optional dependencies:

```sh
python -m pip install -e '.[physics]'
python physics/run_emission_quadrature.py --help
```

Generating trajectories or rerunning a field is unnecessary for rendering.
`prepare_photon_scene.py` can also rebuild the photon display quadrature from
the original six cached trace chunks using `--geometry-cache PATH` and a
fresh `--output PATH`. Those large trace chunks are not required by the
supplied scene and are not silently regenerated.

Recorded protocols preserve their original paths and hashes as provenance.
The new repository layout changes file locations, not their physical
parameters or saved results. The enlarged numerical conductor, finite
injection surface, prescribed particle supply and acceleration, frozen field
rotation and discrepant Vela pulse/spectrum remain explicit limitations.
