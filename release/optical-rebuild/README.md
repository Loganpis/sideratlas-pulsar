# Vela optical reconstruction — research rebuild

The requested product is a physical visible-light camera image of a pulsar: a coupled particle/field simulation, physical optical emission and a supported Vela surface atmosphere, filling a 1080 × 1920 Reel. **That product is not yet available.** The previous emission-map and star/field-guide movies do not satisfy this specification and are not presented here as substitutes.

## Work completed

- Pinned and compiled the real relativistic PIC framework **Entity**, with **Kokkos** on the local CPU
- Passed four unmodified upstream pusher, current-deposition and Maxwell checks
- Completed a coupled three-dimensional periodic pair-plasma smoke run: **16³ cells, 16,384 macroparticles, 55 timesteps**. It has no pulsar boundary, physical Vela scaling, atmosphere or radiation; it is a toolchain check only
- Built and ran an original **48³ oblique pulsar PIC prototype** for approximately one commanded stellar turn. It evolves electron/positron motion and the electromagnetic fields together. An otherwise identical control with current deposition disabled demonstrates real particle-to-field feedback. [Source, configuration and reproduction](../../research/pic-local/README.md)
- Audited native Yee-grid constraints, stable particle lifetimes, charge export, numerical-source drift and physical scale differences. These reveal unresolved source/boundary and scale limitations; they do not validate Vela or a visible-light prediction
- Researched the NASA-linked Brambilla PIC model, Vela optical spectrum and magnetized surface-atmosphere constraints
- Added independently tested radiation, atmosphere-input and spectral-color building blocks. These are not an assembled optical-camera solver or a Vela production run
- Recorded the remaining production prerequisites in [status.json](status.json), with a full-frame camera delivery contract and the requested caption structure

## Why the published illustrations cannot be converted directly

NASA's 2018 image colors represent electrons, positrons and particle energies. Its reduced-scale kinetic model does not supply a literal Vela optical image. In particular, neither a particle's charge color nor the previous GeV curvature-emission map provides its optical spectral radiance toward a camera. Physical pair energies and pitch-angle evolution are needed. The NASA site distributes rendered frames and movies, not the underlying Brambilla model dataset.

The public NSMAXG atmosphere grids cover UV/X-ray energies and store spectral flux. They do not provide local, angle-dependent intensities at 400–700 nm. An X-ray fit temperature or a blackbody with that temperature cannot be relabeled as a calculated magnetized optical atmosphere. Vela's measured optical spectrum is dominated by nonthermal magnetospheric emission; the atmosphere remains a separate component.

Both the magnetosphere and the thin stellar atmosphere remain in scope, as requested. The selected computation route is local. The [surface-atmosphere computation plan](../../research/ATMOSPHERE-SOLVER-PLAN.md) identifies the missing monochromatic magnetic opacity inputs and separate numerical/optical verification needed for a local calculation. No atmosphere solver or surface intensity map has been completed.

At one physical scale, a roughly 13 km stellar radius is only about 0.003 light-cylinder radii. A full magnetosphere view therefore gives a small star. The new camera contract fills the Reel by selecting a physical field of view and crop, without enlarging the star independently or drawing explanatory beams.

## Source and verification records

- [PIC feasibility and source pins](../../research/PIC-REBUILD-FEASIBILITY.md)
- [Local native PIC verification](../../outputs/entity-pic-benchmark/verification.json)
- [PIC benchmark input and reproduction](../../research/pic-benchmark/)
- [Reduced pulsar prototype and measured limitations](../../research/pic-local/README.md)
- [Native prototype evidence manifest](../../outputs/pic-local/native-results.json)
- [Explicit scale and uniform-grid budget](../../outputs/pic-local/physical-scale-budget.json)
- [Observer radiation requirements and primary references](../../research/OBSERVER-RADIATION-REBUILD.md)
- [Vela atmosphere and optical constraints](../../research/VELA-ATMOSPHERE-2026.md)
- [Local surface-atmosphere computation plan](../../research/ATMOSPHERE-SOLVER-PLAN.md)
- [Requested-format caption draft](instagram-description.draft.text)

The draft is explicitly labeled unfinished. There is no new optical Reel or completed Vela release attached to this research rebuild. Unit tests and a successful smoke run do not establish kinetic convergence, a correct particle population or observational agreement.

To regenerate the prerequisite audit:

```sh
python audit_optical_rebuild.py
```

Independent radiation and atmosphere tests run with:

```sh
PYTHONPATH=src python -m pytest -q tests/test_optical_radiation.py tests/test_atmosphere_inputs.py tests/test_optical_color.py
```

The native solver and color datasets retain their own licenses and attributions. No NASA footage, raw simulation dataset or code has been incorporated into a new Vela movie.
