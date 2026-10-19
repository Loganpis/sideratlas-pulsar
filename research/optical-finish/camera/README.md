# Calculated optical camera

This camera assembles an **approximate optical model of a plausible young rotation-powered pulsar**. Vela supplies reference scales; its luminosity, spectrum and pulse are not fitting targets. This is not a reproduction of NASA's global PIC simulation. The photon source, the local reduced PIC experiment and the surface atmosphere are separate calculations. The optical source does not silently rescale reduced PIC particles

Only `@sideratlas`, `THE PULSAR` and `NO.002` appear over a black, 1080 × 1920 image. A fixed observer, one physical image scale, one global exposure and one color transform apply throughout the loop. One 89.371 ms physical rotation is slowed to eight seconds, then repeated four times. Phase zero, image roll and crop are presentation choices, not measured pulse alignment

## Photon source and units

`../particle-emission/load_camera_events.py` reconstructs the compact archive without changing its geometry or radiation. `xyz_lm_sr` is exact-harmonic CIE1931 XYZ **power per physical lepton per steradian** with the factor 683 applied once. Physical source energy is `XYZ * number_weights * source_interval_s`. Counts are linear and incoherent

The prescribed steady co-rotating population has physical SI dipole fields and guiding-centre cooling/focusing. It is not a stored set of native PIC trajectories. The finite source intervals are `P/Nphase`, not the cooling residence intervals already represented by the number weights. Local bulk emission boosts have already been applied; no particle-velocity beaming factor is added in this camera

## Finite spatial and arrival integration

`finite_volume_camera.py` uses the physical injection-coordinate support corners from the photon calculation. `volume_raster.py` fits a uniform affine parallelepiped, intersects observer rays with it, and integrates source/arrival-window overlap along the ray. Its source emissivity is held at the calculated angular/radiative value of the quadrature node. The affine fit is not an exact curved dipolar cell or a solved emissivity distribution. Angular and spatial variation within a support remain uncalibrated; refinement statistics must accompany the image

Rigid motion is integrated by freezing each support at successive source-time midpoint quadrature nodes. Each piece's spatial light-travel-time spread is retained. For supports intersecting the image, the source/arrival-window overlap along each static ray is analytic, with 2×2 subrays per plasma-grid pixel. This is exact for that frozen-motion piece and the chosen pixel quadrature, rather than exact continuous motion. A rotation-invariant bound expands candidate bins across the period boundary

Wholly offscreen supports skip image rasterization. Their cropped-energy timing ledger instead uses 4×4×4 Gauss nodes; this approximates the arrival distribution of light outside the frame and has no effect on image pixels. Covering arrival bins conserve all retained source energy over the period. The output records motion step, affine corner error, raster volume error, crop loss, omitted energy and the full-cycle energy ledger

The spatial support is a finite-volume interpolation assumption, not a broadened radiation cone. It replaces numerical point streaks rather than inventing magnetic field-line light. No random luminous particle trails, educational cones or species colors are added. Up to 1e−7 of each supplied XYZ source-energy component may be omitted computationally, with its actual bound recorded

## Efficient image sampling

The selected Reel uses one physical field of view, 108 × 192 km in the virtual image plane, with a −10° image roll. Its output and stellar surface use 100 m per native pixel. Smooth finite-cell plasma radiance is integrated at 800 m per pixel (135 × 240), then bilinearly resampled onto the same 1080 × 1920 camera. The low-grid centre is shifted so the physical cell centres align exactly with native stellar rays; interpolation conserves integrated XYZ radiance. The original radiation and particle powers are unchanged. Sub-grid plasma morphology is not resolved

Cooling cells are grouped along their ordered cooling path at the same source phase when their nodes and eight support corners remain within a maximum distance of 200 m. Each group sums `XYZ * number_weights` exactly, preserving its already-calculated directional spectral color and source energy without angular or spectral averaging. The grouped geometry is approximate; the report records its actual displacement and timing bounds. The two-pixel motion step is in the low plasma grid (16 native pixels); the piecewise source/arrival overlap remains integrated along each static subray. A native 0.6-pixel PSF is applied after plasma resampling and surface addition

The fixed display divisor is half the maximum positive linear-RGB surface radiance: all source RGB values are divided by this same number before the common tone curve. It enhances the dim stellar disk between magnetospheric pulses without fitting any observed flux or altering the layers' linear light ratios. Bright plasma is compressed by the common tone curve

## Stellar surface

`surface_camera.py` evaluates outward Schwarzschild light-bending and relative travel-time integrals directly. A spherical 1.4-solar-mass, 13 km star is mapped through local SR aberration/Doppler factors. Specific intensity obeys `I_lambda,obs = (g_grav D)^5 I_lambda,em` at emitted wavelength `(g_grav D)*lambda_obs`. Frame dragging, oblateness and magnetospheric propagation of surface photons are absent

`atmosphere_camera.py` accepts an explicit angular optical SI specific-intensity table and checks its numerical-equilibrium flag and full mu=0..1 support. Physical validation is separate; using an unvalidated table requires an explicit approximate-model switch. It never fills optical data using a blackbody, Rayleigh–Jeans or X-ray extrapolation. One homogeneous normal-field column is repeated uniformly over the surface; this is not a dipolar atmosphere library or measured temperature map. The star and plasma use the same metres per image pixel. Stellar subpixels integrate its limb without enlarging the star

Surface GR and flat outer-plasma propagation have different approximations. A uniform surface has no thermal pulse to align, and this bundle makes no prediction of a radio/optical/thermal phase relation

## Color and export

The visible band is explicitly 400–700 nm, not the entire human response. Actual source harmonic wavelengths supply CIE XYZ. The surface spectrum is sampled after its redshift and integrated against the same official CIE table. There is no species/energy palette, source-dependent white balance or frame-dependent normalization

Positive linear sRGB channels share a fixed global exposure and the tone curve `rgb/(1+max(rgb))`, preserving their ratios before encoding. Negative out-of-sRGB-gamut channels are clipped; therefore narrow spectral colors can change under gamut clipping. A fixed 0.6-pixel Gaussian detector/reconstruction PSF is shared by both layers. Display exposure, tone mapping and omitted extinction mean this is not a naked-eye or calibrated telescope image

`reel_export.py` transforms actual sRGB sample values into limited-range Rec.709 YCbCr before encoding SDR H.264 High/yuv420p at 30 fps. The loop is encoded once and repeated without another generation of compression. `verify_export.py` measures the resulting container/stream tags, frame count, duration, lack of audio and decoded equality of repeated cycles. Instagram performs its own further transcoding

## Reproduction

Run from the repository root with Python 3.12, NumPy/SciPy/Pillow and FFmpeg with zscale. Obtain the exact photon archive and surface table listed in the release manifest, or run the provided source calculations. Then:

```sh
python download_optical_assets.py
python research/optical-finish/camera/render_parallel.py \
  --scene scene/optical/physical-pair-events.npz \
  --surface scene/optical/surface-intensity-si.npz \
  --output output/optical --width 1080 --bins 240 --fps 30 --repeats 4 \
  --workers 3 --meters-per-pixel 100 --roll -10 \
  --plasma-scale 8 --motion-arc-pixels 2 --exposure-reference surface
python research/optical-finish/camera/verify_export.py \
  output/optical/sideratlas-no002-optical-32s.mp4
```

The manifest pins the exact inputs, renderer parameters and measured export. Numerical tests establish the declared algorithms' properties; they do not establish Vela agreement
