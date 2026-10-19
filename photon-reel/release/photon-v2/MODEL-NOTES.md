# SiderAtlas NO.002 — calculated outgoing gamma rays

This is an accurate display of rays from the existing **approximate Vela-informed model**, not a newly validated Vela reconstruction. The original spectrum and pulse discrepancies remain. No new field solution, particle trajectories, fitting, pair physics or parameter tuning were introduced. The original v1 files are preserved.

## Why v1 looked like a torus

`prepare_reel_scene.py` deposited angle-integrated curvature-radiation energy at particle emission positions inside 2 light-cylinder radii. `render_reel.py` rotated and projected that volume. Its large toroidal structure represented **where the selected particles radiate**, not a freely propagating gamma-ray beam. A volume render cannot show the lighthouse effect merely by rotating its camera or adjusting its brightness.

The new renderer uses actual emission positions **and actual laboratory photon directions** from the same cached calculation. It displays photon-number-weighted packets moving away from their emission sites. The camera remains external, fixed and distinct from the distant observer direction. No cones, analytic dipole lines, jets, extra emission components or manual peak gains are drawn.

## Preserved physical calculation

- Final 97³ COFFEE force-free snapshot, magnetic inclination 75°, three field rotations, original source and boundary modifications
- 1,536 original streams, 128 azimuths × 6 rings × 2 hemispheres, 1,418,970 active integration cells, all reaching the original outer transport sphere
- Injection at 0.7 R_LC and emission followed to 2 R_LC; no new particle emission is invented outside that region
- Original guiding-center laboratory directions, trajectory curvature and signed acceleration/loss integration; initial Lorentz factor 1,000
- Prescribed acceleration 0.04 cm⁻¹ inside spherical r=R_LC and 0.25 cm⁻¹ outside; these are not self-consistently inferred electric fields
- Original three-node spectral quadrature for 0.1–0.3, 0.3–1, 1–3 and 3–100 GeV curvature photons
- Period 0.08937108859019084 s, R_LC=4,264.203109 km; fixed observer at ζ=65°, original 2° polar interval and 1,000 conservative arrival-phase bins

See the original [model notes](../MODEL-NOTES.md) for the complete physical assumptions, numerical diagnostics, failed refinements, source and boundary limitations, observational discrepancies and references. The same `vela_validation=false` and quantitative release restrictions remain in force. Publishing a visualization is not quantitative model validation.

The photon directions retain the original ultrarelativistic forward-beaming approximation along the calculated particle velocity. A finite angular curvature-radiation kernel is not newly solved. Display glyph width must not be interpreted as the physical 1/γ radiation opening angle.

## Photon geometry and time

Use dimensionless lengths R_LC and time Ωt, so c=1. For an original source-cell midpoint (r_j,n_j,t_j), form the corotating source template

\[
 r_0=R_z(-t_j)r_j,\qquad n_0=R_z(-t_j)n_j
\]

At an emission time t_e the physical source rotates:

\[
 r_e=R_z(t_e)r_0,\qquad n_e=R_z(t_e)n_0
\]

After emission the direction is fixed:

\[
 x(T)=r_e+(T-t_e)n_e
\]

Rotation applies at birth only. Photons in flight are never forced to corotate or follow curved particle or magnetic trajectories. Individual rays are straight; the extended pattern can curve because photons at different distances were emitted at different rotation phases. Its detailed sheets and spirals arise from the saved source positions, directions, weights and propagation times. They are not an imposed beam surface. Free streaming outside the modeled emission region continues the original flat-space propagation assumption, without new emitting particles, absorption, occultation, gravitational ray tracing or pair feedback.

The original distant observer convention is

\[
 \Phi=\frac{t_j-\operatorname{atan2}(n_{jy},n_{jx})-n_j\cdot r_j}{2\pi}\pmod1
\]

For a source oriented so that its emitted photon azimuth is zero, this becomes (t_e−n_e·r_e)/(2π). The observer is a **far-field angular pixel**, not a telescope placed near the pulsar. The crosshair is a scaled bearing at the fixed direction (sin65°,0,cos65°); its drawn location is schematic. The dotted line indicates that bearing and is not a light ray, field line or source of radiation.

For timing display, use a distant-wavefront coordinate D=5.2: a ray crosses n_e·x=D at T_cross=t_e+D−n_e·r_e. Subtracting D recovers precisely the original asymptotic arrival phase. The number D is a display time origin/range coordinate, not a measured Vela distance or a finite point-detector approximation. Each direction inside the original 2° polar pixel uses its own asymptotic wavefront, as in the original calculation. A finite nearby point detector would require different geometry and is not claimed here.

At lab phase p the crosshair flash and the trace use the original full profile at arrival phase p−D/(2π), integrated over one frame. There is no manually adjusted phase alignment. This preserves the relationship between emission direction, flight delay and received pulse. The external diagram's overall brightness is not what this one observer receives.

## Visualization approximations

**Finite packet quadrature.** 12,288 equal-photon-number source representatives are selected deterministically, ordered in 64 polar-direction × 240 arrival-phase strata. Each representative retains an actual source-cell midpoint, its unmodified unit photon direction, original source-cell IDs and four calculated spectral fractions. No angular kernel or finite beam cone is added. Systematic sampling approximates the continuous weighted source. The coarse 16×48 sky histogram has relative L1 error **3.650%** (sum absolute bin differences / total photons). This is a display-sampling diagnostic, not an error bound on every narrow feature. Original source discretization and finite injection assumptions remain.

The sampled band totals differ from the full source by −0.0036%, +0.1271%, −0.3038% and −1.0526%, respectively. No band is renormalized to conceal this. Full-resolution receiver timing and flux use **all original source cells with conservative finite-cell integration**, rather than counts of displayed dots. Individual drawn packets therefore should not be counted as a simulated detector event list. Their ray geometry is exact for the selected source midpoints; the full pulse also retains the original linear within-cell sky reconstruction.

**Steady emission sampling.** Each representative emits 48 times per turn, with staggered phases to avoid synchronized artificial shells. Packets from previous turns are retained for flight ages below 13 Ω⁻¹, approximately 13 R_LC of path. This finite display window limits clutter and cost; it is not photon destruction, a physical nebula boundary or attenuation. The movie repeats a periodic steady source, not a transient turn-on. Source radii stay at 0.7–2 R_LC while photons extend through most of the full vertical frame.

**Rasterization and annotations.** Orthographic camera: polar angle 55°, azimuth −70°, roll −30°. The camera does not orbit. A 2.8-pixel Gaussian glyph footprint at native 1080 width and two samples over the one-frame photon shutter suppress sampling aliasing. These are display glyphs, not physical photon widths, plasma diffusion or a beaming-angle calculation. The shutter spans 1/240 rotation. The short sparse arrows use actual packet positions and directions; their 0.24 R_LC drawn length is a velocity annotation, not additional emitted flux. The tiny 2.1-pixel central marker indicates position only. It is not a resolved neutron-star surface or a new correction to the original enlarged numerical conductor (radius 0.5 R_LC).

**False color and brightness.** Gamma rays are invisible. Explicit band colors are blue (0.1–0.3 GeV), cyan (0.3–1), cream (1–3), orange (3–100). Colors are mixed in linear light in proportion to the source cell's calculated band photon counts. Brightness represents weighted photon number, rather than bolometric energy or a telescope exposure. There is no separate boost for a rare energy band, pole, selected observer ray or bright peak. White/gray arrows, crosshair and trace are diagram annotations, not additional radiation. One exposure and one hue-preserving tone curve apply to every frame; these are display choices and do not imply absolute physical intensity. The supplied frames are sRGB, converted to Rec.709 limited-range YCbCr for video.

**Time and format.** One physical 89.371 ms turn is stretched to 8 seconds (89.51×). Four identical cycles produce 32 seconds, 1080×1920, 30 fps. Smooth endpoint periodicity follows the periodic source and steady photon population. The header retains only @sideratlas, THE PULSAR and NO.002. Scientific explanation is in the accompanying caption and notes.

## Verification before export

The saved scene audit precedes video export:

| Check | Measured result |
|---|---:|
| Recomputed original four-band conservative pulse vs saved original | Exactly identical arrays |
| All source midpoint arrival phases vs saved sky phases | Maximum 8.88×10⁻¹⁶ cycles |
| Unit direction error, all reused active emitting cells | Maximum 4.44×10⁻¹⁶ |
| Displayed birth direction norm error | Maximum 4.44×10⁻¹⁶ |
| Straight-line, c=1 transport error | Maximum 4.44×10⁻¹⁶ R_LC |
| Far observer wavefront phase error | Maximum 2.22×10⁻¹⁶ cycles |
| Wavefront plane crossing error | Maximum 3.55×10⁻¹⁵ R_LC |
| Periodic packet-position endpoint error | Zero |

These verify faithful visualization of the retained model. They do **not** validate that model against Vela. The original development spectrum comparison remains χ²≈32,743 / 6 descriptive degrees of freedom and the pulse is discrepant. No claim of observational agreement is made.

Machine-readable source, calculation, direction/timing and encoded-video checks accompany the delivery. The focused regression suite tests source/event round trips, straight photon motion, far-observer phase, light-travel delay, frame-integrated pulses, conservative quadrature and original radiation transport.

## Reproduce

From the original project root, using its existing geometry cache:

```sh
python prepare_photon_scene.py
PYTHONPATH=src python -m pytest -q tests/test_photon_scene.py tests/test_emission_quadrature.py tests/test_outer_transport.py tests/test_radiation.py
python render_photon_reel.py
python verify_reel.py outputs/sideratlas-pulsar-photons/media
```

The supplied `scene/photons.npz` is sufficient to rerender without any field or particle calculation. Preparation refuses to overwrite an existing scene. Python 3.12+, NumPy, SciPy, Pillow and FFmpeg with zscale/libx264 are required. Model code retains GPL-2.0; font and media rights follow the original project notices.

## Studies and provenance

Spitkovsky (2006), [force-free magnetospheres](https://arxiv.org/abs/astro-ph/0603147); Kalapotharakos, Harding & Kazanas (2014), [particle directions and gamma-ray pulse formation](https://arxiv.org/abs/1310.3545); Barnard et al. (2022), [prescribed acceleration and Vela curvature emission](https://arxiv.org/abs/2111.03405); Smith et al. (2023), [Fermi-LAT 3PC](https://arxiv.org/abs/2307.11132), used for original timing and development comparisons. The force-free solver derives from [COFFEE](https://github.com/fizban007/CoffeeGPU), Yajie Yuan, Alex Chen and contributors, with the original project's full patch chain retained. [Brambilla et al. (2018)](https://arxiv.org/abs/1710.03536) and [NASA SVS 13058](https://svs.gsfc.nasa.gov/13058/) remain visual research references, not sources of photon data or media in this render. No new observational data or NASA footage were used.
