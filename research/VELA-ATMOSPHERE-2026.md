# Vela atmosphere and visible-light rendering: research audit

Research date: 2026-10-17

## Decision relevant to this remake

The requested visible-light camera view is not an X-ray atmosphere image or a
gamma-ray particle-density picture. The visible optical spectrum of Vela is
predominantly nonthermal magnetospheric emission. A hydrogen surface atmosphere
is a plausible interpretation of its soft X-rays, but the accessible NSMAXG
tables do not extend to visible wavelengths or supply a local angle-dependent
optical intensity map. No extrapolation of those tables should be described as a
calculated Vela atmosphere.

Three distinct regions must remain separate:

- **Surface atmosphere:** a thin, hydrostatic, magnetized radiating layer
- **Magnetosphere:** charged particles and electromagnetic fields near the star
- **Pulsar wind nebula:** an extended shocked outflow on much larger scales

An animated glow around magnetic lines is not a surface atmosphere calculation.

## Primary observational constraints

### Optical spectrum: the first quantity to validate in a visible-light render

[Mignani, Zharikov & Caraveo 2007, A&A 473, 891–896](https://arxiv.org/abs/0707.2036)
([DOI](https://doi.org/10.1051/0004-6361:20077774)) measured the VLT/FORS2 spectrum.
The usable 4000–8000 Å continuum is described by
`F_nu proportional to nu^(-alpha)`, with `alpha = -0.04 +/- 0.04` after their
extinction correction. They used `E(B-V)=0.05`, corresponding to `A_V=0.16`, and
reported `V=23.6`. Their interpretation is magnetospheric optical emission;
extrapolating the X-ray power law overpredicts the optical continuum. A
particle-radiation model should therefore be tested at optical wavelengths,
rather than recoloring a gamma-ray calculation. This is a phase-averaged
integrated spectrum, not a resolved spatial emissivity reconstruction.

### Ultraviolet thermal limit is conditional, not an optical atmosphere table

[Romani, Kargaltsev & Pavlov 2005, ApJ 627, 383](https://arxiv.org/abs/astro-ph/0503331)
found complex nonthermal ultraviolet emission through all rotational phases.
Their pulse-minimum Rayleigh–Jeans upper limit implies
`T_eff,infinity < 4.6e5 (d_300 / R_14.2)^2 K`, conditional on their thermal-tail
interpretation and adopted extinction. This constrains a thermal contribution;
it does not provide local optical beaming, a surface temperature map, or a
magnetospheric particle distribution. A blackbody/Rayleigh–Jeans comparison may
be shown as a separate approximation or bound, never called a solved
magnetized atmosphere.

### X-ray atmosphere constraints have model-dependent masses and radii

[Pavlov et al. 2001, ApJ 552, L129](https://arxiv.org/abs/astro-ph/0103171)
fit magnetic hydrogen atmosphere plus nonthermal emission, using `B=3e12 G`,
`M=1.4 M_sun`, `R_infinity=13 km`. They obtained
`T_eff,infinity=0.68 +/- 0.03 MK` and `d=210 +/- 20 pc`. Those values belong
to that historical joint fit: they must not be combined silently with a later
distance and advertised as a new fit. `R_infinity=R/sqrt(1-2GM/Rc^2)` is the
apparent radius, not the physical radius.

[Sanwal et al. 2001, X-ray Spectrum and Pulsations of the Vela Pulsar](https://arxiv.org/abs/astro-ph/0112164)
report the alternative atmosphere fit `kT_infinity=59 +/- 3 eV`,
`R_infinity=15.5 +/- 1.5 km` at 250 pc, alongside a substantially hotter,
smaller-area blackbody alternative. These alternatives demonstrate that a
blackbody temperature is not interchangeable with an atmosphere effective
temperature.

[Ofengeim & Zyuzin 2018, Particles 1, 194–202](https://www.ioffe.ru/astro/Stars/Paper/ofengeim_zyuzin18.pdf)
([DOI](https://doi.org/10.3390/particles1010014)) reanalyzed Chandra using NSMAX
model 1230, a magnetic hydrogen atmosphere with `B=2e12 G`, plus a power law.
They fixed `d=290 pc` and explored mass/radius variation. Their summarized
redshifted temperature is `0.66(+0.04,-0.01) MK`; at the specific
`M=1.4 M_sun, R=13 km` point it is `0.700 +/- 0.005 MK`. They explicitly
identify uncertainty in atmosphere composition, magnetic field, higher
multipoles and polar-cap contributions. This remains an X-ray constraint and
does not license a visible-light atmospheric map.

### Distance: avoid the preprint/journal discrepancy

[Dodson et al. 2003, ApJ 596, 1137](https://doi.org/10.1086/378089)
has a published parallax `3.5 +/- 0.2 mas`, corresponding to
`287(+19,-17) pc`, as also recorded by the
[University of Western Australia author repository](https://research-repository.uwa.edu.au/en/publications/the-vela-pulsars-proper-motion-and-parallax-derived-from-vlbi-obs/).
The current [arXiv abstract](https://arxiv.org/abs/astro-ph/0302374) instead
quotes `3.4 +/- 0.2 mas`, `293(+19,-17) pc`. Cite the published version when
using 287 pc and record which version supplied a parameter.

## Public atmosphere model: what is actually downloadable

[Ho, Potekhin & Chabrier 2008, ApJS 178, 102](https://arxiv.org/abs/0802.2957)
construct coupled, polarized, partially ionized magnetized-hydrogen atmosphere
solutions in radiative and hydrostatic equilibrium. Their internal transfer
grid starts at 0.01 keV, already above visible photon energies; the atmosphere
scale height is about a centimetre. They emphasize anisotropic emission and
the failure of a simple blackbody to represent the atmosphere.

[Ho 2014, NSMAXG](https://arxiv.org/abs/1311.5583)
([DOI](https://doi.org/10.1017/S1743921314000656)) describes the expanded model
implementation. The [current official XSPEC manual](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/node240.html)
documents spectral models depending on surface effective temperature, mass,
radius, distance and a discrete magnetic/composition model. Publicly documented
energy ranges are 0.05–10 keV for single `B,T` models and 0.09–5 keV for dipole
distributions. Neither includes visible 400–700 nm, which is
0.001771–0.003100 keV. Table interpolation is legitimate within its axes;
optical extrapolation is not an atmosphere prediction.

The [legacy NSMAXG page](https://heasarc.gsfc.nasa.gov/docs/software/xspec/models/nsmaxg.html)
still points at older localmodel files. Its model-list/source links returned
404 during this audit. The [current HEASARC model-data instructions](https://heasarc.gsfc.nasa.gov/docs/software/xspec/modeldata.html)
provide individual files and `ftgetmodeldata`; downloading the whole 14 GB
archive is unnecessary.

Verified current file:

```
URL: https://heasarc.gsfc.nasa.gov/FTP/software/xspec/spectral/modelData/nsmaxg.fits
Bytes: 1301760
SHA256: 848ec0c597e0a64627787e7783ba797cf43f52abc01fd406a017e1e2b89dffff
Official index: https://heasarc.gsfc.nasa.gov/FTP/software/xspec/spectral/modelDataFiles.csv
Audit download: /tmp/nsmaxg.fits
```

The downloaded hash matches the official CSV entry. Direct FITS inspection
found `MODELINFO` columns `MODNUM, MODNAME, TEMPS, GRAVS, ENERGIES`; spectral
HDUs contain `FLUX` arrays. Every hydrogen table in this file has energies
from 0.047862999 to 10 keV, with 117 or 233 samples. There are no local
`mu, azimuth` intensity axes. Some integrated dipole models are supplied for
two viewing orientations; those are not a full angular surface-intensity
table. An isotropic intensity reconstructed as `F_E/pi` would be an added
assumption, contrary to the angular physics of the atmosphere. The FITS file
was inspected for feasibility and has not been used to generate an optical
render.

## Implementation requirements for a scientifically defensible optical bundle

These are recommendations inferred from the evidence above, not results of a
completed Vela fit:

1. Compute optical synchrotron/synchro-curvature radiation from the actual
   momentum histories, including pitch angle, local field and radiation losses
2. Select photons that reach the camera, evaluate retarded times and account
   for stellar occultation; all-direction radiation energy is not a camera
   image
3. Integrate the optical spectrum through a stated sensor or CIE response,
   then encode into sRGB; choose and disclose exposure and tone mapping
4. Compare integrated optical slope and flux with observations; a fitted
   normalization is a calibration, not a first-principles luminosity prediction
5. Do not invent glowing magnetic lines, purple cones, or particle-path trails
   in the optical camera output; those belong in optional diagnostic outputs
6. Keep a surface thermal approximation separate from a solved atmosphere
   claim, or leave the optical atmosphere contribution unmodeled and explicitly
   list that requirement as unresolved

A genuine atmosphere match additionally needs visible-band polarized opacity,
equation-of-state and transfer solutions producing
`I_nu(T,B,Theta_B,g,mu,phi)` over the surface, a supported surface temperature
and field map, and relativistic ray tracing. No appropriate public local
optical intensity table was identified in this audit. The absence of a found
table is not a claim that no such author data exist.

## Physical scale and the full-frame request

For the existing `P≈0.08937 s`, `R_LC=cP/(2pi)≈4264 km`.
A 13 km physical radius is only `0.00305 R_LC`. A frame spanning several
light-cylinder radii necessarily resolves the star by only a few pixels.
Filling the Reel can mean cropping/zooming the simulated emitting region, or
a continuous camera approach, but cannot simultaneously make the star large
and keep the complete magnetosphere at a common physical scale. A visible
star enlargement must be disclosed if selected; it cannot be presented as an
optical camera calculation.

## Follow-up availability audit: author archives and runnable source

The user subsequently specified that both magnetospheric plasma and the
stellar surface atmosphere must be calculated locally. This bounded follow-up
looked for an actual magnetized optical atmosphere input or runnable solver;
it did not relax that requirement.

### Ioffe electronic atmosphere archive

The [author archive](https://www.ioffe.ru/astro/Stars/) is reachable through its
legacy read-only CGI entry form. Its [model description](https://www.ioffe.ru/astro/Stars/data1.html)
explicitly states that the publicly presented models are **nonmagnetic** pure
H, He and Fe. The [file-format description](https://www.ioffe.ru/astro/Stars/data3.html)
does supply angle-dependent emergent intensity files (`I*T**G***.data.gz`),
but this does not overcome the field restriction. The inspected selection
inventory has composition, temperature and gravity axes, with no B axis.
These files therefore cannot be relabeled as a Vela-field atmosphere.

### Ioffe magnetic equation-of-state archive

[Potekhin's public magnetic hydrogen tables](https://www.ioffe.ru/astro/NSG/Hmagnet/hmagtab.html)
cover relevant magnetic field and temperature ranges and include partial
ionization. Their documented 14 columns supply thermodynamics, species
fractions and two **Rosseland mean** opacities, parallel and perpendicular to
the field. Those means contain no wavelength axis and cannot supply
`I_lambda(mu,phi)` in the visible band. The interpolation program interpolates
these same tables; it is not the coupled frequency-dependent atmosphere
solver. Published monochromatic-opacity calculations are not equivalent to
a public downloadable opacity grid with verified axes.

### X-PSI and McPHAC

[X-PSI's surface-radiation tutorial](https://xpsi-group.github.io/xpsi/Surface_radiation_field_tools.html)
provides public NSX angular-intensity tables and loading tools, but explicitly
describes the demonstrated H/He atmospheres as **nonmagnetic**. Its general
ability to accept user-defined atmosphere tables does not supply the missing
magnetic input. [McPHAC](https://github.com/McPHAC/McPHAC) is a useful public
nonmagnetic hydrogen atmosphere benchmark, not a direct replacement for
Vela's strongly magnetized surface. Neither was adopted as a Vela atmosphere.

### Public magnetic Julia prototype: an experimental route, not a validated input

[tobiasosborne/neutronstars](https://github.com/tobiasosborne/neutronstars)
contains a public two-mode magnetic transfer implementation. Inspection of
its [magnetic atmosphere module](https://github.com/tobiasosborne/neutronstars/blob/master/src/atmosphere/magnetic_atmosphere.jl)
and [visible rendering script](https://github.com/tobiasosborne/neutronstars/blob/master/scripts/render_rxj1856_visible_magnetic_atmosphere.jl)
identified reasons not to adopt it as the requested rigorous surface model:

- The visible example uses **fully ionized** H and a different neutron star
- The example deliberately uses `B=1e12 G`; its own comments say a roughly
  `3e12 G` case does not satisfy the local temperature-correction tolerance
- The solver documents a heuristic global temperature/flux correction and
  approximate mode-local scattering rather than full inter-mode redistribution
- Output intensity has frequency, emission-cosine and mode axes; an arbitrary
  inclined magnetic field additionally needs verified azimuthal dependence

Those findings describe the inspected version, not every possible future
version. The code was read, not executed or incorporated. Its existence is a
possible research starting point, not evidence of agreement with Vela or a
ready atmosphere component. Extending and independently validating such a
solver would itself be substantial new atmosphere research.

### Result and next concrete requirement

**No verified public atmosphere table or author-provided solver was found that
currently satisfies all required optical coverage, magnetic field, angular
intensity and partial-ionization requirements.** This is an availability result
of the inspected sources, not a statement that suitable author data cannot
exist.

A viable input would be an author-supplied or independently solved grid at the
chosen `T_eff, B, Theta_B, g`, extending at least across the **local**
330.3–578.1 nm band for the current `g=0.8258073` and observer 400–700 nm.
The calculation also needs surrounding frequencies sufficient to enforce
bolometric radiative equilibrium; solving only this narrow optical interval
cannot determine a consistent atmosphere temperature structure. Coupled
polarization modes, partial ionization, hydrostatic equilibrium, transfer
boundary conditions and frequency/depth/angle convergence must be checked.
An external-data request would require explicit user authorization before
contacting authors; no messages were sent.

The optical continuum observations may bound an additional thermal component,
but such a bound **does not meet** the user's requirement to calculate the
stellar atmosphere. A blackbody or Rayleigh–Jeans placeholder remains excluded.
