# Local atmosphere calculations — research only

The user requested plasma **and** the stellar surface atmosphere. This folder
contains actual local atmosphere transfer calculations toward that requirement.
It does not presently provide a validated Vela optical surface. None of its
spectra is renormalized to conceal an energy-balance failure, and no optical
blackbody or Rayleigh–Jeans spectrum is substituted for transfer.

## Final input selected for the disclosed approximate-model Reel

Use `outputs/optical-finish/atmosphere/polished-conservative-64/surface-intensity-si.npz`.
SHA256:
`0915abb9b0e9f4cedfd328e1a2e5e905a42b7c56ac004bd02245b688552d6f37`.
It contains ascending `wavelength_m`, emission cosine `mu` from zero to one,
and summed-mode `i_lambda_si[wavelength,mu]` in `W m^-2 sr^-1 m^-1`.
The explicit local optical nodes span 300–650 nm within a broad equilibrium
frequency grid. The adjacent JSON records physical assumptions and flags.

This **actual 64-node approximate calculation** passes its declared numerical
equation gates: maximum depth flux error `0.0009460` (0.0946%), local relative
energy balance `4.8591e-5`, and conservative energy-ledger error `2.4716e-7`
of target flux. The bounded polish took 57.5 seconds and stopped when the
explicit gates were reached; this is not an optimizer-`gtol` assertion.

Its absolute optical output is still numerically sensitive to depth:
32→64 nodes changes hemispheric 330–580 nm flux by **+4.68%**, median
non-grazing point intensity by 4.70%, and the worst point by **5.88%**.
The one-percent depth-refinement gate **fails**. The thermal nodes retain an
alternating-temperature artifact; the monotone sensitivity calculation was
stopped unfinished after bounded resources. Angular/frequency convergence,
a relaxed bottom-extension comparison and independent magnetic optical
validation have not passed. These must remain disclosed, with
`refinement_pass=false` and `scientific_validation_pass=false`.

The selected input is suitable only for a clearly described **approximate
model visualization**, as authorized for this Reel. It is not a validated
Vela atmosphere or a prediction of the observed surface. The complete
numerical and physical claim boundary is in `CONSERVATIVE-SOLVER-NOTES.md`.
`FINAL-SURFACE-HANDOFF.json` lists the exact minimal source/output copy set,
hashes and omitted external data.

## Actual magnetic calculation

`magnetic_column.py` implements a deliberately approximate, normal-field,
fully ionized hydrogen column at `B=3e12 G`, `g=1.3313e14 cm/s²` and chosen
`T_eff=847655.387 K`. Hydrostatic pressure is `P=g m`; the ideal electron/ion
EOS sets density. The thermal emissivity in each mode is `kappa_abs B_nu/2`.
The two polarization modes exchange radiation through coherent scattering.
The nonlinear temperature iteration has an unscaled target bolometric flux
`sigma T_eff^4`. It covers a broad spectrum; the later linear-source solver
adds actual local 300–650 nm nodes.

The free-free response uses the consistent electron–proton denominator in
[van Adelsberg & Lai 2006, §3.1.2](https://arxiv.org/abs/astro-ph/0607168).
This preserves the low-energy transverse suppression proportional to `E^4`
below **both** cyclotron energies; the older electron-only `E^2` prescription
is excluded. The per-polarization Gaunt factors are prescribed as unity in
the cited magnetic convention. This is a microphysical approximation, not a
claim that true magnetic Gaunt factors are unity. An attempted comparison at
zero field with standard nonmagnetic Kramers `g_ff=1` gave a 0.412 ratio;
those two unity conventions do not establish a matching limit. The magnetic
experiment now rejects `B=0`; the mismatch is retained as an audit finding.

Other approximations are the transverse **linear** X/O basis, no dispersion
or refraction, no partial ionization or bound opacities, no vacuum conversion,
no Compton redistribution, conduction, or external irradiation. These error
terms are not bounded at Vela's chosen parameters. A normal-field column also
cannot be a dipolar whole-star atmosphere library.

The first constant-source-cell integral calculation reached small normalized
local heating residuals but **failed the depth flux gate**. This demonstrated
why a small temperature correction or local residual is insufficient. Its
unaccepted files/log remain under `outputs/optical-finish/atmosphere/`.

`linear_magnetic_column.py` uses analytic linear-source short characteristics,
a thermalized diffusion boundary, and a hybrid temperature iteration: local
LTE energy balance in optically thin layers, constant flux in deeper layers.
Scattering is solved through its three cyclic moments, algebraically reducing
the angular coupled system without discarding mode changing terms. Both flux
and local residuals are still audited independently. `verification.json`
always distinguishes numerical equilibrium from scientific validation and
grid refinement. No result is accepted by virtue of optimizer exit alone.

Run from the project root, with SciPy and NumPy installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/optical-finish/atmosphere/polish_conservative_column.py outputs/optical-finish/atmosphere/checkpoint-refinement-64/input-checkpoint.npz --depths 64 --max-nfev 20 --output /tmp/reproduced-atmosphere-64
python research/optical-finish/atmosphere/export_conservative_surface.py /tmp/reproduced-atmosphere-64 --allow-unvalidated-physical-model
python -m pytest research/optical-finish/atmosphere/test_magnetic_column.py research/optical-finish/atmosphere/test_conservative_column.py -q
```

The seven numerical tests cover polarization completeness, corrected
electron/proton low-energy behavior, rejection of the unvalidated zero-field
limit, scattering detailed balance, an independent exact absorptive slab,
cell-integrated energy conservation, and an isothermal reservoir limit.
These do **not** validate Vela microphysics or a magnetic optical spectrum.

`export_conservative_surface.py` postprocesses the solved radiation field into explicit
angle-dependent `I_lambda` in `W m^-2 sr^-1 m^-1`, summing both modes. It
directly evaluates additional outgoing angles against the solved incoming
quadrature, including grazing rays and `mu=1`. It refuses a numerically
unaccepted input and requires an explicit acknowledgement of the unvalidated
physical model before export. Ten owned numerical/export tests pass.

## Actual nonmagnetic native baseline

As a prior bounded research check, the unmodified public
[McPHAC](https://github.com/McPHAC/McPHAC) repository was
downloaded at commit `ad6df20782bccbacd48dbec443f3f8e3afbccb88`, compiled with
Apple Clang and GNU Fortran 14, and run at `T_eff=1e6 K`, `g=2.43e14 cm/s²`.
The two runs used `(80 depths,64 frequencies,8 mu)` and
`(200 depths,128 frequencies,10 mu)`. They exited successfully in six
iterations; their actual full-transfer surface fluxes were not within our
1% project gate. See `mcphac-200/verification.json` for the independent flux
audit. Their spectral range starts in the ultraviolet, so no visible
intensity is extracted by extrapolation. A successful native run is useful
baseline evidence, **not** a passed independent magnetic atmosphere benchmark.

Build flags were `clang -O2 -fcommon -mmacosx-version-min=15.0`,
`gfortran -O2 -mmacosx-version-min=15.0 -std=legacy -fallow-argument-mismatch`,
plus the upstream BLAS/LAPACK and local GNU Fortran runtime. `-fcommon` keeps
the upstream tentative C globals compatible with the current compiler; no
source edits were made. An OPAL COMMON-block size warning is preserved in
the build log and has not been shown harmless by an independent EOS audit.
No explicit top-level source license was found. The vendor checkout is held
locally and excluded from publication; do not silently redistribute it.

## Actual neutral-fraction audit

`eos_audit.py` reads the author-distributed magnetic-H EOS tables at
`logB=12.4,12.5` and interpolates the actual atomic/molecular fractions and
pressure along a column. The source files are:

- [hmag12_4.dat.gz](http://www.ioffe.ru/astro/NSG/Hmagnet/hmag12_4.dat.gz),
  SHA256 `b9720e3e1929da767730174c2ffd56d5800a6fee84124911ce956e51b2363a71`
- [hmag12_5.dat.gz](http://www.ioffe.ru/astro/NSG/Hmagnet/hmag12_5.dat.gz),
  SHA256 `24a1cf7db1b3b03a1b3cd4b9ffdc7bc57e48cd1b1731b077f18d6c5a000de094`

The [author documentation](https://www.ioffe.ru/astro/NSG/Hmagnet/hmagtab.html)
defines the table axes/columns. The initial provisional column had atomic
fractions reaching about 11%. This is neither a completed atmosphere fit nor
a quantitative bound on the omitted atomic opacity. Rosseland means still
cannot supply missing monochromatic bound-bound/bound-free optical opacity.
The EOS is audited independently, not silently substituted into the fully
ionized transfer model.

The table downloads and vendor source are excluded from the published tree
pending license review. Their exact hashes and URLs remain in the reports.
No authors were contacted. Required future acceptance includes independent
magnetic spectrum comparison, radiative/hydrostatic equilibrium, optical
frequency/depth/angle refinement, boundary sensitivity, and a demonstrated
error budget for the missing microphysics.
