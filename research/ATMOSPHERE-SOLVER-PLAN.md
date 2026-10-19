# Local magnetized-hydrogen atmosphere: staged computation plan

Date: 2026-10-17

## Next step and claim boundary

Build and verify **one plane-parallel local atmosphere column**, initially with
the field normal to the surface. Do not render it as Vela's completed surface.
This makes the first nonlinear transfer problem small enough for local
benchmarking while preserving an upgrade path to inclined-field surface
patches. A validated normal-field patch alone cannot represent the whole star.
The existing redshift/input-audit utilities are useful infrastructure, not an
atmosphere solution.

An explicit trial point is `M=1.4 M_sun`, `R=13 km`,
`g_surface=1.3313e14 cm s^-2`, `B=3e12 G`, and
`T_eff=0.700 MK / 0.8258073 = 0.8477 MK`. These are **chosen model inputs**,
not jointly measured Vela properties. The temperature reflects a particular
X-ray-fit point in [Ofengeim & Zyuzin 2018](https://doi.org/10.3390/particles1010014),
whose mass/radius/composition uncertainties remain relevant.

## Physically minimal column

Use static, pure-H, LTE, plane-parallel hydrostatic and radiative equilibrium,
with no external irradiation for the initial problem. Solve the **coupled
ordinary/extraordinary polarization-mode transfer equations**, including
angle-dependent absorption and mode-changing scattering. Supply partial
ionization, a consistent EOS, magnetic free-free, bound-free/bound-bound and
electron/ion scattering opacities. At Vela temperatures a fully ionized
substitution needs a demonstrated error bound; it is not the default.
[Potekhin & Chabrier 2003](https://arxiv.org/abs/astro-ph/0212062) provide the
opacity framework, including the important low-frequency electron–proton
behavior. Public EOS/Rosseland tables alone do not provide all monochromatic
atomic opacity inputs.

Iterate the temperature/pressure structure with actual local energy-balance
residuals. Include radiative acceleration in hydrostatic balance, or bound
its omission quantitatively. A Planck source term required by LTE and a thermalized deep boundary
are physical transfer ingredients; replacing the emergent spectrum by a
blackbody is excluded. Use zero incident intensity at the outer boundary and
a sufficiently deep diffusion boundary. Implement resonance treatment and
vacuum-mode conversion, or quantify their omission before accepting a case.
[Suleimanov, Potekhin & Werner 2009](https://arxiv.org/abs/0905.3276)
([published DOI](https://doi.org/10.1051/0004-6361/200912121)) give a suitable
coupled-transfer and temperature-correction reference; their energy grid spans
0.001–20 keV and their published convergence criteria include flux and flux
derivative errors. The related conference paper
[arXiv:0909.2997](https://arxiv.org/abs/0909.2997) is a different publication.

For the present compactness, observed 400–700 nm requires **local
330.3–578.1 nm**. Solve a broad spectrum for bolometric equilibrium, explicitly
resolve this local optical interval, and extend energy boundaries until their
effect is negligible. A narrow optical-only equilibrium calculation is
insufficient. Inclined fields later require both emission cosine and azimuth.

## Benchmarks and acceptance gates

The numerical tolerances below are proposed project gates, not claims that
they have already been achieved or universal literature standards.

| Gate | Independent comparison | Required evidence before proceeding |
|---|---|---|
| 1 — Transfer baseline | Public [McPHAC code](https://github.com/McPHAC/McPHAC), [Haakonsen et al. 2012](https://doi.org/10.1088/0004-637X/749/1/52), and analytic absorption-only slab/diffusion tests | Match identical **nonmagnetic** inputs and angular conventions; target ≤1% spectral/angular differences where resolved signal exceeds a declared floor; establish grid convergence separately |
| 2 — Magnetic microphysics | Potekhin–Chabrier opacity equations/figures, [public magnetic EOS tables](https://www.ioffe.ru/astro/NSG/Hmagnet/hmagtab.html) | Verify units, detailed balance, polarization normalization, B→0 behavior and resonance locations; obtain or independently reproduce missing monochromatic atomic inputs before partial-ionization claims |
| 3 — Magnetic equilibrium | [Ho, Potekhin & Chabrier 2008](https://arxiv.org/abs/0802.2957), SPW2009 published cases | Compare the same B, T, gravity, column depth and boundary physics; digitized curves carry digitization uncertainty and cannot support finer accuracy than their resolution |
| 4 — Numeric X-ray spectrum | Official [NSMAXG FITS](https://heasarc.gsfc.nasa.gov/FTP/software/xspec/spectral/modelData/nsmaxg.fits) and [model documentation](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/node240.html) | Compare flux with the **matching discrete table model**, gravity and redshift; integrate local I over solid angle correctly; this tests only overlapping X-ray wavelengths, not optical correctness |
| 5 — Optical prediction | Independent magnetic transfer result or author intensity grid; then optical depth/angle/energy refinement | Require <1% optical-band intensity changes under independent grid/depth refinement, with a stated grazing-ray error criterion; no optical extrapolation from NSMAXG |

Use an **exact NSMAXG benchmark before the Vela trial**: e.g. hydrogen model
1230, `B=2e12 G`, field normal, `log g[cgs]=14.3802`, and an actual tabulated
temperature. The trial column's gravity and field do not exactly match that
table. A comparison with silently changed parameters is not validation.

For every converged case, record `max|F(m)/sigma T_eff^4−1|`, local radiative
equilibrium residual, hydrostatic residual, boundary-depth sensitivity and
iteration history. Initially require flux error <1% and a defined normalized
flux-gradient residual `max|dF/dln m|/(sigma T_eff^4) <1e-4`; tighten the flux target toward 0.1% and demonstrate
that tighter iteration tolerances do not alter the optical result. Uniformly
rescaling an unconverged spectrum does not pass. Preserve failed runs as
failures and never interpolate them into a surface map.

## After the column passes

Compute a converged library over required `T_eff, B, Theta_B, g`; retain
`I_lambda(mu,phi)` and exact source/version hashes. Ray trace a surface map with
gravitational bending/redshift, occultation, rotational shifts and correct
spectral Jacobians. Temperature map, atmosphere composition, field multipoles,
hydrogen column depth and return-current heating remain choices needing
sensitivity studies. Do not couple plasma deposition as heating until its
physical units and energy budget are validated.

Check the combined surface-plus-magnetosphere optical output against the
observed continuum and pulse data. [Mignani et al. 2007](https://arxiv.org/abs/0707.2036)
constrain the integrated optical continuum; they do not provide a resolved
surface map or independently validate a thermal atmosphere. A match in total
brightness can hide compensating errors in the two components.

## Feasibility and unresolved estimates

The immediate runnable milestone is the nonmagnetic transfer comparison plus
a magnetic-opacity audit, **not an optical Reel**. The present hard input gap is
verified monochromatic, partially ionized magnetic opacity/atomic data and an
independent optical intensity reference. Obtaining author data would need
separate authorization before contacting anyone; none has been requested.

A trial `128 depths × 512 energies × 16 mu × 2 modes` intensity array is about
16 MiB; eight azimuths make it about 128 MiB before work arrays/operators.
These are array-size estimates, not measured total memory or runtime. Actual
iteration counts, resonance refinement, stability at `B≈3e12 G`, full library
cost and research time are unknown until a benchmark is profiled. Additional
local computing cannot replace missing opacity physics or independent
validation. If those gates remain unmet, report the surface component as
unfinished and keep it out of a claimed realistic Vela render.
