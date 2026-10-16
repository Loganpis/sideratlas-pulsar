# Radiation upgrade: implemented physics and remaining inputs

Audit date: 2026-10-13. Implementation: `src/pulsar/nonthermal.py`; verification: `tests/test_nonthermal.py`.

The project now includes classical synchrotron spectra and inverse-Compton spectra with Klein–Nishina recoil, in addition to its existing curvature-radiation calculation. These are **tested single-particle kernels**, not a solved Vela particle population or an automatically combined three-component Vela spectrum. The new Maxwell solver does not yet supply the particle energies, pitch distribution, densities, or seed-photon directions needed to use them as a predictive emission model.

## Synchrotron calculation

The input Lorentz factor, magnetic field and pitch angle must describe the same local magnetic frame, where the usual magnetic gyration formula applies. A laboratory magnetic field and pitch cannot be substituted indiscriminately when there is substantial electric drift. Frame conversion and angular emission toward an observer are caller responsibilities.

For a fixed pitch angle $\alpha$, the implementation uses

$$P_{\rm syn}=2\sigma_Tc\frac{B^2}{8\pi}(\gamma^2-1)\sin^2\alpha,
\qquad E_c=\frac{3\hbar eB}{2m_ec}\gamma^2\sin\alpha.$$

The familiar $4/3$ power coefficient emerges only after averaging an isotropic pitch distribution. The photon spectrum uses the normalized broad kernel $F(x)=x\int_x^\infty K_{5/3}(u)\,du$, not a monochromatic approximation. Outputs are angle-integrated photon rates per particle in supplied energy bands. Reusing the existing curvature spectrum's normalized mathematical kernel does not reuse its power or critical energy.

The API requires $\gamma\ge10$, and, for nonzero emission, $\gamma\sin\alpha\ge10$. It rejects $\chi_e=\gamma B\sin\alpha/B_Q>0.1$; this is a regime guard, not a claim of negligible quantum corrections near the threshold. It is not a cyclotron or quantum-synchrotron solver. Zero pitch and zero magnetic field produce zero emission.

## Inverse-Compton calculation

The implemented kernel is the ultrarelativistic isotropic Jones/Blumenthal–Gould expression, also given by [Lefa, Kelner & Aharonian (2012), equations 2–5](https://arxiv.org/html/1205.2929v1). For seed energy $\epsilon$, define $\eta=4\gamma\epsilon/(m_ec^2)$ and $q=E_\gamma/[\eta(\gamma m_ec^2-E_\gamma)]$. The code retains the full recoil term in

$$G(q,\eta)=2q\ln q+(1+2q)(1-q)+\frac{(\eta q)^2(1-q)}{2(1+\eta q)},$$

with $1/(4\gamma^2)\le q\le1$. The scattering rate is integrated over positive seed-energy quadrature weights. Each weight is a photon number density already integrated over its energy interval, not a density per unit energy.

The implementation provides spectra, band photon rates, emitted photon power, intercepted seed power, and **net electron cooling**. The cooling subtracts incident photon energy, so emitted power equals electron energy lost plus seed energy intercepted. Integration in $\log q$ uses the exact band boundaries; photons above the kinematic limit are excluded.

Inputs require $\gamma\ge10$ and seed energies at most 1% of the minimum supplied electron energy. This ultrarelativistic kernel omits the small down-scattering branch; it is not an arbitrary-energy exact Compton redistribution solver. Extreme Klein–Nishina losses are ensemble averages. Continuous cooling alone would not reproduce stochastic large-energy scattering events.

## Seed photons and populations remain explicit

`PhotonField` accepts caller-supplied isotropic photon densities. The Planck helper integrates $10^{-7}\le\epsilon/(kT)\le80$ with supplied temperature and dilution; it has no hidden Vela defaults. The stellar dilution helper uses

$$W(r)=\frac{1-\sqrt{1-(R_\star/r)^2}}{2}.$$

It preserves the flat-space photon energy density of a uniform spherical blackbody. Its name and description explicitly mark the result **isotropized**: real stellar photons propagate outward. It does not include anisotropic Compton geometry, surface temperature maps, atmosphere beaming, general-relativistic propagation, or attenuation.

No pitch-angle distribution, pair multiplicity, injection rate, seed temperature, bulk Doppler factor, cascade, or emission mixture is inferred by these functions. In particular, **do not add synchrotron losses to curvature losses calculated from a fully resolved gyromotion trajectory**: that would count the same acceleration twice. Such trajectories need a consistent total-radiation or synchro-curvature treatment. A separate guiding-centre curvature and unresolved gyro-emission approximation must state its scale separation and evolve the corresponding energy and pitch consistently.

## Verification evidence

All **12 tests passed** on 2026-10-13:

- Fixed-pitch power and critical energy agree with independent Larmor expressions; isotropic pitch integration recovers the $4/3$ coefficient
- Integrating photon energies over broad logarithmic bands recovers bolometric synchrotron power within $5\times10^{-6}$
- Planck number and energy densities agree with their analytic moments within $2\times10^{-9}$; stellar dilution recovers the surface and far-field limits
- In the Thomson limit, scattering rate and net cooling agree with $\sigma_Tcn$ and $(4/3)\sigma_TcU(\gamma^2-1)$ within $3\times10^{-6}$
- A separate integration in emitted photon energy checks Klein–Nishina cooling within $2\times10^{-9}$; density scaling, band additivity, cutoffs and energy accounting are checked independently
- Doubling thermal and scattering quadratures from 96 to 192 nodes changes the tested cooling by less than $3\times10^{-6}$ and band rates by less than 0.6%

These are kernel verification tolerances, not uncertainties on Vela. Reproduce with `PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -p 'test_nonthermal.py' -v` from the project root.

## Claim supported now

“Verified synchrotron and isotropic inverse-Compton building blocks for an evolving pulsar model” is supported. “Vela's multi-process spectrum reconstructed” is not. All previously inspected 3PC energy profiles are now **development data**, including the earlier reserved bands, because their results informed this redesign. Numerical verification and a later fit to those profiles cannot be relabelled independent observational validation.
