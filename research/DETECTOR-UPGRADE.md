# Detector upgrade: calibrated response sensitivity

Audit date: 2026-10-13. Implementation: `src/pulsar/response.py`; verification: `tests/test_response.py`.

The earlier flat within-band detector weighting can now be replaced by the official **P8R3_SOURCE_V3 FRONT and BACK effective area and energy dispersion**. This calculates sensitivity to a real LAT response. It is **not yet a prediction of the published Vela 3PC weighted counts**, an exposure reconstruction, or an absolute Fermi spectral fit.

## Calibration provenance

The local files are pinned to the official [fermi-lat/irfs repository](https://github.com/fermi-lat/irfs/tree/b7cafb0a7a99ec36925cf8984372b4caafb5039c), commit `b7cafb0a7a99ec36925cf8984372b4caafb5039c`, retrieved 2026-10-13. Exact download URLs, sizes and SHA-256 hashes are recorded in `data/response/provenance.json`; the repository licence is preserved alongside the FITS files.

| File | SHA-256 |
| --- | --- |
| `aeff_P8R3_SOURCE_V3_FB.fits` | `2fee015b371b33feda646cf858f71c4c2d6071d8517c704233775cfab75d2adf` |
| `edisp_P8R3_SOURCE_V3_FB.fits` | `a2a4166c17af6a70e86775e86bc32bc1baba82b9c50eee266161d466822b774e` |

The [Fermi response overview](https://fermi.gsfc.nasa.gov/ssc/data/analysis/lat_irfs/irf_overview.html) documents the FRONT/BACK HDUs and IRF naming. The implementation reads effective area from those tables, converts m² to cm² and MeV to GeV, checks the effective-area version and units, and interpolates over logarithmic energy and incidence cosine. It refuses queries outside the calibrated true-energy range or outside $0.2\le\cos\theta\le1$; within the table edges, end-centre values are clamped rather than extrapolated.

## What is folded

For each measured-energy band $b$, the calculation forms

$$\bar A_b(E)=\sum_k w_k\sum_{t\in\{\mathrm{FRONT,BACK}\}}
A_t(E,\theta_k)\,P_t(E_{\rm measured}\in b\mid E,\theta_k).$$

`CurvatureResponse` integrates the actual curvature photon spectrum against $\bar A_b(E)$, over the full calibrated true-energy range. Photons can therefore migrate across measured-energy boundaries. The output has units **photons × cm²**; there is no exposure time or distance factor. This fast folding class currently targets curvature spectra, while the lower-level band-area function can weight another justified photon spectrum.

Energy dispersion uses the FITS scaling polynomial and mixture of two asymmetric exponential-power distributions described in the [official energy-dispersion documentation](https://fermi.gsfc.nasa.gov/ssc/data/analysis/lat_irfs/irf_energy_dispersion.html). Cumulative probabilities are integrated analytically. Interpolation mixes neighbouring **scaled probability distributions**, rather than interpolating their shape parameters; this follows the official [EdispInterpolator implementation](https://github.com/fermi-lat/irfs/blob/b7cafb0a7a99ec36925cf8984372b4caafb5039c/latResponse/src/EdispInterpolator.cxx) and [Edisp3 kernel](https://github.com/fermi-lat/irfs/blob/b7cafb0a7a99ec36925cf8984372b4caafb5039c/latResponse/src/Edisp3.cxx). The fitted distribution's tiny negative-energy tail is not silently reassigned to positive energies.

Finite detector energy resolution means a measured band receives photons from outside its true-energy interval, especially near spectral structure. This is why folding extends beyond the nominal band bounds. [Fermi analysis guidance](https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Pass8_edisp_usage.html)

## Incidence scenarios are assumptions, not observation histories

The caller must supply nonnegative incidence weights. The supplied alternatives are Gaussian quadrature for uniform distributions in $\cos\theta$:

| Scenario | Incidence cosine interval |
| --- | --- |
| `uniform_cosine` | 0.2–1.0 |
| `near_axis` | 0.8–1.0 |
| `off_axis` | 0.2–0.5 |

Weights sum to one. **None is measured Vela livetime**, and the spread between them is a sensitivity comparison, not a statistical confidence interval. The phase histograms alone do not recover the spacecraft pointing and livetime history.

Spatial/PSF acceptance, azimuthal response, livetime-efficiency corrections, event source-probability weights, background selection, and the actual observing history remain absent. In particular, published 3PC counts are weighted by source probabilities: multiplying a model by ordinary effective area cannot reproduce those weights. Separate fitted band normalizations also remove absolute and cross-band flux constraints. Any comparison with the existing profiles must therefore remain a conditional shape diagnostic under the chosen response scenario.

## Verification evidence

All **8 tests passed** on 2026-10-13:

- Downloaded calibration hashes, table units, direct grid values, FRONT+BACK addition, and energy-dependent effective area
- Asymmetric distribution CDF versus independent numerical integration and normalization
- Energy-band partition conservation, real migration across boundaries, and distribution-level interpolation
- Constant-area/no-dispersion limit against the independent emitted-photon calculation, within $5\times10^{-4}$
- Fast spectrum lookup versus direct folding within 0.2% for contributions above $10^{-8}$ of the maximum; doubling energy integration intervals changes the tested direct yields by less than 0.1%
- Explicit angle-weight normalization, domain checks, and zero/invalid-input behaviour

These tolerances verify the implementation on the tested cases; they do not estimate calibration systematics or establish agreement with Fermitools for a complete observation. Reproduce with `PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -p 'test_response.py' -v` from the project root.

## Observational status

All previously inspected 3PC profiles—including the earlier reserved energy bands—are now **development data**, because their discrepancies motivated the model redesign. Re-fitting or response-folding those profiles is not a fresh held-out test. An absolute observational analysis still needs consistent event selections, spacecraft livetime, source/background modelling, timing, and a likelihood or explicitly derived weighted-count response. Independent validation requires data or observables that do not guide further model selection.
