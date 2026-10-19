# Vela multi-band data and validation audit

Audit date: 13 October 2026. This is a proposed validation protocol and a direct inspection of the existing data, not a successful Vela reconstruction.

## Data actually available locally

The retained `data/J0835-4510_3PC_data.fits` contains `GAMMA_LC` (1,600 rows, 16 columns), `BEST_FIT_LC`, and `RADIO_PROFILE`. In every gamma-ray count/error column the second 800 rows repeat the first cycle exactly. Phase edges shift by one cycle. Fit only the first cycle, with edges from 0 to 1 and width 0.00125.

Directly measured totals and empty-bin counts in that first cycle:

| Energy interval | Count column | Total weighted counts | Raw bins with zero error |
| --- | --- | ---: | ---: |
| 50–100 MeV | `50_100_WtCt` | 324,578.55 | 0 |
| 100–300 MeV | `100_300_WtCt` | 1,274,440.36 | 0 |
| 300–1,000 MeV | `300_1000_WtCt` | 1,215,640.96 | 0 |
| 1–3 GeV | `1000_3000_WtCt` | 500,686.99 | 0 |
| 3–100 GeV | `3000_100000_WtCt` | 115,099.84 | 0 |
| 10–100 GeV | `10000_100000_WtCt` | 6,035.98 | 41 |
| Aggregate >100 MeV | `GT100_WtCnt` | 3,105,868.15 | 0 |

The four non-overlapping bands from 100 MeV to 100 GeV sum to the aggregate count profile to floating-point precision (maximum bin difference 0.001465 weighted counts). Do not fit these bands and then describe the aggregate profile as held out. Similarly, 10–100 GeV is contained in 3–100 GeV and cannot validate a fit that already used the latter independently.

The catalog [profile format documentation](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/README_3PC_ProfileData.txt) establishes the column names and extension meanings. The file has no effective-area/exposure table, energy dispersion matrix, event list, background covariance, or joint parameter covariance. These omissions are established by the local header inspection.

## Supplied error bars are not ordinary Poisson counts

[Smith et al. (2023), section 5.1](https://arxiv.org/pdf/2307.11132) defines each bin's error as

\[
\sigma_i^2=\sum_{j\in i}w_j^2+\left(\max_{j\in i} w_j\right)^2.
\]

The maximum-weight term is a small-count correction. The same paper derives a background estimate from sums of weights and squared weights, includes diffuse-background systematics, and describes the empirical pulse templates in section 5.2. `BEST_FIT_LC` is a fit to these observations, not an independent physical prediction.

Consequences for this implementation:

- Fractional weighted counts are not integer Poisson counts; using a Poisson likelihood on the tabulated totals is unjustified.
- Adding supplied variances while combining bins retains one correction per original bin. It is conservative relative to applying the same prescription once to the combined event sample, but is not an exact event-level rebinning.
- The aggregate variance is not the sum of the four band variances because the correction is not additive. In this file their difference ranges from 1.9636 to 2.9584 per raw bin. This does not itself demonstrate correlated noise.
- The unrecorded maximum weights prevent exact reconstruction of the background from count and uncertainty columns alone.
- Disjoint bands use distinct detected photons, conditional on the common source/weight model. Calibration, timing, source-model, and background uncertainties remain shared; treating the bands as independent measurements does not remove these systematic uncertainties.

Rebinning four adjacent bins gives 200 bins per cycle. Every main band still has positive uncertainty. In 3–100 GeV the smallest rebinned weighted count is 13.16 and the minimum diagnostic effective sample size `(sum w)^2 / supplied variance` is 16.60. Eight-bin rebinning gives 100 bins and a minimum effective size of 36.83 there. These are diagnostics, not proofs of a Gaussian likelihood. The 10–100 GeV profile remains sparse even after 16-bin rebinning; retain it for a separate event-level analysis rather than adding an arbitrary error floor.

## Recommended frozen test protocol

1. Use 100–300 MeV and 1–3 GeV for parameter fitting. Reserve 300–1,000 MeV for energy interpolation and 3–100 GeV for high-energy extrapolation. Fix this split before optimization. The previously inspected aggregate Vela morphology means this is a withheld-from-optimization test, not a blind discovery experiment.
2. Fit a shared geometry, phase origin, acceleration prescription, particle injection law, and radiation parameters. Derive peaks from the transport and photon calculation. Do not provide a free brightness coefficient at every observed phase or import the catalog's empirical template into emissivity.
3. Freeze all physical parameters and phase origin before evaluating reserved bands. No independent phase shift, cutoff, sheet width, or peak ratio may be tuned there.
4. Because the current files do not specify the response needed to predict absolute weighted counts, profile out at most a nonnegative normalization and constant separately in each band. Label this **conditional pulse-shape testing**, not spectral or absolute-flux validation. A held-out band's two fitted nuisance parameters still consume information: report its score with `N_bin - 2` as a descriptive denominator, not as an exact calibrated chi-squared distribution after nonlinear model selection.
5. Compare residuals, peak positions and widths, peak ratio, bridge emission, and energy evolution. Report scores separately for fitting and reserved data. Test numerical convergence on the fitted parameters. Do not enlarge error bars merely until the fit becomes acceptable.
6. A failure remains a failure even if the image looks persuasive. If model design changes after seeing held-out residuals, those bands have become development data; record this and obtain a new validation split before claiming predictive performance.

Broad-band shapes themselves depend on energy-dependent detector weighting when the source changes shape with energy. Numerically integrate the photon-number spectrum over each interval and disclose any assumed response. A representative photon energy or an assumed flat response is insufficient to claim a measured phase-resolved spectrum. At best it supports a conditional reduced-model comparison and sensitivity checks.

## Stronger next dataset

The official [photon archive](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/catalogs/3PC/photon/3deg_50MeV/) lists Vela's [individual-event FITS file](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/catalogs/3PC/photon/3deg_50MeV/J0835-4510_3p0deg_50MeV_300000MeV_105deg_128_3.fits), about 428 MB. It has not been downloaded in this audit. The [archive README](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/catalogs/3PC/photon/README) documents photon phase and model-weight columns and the 3-degree, 50 MeV–300 GeV selection. This would permit predefined event splits, a weighted phase likelihood, and exact rebinning. A time split must check valid timing intervals and pulsar glitches rather than assume stationarity.

Even an excellent multi-band profile match cannot uniquely reconstruct a three-dimensional magnetosphere from one viewing direction. A defensible eventual claim is **an observation-constrained Vela emission model**, with its model family, fitted parameters, held-out tests, response treatment, degeneracies, and unresolved physics stated. A self-consistent kinetic reconstruction additionally requires plasma/field evolution and energy accounting, which an emission fit to prescribed fields does not establish.
