# Vela reference observations

Source: NASA Fermi Science Support Center, third gamma-ray pulsar catalog (Smith et al. 2023).

- [Catalog and data](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/)
- [Profile format](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/README_3PC_ProfileData.txt)
- [Corrected archive, 18 July 2025](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/3PC_ProfileData_20250718.tgz)

The archive was retrieved on 12 October 2026. Only its J0835-4510 FITS and ASCII files are retained. These files are unmodified; provenance.json records the source URL and SHA-256 hashes. Calculations use FITS, which avoids ambiguity from printed phase rounding.

GAMMA_LC contains two identical display cycles. Use the first 800 bins (phase 0 to 1), not all 1600 as independent measurements. GT100_WtCnt and Unc_GT100_WtCnt are the supplied >100 MeV weighted counts and uncertainties. They are not raw integer photon counts or calibrated energy flux. The pipeline combines adjacent groups of four, summing counts and variances, before its diagnostic shape comparison.

The model fits a constant component as a nuisance term; the data here are not claimed to be background subtracted. BEST_FIT_LC and RADIO_PROFILE are retained in the original FITS file but are not used to tune this model.
