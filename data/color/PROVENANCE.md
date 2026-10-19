# CIE 1931 color-matching functions

Creator and rights holder: International Commission on Illumination (CIE)

Citation: CIE 2019, *Colour-matching functions of CIE 1931 standard colorimetric observer*, International Commission on Illumination (CIE), Vienna, AT, DOI [10.25039/CIE.DS.xvudnb9b](https://doi.org/10.25039/CIE.DS.xvudnb9b)

Official dataset landing page: https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer

Downloaded 2026-10-17 from the exact official links on that page:

- CSV: https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv
- Metadata: https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv_metadata.json

The CSV and metadata are redistributed **unmodified**. CSV SHA-256: `fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1`; MD5: `17cca777db64b17170f06f67ce9d3ab7`. Both match the official metadata. The CSV has 471 rows, 360–830 nm in 1 nm steps, with wavelength, x̄, ȳ and z̄ columns. Original source: CIE 018:2019, Table 6

## License and attribution

The official metadata identifies **Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)**: https://creativecommons.org/licenses/by-sa/4.0/

Retain this attribution, the dataset DOI, this license link and an indication of any modifications when sharing the dataset. If sharing an adaptation of this dataset, use the same license as required by its terms. The official metadata file is retained alongside the original data. These data have their own license; no claim is made that the CIE dataset is covered by the software repository's code license. CIE does not endorse this project

## Numerical use and scope

`src/pulsar/optical_color.py` uses linear CMF interpolation, as specified by the CIE metadata. The module rejects wavelengths outside the table and requires spectra to cover the entire **400–700 nm** output band. It analytically integrates products of piecewise-linear spectral radiance and CMFs. Per-metre and per-nanometre spectral densities are converted with the wavelength Jacobian

XYZ is photometrically scaled by 683 lm/W. The D65-referenced XYZ-to-linear-sRGB matrix and sRGB transfer function follow the primary W3C specification's [color-conversion section](https://www.w3.org/TR/css-color-4/#color-conversion-code). No source-specific white balance is applied. The same explicit scalar display exposure must be used for every pixel and frame. Negative or excessive linear RGB channels are clipped, with the unclipped values and clipping masks retained; clipping can alter spectral hue/saturation

400–700 nm excludes some human-visible wavelengths. This is a stated display band, not a model of the full human visual system, dark adaptation, a telescope, a detector, Vela's atmosphere or Vela's optical emission
