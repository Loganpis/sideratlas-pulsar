"""Real LAT effective-area/energy-dispersion sensitivity, not Vela exposure.

The supplied official P8R3_SOURCE_V3 FRONT/BACK tables describe instrument
response. Incidence-angle weights MUST be supplied: they are not recoverable
from the 3PC weighted phase histograms. This module omits spatial acceptance,
catalog source-probability weights, azimuth and live-time efficiency corrections.
Its result is area-weighted photon yield, NOT predicted 3PC weighted counts.
"""
from dataclasses import dataclass
from pathlib import Path
from numbers import Integral

import numpy as np
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator
from scipy.special import gammainc, gammaincc

from .radiation import curvature_kernel, GEV_ERG, KERNEL_POWER_INTEGRAL


_DATA = Path(__file__).resolve().parents[2] / 'data' / 'response'


def _bands(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 2 or values.shape[1] != 2
            or np.any(np.isnan(values)) or np.any(values[:, 0] < 0)
            or np.any(values[:, 1] <= values[:, 0])):
        raise ValueError('Require energy bands [minimum,maximum], 0 <= minimum < maximum')
    return values


def _incidence(cos_theta, weights):
    mu, weights = np.asarray(cos_theta, float), np.asarray(weights, float)
    if (mu.ndim != 1 or weights.shape != mu.shape or len(mu) == 0
            or np.any(~np.isfinite(mu)) or np.any(~np.isfinite(weights))
            or np.any(mu < .2) or np.any(mu > 1.) or np.any(weights < 0)
            or weights.sum() <= 0):
        raise ValueError('Supply nonnegative incidence weights and 0.2 <= cos(theta) <= 1')
    return mu, weights / weights.sum()


def incidence_scenarios(n=16):
    """Explicit illustrative angular averages; none is measured Vela livetime.

    Uniform d(cos theta) over [0.2,1], near-axis [0.8,1], off-axis [0.2,0.5].
    These alternatives measure sensitivity, not a confidence interval.
    """
    if isinstance(n, bool) or not isinstance(n, Integral) or n < 2:
        raise ValueError('Use at least two inclination quadrature nodes')
    x, w = np.polynomial.legendre.leggauss(n)
    return {name: ((lo + hi) / 2 + (hi - lo) / 2 * x, w / w.sum())
            for name, lo, hi in [('uniform_cosine', .2, 1.),
                                ('near_axis', .8, 1.), ('off_axis', .2, .5)]}


def _cdf_asymmetric(x, sigma, skew, bias, power):
    """CDF of normalized asymmetric exponential-power kernel (LAT Pass 8).

    Left mass k^2/(1+k^2), right mass 1/(1+k^2). Written using gammaincc
    in the lower tail to avoid subtracting nearly equal numbers.
    """
    u = x - bias
    left_mass = skew**2 / (1 + skew**2)
    left = left_mass * gammaincc(1 / power, (np.maximum(-u, 0) / (skew * sigma))**power)
    right = left_mass + (1 - left_mass) * gammainc(1 / power, (skew * np.maximum(u, 0) / sigma)**power)
    return np.where(u < 0, left, right)


def _mixture_cdf(x, parameters):
    f, s1, k1, b1, b2, s2, k2, p1, p2 = np.moveaxis(parameters, -1, 0)
    return f * _cdf_asymmetric(x, s1, k1, b1, p1) + (1 - f) * _cdf_asymmetric(x, s2, k2, b2, p2)


def _brackets(grid, query):
    q = np.clip(query, grid[0], grid[-1])
    i = np.clip(np.searchsorted(grid, q, side='right') - 1, 0, len(grid) - 2)
    t = (q - grid[i]) / (grid[i + 1] - grid[i])
    return i, t


@dataclass
class LATResponse:
    """Official FRONT/BACK calibration tables, azimuth-averaged sensitivity."""
    area: dict
    dispersion: dict
    energy_limits_gev: tuple
    irf: str = 'P8R3_SOURCE_V3'

    def effective_area_cm2(self, energy_gev, cos_theta, event_type=None):
        energy, mu = np.broadcast_arrays(np.asarray(energy_gev, float), np.asarray(cos_theta, float))
        self._check_domain(energy, mu)
        parts = ('FRONT', 'BACK') if event_type is None else (event_type,)
        total = np.zeros_like(energy)
        for part in parts:
            loge, costh, values = self.area[part]
            xy = np.column_stack((np.clip(mu.ravel(), costh[0], costh[-1]),
                                  np.clip(np.log10(energy.ravel()), loge[0], loge[-1])))
            total += RegularGridInterpolator((costh, loge), values)(xy).reshape(energy.shape)
        return total

    def _check_domain(self, energy, mu):
        lo, hi = self.energy_limits_gev
        if (np.any(~np.isfinite(energy)) or np.any(energy < lo) or np.any(energy > hi)
                or np.any(~np.isfinite(mu)) or np.any(mu < .2) or np.any(mu > 1)):
            raise ValueError(f'Response requires {lo:g} <= true energy <= {hi:g} GeV and 0.2 <= cos(theta) <= 1')

    def band_probability(self, energy_gev, cos_theta, bands_gev, event_type):
        """P(measured energy in band | true energy, inclination, event type).

        Bilinear interpolation is on the scaled probability distributions,
        NOT the fitted shape parameters, following the official EdispInterpolator.
        The parametrization is normalized on all real measured energies; its
        tiny unphysical negative-energy tail is not silently redistributed.
        """
        bands = _bands(bands_gev)
        energy, mu = np.broadcast_arrays(np.asarray(energy_gev, float), np.asarray(cos_theta, float))
        self._check_domain(energy, mu)
        loge, costh, pars, coeff = self.dispersion[event_type]
        lg = np.log10(energy)
        le_mev = lg + 3
        scale = (coeff[0]*le_mev**2 + coeff[1]*mu**2 + coeff[2]*le_mev
                 + coeff[3]*mu + coeff[4]*le_mev*mu + coeff[5])
        if np.any(scale <= 0):
            raise ValueError('Nonpositive energy-dispersion scale')
        ie, te = _brackets(loge, lg)
        im, tm = _brackets(costh, mu)
        xlo = (bands[:, 0] / energy[..., None] - 1) / scale[..., None]
        xhi = (bands[:, 1] / energy[..., None] - 1) / scale[..., None]
        result = np.zeros(energy.shape + (len(bands),))
        for de, we in ((0, 1-te), (1, te)):
            for dm, wm in ((0, 1-tm), (1, tm)):
                p = pars[im+dm, ie+de][..., None, :]
                result += (we * wm)[..., None] * (_mixture_cdf(xhi, p) - _mixture_cdf(xlo, p))
        return np.clip(result, 0, 1)

    def band_area_cm2(self, energy_gev, bands_gev, cos_theta, weights, energy_dispersion=True):
        """Average sum_type A(E,theta,type)*P(E_meas in band) over explicit weights."""
        mu, weights = _incidence(cos_theta, weights)
        bands = _bands(bands_gev)
        energy = np.asarray(energy_gev, float)
        result = np.zeros(energy.shape + (len(bands),))
        for m, w in zip(mu, weights):
            for part in ('FRONT', 'BACK'):
                area = self.effective_area_cm2(energy, m, part)
                if energy_dispersion:
                    probability = self.band_probability(energy, m, bands, part)
                else:
                    probability = ((energy[..., None] >= bands[:, 0]) &
                                   (energy[..., None] < bands[:, 1])).astype(float)
                result += w * area[..., None] * probability
        return result


def load_lat_response(directory=None):
    """Load the version-pinned local official FITS calibration; no network I/O."""
    directory = _DATA if directory is None else Path(directory)
    area, dispersion = {}, {}
    with fits.open(directory / 'aeff_P8R3_SOURCE_V3_FB.fits', checksum=True) as af, fits.open(
            directory / 'edisp_P8R3_SOURCE_V3_FB.fits', checksum=True) as df:
        for part in ('FRONT', 'BACK'):
            h = af['EFFECTIVE AREA_' + part]
            if h.columns['EFFAREA'].unit != 'm2' or h.columns['ENERG_LO'].unit != 'MeV':
                raise ValueError('Unexpected effective-area units')
            if h.header.get('CBD10001') != 'VERSION(P8R3_SOURCE_V3)':
                raise ValueError('Unexpected response version')
            t = h.data[0]
            elo, ehi = np.asarray(t['ENERG_LO'], float)/1000, np.asarray(t['ENERG_HI'], float)/1000
            energy_limits = (float(elo[0]), float(ehi[-1]))
            loge = (np.log10(elo)+np.log10(ehi))/2
            mu = (np.asarray(t['CTHETA_LO'], float)+np.asarray(t['CTHETA_HI'], float))/2
            area[part] = (loge, mu, np.asarray(t['EFFAREA'], float)*10000)
            t = df['ENERGY DISPERSION_' + part].data[0]
            loge = (np.log10(np.asarray(t['ENERG_LO'], float)) +
                    np.log10(np.asarray(t['ENERG_HI'], float)))/2 - 3
            mu = (np.asarray(t['CTHETA_LO'], float)+np.asarray(t['CTHETA_HI'], float))/2
            names = ('F', 'S1', 'K1', 'BIAS', 'BIAS2', 'S2', 'K2', 'PINDEX1', 'PINDEX2')
            pars = np.stack([np.asarray(t[n], float) for n in names], axis=-1)
            scale = np.asarray(df['EDISP_SCALING_PARAMS_' + part].data['EDISPSCALE'][0], float)
            dispersion[part] = (loge, mu, pars, scale)
    return LATResponse(area, dispersion, energy_limits)


def _log_quadrature(lo, hi, breakpoints=(), intervals=256, order=6):
    if intervals < 4 or order < 2:
        raise ValueError('Insufficient energy quadrature')
    edges = np.unique(np.r_[np.linspace(np.log(lo), np.log(hi), intervals+1),
                            np.log([v for v in breakpoints if lo < v < hi])])
    x, w = np.polynomial.legendre.leggauss(order)
    middle, half = (edges[:-1]+edges[1:])/2, np.diff(edges)/2
    return np.exp((middle[:, None]+half[:, None]*x).ravel()), (half[:, None]*w).ravel()


class CurvatureResponse:
    """Fast real-response curvature spectrum fold, output photons * cm².

    There is no distance, exposure, flux normalization, PSF selection, or
    event-probability weighting in this quantity. The lookup is for one
    monoenergetic curvature spectrum per input segment; use spectral nodes
    if gamma varies significantly during a segment. Zero critical energy
    returns zero. No extra Doppler transform belongs here.
    """
    def __init__(self, response, bands_gev, cos_theta, weights,
                 energy_dispersion=True, intervals=256, order=6,
                 critical_energy_nodes=2401):
        self.response = response
        self.bands = _bands(bands_gev)
        self.cos_theta, self.weights = _incidence(cos_theta, weights)
        self.energy_dispersion = bool(energy_dispersion)
        if critical_energy_nodes < 100:
            raise ValueError('Use at least 100 critical-energy lookup nodes')
        lo, hi = response.energy_limits_gev
        # Edisp moves events over measured band edges: integrate the entire
        # calibrated true-energy range rather than truncating at those edges.
        self.energy, self.log_weights = _log_quadrature(lo, hi, self.bands.ravel(), intervals, order)
        self.band_area = response.band_area_cm2(self.energy, self.bands,
                                               self.cos_theta, self.weights, self.energy_dispersion)
        self.log_ec = np.linspace(np.log(1e-7), np.log(1e6), critical_energy_nodes)
        self.ec = np.exp(self.log_ec)
        self.per_erg = self.direct_yield(np.ones_like(self.ec), self.ec)

    def direct_yield(self, radiated_erg, critical_gev):
        """Direct energy quadrature, useful for independent interpolation checks."""
        energy, ec = self._check_inputs(radiated_erg, critical_gev)
        result = np.zeros(energy.shape + (len(self.bands),))
        er, cr, out = energy.ravel(), ec.ravel(), result.reshape(-1, len(self.bands))
        response_weights = self.band_area * self.log_weights[:, None]
        for lo in range(0, er.size, 256):
            hi = min(er.size, lo+256)
            ce = cr[lo:hi]
            safe = np.where(ce > 0, ce, 1.)
            spectra = curvature_kernel(self.energy[None, :] / safe[:, None])
            factors = np.where(ce > 0, er[lo:hi] / safe / GEV_ERG / KERNEL_POWER_INTEGRAL, 0)
            out[lo:hi] = (spectra @ response_weights) * factors[:, None]
        return result

    @staticmethod
    def _check_inputs(radiated_erg, critical_gev):
        energy, ec = np.broadcast_arrays(np.asarray(radiated_erg, float), np.asarray(critical_gev, float))
        if np.any(~np.isfinite(energy)) or np.any(energy < 0) or np.any(~np.isfinite(ec)) or np.any(ec < 0):
            raise ValueError('Require finite nonnegative emitted and critical energies')
        return energy, ec

    def photon_yield(self, radiated_erg, critical_gev):
        energy, ec = self._check_inputs(radiated_erg, critical_gev)
        if np.any(ec > self.ec[-1]):
            raise ValueError('Critical energy exceeds response lookup domain')
        out = np.zeros(ec.shape + (len(self.bands),))
        use = ec >= self.ec[0]
        for i in range(len(self.bands)):
            values = self.per_erg[:, i]
            log_values = np.log(np.maximum(values, np.finfo(float).tiny))
            interp = np.exp(np.interp(np.log(ec[use]), self.log_ec, log_values))
            out[..., i][use] = interp * energy[use]
        return out

