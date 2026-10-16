"""Frozen, conditional pulse-shape comparison to four disjoint Vela bands.

The tabulated weighted counts have no instrument response. Normalization and a
constant are therefore fitted separately in each band, including reserved bands.
This is a shape test conditional on those nuisance fits, not a flux/spectral test.
The phase origin is shared, and is chosen from the two fitting bands only.
"""
from numbers import Integral

import numpy as np
from astropy.io import fits
from scipy.optimize import minimize_scalar

from .observations import bin_average, nonnegative_scale_and_offset


# Freeze this protocol before optimization; aggregate/nested bands are excluded.
BANDS = (
    ((0.1, 0.3), '100_300_WtCt', 'fit'),
    ((0.3, 1.0), '300_1000_WtCt', 'holdout'),
    ((1.0, 3.0), '1000_3000_WtCt', 'fit'),
    ((3.0, 100.0), '3000_100000_WtCt', 'holdout'),
)
ERROR_CAVEAT = (
    'Supplied catalog variances are added when rebinning; this retains a '
    'maximum-weight correction per original bin and is conservative, not an '
    'exact event-level rebinned uncertainty.'
)
VALIDATION_CAVEAT = (
    'Conditional shape-only comparison: nonnegative scale and constant fitted '
    'in every band; one phase origin fitted from bands 0 and 2 only. No '
    'instrument response, absolute flux prediction, or calibrated significance '
    'after physical-parameter/model selection. Reserved bands were withheld '
    'from optimization, not from all earlier inspection of Vela morphology.'
)


def load_bands(path, rebin=8):
    """Load one independent cycle, preserving counts and catalog variances."""
    if isinstance(rebin, bool) or not isinstance(rebin, Integral) or rebin <= 0:
        raise ValueError('Rebin factor must be a positive integer')
    with fits.open(path, checksum=True) as hdul:
        table = hdul['GAMMA_LC'].data
        lo = np.array(table['Ph_Min'], dtype=float)
        hi = np.array(table['Ph_Max'], dtype=float)
        first = (lo >= 0) & (lo < 1) & (hi <= 1 + 1e-6)
        second = (lo >= 1) & (lo < 2) & (hi <= 2 + 1e-6)
        n = int(first.sum())
        if n == 0 or int(second.sum()) != n or n * 2 != len(lo):
            raise ValueError('Unexpected Fermi display-cycle structure')
        if not (np.allclose(lo[second] - 1, lo[first], atol=2e-7, rtol=0)
                and np.allclose(hi[second] - 1, hi[first], atol=2e-7, rtol=0)):
            raise ValueError('Repeated cycle does not have matching phase edges')
        lo, hi = lo[first], hi[first]
        if n % rebin:
            raise ValueError('Rebin factor must divide the independent bin count')
        if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))
                and np.all(hi > lo) and abs(lo[0]) < 1e-7
                and abs(hi[-1] - 1) < 1e-7
                and np.allclose(hi[:-1], lo[1:], atol=1e-7, rtol=0)):
            raise ValueError('Invalid or noncontiguous phase bins')
        lower = lo.reshape(-1, rebin)[:, 0]
        upper = hi.reshape(-1, rebin)[:, -1]
        result = []
        for index, (energy, column, role) in enumerate(BANDS):
            raw = np.array(table[column], dtype=float)
            err = np.array(table['Unc_' + column], dtype=float)
            if not (np.array_equal(raw[first], raw[second])
                    and np.array_equal(err[first], err[second])):
                raise ValueError(f'Nonidentical display cycles in {column}')
            y, sigma = raw[first], err[first]
            if not (np.all(np.isfinite(y)) and np.all(y >= 0)
                    and np.all(np.isfinite(sigma)) and np.all(sigma >= 0)):
                raise ValueError(f'Invalid count/error data in {column}')
            y = y.reshape(-1, rebin).sum(axis=1)
            sigma = np.sqrt((sigma.reshape(-1, rebin) ** 2).sum(axis=1))
            if np.any(sigma <= 0):
                raise ValueError(f'Zero rebinned uncertainty in {column}; no floor is applied')
            centers = (lower + upper) / 2
            result.append({
                'index': index, 'energies_GeV': energy, 'role': role,
                'band': f'{energy[0]:g}–{energy[1]:g} GeV',
                'count_column': column, 'phase_lo': lower.copy(),
                'phase_hi': upper.copy(), 'phase': centers.copy(),
                'phase_centers': centers.copy(), 'counts': y, 'error': sigma,
                'original_independent_bins': n, 'fit_bins': len(y),
                'duplicate_display_cycle_removed': True,
                'rebin_factor': int(rebin), 'error_caveat': ERROR_CAVEAT,
            })
    return result


def _checked_profiles(phase, profiles, data, observation_indices=None):
    phase = np.asarray(phase, dtype=float)
    profiles = np.asarray(profiles, dtype=float)
    if (phase.ndim != 1 or len(phase) < 3 or not np.all(np.isfinite(phase))
            or np.any(np.diff(phase) <= 0) or phase[0] < 0 or phase[-1] >= 1):
        raise ValueError('Model phases must increase strictly within [0,1)')
    if (profiles.shape != (len(data), len(phase))
            or not np.all(np.isfinite(profiles)) or np.any(profiles < 0)):
        raise ValueError('Profiles must be finite nonnegative arrays [band,phase]')
    checked_data = data if observation_indices is None else [data[i] for i in observation_indices]
    for band in checked_data:
        lo, hi, y, err = (np.asarray(band[key], dtype=float)
                          for key in ('phase_lo', 'phase_hi', 'counts', 'error'))
        if (lo.ndim != 1 or len(lo) < 3 or not (lo.shape == hi.shape == y.shape == err.shape)
                or not np.all(np.isfinite([lo, hi, y, err]))
                or np.any(hi <= lo) or np.any(err <= 0) or np.any(y < 0)):
            raise ValueError('Invalid observation arrays or nonpositive errors')
    return phase, profiles


def _score_one(phase, profile, band, shift):
    model = bin_average(phase, profile, band['phase_lo'], band['phase_hi'], shift)
    y, err = np.asarray(band['counts']), np.asarray(band['error'])
    # Scaling does not affect a shape fit and avoids under/overflow for tiny
    # physical emissivities. An exactly dark model has only its constant term.
    maximum = float(np.max(model))
    if maximum > 0:
        normalized = model / maximum
        coef = nonnegative_scale_and_offset(normalized, y, err)
        prediction = coef[0] * normalized + coef[1]
        scale = float(coef[0] / maximum)
    else:
        weight = 1 / err ** 2
        coef = np.array([0., max(0., np.sum(weight * y) / np.sum(weight))])
        prediction = np.full_like(y, coef[1], dtype=float)
        scale = 0.
    residual = (y - prediction) / err
    chi2 = float(residual @ residual)
    dof = len(y) - 2
    return {
        'band': band.get('band', ''), 'energies_GeV': band.get('energies_GeV'),
        'role': band.get('role', ''), 'phase_shift_cycles': float(shift % 1),
        'scale': scale, 'constant_component': float(coef[1]),
        'prediction': prediction, 'residual_sigma': residual, 'chi_squared': chi2,
        'descriptive_degrees_of_freedom': dof, 'degrees_of_freedom': dof,
        'reduced_chi_squared': chi2 / dof,
        'nuisance_parameters_fitted': 2, 'validation_kind': 'conditional shape only',
        'statistical_caveat': VALIDATION_CAVEAT,
    }


def score_profiles(phase, profiles, data, shift):
    """Fit only nonnegative A,C per band for A*profile(phi-shift)+C.

    Returns one dictionary per input band; no geometry or phase is optimized.
    Per-band N-2 scores are descriptive; they omit shared fitted physical/phase
    parameters and must not be interpreted as an exact goodness-of-fit law.
    """
    phase, profiles = _checked_profiles(phase, profiles, data)
    if not np.isfinite(shift):
        raise ValueError('Phase shift must be finite')
    return [_score_one(phase, profile, band, float(shift))
            for profile, band in zip(profiles, data)]


def fit_shared_phase(phase, profiles, data, n_shifts=256, evaluate_holdout=True):
    """Choose one phase using fit bands only, then score the frozen holdouts.

    A periodic grid search followed by local refinement covers phase aliases.
    The energy split is immutable here: roles must match BANDS in their order.
    Set evaluate_holdout=False throughout physical-parameter/model selection;
    reserved observations are then neither inspected nor scored. Evaluate them
    only after selecting and freezing the physical model and shared phase.
    """
    if (len(data) != len(BANDS)
            or [d.get('role') for d in data] != [b[2] for b in BANDS]):
        raise ValueError('Expected frozen band order with fit,holdout,fit,holdout roles')
    if isinstance(n_shifts, bool) or not isinstance(n_shifts, Integral) or n_shifts < 16:
        raise ValueError('Use at least 16 phase-grid points')
    train = [i for i, band in enumerate(data) if band['role'] == 'fit']
    phase, profiles = _checked_profiles(phase, profiles, data,
                                       None if evaluate_holdout else train)

    def objective(shift):
        return sum(_score_one(phase, profiles[i], data[i], shift)['chi_squared']
                   for i in train)

    grid = np.arange(n_shifts) / n_shifts
    values = np.array([objective(shift) for shift in grid])
    best_index = int(np.argmin(values))
    best_shift, best_value = float(grid[best_index]), float(values[best_index])
    minima = np.flatnonzero((values <= np.roll(values, 1))
                            & (values <= np.roll(values, -1)))
    # Refine up to four promising aliases, using a wraparound objective.
    for i in minima[np.argsort(values[minima])[:4]]:
        fit = minimize_scalar(objective, bounds=(grid[i] - 1 / n_shifts,
                                                 grid[i] + 1 / n_shifts),
                              method='bounded', options={'xatol': 1e-10})
        if fit.fun < best_value:
            best_shift, best_value = float(fit.x % 1), float(fit.fun)
    evaluate = range(len(data)) if evaluate_holdout else train
    scores = [_score_one(phase, profiles[i], data[i], best_shift) for i in evaluate]
    return {
        'phase_shift_cycles': best_shift,
        'training_chi_squared': sum(s['chi_squared'] for s in scores if s['role'] == 'fit'),
        'holdout_chi_squared': (sum(s['chi_squared'] for s in scores if s['role'] == 'holdout')
                               if evaluate_holdout else None),
        'bands': scores, 'phase_grid_points': int(n_shifts),
        'training_band_indices': [0, 2], 'holdout_band_indices': [1, 3],
        'holdout_phase_refitted': False, 'holdout_evaluated': bool(evaluate_holdout),
        'validation_kind': 'conditional shape only',
        'statistical_caveat': VALIDATION_CAVEAT,
    }
