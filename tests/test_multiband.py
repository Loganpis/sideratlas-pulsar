"""Independent conservation and phase-recovery tests for the frozen protocol."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
from astropy.io import fits

from pulsar.multiband import BANDS, fit_shared_phase, load_bands, score_profiles

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/J0835-4510_3PC_data.fits'


class MultiBandChecks(unittest.TestCase):
    def test_disjoint_protocol_count_and_variance_conservation(self):
        bands = load_bands(DATA)
        self.assertEqual([b['role'] for b in bands], ['fit', 'holdout', 'fit', 'holdout'])
        self.assertEqual([b['energies_GeV'] for b in bands],
                         [(0.1, 0.3), (0.3, 1.), (1., 3.), (3., 100.)])
        with fits.open(DATA) as hdul:
            table = hdul['GAMMA_LC'].data
            for band in bands:
                column = band['count_column']
                raw = np.asarray(table[column][:800], dtype=float)
                error = np.asarray(table['Unc_' + column][:800], dtype=float)
                self.assertEqual(len(band['counts']), 100)
                self.assertTrue(band['duplicate_display_cycle_removed'])
                self.assertAlmostEqual(band['counts'].sum(), raw.sum(), places=6)
                self.assertAlmostEqual((band['error'] ** 2).sum(), (error ** 2).sum(), places=6)
                self.assertTrue(np.all(band['error'] > 0))
            # Distinct intervals account for the aggregate counts, allowing the
            # independent float32 rounding of the stored aggregate column.
            total = sum(b['counts'].sum() for b in bands)
            aggregate = np.asarray(table['GT100_WtCnt'][:800], dtype=float).sum()
            self.assertLess(abs(total - aggregate), 0.1)

    def test_changed_display_cycle_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'bad.fits'
            with fits.open(DATA) as hdul:
                hdul['GAMMA_LC'].data['100_300_WtCt'][900] += 1
                hdul.writeto(path, checksum=True)
            with self.assertRaisesRegex(ValueError, 'Nonidentical display cycles'):
                load_bands(path)

    def test_bad_rebin_is_rejected(self):
        for factor in (0, -1, 3, 1.5, True):
            with self.subTest(factor=factor), self.assertRaises(ValueError):
                load_bands(DATA, factor)

    @staticmethod
    def synthetic_data(shift):
        phase = np.arange(4096) / 4096
        edges = np.linspace(0, 1, 81)
        lo, hi = edges[:-1], edges[1:]
        profiles, bands = [], []
        for i, (energies, _, role) in enumerate(BANDS):
            a, b = 0.3 + i * .1, .12 + i * .015
            profile = 1 + a * np.cos(2 * np.pi * phase) + b * np.sin(4 * np.pi * phase)
            # Analytical integration independent of the profile interpolator.
            def primitive(x):
                return (x + a * np.sin(2 * np.pi * (x - shift)) / (2 * np.pi)
                        - b * np.cos(4 * np.pi * (x - shift)) / (4 * np.pi))
            average = (primitive(hi) - primitive(lo)) / (hi - lo)
            profiles.append(profile)
            bands.append({'energies_GeV': energies, 'role': role,
                          'phase_lo': lo, 'phase_hi': hi,
                          'counts': (i + 2) * average + (i + 4),
                          'error': np.ones(len(lo)) * .1})
        return phase, np.array(profiles), bands

    def test_one_shared_phase_recovered_with_independent_analytic_integrals(self):
        shift = .987654  # Also exercises local search across phase zero.
        phase, profiles, bands = self.synthetic_data(shift)
        result = fit_shared_phase(phase, profiles, bands, n_shifts=64)
        self.assertAlmostEqual(result['phase_shift_cycles'], shift, places=6)
        self.assertLess(result['training_chi_squared'], 1e-6)
        self.assertLess(result['holdout_chi_squared'], 1e-6)
        self.assertFalse(result['holdout_phase_refitted'])
        for i, score in enumerate(result['bands']):
            self.assertAlmostEqual(score['scale'], i + 2, places=4)
            self.assertAlmostEqual(score['constant_component'], i + 4, places=4)

    def test_holdout_counts_cannot_change_selected_phase(self):
        phase, profiles, bands = self.synthetic_data(.173912)
        before = fit_shared_phase(phase, profiles, bands, n_shifts=64)
        # Deliberately incompatible reserved observations must worsen the score
        # without pulling the shared phase to accommodate them.
        for i in (1, 3):
            bands[i]['counts'] = np.roll(bands[i]['counts'], 17)
        after = fit_shared_phase(phase, profiles, bands, n_shifts=64)
        self.assertEqual(after['phase_shift_cycles'], before['phase_shift_cycles'])
        self.assertEqual(after['training_chi_squared'], before['training_chi_squared'])
        self.assertGreater(after['holdout_chi_squared'], 100)

    def test_dark_model_has_finite_constant_only_fit(self):
        phase, profiles, bands = self.synthetic_data(.1)
        profiles[:] = 0
        for band in bands:
            band['counts'][:] = 7
        scores = score_profiles(phase, profiles, bands, .2)
        for score in scores:
            self.assertEqual(score['scale'], 0)
            self.assertAlmostEqual(score['constant_component'], 7)
            self.assertLess(score['chi_squared'], 1e-20)

    def test_training_search_never_inspects_reserved_counts(self):
        phase, profiles, bands = self.synthetic_data(.271234)
        before = fit_shared_phase(phase, profiles, bands, n_shifts=64,
                                  evaluate_holdout=False)
        for i in (1, 3):
            bands[i]['counts'][:] = np.nan
            bands[i]['error'][:] = 0
        after = fit_shared_phase(phase, profiles, bands, n_shifts=64,
                                 evaluate_holdout=False)
        self.assertEqual(after['phase_shift_cycles'], before['phase_shift_cycles'])
        self.assertEqual(after['training_chi_squared'], before['training_chi_squared'])
        self.assertIsNone(after['holdout_chi_squared'])
        self.assertFalse(after['holdout_evaluated'])
        self.assertEqual([b['role'] for b in after['bands']], ['fit', 'fit'])


if __name__ == '__main__':
    unittest.main()
