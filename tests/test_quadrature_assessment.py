"""Assessments must retain physical phase despite any fitted data alignment."""
from pathlib import Path
import sys
import unittest
import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from assess_emission_quadrature import compare_runs,peaks,comparison_invariants


class AssessmentChecks(unittest.TestCase):
    def test_physical_changes_are_not_quadrature_refinements(self):
        a={'snapshot':{'sha256':'field'},'cap':{'azimuths':128,'rings':6,'seed_radius':.7},
           'transport':{'dt':.004,'r_max':2},'observer_deg':65,
           'acceleration':{'rate':.25},'geometry_cache':'one','source_hashes':{'code':'a'}}
        b={**a,'cap':{**a['cap'],'azimuths':256,'rings':12},
           'transport':{**a['transport'],'dt':.002},'geometry_cache':'two'}
        self.assertEqual(comparison_invariants(a),comparison_invariants(b))
        self.assertNotEqual(comparison_invariants(a),comparison_invariants({**b,'observer_deg':66}))
        self.assertNotEqual(comparison_invariants(a),comparison_invariants({**b,'acceleration':{'rate':.5}}))

    def test_comparison_preserves_unshifted_peak_drift(self):
        edges=np.arange(101)/100
        phase=(edges[:-1]+edges[1:])/2
        pulse=np.exp(-((phase-.225)/.022)**2)+.6*np.exp(-((phase-.725)/.027)**2)
        profile=np.broadcast_to(pulse,(12,100)).copy()
        shifted=np.roll(profile,2,axis=1)
        a={'phase_edges':edges,'profiles_width_2':profile}
        b={'phase_edges':edges,'profiles_width_2':shifted}
        # The observations may prefer cancelling phase shifts, but the
        # convergence comparison must still expose the physical +.02 drift.
        def row(name,az,fitphase):
            return {'name':name,'azimuths':az,'rings':6,'dt':.002,'field_sha256':'same',
                'width_scenarios':{'2.0':{'fitted_common_phase_cycles':fitphase,
                    'spectral_prediction_ph_cm2_s':np.ones(8),
                    'peaks_100_bins':peaks(profile,100,((.15,.4),(.6,.85)))}}}
        answer=compare_runs(row('a',128,0),a,row('b',256,-.02),b,[2.])
        self.assertEqual(answer['changed_parameters'],['azimuths'])
        width=answer['width_scenarios']['2.0']
        self.assertTrue(np.all(np.asarray(width['pulse_unit_mean_L1_100_bins'])>.3))
        assert_allclose(width['unshifted_peak_bin_centre_change_100_bins'],.02,atol=2e-16)
        assert_allclose(width['phase_mean_relative_change_all_bands'],0,atol=2e-15)
        assert_allclose(width['unshifted_peak_separation_change_100_bins'],0,atol=2e-16)


if __name__=='__main__':unittest.main()
