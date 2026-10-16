import unittest
from dataclasses import replace
import numpy as np
from pulsar.model import Parameters, fields, analytic_separation, pulse_profile
from pulsar.validation import field_validation, timing_validation


class PhysicalBenchmarks(unittest.TestCase):
    def test_maxwell_force_free_and_sheet_jump_conditions(self):
        result=field_validation()
        self.assertTrue(result['passed'], str(result))

    def test_arrival_time_caustics_and_repeatability(self):
        result=timing_validation()
        self.assertTrue(result['passed'], str(result))

    def test_no_sheet_crossing_has_no_geometric_double_pulse(self):
        self.assertIsNone(analytic_separation(np.deg2rad(20),np.deg2rad(30)))
        self.assertAlmostEqual(analytic_separation(np.pi/2,np.pi/2),0.5)

    def test_invalid_physical_domain_is_rejected(self):
        for kwargs in ({'r_min':0.2},{'r_max':1.},{'lorentz_factor':2.},
                       {'sheet_width':0.},{'observer_deg':181},{'period_s':float('nan')}):
            with self.assertRaises(ValueError):
                Parameters(**kwargs)

    def test_origin_is_excluded(self):
        with self.assertRaises(ValueError):
            fields(np.zeros((1,3)),0.,1.)

    def test_polar_observer_has_no_rotational_modulation(self):
        p=replace(Parameters(),observer_deg=0.)
        # Same angular sampling at every phase avoids azimuthal grid artifacts.
        f=pulse_profile(p,np.arange(32)/32,(6,48,128))
        self.assertLess(np.ptp(f)/np.mean(f),0.003)


if __name__=='__main__':
    unittest.main()
