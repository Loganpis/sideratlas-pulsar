import unittest
import numpy as np
from pulsar.energetics import PC_CM,C_CM_S,population_power_budget,field_scale_from_source_power
from pulsar.snapshot_fields import SpinDownScale


class PopulationPowerTests(unittest.TestCase):
    def test_numerical_torque_scale_and_point_dipole_limit(self):
        p=.089;scale=SpinDownScale(p,1.25e-13,1e45,75.)
        point_power=1+np.sin(np.deg2rad(75.))**2
        measured=field_scale_from_source_power(p,scale.luminosity_erg_s,point_power)
        self.assertAlmostEqual(measured/scale.dipole_scale_gauss,1)
        half=field_scale_from_source_power(p,scale.luminosity_erg_s,point_power/2)
        self.assertAlmostEqual(half/measured,np.sqrt(2))
        self.assertAlmostEqual(C_CM_S*scale.rlc_cm**2*half**2*point_power/2/scale.luminosity_erg_s,1)

    def test_directional_normalization_and_conservation(self):
        # A synthetic 10-primary/s source with 2 erg emitted by each primary
        # has 20 erg/s luminosity, independent of its angular distribution.
        norm = 10/(2*PC_CM)**2
        result = population_power_budget(norm, 2, 100, 2,
            work_erg_per_primary=3, injected_kinetic_erg_per_primary=.5,
            escaping_kinetic_erg_per_primary=1.5)
        self.assertAlmostEqual(result['primaries_per_second'], 10)
        self.assertAlmostEqual(result['radiation_power_erg_s'], 20)
        self.assertAlmostEqual(result['modeled_supply_fraction_of_spin_down'], .35)
        self.assertAlmostEqual(result['energy_ledger_residual_erg_per_primary'], 0)
        self.assertFalse(result['vela_validation'])

    def test_distance_squared_and_impossible_budget(self):
        a = population_power_budget(1/PC_CM**2, 1, 3, 2)
        b = population_power_budget(1/PC_CM**2, 2, 3, 2)
        self.assertAlmostEqual(b['radiation_power_erg_s']/a['radiation_power_erg_s'], 4)
        self.assertFalse(a['radiation_exceeds_spin_down'])
        self.assertTrue(b['radiation_exceeds_spin_down'])

    def test_signed_work_and_input_validation(self):
        r = population_power_budget(1/PC_CM**2, 1, 10, .5,
            work_erg_per_primary=-.5, injected_kinetic_erg_per_primary=2,
            escaping_kinetic_erg_per_primary=1)
        self.assertAlmostEqual(r['energy_ledger_residual_erg_per_primary'], 0)
        for args in [(-1, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, -1)]:
            with self.assertRaises(ValueError): population_power_budget(*args)
