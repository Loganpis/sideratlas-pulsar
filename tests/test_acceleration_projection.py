"""Independent phase-sum check of the acceleration-scan sky integral."""
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from run_acceleration_scan import observer_cell_fraction
from pulsar.emission_quadrature import project_segment_emission


class DirectionalIntegralTests(unittest.TestCase):
    def test_observer_crossings_equal_sum_over_all_phase_bins(self):
        rng=np.random.default_rng(17);shape=(6,5)
        mu=rng.uniform(-1,1,shape+(2,));phase=rng.uniform(-2,2,shape+(2,))
        # Exact reconstruction provided so no trajectory approximation enters.
        phase[...,1]=phase[...,0]+rng.uniform(-.2,.2,shape)
        y=rng.uniform(0,5,shape+(3,));weight=np.arange(1,6,dtype=float)/15
        bundle={'dt':np.ones(shape),'seed_weight':weight}
        for observer,width in [(65.,2.),(90.,100.),(20.,30.)]:
            f,solid=observer_cell_fraction(mu,observer,width)
            expected=np.einsum('snb,sn,n->b',y,f,weight)/solid
            for bins in [1,7,100]:
                p=project_segment_emission(bundle,y,observer_deg=observer,
                    view_width_deg=width,nphase=bins,sky_endpoints=(mu,phase))
                np.testing.assert_allclose(p['profiles'].mean(axis=1),expected,rtol=2e-13,atol=1e-13)

    def test_stationary_ray_and_reversal(self):
        m=np.cos(np.deg2rad(65.))
        mu=np.array([[m,m],[-1.,1.],[1.,-1.]])
        f,solid=observer_cell_fraction(mu,65.,2.)
        self.assertEqual(f[0],1)
        self.assertAlmostEqual(f[1],f[2])
        self.assertAlmostEqual(f[1],solid/(4*np.pi))
        for angle,width in [(0.,2.),(180.,2.),(65.,0.),(np.nan,2.)]:
            with self.assertRaises(ValueError):observer_cell_fraction(mu,angle,width)
