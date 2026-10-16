import unittest
import numpy as np
from pulsar.emission_projection import project_emission

class ProjectionChecks(unittest.TestCase):
    def test_conservation_and_relative_spectrum(self):
        mu=np.linspace(-.995,.995,200);phase=(np.arange(200)*.137)%1
        direction=np.stack([np.sqrt(1-mu*mu),np.zeros_like(mu),mu],axis=-1)[None]
        b={'direction':direction,'phase':phase[None],'dt':np.ones((1,200)),
           'seed_weight':np.ones(200)/200}
        photons=np.ones((1,200,2))*np.array([2.,7.])
        total=np.zeros(2)
        for low,high in zip(np.arange(0,180,10),np.arange(10,181,10)):
            p=project_emission(b,photons,observer_deg=(low+high)/2,view_width_deg=10,nphase=128)
            total+=p['phase_averaged_yield_per_sr']*p['solid_angle_sr']
            np.testing.assert_allclose(p['profiles'][1],3.5*p['profiles'][0],atol=1e-12)
        np.testing.assert_allclose(total,[2,7],rtol=1e-14)

    def test_phase_wrap_and_no_padding_emission(self):
        b={'direction':np.array([[[1.,0,0]],[[0.,0,0]]]),'phase':np.array([[.999],[np.nan]]),
           'dt':np.array([[1.],[0.]]),'seed_weight':np.ones(1)}
        p=project_emission(b,np.ones((2,1,1)),observer_deg=90,view_width_deg=2,nphase=128)
        self.assertEqual(p['selected_segments'],1)
        self.assertGreater(p['profiles'][0,0],0)
        self.assertEqual(np.count_nonzero(p['profiles']),2)
        self.assertAlmostEqual(p['phase_averaged_yield_per_sr'][0]*p['solid_angle_sr'],1)
