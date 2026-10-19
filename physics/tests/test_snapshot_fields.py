"""Independent interpolation, normalization and pattern-motion checks."""
import unittest
import numpy as np
from numpy.testing import assert_allclose
from pulsar.snapshot_fields import RotatingSnapshot,SpinDownScale,rotate_z,sample_quality


class SnapshotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.axis=np.linspace(-3,3,65)
        x,y,z=np.meshgrid(cls.axis,cls.axis,cls.axis,indexing='ij')
        cls.e=np.array([x,-y,np.zeros_like(z)])
        cls.b=np.array([-y,x,1+0*z])

    def test_linear_affine_field_and_units(self):
        f=RotatingSnapshot((self.axis,)*3,self.e,self.b,amplitude_gauss=123.,reference_time=0.,pattern_omega=0,interpolation='linear')
        q=np.array([[.137,-.549,.345],[-.237,.469,-.577]])
        e,b=f(q,.48)
        assert_allclose(e,123*np.column_stack([q[:,0],-q[:,1],np.zeros(2)]),atol=1e-13)
        assert_allclose(b,123*np.column_stack([-q[:,1],q[:,0],np.ones(2)]),atol=1e-13)
        assert_allclose(sample_quality(f,q,0)['divB_per_RLC'],0,atol=1e-9)

    def test_cubic_smooth_field_error_decreases_under_refinement(self):
        q=np.array([[.17,.31,.23],[-.49,.18,-.4]])
        errors=[]
        for n in (25,49):
            a=np.linspace(-3,3,n);x,y,z=np.meshgrid(a,a,a,indexing='ij')
            b=np.array([np.sin(y),np.cos(z),np.sin(x)]);e=np.zeros_like(b)
            f=RotatingSnapshot((a,)*3,e,b,amplitude_gauss=1,reference_time=0,pattern_omega=0)
            exact=np.column_stack([np.sin(q[:,1]),np.cos(q[:,2]),np.sin(q[:,0])])
            errors.append(np.max(np.abs(f(q,0)[1]-exact)))
        self.assertGreater(errors[0]/errors[1],8)
        self.assertLess(errors[1],2e-6)

    def test_rotating_pattern_and_full_period(self):
        f=RotatingSnapshot((self.axis,)*3,self.e,self.b,amplitude_gauss=2,reference_time=.6,interpolation='linear')
        q=np.array([[.137,-.549,.345]])
        e,b=f(q,.6);angle=.79
        er,br=f(rotate_z(q,angle),.6+angle)
        assert_allclose(er,rotate_z(e,angle),atol=1e-13)
        assert_allclose(br,rotate_z(b,angle),atol=1e-13)
        assert_allclose(f(q,.6+2*np.pi),f(q,.6),atol=1e-13)

    def test_domain_and_release_checks(self):
        f=RotatingSnapshot((self.axis,)*3,self.e,self.b,amplitude_gauss=1,reference_time=0)
        with self.assertRaises(ValueError):f([3.1,0,0],0)
        with self.assertRaises(RuntimeError):f.require_emission_release()
        with self.assertRaises(ValueError):RotatingSnapshot((self.axis,)*3,self.e,self.b,amplitude_gauss=0,reference_time=0)

    def test_spindown_scaling_recovers_independent_torque(self):
        s=SpinDownScale(.089,1.25e-13,1e45,75.)
        mu=s.dipole_scale_gauss*s.rlc_cm**3;omega=2*np.pi/s.period_s
        luminosity=mu**2*omega**4/29979245800.**3*(1+np.sin(np.deg2rad(75))**2)
        expected=1e45*omega*(2*np.pi*1.25e-13/.089**2)
        self.assertAlmostEqual(luminosity/expected,1.,places=14)
        double_I=SpinDownScale(.089,1.25e-13,2e45,75.)
        self.assertAlmostEqual(double_I.dipole_scale_gauss/s.dipole_scale_gauss,np.sqrt(2),places=14)


if __name__=='__main__':unittest.main()
