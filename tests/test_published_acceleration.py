import unittest
import numpy as np
from numpy.testing import assert_allclose
from pulsar.published_acceleration import NearForceFreeFIDO,Barnard2022Work,field_curl
from pulsar.radiation import MEC2_ERG,E_STATC


def helical_field(p,t):
    p=np.asarray(p,float);b=np.stack((-p[...,1],p[...,0],np.ones(p.shape[:-1])),axis=-1)
    return .3*b,b


class AccelerationChecks(unittest.TestCase):
    def test_analytic_curl_conductivity_units_and_screening(self):
        p=np.array([[.4,.2,.1],[1.4,.2,.1]])
        assert_allclose(field_curl(helical_field,p,0,1e-3),[[0,0,2],[0,0,2]],atol=2e-12)
        f=NearForceFreeFIDO(helical_field,sigma_over_omega=30,derivative_step=1e-3)
        e,b=f(p,0)
        assert_allclose(e[0],0,atol=1e-15)
        expected=2/(4*np.pi*30)*b[1]/np.dot(b[1],b[1])
        assert_allclose(e[1],expected,atol=1e-15)
        high=NearForceFreeFIDO(helical_field,sigma_over_omega=60,derivative_step=1e-3)
        assert_allclose(high(p,0)[0],e/2,atol=1e-15)

    def test_signed_parallel_field_not_replaced_by_absolute_value(self):
        def opposite_twist(p,t):
            p=np.asarray(p);b=np.stack((p[...,1],-p[...,0],np.ones(p.shape[:-1])),axis=-1)
            return np.zeros_like(b),b
        f=NearForceFreeFIDO(opposite_twist,derivative_step=1e-3)
        e,b=f([[1.5,.2,.1]],0)
        self.assertLess(float(np.sum(e*b)),0)

    def test_barnard_eq2_rate_is_per_path_length_without_extra_projection(self):
        a=Barnard2022Work();p=np.array([[.5,0,0],[1.5,0,0]])
        b=np.array([[0,0,3],[0,0,3]])
        v=np.array([[.6,0,.8],[.6,0,-.8]])
        expected=np.array([.04,.25])*MEC2_ERG/E_STATC
        assert_allclose(a.work_field(p,v,b),expected)
        # These effective Eq.(2) rates are not a global field-parallel force.
        # Even perpendicular drift does not insert an extra projection factor.
        assert_allclose(a.work_field(p,[[1,0,0],[1,0,0]],b),expected)
        from pulsar.outer_transport import signed_transport_curvature
        energy=signed_transport_curvature([[1e6,1e6]],np.inf,
                                           a.work_field(p,v,b)[None],1000.)
        assert_allclose(energy.gamma[-1],[41000.,251000.],rtol=1e-14)
        assert_allclose(energy.radiated_erg,0.,atol=0.)

    def test_fido_guard_checks_crossing_chords_and_each_stream(self):
        f=NearForceFreeFIDO(helical_field,derivative_step=1e-3)
        # Both endpoints are outside, but their straight stencil crosses LC.
        with self.assertRaisesRegex(ValueError,'sharp screening interface'):
            f.validate_transport_stencil([[-1.2,0,0],[1.2,0,0]])
        with self.assertRaisesRegex(ValueError,'sharp screening interface'):
            f.validate_transport_stencil([[[.8,0,0],[1.2,0,0]],
                                           [[.9,0,0],[.99,0,0]]])
        f.validate_transport_stencil([[[.8,0,0],[1.2,0,0]],
                                      [[.9,0,0],[1.3,0,0]]])


if __name__=='__main__':unittest.main()
