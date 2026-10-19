import unittest
import numpy as np
from numpy.testing import assert_allclose
from pulsar.emission_quadrature import project_segment_emission,reconstruct_sky_endpoints,rebin_phase_averages


def bundle(n=1):
    return {'dt':np.ones((n,1)),'seed_weight':np.ones(1)}


class SegmentProjectionChecks(unittest.TestCase):
    def test_bin_average_rebinning_and_periodic_shift(self):
        edges=np.array([0.,.25,.5,1.]);p=np.array([[4.,0.,2.]])
        assert_allclose(rebin_phase_averages(edges,p,[0,.5],[.5,1]),[[2.,2.]])
        assert_allclose(rebin_phase_averages(edges,p,[0,.25],[.25,.5],shift=.25),[[2.,4.]])
        assert_allclose(rebin_phase_averages(edges,p,[-.2],[1.8],shift=.173),[[2.]])

    def test_analytic_latitude_crossing_and_phase_wrap(self):
        # Linear mu=-1..1 traverses equal solid-angle bins uniformly. A whole
        # rotation in phase concurrently crosses wrap and all phase bins.
        b=bundle();y=np.array([[[2.,7.]]]);mu=np.array([[[-1.,1.]]]);phase=np.array([[[.9,1.9]]])
        total=np.zeros(2)
        for low,high in zip(range(0,180,10),range(10,181,10)):
            out=project_segment_emission(b,y,observer_deg=(low+high)/2,view_width_deg=10,
                                         nphase=100,sky_endpoints=(mu,phase))
            assert_allclose(out['phase_averaged_yield_per_sr'],np.array([2.,7.])/(4*np.pi),rtol=1e-13)
            total+=out['phase_averaged_yield_per_sr']*out['solid_angle_sr']
        assert_allclose(total,[2.,7.],rtol=1e-14)

    def test_uniform_phase_integrates_exactly_in_both_directions(self):
        for ends in ((.37,1.37),(1.37,.37)):
            out=project_segment_emission(bundle(),np.array([[[3.,9.]]]),observer_deg=90,
                view_width_deg=180,nphase=37,sky_endpoints=(np.zeros((1,1,2)),np.array([[ends]])))
            assert_allclose(out['profiles'],np.broadcast_to(np.array([[3.],[9.]])/(4*np.pi),(2,37)),rtol=3e-14)

    def test_cell_subdivision_preserves_exact_linear_map(self):
        n=31;mu=np.linspace(-.8,.7,n+1);phase=np.linspace(.95,1.31,n+1)
        exact=(np.array([[[-.8,.7]]]),np.array([[[.95,1.31]]]))
        coarse=project_segment_emission(bundle(),np.array([[[2.,5.]]]),observer_deg=70,view_width_deg=7,
                                       sky_endpoints=exact)
        fine=project_segment_emission(bundle(n),np.broadcast_to([2./n,5./n],(n,1,2)),observer_deg=70,view_width_deg=7,
             sky_endpoints=(np.stack((mu[:-1],mu[1:]),axis=-1)[:,None],np.stack((phase[:-1],phase[1:]),axis=-1)[:,None]))
        assert_allclose(coarse['profiles'],fine['profiles'],atol=3e-13)

    def test_reconstruction_of_linear_motion_and_partial_batch_weights(self):
        t=(np.arange(8)+.5)*.1;mu=-.2+.3*t;phi=(.95+.4*t)%1
        b={'dt':np.full((8,1),.1),'time':t[:,None],'phase':phi[:,None],
           'direction':np.stack((np.sqrt(1-mu*mu),np.zeros(8),mu),axis=-1)[:,None],
           'seed_weight':np.array([.25])}
        a,p,diag=reconstruct_sky_endpoints(b)
        assert_allclose(a[:,:,0],(-.2+.3*np.arange(8)*.1)[:,None],atol=2e-16)
        assert_allclose(p[:,:,1],(.95+.4*np.arange(1,9)*.1)[:,None],atol=3e-16)
        out=project_segment_emission(b,np.ones((8,1,1)),observer_deg=90,view_width_deg=180)
        assert_allclose(out['phase_averaged_yield_per_sr'],2/(4*np.pi),rtol=2e-15)

    def test_inactive_cells_do_not_emit(self):
        b=bundle(2);b['dt'][1]=0
        out=project_segment_emission(b,np.ones((2,1,1)),observer_deg=90,view_width_deg=180,
             sky_endpoints=(np.zeros((2,1,2)),np.zeros((2,1,2))))
        assert_allclose(out['phase_averaged_yield_per_sr'],1/(4*np.pi),rtol=2e-15)

    def test_undefined_polar_phase_is_rejected_only_if_visible(self):
        b=bundle();phase=np.full((1,1,2),np.nan)
        out=project_segment_emission(b,np.ones((1,1,1)),observer_deg=90,view_width_deg=2,
                                     sky_endpoints=(np.full((1,1,2),.99),phase))
        assert_allclose(out['profiles'],0.)
        with self.assertRaisesRegex(ValueError,'reaches observer'):
            project_segment_emission(b,np.ones((1,1,1)),observer_deg=90,view_width_deg=2,
                                     sky_endpoints=(np.zeros((1,1,2)),phase))

    def test_nonlinear_latitude_reconstruction_converges_to_analytic_profile(self):
        # mu(t)=2t^2-1, phase=.8+.6t, and unit photon rate. The observer's
        # entry/exit times have exact roots, independently defining bin mass.
        low=np.cos(np.deg2rad(65));high=np.cos(np.deg2rad(55))
        tlo=np.sqrt((low+1)/2);thi=np.sqrt((high+1)/2);solid=2*np.pi*(high-low)
        a=.8+.6*tlo;b=.8+.6*thi;nphase=128;edges=np.arange(nphase+1)/nphase
        mass=np.zeros(nphase)
        for period in (0,1):
            mass+=np.maximum(0,np.minimum(edges[1:]+period,b)-np.maximum(edges[:-1]+period,a))/.6
        exact=mass*nphase/solid;errors=[]
        for n in (20,40,80,160):
            t=(np.arange(n)+.5)/n;mu=2*t*t-1
            stream={'dt':np.full((n,1),1/n),'time':t[:,None],'phase':((.8+.6*t)%1)[:,None],
                    'direction':np.stack((np.sqrt(1-mu*mu),np.zeros(n),mu),axis=-1)[:,None],
                    'seed_weight':np.ones(1)}
            actual=project_segment_emission(stream,np.full((n,1,1),1/n),observer_deg=60,
                                            view_width_deg=10,nphase=nphase)['profiles'][0]
            errors.append(np.sum(abs(actual-exact))/np.sum(exact))
        # Boundary alignment makes the leading coefficient oscillate; the
        # envelope, rather than every adjacent ratio, must converge as h^2.
        for n,error in zip((20,40,80,160),errors):self.assertLess(error,8/n**2)
        self.assertLess(errors[-1],errors[0]/16)


if __name__=='__main__':unittest.main()
