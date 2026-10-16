import sys
from pathlib import Path
import unittest
import numpy as np
from numpy.testing import assert_allclose
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit_curl_source_independent import analytic_finite_dipole,surface_diagnostics,compatible_cube_flux,calibration,fit_source
from run_coffee_assessment import curls,rotation
from run_coffee_curl_audit import source_basis


class SourceDiagnosticChecks(unittest.TestCase):
    def test_finite_source_is_regular_at_origin(self):
        moment=np.array([.3,-.4,np.sqrt(.75)])
        with np.errstate(all='raise'):
            field=analytic_finite_dipole(np.zeros((1,3)),moment)
        assert_allclose(field,2*(1+1.5+1.875+2.1875+2.4609375)/.5**3*moment[None,:])

    def test_analytic_dipole_surface_moment_signed_flux_and_normalization(self):
        for moment in ([0.,0.,1.],[.3,-.4,np.sqrt(.75)],[0.,-2.,0.]):
            m=np.array(moment)
            for radius in (.25,.5,1.,2.):
                provider=lambda p:(np.zeros_like(p),analytic_finite_dipole(p,m))
                q=surface_diagnostics(provider,radius)
                assert_allclose(q['radial_dipole_harmonic_moment'],m,atol=2e-14)
                self.assertLess(abs(q['net_over_unsigned_flux']),2e-15)
                self.assertEqual(q['luminosity_mu2_omega4_over_c3'],0)
                if radius>=.5:self.assertLess(abs(q['half_unsigned_magnetic_flux']/(2*np.pi*np.linalg.norm(m)/radius)-1),3e-4)

    def test_known_wind_power_and_net_flux(self):
        def wind(p):
            r=np.linalg.norm(p,axis=-1);n=p/r[...,None];sin=np.sqrt(n[...,0]**2+n[...,1]**2)
            az=np.stack((-n[...,1],n[...,0],np.zeros_like(r)),axis=-1)/sin[...,None]
            theta=np.cross(az,n);et=-sin/r
            return et[...,None]*theta,n/r[...,None]**2+et[...,None]*az
        q=surface_diagnostics(wind,1.7)
        self.assertAlmostEqual(q['signed_magnetic_flux'],4*np.pi,places=12)
        self.assertAlmostEqual(q['luminosity_mu2_omega4_over_c3'],2/3,places=13)

    def test_source_basis_converges_to_independent_analytic_derivative(self):
        rows=calibration()['analytic_vs_D4'];errors=[r['D4_field_relative_error_to_analytic'] for r in rows]
        for previous,next_row,e0,e1 in zip(rows[:-1],rows[1:],errors[:-1],errors[1:]):
            order=np.log(e0/e1)/np.log(previous['h']/next_row['h'])
            self.assertGreater(order,3.5)
        self.assertLess(rows[-1]['fit_analytic_field_moment_error'],1e-4)

    def test_source_fit_phase_amplitude_and_conductor_E(self):
        axis=np.linspace(-.6,.6,33);h=axis[1]-axis[0]
        xyz=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'));mask=np.linalg.norm(xyz,axis=0)<.35
        initial=np.array([.3,-.4,np.sqrt(.75)]);phase=1.2
        # Measured lag and amplitude are prescribed independently of expected phase.
        moment=.97*rotation(initial,phase-.03)
        b=np.zeros_like(xyz);b[:,mask]=(source_basis(xyz[:,mask].T,h,.5)@moment).reshape(-1,3).T
        beta=np.array([-xyz[1],xyz[0],np.zeros_like(xyz[0])])
        e=-np.cross(beta,b,axisa=0,axisb=0,axisc=0)
        d={'axes':(axis,)*3,'e':e,'b':b,'time':phase,'metadata':{'star_radius':.5,'magnetic_phase_rad':phase,
            'config':dict(zip(('p1','p2','p3'),initial))}}
        q=fit_source(d)
        assert_allclose(q['moment'],moment,atol=1e-15)
        self.assertAlmostEqual(q['magnitude'],.97,places=14)
        self.assertAlmostEqual(q['phase_error_degrees'],np.degrees(-.03),places=12)
        self.assertLess(q['relative_fit_residual'],1e-14)
        self.assertEqual(q['conductor_E_relative_L2'],0.)

    def test_native_cube_flux_equals_D4_volume_divergence(self):
        n=33;a=np.linspace(-2,2,n);h=a[1]-a[0];x,y,z=np.meshgrid(a,a,a,indexing='ij')
        # Nonzero known divergence catches a spurious always-zero flux check.
        b=np.array([x,2*y,3*z]);q=compatible_cube_flux(b,(a,)*3,.7)
        self.assertAlmostEqual(q['signed_discrete_flux'],q['volume_integral_D4_divB'],places=12)
        self.assertGreater(q['signed_discrete_flux'],10)
        rng=np.random.default_rng(77);solenoidal=curls(rng.normal(size=b.shape),h)
        q=compatible_cube_flux(solenoidal,(a,)*3,.7)
        self.assertLess(abs(q['signed_discrete_flux']),1e-12)
        self.assertLess(abs(q['closure_error']),1e-12)


if __name__=='__main__':unittest.main()
