"""Independent identities and analytic GLM modes, not emission release tests."""
from pathlib import Path
import importlib.util
import unittest
import numpy as np

path=Path(__file__).resolve().parents[1]/'run_coffee_divergence_audit.py'
spec=importlib.util.spec_from_file_location('coffee_divergence_audit',path)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


class DivergenceAuditTests(unittest.TestCase):
    def test_compatible_curl_has_zero_discrete_divergence(self):
        field=np.random.default_rng(11).normal(size=(3,17,19,21))
        residual=audit.divergence(audit.curl(field,.1),.1)
        self.assertLess(np.max(np.abs(residual)),1e-12)

    def test_regular_potential_includes_origin_and_has_solenoidal_curl(self):
        axis=np.linspace(-1.5,1.5,33);h=axis[1]-axis[0]
        xyz=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'))
        a=audit.regular_dipole_potential(xyz,np.array([.8,0.,.6]))
        self.assertTrue(np.isfinite(a).all())
        b=audit.curl(a,h)
        self.assertTrue(np.isfinite(b).all())
        self.assertLess(np.max(np.abs(audit.divergence(b,h))),5e-12)

    def test_field_blending_injects_divergence_but_curl_blending_does_not(self):
        # Both initial states are exact discrete curls. Multiplying fields
        # by a spatial weight still violates the constraint.
        n=32;h=2*np.pi/n
        x,y,z=np.meshgrid(*(np.arange(n)*h for _ in range(3)),indexing='ij')
        zero=np.zeros_like(x)
        a0=np.stack((zero,zero,np.sin(x)*np.sin(y)))
        a1=np.stack((zero,zero,np.cos(x)*np.sin(y)))
        weight=.5+.2*np.cos(x)
        blended_b=weight*audit.curl(a1,h)+(1-weight)*audit.curl(a0,h)
        curl_blended_a=audit.curl(weight*a1+(1-weight)*a0,h)
        self.assertGreater(np.max(np.abs(audit.divergence(blended_b,h))),.1)
        self.assertLess(np.max(np.abs(audit.divergence(curl_blended_a,h))),1e-12)

    def test_grid_glm_matches_exact_semidiscrete_telegraph_mode(self):
        result=audit.glm_mode(n=48,dt=.02,duration=2.,tau_steps=2.)
        self.assertLess(result['max_abs_B_error'],1e-10)
        self.assertLess(result['max_abs_P_error'],1e-10)

    def test_fixed_physical_damping_has_fourth_order_time_convergence(self):
        coarse=audit.glm_mode(n=48,dt=.08,duration=2.,physical_tau=1.)
        fine=audit.glm_mode(n=48,dt=.04,duration=2.,physical_tau=1.)
        self.assertGreater(coarse['max_abs_B_error']/fine['max_abs_B_error'],12.)

    def test_fixed_timestep_tau_overdamps_long_waves_more_under_refinement(self):
        coarse=audit.glm_mode(n=48,dt=.04,duration=2.,tau_steps=2.)
        fine=audit.glm_mode(n=48,dt=.02,duration=2.,tau_steps=2.)
        self.assertAlmostEqual(coarse['tau_physical'],.08)
        self.assertAlmostEqual(fine['tau_physical'],.04)
        self.assertGreater(fine['divergence_amplitude_ratio'],coarse['divergence_amplitude_ratio'])


if __name__=='__main__':unittest.main()
