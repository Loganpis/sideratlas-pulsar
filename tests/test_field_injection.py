import unittest
from unittest.mock import patch
from dataclasses import replace
import numpy as np
from pulsar.field_injection import CapConfig,_classify,trace_cap_rim,snapshot_injection

def static_dipole(p,t):
    p=np.asarray(p);r=np.linalg.norm(p,axis=-1);n=p/r[...,None]
    moment=np.array([0.,0.,1.])
    b=(3*n*np.sum(n*moment,axis=-1)[...,None]-moment)/r[...,None]**3
    return np.zeros_like(b),b

class CapChecks(unittest.TestCase):
    def test_polar_split_monopole_can_escape_outward_without_cylinder_crossing(self):
        # A radial polar line has R_cyl=0 at every radius. It is a valid
        # outward finite-domain escape, not evidence of a cylinder crossing.
        def split_monopole(p,t):
            r=np.linalg.norm(p,axis=-1)
            b=np.sign(p[...,2])[...,None]*p/r[...,None]**3
            return np.zeros_like(b),b
        cfg=CapConfig(inclination_deg=0.,seed_radius=.2,inner_radius=.19,
                      step=.02,max_length=4.,rim_boundary='outgoing_sphere',azimuths=4)
        for pole in (1,-1):
            np.testing.assert_array_equal(_classify(split_monopole,np.array([0.]),np.array([0.]),pole,cfg),[2])

    def test_outgoing_boundary_rim_is_opt_in_and_keeps_selected_seed_boundary_gate(self):
        cfg=CapConfig(inclination_deg=0.,seed_radius=.2,inner_radius=.19,azimuths=4)
        def sphere_open(field,theta,az,pole,cfg):return np.where(theta<.6,2,-1)
        with patch('pulsar.field_injection._classify',side_effect=sphere_open):
            with self.assertRaisesRegex(RuntimeError,'light-cylinder'):
                trace_cap_rim(None,np.arange(4)*np.pi/2,1,cfg)
            finite=replace(cfg,rim_boundary='outgoing_sphere')
            rim=trace_cap_rim(None,np.arange(4)*np.pi/2,1,finite)
            np.testing.assert_allclose(rim,.6,atol=1e-5,rtol=0.)
        with patch('pulsar.field_injection.trace_cap_rim',return_value=np.full(4,.6)):
            with patch('pulsar.field_injection._classify',return_value=np.array([1,2,-1,1]*3)):
                with self.assertRaisesRegex(RuntimeError,'Selected annular seeds'):
                    snapshot_injection(None,finite)

    def test_return_to_launch_sphere_closes_line_above_inner_transport_cutoff(self):
        # B=(z,0,xc-x) is divergence-free. Its x/z field lines are exact
        # circles around (xc,0,0), with every point on this circle above .695.
        # The chosen intersection starts outward from r=.7, then returns to
        # r=.7; requiring r=.68 would leave it unclassified forever.
        radius,xc,orbit_radius=.7,.705,.01
        x=(radius**2+xc**2-orbit_radius**2)/(2*xc)
        z=np.sqrt(radius**2-x*x)
        theta=np.array([np.arctan2(x,z)])
        def circular_field(p,t):
            p=np.asarray(p)
            b=np.stack((p[...,2],np.zeros(p.shape[:-1]),xc-p[...,0]),axis=-1)
            return np.zeros_like(b),b
        self.assertGreater(np.dot(circular_field(np.array([x,0.,z]),0)[1],np.array([x,0.,z])),0.)
        for inner in (.6,.68,.69):
            self.assertGreater(xc-orbit_radius,inner)
            cfg=CapConfig(inclination_deg=0.,seed_radius=radius,inner_radius=inner,
                          step=.001,max_length=3.,azimuths=4)
            np.testing.assert_array_equal(_classify(circular_field,theta,np.array([0.]),1,cfg),[-1])

    def test_tangent_orbit_without_surface_reentry_remains_unclassified(self):
        # A circular line on the launch sphere has no inward crossing and
        # never reaches the light cylinder. A timeout must not become open.
        def tangent_field(p,t):
            p=np.asarray(p)
            b=np.stack((-p[...,1],p[...,0],np.zeros(p.shape[:-1])),axis=-1)
            return np.zeros_like(b),b
        cfg=CapConfig(inclination_deg=0.,seed_radius=.7,inner_radius=.68,
                      step=.01,max_length=3.,azimuths=4)
        with self.assertRaisesRegex(RuntimeError,'Unclassified'):
            _classify(tangent_field,np.array([np.pi/2]),np.array([0.]),1,cfg)

    def test_null_and_nonmonotonic_topology_are_still_rejected(self):
        cfg=CapConfig(inclination_deg=0.,seed_radius=.7,inner_radius=.68,azimuths=4)
        def null_field(p,t):return np.zeros_like(p),np.zeros_like(p)
        with self.assertRaisesRegex(ValueError,'Magnetic null'):
            _classify(null_field,np.array([.1]),np.array([0.]),1,cfg)
        def disconnected(field,theta,az,pole,cfg):
            return np.where((theta<.4)|((theta>.7)&(theta<1.)),1,-1)
        with patch('pulsar.field_injection._classify',side_effect=disconnected):
            with self.assertRaisesRegex(RuntimeError,'not star-shaped'):
                trace_cap_rim(null_field,np.arange(4)*np.pi/2,1,cfg)

    def test_full_spoke_bracket_finds_connected_cap_beyond_old_upper_angle(self):
        cfg=CapConfig(inclination_deg=0.,seed_radius=.7,inner_radius=.68,azimuths=4)
        def wide_cap(field,theta,az,pole,cfg):return np.where(theta<1.65,1,-1)
        with patch('pulsar.field_injection._classify',side_effect=wide_cap):
            rim=trace_cap_rim(None,np.arange(4)*np.pi/2,1,cfg)
        np.testing.assert_allclose(rim,1.65,atol=1e-5,rtol=0.)

    def test_dipole_exact_last_cylinder_crossing_both_poles(self):
        cfg=CapConfig(inclination_deg=0,seed_radius=.2,inner_radius=.19,step=.018,
                      azimuths=4,rings=1,bisections=13,max_length=10)
        az=np.arange(4)*np.pi/2
        expected=np.arcsin(np.sqrt(cfg.seed_radius))
        for pole in (1,-1):
            rim=trace_cap_rim(static_dipole,az,pole,cfg)
            np.testing.assert_allclose(rim,expected,atol=2e-4,rtol=0)

    def test_weights_and_seeds_retain_pole_symmetry_without_copying_rim(self):
        cfg=CapConfig(inclination_deg=0,seed_radius=.2,inner_radius=.19,step=.018,
                      azimuths=4,rings=2,bisections=10,max_length=10)
        bundle=snapshot_injection(static_dipole,cfg)
        self.assertEqual(bundle['position'].shape,(16,3))
        np.testing.assert_allclose(bundle['position'][:8],-bundle['position'][8:],atol=1e-14)
        np.testing.assert_allclose(bundle['seed_weight'][:8],bundle['seed_weight'][8:],atol=1e-14)
        self.assertAlmostEqual(bundle['seed_weight'].sum(),1)
        self.assertFalse(bundle['population_is_predicted'])
