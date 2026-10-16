"""Independent limits and quadrature checks for primary seed loading."""
import unittest
from dataclasses import asdict, replace
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from pulsar.injection import _corotation_charge, seed_injection_weights
from pulsar.trajectories import Geometry, dipole_field, seed_positions, velocity


def bundle_quadrature(g, nphi=16, rings=3):
    az = (np.arange(nphi)+.5)*2*np.pi/nphi
    rim = np.full(nphi,np.arcsin(np.sqrt(g.seed_radius)))
    # This helper needs no expensive field-line tracing: injection only uses
    # the supplied rim quadrature, whose aligned limit is independently known.
    u = g.annulus_inner+(np.arange(rings)+.5)/rings*(g.annulus_outer-g.annulus_inner)
    theta = (u[:,None]*rim).ravel()
    phi = np.tile(az,rings)
    q = seed_positions(theta,phi,g.chi,g.seed_radius)
    field = dipole_field(q,g.chi)
    darea = g.seed_radius**2*np.sin(theta)*np.tile(rim,rings)*(g.annulus_outer-g.annulus_inner)/rings*2*np.pi/nphi
    weight = np.abs(np.sum(field*q/g.seed_radius,axis=-1))*darea
    weights = np.r_[weight,weight]
    weights /= weights.sum()
    bundle = dict(seed_weight=weights,cap_azimuth=az,cap_rim=rim,
                  n_rings=rings,n_azimuth=nphi,geometry=asdict(g),
                  pole=np.r_[np.ones(len(q)),-np.ones(len(q))])
    return bundle, theta, phi


class InjectionTests(unittest.TestCase):
    def test_uniform_flux_preserves_existing_weights_without_aliasing(self):
        g = Geometry()
        bundle,_,_ = bundle_quadrature(g)
        result = seed_injection_weights(bundle,g,"uniform_flux")
        assert_array_equal(result,bundle['seed_weight'])
        result[0] = 0.
        self.assertGreater(bundle['seed_weight'][0],0.)

    def test_aligned_static_dipole_matches_analytic_surface_current(self):
        g = Geometry(inclination_deg=0.)
        bundle,theta,_ = bundle_quadrature(g)
        result = seed_injection_weights(bundle,g,"goldreich_julian")
        # Static dipole Br=2cos(theta)/r^3, Bz=(3cos²(theta)-1)/r³.
        # Corotation is azimuthal; the unit-speed parallel component has
        # magnitude sqrt(1-r²sin²(theta)) and radial projection Br/|B|.
        co, si, r = np.cos(theta),np.sin(theta),g.seed_radius
        radial_speed = np.sqrt(1-r*r*si*si)*2*co/np.sqrt(1+3*co*co)
        expected_north = np.abs(3*co*co-1)*radial_speed*si
        expected = np.r_[expected_north,expected_north]
        expected /= expected.sum()
        assert_allclose(result,expected,rtol=2e-15,atol=1e-17)

    def test_oblique_pole_parity_normalization_and_nonuniform_loading(self):
        g = Geometry(inclination_deg=75.)
        bundle,_,_ = bundle_quadrature(g)
        result = seed_injection_weights(bundle,g,"goldreich_julian")
        half = len(result)//2
        assert_allclose(result[:half],result[half:],rtol=3e-15)
        self.assertAlmostEqual(result.sum(),1.,places=15)
        self.assertTrue(np.all(result >= 0))
        self.assertGreater(np.ptp(result[:16]),.003)
        self.assertGreater(np.linalg.norm(result-bundle['seed_weight']),.02)

    def test_reweighting_matches_flux_multiplier_at_true_seeds(self):
        g = Geometry()
        bundle,theta,phi = bundle_quadrature(g)
        q = seed_positions(theta,phi,g.chi,g.seed_radius)
        b = dipole_field(q,g.chi)
        v = velocity(q,g.chi)[0]
        radial = q/g.seed_radius
        multiplier = np.abs(b[:,2])*np.abs(np.sum(v*radial,axis=-1))/np.abs(np.sum(b*radial,axis=-1))
        expected = bundle['seed_weight']*np.r_[multiplier,multiplier]
        expected /= expected.sum()
        actual = seed_injection_weights(bundle,g,"goldreich_julian")
        assert_allclose(actual,expected,rtol=2e-15,atol=1e-17)
        # Trajectory midpoints are intentionally irrelevant to seed loading.
        bundle['position'] = np.full((1,len(actual),3),np.nan)
        assert_array_equal(seed_injection_weights(bundle,g,"goldreich_julian"),actual)
        bundle['seed_weight'] = actual
        assert_array_equal(seed_injection_weights(bundle,g,"goldreich_julian"),actual)

    def test_narrow_aligned_cap_approaches_uniform_flux_loading(self):
        g = Geometry(inclination_deg=0.,seed_radius=1e-5)
        bundle,_,_ = bundle_quadrature(g)
        actual = seed_injection_weights(bundle,g,"goldreich_julian")
        assert_allclose(actual,bundle['seed_weight'],rtol=4e-7)

    def test_corotation_charge_against_independent_electric_divergence(self):
        g = Geometry()
        bundle,theta,phi = bundle_quadrature(g)
        north = seed_positions(theta,phi,g.chi,g.seed_radius)
        q = np.concatenate((north,-north))
        def electric(x):
            beta = np.column_stack((-x[:,1],x[:,0],np.zeros(len(x))))
            return -np.cross(beta,dipole_field(x,g.chi))
        density = _corotation_charge(q,g.chi,1e-5)
        # Independent central derivative of the complete E field, rather
        # than repeating the code's vector identity and fourth-order curl.
        errors = []
        references = []
        for step in (4e-5,2e-5,1e-5):
            divergence = np.zeros(len(q))
            for axis in range(3):
                delta = np.zeros(3)
                delta[axis] = step
                divergence += (electric(q+delta)[:,axis]-electric(q-delta)[:,axis])/(2*step)
            references.append(divergence)
            errors.append(np.linalg.norm(divergence-density)/np.linalg.norm(density))
        self.assertLess(errors[-1],1e-7)
        self.assertGreater(errors[0]/errors[1],3.8)
        self.assertGreater(errors[1]/errors[2],3.8)
        v = velocity(q,g.chi,bundle['pole'])[0]
        radial = q/g.seed_radius
        field = dipole_field(q,g.chi)
        multiplier = np.abs(references[-1])*np.abs(np.sum(v*radial,axis=-1))/np.abs(np.sum(field*radial,axis=-1))
        expected = bundle['seed_weight']*multiplier
        expected /= expected.sum()
        assert_allclose(seed_injection_weights(bundle,g,'corotation_charge'),expected,rtol=3e-6,atol=1e-10)

    def test_aligned_corotation_charge_reduces_to_leading_gj_loading(self):
        g = Geometry(inclination_deg=0.)
        bundle,_,_ = bundle_quadrature(g)
        full = seed_injection_weights(bundle,g,'corotation_charge')
        leading = seed_injection_weights(bundle,g,'goldreich_julian')
        assert_allclose(full,leading,rtol=3e-12,atol=1e-15)

    def test_corotation_charge_pole_parity_normalization_and_step_refinement(self):
        g = Geometry()
        bundle,theta,phi = bundle_quadrature(g)
        result = seed_injection_weights(bundle,g,'corotation_charge')
        half = len(result)//2
        assert_allclose(result[:half],result[half:],rtol=5e-12,atol=1e-15)
        self.assertAlmostEqual(result.sum(),1.,places=15)
        self.assertTrue(np.all(result >= 0))
        q = seed_positions(theta,phi,g.chi,g.seed_radius)
        coarse = _corotation_charge(q,g.chi,1e-5)
        fine = _corotation_charge(q,g.chi,5e-6)
        assert_allclose(coarse,fine,rtol=1e-9,atol=1e-8)

    def test_tangent_field_and_zero_gj_density_are_rejected(self):
        g = Geometry()
        bundle,_,_ = bundle_quadrature(g)
        def tangent(q,chi):
            return np.column_stack((-q[:,1],q[:,0],np.zeros(len(q))))
        with patch('pulsar.injection.dipole_field',side_effect=tangent):
            with self.assertRaisesRegex(ValueError,'tangent'):
                seed_injection_weights(bundle,g,"goldreich_julian")
        def zero_bz(q,chi):
            return np.column_stack((q[:,0],q[:,1],np.zeros(len(q))))
        with patch('pulsar.injection.dipole_field',side_effect=zero_bz):
            with self.assertRaisesRegex(ValueError,'normalization'):
                seed_injection_weights(bundle,g,"goldreich_julian")

    def test_mismatched_geometry_and_malformed_quadrature_are_rejected(self):
        g = Geometry()
        bundle,_,_ = bundle_quadrature(g)
        with self.assertRaisesRegex(ValueError,'geometry differs'):
            seed_injection_weights(bundle,replace(g,seed_radius=.1),"goldreich_julian")
        with self.assertRaises(ValueError):
            seed_injection_weights(bundle,g,"fitted_per_phase")
        broken = dict(bundle,cap_azimuth=bundle['cap_azimuth'].copy())
        broken['cap_azimuth'][0] += .01
        with self.assertRaisesRegex(ValueError,'uniform'):
            seed_injection_weights(broken,g,"goldreich_julian")


if __name__ == '__main__':
    unittest.main()
