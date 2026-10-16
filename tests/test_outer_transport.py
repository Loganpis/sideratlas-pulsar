"""Analytic and independently integrated checks of outer transport.

These validate a prescribed-velocity model, not a physical plasma solution.
"""
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad, solve_ivp

from pulsar.model import fields as monopole_fields
from pulsar.outer_transport import (
    OuterTraceConfig, _advance_signed, _cooling_integral, arrival_phase,
    inertial_curvature, prescribed_velocity, select_branch,
    signed_transport_curvature, sky_phase, trace_outer,
)
from pulsar.radiation import E_STATC, MEC2_ERG, transport_curvature
from pulsar.published_acceleration import NearForceFreeFIDO
from pulsar.trajectories import rotate_z


def uniform(e, b):
    def provider(position, time):
        return np.broadcast_to(e, np.shape(position)), np.broadcast_to(b, np.shape(position))
    return provider


class VelocityTests(unittest.TestCase):
    def test_uniform_crossed_parallel_and_electric_dominated_fields(self):
        # Crossed magnetically dominated case: drift x parallel y = z.
        e, b = np.array([.6, 0., 0.]), np.array([0., 1., 0.])
        assert_allclose(prescribed_velocity(e, b, 1), [0., .8, .6], atol=2e-15)
        assert_allclose(prescribed_velocity(e, b, -1), [0., -.8, .6], atol=2e-15)
        assert_allclose(prescribed_velocity([0., 4., 0.], b, -1), [0., -1., 0.])
        # E0²=E²-B² regularizes drift without rescaling physical E or B.
        assert_allclose(prescribed_velocity([2., 0., 0.], b), [0., np.sqrt(.75), .5])
        electric = np.array([[.2, .4, 1.1], [3., .4, -.2], [0., 0., 0.]])
        magnetic = np.array([[1., .1, .2], [.1, .2, .1], [1., 2., 3.]])
        velocity = prescribed_velocity(electric, magnetic, [1, -1, 1])
        assert_allclose(np.linalg.norm(velocity, axis=-1), 1., atol=4e-15)
        assert_allclose(prescribed_velocity(electric*1e120, magnetic*1e120, [1, -1, 1]), velocity, atol=3e-15)
        assert_allclose(prescribed_velocity(electric*1e-120, magnetic*1e-120, [1, -1, 1]), velocity, atol=3e-15)

    def test_invalid_fields_and_signs_raise(self):
        for e, b, sign in [([0, 0, 0], [0, 0, 0], 1),
                           ([1, 0, 0], [0, 0, 0], 1),
                           ([0, 0, 0], [1, 0, 0], 0),
                           ([np.nan, 0, 0], [1, 0, 0], 1)]:
            with self.assertRaises(ValueError):
                prescribed_velocity(e, b, sign)

    def test_select_outflow_inflow_and_unavailable_inflow(self):
        provider = uniform([0., 0., 0.], [1., 0., 0.])
        points = np.array([[2., 0., 0.], [-2., 0., 0.]])
        assert_allclose(select_branch(provider, points, flow=1), [1., -1.])
        assert_allclose(select_branch(provider, points, flow=-1), [-1., 1.])
        wind = lambda p, t: monopole_fields(p, t, .4, polarity=1)
        with self.assertRaisesRegex(ValueError, "unavailable"):
            select_branch(wind, np.array([3., 0., 0.]), flow=-1)

    def test_signed_charge_and_branch_do_not_force_positive_work(self):
        provider = uniform([2., 0., 0.], [3., 0., 0.])
        result = trace_outer(provider, np.array([[.5, .1, 0.]]*4),
                             branch=[1, 1, -1, -1], charge_sign=[1, -1, 1, -1],
                             config=OuterTraceConfig(dt=.01, max_time=.02, r_min=.1))
        assert_allclose(result["e_work"], [[2., -2., -2., 2.]]*2)
        assert_allclose(result["direction"][0, :, 0], [1., 1., -1., -1.])


class GeometryTests(unittest.TestCase):
    def test_sharp_fido_interface_rejects_artificial_curvature_and_rk_crossings(self):
        # Smooth, divergence-free B with nonzero parallel curl and magnetic
        # dominance. Only FIDO's screening switch makes the velocity jump.
        def smooth_base(p,t):
            p=np.asarray(p,float)
            b=np.stack((np.ones(p.shape[:-1]),20*(p[...,0]-1),
                        np.ones(p.shape[:-1])),axis=-1)
            e=np.zeros_like(b);e[...,1]=1.2
            return e,b
        f=NearForceFreeFIDO(smooth_base,derivative_step=.001)
        boundary=np.array([1.,0.,0.])
        for h in (1e-3,1e-4,1e-5,1e-6,1e-7):
            with self.assertRaisesRegex(ValueError,'sharp screening interface'):
                inertial_curvature(f,boundary,h=h)
        # A step must fail even when the curvature midpoint misses the switch.
        with self.assertRaisesRegex(ValueError,'sharp screening interface'):
            trace_outer(f,[[.999,0,0]],branch=1,
                        config=OuterTraceConfig(dt=.01,max_time=.02,r_min=.5,r_max=2))
        # Smooth-side derivative refinement remains usable on either side.
        for x in (.99,1.01):
            kappas=[1/inertial_curvature(f,[x,0,0],h=h)
                    for h in (1e-4,1e-5,1e-6)]
            assert_allclose(kappas[-1],kappas[-2],rtol=2e-8)
            self.assertLess(abs(kappas[-1]-kappas[-2]),
                            abs(kappas[-2]-kappas[-3])/20)

    def test_split_monopole_outward_is_radial_and_has_zero_curvature(self):
        points = np.array([[1.2, .3, .4], [2., -.4, -.2], [-2.3, .7, .2]])
        for polarity in (-1, 1):
            provider = lambda p, t: monopole_fields(p, t, .7, polarity=polarity)
            branch = select_branch(provider, points)
            direction = points/np.linalg.norm(points, axis=-1)[:, None]
            assert_allclose(prescribed_velocity(*provider(points, .2), branch), direction, atol=4e-15)
            self.assertTrue(np.all(np.isinf(inertial_curvature(provider, points, .2, branch))))
            trace = trace_outer(provider, points, branch, charge_sign=-1,
                                config=OuterTraceConfig(dt=.08, max_time=3., r_max=3.))
            self.assertTrue(trace["all_escaped"])
            self.assertTrue(np.all(trace["status"] == "r_max"))
            assert_allclose(trace["end_position"], 3*direction, atol=3e-11)
            assert_allclose(trace["dt"].sum(axis=0), 3-np.linalg.norm(points, axis=-1), atol=3e-11)
            active = trace["dt"] > 0
            self.assertTrue(np.all(np.isinf(trace["rho"][active])))
            assert_allclose(trace["e_work"][active], 0., atol=5e-16)
            # An outward relativistic radial wind has no curvature power,
            # despite its spiral B lines. This guards against line curvature.
            energy = signed_transport_curvature(trace["dt"]*1e8, trace["rho"]*1e8,
                                                0., np.full(len(points), 1e7),
                                                energy_bands_gev=[[.1, 10.]])
            assert_allclose(energy.radiated_erg, 0., atol=0.)
            assert_allclose(energy.photon_yield, 0., atol=0.)

    def test_curvature_of_time_dependent_field_is_not_frozen_field_curvature(self):
        def rotating(position, time):
            t = np.broadcast_to(time, np.shape(position)[:-1])
            b = np.stack((np.cos(t), np.sin(t), np.zeros_like(t)), axis=-1)
            return np.zeros_like(b), b
        points = np.array([[.6, .2, .1], [1.4, -.5, .1]])
        assert_allclose(inertial_curvature(rotating, points, [.3, .7]), 1., rtol=2e-9)
        result = trace_outer(rotating, points, config=OuterTraceConfig(dt=.04, max_time=.8, r_max=5.), start_time=.3)
        exact = points+np.array([np.sin(1.1)-np.sin(.3), np.cos(.3)-np.cos(1.1), 0.])
        assert_allclose(result["end_position"], exact, atol=1e-9)
        self.assertTrue(np.all(result["status"] == "horizon"))
        self.assertFalse(result["all_escaped"])
        assert_allclose(result["rho"], 1., rtol=2e-9)

    def test_curvature_against_independently_integrated_lab_trajectory(self):
        def provider(p, t):
            x = np.asarray(p)
            t = np.broadcast_to(t, x.shape[:-1])
            b = np.stack((-x[..., 1]+.3*np.cos(t), x[..., 0]+.2*np.sin(t),
                          np.full_like(t, .4)), axis=-1)
            e = np.stack((.05*np.sin(t), .03*np.cos(t), np.full_like(t, .02)), axis=-1)
            return e, b
        # Independent formula and ODE; no production velocity or curvature.
        def reference(t, p):
            e, b = provider(p, t)
            b2, e2, eb = np.dot(b, b), np.dot(e, e), np.dot(e, b)
            e02 = .5*(np.sqrt((b2-e2)**2+4*eb*eb)-(b2-e2))
            drift = np.cross(e, b)/(b2+e02)
            return drift+np.sqrt(1-np.dot(drift, drift))*b/np.sqrt(b2)
        p, t = np.array([1.3, -.4, .2]), .7
        accelerations = []
        for h in (.002, .001, .0005):
            plus = solve_ivp(reference, (t, t+h), p, rtol=1e-12, atol=1e-14).y[:, -1]
            minus = solve_ivp(reference, (t, t-h), p, rtol=1e-12, atol=1e-14).y[:, -1]
            a = (reference(t+h, plus)-reference(t-h, minus))/(2*h)
            accelerations.append(np.linalg.norm(np.cross(reference(t, p), a)))
        actual = 1/inertial_curvature(provider, p, t, h=1e-5)
        assert_allclose(actual, accelerations[-1], rtol=2e-7)
        self.assertLess(abs(accelerations[-1]-actual), abs(accelerations[0]-actual)/10)

    def test_inflow_exit_clips_time_and_padding(self):
        provider = uniform([0., 0., 0.], [1., 0., 0.])
        result = trace_outer(provider, [[.5, 0., 0.], [.8, 0., 0.]], -1,
                             config=OuterTraceConfig(dt=.07, max_time=2., r_min=.2))
        assert_allclose(result["end_position"][:, 0], .2, atol=2e-11)
        assert_allclose(result["dt"].sum(axis=0), [.3, .6], atol=2e-11)
        self.assertTrue(np.all(result["status"] == "r_min"))
        active = result["dt"] > 0
        self.assertTrue(np.all(np.diff(active.astype(int), axis=0) <= 0))
        assert_allclose(result["direction"][~active], 0.)

    def test_fixed_and_rotating_pattern_arrival_phases(self):
        p = np.array([[1.2, -.7, .1], [-.2, 1.3, .4]])
        t = np.array([.4, 2.1])
        v = np.array([[.8, .6, 0.], [.3, .4, np.sqrt(.75)]])
        azimuth = np.arctan2(v[:, 1], v[:, 0])
        observer = rotate_z(v, -azimuth)
        rotated = rotate_z(p, -azimuth)
        explicit = arrival_phase(rotated, t-azimuth, observer)
        assert_allclose(sky_phase(p, t, v), explicit, atol=3e-16)
        assert_allclose(arrival_phase(p, t, v), (t-np.sum(v*p, axis=-1))/(2*np.pi) % 1)
        # A simultaneous rigid rotation/time shift leaves the sky phase fixed.
        assert_allclose(sky_phase(rotate_z(p, .6), t+.6, rotate_z(v, .6)), explicit, atol=3e-16)
        with self.assertRaises(ValueError):
            sky_phase([1., 0., 0.], 0., [0., 0., 1.])


class SignedEnergyTests(unittest.TestCase):
    def test_signed_integral_against_independent_quadrature(self):
        for lo, hi in [(1e-8, 1e-4), (.01, .9), (.4, 3.), (2., 20.), (1e3, 2e4)]:
            reference = quad(lambda x: 1/(1+x**4), lo, hi, epsabs=1e-30, epsrel=2e-12)[0]
            assert_allclose(_cooling_integral(lo, hi), reference, rtol=2e-12)

    def test_signed_cooling_against_independent_adaptive_ode(self):
        for g, a, b, h in [(100., -2., .00001, .3), (1e7, -.03, 2e-29, 2e6),
                            (20., -2., 1e-10, 1.), (500., -1e-5, 1e-3, .002)]:
            reference = solve_ivp(lambda t, state: a-b*state**4, (0., h), [g],
                                  method="DOP853", rtol=3e-12, atol=1e-12).y[0, -1]
            assert_allclose(_advance_signed(np.array(g), a, b, h, 1.), reference, rtol=2e-11)

    def test_positive_reduces_to_existing_transport_and_energy_balances(self):
        ds = np.full((6, 2), 3e5)
        rho = np.array([[1e8, 4e8]]*6)
        work = np.array([[15., 40.]]*6)
        gamma0 = np.array([2e7, 3e7])
        old = transport_curvature(ds, rho, work, gamma0, [[.1, 1.], [1., 10.]], 3)
        new = signed_transport_curvature(ds, rho, work, gamma0, [[.1, 1.], [1., 10.]], 3)
        assert_allclose(new.gamma, old.gamma, rtol=0., atol=0.)
        assert_allclose(new.photon_yield, old.photon_yield, rtol=1e-13)
        signed = signed_transport_curvature(ds, rho, -work, gamma0, [[.1, 10.]], 3)
        self.assertTrue(np.all(signed.work_erg < 0))
        self.assertTrue(np.all(signed.radiated_erg > 0))
        self.assertTrue(np.all(np.diff(signed.gamma, axis=0) < 0))
        assert_allclose(signed.energy_balance_error_erg, 0., atol=5e-15)
        self.assertTrue(np.all(signed.photon_yield > 0))

    def test_straight_signed_work_and_inactive_padding(self):
        field = MEC2_ERG/E_STATC
        result = signed_transport_curvature([[2., 2.], [0., 0.]], np.inf,
                                            [[field, -field], [field, -field]], [100., 100.],
                                            [[.1, 100.]])
        assert_allclose(result.gamma, [[100., 100.], [102., 98.], [102., 98.]])
        assert_allclose(result.radiated_erg, 0., atol=0.)
        assert_allclose(result.photon_yield, 0., atol=0.)
        assert_allclose(result.energy_balance_error_erg, 0., atol=1e-19)

    def test_gamma_floor_is_not_silently_clamped(self):
        field = -MEC2_ERG/E_STATC
        with self.assertRaisesRegex(ValueError, "gamma_floor"):
            signed_transport_curvature([30.], np.inf, field, 20.)
        with self.assertRaisesRegex(ValueError, "gamma_floor"):
            _advance_signed(np.array(20.), -3., .001, 20., 10.)

    def test_negative_energy_near_floor_and_first_crossing(self):
        # Independent quadrature fixes the elapsed arclength needed to reach
        # a specified Lorentz factor. Exercise both tiny and large gamma/scale
        # branches of the analytic primitive, just above the validity floor.
        floor = 10.
        target = floor*(1+1e-8)
        for initial, acceleration, loss in [(1000., -1., 1e-28),
                                             (100., -1., 1e-4),
                                             (100., -1e-8, 1e-3)]:
            integrand = lambda g: 1/(-acceleration+loss*g**4)
            duration = quad(integrand, target, initial, epsabs=1e-22, epsrel=3e-13)[0]
            actual = _advance_signed(np.array(initial), acceleration, loss, duration, floor)
            assert_allclose(actual, target, rtol=2e-12)
            self.assertGreater(float(actual), floor)
            crossing = quad(integrand, floor, initial, epsabs=1e-22, epsrel=3e-13)[0]
            with self.assertRaisesRegex(ValueError, "gamma_floor"):
                _advance_signed(np.array(initial), acceleration, loss, crossing*(1+1e-7), floor)


if __name__ == "__main__":
    unittest.main()
