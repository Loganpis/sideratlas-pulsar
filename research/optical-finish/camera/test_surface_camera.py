"""Independent geodesic/relativity limits of the calculated surface camera."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import roots_legendre

spec = importlib.util.spec_from_file_location('surface_camera', Path(__file__).with_name('surface_camera.py'))
camera = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = camera
spec.loader.exec_module(camera)


def test_flat_space_surface_geometry_and_travel_time_including_limb():
    star = camera.Star(0, 13000)
    q = np.array([0, .1, .5, .9, .999999, 1])
    psi, mu, delay = camera.outward_geodesics(q*star.radius_m, star)
    np.testing.assert_allclose(psi, np.arcsin(q), rtol=0, atol=2e-10)
    np.testing.assert_allclose(mu, np.sqrt(1-q*q), atol=1e-14)
    np.testing.assert_allclose(delay, star.radius_m/camera.C_M_S*(1-mu), atol=2e-14)
    rays = camera.surface_raymap(q*star.radius_m, np.zeros_like(q), star=star)
    np.testing.assert_allclose(rays.normal_static[:, 2], mu, atol=2e-10)
    np.testing.assert_allclose(rays.normal_static[:, 0], q, atol=2e-10)
    np.testing.assert_allclose(rays.photon_static, np.broadcast_to([0, 0, 1], (len(q), 3)), atol=2e-10)
    np.testing.assert_allclose(rays.arrival_offset_s, -star.radius_m*mu/camera.C_M_S, atol=2e-14)


def test_bending_and_delay_match_independent_adaptive_radial_integrals():
    star = camera.Star()
    for fraction in (.1, .7, .98):
        b = fraction*star.apparent_radius_m
        k, u = b/star.radius_m, star.compactness
        root = lambda x: np.sqrt(1-k*k*x*x*(1-u*x))
        psi_ref = quad(lambda x: k/root(x), 0, 1, epsabs=1e-12)[0]
        # Independent original difference formula; avoid x=0 endpoint explicitly.
        f = lambda x: (1/root(x)-1)/(x*x*(1-u*x)) if x > 1e-6 else k*k/2
        delay_ref = star.radius_m/camera.C_M_S*quad(f, 0, 1, epsabs=1e-11)[0]
        psi, _, delay = camera.outward_geodesics(b, star)
        np.testing.assert_allclose(psi, psi_ref, atol=3e-12)
        np.testing.assert_allclose(delay, delay_ref, atol=3e-13)


def test_compact_star_reveals_far_side_and_quad_converges():
    star = camera.Star()
    b = np.array([0, .5, .99, 1])*star.apparent_radius_m
    a = camera.outward_geodesics(b, star, order=96)
    z = camera.outward_geodesics(b, star, order=384)
    assert a[0][-1] > np.pi/2
    assert 1 < star.apparent_radius_m/star.radius_m < 1.3
    np.testing.assert_allclose(a[0], z[0], atol=2e-10)
    np.testing.assert_allclose(a[2], z[2], atol=2e-14)


def test_aberration_preserves_null_direction_and_emission_mu():
    beta = np.array([.1, -.05, 0.])
    k = np.array([.3, .4, np.sqrt(.75)])
    comoving, doppler = camera.aberrate_to_comoving(k, beta)
    np.testing.assert_allclose(np.dot(comoving, comoving), 1, atol=1e-14)
    np.testing.assert_allclose(comoving[2], doppler*k[2], atol=1e-14)
    _, forward = camera.aberrate_to_comoving([1, 0, 0], [.2, 0, 0])
    _, backward = camera.aberrate_to_comoving([-1, 0, 0], [.2, 0, 0])
    np.testing.assert_allclose(forward/backward, 1.5, atol=1e-14)


def test_uniform_surface_flux_has_correct_lensed_area_and_redshift():
    star = camera.Star()
    # Polar quadrature integrates physical image area independently of ray map.
    z, w = roots_legendre(16)
    radii2 = (z+1)/2*star.apparent_radius_m**2
    phi = (np.arange(32)+.5)*2*np.pi/32
    x = np.sqrt(radii2)[:, None]*np.cos(phi)
    y = np.sqrt(radii2)[:, None]*np.sin(phi)
    area = np.broadcast_to((w/2*np.pi*star.apparent_radius_m**2/32)[:, None], x.shape)
    rays = camera.surface_raymap(x, y, star=star)
    wavelength = np.linspace(400, 700, 11)
    # Constant synthetic spectral radiance fixture, NOT an atmosphere model.
    callback = lambda **kw: np.ones_like(kw['wavelength_emitted_nm'])*3.
    image = camera.spectral_surface_image(rays, 0, wavelength, callback)
    distance = 1e18
    flux = camera.integrate_surface_flux(image, area, distance)
    expected = 3*star.frequency_ratio_g**5*np.pi*star.apparent_radius_m**2/distance**2
    np.testing.assert_allclose(flux, expected, rtol=2e-14)


def test_body_phase_uses_retarded_time_and_spectral_coordinate_shift():
    star = camera.Star()
    rays = camera.surface_raymap(np.array([0, .8*star.apparent_radius_m]), [0, 0], star=star)
    seen = {}
    def callback(**kw):
        seen.update(kw)
        return np.ones_like(kw['wavelength_emitted_nm'])
    omega, t = 70., .01
    wave = np.array([400., 700.])
    camera.spectral_surface_image(rays, t, wave, callback, omega_rad_s=omega, spin_axis=[0, 1, 0])
    np.testing.assert_allclose(seen['emission_time_s'], t-rays.arrival_offset_s)
    expected = camera.rotate_vectors(rays.normal_static, [0, 1, 0], -omega*seen['emission_time_s'])
    np.testing.assert_allclose(seen['normal_body'], expected, atol=1e-14)
    assert np.all(seen['wavelength_emitted_nm'] < wave)


def test_missed_pixels_stay_black_and_do_not_request_surface_data():
    rays = camera.surface_raymap([1e6], [0])
    def forbidden(**kw):
        raise AssertionError('Missed rays have no atmospheric emission')
    result = camera.spectral_surface_image(rays, 0, [400, 700], forbidden)
    assert np.all(result == 0)


def test_invalid_model_inputs_are_rejected():
    with pytest.raises(ValueError, match='photon sphere'):
        camera.Star(1.4, 5000)
    with pytest.raises(ValueError):
        camera.outward_geodesics(1e6)
    with pytest.raises(ValueError, match='subluminal'):
        camera.aberrate_to_comoving([1, 0, 0], [1, 0, 0])
    with pytest.raises(ValueError, match='unit'):
        camera.aberrate_to_comoving([2, 0, 0], [0, 0, 0])
