"""Geometry and presentation checks; no atmosphere or plasma solve is run."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import render_pulsar_3d as render
from surface_camera import Star


@pytest.fixture
def scene(monkeypatch):
    def synthetic_surface(self, cache):
        # The real atmosphere has separate validation. This keeps the geometry
        # tests independent of its table, transfer quadrature and cache files.
        yy, xx = np.mgrid[:self.height, :self.width]
        r = np.hypot(xx-self.center[0], yy-self.center[1])/self.scale
        disk = r < self.apparent_radius
        rgb = disk[..., None]*np.array([.6, .8, 1.], dtype='f4')
        return rgb, {'test_surface': True, 'mass_solar': 1.4, 'radius_m': 13000.}
    monkeypatch.setattr(render.PulsarScene, 'surface', synthetic_surface)
    return render.PulsarScene(width=180)


def test_closed_guides_are_tangent_to_the_adopted_dipole():
    """Independent B-vector test catches an incorrect obliquity/curve transform."""
    axis = render.magnetic_to_spin([0., 0., 1.])
    for line in render.dipole_loops(samples=600):
        r = np.linalg.norm(line, axis=-1)
        n = line/r[:, None]
        field = 3*n*(n@axis)[:, None]-axis
        tangent = line[2:]-line[:-2]
        b = field[1:-1]
        sine = np.linalg.norm(np.cross(tangent, b), axis=-1)
        sine /= np.linalg.norm(tangent, axis=-1)*np.linalg.norm(b, axis=-1)
        assert np.max(sine) < 5e-5
        np.testing.assert_allclose(r[[0, -1]], 1.3, atol=1e-12)


def test_nominal_polar_guides_match_last_open_dipole_surface():
    axis = render.magnetic_to_spin([0., 0., 1.])
    for line in render.polar_guides():
        r = np.linalg.norm(line, axis=-1)
        cos_theta = line@axis/r
        expected = (render.LIGHT_CYLINDER_M/render.RADIUS_M)*(1-cos_theta**2)
        np.testing.assert_allclose(r, expected, rtol=1e-12, atol=2e-12)
        assert r.max()*render.RADIUS_M/render.LIGHT_CYLINDER_M < .034


def test_rotation_preserves_obliquity_length_and_camera_handedness(scene):
    np.testing.assert_allclose(scene.basis@scene.basis.T, np.eye(3), atol=1e-14)
    assert np.linalg.det(scene.basis) == pytest.approx(1.)
    axis = scene.axis
    assert np.linalg.norm(axis) == pytest.approx(1.)
    assert np.rad2deg(np.arccos(axis[2])) == pytest.approx(75.)
    sample = scene.lines[3]
    for phase in (0, .17, .5, .81):
        np.testing.assert_allclose(np.linalg.norm(scene.camera(sample, phase), axis=-1),
                                   np.linalg.norm(sample, axis=-1), atol=3e-15)
        np.testing.assert_allclose(render.rotate_z(axis, 2*np.pi*phase)[2],
                                   np.cos(render.ALPHA), atol=1e-15)


def test_opaque_star_hides_back_and_interior_guides():
    # Positive camera z points toward the viewer. The silhouette includes its
    # exact limb; an exterior ray is not occulted even far behind the star.
    points = np.array([[0, 0, 2], [0, 0, -2], [0, 0, .5],
                       [.6, 0, .81], [.6, 0, .79], [1.01, 0, -2], [1, 0, 0]])
    np.testing.assert_array_equal(render.visible(points, 1.),
                                  [True, False, False, True, False, True, True])


def test_camera_and_composited_movie_are_exactly_periodic(scene):
    points = np.vstack((scene.lines[0], scene.axis))
    np.testing.assert_array_equal(scene.camera(points, 0), scene.camera(points, 1))
    np.testing.assert_allclose(scene.camera(points, .231), scene.camera(points, 1.231), atol=1e-14)
    frame0 = np.asarray(scene.frame(0))
    np.testing.assert_array_equal(frame0, np.asarray(scene.frame(1)))
    assert frame0.shape == (320, 180, 3)
    assert frame0.dtype == np.uint8
    # The endpoint derivative must join smoothly too, not merely the position.
    eps = 1e-5
    left = (scene.camera(points, 0)-scene.camera(points, -eps))/eps
    right = (scene.camera(points, eps)-scene.camera(points, 0))/eps
    np.testing.assert_allclose(left, right, rtol=0, atol=.003)


def test_guides_and_beams_can_be_removed_without_replacing_surface(scene):
    phase = .31
    original_star = scene.star.copy()
    fields = scene.field_light(phase)
    beams = scene.beam_light(phase)
    assert np.isfinite(fields).all() and np.isfinite(beams).all()
    assert np.max(fields) > 0 and np.max(beams) > 0
    assert np.mean(fields[..., 2]) > np.mean(fields[..., 0])
    scene.guides = scene.beams = False
    assert not np.any(scene.field_light(phase))
    assert not np.any(scene.beam_light(phase))
    np.testing.assert_array_equal(scene.star, original_star)
    np.testing.assert_array_equal(scene.frame(.1), scene.frame(.6))


def test_scene_parameters_share_the_existing_physical_model(scene):
    path = ROOT/'research/optical-finish/particle-emission/physical_synchrotron.py'
    spec = importlib.util.spec_from_file_location('pulsar_3d_pair_reference', path)
    pair_module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = pair_module
    spec.loader.exec_module(pair_module)
    pair, star = pair_module.PairModel(), Star()
    meta = scene.metadata()
    assert render.PERIOD == pair.period_s
    assert render.RADIUS_M == pair.star_radius_m == star.radius_m
    assert meta['parameters']['mass_solar'] == star.mass_solar
    assert meta['parameters']['polar_field_gauss'] == pair.polar_field_t*1e4
    assert meta['parameters']['magnetic_obliquity_deg'] == pair.inclination_deg
    assert render.LIGHT_CYLINDER_M == pytest.approx(pair.rlc_m)
    assert scene.apparent_radius == pytest.approx(star.apparent_radius_m/star.radius_m)
    assert 0 < scene.beam_angle < np.deg2rad(15.)
    assert meta['beams']['reference_height_assumed']
    assert not meta['calculated_optical_particle_archive_displayed']
    assert not meta['physical_scale']['star_separately_enlarged']
