import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import event_camera as camera


def test_periodic_arrival_intervals_conserve_every_event_across_wraps():
    left = np.array([-.1, .1, .9, 1.9, 80.1])
    right = left+np.array([.2, .7, .2, 1., .4])
    ids, bins, weights = camera.periodic_interval_weights(left, right, 1., 10)
    np.testing.assert_allclose(np.bincount(ids, weights=weights, minlength=5), 1, atol=2e-14)
    assert np.all((bins >= 0) & (bins < 10))
    np.testing.assert_allclose(np.bincount(bins[ids == 0], weights=weights[ids == 0], minlength=10)[[0, 9]], [.5, .5])


def test_static_source_steady_intervals_give_steady_flux():
    count = 64
    t = (np.arange(count)+.5)/count
    left, right = camera.rotating_pattern_arrival_edges(t, np.full(count, 1/count), np.zeros((count, 3)), [0, 0, 1], 0)
    ids, bins, fraction = camera.periodic_interval_weights(left, right, 1, 87)
    energy = camera.source_xyz_energy(np.ones((count, 3))*[1, 2, 3], np.ones(count), np.full(count, 1/count))
    received = np.zeros((87, 3)); np.add.at(received, bins, energy[ids]*fraction[:, None])
    np.testing.assert_allclose(received*87, np.broadcast_to([1, 2, 3], (87, 3)), atol=2e-13)


def test_rotating_edges_include_photon_arrival_compression():
    radius = .1*camera.C_M_S
    dt = np.array([1e-5])
    left, right = camera.rotating_pattern_arrival_edges([0.], dt, [[radius, 0, 0]], [0, 1, 0], 1.)
    np.testing.assert_allclose((right-left)/dt, .9, atol=1e-11)


def test_incoherent_particle_splitting_preserves_flux():
    p = np.array([[1., 3., 2.]])
    one = camera.source_xyz_energy(p, [100.], [.2])
    split = camera.source_xyz_energy(np.repeat(p, 10, axis=0), np.full(10, 10.), np.full(10, .2))
    np.testing.assert_allclose(one[0], split.sum(axis=0), atol=1e-14)


def test_rotating_arrival_inverse_matches_exact_retarded_geometry():
    rng = np.random.default_rng(31)
    center = rng.uniform(-.1, .2, 100)
    position = rng.normal(size=(100, 3))
    position *= (.2*camera.C_M_S/np.linalg.norm(position, axis=-1))[:, None]
    source = center+rng.uniform(-.01, .01, 100)
    actual = camera.rotate_vectors(position, [0, 0, 1], source-center)
    n = camera.unit_vector([1, 0, .4])
    arrival = source-np.sum(actual*n, axis=-1)/camera.C_M_S
    recovered = camera.inverse_rotating_arrival(arrival, center, position, n, 1.)
    np.testing.assert_allclose(recovered, source, atol=3e-16)


def test_pixel_projection_and_crop_keep_energy_ledger():
    pixel = np.array([[3.2, 4.5], [0, 0], [-.25, 3.], [8.7, 9.9], [100, 100]])
    energy = np.arange(15).reshape(5, 3)+1.
    image, lost = camera.splat_xyz_samples(pixel, energy, 10, 10)
    np.testing.assert_allclose(image.sum(axis=(0, 1))+lost, energy.sum(axis=0), atol=1e-13)
    assert np.all(lost > 0)
    projected = camera.project_m_to_pixels([[0, 0, 0]], [0, 0, 1], meters_per_pixel=2, width=10, height=20, roll_deg=0)
    np.testing.assert_allclose(projected, [[5, 10.6]])


def test_invalid_intervals_and_superluminal_patterns_are_rejected():
    with pytest.raises(ValueError):
        camera.periodic_interval_weights([0], [-1], 1., 10)
    with pytest.raises(ValueError):
        camera.periodic_interval_weights([0], [2], 1., 10)
    with pytest.raises(ValueError, match='light cylinder'):
        camera.rotating_pattern_arrival_edges([0], [1], [[camera.C_M_S, 0, 0]], [0, 1, 0], 1)


def test_empty_overlap_keeps_requested_schema():
    result=camera.periodic_interval_weights(np.array([]),np.array([]),1.,10,return_edges=True)
    assert len(result)==5
    assert all(len(a)==0 for a in result)
