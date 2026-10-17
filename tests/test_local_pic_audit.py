"""Manufactured Yee fields and corrupt particle logs test the native audit."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest


spec = importlib.util.spec_from_file_location(
    'local_pic_audit', Path(__file__).parents[1]/'research/pic-local/audit_native.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_discrete_curl_has_no_magnetic_monopoles():
    rng = np.random.default_rng(413)
    potential = rng.normal(size=(11, 13, 15, 3))
    # Edge-centred A -> face-centred B -> cell-centred divB. For example,
    # A_z(i,j,k+.5) and A_y(i,j+.5,k) curl to Bx(i,j+.5,k+.5).
    forward = lambda a, axis: np.roll(a, -1, axis=axis)-a
    ax, ay, az = np.moveaxis(potential, -1, 0)
    magnetic = np.stack((forward(az, 1)-forward(ay, 2),
                         forward(ax, 2)-forward(az, 0),
                         forward(ay, 0)-forward(ax, 1)), axis=-1)
    fields = np.zeros((11, 13, 15, 10))
    fields[..., 3:6] = magnetic
    divb, _ = audit.staggered_divergences(fields)
    np.testing.assert_allclose(divb, 0, atol=3e-15)


def test_linear_electric_field_has_known_uniform_charge():
    i, j, k = np.meshgrid(np.arange(11), np.arange(13), np.arange(15), indexing='ij')
    fields = np.zeros((11, 13, 15, 10))
    fields[..., 0], fields[..., 1], fields[..., 2] = 2*(i+.5), -3*(j+.5), 5*(k+.5)
    _, dive = audit.staggered_divergences(fields)
    # Exclude the synthetic periodic wrap, which is not part of the linear test.
    np.testing.assert_allclose(dive[1:, 1:, 1:], 4, atol=1e-15)


def event(kind, step=1, particle=2):
    return {'event': kind, 'completed_step': str(step), 'species': '1', 'id': str(particle)}


def test_event_logs_preserve_unique_particle_lifetimes():
    birth, death = audit.event_inventory([event('birth'), event('star_death', 2)])
    assert birth == {(1, 2): 1} and death == {(1, 2): 2}
    with pytest.raises(ValueError, match='Duplicate'):
        audit.event_inventory([event('birth'), event('birth', 2)])
    with pytest.raises(ValueError, match='Duplicate'):
        audit.event_inventory([event('birth'), event('star_death', 2), event('boundary_death', 3)])


def test_event_logs_reject_impossible_deaths_and_unknown_events():
    with pytest.raises(ValueError, match='preceding birth'):
        audit.event_inventory([event('star_death')])
    with pytest.raises(ValueError, match='preceding birth'):
        audit.event_inventory([event('birth', 3), event('star_death', 2)])
    with pytest.raises(ValueError, match='Unknown'):
        audit.event_inventory([event('teleport')])
