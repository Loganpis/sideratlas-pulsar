import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from compact_render_geometry import compact_reference_geometry, compact_render_events
from volume_raster import SIGNS


def inputs():
    reference = np.array([[100., 20., 10.], [100.01, 20.01, 10.01], [300., 20., 10.]])
    vertices = reference[:, None, :] + SIGNS[None, :, :] * 2
    return dict(reference_positions_m=reference, reference_vertices_m=vertices,
                cell_index=np.array([0, 1, 2, 0, 1]),
                source_time_s=np.array([0., 0., 0., .5, .5]),
                source_interval_s=np.full(5, .2),
                number_weights=np.array([2., 3., 4., 5., 7.]),
                xyz_lm_sr=np.arange(15, dtype=float).reshape(5, 3) + 1,
                phase_index=np.array([0, 0, 0, 1, 1]), grid_m=.1, period_s=1.)


def test_xyz_population_energy_is_summed_without_color_averaging():
    data = inputs()
    arrays, report = compact_render_events(**data)
    assert report['input_events'] == 5 and report['output_events'] == 3
    assert np.array_equal(arrays['number_weights'], np.ones(3))
    before = (data['xyz_lm_sr'] * data['number_weights'][:, None]
              * data['source_interval_s'][:, None]).sum(axis=0)
    after = (arrays['xyz_lm_sr'] * arrays['source_interval_s'][:, None]).sum(axis=0)
    np.testing.assert_allclose(before, after, rtol=1e-15)
    first = np.flatnonzero((arrays['source_time_s'] == 0) & (arrays['positions_m'][:, 0] < 200))[0]
    np.testing.assert_array_equal(arrays['xyz_lm_sr'][first],
                                  data['xyz_lm_sr'][0] * 2 + data['xyz_lm_sr'][1] * 3)


def test_reference_geometry_error_and_retarded_time_bound():
    data = inputs()
    _, _, mapping, report = compact_reference_geometry(
        data['reference_positions_m'], data['reference_vertices_m'], data['grid_m'])
    assert mapping[0] == mapping[1] and mapping[0] != mapping[2]
    assert 0 < report['maximum_node_displacement_m'] <= data['grid_m']
    assert report['maximum_raw_corner_displacement_m'] <= data['grid_m']
    assert report['maximum_affine_corner_displacement_m'] <= report['affine_geometry_displacement_bound_m']
    assert report['maximum_node_arrival_shift_bound_s'] == report['maximum_node_displacement_m'] / 299792458.


def test_different_source_phases_never_merge_even_at_identical_time():
    data = inputs()
    data['source_time_s'][:] = 0
    arrays, _ = compact_render_events(**data)
    assert len(arrays['source_time_s']) == 3
    assert set(arrays['phase_index']) == {0, 1}


@pytest.mark.parametrize('field', ['source_time_s', 'source_interval_s'])
def test_mismatched_source_clock_within_phase_is_rejected(field):
    data = inputs()
    data[field][1] += .001
    with pytest.raises(ValueError, match='exact source midpoint and interval'):
        compact_render_events(**data)


def test_position_reconstruction_preserves_rotation_and_rejects_other_motion():
    data = inputs()
    reference = data['reference_positions_m'][data['cell_index']]
    angle = 2 * np.pi * data['source_time_s']
    co, si = np.cos(angle), np.sin(angle)
    data['positions_m'] = np.stack((co * reference[:, 0] - si * reference[:, 1],
                                   si * reference[:, 0] + co * reference[:, 1], reference[:, 2]), axis=1)
    arrays, report = compact_render_events(**data)
    assert arrays['positions_m'][-1, 0] < 0
    assert report['source_phases_averaged'] is False
    data['positions_m'][0, 0] += 1
    with pytest.raises(ValueError, match='rigidly rotating'):
        compact_render_events(**data)


def test_unstated_phase_identifiers_use_exact_source_times():
    data = inputs()
    del data['phase_index']
    arrays, _ = compact_render_events(**data)
    assert len(arrays['source_time_s']) == 3


def test_invalid_grid_is_rejected():
    data = inputs()
    data['grid_m'] = 0
    with pytest.raises(ValueError, match='positive geometry grid'):
        compact_render_events(**data)


def test_ordered_path_bounds_corners_even_when_motion_changes_direction():
    position = np.c_[np.sin(np.arange(100) * .1), np.cos(np.arange(100) * .1), np.zeros(100)]
    vertices = position[:, None, :] + SIGNS[None, :, :] * 2
    # Independent corner motion must contribute to the path metric too.
    vertices[:, 0, 2] += .08 * np.sin(np.arange(100) * .3)
    r, v, mapping, report = compact_reference_geometry(position, vertices, .4)
    assert len(r) < len(position)
    assert np.max(np.linalg.norm(position - r[mapping], axis=-1)) <= .4
    assert np.max(np.linalg.norm(vertices - v[mapping], axis=-1)) <= .4
    assert report['raw_geometry_displacement_bound_m'] == .4
    # A repeated physical location reached much later along the path is not
    # silently merged across a source-path excursion.
    assert mapping[0] != mapping[63]
