"""Conservative radiance aggregation for geometrically subpixel cooling cells.

This is a rendering approximation, not a new particle or radiation model.
Cells may share a support only along a bounded path through their ordered
reference node and eight corners. Actual directional XYZ powers
times physical lepton counts are summed at the same source time and interval.
No source phases, spectra, directions, or source clocks are averaged.
"""
import numpy as np


C_M_S = 299792458.
_SIGNS = np.array([[x, y, z] for x in (-1., 1.)
                   for y in (-1., 1.) for z in (-1., 1.)])


def _geometry(positions, vertices):
    positions = np.asarray(positions, dtype=np.float64)
    vertices = np.asarray(vertices, dtype=np.float64)
    if (positions.ndim != 2 or positions.shape[1:] != (3,)
            or vertices.shape != (len(positions), 8, 3)
            or not np.isfinite(positions).all() or not np.isfinite(vertices).all()):
        raise ValueError('Finite reference positions [cells,3] and corners [cells,8,3] required')
    return positions, vertices


def _affine_corners(vertices):
    center = vertices.mean(axis=1)
    half = np.einsum('vi,nvj->nji', _SIGNS, vertices - center[:, None, :]) / 8
    return center[:, None, :] + np.einsum('vi,nji->nvj', _SIGNS, half)


def compact_reference_geometry(reference_positions_m, reference_vertices_m, grid_m):
    """Return representative geometry, input-to-group mapping, and error ledger.

    Representatives are existing cells (the first input cell in each group),
    rather than averaged synthetic cells. A path step is the greatest of the
    reference-node displacement and the eight ordered-corner displacements.
    Cumulative path-length bins of width grid_m define consecutive groups.
    The triangle inequality bounds each raw node/corner movement by grid_m.
    """
    position, vertices = _geometry(reference_positions_m, reference_vertices_m)
    if not np.isfinite(grid_m) or grid_m <= 0:
        raise ValueError('A finite positive geometry grid is required')
    if len(position):
        step = np.maximum(np.linalg.norm(np.diff(position, axis=0), axis=-1),
                          np.linalg.norm(np.diff(vertices, axis=0), axis=-1).max(axis=1))
        # Extended-precision accumulation keeps small cooling steps resolved
        # after the much larger jumps between separate injection paths.
        path = np.r_[np.longdouble(0), np.cumsum(step, dtype=np.longdouble)]
        scaled = path / np.longdouble(grid_m)
        if not np.isfinite(scaled).all() or np.max(scaled) >= 9e18:
            raise ValueError('Geometry path/grid ratio cannot be represented by int64 bins')
        labels = np.floor(scaled).astype(np.int64)
        boundary = np.r_[True, np.diff(labels) != 0]
        representative = np.flatnonzero(boundary)
        inverse = np.cumsum(boundary, dtype=np.int64) - 1
        path_length = float(path[-1])
    else:
        representative = inverse = np.array([], dtype=np.int64)
        path_length = 0.
    compact_position, compact_vertices = position[representative], vertices[representative]
    if len(position):
        node_error = float(np.max(np.linalg.norm(position - compact_position[inverse], axis=-1)))
        corner_error = float(np.max(np.linalg.norm(vertices - compact_vertices[inverse], axis=-1)))
        fitted = _affine_corners(vertices)
        affine_error = float(np.max(np.linalg.norm(fitted - fitted[representative][inverse], axis=-1)))
        if max(node_error, corner_error) > grid_m * (1 + 1e-10):
            raise ValueError('Ordered-path grouping exceeded its raw geometry displacement bound')
    else:
        node_error = corner_error = affine_error = 0.
    report = {'input_reference_cells': len(position), 'output_reference_cells': len(representative),
              'geometry_grid_m': float(grid_m),
              'membership': 'Consecutive cells in one floor(cumulative_path/grid_m) bin; path steps are max(node movement, all eight ordered-corner movements)',
              'ordered_reference_path_length_m': path_length,
              'representative': 'First existing reference cell in each bounded-path group; no synthetic geometry average',
              'maximum_node_displacement_m': node_error,
              'maximum_raw_corner_displacement_m': corner_error,
              'maximum_affine_corner_displacement_m': affine_error,
              'raw_geometry_displacement_bound_m': float(grid_m),
              'affine_geometry_displacement_bound_m': float(1.5 * grid_m),
              'maximum_node_arrival_shift_bound_s': node_error / C_M_S,
              'maximum_affine_corner_arrival_shift_bound_s': affine_error / C_M_S,
              'scope': 'Subpixel image support aggregation; source particle/radiation calculation unchanged'}
    return compact_position, compact_vertices, inverse.astype(np.int32), report


def compact_render_events(*, reference_positions_m, reference_vertices_m,
                          cell_index, source_time_s, source_interval_s,
                          number_weights, xyz_lm_sr, grid_m, period_s,
                          phase_index=None, positions_m=None):
    """Aggregate event powers for shared source-phase and subpixel supports.

    Return ``(arrays, report)`` with the ordinary event-camera array names.
    ``xyz_lm_sr`` in the returned arrays is total directional XYZ power for the
    grouped physical population, and its internal number weight is exactly 1.
    If phase identifiers are supplied, each identifier must have a single
    *exactly identical* source midpoint and source interval. Different phases
    never merge, even if their positions or source times happen to coincide.
    """
    reference, vertices = _geometry(reference_positions_m, reference_vertices_m)
    ci = np.asarray(cell_index)
    time = np.asarray(source_time_s, dtype=np.float64)
    dt = np.asarray(source_interval_s, dtype=np.float64)
    number = np.asarray(number_weights, dtype=np.float64)
    power = np.asarray(xyz_lm_sr, dtype=np.float64)
    count = len(time)
    if (time.ndim != 1 or dt.shape != (count,) or number.shape != (count,)
            or power.shape != (count, 3) or ci.shape != (count,)
            or not np.issubdtype(ci.dtype, np.integer)
            or np.any(ci < 0) or np.any(ci >= len(reference))
            or not all(np.isfinite(x).all() for x in (time, dt, number, power))
            or np.any(dt <= 0) or np.any(number < 0) or np.any(power < 0)
            or not np.isfinite(period_s) or period_s <= 0):
        raise ValueError('Invalid event arrays, cell identifiers, or source period')
    if phase_index is None:
        _, pi = np.unique(time, return_inverse=True)
    else:
        pi = np.asarray(phase_index)
        if pi.shape != (count,) or not np.issubdtype(pi.dtype, np.integer):
            raise ValueError('Integer phase identifiers must match the event count')
    _, phase_first, phase_inverse = np.unique(pi, return_index=True, return_inverse=True)
    if (np.any(time != time[phase_first][phase_inverse])
            or np.any(dt != dt[phase_first][phase_inverse])):
        raise ValueError('Every phase must have one exact source midpoint and interval; clocks cannot be averaged')
    omega = 2 * np.pi / period_s
    angle = omega * time
    co, si = np.cos(angle), np.sin(angle)
    expected = np.stack((co * reference[ci, 0] - si * reference[ci, 1],
                         si * reference[ci, 0] + co * reference[ci, 1], reference[ci, 2]), axis=1)
    if positions_m is not None:
        given = np.asarray(positions_m, dtype=np.float64)
        tolerance = 32 * np.finfo(float).eps * np.maximum(np.linalg.norm(expected, axis=-1), 1.)
        if (given.shape != expected.shape or not np.isfinite(given).all()
                or np.any(np.linalg.norm(given - expected, axis=-1) > tolerance)):
            raise ValueError('Only the prescribed rigidly rotating reference-cell geometry may be compacted')
    compact_reference, compact_vertices, mapping, geometry_report = compact_reference_geometry(
        reference, vertices, grid_m)
    # Dense phase indices avoid dtype promotion (e.g. uint64 identifiers mixed
    # with int32 cell indices) changing the exact phase grouping.
    event_keys = np.column_stack((phase_inverse, mapping[ci]))
    _, event_first, event_inverse = np.unique(event_keys, axis=0, return_index=True, return_inverse=True)
    grouped = np.zeros((len(event_first), 3), dtype=np.float64)
    weighted = power * number[:, None]
    np.add.at(grouped, event_inverse, weighted)
    output_ci = mapping[ci[event_first]]
    output_time = time[event_first]
    output_angle = omega * output_time
    co, si = np.cos(output_angle), np.sin(output_angle)
    r = compact_reference[output_ci]
    output_position = np.stack((co * r[:, 0] - si * r[:, 1],
                                si * r[:, 0] + co * r[:, 1], r[:, 2]), axis=1)
    before = (weighted * dt[:, None]).sum(axis=0)
    after = (grouped * dt[event_first, None]).sum(axis=0)
    absolute_error = np.abs(after - before)
    fractional_error = np.divide(absolute_error, np.abs(before), out=np.zeros(3), where=before != 0)
    if np.any(fractional_error > 1e-10) or not np.isfinite(grouped).all():
        raise ValueError('XYZ source-energy aggregation failed its roundoff-level conservation check')
    arrays = {'reference_cell_positions_m': compact_reference,
              'cell_footprint_vertices_reference_m': compact_vertices,
              'cell_index': output_ci.astype(np.int32), 'phase_index': pi[event_first],
              'source_time_s': output_time, 'source_interval_s': dt[event_first],
              'number_weights': np.ones(len(event_first)), 'xyz_lm_sr': grouped,
              'positions_m': output_position}
    report = {**geometry_report, 'input_events': count, 'output_events': len(event_first),
              'input_source_xyz_energy': before.tolist(), 'output_source_xyz_energy': after.tolist(),
              'source_xyz_energy_absolute_error': absolute_error.tolist(),
              'source_xyz_energy_relative_error': fractional_error.tolist(),
              'event_aggregation': 'Sum of each original directional XYZ power times physical number weight, at exactly the same source phase/time/interval',
              'internal_weight_convention': 'Grouped power already contains the physical population; number_weights equals one',
              'source_phases_averaged': False, 'directional_powers_recomputed': False,
              'reference_geometry': geometry_report}
    return arrays, report
