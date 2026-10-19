"""Conservative optical-event arrival bins and physical image-plane samples.

No radiation, beaming, field lines or species colors are invented here. The
input is already calculated observer-directed XYZ power per physical lepton.
Finite source intervals and physical number weights become detector-bin energy
and then radiance at one shared physical scale. Spatial quadrature and camera
sampling remain separate convergence questions.
"""
import numpy as np

from surface_camera import C_M_S, camera_basis, rotate_vectors, unit_vector


def rotating_pattern_arrival_edges(source_time_s, source_interval_s, position_m,
                                  observer_direction, omega_rad_s, spin_axis=(0, 0, 1)):
    """Map the ends of a corotating-pattern cell interval to arrival time.

    This is a prescribed rotating emissivity-cell world tube, not an electron
    trajectory. Its lab-frame power must already include any radiative boosts.
    Plasma clocks/particle velocity cannot be replaced by this pattern speed.
    Flat-spacetime far-observer propagation is explicit: tau=t-n.r/c.
    """
    t, dt = np.broadcast_arrays(np.asarray(source_time_s, float), np.asarray(source_interval_s, float))
    position = np.asarray(position_m, float)
    if position.shape != t.shape+(3,) or np.any(~np.isfinite(position)):
        raise ValueError('One finite source position is required per interval')
    if np.any(~np.isfinite(t)) or np.any(~np.isfinite(dt)) or np.any(dt <= 0):
        raise ValueError('Source intervals must be positive and finite')
    if not np.isfinite(omega_rad_s):
        raise ValueError('Nonfinite pattern angular speed')
    n, axis = unit_vector(observer_direction), unit_vector(spin_axis)
    velocity = omega_rad_s*np.cross(axis, position)
    if np.any(np.linalg.norm(velocity, axis=-1) >= C_M_S):
        raise ValueError('Rotating cells must remain inside their light cylinder')
    left_r = rotate_vectors(position, axis, -omega_rad_s*dt/2)
    right_r = rotate_vectors(position, axis, omega_rad_s*dt/2)
    left = t-dt/2-np.sum(left_r*n, axis=-1)/C_M_S
    right = t+dt/2-np.sum(right_r*n, axis=-1)/C_M_S
    if np.any(right <= left):
        raise ValueError('Nonmonotonic arrival interval')
    return left, right


def periodic_interval_weights(left_s, right_s, period_s, bins, *, return_edges=False):
    """Sparse interval/bin overlaps; weights sum to one per source event.

    Positive intervals may straddle the cycle origin. The fractions distribute
    the event's already source-time-integrated energy. Dividing by the detector
    bin duration occurs once, after accumulation; no additional arrival Jacobian
    is applied. Power is approximated as constant over each supplied interval.
    """
    left, right = np.broadcast_arrays(np.asarray(left_s, float), np.asarray(right_s, float))
    if left.ndim != 1 or not np.isfinite(period_s) or period_s <= 0:
        raise ValueError('Expected one-dimensional intervals and a positive period')
    if not isinstance(bins, int) or bins < 1:
        raise ValueError('Positive integer detector bin count required')
    length = right-left
    if np.any(~np.isfinite(left)) or np.any(~np.isfinite(right)) or np.any(length <= 0) or np.any(length > period_s*(1+1e-12)):
        raise ValueError('Each interval must have length in (0,period]')
    width = period_s/bins
    # Bring the left edge into one cycle before integer conversion.
    phase_left = np.remainder(left, period_s)
    phase_right = phase_left+length
    first = np.floor(phase_left/width).astype(np.int64)
    last = np.floor(np.nextafter(phase_right, -np.inf)/width).astype(np.int64)
    event_ids, target_bins, weights, overlap_left, overlap_right = [], [], [], [], []
    span = int(np.max(last-first, initial=-1))+1
    for step in range(span):
        index = first+step
        valid = index <= last
        overlap = np.maximum(0., np.minimum(phase_right, (index+1)*width)-np.maximum(phase_left, index*width))
        valid &= overlap > 0
        ids = np.flatnonzero(valid)
        event_ids.append(ids)
        target_bins.append(np.remainder(index[valid], bins))
        weights.append(overlap[valid]/length[valid])
        if return_edges:
            overlap_left.append(left[valid]+np.maximum(phase_left[valid], index[valid]*width)-phase_left[valid])
            overlap_right.append(left[valid]+np.minimum(phase_right[valid], (index[valid]+1)*width)-phase_left[valid])
    if not event_ids:
        result = np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([], dtype=float)
        return result+(np.array([], dtype=float), np.array([], dtype=float)) if return_edges else result
    result = np.concatenate(event_ids), np.concatenate(target_bins), np.concatenate(weights)
    if return_edges:
        return result+(np.concatenate(overlap_left), np.concatenate(overlap_right))
    return result


def inverse_rotating_arrival(arrival_s, source_center_s, position_center_m,
                            observer_direction, omega_rad_s, spin_axis=(0, 0, 1)):
    """Invert a subluminal rotating-pattern arrival map with Newton's method."""
    arrival, center = np.broadcast_arrays(np.asarray(arrival_s, float), np.asarray(source_center_s, float))
    position = np.asarray(position_center_m, float)
    if position.shape != arrival.shape+(3,) or np.any(~np.isfinite(arrival)) or np.any(~np.isfinite(center)):
        raise ValueError('Invalid source/arrival shapes or coordinates')
    n, axis = unit_vector(observer_direction), unit_vector(spin_axis)
    velocity = omega_rad_s*np.cross(axis, position)
    if not np.isfinite(omega_rad_s) or np.any(np.linalg.norm(velocity, axis=-1) >= C_M_S):
        raise ValueError('Pattern must remain subluminal')
    estimate = center+(arrival-center+np.sum(position*n, axis=-1)/C_M_S)/(1-np.sum(velocity*n, axis=-1)/C_M_S)
    for _ in range(6):
        actual = rotate_vectors(position, axis, omega_rad_s*(estimate-center))
        derivative = 1-np.sum(omega_rad_s*np.cross(axis, actual)*n, axis=-1)/C_M_S
        estimate -= (estimate-np.sum(actual*n, axis=-1)/C_M_S-arrival)/derivative
    actual = rotate_vectors(position, axis, omega_rad_s*(estimate-center))
    residual = estimate-np.sum(actual*n, axis=-1)/C_M_S-arrival
    if np.any(np.abs(residual) > 2e-12*np.maximum(1, np.abs(arrival))):
        raise ValueError('Arrival inverse did not converge')
    return estimate


def source_xyz_energy(power_xyz_per_sr, number_weights, source_interval_s):
    power = np.asarray(power_xyz_per_sr, float)
    weights, dt = np.broadcast_arrays(np.asarray(number_weights, float), np.asarray(source_interval_s, float))
    if power.shape != weights.shape+(3,) or np.any(~np.isfinite(power)) or np.any(power < 0):
        raise ValueError('Expected nonnegative observer-directed physical XYZ power')
    if np.any(~np.isfinite(weights)) or np.any(weights < 0) or np.any(~np.isfinite(dt)) or np.any(dt <= 0):
        raise ValueError('Invalid physical number weights or source intervals')
    return power*weights[..., None]*dt[..., None]


def project_m_to_pixels(position_m, observer_direction, *, meters_per_pixel,
                        width, height, roll_deg=-25., center_fraction=(.5, .53)):
    position = np.asarray(position_m, float)
    if position.shape[-1:] != (3,) or np.any(~np.isfinite(position)):
        raise ValueError('Physical positions must be finite three-vectors')
    if not np.isfinite(meters_per_pixel) or meters_per_pixel <= 0:
        raise ValueError('Positive physical pixel size required')
    if width < 1 or height < 1 or not np.isfinite(roll_deg):
        raise ValueError('Invalid camera dimensions or roll')
    ex, ey, _ = camera_basis(observer_direction, (0, 0, 1))
    x, y = position@ex, position@ey
    angle = np.deg2rad(roll_deg)
    x, y = np.cos(angle)*x-np.sin(angle)*y, np.sin(angle)*x+np.cos(angle)*y
    return np.stack((width*center_fraction[0]+x/meters_per_pixel,
                     height*center_fraction[1]-y/meters_per_pixel), axis=-1)


def splat_xyz_samples(pixel_xy, xyz_energy, width, height):
    """Conservative CIC rasterization of point quadrature samples.

    Returns the in-frame image and lost XYZ energy. CIC is a pixel reconstruction
    footprint, not particle diffusion or a physical radiation-cone width. It
    cannot make an undersampled source field scientifically converged.
    """
    pixel = np.asarray(pixel_xy, float)
    energy = np.asarray(xyz_energy, float)
    if pixel.ndim != 2 or pixel.shape[1] != 2 or energy.shape != (len(pixel), 3):
        raise ValueError('Invalid pixel/energy dimensions')
    if np.any(~np.isfinite(pixel)) or np.any(~np.isfinite(energy)) or np.any(energy < 0):
        raise ValueError('Nonfinite or negative physical samples')
    base = np.floor(pixel).astype(np.int64)
    fraction = pixel-base
    output = np.zeros((height*width, 3), float)
    lost = np.zeros(3)
    for dx, dy in ((0, 0), (1, 0), (0, 1), (1, 1)):
        ij = base+np.array([dx, dy])
        factor = (fraction[:, 0] if dx else 1-fraction[:, 0])*(fraction[:, 1] if dy else 1-fraction[:, 1])
        valid = (ij[:, 0] >= 0) & (ij[:, 0] < width) & (ij[:, 1] >= 0) & (ij[:, 1] < height)
        np.add.at(output, ij[valid, 1]*width+ij[valid, 0], factor[valid, None]*energy[valid])
        lost += np.sum(factor[~valid, None]*energy[~valid], axis=0)
    return output.reshape(height, width, 3), lost
