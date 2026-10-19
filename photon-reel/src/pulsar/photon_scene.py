"""Display transport of photons from the existing rotating emission model.

All lengths are R_LC and all times are Omega*t, so c=1. A rotating
source is sampled at emission, then each photon direction is held fixed.
These functions do not solve or modify the magnetosphere or particle model.
"""
import numpy as np
from .snapshot_fields import rotate_z
from .emission_quadrature import rebin_phase_averages

TAU = 2*np.pi


def source_frame(position, direction, emission_time):
    return rotate_z(position, -emission_time), rotate_z(direction, -emission_time)


def emission_event(position0, direction0, emission_time):
    return rotate_z(position0, emission_time), rotate_z(direction0, emission_time)


def propagate(position, direction, age):
    """Straight free-streaming rays; no post-emission corotation."""
    return np.asarray(position) + np.asarray(age)[..., None]*np.asarray(direction)


def observer_event(position0, direction0):
    """Rotate a source event until its photon azimuth is the observer's zero.

    The existing observer integrates a finite polar bin about zeta=65deg.
    Each direction in that bin uses the asymptotic phase t - n.r, exactly
    as in the original model, rather than a nearby point-detector distance.
    """
    te = -np.arctan2(direction0[..., 1], direction0[..., 0])
    r, n = emission_event(position0, direction0, te)
    return te, r, n, (te-np.sum(n*r, axis=-1))/TAU % 1


def wavefront_crossing(position, direction, emission_time, distance):
    """Cross n.x=D; subtracting D gives the model's distant arrival time."""
    age = distance-np.sum(position*direction, axis=-1)
    return propagate(position, direction, age), emission_time+age


def received_flux(edges, profile, lab_phase, distance, shutter_phase):
    """Exact frame average of saved phase bins, including flight-time delay."""
    center = np.asarray(lab_phase)-distance/TAU
    return rebin_phase_averages(edges, profile,
        np.atleast_1d(center-shutter_phase/2), np.atleast_1d(center+shutter_phase/2))


def systematic_sample(weights, count, offset=.5):
    """Deterministic equal-weight packets, without perturbing ray directions."""
    weights=np.asarray(weights, float)
    if np.any(weights<0) or not np.isfinite(weights).all() or weights.sum()<=0:
        raise ValueError('Finite positive total photon weight required')
    if count<1 or not 0<=offset<1:raise ValueError('Invalid packet sampling')
    cumulative=np.cumsum(weights);total=cumulative[-1]
    return np.searchsorted(cumulative, (np.arange(count)+offset)*total/count), total/count
