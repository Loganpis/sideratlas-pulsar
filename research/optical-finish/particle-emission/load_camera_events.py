"""Reader for either exact full or compact physical photon event archives."""
from pathlib import Path

import numpy as np
from scipy.constants import c


def load_camera_events(path):
    """Return ordinary ndarray dictionary compatible with the original schema.

    Compact XYZ has float32 storage but is promoted here; its measured weighted
    error bound lives in the model report. Geometric/time reconstruction is
    double precision and deterministic, using the same Rz expressions as the
    producer. No spatial, spectral, angular or time smoothing is introduced.
    Per-event wavelength-bin spectra are intentionally absent in a compact
    file; aggregate per-phase spectra and each event's total band power remain.
    """
    with np.load(Path(path), allow_pickle=False) as saved:
        d = {key: saved[key] for key in saved.files}
    if not bool(d.get('compact_camera', False)):
        return d
    ci, pi = d['cell_index'], d['phase_index']
    reference = d['reference_cell_positions_m'][ci]
    angle = d['phase_angles_rad'][pi]
    co, si = np.cos(angle), np.sin(angle)
    position = np.stack((co*reference[:, 0]-si*reference[:, 1],
                         si*reference[:, 0]+co*reference[:, 1], reference[:, 2]), axis=1)
    period = float(d['period_s'])
    omega = 2*np.pi/period
    d['positions_m'] = position
    d['source_time_s'] = angle/omega
    d['source_interval_s'] = np.full(len(ci), float(d['source_interval_scalar_s']))
    d['number_weights'] = d['reference_cell_number_weights'][ci]
    d['pattern_beta'] = omega*np.stack((-position[:, 1], position[:, 0], np.zeros(len(ci))), axis=1)/c
    d['xyz_lm_sr'] = d['xyz_lm_sr'].astype(np.float64)
    return d
