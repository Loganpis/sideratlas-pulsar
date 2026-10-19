#!/usr/bin/env python3
"""Independent diagnostics of raw, staggered local PIC exports.

These are constraints and data-integrity measurements, not a Vela validation,
optical camera, kinetic-convergence certificate or complete energy ledger.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def staggered_divergences(fields):
    """Yee divB at cell centres and divE at vertices; ghost border discarded.

    Raw contravariant E/B have a common physical factor B0*dx. Derivatives
    with respect to physical code positions cancel dx, so returned differences
    are divergence divided by B0, not derivatives of collocated vectors.
    """
    b = fields[..., 3:6]
    e = fields[..., :3]
    divb = sum(np.roll(b[..., i], -1, axis=i)-b[..., i] for i in range(3))
    dive = sum(e[..., i]-np.roll(e[..., i], 1, axis=i) for i in range(3))
    return divb, dive


def csv_rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def inspect_snapshot(meta_path):
    meta_path = Path(meta_path)
    m = json.loads(meta_path.read_text())
    if m['storage_order'] != 'C: i,j,k,component' or m['dtype'] != 'little-endian float64':
        raise ValueError('Unsupported explicit native array format')
    path = meta_path.with_name(m['fields_file'])
    fields = np.fromfile(path, dtype='<f8').reshape(m['field_shape'])
    if not np.isfinite(fields).all():
        raise ValueError('Nonfinite native field values')
    g, n = int(m['ghost_cells_each_side']), np.array(m['active_resolution'], int)
    h = float(m['dx'])
    if g < 1 or not np.all(n > 4) or h <= 0:
        raise ValueError('Invalid native mesh metadata')
    interior = tuple(slice(g, g+int(v)) for v in n)
    axes = [(np.arange(int(v))+offset)*h+origin for v, origin, offset in
            zip(n, m['origin'], (0., 0., 0.))]
    xyz_vertex = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
    # Vertex and cell-centre constraints live at DIFFERENT positions.
    rv = np.linalg.norm(xyz_vertex, axis=-1)
    rb = np.linalg.norm(xyz_vertex+h/2, axis=-1)
    radius = float(m['stellar_radius'])
    maskb = (rb > radius+3*h) & (rb < 1.8)
    maske = (rv > radius+3*h) & (rv < 1.8)
    if not maskb.any() or not maske.any():
        raise ValueError('No cells in the declared diagnostic shell')
    divb, dive = [x[interior] for x in staggered_divergences(fields)]
    rho = fields[interior+(9,)]
    db = divb[maskb]
    de, charge = dive[maske], rho[maske]
    b_raw = fields[interior+(slice(3, 6),)]
    b_ref = np.sqrt(np.mean(np.sum((h*b_raw[maskb])**2, axis=-1)))
    gauss = de-charge
    denom = np.linalg.norm(de)+np.linalg.norm(charge)
    reference_rms = denom/np.sqrt(len(de))
    # A ratio with a roundoff-sized denominator is meaningless (e.g. before
    # the injected charge reaches this shell). Retain absolute residuals.
    roundoff_floor = 128*np.finfo(float).eps*max(1., float(np.max(np.abs(fields[..., :6]))))
    relative_defined = reference_rms > roundoff_floor
    particles_path = meta_path.with_name(m['particles_file'])
    rows = csv_rows(particles_path)
    ids = [(int(r['species']), int(r['id'])) for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate stable species/particle IDs')
    values = np.array([[float(r[c]) for c in ('x', 'y', 'z', 'ux', 'uy', 'uz', 'weight')]
                       for r in rows], float).reshape(-1, 7)
    if not np.isfinite(values).all() or np.any(values[:, -1] < 0):
        raise ValueError('Invalid particle coordinates, momenta or weights')
    gamma = np.sqrt(1+np.sum(values[:, 3:6]**2, axis=-1))
    counts = {str(s): sum(a == s for a, _ in ids) for s in (1, 2)}
    currents = fields[interior+(slice(6, 9),)]
    diagnostics = {
        'step': m['completed_step'], 'position_and_E_time': m['position_and_E_time'],
        'momentum_B_and_J_time': m['momentum_B_and_J_time'],
        'mesh': m['active_resolution'], 'dx': h,
        'alive_particles': len(rows), 'species_counts': counts,
        'gamma_min': float(gamma.min()) if len(gamma) else None,
        'gamma_max': float(gamma.max()) if len(gamma) else None,
        'gamma_mean': float(gamma.mean()) if len(gamma) else None,
        'divB_rms_over_B0': float(np.sqrt(np.mean(db*db))),
        'divB_rms_over_rms_B_per_rlc': float(np.sqrt(np.mean(db*db))/b_ref) if b_ref else None,
        'gauss_rms_over_B0': float(np.sqrt(np.mean(gauss*gauss))),
        'gauss_reference_rms_over_B0': float(reference_rms),
        'relative_gauss_reference_roundoff_floor_over_B0': float(roundoff_floor),
        'gauss_relative_l2_symmetric': float(np.linalg.norm(gauss)/denom) if relative_defined else None,
        'gauss_relative_defined_above_roundoff': bool(relative_defined),
        'charge_rms_over_B0': float(np.sqrt(np.mean(charge*charge))),
        'raw_current_l2': float(np.linalg.norm(currents)),
        'no_surface_inside_alive_particles': bool(np.all(np.linalg.norm(values[:, :3], axis=-1) >= radius)),
        'field_sha256': sha(path), 'particles_sha256': sha(particles_path),
        'metadata_sha256': sha(meta_path),
    }
    return diagnostics, ids


def event_inventory(events):
    """Require unique births/deaths and reject silently overwritten event IDs."""
    birth, death = {}, {}
    for row in events:
        key = int(row['species']), int(row['id'])
        step = int(row['completed_step'])
        kind = row['event']
        if key[0] not in (1, 2) or key[1] < 0 or step < 0:
            raise ValueError('Invalid stable particle event key or step')
        target = birth if kind == 'birth' else death if kind in ('star_death', 'boundary_death') else None
        if target is None:
            raise ValueError('Unknown particle event type')
        if key in target:
            raise ValueError('Duplicate birth/death event for stable particle ID')
        target[key] = step
    if any(key not in birth or step < birth[key] for key, step in death.items()):
        raise ValueError('Particle death has no preceding birth')
    return birth, death


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('native', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    paths = sorted(a.native.glob('state-*.json'))
    if not paths:
        raise FileNotFoundError('No native state exports exist yet')
    results, inventories = [], []
    for path in paths:
        measured, ids = inspect_snapshot(path)
        results.append(measured)
        inventories.append(set(ids))
    events_path = a.native/'events.csv'
    events = csv_rows(events_path)
    birth, death = event_inventory(events)
    event_checks = []
    for report, recorded in zip(results, inventories):
        step = report['step']
        expected = {key for key, born in birth.items() if born <= step and death.get(key, step+1) > step}
        event_checks.append(expected == recorded)
    output = {
        'scope': 'Raw export integrity and discrete constraints of a reduced local kinetic prototype; not a Vela model or optical image',
        'finite_native_exports': True, 'snapshots': results,
        'relative_gauss_definition': 'norm(divE-rho)/(norm(divE)+norm(rho)); null when the RMS reference is below an explicitly recorded roundoff floor; not a convergence tolerance',
        'particle_inventory_matches_birth_death_events': all(event_checks),
        'particle_inventory_checks_per_snapshot': event_checks,
        'events_sha256': sha(events_path), 'script_sha256': sha(__file__),
        'source_and_boundary_interventions_require_separate_energy_accounting': True,
        'complete_energy_conservation_test': False, 'kinetic_convergence_established': False,
        'physical_vela_scaling_established': False, 'optical_prediction': False,
    }
    destination = a.output or a.native/'independent-audit.json'
    destination.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))
    if not all(event_checks):
        raise SystemExit('Stable particle inventory does not match recorded events')


if __name__ == '__main__':
    main()
