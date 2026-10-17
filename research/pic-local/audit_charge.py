#!/usr/bin/env python3
"""Independently reconstruct exported CIC vertex charge from particle CSVs.

This checks export consistency and charge bookkeeping, not PIC continuity,
Gauss's law, convergence, physical Vela normalization, or optical emission.
It reproduces the declared first-order CIC vertex convention using only
positions, species charges, number weights and recorded normalization factors.
No native particle/field files or producer source are changed.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import tempfile
import tomllib
import unittest

import numpy as np


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 16), b''):
            digest.update(block)
    return digest.hexdigest()


def reconstruct_vertex_charge(positions, signed_number_weights, origin, dx,
                              stored_shape, ghost_cells, coefficient):
    """Return CIC density and explicit lost-support accounting.

    The vertex corresponding to local cell index i is stored at i+ghost_cells.
    Position xi=(x-origin)/dx contributes prod(1-f) or prod(f) to the eight
    vertices at floor(xi)+corner.  coefficient=q0/(B0*V0), so the returned
    density has exactly the export's rho_vertex_gauss_units convention.
    Signed number weights already include the dimensionless species charge.
    """
    positions = np.asarray(positions, float)
    signed = np.asarray(signed_number_weights, float)
    origin = np.asarray(origin, float)
    shape = np.asarray(stored_shape, int)
    if (positions.ndim != 2 or positions.shape[1] != 3 or signed.shape != (len(positions),)
            or origin.shape != (3,) or shape.shape != (3,) or np.any(shape < 2)
            or not isinstance(ghost_cells, (int, np.integer)) or ghost_cells < 0
            or not np.isfinite(dx) or dx <= 0 or not np.isfinite(coefficient) or coefficient <= 0
            or np.any(~np.isfinite(positions)) or np.any(~np.isfinite(signed))
            or np.any(~np.isfinite(origin))):
        raise ValueError('Invalid finite particle, grid or charge-normalization inputs')
    xi = (positions-origin)/dx
    if np.any(np.abs(xi) > np.iinfo(np.int64).max/2):
        raise ValueError('Particle coordinates cannot be represented as grid indices')
    cells = np.floor(xi).astype(np.int64)
    fraction = xi-cells
    rho = np.zeros(tuple(shape), dtype=float)
    lost_signed, lost_absolute = 0., 0.
    unsupported = np.zeros(len(positions), dtype=bool)
    for ax in (0, 1):
        for ay in (0, 1):
            for az in (0, 1):
                corner = np.array([ax, ay, az])
                index = cells+corner+ghost_cells
                support = np.prod(np.where(corner, fraction, 1-fraction), axis=1)
                contribution = coefficient*signed*support
                valid = np.all((index >= 0) & (index < shape), axis=1)
                np.add.at(rho, tuple(index[valid].T), contribution[valid])
                lost_signed += float(np.sum(contribution[~valid]))
                lost_absolute += float(np.sum(np.abs(contribution[~valid])))
                unsupported |= (~valid) & (support > 0) & (signed != 0)
    expected = coefficient*float(np.sum(signed))
    return rho, {
        'expected_full_density_sum_from_particles': expected,
        'reconstructed_retained_density_sum': float(np.sum(rho)),
        'lost_support_signed_density_sum': lost_signed,
        'lost_support_absolute_density_sum': lost_absolute,
        'particles_with_nonzero_charge_support_outside_storage': int(np.count_nonzero(unsupported)),
        'all_nonzero_charge_shape_support_inside_storage': not bool(np.any(unsupported)),
        'retained_plus_lost_minus_expected': float(np.sum(rho))+lost_signed-expected,
        'absolute_particle_charge_density_sum': coefficient*float(np.sum(np.abs(signed))),
    }


def comparison_stats(reconstructed, exported, selection, atol, rtol):
    computed = reconstructed[selection]
    recorded = exported[selection]
    if not computed.size:
        return {'vertices': 0, 'comparison_available': False}
    difference = computed-recorded
    error_l2 = float(np.linalg.norm(difference.ravel()))
    reference_l2 = float(np.linalg.norm(recorded.ravel()))
    symmetric_l2 = reference_l2+float(np.linalg.norm(computed.ravel()))
    return {
        'vertices': int(computed.size), 'comparison_available': True,
        'maximum_absolute_error': float(np.max(np.abs(difference))),
        'rms_error': float(np.sqrt(np.mean(difference*difference))),
        'relative_l2_to_exported': error_l2/reference_l2 if reference_l2 else (0. if not error_l2 else None),
        'symmetric_relative_l2': error_l2/symmetric_l2 if symmetric_l2 else 0.,
        'reconstructed_density_sum': float(np.sum(computed)),
        'exported_density_sum': float(np.sum(recorded)),
        'maximum_absolute_exported_density': float(np.max(np.abs(recorded))),
        'pointwise_within_recorded_tolerance': bool(np.allclose(computed, recorded, atol=atol, rtol=rtol)),
    }


def read_particles(path, species_charge):
    position, signed, counts = [], [], {}
    with Path(path).open() as stream:
        reader = csv.DictReader(stream)
        required = {'species', 'id', 'x', 'y', 'z', 'weight'}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError('Particle CSV lacks species, IDs, positions or number weights')
        identifiers = set()
        for row in reader:
            species, identifier = int(row['species']), int(row['id'])
            key = (species, identifier)
            if key in identifiers:
                raise ValueError('Duplicate stable particle ID in snapshot')
            identifiers.add(key)
            if species not in species_charge:
                raise ValueError(f'Species {species} has no recorded charge mapping')
            weight = float(row['weight'])
            if not np.isfinite(weight) or weight < 0:
                raise ValueError('Particle number weight must be finite and nonnegative')
            position.append([float(row[key]) for key in ('x', 'y', 'z')])
            signed.append(weight*species_charge[species])
            counts[species] = counts.get(species, 0)+1
    return np.asarray(position, float).reshape(-1, 3), np.asarray(signed, float), counts


def inspect_charge_snapshot(meta_path, species_charge, *, atol=1e-10, rtol=1e-9,
                            maximum_stored_vertices=5_000_000):
    meta_path = Path(meta_path)
    meta = json.loads(meta_path.read_text())
    if (meta['storage_order'] != 'C: i,j,k,component'
            or meta['dtype'] != 'little-endian float64'
            or meta['field_components'][9] != 'rho_vertex_gauss_units'
            or not meta['rho_definition'].startswith('CIC vertex')):
        raise ValueError('Unsupported explicit native CIC charge/array convention')
    shape = np.asarray(meta['field_shape'], int)
    active = np.asarray(meta['active_resolution'], int)
    ghosts = int(meta['ghost_cells_each_side'])
    if (shape.shape != (4,) or shape[-1] != 10 or active.shape != (3,) or ghosts < 0
            or np.any(active < 2) or not np.all(shape[:3] == active+2*ghosts)):
        raise ValueError('Native array dimensions do not match active grid and ghosts')
    if int(np.prod(shape[:3])) > maximum_stored_vertices:
        raise ValueError('Native grid exceeds this bounded audit size; explicitly raise --maximum-stored-vertices')
    for key in ('dx', 'q0', 'B0', 'V0'):
        if not np.isfinite(meta[key]) or meta[key] <= 0:
            raise ValueError(f'Invalid positive native normalization: {key}')
    fields_path = meta_path.with_name(meta['fields_file'])
    particles_path = meta_path.with_name(meta['particles_file'])
    if fields_path.stat().st_size != int(np.prod(shape))*8:
        raise ValueError('Binary field size disagrees with metadata')
    fields = np.memmap(fields_path, dtype='<f8', mode='r', shape=tuple(shape))
    exported = fields[..., 9]
    if not np.isfinite(exported).all():
        raise ValueError('Nonfinite exported charge density')
    positions, signed, counts = read_particles(particles_path, species_charge)
    coefficient = meta['q0']/(meta['B0']*meta['V0'])
    reconstructed, support = reconstruct_vertex_charge(
        positions, signed, meta['origin'], meta['dx'], shape[:3], ghosts, coefficient)
    active_selection = tuple(slice(ghosts, ghosts+int(n)) for n in active)
    strict_selection = tuple(slice(ghosts+1, ghosts+int(n)-1) for n in active)
    regions = {
        'full_stored_vertices_including_ghosts': comparison_stats(
            reconstructed, exported, np.s_[:, :, :], atol, rtol),
        'all_active_vertices_excluding_ghosts': comparison_stats(
            reconstructed, exported, active_selection, atol, rtol),
        'strict_interior_active_vertices_one_vertex_margin': comparison_stats(
            reconstructed, exported, strict_selection, atol, rtol),
    }
    # The same source-excluded shell as audit_native.py, on VERTEX positions.
    axes = [np.arange(int(n))*meta['dx']+origin for n, origin in zip(active, meta['origin'])]
    xv, yv, zv = np.meshgrid(*axes, indexing='ij', sparse=True)
    radius2 = xv*xv+yv*yv+zv*zv
    lower = meta['stellar_radius']+3*meta['dx']
    upper = 1.8/float(meta.get('final_omega', 1.))
    shell = (radius2 > lower*lower) & (radius2 < upper*upper)
    regions['exterior_diagnostic_shell_active_vertices'] = comparison_stats(
        reconstructed[active_selection], exported[active_selection], shell, atol, rtol)
    support.update({
        'particle_signed_number_weight_sum': float(np.sum(signed)),
        'particle_total_charge_code_units': meta['q0']*float(np.sum(signed)),
        'exported_total_charge_code_units': meta['B0']*meta['V0']*float(np.sum(exported)),
        'exported_full_density_sum': float(np.sum(exported)),
        'exported_plus_lost_minus_particle_density_sum':
            float(np.sum(exported))+support['lost_support_signed_density_sum']
            - support['expected_full_density_sum_from_particles'],
    })
    return {
        'completed_step': meta['completed_step'],
        'position_time': meta['position_and_E_time'],
        'particles': len(positions), 'species_counts': counts,
        'coefficient_q0_over_B0_V0': coefficient,
        'comparison_regions': regions, 'shape_support_and_global_charge': support,
        'shell_vertex_radii_code': [lower, upper],
        'metadata_sha256': file_hash(meta_path),
        'fields_sha256': file_hash(fields_path), 'particles_sha256': file_hash(particles_path),
    }


def select_configuration(native, requested=None):
    if requested is not None:
        return Path(requested)
    choices = [path for name in ('startup.toml', 'one-rotation.toml')
               if (path := Path(native).parent/name).is_file()]
    if len(choices) != 1:
        raise ValueError('Select the actual run configuration with --configuration')
    return choices[0]


def audit(native, configuration, *, atol=1e-10, rtol=1e-9, maximum_stored_vertices=5_000_000):
    native, configuration = Path(native), Path(configuration)
    data = tomllib.loads(configuration.read_text())
    species_charge = {index: float(species['charge'])
                      for index, species in enumerate(data['particles']['species'], start=1)}
    if not species_charge or not all(np.isfinite(value) for value in species_charge.values()):
        raise ValueError('Run configuration does not define finite species charges')
    paths = sorted(native.glob('state-*.json'))
    if not paths:
        raise FileNotFoundError('No native snapshots exist')
    reports = [inspect_charge_snapshot(path, species_charge, atol=atol, rtol=rtol,
                maximum_stored_vertices=maximum_stored_vertices) for path in paths]
    producer = native.parent/'pgen.hpp'
    source_warning = None
    if producer.is_file():
        source_warning = bool(re.search(r'sp\.i[123]\(p\)\s*\+\s*sp\.dx[123]\(p\)', producer.read_text()))
    return {
        'scope': 'Independent first-order CIC vertex export reconstruction and total-charge bookkeeping; not a charge-continuity, Gauss-law, physical-Vela or optical validation',
        'native_directory': str(native.resolve()),
        'configuration': str(configuration.resolve()), 'configuration_sha256': file_hash(configuration),
        'script_sha256': file_hash(__file__),
        'producer_sha256': file_hash(producer) if producer.is_file() else None,
        'potential_float32_coordinate_addition_in_recorded_producer': source_warning,
        'coordinate_precision_note': 'If int particle indices are added to float32 subcell offsets before promotion to double, exported Cartesian positions lose subcell bits and pointwise reconstructed charge can differ; this is an export-consistency issue, not evidence of a physical charge source',
        'species_charge_from_recorded_configuration': species_charge,
        'charge_convention': 'rho_vertex_gauss_units = CIC sum(q_species*number_weight)*q0/(B0*V0)',
        'density_units': 'Native dimensionless Gauss units; no conversion to coulombs per cubic metre',
        'total_charge_code_units': 'q0*sum(q_species*number_weight) = B0*V0*sum(rho) when all shape support is retained',
        'active_vertex_convention': 'Exactly indices ghost..ghost+n_active-1; upper-face vertices can lie in the stored ghost region',
        'required_kernel_contract': 'First-order vertex CIC, matching the actual SHAPE_ORDER=1 build; this audit does not authenticate arbitrary higher-order builds',
        'tolerance': {'absolute_native_density': atol, 'relative_pointwise': rtol},
        'snapshots': reports,
        'all_snapshots_pointwise_match_full_storage': all(
            row['comparison_regions']['full_stored_vertices_including_ghosts']['pointwise_within_recorded_tolerance']
            for row in reports),
        'all_snapshots_have_retained_full_charge_shape_support': all(
            row['shape_support_and_global_charge']['all_nonzero_charge_shape_support_inside_storage']
            for row in reports),
        'physical_charge_conservation_or_convergence_established': False,
    }


class ChargeAuditTests(unittest.TestCase):
    def test_one_particle_cic_vertex_and_charge_sum(self):
        rho, support = reconstruct_vertex_charge([[1.25, 2.5, 3.75]], [-4.], [0, 0, 0],
                                                  1., [9, 9, 9], 2, .3125)
        self.assertAlmostEqual(float(rho.sum()), -1.25)
        self.assertEqual(np.count_nonzero(rho), 8)
        self.assertAlmostEqual(float(rho[3, 4, 5]), -1.25*.75*.5*.25)
        self.assertTrue(support['all_nonzero_charge_shape_support_inside_storage'])

    def test_opposite_species_cancel_and_number_weights_are_linear(self):
        position = [[1.2, 2.3, 3.4], [1.2, 2.3, 3.4]]
        rho, _ = reconstruct_vertex_charge(position, [-17., 17.], [0, 0, 0], 1., [9]*3, 2, .25)
        np.testing.assert_array_equal(rho, 0.)
        one, _ = reconstruct_vertex_charge(position[:1], [-17.], [0, 0, 0], 1., [9]*3, 2, .25)
        split, _ = reconstruct_vertex_charge(position*2, [-4.25]*4, [0, 0, 0], 1., [9]*3, 2, .25)
        np.testing.assert_allclose(one, split, rtol=1e-15, atol=0.)

    def test_ghost_support_and_lost_support_are_accounted_separately(self):
        rho, support = reconstruct_vertex_charge([[-.2, 1., 1.]], [2.], [0, 0, 0], 1., [8]*3, 2, 1.)
        self.assertAlmostEqual(float(rho.sum()), 2.)
        self.assertAlmostEqual(float(rho[2:6, 2:6, 2:6].sum()), 1.6)
        self.assertTrue(support['all_nonzero_charge_shape_support_inside_storage'])
        rho, support = reconstruct_vertex_charge([[-2.2, 1., 1.]], [-2.], [0, 0, 0], 1., [8]*3, 2, 1.)
        self.assertAlmostEqual(float(rho.sum()), -1.6)
        self.assertAlmostEqual(support['lost_support_signed_density_sum'], -.4)
        self.assertAlmostEqual(support['lost_support_absolute_density_sum'], .4)
        self.assertAlmostEqual(support['retained_plus_lost_minus_expected'], 0.)

    def test_native_coordinates_and_normalization_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            position = np.array([[.1, .2, .3]])
            rho, _ = reconstruct_vertex_charge(position, [-3.], [-1.]*3, .25, [12]*3, 2, .125)
            fields = np.zeros((12, 12, 12, 10)); fields[..., 9] = rho
            fields.astype('<f8').tofile(root/'fields.bin')
            (root/'particles.csv').write_text('species,id,x,y,z,ux,uy,uz,weight\n1,7,.1,.2,.3,0,0,0,3\n')
            meta = {'storage_order': 'C: i,j,k,component', 'dtype': 'little-endian float64',
                'field_shape': [12, 12, 12, 10], 'active_resolution': [8]*3, 'ghost_cells_each_side': 2,
                'field_components': ['x']*9+['rho_vertex_gauss_units'], 'rho_definition': 'CIC vertex fixture',
                'dx': .25, 'q0': .25, 'B0': 2., 'V0': 1., 'origin': [-1.]*3,
                'fields_file': 'fields.bin', 'particles_file': 'particles.csv', 'stellar_radius': .1,
                'completed_step': 1, 'position_and_E_time': .2, 'final_omega': 1.}
            path = root/'state-000001.json'; path.write_text(json.dumps(meta))
            result = inspect_charge_snapshot(path, {1: -1.})
            self.assertTrue(result['comparison_regions']['full_stored_vertices_including_ghosts']
                            ['pointwise_within_recorded_tolerance'])
            self.assertAlmostEqual(result['shape_support_and_global_charge']['particle_total_charge_code_units'], -.75)

    def test_invalid_inputs_and_empty_arrays(self):
        with self.assertRaises(ValueError):
            reconstruct_vertex_charge([[np.nan, 0, 0]], [1.], [0]*3, 1., [9]*3, 2, 1.)
        rho, support = reconstruct_vertex_charge(np.empty((0, 3)), [], [0]*3, 1., [9]*3, 2, 1.)
        self.assertEqual(float(rho.sum()), 0.)
        self.assertTrue(support['all_nonzero_charge_shape_support_inside_storage'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native', type=Path, nargs='?')
    parser.add_argument('--configuration', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--atol', type=float, default=1e-10)
    parser.add_argument('--rtol', type=float, default=1e-9)
    parser.add_argument('--maximum-stored-vertices', type=int, default=5_000_000)
    parser.add_argument('--require-match', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.atol < 0 or args.rtol < 0 or not np.isfinite(args.atol+args.rtol):
        parser.error('Comparison tolerances must be finite and nonnegative')
    if args.self_test:
        result = unittest.TextTestRunner(verbosity=2).run(
            unittest.defaultTestLoader.loadTestsFromTestCase(ChargeAuditTests))
        if not result.wasSuccessful():
            raise SystemExit(1)
    if args.native is None:
        if args.self_test:
            return
        parser.error('Provide the native snapshot directory')
    configuration = select_configuration(args.native, args.configuration)
    report = audit(args.native, configuration, atol=args.atol, rtol=args.rtol,
                   maximum_stored_vertices=args.maximum_stored_vertices)
    output = args.output or args.native.parent/'independent-charge-audit.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'snapshots': len(report['snapshots']),
        'all_snapshots_pointwise_match_full_storage': report['all_snapshots_pointwise_match_full_storage'],
        'all_shape_support_inside_storage': report['all_snapshots_have_retained_full_charge_shape_support'],
        'potential_float32_coordinate_addition_in_recorded_producer':
            report['potential_float32_coordinate_addition_in_recorded_producer'],
        'largest_pointwise_native_density_error': max(
            row['comparison_regions']['full_stored_vertices_including_ghosts']['maximum_absolute_error']
            for row in report['snapshots'])}, indent=2))
    if args.require_match and not report['all_snapshots_pointwise_match_full_storage']:
        raise SystemExit('Exported charge and CSV reconstruction differ beyond the recorded tolerance')


if __name__ == '__main__':
    main()
