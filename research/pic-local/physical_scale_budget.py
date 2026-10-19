#!/usr/bin/env python3
"""Estimate explicit uniform PIC scale budgets without allocating a grid.

This is dimensional analysis, not a simulation, hardware measurement, or claim
that every numerical method must use a uniform grid.  The default numbers are
explicit Vela reference/model inputs, not new measurements.  AMR, implicit
methods, guiding centres and multiscale models have different requirements.
"""
from dataclasses import asdict, dataclass
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import tomllib
import unittest


C_M_S = 299_792_458.0
E_C = 1.602_176_634e-19
EPSILON_0_F_M = 8.854_187_8128e-12
M_E_KG = 9.109_383_7015e-31


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a finite positive number')
    return value


def positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return value


def gauss_to_tesla(magnetic_field_gauss):
    return positive(magnetic_field_gauss, 'Surface magnetic field')*1e-4


def goldreich_julian_density_m3(omega_rad_s, magnetic_field_t):
    """Leading aligned reference |rho_GJ|/e = 2 epsilon_0 Omega B/e.

    This treats the local field as parallel to rotation, omits GR and current
    corrections, and is not an actual density solution near the light cylinder.
    """
    return (2*EPSILON_0_F_M*positive(omega_rad_s, 'Omega')
            * positive(magnetic_field_t, 'Magnetic field')/E_C)


def cold_pair_plasma_scales(total_lepton_density_m3):
    """Cold pair-plasma omega_p and c/omega_p for total e- plus e+ density.

    Both species contribute q^2/m.  Density is their sum, not the density of
    pairs multiplied by two a second time.  Relativistic anisotropic response
    cannot generally be represented by replacing m_e with one mean gamma*m_e.
    """
    density = positive(total_lepton_density_m3, 'Total lepton number density')
    omega = math.sqrt(density*E_C**2/(EPSILON_0_F_M*M_E_KG))
    return {'omega_p_rad_s': omega, 'frequency_p_hz': omega/(2*math.pi),
            'cold_skin_depth_m': C_M_S/omega}


@dataclass(frozen=True)
class Inputs:
    period_s: float = .08937108859019084
    stellar_radius_km: float = 13.
    surface_reference_field_gauss: float = 3e12
    total_lepton_multiplicity: float = 6000.
    illustrative_gamma: float = 1.
    domain_width_rlc: float = 4.
    cells_per_stellar_radius: int = 10
    cells_per_skin_depth: int = 2
    illustrative_particles_per_cell_total: int = 8
    local_grid_cells_per_axis: int = 48
    local_stellar_radius_rlc: float = .5
    local_domain_width_rlc: float = 5.

    def validate(self):
        for name in ('period_s', 'stellar_radius_km', 'surface_reference_field_gauss',
                     'total_lepton_multiplicity', 'domain_width_rlc',
                     'local_stellar_radius_rlc', 'local_domain_width_rlc'):
            positive(getattr(self, name), name)
        if (isinstance(self.illustrative_gamma, bool)
                or not math.isfinite(self.illustrative_gamma) or self.illustrative_gamma < 1):
            raise ValueError('Illustrative gamma must be finite and at least one')
        for name in ('cells_per_stellar_radius', 'cells_per_skin_depth',
                     'illustrative_particles_per_cell_total', 'local_grid_cells_per_axis'):
            positive_integer(getattr(self, name), name)
        if 2*self.local_stellar_radius_rlc >= self.local_domain_width_rlc:
            raise ValueError('Local numerical star must fit inside the comparison domain')


def uniform_grid_budget(domain_width_m, maximum_cell_size_m, particles_per_cell_total=8):
    """Illustrative memory lower bound for one cubic uniform grid, no allocation.

    Fields: six E/B plus three current components, each float64, 72 bytes/cell.
    Particles: six phase-space coordinates, one number weight and one uint64 ID,
    64 bytes/particle.  Actual Entity storage includes additional arrays.
    """
    width = positive(domain_width_m, 'Domain width')
    requested_dx = positive(maximum_cell_size_m, 'Maximum cell size')
    positive_integer(particles_per_cell_total, 'Illustrative particles per cell')
    n = math.ceil(width/requested_dx)
    cells = n**3
    actual_dx = width/n
    fields = cells*9*8
    particle_count = cells*particles_per_cell_total
    particles = particle_count*(7*8+8)
    return {
        'maximum_requested_cell_size_m': requested_dx,
        'actual_uniform_cell_size_m': actual_dx,
        'cells_per_axis': n, 'cells_total': cells,
        'field_and_current_bytes_lower_bound': fields,
        'illustrative_particles_total': particle_count,
        'particle_phase_space_bytes_lower_bound': particles,
        'combined_bytes_lower_bound': fields+particles,
        'combined_decimal_terabytes_lower_bound': (fields+particles)/1e12,
        'vacuum_3d_yee_cfl_timestep_s': actual_dx/(math.sqrt(3)*C_M_S),
    }


def calculate_budget(inputs=Inputs()):
    inputs.validate()
    omega = 2*math.pi/inputs.period_s
    rlc = C_M_S/omega
    radius = inputs.stellar_radius_km*1e3
    field = gauss_to_tesla(inputs.surface_reference_field_gauss)
    ratio = radius/rlc
    # Same-axis vacuum-dipole r^-3 scaling, only a reference magnitude.
    # For a polar B_star the vacuum equatorial B at RLC is half this value.
    lc_field = field*ratio**3
    star_gj = goldreich_julian_density_m3(omega, field)
    lc_gj = goldreich_julian_density_m3(omega, lc_field)
    star_density = inputs.total_lepton_multiplicity*star_gj
    lc_density = inputs.total_lepton_multiplicity*lc_gj
    star_scales = cold_pair_plasma_scales(star_density)
    lc_scales = cold_pair_plasma_scales(lc_density)
    domain_width = inputs.domain_width_rlc*rlc
    star_dx = radius/inputs.cells_per_stellar_radius
    lc_skin_dx = lc_scales['cold_skin_depth_m']/inputs.cells_per_skin_depth
    surface_skin_dx = star_scales['cold_skin_depth_m']/inputs.cells_per_skin_depth
    grids = {
        'resolve_stellar_radius_only': uniform_grid_budget(
            domain_width, star_dx, inputs.illustrative_particles_per_cell_total),
        'resolve_reference_lc_skin_depth': uniform_grid_budget(
            domain_width, min(star_dx, lc_skin_dx), inputs.illustrative_particles_per_cell_total),
        'resolve_reference_surface_skin_depth': uniform_grid_budget(
            domain_width, min(star_dx, surface_skin_dx), inputs.illustrative_particles_per_cell_total),
    }
    for grid in grids.values():
        grid['minimum_yee_steps_per_rotation'] = math.ceil(
            inputs.period_s/grid['vacuum_3d_yee_cfl_timestep_s'])
    local_dx = inputs.local_domain_width_rlc*rlc/inputs.local_grid_cells_per_axis
    local_radius = inputs.local_stellar_radius_rlc*rlc
    local = {
        'grid_cells_per_axis': inputs.local_grid_cells_per_axis,
        'comparison_domain_width_rlc': inputs.local_domain_width_rlc,
        'numerical_stellar_radius_rlc': inputs.local_stellar_radius_rlc,
        'numerical_radius_if_period_mapped_to_vela_km': local_radius/1e3,
        'radius_enlargement_relative_to_assumed_physical_star': local_radius/radius,
        'cell_size_if_period_mapped_to_vela_m': local_dx,
        'numerical_stellar_radius_in_cells': local_radius/local_dx,
        'physical_stellar_radius_in_same_cells': radius/local_dx,
        'cell_size_over_reference_lc_cold_skin_depth': local_dx/lc_scales['cold_skin_depth_m'],
        'cell_size_over_reference_surface_cold_skin_depth': local_dx/star_scales['cold_skin_depth_m'],
        'physical_omega_p_dt_not_inferred': True,
        'physical_magnetic_field_and_density_mapping_validated': False,
        'claim': 'A reduced numerical pulsar can be a real coupled PIC calculation while failing to reproduce the physical Vela scale ratios',
    }
    return {
        'scope': 'Dimensional scale and illustrative explicit uniform-grid budget; not a run or a hardware survey',
        'inputs': asdict(inputs),
        'input_status': 'Period reference and explicit model choices; R, B, multiplicity and gamma are not newly measured or independently inferred',
        'constants_si': {'c_m_s': C_M_S, 'elementary_charge_c': E_C,
                         'epsilon0_f_m': EPSILON_0_F_M, 'electron_mass_kg': M_E_KG},
        'equations': {
            'omega': 'Omega = 2*pi/P', 'light_cylinder': 'R_LC = c/Omega',
            'aligned_reference_number_density': 'n_GJ = 2*epsilon_0*Omega*B/e',
            'total_lepton_density': 'n_total(e- plus e+) = kappa_total*n_GJ',
            'cold_plasma_frequency': 'omega_p = sqrt(n_total*e^2/(epsilon_0*m_e))',
            'cold_skin_depth': 'lambda_skin = c/omega_p',
            'lc_reference_field': 'B_LC,reference = B_star*(R_star/R_LC)^3',
            'uniform_grid': 'N_axis = ceil(L_domain/dx_max); N_cells = N_axis^3',
        },
        'physical_reference_scales': {
            'omega_rad_s': omega, 'light_cylinder_radius_m': rlc,
            'light_cylinder_radius_km': rlc/1e3, 'stellar_radius_over_rlc': ratio,
            'surface_reference_field_t': field, 'lc_dipole_reference_field_t': lc_field,
            'surface_cold_cyclotron_frequency_rad_s': E_C*field/M_E_KG,
            'surface_cold_cyclotron_to_rotation_frequency_ratio': E_C*field/(M_E_KG*omega),
            'surface_aligned_gj_reference_density_m3': star_gj,
            'lc_formal_aligned_gj_reference_density_m3': lc_gj,
            'surface_total_lepton_reference_density_m3': star_density,
            'lc_total_lepton_reference_density_m3': lc_density,
            'surface_cold_pair_plasma': star_scales,
            'lc_cold_pair_plasma': lc_scales,
        },
        'gamma_illustration_only': {
            'gamma': inputs.illustrative_gamma,
            'surface_skin_if_effective_inertia_is_gamma_m_e_m':
                star_scales['cold_skin_depth_m']*math.sqrt(inputs.illustrative_gamma),
            'lc_skin_if_effective_inertia_is_gamma_m_e_m':
                lc_scales['cold_skin_depth_m']*math.sqrt(inputs.illustrative_gamma),
            'used_in_uniform_grid_budgets': False,
            'qualification': 'sqrt(gamma) scaling assumes a specific transverse effective inertia; anisotropic relativistic distributions and parallel streaming have different response, so cold budgets remain explicitly gamma=1',
        },
        'uniform_grid_budgets': grids,
        'memory_assumptions': {
            'field_components_float64': 9, 'bytes_per_float64': 8,
            'particle_coordinates_and_weight_float64': 7, 'id_bytes': 8,
            'particles_per_cell_total': inputs.illustrative_particles_per_cell_total,
            'excluded': ['ghost cells', 'field backups and auxiliary buffers', 'particle previous positions and tags',
                         'sorting and deposition buffers', 'communication', 'output', 'runtime overhead'],
            'qualification': 'Illustrative storage lower bounds for the specified uniform layout; no actual machine RAM or runtime prediction',
        },
        'local_reduced_model_comparison': local,
        'limitations': [
            'Aligned GJ magnitude is a reference charge-density scale, not the self-consistent oblique density or multiplicity',
            'The r^-3 LC reference is not an evolved force-free/PIC field; polar versus equatorial dipole conventions add a factor of two',
            'Cold pair skin depths use total lepton density and do not describe a solved relativistic momentum distribution',
            'Two cells per skin depth and ten cells per stellar radius are illustrative minimum choices, not convergence certificates',
            'Uniform-grid and vacuum-Yee estimates do not describe AMR, implicit, guiding-centre or multiscale algorithms',
            'Matching period alone cannot restore the physical star-to-LC, gyro, plasma, radiative-cooling and optical-emission scales',
            'A computed reduced PIC prototype must not be described as a validated optical view or atmosphere of Vela',
        ],
        'references': [
            {'citation': 'Goldreich & Julian 1969, Pulsar Electrodynamics',
             'url': 'https://doi.org/10.1086/150119'},
            {'citation': 'Brambilla et al. 2018, Electron–Positron Pair Flow and Current Composition in the Pulsar Magnetosphere',
             'url': 'https://arxiv.org/abs/1710.03536'},
            {'citation': 'Kalapotharakos et al. 2018, 3D Kinetic Pulsar Magnetosphere Models',
             'url': 'https://arxiv.org/abs/1710.03170'},
        ],
    }


def file_hash(path):
    digest = sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4096), b''):
            digest.update(chunk)
    return digest.hexdigest()


def native_reference_from_config(path, inputs):
    """Record actual native parameters; this does not assert an SI mapping.

    Entity's checked source defines omegaB0=1/larmor0.  Its SR pusher uses
    (charge/mass)*omegaB0 multiplying the dimensionless magnetic field.  The
    reported gyro ratio uses polar_field as a reference amplitude, not a
    measurement of the smoothed grid's field or a relativistic trajectory.
    """
    path = Path(path).resolve()
    data = tomllib.loads(path.read_text())
    omega = positive(data['setup']['omega'], 'Native target Omega')
    larmor0 = positive(data['scales']['larmor0'], 'Native larmor0')
    polar = positive(data['setup']['polar_field'], 'Native polar field')
    star = positive(data['setup']['star_radius'], 'Native stellar radius')
    widths = [positive(upper-lower, 'Native domain width')
              for lower, upper in data['grid']['extent']]
    resolution = data['grid']['resolution']
    gyro = []
    for species in data['particles']['species']:
        mass = positive(species['mass'], 'Native species mass')
        charge = abs(species['charge'])
        gyro.append({'species': species['label'],
                     'cold_cyclotron_to_target_rotation_reference': charge/mass*polar/(larmor0*omega)})
    return {
        'config_path': str(path), 'config_sha256': file_hash(path),
        'larmor0_code': larmor0, 'omegaB0_code': 1/larmor0,
        'target_omega_code': omega, 'polar_field_reference_code': polar,
        'stellar_radius_rlc': star*omega,
        'domain_widths_rlc': [width*omega for width in widths],
        'resolution': resolution, 'species_reference_ratios': gyro,
        'comparison_geometry_matches_recorded_config': (
            all(value == inputs.local_grid_cells_per_axis for value in resolution)
            and all(math.isclose(width*omega, inputs.local_domain_width_rlc) for width in widths)
            and math.isclose(star*omega, inputs.local_stellar_radius_rlc)),
        'qualification': 'Reference unit-amplitude cold gyro ratio at final target rotation; spin-up Omega varies, local B and particle gamma vary, and no physical magnetic-field mapping is inferred',
    }


class ScaleBudgetTests(unittest.TestCase):
    """Small independent unit and scaling checks; no grid arrays are allocated."""
    def test_gauss_tesla_and_light_cylinder(self):
        self.assertEqual(gauss_to_tesla(3e12), 3e8)
        scales = calculate_budget()['physical_reference_scales']
        self.assertAlmostEqual(scales['light_cylinder_radius_km'], 4264.2, delta=.1)
        self.assertAlmostEqual(scales['stellar_radius_over_rlc'], 13/4264.2, delta=1e-7)

    def test_gj_si_against_independent_cgs_expression(self):
        period, gauss = .08937108859019084, 3e12
        # In Gaussian units n_GJ=B/(c*e*P), in cm^-3; convert to m^-3.
        c_cgs, e_statc = 2.99792458e10, 4.803204712570263e-10
        cgs = gauss/(c_cgs*e_statc*period)*1e6
        si = goldreich_julian_density_m3(2*math.pi/period, gauss_to_tesla(gauss))
        self.assertAlmostEqual(si/cgs, 1., delta=1e-9)

    def test_cold_plasma_units_and_density_scaling(self):
        one_per_cm3 = cold_pair_plasma_scales(1e6)
        self.assertAlmostEqual(one_per_cm3['frequency_p_hz'], 8978.66, delta=.01)
        high = cold_pair_plasma_scales(16e6)
        self.assertAlmostEqual(high['cold_skin_depth_m']/one_per_cm3['cold_skin_depth_m'], .25)

    def test_grid_ceil_memory_and_no_particle_double_count(self):
        grid = uniform_grid_budget(100., 10., particles_per_cell_total=8)
        self.assertEqual(grid['cells_per_axis'], 10)
        self.assertEqual(grid['cells_total'], 1000)
        self.assertEqual(grid['field_and_current_bytes_lower_bound'], 72000)
        self.assertEqual(grid['particle_phase_space_bytes_lower_bound'], 512000)
        self.assertEqual(grid['combined_bytes_lower_bound'], 584000)
        self.assertEqual(uniform_grid_budget(100., 9.9)['cells_per_axis'], 11)

    def test_local_enlargement_is_explicit(self):
        result = calculate_budget()
        local = result['local_reduced_model_comparison']
        self.assertAlmostEqual(local['numerical_radius_if_period_mapped_to_vela_km'], 2132.1, delta=.1)
        self.assertAlmostEqual(local['numerical_stellar_radius_in_cells'], 4.8)
        self.assertLess(local['physical_stellar_radius_in_same_cells'], .04)
        self.assertFalse(local['physical_magnetic_field_and_density_mapping_validated'])

    def test_multiplicity_and_gamma_are_not_hidden_pair_factors(self):
        first = calculate_budget(Inputs(total_lepton_multiplicity=1.))
        four = calculate_budget(Inputs(total_lepton_multiplicity=4., illustrative_gamma=9.))
        s1 = first['physical_reference_scales']['lc_cold_pair_plasma']['cold_skin_depth_m']
        s4 = four['physical_reference_scales']['lc_cold_pair_plasma']['cold_skin_depth_m']
        self.assertAlmostEqual(s4/s1, .5)
        illustrative = four['gamma_illustration_only']['lc_skin_if_effective_inertia_is_gamma_m_e_m']
        self.assertAlmostEqual(illustrative/s4, 3.)
        self.assertFalse(four['gamma_illustration_only']['used_in_uniform_grid_budgets'])

    def test_invalid_inputs_are_rejected(self):
        for inputs in (Inputs(period_s=0), Inputs(stellar_radius_km=-1),
                       Inputs(surface_reference_field_gauss=math.nan),
                       Inputs(total_lepton_multiplicity=0), Inputs(illustrative_gamma=.9),
                       Inputs(local_grid_cells_per_axis=1.5)):
            with self.assertRaises(ValueError):
                calculate_budget(inputs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--period-s', type=float, default=Inputs.period_s)
    parser.add_argument('--radius-km', type=float, default=Inputs.stellar_radius_km)
    parser.add_argument('--field-gauss', type=float, default=Inputs.surface_reference_field_gauss)
    parser.add_argument('--multiplicity', type=float, default=Inputs.total_lepton_multiplicity)
    parser.add_argument('--gamma', type=float, default=Inputs.illustrative_gamma)
    parser.add_argument('--domain-width-rlc', type=float, default=Inputs.domain_width_rlc)
    parser.add_argument('--local-domain-width-rlc', type=float, default=Inputs.local_domain_width_rlc)
    parser.add_argument('--local-grid', type=int, default=Inputs.local_grid_cells_per_axis)
    parser.add_argument('--local-star-rlc', type=float, default=Inputs.local_stellar_radius_rlc)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[2]/'outputs/pic-local/physical-scale-budget.json')
    parser.add_argument('--native-config', type=Path,
                        default=Path(__file__).resolve().parent/'startup.toml')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ScaleBudgetTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            raise SystemExit(1)
        return
    inputs = Inputs(period_s=args.period_s, stellar_radius_km=args.radius_km,
        surface_reference_field_gauss=args.field_gauss,
        total_lepton_multiplicity=args.multiplicity, illustrative_gamma=args.gamma,
        domain_width_rlc=args.domain_width_rlc,
        local_grid_cells_per_axis=args.local_grid, local_stellar_radius_rlc=args.local_star_rlc,
        local_domain_width_rlc=args.local_domain_width_rlc)
    report = calculate_budget(inputs)
    root = Path(__file__).resolve().parents[2]
    sources = [Path(__file__), root/'third_party/entity-research/src/framework/parameters/parameters.cpp',
               root/'third_party/entity-research/src/kernels/pushers/sr.hpp']
    report['analysis_source_provenance'] = {
        str(path.relative_to(root)): file_hash(path) for path in sources if path.is_file()}
    if args.native_config.is_file():
        native = native_reference_from_config(args.native_config, inputs)
        physical_ratio = report['physical_reference_scales']['surface_cold_cyclotron_to_rotation_frequency_ratio']
        for species in native['species_reference_ratios']:
            reference = species['cold_cyclotron_to_target_rotation_reference']
            species['physical_reference_ratio_over_native_reference_ratio'] = (
                physical_ratio/reference if reference > 0 else None)
        report['recorded_native_configuration_reference'] = native
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    scales = report['physical_reference_scales']
    print(json.dumps({'output': str(args.output),
                      'R_LC_km': scales['light_cylinder_radius_km'],
                      'R_star_over_R_LC': scales['stellar_radius_over_rlc'],
                      'reference_LC_cold_skin_depth_m': scales['lc_cold_pair_plasma']['cold_skin_depth_m'],
                      'numerical_star_if_period_mapped_km': report['local_reduced_model_comparison']['numerical_radius_if_period_mapped_to_vela_km'],
                      'claim': 'Scale budget only; not a Vela reconstruction'}, indent=2))


if __name__ == '__main__':
    main()
