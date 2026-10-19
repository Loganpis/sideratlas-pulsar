#!/usr/bin/env python3
"""Independent output/observational audit; no model fitting or renormalization."""
import argparse
import hashlib
import json
from pathlib import Path
from contextlib import contextmanager

import numpy as np
from scipy.constants import c, parsec, m_e
from load_camera_events import load_camera_events


@contextmanager
def _loaded_camera_archive(path):
    yield load_camera_events(path)


def audit(folder, distance_pc=290.):
    folder = Path(folder)
    path = folder/'physical-pair-events.npz'
    model_report = json.loads((folder/'physical-pair-report.json').read_text())
    actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual_hash != model_report['npz_sha256']:
        raise ValueError('Native event archive does not match its recorded SHA256')
    with _loaded_camera_archive(path) as d:
        required = ('positions_m', 'source_time_s', 'source_interval_s', 'number_weights',
                    'xyz_lm_sr', 'wavelength_edges_m',
                    'pattern_beta', 'observer_direction', 'period_s', 'phase_index')
        for key in required:
            if key not in d or np.any(~np.isfinite(d[key])):
                raise ValueError(f'Missing/nonfinite {key}')
        weight = d['number_weights']
        duration = d['source_interval_s']
        powers = d.get('power_bins_w_sr', d.get('event_band_power_w_sr'))
        if np.any(weight < 0) or np.any(duration <= 0) or np.any(powers < 0) or np.any(d['xyz_lm_sr'] < 0):
            raise ValueError('Nonphysical weights, times or powers')
        if np.any(np.linalg.norm(d['pattern_beta'], axis=1) >= 1):
            raise ValueError('Superluminal rotating pattern')
        n = d['observer_direction']
        period = float(d['period_s'])
        xyz = d['positions_m']
        arrival = np.remainder(d['source_time_s']-xyz@n/c, period)
        event_energy = (powers.sum(axis=1) if powers.ndim == 2 else powers)*weight*duration
        received, _ = np.histogram(arrival, np.linspace(0, period, model_report['model']['phase_nodes']+1), weights=event_energy)
        xyz_energy = (d['xyz_lm_sr']*(weight*duration)[:, None]).sum(axis=0)
        spectral_energy = ((powers*(weight*duration)[:, None]).sum(axis=0)
                           if powers.ndim == 2 else d['phase_spectral_power_w_sr'].sum(axis=0)*float(d['source_interval_scalar_s']))
        edges = d['wavelength_edges_m']
        frequency_width = c/edges[:-1]-c/edges[1:]
        frequency = c/np.sqrt(edges[:-1]*edges[1:])
        # Flat-space distant-observer geometric dilution. No extinction is
        # applied: compare to intrinsic/dereddened source estimates only.
        distance = distance_pc*parsec
        flux_density_si = spectral_energy/period/frequency_width/distance**2
        flux_jy = flux_density_si/1e-26
        wavelength_nm = np.sqrt(edges[:-1]*edges[1:])*1e9
        finite = flux_density_si > 0
        slope = np.polyfit(np.log(frequency[finite]), np.log(flux_density_si[finite]), 1)[0]
        f550_jy = np.interp(550., wavelength_nm, flux_jy)
        approximate_mag_ab = -2.5*np.log10(f550_jy/3631.)
        # Mignani2007 gives V~23.6; adopt A_V=.16 only for a rough intrinsic
        # normalization check. Johnson V/AB bandpass offsets are not modeled.
        reference_v, assumed_av = 23.6, .16
        approximate_reference_dereddened_jy = 3631*10**(-.4*(reference_v-assumed_av))
        observed_alpha = -.04
        nu0 = c/(550e-9)
        nu_low, nu_high = c/(700e-9), c/(400e-9)
        reference_band_flux = (approximate_reference_dereddened_jy*1e-26*nu0**observed_alpha
                               *(nu_high**(1-observed_alpha)-nu_low**(1-observed_alpha))/(1-observed_alpha))
        predicted_band_flux = float(event_energy.sum()/period/distance**2)
        result = dict(
            scientific_status='Prediction audit for an explicitly assumed population; NOT a fit, calibrated Vela reconstruction, detector simulation or validation',
            event_archive_sha256=actual_hash,
            event_count=len(weight), all_saved_values_finite=True,
            all_weights_and_spectral_powers_nonnegative=True,
            independent_source_band_energy_j_sr=float(event_energy.sum()),
            independent_retarded_band_energy_j_sr=float(received.sum()),
            relative_source_report_error=float(abs(event_energy.sum()-model_report['total_source_frame_band_energy_j_sr'])/event_energy.sum()),
            relative_retarded_energy_error=float(abs(received.sum()-event_energy.sum())/event_energy.sum()),
            source_xyz_luminous_energy_lm_s_sr=xyz_energy.tolist(),
            chosen_distance_pc=distance_pc,
            predicted_unextincted_integrated_400_700nm_flux_w_m2=predicted_band_flux,
            approximate_intrinsic_reference_400_700nm_flux_w_m2=float(reference_band_flux),
            approximate_model_to_reference_400_700nm_band_flux_ratio=predicted_band_flux/reference_band_flux,
            wavelength_bin_centres_nm=wavelength_nm.tolist(),
            predicted_unextincted_bin_mean_fnu_jy=flux_jy.tolist(),
            predicted_approximate_fnu_at_550nm_jy=float(f550_jy),
            predicted_AB_like_magnitude_at_550nm_not_Johnson_V=float(approximate_mag_ab),
            predicted_alpha_from_binned_band_fnu=float(-slope),
            spectral_convention='F_nu proportional nu^(-alpha); unweighted logarithmic slope of bin-mean model flux over400-700nm; not a data likelihood',
            observed_reference=dict(
                paper='Mignani, Zharikov & Caraveo2007, https://arxiv.org/abs/0707.2036, doi:10.1051/0004-6361:20077774',
                reported_alpha=-.04, reported_alpha_uncertainty=.04,
                reference_spectral_interval_nm=[400, 800], approximate_V_magnitude=reference_v,
                adopted_AV_for_rough_flux_comparison=assumed_av,
                approximate_dereddened_550nm_reference_fnu_jy=approximate_reference_dereddened_jy),
            approximate_model_to_reference_550nm_flux_ratio=float(f550_jy/approximate_reference_dereddened_jy),
            magnitude_comparison_limit='Reference V and monochromatic AB-like model quantity have different bandpass conventions; ratio is an approximate order-of-magnitude check only',
            no_model_parameters_fit_to_optical_data=True,
            no_display_exposure_used_in_this_audit=True,
            conservative_residual_perpendicular_energy_bounds_j=[
                .5*column['final_gamma']*column['final_transverse_momentum']**2*m_e*c*c
                for column in model_report['column_diagnostics']],
            residual_energy_bound_scope='For continued outward motion in decreasing B without heating: -dgamma/dt=a*p_perp^2 and -dp_perp/dt>=a*p_perp*(1+p_perp^2)/gamma imply remaining radiation <= gamma_final*p_perp_final^2*m_e*c^2/2. The producer uniform-field residual is a slightly smaller reference estimate.',
            optical_validation=False,
            assumptions_with_potential_large_effect=[
                'Monoenergetic gamma300/pitch0.1 activation instead of a calculated cascade/radio-heated distribution',
                'Prescribed uniform activation-radius probability and cap supply, not measured spatial emissivity',
                'Dipole/corotation geometry instead of a force-free or PIC-derived physical field',
                'Optically thin propagation; cyclotron/synchrotron absorption and plasma response not calculated',
                'Numerical quadrature and virtual imaging need separate convergence assessment'])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder', type=Path)
    p.add_argument('--distance-pc', type=float, default=290.)
    args = p.parse_args()
    result = audit(args.folder, args.distance_pc)
    out = args.folder/'independent-optical-audit.json'
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in (
        'predicted_unextincted_integrated_400_700nm_flux_w_m2',
        'predicted_approximate_fnu_at_550nm_jy',
        'predicted_AB_like_magnitude_at_550nm_not_Johnson_V',
        'predicted_alpha_from_binned_band_fnu',
        'approximate_model_to_reference_550nm_flux_ratio',
        'relative_retarded_energy_error')}, indent=2))


if __name__ == '__main__':
    main()
