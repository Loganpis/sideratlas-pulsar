#!/usr/bin/env python3
"""Compute a disclosed physical pair-radiation component, never a PIC rescale."""
import argparse
import json
from pathlib import Path

from physical_synchrotron import PairModel, generate_events


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phases', type=int, default=128)
    parser.add_argument('--azimuths', type=int, default=64)
    parser.add_argument('--cooling', type=int, default=32)
    parser.add_argument('--radii', type=int, default=3)
    parser.add_argument('--annuli', type=int, default=1)
    parser.add_argument('--wavelength-bins', type=int, default=64)
    parser.add_argument('--no-aberration', action='store_true')
    parser.add_argument('--compact-camera', action='store_true', help='Preserve exact line-color calculation; store reconstructible geometry and measured-error float32 XYZ without per-event spectra')
    args = parser.parse_args()
    model = PairModel(phase_nodes=args.phases, azimuth_nodes=args.azimuths,
                      cooling_nodes=args.cooling, radial_nodes=args.radii,
                      annulus_nodes=args.annuli, wavelength_nodes=args.wavelength_bins,
                      aberration=not args.no_aberration)
    args.output.mkdir(parents=True, exist_ok=True)
    report = generate_events(model, args.output/'physical-pair-events.npz',
                             args.output/'physical-pair-report.json', compact_camera=args.compact_camera)
    print(json.dumps({k: report[k] for k in (
        'scientific_status', 'source_event_count', 'injected_total_leptons_per_local_second',
        'local_synchrotron_loss_luminosity_w', 'total_source_frame_band_energy_j_sr',
        'relative_retarded_energy_histogram_error', 'phases_above_one_percent_peak',
        'max_corotation_beta')}, indent=2))


if __name__ == '__main__':
    main()
