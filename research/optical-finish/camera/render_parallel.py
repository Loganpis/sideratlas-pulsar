#!/usr/bin/env python3
"""Render the calculated optical scene with a memory-sharing worker pool.

Only image integration is parallelized. Each worker inherits the same physical
camera, atmosphere image, fixed display exposure and source event data. Source
energies, SI units, arrival intervals and the full-cycle energy audit are those
of FiniteVolumeCamera; the worker pool introduces no new emitting structures.
"""
import os

# Keep numerical libraries from starting a thread pool inside every process.
# These assignments must precede NumPy/SciPy imports, including indirect ones.
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'

import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from PIL import Image

from multiresolution_camera import MultiResolutionCamera
from atmosphere_camera import TabulatedAtmosphere, SurfaceComponent
from render_optical import ROOT, XYZ_TO_LINEAR_SRGB, encoder, encoder_metadata


_CAMERA = None
_FRAME_DIR = None


def log(event, **values):
    print(json.dumps({'event': event, **values}), flush=True)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def records(index):
    return {'index': int(index), 'crop': _CAMERA.crop_records[int(index)],
            'volume': _CAMERA.volume_records[int(index)]}


def exposure_job(index):
    started = time.monotonic()
    log('exposure_job_started', bin=int(index), worker=os.getpid())
    radiance = _CAMERA.radiance(index)
    peak = np.maximum(radiance @ XYZ_TO_LINEAR_SRGB.T, 0).max(axis=-1)
    samples = peak[peak > 0].astype(np.float32)
    result = records(index)
    result.update(samples=samples, seconds=time.monotonic() - started)
    log('exposure_job_finished', bin=int(index), seconds=result['seconds'],
        positive_pixels=int(samples.size), worker=os.getpid())
    return result


def frame_job(job):
    frame_index, bin_index = job
    started = time.monotonic()
    log('frame_job_started', frame=int(frame_index), bin=int(bin_index), worker=os.getpid())
    _CAMERA.frame(bin_index).save(_FRAME_DIR / f'frame-{frame_index:04d}.png',
                                compress_level=3)
    result = records(bin_index)
    result.update(frame=int(frame_index), seconds=time.monotonic() - started)
    log('frame_job_finished', frame=int(frame_index), bin=int(bin_index),
        seconds=result['seconds'], worker=os.getpid())
    return result


def merge_record(camera, result):
    camera.crop_records[result['index']] = result['crop']
    camera.volume_records[result['index']] = result['volume']


def indices_list(value):
    try:
        result = [int(x) for x in value.split(',')]
    except ValueError as exc:
        raise argparse.ArgumentTypeError('Use comma-separated integer detector bins') from exc
    if not result:
        raise argparse.ArgumentTypeError('At least one detector bin is required')
    return result


def main():
    global _CAMERA, _FRAME_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--surface', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--width', type=int, default=1080)
    parser.add_argument('--bins', type=int, default=240)
    parser.add_argument('--fps', type=int, default=30)
    parser.add_argument('--repeats', type=int, default=4)
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--motion-arc-pixels', type=float, default=2.)
    parser.add_argument('--plasma-scale', type=int, default=1,
                        help='Integer plasma-raster reduction; surface and labels retain native resolution')
    parser.add_argument('--exposure-reference', choices=['auto', 'surface'], default='auto',
                        help='Calibrate one display divisor from sampled scene pixels or the native stellar surface')
    parser.add_argument('--meters-per-pixel', type=float)
    parser.add_argument('--roll', type=float, default=-55.)
    parser.add_argument('--phase-offset', type=int)
    parser.add_argument('--previews-only', action='store_true')
    parser.add_argument('--preview-indices', type=indices_list,
                        help='Absolute detector bins for preview PNGs; the selected exposure reference stays fixed')
    parser.add_argument('--exposure-from', type=Path,
                        help='Reuse exposure from this wrapper’s matching render-metadata.json')
    args = parser.parse_args()
    if (args.width < 2 or args.width % 2 or args.bins < 1 or args.fps < 1
            or args.repeats < 1 or not 1 <= args.workers <= 3):
        parser.error('Use an even positive width, positive counts, and one to three workers')
    if 'fork' not in mp.get_all_start_methods():
        parser.error('A fork-capable platform is required for memory-sharing rendering')
    if args.preview_indices is not None and not args.previews_only:
        parser.error('--preview-indices is only valid with --previews-only')
    args.output.mkdir(parents=True, exist_ok=True)
    _FRAME_DIR = args.output / 'frames'
    _FRAME_DIR.mkdir(exist_ok=True)
    started = time.monotonic()
    log('camera_initializing', scene=str(args.scene), width=args.width, bins=args.bins)
    _CAMERA = camera = MultiResolutionCamera(
        args.scene, width=args.width, arrival_bins=args.bins,
        meters_per_pixel=args.meters_per_pixel, roll_deg=args.roll,
        motion_arc_pixels=args.motion_arc_pixels, plasma_scale=args.plasma_scale)
    log('camera_initialized', seconds=time.monotonic() - started,
        retained_events=len(camera.time), meters_per_pixel=camera.pixel_m)
    atmosphere = TabulatedAtmosphere(args.surface, allow_unvalidated_model=True)
    surface = SurfaceComponent(atmosphere, width=camera.width, height=camera.height,
                               meters_per_pixel=camera.pixel_m,
                               observer_direction=camera.observer,
                               period_s=camera.period, roll_deg=camera.roll)
    camera.surface_xyz = surface
    # The uniform surface is independent of arrival phase. Calculate it once
    # before fork, so all workers share its immutable cached image.
    log('surface_initializing')
    surface(0.)
    log('surface_initialized', seconds=time.monotonic() - started)
    identity = {'scene_sha256': sha(args.scene), 'surface_sha256': atmosphere.sha256,
                'width': camera.width, 'height': camera.height, 'arrival_bins': camera.bins,
                'meters_per_pixel': camera.pixel_m, 'roll_deg': camera.roll,
                'motion_quadrature_max_arc_pixels': camera.motion_arc,
                'motion_arc_pixel_grid': 'Plasma raster pixels',
                'plasma_scale': camera.plasma_scale,
                'plasma_center_fraction': list(camera.low_center_fraction),
                'exposure_reference': args.exposure_reference,
                'psf_sigma_pixels': camera.psf,
                'renderer_sha256': sha(__file__),
                'optical_camera_sha256': sha(Path(__file__).with_name('render_optical.py')),
                'finite_volume_camera_sha256': sha(Path(__file__).with_name('finite_volume_camera.py')),
                'geometry_compaction_sha256': sha(Path(__file__).with_name('compact_render_geometry.py')),
                'multiresolution_camera_sha256': sha(Path(__file__).with_name('multiresolution_camera.py')),
                'volume_raster_sha256': sha(Path(__file__).with_name('volume_raster.py')),
                'surface_adapter_sha256': sha(Path(__file__).with_name('atmosphere_camera.py')),
                'surface_camera_sha256': sha(Path(__file__).with_name('surface_camera.py')),
                'optical_color_sha256': sha(ROOT / 'src/pulsar/optical_color.py'),
                'cie_csv_sha256': sha(ROOT / 'data/color/CIE_xyz_1931_2deg.csv')}
    calibration_bins = (np.array([], dtype=int) if args.exposure_reference == 'surface' else
                        np.unique(np.r_[np.linspace(0, camera.bins, 24, endpoint=False).astype(int),
                                         np.argsort(camera.bin_xyz_energy[:, 1])[-8:]]))
    ctx = mp.get_context('fork')
    if args.exposure_from is not None:
        previous = json.loads(args.exposure_from.read_text())
        if previous.get('render_identity') != identity:
            raise ValueError('Exposure metadata does not match this scene, camera and renderer')
        camera.exposure = float(previous['exposure'])
        if not np.isfinite(camera.exposure) or camera.exposure <= 0:
            raise ValueError('Reused exposure must be finite and positive')
        log('fixed_exposure_reused', exposure=camera.exposure, source=str(args.exposure_from))
    elif args.exposure_reference == 'surface':
        surface_rgb = np.maximum(surface(0.) @ XYZ_TO_LINEAR_SRGB.T, 0)
        camera.exposure = .5 * float(surface_rgb.max())
        if not np.isfinite(camera.exposure) or camera.exposure <= 0:
            raise ValueError('No finite positive native stellar surface radiance for display exposure')
        del surface_rgb
        calibration_bins = np.array([], dtype=int)
        log('surface_fixed_exposure_selected', exposure=camera.exposure)
    else:
        log('exposure_calibration_started', bins=calibration_bins.tolist(), workers=args.workers)
        samples = []
        with ctx.Pool(args.workers) as pool:
            for result in pool.imap_unordered(exposure_job, calibration_bins.tolist(), chunksize=1):
                merge_record(camera, result)
                samples.append(result['samples'])
                log('exposure_progress', completed=len(samples), total=len(calibration_bins),
                    bin=result['index'], seconds=result['seconds'])
        values = np.concatenate(samples)
        if values.size == 0:
            raise ValueError('No optical photons available for fixed exposure')
        camera.exposure = float(np.quantile(values, .99)) * .35
        del samples, values
        if not np.isfinite(camera.exposure) or camera.exposure <= 0:
            raise ValueError('No finite positive fixed exposure')
        log('fixed_exposure_selected', exposure=camera.exposure)
    # Preserve the expensive global calibration even if a later render is
    # interrupted. --exposure-from accepts this receipt as well as a complete
    # render-metadata.json, after checking every camera/source identity field.
    exposure_receipt = {'render_identity': identity, 'exposure': camera.exposure,
                        'exposure_calibration_bins': calibration_bins.tolist(),
                        'scope': 'Fixed display exposure only; no physical luminosity fit'}
    (args.output / 'exposure-metadata.json').write_text(
        json.dumps(exposure_receipt, indent=2) + '\n')
    phase_offset = (int(np.argmax(camera.bin_xyz_energy[:, 1]))
                    if args.phase_offset is None else args.phase_offset) % camera.bins
    if args.previews_only:
        frame_bins = (args.preview_indices if args.preview_indices is not None
                      else [round(i * camera.bins / 8) + phase_offset for i in range(8)])
        jobs = [(i, int(index) % camera.bins) for i, index in enumerate(frame_bins)]
    else:
        jobs = [(i, (i + phase_offset) % camera.bins) for i in range(camera.bins)]
    log('frame_rendering_started', total=len(jobs), workers=args.workers,
        exposure=camera.exposure, phase_offset=phase_offset)
    # Recreate the pool after setting the global exposure. Newly forked workers
    # inherit that single value instead of calculating private per-frame gains.
    with (args.output / 'frame-ledger.jsonl').open('w') as ledger:
        ledger.write(json.dumps({'render_identity': identity, 'exposure': camera.exposure}) + '\n')
        ledger.flush()
        with ctx.Pool(args.workers) as pool:
            for completed, result in enumerate(pool.imap_unordered(frame_job, jobs, chunksize=1), 1):
                merge_record(camera, result)
                ledger.write(json.dumps(result) + '\n')
                ledger.flush()
                log('render_progress', completed=completed, total=len(jobs), frame=result['frame'],
                    bin=result['index'], seconds=result['seconds'], elapsed_seconds=time.monotonic() - started)
    if args.previews_only:
        for i, _ in jobs:
            shutil.copyfile(_FRAME_DIR / f'frame-{i:04d}.png', args.output / f'preview-{i:02d}.png')
    else:
        for i in range(8):
            frame_index = min(camera.bins - 1, round(i * camera.bins / 8))
            shutil.copyfile(_FRAME_DIR / f'frame-{frame_index:04d}.png', args.output / f'preview-{i:02d}.png')
    spatial = camera.volume_metadata()
    # This check addresses the complete source/arrival integral, including
    # energy falling outside the image crop. Cropping is a separate ledger.
    cycle_error = spatial.get('full_cycle_relative_energy_error')
    if cycle_error is not None and cycle_error > 1e-6:
        raise ValueError(f'Finite-support full-cycle energy audit failed: {cycle_error}')
    export = None
    if not args.previews_only:
        loop = args.output / 'optical-loop-8s.mp4'
        movie = args.output / 'sideratlas-no002-optical-32s.mp4'
        log('encoding_started', loop=str(loop), frames=camera.bins)
        process = encoder(loop, camera.width, camera.height, args.fps)
        try:
            for frame_index in range(camera.bins):
                with Image.open(_FRAME_DIR / f'frame-{frame_index:04d}.png') as image:
                    process.stdin.write(image.convert('RGB').tobytes())
                if frame_index % args.fps == 0:
                    log('encoding_progress', frames=frame_index, total=camera.bins)
        except BaseException:
            process.stdin.close()
            process.wait()
            raise
        process.stdin.close()
        if process.wait() != 0:
            raise RuntimeError('Optical loop encoding failed')
        ffmpeg = shutil.which('ffmpeg')
        if ffmpeg is None:
            raise RuntimeError('FFmpeg is required for repeat concatenation')
        subprocess.run([ffmpeg, '-hide_banner', '-loglevel', 'error', '-y',
                        '-stream_loop', str(args.repeats - 1), '-i', str(loop),
                        '-map', '0:v:0', '-c', 'copy', '-movflags', '+faststart', str(movie)], check=True)
        export = {**encoder_metadata(), 'loop_frames': camera.bins, 'repeats': args.repeats,
                  'frames': camera.bins * args.repeats,
                  'duration_seconds': camera.bins * args.repeats / args.fps,
                  'loop_sha256': sha(loop), 'movie_sha256': sha(movie), 'movie': movie.name}
        log('encoding_finished', movie=str(movie), bytes=movie.stat().st_size)
    metadata = {'scope': 'Calculated approximate optical visualization of a plausible rotation-powered pulsar; not a Vela observational reconstruction',
                'physical_model': 'Prescribed near-zone dipole, pair activation and physical harmonic radiation; optical emission is separate from the reduced PIC dynamics',
                'render_identity': identity, 'scene_sha256': identity['scene_sha256'],
                'width': camera.width, 'height': camera.height,
                'meters_per_pixel': camera.pixel_m, 'roll_deg': camera.roll,
                'exposure': camera.exposure, 'exposure_calibration_bins': calibration_bins.tolist(),
                'exposure_reference': args.exposure_reference,
                'exposure_definition': ('Half the maximum positive linear-sRGB component of the native calculated stellar surface; one divisor for all frames'
                                        if args.exposure_reference == 'surface' else
                                        '0.35 times the 99th percentile of positive sampled scene RGB peaks; one divisor for all frames'),
                'exposure_reused_from': str(args.exposure_from) if args.exposure_from else None,
                'source_period_s': camera.period, 'arrival_bins': camera.bins,
                'temporal_energy_relative_error': camera.temporal_energy_error,
                'psf_sigma_pixels': camera.psf,
                'spectral_color': 'CIE1931 XYZ from calculated optical harmonic lines; 683 applied once; fixed exposure to sRGB',
                'tone_map': 'Positive in-gamut RGB / (1+maximumRGB); negative out-of-gamut channels clipped',
                'surface_included': True, 'surface': surface.metadata(),
                'phase_offset_bins': phase_offset,
                'rendered_frame_detector_bins': [index for _, index in jobs],
                'spatial_camera': spatial, 'crop': camera.crop_records,
                'time_mapping': f'One physical rotation slowed to {camera.bins / args.fps:g} seconds; constructed loop repeated {args.repeats} times',
                'worker_processes': args.workers, 'export': export,
                'elapsed_seconds': time.monotonic() - started}
    (args.output / 'render-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    log('render_finished', metadata=str(args.output / 'render-metadata.json'),
        elapsed_seconds=metadata['elapsed_seconds'], full_cycle_relative_energy_error=cycle_error)


if __name__ == '__main__':
    main()
