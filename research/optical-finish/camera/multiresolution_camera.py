"""Shared physical camera with a coarse plasma raster and native stellar image.

Only raster sampling is changed. The finite-support source/arrival calculation
retains its original SI energies; integer-factor bilinear interpolation is
audited for integrated radiance conservation. The star is evaluated on the
native physical pixel grid and receives the same final optical sampling PSF.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom

from finite_volume_camera import FiniteVolumeCamera
from render_optical import OpticalCamera


def upsample_radiance(low, scale):
    """Interpolate cell-centred radiance without altering integrated XYZ flux."""
    low = np.asarray(low, dtype=float)
    if (low.ndim != 3 or low.shape[-1] != 3 or not np.isfinite(low).all()
            or not isinstance(scale, int) or scale < 1):
        raise ValueError('Finite XYZ image and positive integer resampling scale required')
    high = (low.copy() if scale == 1 else
            zoom(low, (scale, scale, 1), order=1, mode='nearest',
                 prefilter=False, grid_mode=True))
    # Physical pixel area shrinks by scale²; never renormalize an image to
    # disguise interpolation loss. Integer grid-mode bilinear interpolation
    # with nearest edge continuation preserves this integral to roundoff.
    before = low.sum(axis=(0, 1))
    after = high.sum(axis=(0, 1)) / scale ** 2
    absolute = np.abs(after - before)
    relative = np.divide(absolute, np.abs(before), out=np.zeros(3), where=before != 0)
    if np.any(relative > 1e-10):
        raise ValueError('Bilinear plasma resampling failed integrated-radiance conservation')
    return high, relative


class MultiResolutionCamera:
    """FiniteVolumeCamera proxy exposing native output dimensions and scale."""

    def __init__(self, scene, *, width=1080, plasma_scale=1, meters_per_pixel=None,
                 psf_sigma_pixels=.6, motion_arc_pixels=2., **kwargs):
        height = round(width * 16 / 9) // 2 * 2
        if (not isinstance(plasma_scale, int) or plasma_scale < 1
                or width % plasma_scale or height % plasma_scale):
            raise ValueError('Plasma scale must divide both native image dimensions')
        if not np.isfinite(psf_sigma_pixels) or psf_sigma_pixels < 0:
            raise ValueError('Finite nonnegative final sampling PSF required')
        self.width, self.height = width, height
        self.plasma_scale = plasma_scale
        self.psf = float(psf_sigma_pixels)
        low_width, low_height = width // plasma_scale, height // plasma_scale
        if round(low_width * 16 / 9) // 2 * 2 != low_height:
            raise ValueError('Low raster must retain the exact native aspect ratio')
        # ndimage.zoom(grid_mode=True) maps low pixel j to native coordinate
        # scale*j+(scale-1)/2. Shift its physical camera centre accordingly.
        low_cx = (width * .5 - (plasma_scale - 1) / 2) / plasma_scale
        low_cy = (height * .53 - (plasma_scale - 1) / 2) / plasma_scale
        self.low_center_fraction = (low_cx / low_width, low_cy / low_height)
        self.inner = FiniteVolumeCamera(
            scene, width=low_width,
            meters_per_pixel=(None if meters_per_pixel is None else meters_per_pixel * plasma_scale),
            psf_sigma_pixels=0, motion_arc_pixels=motion_arc_pixels,
            center_fraction=self.low_center_fraction, **kwargs)
        self.pixel_m = self.inner.pixel_m / plasma_scale
        self.surface_xyz = None
        self.exposure = None
        self.overlay = OpticalCamera.make_overlay(self)
        # Parent and child record merging uses the same dict objects as the
        # inner finite-volume ledger; no separate energy ledger is invented.
        self.volume_records = self.inner.volume_records
        self.crop_records = self.inner.crop_records

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, 'inner'), name)

    choose_exposure = OpticalCamera.choose_exposure
    frame = OpticalCamera.frame
    frame_from_radiance = OpticalCamera.frame_from_radiance

    def radiance(self, bin_index):
        index = int(bin_index) % self.bins
        low = self.inner.radiance(index)
        frame, relative = upsample_radiance(low, self.plasma_scale)
        if frame.shape != (self.height, self.width, 3):
            raise ValueError('Resampled plasma grid does not match native detector dimensions')
        record = self.volume_records[index]
        record['plasma_resampling_xyz_relative_flux_error'] = relative.tolist()
        if self.surface_xyz is not None:
            frame += self.surface_xyz((index + .5) * self.period / self.bins)
        if self.psf > 0:
            before = frame.sum(axis=(0, 1))
            frame = gaussian_filter(frame, (self.psf, self.psf, 0), mode='constant', truncate=4)
            record['psf_edge_lost_xyz_radiance_sum'] = (before - frame.sum(axis=(0, 1))).tolist()
        return frame

    def volume_metadata(self):
        out = self.inner.volume_metadata()
        errors = [x['plasma_resampling_xyz_relative_flux_error']
                  for x in self.volume_records.values()
                  if 'plasma_resampling_xyz_relative_flux_error' in x]
        out['multiresolution_camera'] = {
            'native_width': self.width, 'native_height': self.height,
            'native_meters_per_pixel': self.pixel_m,
            'plasma_width': self.inner.width, 'plasma_height': self.inner.height,
            'plasma_meters_per_pixel': self.inner.pixel_m,
            'integer_plasma_scale': self.plasma_scale,
            'plasma_center_fraction': list(self.low_center_fraction),
            'native_center_fraction': [.5, .53],
            'resampling': 'Integer-factor, cell-centred bilinear interpolation; grid_mode=True, nearest edge continuation',
            'maximum_resampling_xyz_relative_flux_error': np.max(errors, axis=0).tolist() if errors else None,
            'motion_quadrature_max_arc_plasma_pixels': self.inner.motion_arc,
            'motion_quadrature_max_arc_native_pixels': self.inner.motion_arc * self.plasma_scale,
            'surface_raster': 'Native output grid at the same physical camera scale',
            'plasma_psf_sigma_pixels': 0,
            'shared_final_psf_sigma_native_pixels': self.psf,
            'limitation': 'Sub-plasma-pixel morphology is not resolved; interpolation changes raster sampling, not source radiation or display exposure'}
        return out
