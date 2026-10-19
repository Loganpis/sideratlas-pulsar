#!/usr/bin/env python3
"""Render calculated observer-directed optical photons at a shared physical scale.

This renderer cannot validate particle supply, field dynamics or an atmosphere.
No field-guide emitters, diagnostic species colors or decorative stellar sphere
are added. A surface is included only when an explicit transfer table is given.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'research/optical-finish/particle-emission'))
from pulsar.optical_color import XYZ_TO_LINEAR_SRGB
from reel_export import encoder, metadata as encoder_metadata
from surface_camera import rotate_vectors
from atmosphere_camera import TabulatedAtmosphere, SurfaceComponent
from load_camera_events import load_camera_events
from event_camera import (rotating_pattern_arrival_edges, periodic_interval_weights,
                          inverse_rotating_arrival, project_m_to_pixels, splat_xyz_samples)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def weighted_quantile(value, weights, fraction):
    order = np.argsort(value)
    cdf = np.cumsum(weights[order])
    if cdf[-1] <= 0:
        raise ValueError('Scene contains no observer-directed optical light')
    return np.interp(fraction*cdf[-1], cdf, value[order])


class OpticalCamera:
    def __init__(self, scene, *, width=1080, arrival_bins=240, meters_per_pixel=None,
                 roll_deg=-55., psf_sigma_pixels=.6, fill_scale=1.1):
        self.scene = Path(scene)
        self.width, self.height = width, round(width*16/9)//2*2
        self.bins = arrival_bins
        self.roll = roll_deg
        self.psf = float(psf_sigma_pixels)
        self.surface_xyz = None
        with np.load(scene, allow_pickle=False) as saved:
            compact = bool(saved['compact_camera']) if 'compact_camera' in saved else False
        data = load_camera_events(scene) if compact else np.load(scene,allow_pickle=False)
        try:
            self.position = data['positions_m'].astype(float)
            self.time = data['source_time_s'].astype(float)
            self.dt = data['source_interval_s'].astype(float)
            self.number = data['number_weights'].astype(float)
            self.power = data['xyz_lm_sr'].astype(float)
            self.observer = data['observer_direction'].astype(float)
            self.period = float(data['period_s'])
            self.cell_index = data['cell_index'].astype(np.int32)
            self.reference_vertices = data['cell_footprint_vertices_reference_m'].astype(float)
        finally:
            if not compact:
                data.close()
        if (self.power.shape != self.position.shape or self.power.shape[-1] != 3
                or np.any(self.power < 0) or not np.isfinite(self.power).all()
                or np.any(self.number < 0) or not np.isfinite(self.number).all()
                or not np.isfinite(self.period) or self.period <= 0):
            raise ValueError('Invalid physical optical event data')
        self.omega = 2*np.pi/self.period
        left, right = rotating_pattern_arrival_edges(self.time, self.dt, self.position, self.observer, self.omega)
        ids, bins, _, arrival_left, arrival_right = periodic_interval_weights(left, right, self.period, self.bins, return_edges=True)
        source_left = inverse_rotating_arrival(arrival_left, self.time[ids], self.position[ids], self.observer, self.omega)
        source_right = inverse_rotating_arrival(arrival_right, self.time[ids], self.position[ids], self.observer, self.omega)
        self.record_order = np.argsort(bins, kind='stable')
        self.event_ids = ids[self.record_order]
        self.source_left, self.source_right = source_left[self.record_order], source_right[self.record_order]
        self.offset = np.r_[0, np.cumsum(np.bincount(bins, minlength=self.bins))]
        # Energy is integrated in SOURCE time. Arrival inversion supplies the
        # exact overlap duration for the prescribed pattern, not a second boost.
        deposited = np.bincount(ids, weights=source_right-source_left, minlength=len(self.time))
        energy = self.power*self.number[:, None]*self.dt[:, None]
        recovered = self.power*self.number[:, None]*deposited[:, None]
        self.temporal_energy_error = float(np.linalg.norm(recovered.sum(axis=0)-energy.sum(axis=0))/np.linalg.norm(energy.sum(axis=0)))
        if self.temporal_energy_error > 1e-9:
            raise ValueError('Optical event arrival integration lost energy')
        if meters_per_pixel is None:
            # One camera crop over a whole rotation, weighted by calculated Y.
            q = project_m_to_pixels(self.position, self.observer, meters_per_pixel=1,
                                    width=1, height=1, roll_deg=self.roll, center_fraction=(0, 0))
            weights = energy[:, 1]
            sx = weighted_quantile(np.abs(q[:, 0]), weights, .995)/(self.width*.5)
            sy = weighted_quantile(np.abs(q[:, 1]), weights, .995)/(self.height*.47)
            meters_per_pixel = max(sx, sy)/fill_scale
        if not np.isfinite(meters_per_pixel) or meters_per_pixel <= 0:
            raise ValueError('Invalid common physical pixel scale')
        self.pixel_m = float(meters_per_pixel)
        self.bin_xyz_energy = np.zeros((self.bins, 3))
        np.add.at(self.bin_xyz_energy, bins, self.power[ids]*self.number[ids, None]*(source_right-source_left)[:, None])
        self.crop_records = {}
        self.overlay = self.make_overlay()
        self.exposure = None

    def make_overlay(self):
        s = self.width/1080
        im = Image.new('RGBA', (self.width, self.height))
        draw = ImageDraw.Draw(im)
        def label(x, y, text, size, font='lmroman10-regular.otf', anchor='lt'):
            f = ImageFont.truetype(str(ROOT/'fonts'/font), round(size*s))
            draw.text((x*s, y*s), text, font=f, fill=(246,246,246), anchor=anchor,
                      stroke_width=max(1, round(1.5*s)), stroke_fill=(0,0,0,230))
        label(72,258,'@sideratlas',34)
        label(965,265,'NO.002',25,anchor='rt')
        label(65,320,'THE PULSAR',106,'lmroman12-bold.otf')
        return im

    def radiance(self, bin_index):
        index = int(bin_index)%self.bins
        sl = slice(self.offset[index], self.offset[index+1])
        ids = self.event_ids[sl]
        t0, t1 = self.source_left[sl], self.source_right[sl]
        if len(ids):
            # Integrate moving pattern cells during each detector exposure.
            # The upper bound on arc length keeps quadrature movement <0.6px.
            radius_xy = np.linalg.norm(self.position[ids, :2], axis=-1)
            pieces = np.maximum(1, np.ceil(self.omega*radius_xy*(t1-t0)/(self.pixel_m*.6)).astype(int))
            if np.any(pieces > 8192):
                raise ValueError('Exposure quadrature too coarse; refine source phase sampling')
            row = np.repeat(np.arange(len(ids)), pieces)
            begin = np.repeat(np.cumsum(pieces)-pieces, pieces)
            local = (np.arange(len(row))-begin+.5)/pieces[row]
            source_t = t0[row]+local*(t1[row]-t0[row])
            event = ids[row]
            position = rotate_vectors(self.position[event], [0,0,1], self.omega*(source_t-self.time[event]))
            xyz = (self.power[event]*self.number[event,None]
                   *((t1-t0)/pieces)[row,None])
            pixel = project_m_to_pixels(position, self.observer, meters_per_pixel=self.pixel_m,
                                        width=self.width, height=self.height, roll_deg=self.roll)
            image, lost = splat_xyz_samples(pixel, xyz, self.width, self.height)
        else:
            image = np.zeros((self.height,self.width,3)); lost = np.zeros(3)
        total = self.bin_xyz_energy[index]
        self.crop_records[index] = {'source_Y_energy':float(total[1]),'cropped_Y_energy':float(lost[1]),
                                    'retained_fraction':float(1-lost[1]/total[1]) if total[1]>0 else None}
        image /= self.pixel_m**2*(self.period/self.bins)
        if self.surface_xyz is not None:
            image += self.surface_xyz(index*self.period/self.bins)
        if self.psf > 0:
            image = gaussian_filter(image, (self.psf,self.psf,0), mode='constant', truncate=4)
        return image

    def choose_exposure(self):
        bins = np.unique(np.r_[np.linspace(0,self.bins,24,endpoint=False).astype(int),
                               np.argsort(self.bin_xyz_energy[:,1])[-8:]])
        samples = []
        for index in bins:
            rgb = np.maximum(self.radiance(index)@XYZ_TO_LINEAR_SRGB.T,0)
            peak = rgb.max(axis=-1)
            samples.append(peak[peak>0].astype('f4'))
        values = np.concatenate(samples)
        self.exposure = float(np.quantile(values,.99))*.35
        if not np.isfinite(self.exposure) or self.exposure<=0:
            raise ValueError('No optical photons available for fixed exposure')

    def frame(self, index):
        if self.exposure is None:
            self.choose_exposure()
        return self.frame_from_radiance(self.radiance(index))

    def frame_from_radiance(self, radiance):
        if self.exposure is None:
            raise ValueError('Set one fixed exposure before displaying radiance')
        rgb = np.maximum(np.asarray(radiance)@XYZ_TO_LINEAR_SRGB.T,0)/self.exposure
        linear = rgb/(1+rgb.max(axis=-1,keepdims=True))
        srgb = np.where(linear<=.0031308,12.92*linear,1.055*linear**(1/2.4)-.055)
        im = Image.fromarray(np.uint8(np.clip(srgb*255+.5,0,255))).convert('RGBA')
        return Image.alpha_composite(im,self.overlay).convert('RGB')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scene',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--width',type=int,default=1080)
    parser.add_argument('--bins',type=int,default=240)
    parser.add_argument('--meters-per-pixel',type=float)
    parser.add_argument('--roll',type=float,default=-55.)
    parser.add_argument('--surface',type=Path)
    parser.add_argument('--allow-approximate-model',action='store_true')
    parser.add_argument('--fps',type=int,default=30)
    parser.add_argument('--repeats',type=int,default=4)
    parser.add_argument('--phase-offset',type=int)
    parser.add_argument('--point-samples',action='store_true',help='Diagnostic point quadrature instead of finite supports')
    parser.add_argument('--motion-arc-pixels',type=float,default=2.)
    parser.add_argument('--previews-only',action='store_true')
    parser.add_argument('--allow-magnetosphere-only-preview',action='store_true')
    args=parser.parse_args()
    if args.surface is None and not (args.allow_magnetosphere_only_preview and args.previews_only):
        parser.error('A numerically accepted atmosphere table is required for combined export')
    if args.fps<1 or args.repeats<1:
        parser.error('Positive frame rate and repeat count required')
    args.output.mkdir(parents=True,exist_ok=True)
    cls=OpticalCamera
    extra={}
    if not args.point_samples:
        from finite_volume_camera import FiniteVolumeCamera
        cls=FiniteVolumeCamera
        extra={'motion_arc_pixels':args.motion_arc_pixels}
    camera=cls(args.scene,width=args.width,arrival_bins=args.bins,
               meters_per_pixel=args.meters_per_pixel,roll_deg=args.roll,**extra)
    if args.surface is not None:
        atmosphere=TabulatedAtmosphere(args.surface,allow_unvalidated_model=args.allow_approximate_model)
        surface=SurfaceComponent(atmosphere,width=camera.width,height=camera.height,
                                 meters_per_pixel=camera.pixel_m,observer_direction=camera.observer,
                                 period_s=camera.period,roll_deg=camera.roll)
        camera.surface_xyz=surface
    started=time.monotonic()
    camera.choose_exposure()
    phase_offset=int(np.argmax(camera.bin_xyz_energy[:,1])) if args.phase_offset is None else args.phase_offset
    for i in range(8):
        camera.frame(round(i*camera.bins/8)+phase_offset).save(args.output/f'preview-{i:02d}.png')
    export=None
    if not args.previews_only:
        movie=args.output/'sideratlas-no002-optical-32s.mp4'
        loop=args.output/'optical-loop-8s.mp4'
        process=encoder(loop,camera.width,camera.height,args.fps)
        try:
            for index in range(camera.bins):
                im=camera.frame(index+phase_offset)
                process.stdin.write(im.tobytes())
                if index%30==0:
                    print(json.dumps({'rendered_frames':index,'total_frames':camera.bins,'seconds':time.monotonic()-started}),flush=True)
        finally:
            process.stdin.close()
        if process.wait()!=0:
            raise RuntimeError('Optical loop encoding failed')
        subprocess.run(['ffmpeg','-hide_banner','-loglevel','error','-y','-stream_loop',str(args.repeats-1),
                        '-i',str(loop),'-map','0:v:0','-c','copy','-movflags','+faststart',str(movie)],check=True)
        export={**encoder_metadata(),'loop_frames':camera.bins,'repeats':args.repeats,
                'frames':camera.bins*args.repeats,'duration_seconds':camera.bins*args.repeats/args.fps,
                'loop_sha256':sha(loop),'movie_sha256':sha(movie),'movie':movie.name}
    metadata={'scope':('Calculated approximate Vela-informed optical model; not an observational reconstruction'
                       if args.surface else 'Physical magnetosphere component preview; stellar atmosphere not included'),
              'physical_model':'Prescribed near-zone dipole, pair activation and physical harmonic radiation; not PIC-derived Vela emission',
              'scene_sha256':sha(args.scene),'width':camera.width,'height':camera.height,
              'meters_per_pixel':camera.pixel_m,'roll_deg':camera.roll,'exposure':camera.exposure,
              'source_period_s':camera.period,'arrival_bins':camera.bins,
              'temporal_energy_relative_error':camera.temporal_energy_error,
              'motion_quadrature_maximum_arc_pixels':.6,'psf_sigma_pixels':camera.psf,
              'spectral_color':'CIE1931 XYZ from actual optical harmonic lines,683 applied once; fixed exposure to sRGB',
              'tone_map':'Positive in-gamut RGB / (1+maximumRGB); negative out-of-gamut channels clipped',
              'surface_included':args.surface is not None,
              'surface':surface.metadata() if args.surface else None,
              'phase_offset_bins':phase_offset,
              'spatial_camera':camera.volume_metadata() if hasattr(camera,'volume_metadata') else {'mode':'Point quadrature diagnostic'},
              'time_mapping':'One physical rotation slowed to 8 seconds at 240 bins / 30 fps; exact constructed loop repeated four times',
              'export':export,'crop':camera.crop_records,'elapsed_seconds':time.monotonic()-started}
    (args.output/'render-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps({k:metadata[k] for k in ('scope','meters_per_pixel','exposure','temporal_energy_relative_error')},indent=2))


if __name__=='__main__':
    main()
