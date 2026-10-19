#!/usr/bin/env python3
"""Close-up pulsar geometry with a calculated surface and illustrative guides.

The dipole integral curves are mathematical field geometry. Soft polar volumes
are emission-direction envelopes, never particle tracks or optical intensities.
Existing particle calculations are preserved separately, not recolored here.
"""
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing
import subprocess
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'research/optical-finish/camera')]
from atmosphere_camera import TabulatedAtmosphere, SurfaceComponent
from pulsar.optical_color import XYZ_TO_LINEAR_SRGB
from reel_export import encoder, metadata as encoding_metadata

PERIOD = .08937108859
RADIUS_M = 13000.
LIGHT_CYLINDER_M = 299792458.*PERIOD/(2*np.pi)
ALPHA = np.deg2rad(75.)


def rotate_z(p, angle):
    p = np.asarray(p)
    c, s = np.cos(angle), np.sin(angle)
    return np.stack((c*p[..., 0]-s*p[..., 1], s*p[..., 0]+c*p[..., 1], p[..., 2]), -1)


def magnetic_to_spin(p):
    """Rotate magnetic z onto the adopted oblique magnetic moment."""
    c, s = np.cos(ALPHA), np.sin(ALPHA)
    return np.asarray(p) @ np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])


def dipole_loops(shells=(2., 3.2, 4.7, 5.6), azimuths=10, samples=300):
    """Exact static dipole curves r=L sin²theta, expressed in stellar radii.

    Fade/crop near the surface avoids suggesting the Euclidean guides were
    ray-traced through the atmosphere's Schwarzschild image.
    """
    lines = []
    for j, shell in enumerate(shells):
        theta0 = np.arcsin(np.sqrt(1.3/shell))
        theta = np.linspace(theta0, np.pi-theta0, samples)
        r = shell*np.sin(theta)**2
        for phi in np.arange(azimuths)*2*np.pi/azimuths+(j%2)*np.pi/azimuths:
            q = np.stack((r*np.sin(theta)*np.cos(phi),
                          r*np.sin(theta)*np.sin(phi), r*np.cos(theta)), -1)
            lines.append(magnetic_to_spin(q))
    return lines


def polar_guides(samples=200):
    """Local segments of the nominal last-open dipole surface, not global topology."""
    r = np.linspace(1.3, 11., samples)
    theta = np.arcsin(np.sqrt(r/(LIGHT_CYLINDER_M/RADIUS_M)))
    lines = []
    for sign in (-1, 1):
        for phi in np.arange(6)*np.pi/3:
            q = np.stack((r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),
                          sign*r*np.cos(theta)), -1)
            lines.append(magnetic_to_spin(q))
    return lines


def visible(p, radius):
    transverse = np.sum(p[..., :2]**2, axis=-1)
    return (transverse >= radius*radius) | (p[..., 2] >= np.sqrt(np.maximum(radius*radius-transverse, 0)))


def pieces(mask):
    edge = np.diff(np.r_[False, mask, False].astype('i1'))
    return zip(np.flatnonzero(edge == 1), np.flatnonzero(edge == -1))


def smoothstep(a, b, x):
    u = np.clip((x-a)/(b-a), 0, 1)
    return u*u*(3-2*u)


class PulsarScene:
    def __init__(self, width=1080, surface_cache=None, guides=True, beams=True):
        if width < 180 or width%2:
            raise ValueError('Use an even width of at least 180 pixels')
        self.width, self.height = width, round(width*16/9)//2*2
        self.s = width/1080.
        self.scale = 105*self.s
        self.center = np.array([width*.5, self.height*.55])
        self.guides, self.beams = guides, beams
        self.lines = dipole_loops()
        self.polar = polar_guides()
        self.axis = magnetic_to_spin([0., 0., 1.])
        zeta, roll = np.deg2rad([65., -32.])
        right = np.array([0., 1., 0.])
        up = np.array([-np.cos(zeta), 0., np.sin(zeta)])
        eye = np.array([np.sin(zeta), 0., np.cos(zeta)])
        self.basis = np.stack((np.cos(roll)*right-np.sin(roll)*up,
                               np.sin(roll)*right+np.cos(roll)*up, eye))
        self.apparent_radius = 1/np.sqrt(1-2*6.67430e-11*1.4*1.988409870698051e30/(RADIUS_M*299792458.**2))
        # Dipole tangent angle at a declared reference height; this is a diagram
        # convention, not an inferred emission height or a radiation solution.
        theta = np.arcsin(np.sqrt(6*RADIUS_M/LIGHT_CYLINDER_M))
        self.beam_angle = theta+np.arctan(.5*np.tan(theta))
        self.overlay = self.labels()
        self.star, self.surface_metadata = self.surface(surface_cache)
        # Fixed diagram exposure: only the stellar layer has physical radiance.
        self.star = self.star/max(float(self.star.max()), 1e-30)*4.5
        self.star_halo = gaussian_filter(np.max(self.star, -1), 8*self.s)
        self.bw = max(180, width//3)
        self.bh = round(self.bw*self.height/width)
        xx = (np.arange(self.bw)+.5)*width/self.bw
        yy = (np.arange(self.bh)+.5)*self.height/self.bh
        self.rx = ((xx-self.center[0])/self.scale)[None, :].astype('f4')
        self.ry = ((self.center[1]-yy)/self.scale)[:, None].astype('f4')
        self.r2 = self.rx**2+self.ry**2
        self.disk = self.r2 < self.apparent_radius**2
        self.front = np.sqrt(np.maximum(self.apparent_radius**2-self.r2, 0))

    def surface(self, cache):
        if cache and Path(cache).exists():
            with np.load(cache, allow_pickle=False) as d:
                if tuple(d['rgb'].shape) != (self.height, self.width, 3):
                    raise ValueError('Cached surface dimensions differ')
                return d['rgb'].copy(), json.loads(str(d['metadata']))
        table = ROOT/'scene/optical/surface-intensity-si.npz'
        atmosphere = TabulatedAtmosphere(table, allow_unvalidated_model=True)
        # SurfaceComponent uses center .53H. Translate its bounding box onto
        # this composition; no independent stellar enlargement is performed.
        component = SurfaceComponent(atmosphere, width=self.width, height=self.height,
            meters_per_pixel=RADIUS_M/self.scale,
            observer_direction=(np.sin(np.deg2rad(65)), 0, np.cos(np.deg2rad(65))),
            period_s=PERIOD, roll_deg=-32., subpixels=2)
        xyz = component(0.)
        rgb = np.maximum(xyz @ XYZ_TO_LINEAR_SRGB.T, 0).astype('f4')
        dy = round(self.center[1]-self.height*.53)
        rgb = np.roll(rgb, dy, axis=0)
        if dy>0:
            rgb[:dy] = 0
        meta = component.metadata()
        meta['composition_translation_pixels_y'] = dy
        if cache:
            Path(cache).parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(cache, rgb=rgb, metadata=json.dumps(meta))
        return rgb, meta

    def labels(self):
        image = Image.new('RGBA', (self.width, self.height))
        draw = ImageDraw.Draw(image)
        def label(x, y, text, size, bold=False, anchor='lt'):
            font = ImageFont.truetype(str(ROOT/'fonts'/('lmroman12-bold.otf' if bold else 'lmroman10-regular.otf')), round(size*self.s))
            draw.text((x*self.s, y*self.s), text, font=font, fill=(247, 248, 251),
                      stroke_width=max(1, round(1.5*self.s)), stroke_fill=(0, 0, 0, 220), anchor=anchor)
        label(72, 258, '@sideratlas', 34)
        label(965, 265, 'NO.002', 25, anchor='rt')
        label(65, 320, 'THE PULSAR', 106, True)
        return image

    def camera(self, points, phase):
        phase = phase%1
        azimuth = np.deg2rad(12)*np.sin(2*np.pi*phase)
        return rotate_z(points, 2*np.pi*(phase+.25)-azimuth) @ self.basis.T

    def pixels(self, p):
        return self.center+p[..., :2]*np.array([1., -1.])*self.scale

    def field_light(self, phase):
        image = Image.new('L', (self.width*2, self.height*2))
        draw = ImageDraw.Draw(image)
        if self.guides:
            for i, line in enumerate(self.lines+self.polar):
                cam = self.camera(line, phase)
                pixels = self.pixels(cam)*2
                mask = visible(cam, self.apparent_radius*1.003)
                # Piecewise shading gives a continuous, depth-readable curve.
                for a, b in pieces(mask):
                    for u in range(a, b-1, 5):
                        v = min(u+6, b)
                        depth = cam[u:v, 2].mean()
                        r = np.linalg.norm(line[u:v], axis=-1).mean()
                        fade = smoothstep(1.3, 1.7, r)
                        if i >= len(self.lines):
                            fade *= 1-smoothstep(7, 11, r)
                        value = round((130+80*np.tanh(depth/2))*fade)
                        draw.line([tuple(q) for q in pixels[u:v]], fill=value,
                                  width=max(1, round(2.3*self.s)), joint='curve')
        mask = np.asarray(image.resize((self.width, self.height), Image.Resampling.LANCZOS), dtype='f4')/255
        glow = gaussian_filter(mask, max(.5, 2.5*self.s))
        haze = gaussian_filter(mask, max(.8, 10*self.s))
        return (mask[..., None]*np.array([.28, 1.0, 2.2])+
                glow[..., None]*np.array([.18, .6, 1.4])*2.0+
                haze[..., None]*np.array([.08, .25, .6])*1.2).astype('f4')

    def beam_light(self, phase):
        broad = np.zeros((self.bh, self.bw), 'f4')
        core = broad.copy()
        if self.beams:
            m = self.camera(self.axis, phase)
            dz = .32
            for depth in np.arange(-14+dz/2, 14, dz):
                along = abs(m[0]*self.rx+m[1]*self.ry+m[2]*depth)
                perp2 = np.maximum(self.r2+depth**2-along**2, 0)
                sigma = .055+np.tan(self.beam_angle)*np.maximum(along-1, 0)
                envelope = smoothstep(1., 1.7, along)*(1-smoothstep(10, 14, along))
                envelope *= np.exp(-along/10)/(1+.06*along)
                base = envelope*((~self.disk)|(depth >= self.front))*dz
                broad += np.exp(-.5*perp2/sigma**2)*base
                core += np.exp(-.5*perp2/(sigma*.3)**2)*base
        # A coherent blue-white palette. No diagnostic rainbow or species color.
        light = (broad[..., None]*np.array([.13, .34, .72])+
                 core[..., None]*np.array([.60, 1.0, 1.4]))
        channels = [np.asarray(Image.fromarray(light[..., i].astype('f4')).resize(
            (self.width, self.height), Image.Resampling.BICUBIC)) for i in range(3)]
        return np.maximum(np.stack(channels, -1), 0)

    def frame(self, phase):
        phase = phase%1
        light = self.star+self.field_light(phase)+self.beam_light(phase)
        # Local halo is an explicit screen bloom, not magnetospheric radiation.
        light += self.star_halo[..., None]*np.array([.06, .13, .25])
        # A fixed per-channel shoulder lets the intensely exposed stellar
        # surface approach white. This display map is not a spectral result.
        value = light/(1+light)
        rgb = np.where(value <= .0031308, 12.92*value, 1.055*value**(1/2.4)-.055)
        im = Image.fromarray(np.uint8(np.clip(rgb*255+.5, 0, 255))).convert('RGBA')
        return Image.alpha_composite(im, self.overlay).convert('RGB')

    def metadata(self):
        return dict(scope='Conceptual 3D close-up, not a telescope image or global PIC reconstruction',
            parameters=dict(period_s=PERIOD, mass_solar=1.4, radius_m=RADIUS_M,
                polar_field_gauss=3e12, magnetic_obliquity_deg=75, temperature_infinity_K=700000),
            physical_scale=dict(pixels_per_stellar_radius_at_1080=105,
                field_apex_radii_Rstar=[2, 3.2, 4.7, 5.6], polar_guides_max_radius_Rstar=11,
                light_cylinder_m=LIGHT_CYLINDER_M, star_separately_enlarged=False),
            geometry=dict(field='Analytic near-zone dipole: r=L sin^2(theta), rigidly rotated',
                polar_guides='Local nominal last-open dipole surface, r=R_LC sin^2(theta); no global openness solution',
                field_projection='Euclidean orthographic; not lensed field lines',
                line_sampling='Representative, not a magnetic flux density or plasma density map'),
            beams=dict(role='Illustrative emission-direction envelopes; not optical intensity, photon paths or material jets',
                reference_radius_Rstar=6, angular_width_rad=float(self.beam_angle),
                width_rule='theta+atan(tan(theta)/2), theta=asin(sqrt(6R*/R_LC))',
                reference_height_assumed=True, gaussian_profile_assumed=True,
                graphical_axial_fade_Rstar=[10, 14], two_poles_symmetric=True),
            surface=self.surface_metadata,
            display=dict(stellar_exposure_normalization=4.5, guides_palette='Illustrative blue-white',
                common_absolute_radiometry=False, bloom='Screen-space presentation effect',
                tone_map='Fixed per-channel x/(1+x), then sRGB; compresses color toward white at high exposure'),
            camera=dict(inclination_deg=65, roll_deg=-32, azimuth_excursion_deg=12,
                center_height_fraction=.55, projection='Orthographic'),
            time='One 89.371ms stellar rotation displayed in 8 seconds; smooth periodic virtual-camera azimuth; repeated four times',
            calculated_optical_particle_archive_displayed=False,
            existing_physics_recalculated=False, labels=['@sideratlas', 'NO.002', 'THE PULSAR'])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


_WORKER_SCENE = None


def render_worker(index):
    return index, _WORKER_SCENE.frame(index/240).tobytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--width', type=int, default=1080)
    parser.add_argument('--preview-only', action='store_true')
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    scene = PulsarScene(args.width, args.output/f'surface-cache-{args.width}.npz')
    phases = [0, .125, .25, .375, .5, .625, .75, .875]
    previews = []
    for i, p in enumerate(phases):
        frame = scene.frame(p)
        frame.save(args.output/f'preview-{i:02d}.png')
        previews.append(frame.resize((270, 480), Image.Resampling.LANCZOS))
    sheet = Image.new('RGB', (270*4, 480*2))
    for i, im in enumerate(previews):
        sheet.paste(im, ((i%4)*270, (i//4)*480))
    sheet.save(args.output/'contact-sheet.jpg', quality=95)
    if args.preview_only:
        print(args.output/'contact-sheet.jpg', flush=True)
        return
    if args.workers not in (1, 2, 3):
        raise ValueError('Use one to three rendering workers')
    loop = args.output/'sideratlas-no002-pulsar-3d-loop.mp4'
    movie = args.output/'sideratlas-no002-pulsar-3d-32s.mp4'
    process = encoder(loop, scene.width, scene.height, 30)
    global _WORKER_SCENE
    _WORKER_SCENE = scene
    # A small bounded queue avoids retaining hundreds of full-size RGB frames.
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('fork')) as pool:
        pending = {i: pool.submit(render_worker, i) for i in range(args.workers*2)}
        for i in range(240):
            _, pixels = pending.pop(i).result()
            process.stdin.write(pixels)
            nxt = i+args.workers*2
            if nxt < 240:
                pending[nxt] = pool.submit(render_worker, nxt)
            if i%30 == 0:
                print(f'Frame {i}/240; elapsed {time.monotonic()-started:.1f}s', flush=True)
    process.stdin.close()
    if process.wait():
        raise RuntimeError('Video encoder failed')
    subprocess.run(['ffmpeg', '-hide_banner', '-loglevel', 'error', '-y', '-stream_loop', '3',
                    '-i', str(loop), '-map', '0:v:0', '-c', 'copy', '-movflags', '+faststart', str(movie)], check=True)
    meta = scene.metadata()
    meta.update(width=scene.width, height=scene.height, fps=30, frames=960, duration_seconds=32,
                export=encoding_metadata(), elapsed_seconds=round(time.monotonic()-started, 2),
                source_sha256=sha(__file__), video_sha256=sha(movie))
    (args.output/'render-metadata.json').write_text(json.dumps(meta, indent=2)+'\n')
    print(movie, flush=True)


if __name__ == '__main__':
    main()
