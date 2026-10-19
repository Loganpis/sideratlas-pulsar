"""Calculated atmosphere table adapter and resolved spherical surface camera.

A single normal-field atmosphere column is applied uniformly over the star.
This is an explicit approximation, not a dipolar field atmosphere library.
The table must pass its numerical equilibrium gate; physical validation is
separate and never inferred from that gate. No blackbody extrapolation occurs.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from surface_camera import Star, SurfaceRayMap, surface_raymap, spectral_surface_image
from pulsar.optical_color import spectrum_to_xyz


class TabulatedAtmosphere:
    def __init__(self, path, *, allow_unvalidated_model=False):
        self.path = Path(path)
        self.sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()
        with np.load(path, allow_pickle=False) as data:
            wavelength = data['wavelength_m'].astype(float)*1e9
            mu = data['mu'].astype(float)
            intensity = data['i_lambda_si'].astype(float)/1e9
            numerical = bool(data['numerical_equilibrium_pass']) if 'numerical_equilibrium_pass' in data else False
            physical = bool(data['scientific_validation_pass'])
            units = str(data['units'])
        if units != 'W m^-2 sr^-1 m^-1':
            raise ValueError('Atmosphere intensity units must be explicit SI per metre')
        if (intensity.shape != (len(wavelength),len(mu)) or len(mu)<2
                or not np.isfinite(intensity).all() or np.any(intensity<0)
                or not np.isfinite(wavelength).all() or np.any(wavelength<=0)
                or np.any(np.diff(wavelength)<=0) or np.any(np.diff(mu)<=0)
                or mu[0] != 0 or mu[-1] != 1):
            raise ValueError('Invalid atmosphere table; full [0,1] angular support required')
        if not numerical:
            raise ValueError('Atmosphere has not passed numerical equilibrium')
        if not physical and not allow_unvalidated_model:
            raise ValueError('Explicit approximate-model opt-in required for unvalidated atmosphere')
        self.numerical_equilibrium_pass, self.scientific_validation_pass = numerical, physical
        self.wave, self.mu = wavelength, mu
        self.interpolator = RegularGridInterpolator((wavelength,mu),intensity,bounds_error=True)

    def __call__(self, *, wavelength_emitted_nm, cos_emission_angle, **unused):
        wave=np.asarray(wavelength_emitted_nm,float)
        mu=np.broadcast_to(np.asarray(cos_emission_angle,float)[:,None],wave.shape)
        # Bounds errors must propagate: no X-ray extension or grazing-angle fill.
        return self.interpolator(np.stack((wave,mu),axis=-1))


class SurfaceComponent:
    def __init__(self, atmosphere, *, width,height,meters_per_pixel,observer_direction,
                 period_s,roll_deg=-55.,subpixels=3):
        if not isinstance(subpixels,int) or subpixels<1:
            raise ValueError('Positive subpixel quadrature required')
        self.atmosphere=atmosphere
        self.width,self.height,self.pixel_m=width,height,meters_per_pixel
        self.period=period_s
        self.star=Star()
        self.subpixels=subpixels
        self.wavelength=np.arange(400.,701.,2.)
        # Evaluate only the resolved stellar bounding box, on the same physical
        # pixel grid as the magnetosphere. The star is never enlarged separately.
        radius=self.star.apparent_radius_m/meters_per_pixel
        cx,cy=width*.5,height*.53
        self.x0=max(0,int(np.floor(cx-radius-1)))
        self.x1=min(width,int(np.ceil(cx+radius+1)))
        self.y0=max(0,int(np.floor(cy-radius-1)))
        self.y1=min(height,int(np.ceil(cy+radius+1)))
        offset=(np.arange(subpixels)+.5)/subpixels-.5
        xx=np.arange(self.x0,self.x1)[:,None]+offset
        yy=np.arange(self.y0,self.y1)[:,None]+offset
        # Order h,sub_y,w,sub_x so averaging the two subpixel axes is exact.
        x=np.broadcast_to((xx-cx)[None,None,:,:],(len(yy),subpixels,len(xx),subpixels))*meters_per_pixel
        y=np.broadcast_to(-(yy-cy)[:,:,None,None],x.shape)*meters_per_pixel
        angle=np.deg2rad(roll_deg)
        # Undo the image roll before building the observer's physical ray map.
        raw_x=np.cos(angle)*x+np.sin(angle)*y
        raw_y=-np.sin(angle)*x+np.cos(angle)*y
        self.rays=surface_raymap(raw_x,raw_y,star=self.star,
                                  observer_direction=observer_direction,up_hint=(0,0,1))
        self.cached=None

    def __call__(self, arrival_time_s):
        # One uniform, azimuth-averaged normal-field column is axisymmetric;
        # image is time independent even with local rotational Doppler shifts.
        if self.cached is None:
            # Batches retain the same rays and spectra while avoiding a full
            # subpixel-by-wavelength cube for a large, resolved stellar image.
            rays=self.rays
            hit=np.flatnonzero(rays.hit.ravel())
            flat_xyz=np.zeros((rays.hit.size,3))
            for start in range(0,len(hit),4096):
                take=hit[start:start+4096]
                subset=SurfaceRayMap(rays.star,np.ones(len(take),dtype=bool),
                    rays.normal_static.reshape(-1,3)[take],rays.photon_static.reshape(-1,3)[take],
                    rays.cos_alpha.ravel()[take],rays.relative_delay_s.ravel()[take],
                    rays.arrival_offset_s.ravel()[take],rays.observer_direction,rays.reference_radius_m)
                spectral=spectral_surface_image(subset,arrival_time_s,self.wavelength,
                    self.atmosphere,omega_rad_s=2*np.pi/self.period)
                flat_xyz[take]=spectrum_to_xyz(self.wavelength,spectral)
            xyz=flat_xyz.reshape(rays.hit.shape+(3,)).mean(axis=(1,3))
            self.cached=np.zeros((self.height,self.width,3))
            self.cached[self.y0:self.y1,self.x0:self.x1]=xyz
        return self.cached

    def metadata(self):
        return {'surface_table':str(self.atmosphere.path),'surface_table_sha256':self.atmosphere.sha256,
                'numerical_equilibrium_pass':self.atmosphere.numerical_equilibrium_pass,
                'scientific_validation_pass':self.atmosphere.scientific_validation_pass,
                'mass_solar':self.star.mass_solar,'radius_m':self.star.radius_m,
                'frequency_ratio_gravity':self.star.frequency_ratio_g,
                'apparent_radius_pixels':self.star.apparent_radius_m/self.pixel_m,
                'surface_subpixels_per_axis':self.subpixels,'observer_wavelength_step_nm':2,
                'surface_assumption':'Uniform normal-B atmosphere column; no surface temperature or field map',
                'geodesics':'Direct outward Schwarzschild integrals with local SR Doppler and aberration',
                'phase_scope':'Uniform surface has no temperature pattern; common GR-vs-flat outer propagation offset is not an aligned pulse prediction'}
