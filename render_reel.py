#!/usr/bin/env python3
"""Render the original Sideratlas Vela-informed scientific illustration.

Volume emissivity and the received pulse are distinct calculated quantities.
No observer-beaming brightness is invented for the all-direction volume.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time,subprocess
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.fft import rfft2,irfft2,next_fast_len
from PIL import Image,ImageDraw,ImageFont
from pulsar.snapshot_fields import rotate_z
from reel_export import encoder,metadata as encoding_metadata


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class PulsarReel:
    def __init__(self,scene,width=1080):
        if width<270 or width%2:raise ValueError('Even width at least270 required')
        self.width=width;self.height=round(width*16/9)//2*2;self.scale=width/1080
        with np.load(scene,allow_pickle=False) as d:
            self.points=d['points'].astype(float);self.energy=d['energy'];self.loge=d['log10_critical_energy_gev']
            self.tracks=d['track_points'].astype(float);self.offsets=d['track_offsets']
            self.edges=d['phase_edges'];self.pulse=d['pulse_photons_per_sr']
            self.meta=json.loads(str(d['metadata']))
        if np.any(~np.isfinite(self.points)) or np.any(self.energy<0):raise ValueError('Invalid scene')
        self.hero_y=0;self.hero_h=self.height
        self.hero_w=self.width;self.center=np.array([self.width/2,self.hero_h/2])
        self.center[1]=self.height*.53
        self.worldscale=480*self.scale
        self.roll=np.deg2rad(-28.)
        voxel_width=2*self.meta['voxel_extent_rlc']/(self.meta['voxel_grid']-1)
        # Project the trilinear voxel footprint, rather than treating each
        # voxel as a delta-function light. Its per-axis variance is h²/6.
        azimuths=self.meta['physical_protocol']['cap']['azimuths']
        cylindrical_second_moment=float(np.average(np.sum(self.points[:,:2]**2,axis=-1),weights=self.energy))
        # Record the finite cap-cell scale as a display-resolution diagnostic.
        # The final map retains its sampled emission structures; only the
        # voxel reconstruction footprint is applied to its rasterization.
        self.cap_display_sigma_world=2*np.pi/azimuths*np.sqrt(cylindrical_second_moment/6)
        self.display_sigma=float(np.sqrt(voxel_width**2/6*self.worldscale**2+(.35*self.scale)**2))
        self.voxel_width=voxel_width
        self.pad=max(16,int(np.ceil(3*voxel_width*self.worldscale)))
        self.fft_h=next_fast_len(self.hero_h+2*self.pad)
        self.fft_w=next_fast_len(self.hero_w+2*self.pad)
        self.fy=np.fft.fftfreq(self.fft_h).astype('f4')[:,None]
        self.fx=np.fft.rfftfreq(self.fft_w).astype('f4')[None,:]
        self.zeta=np.deg2rad(self.meta['observer_deg'])
        self.palette_stops=np.log10([.03,.15,.7,3.])
        self.colors=np.stack([np.interp(self.loge,self.palette_stops,c) for c in
            ([25,35,255,255],[75,204,177,245],[230,244,70,225])],axis=-1)/255
        # Convert display palette to linear light before additive integration.
        self.linear_colors=np.where(self.colors<=.04045,self.colors/12.92,((self.colors+.055)/1.055)**2.4)
        self.channels=self.energy[:,None]*self.linear_colors
        self.overlay=self.make_overlay()
        # A fixed exposure for every frame, measured over a whole turn.
        luminance=np.concatenate([np.max(self.volume(i/12),axis=-1).astype('f4').ravel() for i in range(12)])
        positive=luminance[luminance>0]
        self.exposure=float(np.quantile(positive,.997))*.5 if len(positive) else 1.
        self.max_curve=float(self.pulse.max())
        if self.max_curve<=0:raise ValueError('No calculated pulse')
        self.extent_check=self.check_bounds()

    def font(self,size,role='regular'):
        file={'regular':'lmroman10-regular.otf','bold':'lmroman12-bold.otf','math':'STIXTwoMath-Regular.otf'}[role]
        return ImageFont.truetype(str(ROOT/'fonts'/file),max(8,round(size*self.scale)))

    def project(self,points,phase):
        p=rotate_z(points,2*np.pi*(phase%1))
        x=p[:,1];y=-np.cos(self.zeta)*p[:,0]+np.sin(self.zeta)*p[:,2]
        xx=np.cos(self.roll)*x-np.sin(self.roll)*y
        yy=np.sin(self.roll)*x+np.cos(self.roll)*y
        return self.center+np.stack((xx,-yy),axis=-1)*self.worldscale

    def check_bounds(self):
        minimum=np.full(2,np.inf);maximum=np.full(2,-np.inf);fractions=[]
        for i in range(72):
            q=self.project(self.points,i/72);minimum=np.minimum(minimum,q.min(axis=0));maximum=np.maximum(maximum,q.max(axis=0))
            visible=(q[:,0]>=0)&(q[:,0]<self.width)&(q[:,1]>=0)&(q[:,1]<self.height)
            fractions.append(float(self.energy[visible].sum()/self.energy.sum()))
        return {'minimum_pixel':minimum.tolist(),'maximum_pixel':maximum.tolist(),
            'intentional_full_bleed_crop':True,'camera_roll_deg':float(np.degrees(self.roll)),
            'fraction_of_display_emission_centres_in_frame_range':[min(fractions),max(fractions)],'phases_checked':72}

    def volume(self,phase):
        q=self.project(self.points,phase)+self.pad;base=np.floor(q).astype(int);f=q-base
        out=np.zeros((self.fft_h*self.fft_w,3),'f4')
        for dx,dy in ((0,0),(0,1),(1,0),(1,1)):
            ij=base+np.array([dx,dy]);w=np.where(dx,f[:,0],1-f[:,0])*np.where(dy,f[:,1],1-f[:,1])
            valid=(ij[:,0]>=0)&(ij[:,0]<self.fft_w)&(ij[:,1]>=0)&(ij[:,1]<self.fft_h)
            ids=ij[valid,1]*self.fft_w+ij[valid,0]
            for c in range(3):out[:,c]+=np.bincount(ids,weights=w[valid]*self.channels[valid,c],minlength=len(out))
        out=out.reshape(self.fft_h,self.fft_w,3)
        # Exact Fourier projection of a normalized 3D trilinear voxel basis:
        # each Cartesian tent contributes sinc² of its projected frequency.
        # This removes visible voxel lattices without blurring cap structure.
        origin=self.project(np.zeros((1,3)),phase)[0]
        matrix=(self.project(np.eye(3),phase)-origin).T
        transfer=np.ones((self.fft_h,self.fft_w//2+1),'f4')
        for j in range(3):
            transfer*=np.sinc(self.voxel_width*(self.fx*matrix[0,j]+self.fy*matrix[1,j]))**2
        transfer*=np.exp(-2*np.pi**2*(.35*self.scale)**2*(self.fx**2+self.fy**2))
        projected=irfft2(rfft2(out,axes=(0,1))*transfer[...,None],s=(self.fft_h,self.fft_w),axes=(0,1))
        return np.maximum(projected[self.pad:self.pad+self.hero_h,self.pad:self.pad+self.hero_w],0)

    def make_overlay(self):
        im=Image.new('RGBA',(self.width,self.height));d=ImageDraw.Draw(im);s=self.scale
        white=(246,246,246);soft=(221,225,229)
        def text(x,y,v,size=28,color=white,role='regular',anchor='lt'):
            d.text((x*s,y*s),v,font=self.font(size,role),fill=color,anchor=anchor,
                stroke_width=max(1,round(1.5*s)),stroke_fill=(0,0,0,220))
        text(72,258,'@sideratlas',34)
        text(965,265,'NO.002',25,soft,anchor='rt')
        text(65,320,'THE PULSAR',106,role='bold')
        return im

    def frame(self,phase):
        p=phase%1;light=self.volume(p)/self.exposure
        # One fixed hue-preserving filmic display curve, no frame normalization.
        value=light/(1+np.max(light,axis=-1,keepdims=True))
        rgb=np.where(value<=.0031308,12.92*value,1.055*value**(1/2.4)-.055)
        hero=Image.fromarray(np.uint8(np.clip(rgb*255+.5,0,255)))
        im=Image.new('RGB',(self.width,self.height),'black');im.paste(hero,(0,self.hero_y));d=ImageDraw.Draw(im);s=self.scale
        im=im.convert('RGBA')
        d=ImageDraw.Draw(im)
        cx,cy=self.center
        r=max(2,round(2.6*s));d.ellipse((cx-r,cy-r,cx+r,cy+r),fill=(230,235,240))
        im=Image.alpha_composite(im,self.overlay)
        return im.convert('RGB')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('scene',type=Path)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--width',type=int,default=1080)
    ap.add_argument('--fps',type=int,default=30);ap.add_argument('--repeats',type=int,default=4)
    ap.add_argument('--stills-only',action='store_true');a=ap.parse_args()
    if a.repeats<1 or a.fps<1:raise ValueError('Positive repeats and frame rate required')
    a.output.mkdir(parents=True,exist_ok=True);started=time.monotonic();renderer=PulsarReel(a.scene,a.width)
    for i,phase in enumerate((0.,.125,.25,.375,.5,.625,.75,.875)):
        renderer.frame(phase).save(a.output/f'preview-{i:02d}.png')
    if not a.stills_only:
        loop=a.output/'sideratlas-pulsar-loop.mp4';process=encoder(loop,renderer.width,renderer.height,a.fps)
        try:
            for index in range(8*a.fps):
                process.stdin.write(renderer.frame(index/(8*a.fps)).tobytes())
                if index%a.fps==0:print(json.dumps({'frame':index,'frames':8*a.fps,'elapsed_s':time.monotonic()-started}),flush=True)
            process.stdin.close()
            if process.wait()!=0:raise RuntimeError('Video encoding failed')
        except BaseException:
            process.kill();raise
        concat=a.output/'repeat-inputs.txt';concat.write_text(''.join("file 'sideratlas-pulsar-loop.mp4'\n" for _ in range(a.repeats)))
        subprocess.run(['ffmpeg','-v','error','-y','-f','concat','-safe','0','-i',str(concat),'-c','copy','-movflags','+faststart',
            str(a.output/'sideratlas-002-pulsar-reel.mp4')],check=True)
    meta=dict(scope=renderer.meta['scope'],scientific_scene=renderer.meta,
        scene_path=str(a.scene.resolve()),scene_sha256=sha(a.scene),width=renderer.width,height=renderer.height,
        fps=a.fps,loop_seconds=8,repeats=a.repeats,duration_seconds=8*a.repeats,
        time_scale=8/renderer.meta['period_s'],camera='Fixed orthographic observer at declared viewing inclination',
        palette='Illustrative characteristic photon energies, not visible-light true color',
        palette_stops_gev=[.03,.15,.7,3.],
        palette_stops_srgb_255=[[25,75,230],[35,204,244],[255,177,70],[255,245,225]],
        exposure=renderer.exposure,exposure_selection='Half the99.7th percentile of nonzero linear-channel pixel peaks across12uniform phases; fixed for whole movie',
        tone_map='Hue-preserving x/(1+max_channel) in linear light, then sRGB encoding',
        display_pixel_footprint_equivalent_sigma=renderer.display_sigma,
        display_cap_footprint_sigma_rlc=renderer.cap_display_sigma_world,
        display_pixel_footprint='Fourier-slice projection of the exact trilinear3Dvoxel basis (product of3projected sinc² factors), with0.35-reference-pixel Gaussian raster footprint; finite padding before deliberate frame crop; no cap smoothing or plasma diffusion model',
        particle_paths='Not overlaid; visible brightness is only the calculated all-direction emission scene',
        field_of_view_check=renderer.extent_check,loop_endpoint_max_pixel_difference=int(np.max(abs(np.asarray(renderer.frame(0),int)-np.asarray(renderer.frame(1),int)))),
        overlay_labels=['@sideratlas','THE PULSAR','NO.002'],
        text_safe_region_1080p={'left':65,'right':965,'top':258,'bottom':400},
        export=encoding_metadata(),quantitative_vela_validation=False,elapsed_seconds=time.monotonic()-started,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'reel_export.py',ROOT/'prepare_reel_scene.py']})
    (a.output/'render-metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(a.output,flush=True)


if __name__=='__main__':main()
