#!/usr/bin/env python3
"""External 3D diagram of actual gamma-ray rays, with a delayed receiver.

The camera sees diagram glyphs, not scattered optical light. Orthographic
projection preserves geometry; only display size, sampling and color change.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time,subprocess
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.ndimage import gaussian_filter
from PIL import Image,ImageDraw,ImageFont
from pulsar.photon_scene import TAU,emission_event,propagate,observer_event,wavefront_crossing,received_flux
from reel_export import encoder,metadata as encoding_metadata

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class PhotonReel:
    def __init__(self,scene,width=1080,births=48):
        self.width=width;self.height=round(width*16/9)//2*2;self.scale=width/1080
        with np.load(scene) as d:
            self.r0=d['position0'];self.n0=d['direction0'];self.fractions=d['band_fraction']
            self.weights=d['packet_weight'];self.edges=d['phase_edges'];self.pulse=d['pulse_bands'].sum(axis=0)
            self.meta=json.loads(str(d['metadata']))
        self.births=births;self.worldscale=91*self.scale
        self.center=np.array([width*.5,self.height*.53]);self.max_age=13.
        self.distance=5.2;self.frames=240;self.shutter=1/self.frames
        # Fixed camera, distinct from the fixed distant observer's direction.
        az,polar,roll=np.deg2rad([-70.,55.,-30.])
        right=np.array([-np.sin(az),np.cos(az),0.])
        up=np.array([-np.cos(polar)*np.cos(az),-np.cos(polar)*np.sin(az),np.sin(polar)])
        self.camera=np.stack((np.cos(roll)*right-np.sin(roll)*up,
            np.sin(roll)*right+np.cos(roll)*up))
        self.zeta=np.deg2rad(self.meta['observer_deg'])
        self.observer=np.array([np.sin(self.zeta),0.,np.cos(self.zeta)])
        self.receiver_pixel=self.project(self.distance*self.observer)
        self.palette=np.array([[65,134,255],[93,221,255],[255,231,173],[255,137,58]])/255
        linear=np.where(self.palette<=.04045,self.palette/12.92,((self.palette+.055)/1.055)**2.4)
        self.source_colors=self.fractions@linear
        # Emission cadence is staggered among packets. Every point is a real
        # direction from the model, rotated only at its birth time.
        offset=(np.arange(len(self.r0))*.6180339887498949)%1
        birth=(np.arange(births)[:,None]+offset[None,:])*TAU/births
        self.birth=birth.ravel()
        self.birth_r,self.birth_n=emission_event(np.broadcast_to(self.r0,(births,)+self.r0.shape).reshape(-1,3),
            np.broadcast_to(self.n0,(births,)+self.n0.shape).reshape(-1,3),self.birth)
        self.colors=np.tile(self.source_colors,(births,1)).astype('f4')
        self.vector_ids=np.linspace(0,len(self.birth)-1,160,dtype=int)
        self.overlay=self.make_overlay()
        # Fixed common exposure from a turn; never independently rescale peaks.
        samples=[]
        for phase in np.arange(12)/12:
            image=self.radiation(phase)
            values=image.max(axis=-1);samples.append(values[values>0])
        self.exposure=float(np.quantile(np.concatenate(samples),.995))*2.0
        self.flux_peak=float(received_flux(self.edges,self.pulse,np.arange(240)/240,self.distance,self.shutter).max())

    def project(self,points):
        q=np.asarray(points)@self.camera.T
        return self.center+q*np.array([self.worldscale,-self.worldscale])

    def font(self,size,bold=False):
        return ImageFont.truetype(str(ROOT/'fonts'/('lmroman12-bold.otf' if bold else 'lmroman10-regular.otf')),max(8,round(size*self.scale)))

    def make_overlay(self):
        im=Image.new('RGBA',(self.width,self.height));d=ImageDraw.Draw(im);s=self.scale
        for x,y,label,size,bold,anchor in [(72,258,'@sideratlas',34,False,'lt'),(965,265,'NO.002',25,False,'rt'),(65,320,'THE PULSAR',106,True,'lt')]:
            d.text((x*s,y*s),label,font=self.font(size,bold),fill=(246,246,246),anchor=anchor,
                stroke_width=max(1,round(1.5*s)),stroke_fill=(0,0,0,235))
        return im

    def packet_positions(self,phase):
        age=(TAU*(phase%1)-self.birth)%TAU
        positions=[];directions=[];colors=[]
        for previous_turn in range(3):
            a=age+previous_turn*TAU;valid=a<self.max_age
            positions.append(propagate(self.birth_r[valid],self.birth_n[valid],a[valid]))
            directions.append(self.birth_n[valid]);colors.append(self.colors[valid])
        return np.concatenate(positions),np.concatenate(directions),np.concatenate(colors)

    def radiation(self,phase):
        points,directions,colors=self.packet_positions(phase)
        q=self.project(points)
        # Integrate a one-frame physical shutter along each straight ray.
        # A compact pixel footprint is a glyph, not a gamma-ray beam opening.
        dq=(directions@self.camera.T)*np.array([self.worldscale,-self.worldscale])*TAU*self.shutter
        out=np.zeros((self.height*self.width,3),'f4')
        for offset in (-.25,.25):
            p=q+offset*dq;ij=np.floor(p+.5).astype(int)
            valid=(ij[:,0]>=0)&(ij[:,0]<self.width)&(ij[:,1]>=0)&(ij[:,1]<self.height)
            ids=ij[valid,1]*self.width+ij[valid,0]
            for channel in range(3):
                out[:,channel]+=np.bincount(ids,weights=colors[valid,channel]*.5,minlength=len(out))
        image=out.reshape(self.height,self.width,3)
        return gaussian_filter(image,sigma=(max(.6,2.8*self.scale),max(.6,2.8*self.scale),0),mode='constant')

    def frame(self,phase):
        phase=phase%1
        light=self.radiation(phase)/self.exposure
        v=light/(1+light.max(axis=-1,keepdims=True))
        rgb=np.where(v<=.0031308,12.92*v,1.055*v**(1/2.4)-.055)
        im=Image.fromarray(np.uint8(np.clip(rgb*255+.5,0,255))).convert('RGBA')
        # Sparse velocity glyphs distinguish photon flow from a glowing gas.
        # Their anchors and axes are actual packet positions and directions;
        # their drawn lengths are annotations, not a radiation opening angle.
        vector_layer=Image.new('RGBA',im.size);vd=ImageDraw.Draw(vector_layer)
        ids=self.vector_ids;age=(TAU*phase-self.birth[ids])%TAU
        age+=TAU*(np.arange(len(ids))%2)
        valid=age<self.max_age
        points=propagate(self.birth_r[ids][valid],self.birth_n[ids][valid],age[valid])
        ends=self.project(points+.12*self.birth_n[ids][valid])
        starts=self.project(points-.12*self.birth_n[ids][valid])
        for start,end in zip(starts,ends):
            if not (0<=end[0]<self.width and 0<=end[1]<self.height):continue
            delta=end-start;length=np.linalg.norm(delta)
            if length<3*self.scale:continue
            unit=delta/length;side=np.array([-unit[1],unit[0]])
            vd.line((*start,*end),fill=(210,226,242,85),width=max(1,round(self.scale)))
            vd.line([tuple(end-unit*4*self.scale+side*2*self.scale),tuple(end),tuple(end-unit*4*self.scale-side*2*self.scale)],
                fill=(210,226,242,100),width=max(1,round(self.scale)))
        im=Image.alpha_composite(im,vector_layer)
        d=ImageDraw.Draw(im);s=self.scale;cx,cy=self.center
        # Central marker only, deliberately much smaller than the emission region.
        r=max(1.,2.1*s);d.ellipse((cx-r,cy-r,cx+r,cy+r),fill=(245,248,255))
        # Fixed observer bearing, not a field line or an imposed emission cone.
        receiver=self.receiver_pixel;delta=receiver-self.center;unit=delta/np.linalg.norm(delta)
        for fraction in np.arange(.13,.92,.055):
            p=self.center+fraction*delta;p2=p+unit*4*s
            d.line((*p,*p2),fill=(118,128,140,120),width=max(1,round(s)))
        flux=float(received_flux(self.edges,self.pulse,phase,self.distance,self.shutter)[0])/self.flux_peak
        # Telescope/crosshair schematic and its calculated delayed pulse.
        x,y=receiver;radius=(8+12*flux)*s;white=round(125+130*flux)
        d.ellipse((x-radius,y-radius,x+radius,y+radius),outline=(white,white,white),width=max(1,round(2*s)))
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            d.line((x+dx*11*s,y+dy*11*s,x+dx*18*s,y+dy*18*s),fill=(220,225,230),width=max(1,round(s)))
        if flux>0:
            r=(2+4*flux)*s;d.ellipse((x-r,y-r,x+r,y+r),fill=(white,white,white))
        # Small scientific trace under the receiver. No new text labels.
        tx=x-86*s;ty=y+48*s;tw=172*s;th=38*s
        phases=np.arange(240)/240
        curve=received_flux(self.edges,self.pulse,phases,self.distance,self.shutter)/self.flux_peak
        d.line([(tx+i/239*tw,ty-float(f)*th) for i,f in enumerate(curve)],fill=(160,169,182),width=max(1,round(s)))
        xx=tx+phase*tw;yy=ty-flux*th;rr=2.5*s
        d.ellipse((xx-rr,yy-rr,xx+rr,yy+rr),fill=(255,255,255))
        return Image.alpha_composite(im,self.overlay).convert('RGB')

    def verify(self):
        te,r,n,phase=observer_event(self.r0,self.n0)
        crossing,arrival=wavefront_crossing(r,n,te,self.distance)
        residual=((arrival-self.distance)/TAU-phase+.5)%1-.5
        straight=propagate(self.birth_r[:1000],self.birth_n[:1000],2.)-propagate(self.birth_r[:1000],self.birth_n[:1000],1.)
        checks=dict(unit_direction_max_error=float(abs(np.linalg.norm(self.birth_n,axis=-1)-1).max()),
            straight_ray_c1_max_error=float(abs(straight-self.birth_n[:1000]).max()),
            observer_wavefront_phase_max_error_cycles=float(abs(residual).max()),
            observer_wavefront_plane_max_error_rlc=float(abs(np.sum(crossing*n,axis=-1)-self.distance).max()),
            packet_birth_position_max_radius_rlc=float(np.linalg.norm(self.birth_r,axis=-1).max()),
            loop_endpoint_position_max_error_rlc=float(abs(self.packet_positions(0)[0]-self.packet_positions(1)[0]).max()))
        passed=all(checks[k]<1e-10 for k in checks if 'radius' not in k)
        if not passed:raise ValueError(checks)
        return dict(passed=passed,checks=checks,
            observer_scope='Asymptotic angular observer, original 2deg polar pixel centered at zeta=65deg and azimuth=0; symbol position is a scaled bearing',
            receiver_equation='t_cross=t_emit+D-n.r_emit; displayed pulse uses Phi=(t_cross-D)/(2pi)',
            finite_packets='Weighted source midpoint display quadrature; receiver instead retains original full cell integration')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scene',type=Path,default=ROOT/'outputs/sideratlas-pulsar-photons/scene/photons.npz')
    p.add_argument('--output',type=Path,default=ROOT/'outputs/sideratlas-pulsar-photons/media')
    p.add_argument('--width',type=int,default=1080);p.add_argument('--births',type=int,default=48)
    p.add_argument('--stills-only',action='store_true');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    renderer=PhotonReel(a.scene,a.width,a.births)
    audit=renderer.verify();(a.output/'photon-direction-timing-verification.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)
    for i in range(8):renderer.frame(i/8).save(a.output/f'preview-{i:02d}.png')
    endpoint=int(abs(np.asarray(renderer.frame(0),int)-np.asarray(renderer.frame(1),int)).max())
    if endpoint:raise ValueError('Nonperiodic displayed loop')
    if not a.stills_only:
        loop=a.output/'sideratlas-pulsar-loop.mp4';process=encoder(loop,renderer.width,renderer.height,30)
        try:
            for i in range(240):
                process.stdin.write(renderer.frame(i/240).tobytes())
                if i%30==0:print(json.dumps(dict(frame=i,frames=240,elapsed_seconds=time.monotonic()-start)),flush=True)
            process.stdin.close()
            if process.wait():raise RuntimeError('Encoding failed')
        except BaseException:process.kill();raise
        concat=a.output/'repeat-inputs.txt';concat.write_text("file 'sideratlas-pulsar-loop.mp4'\n"*4)
        subprocess.run(['ffmpeg','-v','error','-y','-f','concat','-safe','0','-i',str(concat),'-c','copy','-movflags','+faststart',str(a.output/'sideratlas-002-pulsar-reel.mp4')],check=True)
    meta=dict(scientific_scene=renderer.meta,width=renderer.width,height=renderer.height,fps=30,loop_seconds=8,repeats=4,duration_seconds=32,
        scene_sha256=sha(a.scene),camera=dict(projection='Fixed external orthographic',azimuth_deg=-70,polar_deg=55,roll_deg=-30),
        propagation='x(T)=Rz(te)r0+(T-te)Rz(te)n0; c=1; no rotation of emitted photons',
        display_extent_max_age_rlc=renderer.max_age,worldscale_pixels_per_rlc=renderer.worldscale,
        births_per_packet_per_turn=a.births,shutter_phase=renderer.shutter,raster_sigma_reference_pixels=2.8,
        display_weight='Each packet equal integrated 0.1–100GeV photon number; band fractions color mixed in linear light',
        velocity_annotations='160 sparse sampled packet velocity arrows; calculated positions/directions, display-only length 0.24RLC, not field lines or extra emissivity',
        palette_bands_gev=renderer.meta['bands_gev'],palette_srgb_255=(renderer.palette*255).astype(int).tolist(),
        color='Explicit false color: increasing gamma-ray energy from blue through cyan and cream to orange; not visible light',
        exposure=renderer.exposure,exposure_selection='2 times 99.5 percentile nonzero pixel peaks over 12 phases, fixed',
        receiver_distance_display_rlc=renderer.distance,receiver='Fixed bearing glyph; true far-observer time convention, original finite 2deg polar bin; pulse trace is full original profile, not finite dot counts',
        star_marker='2.1-reference-pixel position marker; not a resolved stellar surface',
        time_scale=8/renderer.meta['period_s'],overlay_labels=['@sideratlas','THE PULSAR','NO.002'],
        loop_endpoint_max_pixel_difference=endpoint,export=encoding_metadata(),quantitative_vela_validation=False,
        source_hashes={p.name:sha(p) for p in [Path(__file__),ROOT/'prepare_photon_scene.py',ROOT/'src/pulsar/photon_scene.py']},
        elapsed_seconds=time.monotonic()-start)
    (a.output/'render-metadata.json').write_text(json.dumps(meta,indent=2)+'\n');print(a.output,flush=True)

if __name__=='__main__':main()
