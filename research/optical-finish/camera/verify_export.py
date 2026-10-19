"""Measure the delivered H.264 file rather than trusting requested settings."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np


def probe(path):
    return json.loads(subprocess.run(['ffprobe','-v','error','-count_frames','-show_streams','-show_format',
                                     '-of','json',str(path)],check=True,capture_output=True,text=True).stdout)


def decode_frames(path,indices,width,height):
    expression='+'.join(f'eq(n,{i})' for i in indices)
    p=subprocess.run(['ffmpeg','-hide_banner','-loglevel','error','-i',str(path),
                      '-vf',f"select='{expression}'",'-fps_mode','passthrough','-f','rawvideo',
                      '-pix_fmt','rgb24','-'],check=True,capture_output=True)
    return np.frombuffer(p.stdout,np.uint8).reshape(-1,height,width,3)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('movie',type=Path)
    parser.add_argument('--loop-frames',type=int,default=240)
    parser.add_argument('--repeats',type=int,default=4)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    p=probe(args.movie)
    video=[x for x in p['streams'] if x['codec_type']=='video']
    if len(video)!=1:raise ValueError('Exactly one video stream required')
    s=video[0]
    expected={'codec_name':'h264','pix_fmt':'yuv420p','width':1080,'height':1920,
              'r_frame_rate':'30/1','avg_frame_rate':'30/1','color_space':'bt709',
              'color_transfer':'bt709','color_primaries':'bt709','color_range':'tv'}
    checks={key:s.get(key)==value for key,value in expected.items()}
    checks['frames']=int(s.get('nb_read_frames',0))==args.loop_frames*args.repeats
    checks['duration']=abs(float(p['format']['duration'])-args.loop_frames*args.repeats/30)<1/300
    checks['no_audio']=all(x['codec_type']!='audio' for x in p['streams'])
    indices=[j*args.loop_frames+i for j in range(args.repeats) for i in (0,args.loop_frames//2,args.loop_frames-1)]
    frames=decode_frames(args.movie,indices,1080,1920)
    # Selected frames are sorted by decoder order, three per repeated cycle.
    repeats=frames.reshape(args.repeats,3,1920,1080,3)
    checks['decoded_repeats_identical']=bool(np.array_equal(repeats,np.broadcast_to(repeats[:1],repeats.shape)))
    seam=np.sqrt(np.mean((repeats[0,-1].astype(float)-repeats[0,0])**2))
    report={'file':str(args.movie),'sha256':hashlib.sha256(args.movie.read_bytes()).hexdigest(),
            'checks':checks,'all_export_checks_pass':all(checks.values()),'measured_stream':s,
            'duration_seconds':float(p['format']['duration']),
            'decoded_cycle_last_to_first_rgb_rmse':float(seam),
            'seam_note':'Adjacent phases of a periodic source are not identical frames; repeated cycles must be identical',
            'scientific_validation':'This verifies delivery encoding only, not the Vela physical model'}
    target=args.output or args.movie.with_suffix('.verification.json')
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_export_checks_pass':report['all_export_checks_pass'],'checks':checks},indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
