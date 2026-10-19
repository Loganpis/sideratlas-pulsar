#!/usr/bin/env python3
"""Inspect final encoded bytes, loop repeats and color-managed frame quality."""
from pathlib import Path
import argparse,hashlib,json,subprocess,struct
import numpy as np
from PIL import Image
from reel_export import BT709_TO_SRGB_FILTER


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path)
    p.add_argument('--video-name',default='sideratlas-002-pulsar-reel.mp4');a=p.parse_args()
    folder=a.folder;video=folder/a.video_name
    metadata=json.loads((folder/'render-metadata.json').read_text())
    info=json.loads(subprocess.check_output(['ffprobe','-v','error','-count_frames','-show_streams','-show_format','-of','json',str(video)]))
    streams=info['streams'];stream=next(x for x in streams if x['codec_type']=='video')
    w,h=stream['width'],stream['height'];fps=metadata['fps'];loop_frames=8*fps;repeats=metadata['repeats']
    ids=[i*loop_frames for i in range(repeats)]
    select='select='+"'"+'+'.join(f'eq(n,{i})' for i in ids)+"'"
    raw=subprocess.check_output(['ffmpeg','-v','error','-i',str(video),'-vf',select+','+BT709_TO_SRGB_FILTER,
        '-fps_mode','vfr','-f','rawvideo','-pix_fmt','rgb24','-'])
    frames=np.frombuffer(raw,np.uint8).reshape(repeats,h,w,3)
    original=np.asarray(Image.open(folder/'preview-00.png').convert('RGB'))
    error=frames[0].astype(float)-original.astype(float);mse=float(np.mean(error**2))
    psnr=float(10*np.log10(255**2/mse)) if mse else 99.
    atoms=[]
    with video.open('rb') as f:
        while head:=f.read(8):
            size,name=struct.unpack('>I4s',head);header=8
            if size==1:size=struct.unpack('>Q',f.read(8))[0];header=16
            atoms.append(name.decode('ascii'))
            if size==0:break
            if size<header:raise ValueError('Malformed MP4 atom')
            f.seek(size-header,1)
    checks=dict(native1080x1920=(w,h)==(1080,1920),h264_high_420=stream['codec_name']=='h264' and stream['profile']=='High' and stream['pix_fmt']=='yuv420p',
        duration32=abs(float(info['format']['duration'])-32)<.01,
        fps30=stream['avg_frame_rate']=='30/1',frame_count960=int(stream['nb_read_frames'])==960,
        rec709_limited=all(stream.get(k)==v for k,v in {'color_space':'bt709','color_transfer':'bt709','color_primaries':'bt709','color_range':'tv','chroma_location':'left'}.items()),
        four_exact_repeats=bool(np.all(frames==frames[0])),faststart='moov' in atoms and atoms.index('moov')<atoms.index('mdat'),
        no_audio=all(s['codec_type']!='audio' for s in streams),no_b_frames=stream['has_b_frames']==0,
        three_labels_only=metadata['overlay_labels']==['@sideratlas','THE PULSAR','NO.002'],
        mathematical_loop_endpoint=metadata['loop_endpoint_max_pixel_difference']==0,
        source_frame_psnr_above35=psnr>35,
        quantitative_claim_disabled=not metadata['quantitative_vela_validation'])
    report=dict(passed=all(checks.values()),checks=checks,scope='Local export and repeated-loop checks; not an Instagram upload or scientific model validation',
        source_frame_psnr_db=psnr,source_frame_mean_absolute_error_codes=float(abs(error).mean()),
        final_sha256=hashlib.sha256(video.read_bytes()).hexdigest(),file_bytes=video.stat().st_size,
        stream={k:stream.get(k) for k in ['codec_name','profile','pix_fmt','width','height','avg_frame_rate','nb_read_frames','bit_rate','color_space','color_transfer','color_primaries','color_range','chroma_location']},
        mp4_atoms=atoms)
    (folder/'export-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    Image.fromarray(frames[0]).save(folder/'encoded-frame-00.png')
    print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit('Final export verification failed')


if __name__=='__main__':main()
