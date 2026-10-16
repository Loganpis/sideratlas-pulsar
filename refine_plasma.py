"""Continue a coarse periodic candidate on a finer grid; no Vela fitting."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from dataclasses import replace
from scipy.interpolate import RegularGridInterpolator
from run_plasma_upgrade import ROOT,OUT,dump,spherical_flux,compare_fields
from pulsar.plasma import Magnetosphere,PlasmaConfig,energy


def main():
    source=OUT/'resistive-n96.npz'; data=np.load(source)
    config=replace(PlasmaConfig(**json.loads(str(data['config']))),n=128)
    dest=OUT/'resistive-n128-restarted.npz'; hp=OUT/'resistive-n128-restarted-history.json'
    report=OUT/'refinement-validation.json'
    checkpoint=OUT/'refinement-checkpoint.npz'
    input_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    manifest_path=OUT/'refinement-generation.json'
    if dest.exists() and hp.exists() and report.exists():
        prior=json.loads(report.read_text())
        if (prior.get('initial_snapshot_sha256')!=input_hash
                or prior.get('output_snapshot_sha256')!=hashlib.sha256(dest.read_bytes()).hexdigest()):
            raise ValueError('Refinement cache lacks matching input/output provenance; preserve it and audit before reuse')
        print('Reusing completed refinement with its recorded original generation provenance',flush=True)
        return
    generation={'initial_snapshot_sha256':input_hash,
                'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                    [ROOT/'refine_plasma.py',ROOT/'src/pulsar/plasma.py']}}
    if checkpoint.exists():
        if not manifest_path.exists() or json.loads(manifest_path.read_text())!=generation:
            raise ValueError('Checkpoint input or generating source differs; do not silently combine solver versions')
    else:
        dump(manifest_path,generation)
    s=Magnetosphere(config)
    start=float(data['time']); end=start+2*np.pi; history=[]; wall=time.monotonic()
    if checkpoint.exists():
        saved=np.load(checkpoint)
        if json.loads(str(saved['config']))!=config.__dict__:
            raise ValueError('Checkpoint configuration differs from requested refinement')
        s.a=saved['a']; s.e=saved['e']; s.time=float(saved['time'])
        s.initial_energy=float(saved['initial_energy'])
        diag=json.loads(str(saved['diagnostics']))
        s.step_number=diag['steps']; s.ledger=diag['ledger']
        history=[v for v in json.loads(hp.read_text()) if v['time']<=s.time]
        print(f'Resuming saved refinement at extra rotations={(s.time-start)/(2*np.pi):.2f}',flush=True)
    else:
        points=s.xyz.reshape(3,-1).T
        # Extrapolation occurs only in the damped outermost half cell.
        # The physical comparison region is wholly inside both meshes.
        for key in ('a','e'):
            interp=RegularGridInterpolator((data['axis'],)*3,data[key].transpose(1,2,3,0),
                                           bounds_error=False,fill_value=None)
            setattr(s,key,interp(points).T.reshape(s.xyz.shape))
        s.time=start; s.initial_energy=energy(s.e,s.b,s.h)
        del points
        print('Starting 128-cubed refinement from completed 96-cubed fields',flush=True)
    interval=round(.1*2*np.pi/(s.h*config.courant))
    dump(OUT/'refinement-protocol.json',{'initial_snapshot':str(source.relative_to(ROOT)),
         'config':json.loads(json.dumps(config.__dict__)),'relaxation_rotations':1,
         'method':'Interpolate A and E, reconstruct B=curlA, evolve one more full rotation',
         'not_equivalent_to':'A fresh independent high-resolution run or resolved current-sheet microphysics'})
    while s.time<end:
        s.step(min(s.h*config.courant,end-s.time))
        if s.step_number%interval==0 or s.time>=end:
            d=s.diagnostics(); d['poynting_flux_r_1_1p5_2']=spherical_flux(s.axis,s.e,s.b)
            history.append(d); dump(hp,history)
            # Atomic, uncompressed checkpoints avoid losing the computation
            # when a desktop session is interrupted or the computer sleeps.
            temp=OUT/'refinement-checkpoint-next.npz'
            s.save(temp,compressed=False); temp.replace(checkpoint)
            print(f"128³ restart: extra rotations={(s.time-start)/(2*np.pi):.2f}, L(1.5)={d['poynting_flux_r_1_1p5_2'][1]:.3f}, elapsed={time.monotonic()-wall:.0f}s",flush=True)
    s.save(dest)
    result={'96_to_128_restarted':compare_fields(source,dest),
            'last_diagnostics':history[-1],
            'same_rotation_phase':True,'field_refinement_target':.10,
            'separate_high_resolution_startup_run':False,
            'validated_for_Vela':False,
            'generation_source_sha256':generation['source_sha256'],
            'initial_snapshot_sha256':input_hash,
            'output_snapshot_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
            'comparison_kind':'Combined grid refinement and one extra rotation of evolution; use matched-time control for spatial comparison'}
    result['combined_change_below_target']=max(result['96_to_128_restarted'].values())<.10
    result['field_refinement_passed']=False
    dump(OUT/'refinement-validation.json',result)
    checkpoint.unlink(missing_ok=True)


if __name__=='__main__': main()
