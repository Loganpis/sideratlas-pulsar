"""Same-duration coarse continuation to separate relaxation from grid effects."""
import json,time,hashlib
import numpy as np
from run_plasma_upgrade import ROOT,OUT,dump,compare_fields,spherical_flux
from pulsar.plasma import Magnetosphere,PlasmaConfig,energy,SOLVER_SOURCE_SHA256


def main():
    source=OUT/'resistive-n96.npz'
    dest=OUT/'resistive-n96-continued.npz'
    hp=OUT/'resistive-n96-continued-history.json'
    data=np.load(source)
    generation={'solver_sha256':SOLVER_SOURCE_SHA256,
                'runner_sha256':hashlib.sha256((ROOT/'check_refinement.py').read_bytes()).hexdigest(),
                'initial_snapshot_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
    if not dest.exists():
        s=Magnetosphere(PlasmaConfig(**json.loads(str(data['config']))))
        s.a=data['a'];s.e=data['e'];s.time=float(data['time'])
        s.initial_energy=energy(s.e,s.b,s.h)
        start=s.time;end=start+2*np.pi;history=[];wall=time.monotonic()
        every=max(1,round(.25*2*np.pi/(s.h*s.config.courant)))
        while s.time<end:
            s.step(min(s.h*s.config.courant,end-s.time))
            if s.step_number%every==0 or s.time>=end:
                d=s.diagnostics();d['poynting_flux_r_1_1p5_2']=spherical_flux(s.axis,s.e,s.b)
                history.append(d);dump(hp,history)
                print(f'96-cubed control: extra rotation={(s.time-start)/(2*np.pi):.2f}, elapsed={time.monotonic()-wall:.0f}s',flush=True)
        s.save(dest)
        dump(OUT/'coarse-continuation-generation.json',generation)
    else:
        generation=json.loads((OUT/'coarse-continuation-generation.json').read_text())
        if generation['initial_snapshot_sha256']!=hashlib.sha256(source.read_bytes()).hexdigest():
            raise ValueError('Initial snapshot changed after coarse continuation')
    fine=OUT/'resistive-n128-restarted.npz'
    if fine.exists():
        result={'same_time_96_to_128':compare_fields(dest,fine),
                'coarse_3_to_4_periods':compare_fields(source,dest),
                'same_final_time':bool(np.isclose(np.load(dest)['time'],np.load(fine)['time'])),
                'fine_run_started_from_interpolated_coarse_state':True,
                'independent_fine_startup':False,
                'validated_for_Vela':False,'field_refinement_target':.10,
                'coarse_generation':generation,
                'input_snapshot_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,dest,fine)}}
        result['same_time_field_change_below_target']=max(result['same_time_96_to_128'].values())<.10
        dump(OUT/'matched-time-validation.json',result)
        print(json.dumps(result,indent=2))
    else:
        print('Control saved; run again after the fine snapshot finishes to compare the fields',flush=True)


if __name__=='__main__':main()
