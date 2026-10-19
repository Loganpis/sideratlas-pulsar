#!/usr/bin/env python3
"""Aggregate the bounded source tests; never mark an optical model ready."""
from pathlib import Path
import hashlib
import json
import math

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs'
CASES={
    'original48':('local-pic-one-rotation-v2','source-audit.json'),
    'compatible48_core80':('pic-source-repair-imported-one-turn-v2','source-audit.json'),
    'compatible48_core95':('pic-source-repair-core95-one-turn','source-audit.json'),
    'compatible48_no_feedback':('pic-source-repair-no-feedback-one-turn','source-audit.json'),
    'compatible96_core80':('pic-source-repair-resolution96-one-turn','source-audit.json'),
}


def read(path): return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def assemble():
    result={
      'scope':'Bounded rotating-conductor source repair tests of a reduced local kinetic prototype',
      'optical_prediction':False,'physical_Vela_scaling':False,
      'validated_maintained_stellar_surface':False,'kinetic_convergence':False,
      'source_change':'Only conductor E inside the star changes. Compatible temporal edge-potential difference plus compatible vertex-potential gradient; B never reset; native Faraday/Ampere/plasma currents remain active',
      'cases':{},'files':{},
      'limits':[
       'The guaranteed core matches prescribed magnetization by construction; this is a boundary implementation check, not a prediction of Vela magnetic geometry',
       'Same-mask 48-cubed comparisons use the original r<Rstar-dx staggered interior fit mask',
       'Each guaranteed-core fit has its own documented radius; those radii must not be confused with the stellar surface',
       'Surface radial/tangential harmonic fits include non-dipolar distortion and interpolation, and retain substantial residuals',
       '96-cubed resolution control preserves reduced-model source parameters and normalized injection rate, but is not a completed resolution/particle/statistical convergence study',
       'Higher numerical resolution does not establish a mapping to Vela microscopic scales or optical photons',
       'Outer damping, startup transient, prescribed supply, finite enlarged star and incomplete source/boundary energy accounting remain',
       'No atmosphere, physical emission/radiation reaction, camera view or new Reel is produced',
      ]}
    for key,(directory,source_name) in CASES.items():
        folder=OUT/directory
        s=read(folder/source_name)
        v=read(folder/'verification.json')
        first=read(folder/'native'/v['states'][0]['metadata'])
        final=read(folder/'native'/v['states'][-1]['metadata'])
        row={
          'run_directory':str(folder.relative_to(ROOT)),
          'exit_code':v['exit_code'],'wall_seconds':v['wall_seconds'],
          'N':first['active_resolution'][0],'dx':first['dx'],'dt':first['dt'],
          'final_position_E_time':final['position_and_E_time'],
          'final_B_time':final['momentum_B_and_J_time'],
          'normalized_injected_number_flux_per_species':None,
          'guaranteed_or_original_interior_fit':s['summary'],
          'final_surface_radial_fit':s['states'][-1]['spherical_shell_fits'][1]['radial'],
          'final_surface_tangential_fit':s['states'][-1]['spherical_shell_fits'][1]['tangential'],
          'external_pgen_sha256':v['external_pgen_sha256'],
          'native_executable_sha256':v['native_executable_sha256'],
        }
        import tomllib
        config_candidates=[p for p in folder.glob('*.toml')]
        if len(config_candidates)!=1:raise ValueError('Ambiguous run input')
        cfg=tomllib.loads(config_candidates[0].read_text())
        pairs=cfg['setup']['pairs_per_step']
        row['pairs_per_step']=pairs
        row['normalized_injected_number_flux_per_species']=first['q0']*pairs/first['dt']
        row['anchor_core_fraction']=first.get('anchor_core_fraction')
        for report_name,slot in [('original-body-mask-source-audit.json','original_48_style_interior_fit'),('independent-audit.json','native_export_audit'),('independent-charge-audit.json','charge_export_audit')]:
          p=folder/report_name
          if p.exists():
            a=read(p)
            if slot=='original_48_style_interior_fit':row[slot]=a['summary']
            elif slot=='native_export_audit':row[slot]={'inventory_matches':a['particle_inventory_matches_birth_death_events'],'final_snapshot':a['snapshots'][-1]}
            else:row[slot]={'all_snapshots_pointwise_match_full_storage':a['all_snapshots_pointwise_match_full_storage'],'snapshot_count':len(a['snapshots'])}
        result['cases'][key]=row
        for filename in (source_name,'verification.json','original-body-mask-source-audit.json','independent-audit.json','independent-charge-audit.json','native/ledger.csv'):
          path=folder/filename
          if path.is_file():result['files'][str(path.relative_to(ROOT))]={'sha256':sha(path),'bytes':path.stat().st_size}
    c48=result['cases']['compatible48_core80'];c96=result['cases']['compatible96_core80']
    f48=c48['normalized_injected_number_flux_per_species'];f96=c96['normalized_injected_number_flux_per_species']
    assert math.isclose(f48,f96,rel_tol=1e-14)
    result['resolution_control']={
       'normalized_injection_flux_ratio_96_to_48':f96/f48,
       'matched_final_position_E_time_absolute_difference':abs(c48['final_position_E_time']-c96['final_position_E_time']),
       'why_pairs_per_step_4x':'dt halves and q0 scales as dx^3; q0*pairs/dt remains equal',
       'physical_Vela_number_flux_assigned':False,
       'initial_batch_difference':'The dt-proportional initial batch is half as much normalized charge-equivalent number at 96^3; it is not an identical discrete particle realization',
       'other_resolution_changes':'Injection shell Rstar+0.7dx approaches the surface with dx; more macro particles reduce sampling noise; surface staircasing and kinetic resolution both change',
    }
    result['feedback_control']=read(OUT/'pic-source-repair-imported-one-turn-v2/feedback-comparison.json')
    result['algebra_checks']=read(OUT/'pic-source-repair-algebra.json')
    main_receipt=read(OUT/'pic-source-repair-imported-one-turn-v2/verification.json')
    result['build_provenance']={k:main_receipt[k] for k in ['external_pgen_sha256','native_executable_sha256','source_pins','build_commands','imported_kokkos_package','upstream_receipt']}
    result['authored_files']={str(p.relative_to(ROOT)):{'sha256':sha(p),'bytes':p.stat().st_size}for p in (ROOT/'research/pic-source-repair').glob('*') if p.is_file()}
    result['incomplete_attempts']={
       'outputs/pic-source-repair-ct-one-turn':'Fresh full dependency build was interrupted during Kokkos compilation; no simulation run or success receipt',
       'outputs/pic-source-repair-imported-one-turn':'Initial compiler attempt failed because C-array copy needed explicit element copy; preserved build.log; no simulation success',
       'outputs/pic-source-repair-one-turn':'First source-scan attempt was interrupted before any output directory was created',
    }
    (OUT/'pic-source-repair-report.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    x=assemble()
    print(json.dumps({k:{'surface_amplitude_ratio':v['final_surface_radial_fit']['magnitude_relative_to_initial_seed'],'surface_phase_error_deg':v['final_surface_radial_fit']['phase_error_deg_vs_rigid_rotation_at_B_time'],'surface_fit_residual':v['final_surface_radial_fit']['relative_fit_residual_L2']}for k,v in x['cases'].items()},indent=2))
