#!/usr/bin/env python3
"""Record optical production prerequisites without promoting tests into physics.

This does not render an image or certify a Vela model. Local kernel checks,
a periodic pair-plasma smoke run and a reduced oblique pulsar prototype are
separate from physical Vela data, emissivity and a stellar atmosphere.
"""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parent

REQUIRED_EVIDENCE = {
    'pulsar_pic_dataset': 'An evolved oblique 3D pulsar PIC dataset with particle/field feedback, IDs, charges, weights, times, momenta and physical normalization',
    'kinetic_convergence': 'Grid, timestep, particle number, domain, source-size and runtime convergence for the target kinetic observables',
    'physical_optical_pairs': 'Physical pair energies, pitch angles and their evolution; a justified mapping from any reduced PIC scales',
    'directional_optical_radiation': 'Observer-angle-dependent 400–700 nm emission, retarded arrival times and independently checked radiation normalization',
    'optical_surface_atmosphere': 'Optical, angle-dependent magnetized atmosphere intensities with supported mass/radius/temperature/composition; no X-ray-table extrapolation',
    'surface_ray_transport': 'Surface gravitational redshift, ray bending, occultation and consistent matching to the outer camera',
    'vela_comparison': 'Separate recorded optical spectrum, normalization and pulse comparisons, identifying fitted and independently assessed data',
    'full_frame_camera': 'A calculated 1080×1920 optical image sequence with a single physical pixel scale, no letterboxing, fixed recorded exposure and no schematic emitters',
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inspect_evidence(manifest, root):
    """Check presence and byte identity only; never infer scientific validity."""
    records = {}
    for key, description in REQUIRED_EVIDENCE.items():
        entry = manifest.get(key)
        record = {'requirement': description, 'artifact_present_and_hashed': False,
                  'scientific_review': 'not_reviewed'}
        if isinstance(entry, dict):
            relative = entry.get('path', '')
            expected = entry.get('sha256', '')
            path = (root/relative).resolve()
            if relative and path.is_relative_to(root.resolve()) and path.is_file():
                actual = sha256(path)
                record.update(path=relative, actual_sha256=actual,
                              artifact_present_and_hashed=actual == expected)
        records[key] = record
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--production-evidence', type=Path,
                   default=ROOT/'data/optical/production-evidence.json')
    p.add_argument('--output', type=Path, default=ROOT/'release/optical-rebuild/status.json')
    a = p.parse_args()
    manifest = json.loads(a.production_evidence.read_text()) if a.production_evidence.exists() else {}
    evidence = inspect_evidence(manifest, ROOT)
    benchmark_path = ROOT/'outputs/entity-pic-benchmark/verification.json'
    benchmark = json.loads(benchmark_path.read_text()) if benchmark_path.exists() else {}
    benchmark_ok = (benchmark.get('status') == 'local_pic_foundation_checks_passed'
                    and benchmark.get('coupled_pic_smoke_test', {}).get('exit_code') == 0
                    and len(benchmark.get('upstream_kernel_checks', [])) == 4
                    and all(x.get('exit_code') == 0 for x in benchmark.get('upstream_kernel_checks', [])))
    local_manifest = ROOT/'outputs/pic-local/native-results.json'
    local_records = {}
    if local_manifest.exists():
        local = json.loads(local_manifest.read_text())
        for relative, record in local.get('files', {}).items():
            path = (ROOT/relative).resolve()
            local_records[relative] = bool(path.is_relative_to(ROOT) and path.is_file()
                                          and sha256(path) == record.get('sha256'))
    available = all(r['artifact_present_and_hashed'] for r in evidence.values())
    report = {
        'status': 'requires_scientific_review' if available else 'production_inputs_missing',
        'scope': 'Prerequisite and provenance audit, not a physical validation or rendering certificate',
        'requested_product': 'Vela optical camera reconstruction with a physical surface atmosphere, full-frame Reel and published source',
        'requested_computation': 'Local computation',
        'requested_components': ['Magnetospheric plasma', 'Thin stellar surface atmosphere'],
        'local_pic_foundation_checks_passed': benchmark_ok,
        'local_pic_benchmark_is_production_input': False,
        'local_reduced_pulsar_evidence_manifest': str(local_manifest.relative_to(ROOT)),
        'local_reduced_pulsar_artifacts_present_and_hashed': bool(local_records) and all(local_records.values()),
        'local_reduced_pulsar_file_checks': local_records,
        'local_reduced_pulsar_is_production_input': False,
        'production_vela_simulation_available': False,
        'optical_camera_movie_available': False,
        'validated_vela_atmosphere_available': False,
        'publish_as_completed_vela_reconstruction': False,
        'previous_illustration_qualifies_as_optical_camera': False,
        'required_evidence': evidence,
        'missing_requirements': [key for key, r in evidence.items() if not r['artifact_present_and_hashed']],
        'camera_delivery_contract': {
            'width': 1080, 'height': 1920, 'fps': 30,
            'target_duration_seconds': 32,
            'wavelength_band_nm': [400, 700],
            'fill_entire_canvas': True, 'letterboxing': False,
            'spatial_scale': 'One physical pixel scale shared by star and magnetosphere; crop based on calculated emission, no independent stellar enlargement',
            'overlays': ['@sideratlas', 'THE PULSAR', 'NO.002'],
            'bottom_text': False,
            'optical_layers': 'Calculated photon radiance only; no luminous field guides, species-colored tracks, artificial cones or arbitrary surface shading',
            'color': 'CIE 1931 spectral integration; fixed recorded exposure to sRGB, then color-managed Rec.709 export',
            'timeline': 'Must record actual simulated physical duration and slowdown; any periodic construction must be disclosed',
        },
        'known_atmosphere_limitation': 'Official NSMAXG hydrogen grids start at 0.047863 keV (~25.9 nm), outside 400–700 nm, and do not contain the required local angle-dependent optical intensity',
        'benchmark_evidence': str(benchmark_path.relative_to(ROOT)),
        'benchmark_sha256': sha256(benchmark_path) if benchmark_path.exists() else None,
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'local_pic_foundation_checks_passed',
          'optical_camera_movie_available', 'publish_as_completed_vela_reconstruction', 'missing_requirements')}, indent=2))


if __name__ == '__main__':
    main()
