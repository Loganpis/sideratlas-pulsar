"""Bundle the ready photon scene, renderer, media and recorded checks."""
from pathlib import Path
import ast,hashlib,json,shutil,zipfile

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'outputs/sideratlas-pulsar-photons'
MEDIA=BASE/'media'

def main():
    for name in ['instagram-description.text','MODEL-NOTES.md']:
        shutil.copy2(ROOT/'release/photon-v2'/name,MEDIA/name)
    checks=json.loads((BASE/'verification-summary.json').read_text())
    files=[ROOT/n for n in ['prepare_photon_scene.py','render_photon_reel.py','reel_export.py','verify_reel.py','LICENSE','NOTICE']]
    files+=list((ROOT/'release/photon-v2').glob('*'))+[ROOT/'release/MODEL-NOTES.md']
    files+=[ROOT/'fonts'/n for n in ['FONT-SOURCES.md','GUST-FONT-LICENSE.txt','lmroman10-regular.otf','lmroman12-bold.otf']]
    modules={'photon_scene','snapshot_fields','emission_quadrature','outer_transport','radiation','model','published_acceleration','trajectories'}
    pending=list(modules)
    while pending:
        name=pending.pop();path=ROOT/'src/pulsar'/f'{name}.py'
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,ast.ImportFrom) and node.level==1 and node.module:
                dep=node.module.split('.')[0]
                if (ROOT/'src/pulsar'/f'{dep}.py').exists() and dep not in modules:
                    modules.add(dep);pending.append(dep)
    files+=[ROOT/'src/pulsar/__init__.py']+[ROOT/'src/pulsar'/f'{m}.py' for m in sorted(modules)]
    files+=[ROOT/'tests'/f'test_{m}.py' for m in ['photon_scene','emission_quadrature','outer_transport','radiation']]
    files+=list((BASE/'scene').glob('*'))
    files+=[MEDIA/n for n in ['sideratlas-002-pulsar-reel.mp4','preview-00.png','contact-sheet.jpg','pulse-timing-check.png','observer-frame-flux.csv','render-metadata.json','export-verification.json','photon-direction-timing-verification.json']]
    files+=[ROOT/'outputs/emission-quadrature/n97-final-illustration-a128-r6'/n for n in ['protocol.json','results.json','profiles.npz']]
    files+=[BASE/'verification-summary.json']
    mapping={str(p.relative_to(ROOT)):p for p in set(files) if p.is_file()}
    mapping['fonts/LPPL-1.3c.txt']=ROOT/'release/photon-v2/LPPL-1.3c.txt'
    mapping['README.md']=ROOT/'release/photon-v2/README.md'
    manifest=dict(scope='Original gamma-ray physics; new outgoing-photon visualization; no remote publication performed',
        all_verifications_passed=all(checks[k]['passed'] for k in ['calculation','directions_and_timing','video','focused_regression_tests']),
        video=str((MEDIA/'sideratlas-002-pulsar-reel.mp4').relative_to(ROOT)),
        caption='release/photon-v2/instagram-description.text',model_notes='release/photon-v2/MODEL-NOTES.md',
        files={name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for name,p in sorted(mapping.items())})
    (BASE/'DELIVERY.json').write_text(json.dumps(manifest,indent=2)+'\n')
    mapping[str((BASE/'DELIVERY.json').relative_to(ROOT))]=BASE/'DELIVERY.json'
    archive=BASE/'sideratlas-002-photon-package.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,p in sorted(mapping.items()):
            compression=zipfile.ZIP_STORED if p.suffix.lower() in {'.mp4','.png','.jpg','.npz'} else zipfile.ZIP_DEFLATED
            z.write(p,Path('SiderAtlas-NO002-Photon-Reel')/name,compress_type=compression)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('Archive CRC verification failed')
    print(json.dumps(dict(all_verifications_passed=manifest['all_verifications_passed'],archive=str(archive),archive_bytes=archive.stat().st_size,code_modules=sorted(modules)),indent=2))

if __name__=='__main__':main()
