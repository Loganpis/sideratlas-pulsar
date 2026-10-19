"""Package only the current gamma-ray renderer, inputs, notes and media."""
from pathlib import Path
import hashlib,json,zipfile

ROOT=Path(__file__).resolve().parent
VERSION='photon-v1.1.0'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    verification=json.loads((ROOT/'verification.json').read_text())
    export=json.loads((ROOT/'media/export-verification.json').read_text())
    if not verification['passed'] or not export['passed']:
        raise ValueError('Verification must pass before packaging')
    required=['README.md','instagram-description.text','docs/MODEL-NOTES.md',
        'scene/photons.npz','media/sideratlas-002-pulsar-reel.mp4',
        'media/render-metadata.json','media/export-verification.json']
    for name in required:
        if not (ROOT/name).is_file():raise FileNotFoundError(name)
    files={}
    for p in ROOT.iterdir():
        if p.is_file() and (p.suffix in {'.py','.md','.txt','.toml','.json'} or p.name in {'LICENSE','NOTICE','.gitignore','.gitattributes'}):
            files[str(p.relative_to(ROOT))]=p
    for folder in ['src','tests','docs','scene','fonts','physics']:
        for p in (ROOT/folder).rglob('*'):
            if not p.is_file() or any(part in {'__pycache__','.pytest_cache','.git','native-work'} for part in p.parts):continue
            if str(p.relative_to(ROOT)).startswith('physics/outputs/'):continue
            files[str(p.relative_to(ROOT))]=p
    for name in ['sideratlas-002-pulsar-reel.mp4','preview-00.png','contact-sheet.jpg',
        'pulse-timing-check.png','observer-frame-flux.csv','render-metadata.json',
        'export-verification.json','photon-direction-timing-verification.json']:
        p=ROOT/'media'/name
        if not p.exists():raise FileNotFoundError(p)
        files[str(p.relative_to(ROOT))]=p
    manifest=dict(version=VERSION,scope='Current calculated gamma-ray photon animation only',
        verification_passed=True,scene_sha256=sha(ROOT/'scene/photons.npz'),
        video='media/sideratlas-002-pulsar-reel.mp4',caption='instagram-description.text',
        files={name:dict(bytes=p.stat().st_size,sha256=sha(p)) for name,p in sorted(files.items())})
    release=ROOT/'release';release.mkdir(exist_ok=True)
    path=release/'manifest.json';path.write_text(json.dumps(manifest,indent=2)+'\n')
    files['release/manifest.json']=path
    dist=ROOT/'dist';dist.mkdir(exist_ok=True);archive=dist/'sideratlas-002-photon-package.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,p in sorted(files.items()):
            method=zipfile.ZIP_STORED if p.suffix in {'.mp4','.png','.jpg','.npz'} else zipfile.ZIP_DEFLATED
            z.write(p,Path('SiderAtlas-NO002-Photon-Reel')/name,compress_type=method)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('Invalid archive CRC')
    assets=[archive,ROOT/'media/sideratlas-002-pulsar-reel.mp4',ROOT/'instagram-description.text',
        ROOT/'docs/MODEL-NOTES.md',ROOT/'verification.json']
    (release/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in assets))
    print(json.dumps(dict(version=VERSION,archive=str(archive),bytes=archive.stat().st_size,
        sha256=sha(archive),packaged_files=len(files)),indent=2))

if __name__=='__main__':main()
