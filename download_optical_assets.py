"""Fetch the exact large optical event archive from this repository's release."""
import argparse
import hashlib
from pathlib import Path
import tempfile
import urllib.request

URL='https://github.com/Loganpis/sideratlas-pulsar/releases/download/optical-v1.0.0/physical-pair-events.npz'
SHA256='c2044bb43659951d4a7f78830f8be4b71517f462683fcb2f00cf1e3be92451f8'
SIZE=124646782


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,default=Path(__file__).resolve().parent/'scene/optical')
    args=p.parse_args()
    args.directory.mkdir(parents=True,exist_ok=True)
    target=args.directory/'physical-pair-events.npz'
    if target.exists():
        if target.stat().st_size==SIZE and digest(target)==SHA256:
            print(f'Exact optical archive already present: {target}');return
        raise SystemExit(f'Existing file differs; preserve or move it before downloading: {target}')
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(dir=args.directory,prefix='.optical-',delete=False) as output:
            temporary=Path(output.name)
            with urllib.request.urlopen(URL,timeout=120) as incoming:
                size=0
                for block in iter(lambda:incoming.read(1<<20),b''):
                    size+=len(block)
                    if size>SIZE:raise ValueError('Download exceeds recorded input size')
                    output.write(block)
        if temporary.stat().st_size!=SIZE or digest(temporary)!=SHA256:
            raise ValueError('Downloaded optical archive does not match release hash/size')
        if target.exists():raise FileExistsError(target)
        temporary.rename(target);temporary=None
        print(f'Verified optical archive: {target}')
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)


if __name__=='__main__':main()
