#!/usr/bin/env python3
"""Prepare independent pinned COFFEE variants from the bundled pristine source.

Requires Python, patch, CMake, a C++17 compiler, OpenMPI and parallel HDF5.
Use --prepare-only to inspect the complete patched corresponding source without
compiling. No source tree or existing build is overwritten.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
CHAINS = {
    "portable": ["portable"],
    "surface-anchored": ["portable", "surface-anchored"],
    "curl-conductor": ["portable", "surface-anchored", "curl-conductor"],
    "curl-ohmic": ["portable", "surface-anchored", "curl-conductor", "curl-ohmic"],
    "curl-radiative": ["portable", "surface-anchored", "curl-conductor", "curl-radiative"],
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("variant", choices=CHAINS)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--output", type=Path, help="Fresh work directory; defaults to native-work/VARIANT")
    ap.add_argument("--local-layout", action="store_true", help="Prepare third_party/coffee-VARIANT for the historical run entry points")
    ap.add_argument("--dependency-prefix", type=Path, help="macOS MPI/HDF5 prefix; omit for system package discovery")
    ap.add_argument("--cmake", default="cmake")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    default = ROOT / "third_party" / f"coffee-{args.variant}" if args.local_layout else ROOT / "native-work" / args.variant
    output = (args.output or default).absolute()
    source = output if args.local_layout else output / "source"
    if output.exists():
        raise FileExistsError(f"Use a new work directory: {output}")
    pristine = ROOT / "third_party/coffee-src"
    if not (pristine / "LICENSE").is_file():
        raise FileNotFoundError("The release must include third_party/coffee-src and its unchanged license")
    source.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(pristine, source, ignore=shutil.ignore_patterns(".git", ".vscode", "__pycache__", "*.pyc"))
    patches = []
    for variant in CHAINS[args.variant]:
        patch = ROOT / f"third_party/coffee-{variant}.patch"
        with patch.open("rb") as stream:
            subprocess.run(["patch", "--batch", "--forward", "-p1"], cwd=source, stdin=stream, check=True)
        patches.append(dict(path=str(patch.relative_to(ROOT)), sha256=sha(patch)))
    record = dict(upstream_commit="9986522d08b4d887e5ef9a53e12cecae4e376d22", variant=args.variant,
                  patches=patches, source=str(source), prepare_only=args.prepare_only)
    if not args.prepare_only:
        command = [args.cmake, "-S", str(source), "-B", str(output / "build"),
            "-DCMAKE_POLICY_VERSION_MINIMUM=3.5", "-DCMAKE_BUILD_TYPE=Release",
            "-Duse_cuda=0", "-Duse_double=1", "-Dportable_scalar=ON", "-Dbuild_tests=OFF"]
        if args.dependency_prefix:
            command.append(f"-DCOFFEE_DEPENDENCY_PREFIX={args.dependency_prefix.absolute()}")
        subprocess.run(command, check=True)
        subprocess.run([args.cmake, "--build", str(output / "build"), "--target", "coffee-3d", "coffee-verification", "-j", str(args.jobs)], check=True)
        record["configure_command"] = command
        record["binaries"] = {name: sha(source / "bin" / name) for name in ("coffee-3d", "coffee-verification")}
    (output / "source-manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
