#!/usr/bin/env python3
"""Build/run the scaled oblique 3D kinetic prototype in a new output directory."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BUILD = ROOT / "outputs/pic-source-repair-build"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_build_receipt(path, sources, executable, pgen):
    """Reuse exactly a previously built, successfully executed native binary.

    This authenticates byte identity against local build evidence, not a fresh
    upstream working-tree scan or the scientific validity of the model.
    """
    receipt = json.loads(path.read_text())
    if (receipt.get('exit_code') != 0 or not receipt.get('build_commands')
            or receipt.get('reused_verified_build')
            or receipt.get('source_pins') != json.loads(json.dumps(sources))):
        raise ValueError('Receipt must record a direct successful build of the exact pinned sources')
    if receipt.get('external_pgen_sha256') != digest(pgen):
        raise ValueError('Problem generator differs from the verified build receipt')
    if receipt.get('native_executable_sha256') != digest(executable):
        raise ValueError('Native binary differs from the verified build receipt')
    return {'path': str(path.resolve()), 'sha256': digest(path),
            'meaning': 'Exact prior binary and problem-generator bytes; no fresh upstream tree scan'}


def main():
    global BUILD
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--build-directory",type=Path,default=BUILD)
    parser.add_argument("--kokkos-package",type=Path,help="Use an existing KokkosConfig.cmake package; records imported library hashes")
    parser.add_argument("--upstream-receipt",type=Path,help="Use a prior direct clean-source build receipt and verify current HEAD pins; records that worktree bytes were not freshly rescanned")
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--config", type=Path, default=HERE / "one-rotation.toml")
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument('--reuse-verified-build', type=Path,
                        help='Reuse the exact binary from a successful direct-build verification.json; refuses changed bytes')
    args = parser.parse_args()
    BUILD=args.build_directory.resolve()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if sys.byteorder != "little":
        parser.error("The native binary exporter currently requires a little-endian CPU")
    config = args.config.resolve()
    if not config.is_file():
        parser.error("Configuration file is missing")
    spec = importlib.util.spec_from_file_location("pic_foundation_builder",
                                                 HERE.parent / "pic-benchmark/build_local_benchmark.py")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    executable = BUILD / 'src/entity.xc'
    reused = None
    if args.reuse_verified_build:
        if args.build or args.fetch:
            parser.error('Reusing a verified binary cannot be combined with --build or --fetch')
        try:
            reused = verify_build_receipt(args.reuse_verified_build.resolve(), helper.SOURCES,
                                           executable, HERE/'pgen.hpp')
        except (ValueError, OSError, KeyError) as error:
            parser.error(str(error))
    elif args.upstream_receipt:
        prior = json.loads(args.upstream_receipt.read_text())
        if prior.get('exit_code') != 0 or not prior.get('build_commands') or prior.get('source_pins') != json.loads(json.dumps(helper.SOURCES)):
            parser.error('Upstream receipt must be a successful direct exact-pinned build')
        for source in helper.SOURCES:
            path=helper.THIRD_PARTY/source['directory']
            if path.is_symlink() or helper.git_output(path,'rev-parse','HEAD') != source['commit']:
                parser.error('Current upstream HEAD differs from prior receipt')
    else:
        for source in helper.SOURCES:
            path = helper.THIRD_PARTY / source["directory"]
            if args.fetch:
                helper.fetch_missing(source, path)
            if helper.check_source(source, path) != "verified":
                parser.error("Pinned source is missing; provide --fetch")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_directory = (args.run_dir or ROOT / "outputs" / f"local-pic-{stamp}").resolve()
    if run_directory.exists():
        parser.error("Run directory already exists; use a new path to preserve previous results")
    run_directory.mkdir(parents=True)
    build_commands = []
    if args.build:
        cmake = helper.checked_cmake(args.cmake)
        command = helper.configure_arguments(cmake)
        command[command.index("-B") + 1] = str(BUILD)
        if args.kokkos_package:
            command.append("-DKokkos_DIR="+str(args.kokkos_package.resolve()))
        command = [f"-Dpgen={HERE}" if x.startswith("-Dpgen=") else
                   "-DTESTS=OFF" if x == "-DTESTS=ON" else x for x in command]
        build_commands.append(command)
        build_commands.append([cmake, "--build", str(BUILD), "--target", "entity.xc",
                               "--parallel", str(args.jobs)])
        with (run_directory / "build.log").open("w") as log:
            for command in build_commands:
                helper.run(command, log=log)
        (BUILD / "pgen-source-sha256.txt").write_text(digest(HERE / "pgen.hpp") + "\n")
    if not executable.is_file():
        parser.error("Native executable is missing; provide --build")
    stamp_path = BUILD / "pgen-source-sha256.txt"
    if not stamp_path.exists() or stamp_path.read_text().strip() != digest(HERE / "pgen.hpp"):
        parser.error("Executable source hash is unrecorded or stale; provide --build")
    for filename in ("pgen.hpp", config.name):
        source = HERE / filename if filename == "pgen.hpp" else config
        (run_directory / filename).write_bytes(source.read_bytes())
    command = [str(executable), "-input", str(config)]
    started = time.perf_counter()
    with (run_directory / "run.log").open("w") as log:
        result = subprocess.run(command, cwd=run_directory, stdout=log,
                                stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - started
    states = []
    for path in sorted((run_directory / "native").glob("state-*.json")):
        state = json.loads(path.read_text())
        binary = path.with_name(state["fields_file"])
        particles = path.with_name(state["particles_file"])
        expected = 8
        for size in state["field_shape"]:
            expected *= size
        if binary.stat().st_size != expected:
            raise RuntimeError(f"Incomplete native field export: {binary}")
        states.append({"completed_step": state["completed_step"],
                       "position_and_E_time": state["position_and_E_time"],
                       "metadata": path.name, "metadata_sha256": digest(path),
                       "fields_sha256": digest(binary),
                       "particles_sha256": digest(particles)})
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "imported_kokkos_package": {"path":str(args.kokkos_package.resolve()),"libraries_sha256":{str(p.relative_to(args.kokkos_package)):digest(p) for p in args.kokkos_package.rglob("libkokkos*.a")}} if args.kokkos_package else None,
        "upstream_receipt": {"path":str(args.upstream_receipt.resolve()),"sha256":digest(args.upstream_receipt),"fresh_worktree_scan":False} if args.upstream_receipt else None,
        "scope": "scaled local oblique 3D PIC startup prototype",
        "production_vela_dataset": False, "optical_prediction": False,
        "stellar_atmosphere_model": False, "convergence_claim": False,
        "exit_code": result.returncode, "wall_seconds": elapsed,
        "run_command": command, "build_commands": build_commands,
        "source_pins": helper.SOURCES,
        "external_pgen_sha256": digest(HERE / "pgen.hpp"),
        "configuration_sha256": digest(config),
        "native_executable_sha256": digest(executable),
        "run_log_sha256": digest(run_directory / "run.log"), "states": states,
        'reused_verified_build': reused,
    }
    (run_directory / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"exit_code": result.returncode, "wall_seconds": elapsed,
                      "snapshots": len(states), "run_directory": str(run_directory)}))
    if result.returncode:
        raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
