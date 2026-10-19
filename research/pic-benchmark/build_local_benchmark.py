#!/usr/bin/env python3
"""Prepare or explicitly fetch/build the pinned Entity CPU PIC benchmark.

No flags: inspect existing checkouts and print the plan; no network or build.
--fetch: fetch missing exact commits into new directories, never overwrite.
--build: configure/build the five known local targets; does not run tests.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
THIRD_PARTY = ROOT / "third_party"
BUILD = THIRD_PARTY / "entity-build-research"
SOURCES = (
    {"name": "entity", "directory": "entity-research",
     "repository": "https://github.com/entity-toolkit/entity.git",
     "commit": "62ab34fa0a0d1d6965918abb16d1d0c4b4740cf3",
     "license": "BSD-3-Clause"},
    {"name": "kokkos", "directory": "kokkos-research",
     "repository": "https://github.com/kokkos/kokkos.git",
     "commit": "89f0b12e43583d02bd19739810ceadcedf1ccd56",
     "license": "Apache-2.0 WITH LLVM-exception"},
)
TARGETS = ("entity.xc", "test-kernels-pusher.xc", "test-kernels-deposit.xc",
           "test-kernels-faraday_mink.xc", "test-kernels-ampere_mink.xc")


def git_output(path, *arguments):
    return subprocess.check_output(["git", "-C", str(path), *arguments],
                                   text=True, stderr=subprocess.STDOUT).strip()


def check_source(source, path):
    if path.is_symlink():
        raise RuntimeError(f"Refusing a symbolic-link source directory: {path}")
    if not path.exists():
        return "missing"
    if not path.is_dir() or not (path / ".git").exists():
        raise RuntimeError(f"Refusing existing non-checkout path: {path}")
    actual = git_output(path, "rev-parse", "HEAD")
    if actual != source["commit"]:
        raise RuntimeError(f"Refusing non-pinned {source['name']} checkout: {actual}")
    changes = git_output(path, "status", "--porcelain", "--untracked-files=all")
    if changes:
        raise RuntimeError(f"Refusing modified or untracked source files in {path}:\n{changes}")
    return "verified"


def run(command, *, cwd=ROOT, log=None, env=None):
    print(shlex.join(map(str, command)), flush=True)
    if log is None:
        subprocess.run(list(map(str, command)), cwd=cwd, env=env, check=True)
    else:
        log.write("\n$ " + shlex.join(map(str, command)) + "\n")
        log.flush()
        subprocess.run(list(map(str, command)), cwd=cwd, env=env, check=True,
                       stdout=log, stderr=subprocess.STDOUT)


def fetch_missing(source, path):
    if check_source(source, path) == "verified":
        return
    THIRD_PARTY.mkdir(parents=True, exist_ok=True)
    # Prepare a uniquely named sibling, then move it into the absent target.
    # A failed fetch leaves its directory for inspection; nothing is deleted.
    temporary = Path(tempfile.mkdtemp(prefix=f".{source['name']}-fetch-",
                                      dir=THIRD_PARTY))
    environment = os.environ.copy()
    environment["GIT_LFS_SKIP_SMUDGE"] = "1"  # Images are not required to build.
    try:
        run(["git", "init", temporary], env=environment)
        run(["git", "-C", temporary, "remote", "add", "origin",
             source["repository"]], env=environment)
        run(["git", "-C", temporary, "fetch", "--depth", "1", "origin",
             source["commit"]], env=environment)
        run(["git", "-C", temporary, "checkout", "--detach", "FETCH_HEAD"],
            env=environment)
        check_source(source, temporary)
        if path.exists() or path.is_symlink():
            raise RuntimeError(f"Destination appeared during fetch; refusing to replace {path}")
        temporary.rename(path)
    except Exception:
        print(f"Fetch did not finish; inspection directory retained: {temporary}", flush=True)
        raise


def checked_cmake(value):
    candidate = shutil.which(value)
    if candidate is None:
        path = Path(value).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            candidate = str(path.resolve())
    if candidate is None:
        raise RuntimeError("CMake not found; provide --cmake /absolute/path/to/cmake")
    return str(Path(candidate).resolve())


def configure_arguments(cmake):
    result = [cmake, "-S", str(THIRD_PARTY / "entity-research"), "-B", str(BUILD),
              "-Doutput=OFF", "-Dmpi=OFF", "-Dpgen=streaming",
              "-Dprecision=double", "-DTESTS=ON",
              "-DKokkos_ENABLE_SERIAL=ON", "-DKokkos_ENABLE_OPENMP=OFF",
              "-DKokkos_VERSION=5.2.1",
              "-DFETCHCONTENT_SOURCE_DIR_KOKKOS=" + str(THIRD_PARTY / "kokkos-research"),
              "-DCONNECTION_CHECKED=TRUE", "-DFETCHCONTENT_FULLY_DISCONNECTED=OFF"]
    if platform.machine().lower() in ("arm64", "aarch64"):
        result.append("-DKokkos_ARCH_ARMV8=ON")
    if platform.system() == "Darwin":
        result += ["-DCMAKE_CXX_COMPILER=/usr/bin/clang++",
                   "-DCMAKE_C_COMPILER=/usr/bin/clang"]
    return result


def check_build_directory():
    if BUILD.is_symlink():
        raise RuntimeError(f"Refusing symbolic-link build directory: {BUILD}")
    if not BUILD.exists():
        return
    if not BUILD.is_dir():
        raise RuntimeError(f"Build path exists and is not a directory: {BUILD}")
    cache = BUILD / "CMakeCache.txt"
    if not cache.exists():
        if any(BUILD.iterdir()):
            raise RuntimeError("Existing nonempty build directory has no CMake cache; "
                               "move it aside explicitly before building")
        return
    entries = dict(line.split("=", 1) for line in cache.read_text().splitlines()
                   if "=" in line and not line.startswith(("#", "//")))
    home = entries.get("CMAKE_HOME_DIRECTORY:INTERNAL")
    if home is None or Path(home).resolve() != (THIRD_PARTY / "entity-research").resolve():
        raise RuntimeError("Existing build cache belongs to another source directory")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true",
                        help="Fetch missing exact source commits using HTTPS")
    parser.add_argument("--build", action="store_true",
                        help="Configure and build only the five benchmark targets")
    parser.add_argument("--cmake", default="cmake", help="CMake executable path")
    parser.add_argument("--jobs", type=int, default=2,
                        help="Parallel compiler jobs (default: 2)")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")

    states = [{**source, "status": check_source(source, THIRD_PARTY / source["directory"])}
              for source in SOURCES]
    print(json.dumps({"sources": states, "build_directory": str(BUILD),
                      "targets": list(TARGETS), "jobs": args.jobs,
                      "actions": {"fetch": args.fetch, "build": args.build}}, indent=2))
    if args.fetch:
        for source in SOURCES:
            fetch_missing(source, THIRD_PARTY / source["directory"])
    if not args.build:
        print("No build requested; no executable or benchmark result was created.")
        return
    for source in SOURCES:
        if check_source(source, THIRD_PARTY / source["directory"]) != "verified":
            raise RuntimeError("Pinned sources are missing; rerun with --fetch --build")
    check_build_directory()
    cmake = checked_cmake(args.cmake)
    output = ROOT / "outputs/entity-pic-benchmark"
    output.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    descriptor, name = tempfile.mkstemp(prefix=f"build-{stamp}-", suffix=".log", dir=output)
    with os.fdopen(descriptor, "w") as log:
        run(configure_arguments(cmake), log=log)
        run([cmake, "--build", BUILD, "--target", *TARGETS, "--parallel", args.jobs], log=log)
    for source in SOURCES:
        check_source(source, THIRD_PARTY / source["directory"])
    print(f"Build complete; unique build log: {name}")
    print("Run research/pic-benchmark/run_local_benchmark.py to record kernel checks "
          "and the bounded 3D pair-plasma smoke test")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error
