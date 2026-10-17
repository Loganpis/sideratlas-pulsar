#!/usr/bin/env python3
"""Run the pinned Entity CPU checks and 3D PIC smoke test, not a Vela model.

Prerequisite checkouts and build instructions are in PIC-REBUILD-FEASIBILITY.md.
This script neither changes upstream source nor downloads anything.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENTITY = ROOT / "third_party/entity-research"
KOKKOS = ROOT / "third_party/kokkos-research"
BUILD = ROOT / "third_party/entity-build-research"
OUT = ROOT / "outputs/entity-pic-benchmark"
CONFIG = Path(__file__).with_name("periodic-pairs-3d.toml")
PINS = {
    "entity": "62ab34fa0a0d1d6965918abb16d1d0c4b4740cf3",
    "kokkos": "89f0b12e43583d02bd19739810ceadcedf1ccd56",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args: str, cwd: Path) -> str:
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, path in (("entity", ENTITY), ("kokkos", KOKKOS)):
        actual = git("rev-parse", "HEAD", cwd=path)
        if actual != PINS[name]:
            raise RuntimeError(f"{name} source pin changed: {actual}")
        if git("status", "--porcelain", "--untracked-files=no", cwd=path):
            raise RuntimeError(f"{name} tracked source is modified")

    records = []
    for name in ("pusher", "deposit", "faraday_mink", "ampere_mink"):
        executable = BUILD / "tests/kernels" / f"test-kernels-{name}.xc"
        t0 = time.perf_counter()
        result = subprocess.run([str(executable)], cwd=OUT, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = OUT / f"test-{name}.log"
        log.write_text(result.stdout)
        records.append({"name": name, "exit_code": result.returncode,
                        "wall_seconds": time.perf_counter() - t0,
                        "executable_sha256": sha256(executable),
                        "log_sha256": sha256(log)})
        if result.returncode:
            raise RuntimeError(f"Upstream {name} check failed; see {log}")

    executable = BUILD / "src/entity.xc"
    t0 = time.perf_counter()
    result = subprocess.run([str(executable), "-input", str(CONFIG)], cwd=OUT,
                            text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    wall = time.perf_counter() - t0
    log = OUT / "run.log"
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"3D PIC run failed; see {log}")
    if "Resolution: [16, 16, 16]" not in result.stdout:
        raise RuntimeError("Unexpected 3D mesh")
    if "Runtime: 8.000e+00 [55 steps]" not in result.stdout:
        raise RuntimeError("Unexpected run duration or timestep count")
    if "CurrentDeposit" not in result.stdout or "ParticlePusher" not in result.stdout:
        raise RuntimeError("Required PIC stages absent from timing report")
    if re.search(r"\b(?:nan|inf|fatal)\b", result.stdout, flags=re.IGNORECASE):
        raise RuntimeError("Nonfinite/fatal marker in benchmark log")
    counts = [int(x) for x in re.findall(
        r"species\s+\d+\s+\([^\n]+?\)\.+(\d+)\s+\d+%", result.stdout)]
    if not counts or set(counts) != {4096}:
        raise RuntimeError(f"Unexpected particle counts: {counts}")

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "local_pic_foundation_checks_passed",
        "production_vela_dataset": False,
        "visible_light_prediction": False,
        "pulsar_magnetosphere_validation": False,
        "sources": {
            "entity": {"repository": "https://github.com/entity-toolkit/entity",
                       "commit": PINS["entity"], "license": "BSD-3-Clause",
                       "license_sha256": sha256(ENTITY / "LICENSE")},
            "kokkos": {"repository": "https://github.com/kokkos/kokkos",
                       "commit": PINS["kokkos"], "tag": "5.2.1",
                       "license": "Apache-2.0 WITH LLVM-exception",
                       "license_sha256": sha256(KOKKOS / "LICENSE")},
        },
        "build": {"backend": "Serial CPU", "precision": "double",
                  "compiler": "AppleClang 21.0.0.21000334 / C++20",
                  "mpi": False, "adios2_output": False,
                  "upstream_tracked_sources_modified": False},
        "upstream_kernel_checks": records,
        "coupled_pic_smoke_test": {
            "exit_code": result.returncode, "wall_seconds": wall,
            "dimension": 3, "mesh": [16, 16, 16], "boundary": "periodic",
            "input_sha256": sha256(CONFIG),
            "input_path": str(CONFIG.relative_to(ROOT)),
            "executable_sha256": sha256(executable),
            "log_sha256": sha256(log),
            "species": 4, "particles_per_species": 4096,
            "total_macroparticles": 16384, "steps": 55,
            "dt": 0.5 / (3 ** 0.5) * 0.5,
            "requested_runtime": 8.0,
            "simulated_time_at_55_steps": 55 * 0.5 / (3 ** 0.5) * 0.5,
            "field_solver": "enabled", "current_deposition": "enabled",
            "particle_pusher": "Boris", "radiation_drag": "none",
            "emission_policy": "none",
        },
        "limits": [
            "The 3D smoke test establishes that the actual coupled PIC execution path runs locally",
            "It does not measure instability growth, energy conservation or kinetic convergence",
            "The upstream pusher and Maxwell/deposition regression checks test separate algorithms",
            "No pulsar conductor, atmosphere, open boundary, pair cascade or physical Vela scaling",
            "No spectral emissivity, optical radiative transfer or observer camera prediction",
            "No benchmark data are suitable as final Reel source data",
        ],
    }
    (OUT / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": summary["status"], "kernel_checks": len(records),
                      "pic_steps": 55, "pic_particles": 16384,
                      "production_vela_dataset": False, "wall_seconds": wall}))


if __name__ == "__main__":
    main()
