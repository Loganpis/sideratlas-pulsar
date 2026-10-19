#!/usr/bin/env python3
"""Compare the coupled startup with its otherwise identical no-deposit control.

This is a causal feedback check, not a convergence or Vela validation test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def state(directory, step):
    metadata = json.loads((directory / f"state-{step:06d}.json").read_text())
    array = np.fromfile(directory / metadata["fields_file"], dtype="<f8")
    return metadata, array.reshape(metadata["field_shape"])


def compare(feedback, control, step):
    s, full = state(feedback, step)
    c, empty = state(control, step)
    for key in ("completed_step", "position_and_E_time", "dt", "field_shape",
                "ghost_cells_each_side", "active_resolution", "origin", "dx", "B0"):
        if s[key] != c[key]:
            raise ValueError(f"Control geometry/timing mismatch: {key}")
    g, n, h = s["ghost_cells_each_side"], s["active_resolution"][0], s["dx"]
    active = (slice(g, g + n),) * 3
    x = s["origin"][0] + np.arange(n) * h
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    radius = np.sqrt(xx ** 2 + yy ** 2 + zz ** 2)
    mask = (radius > s["stellar_radius"] + 3 * h) & (radius < 1.8)
    f, c = full[active], empty[active]
    initial_full = state(feedback, 0)[1]
    initial_empty = state(control, 0)[1]
    csv_name = "state-000000-particles.csv"
    csv_hashes = [hashlib.sha256((d / csv_name).read_bytes()).hexdigest()
                  for d in (feedback, control)]
    if not np.array_equal(initial_full, initial_empty) or csv_hashes[0] != csv_hashes[1]:
        raise ValueError("Control and coupled simulation have different initial conditions")
    result = {
        "scope": "causal field-feedback check for the scaled local PIC startup; not Vela validation",
        "time": s["position_and_E_time"], "completed_step": step,
        "initial_field_arrays_bit_identical": True,
        "initial_particle_csv_sha256": csv_hashes,
        "changed_input": "algorithms.deposit.enable=false in the separate control; geometry/injection/solver/source otherwise identical",
        "mask": "Rstar+3*dx < vertex-radius < 1.8",
        "measured_vertices": int(mask.sum()),
        "claim": "Disabling particle current deposition changes the electromagnetic fields; the enabled run is not a prescribed-field test-particle replay",
    }
    for name, components in (("E", slice(0, 3)), ("B", slice(3, 6))):
        a, b = f[..., components][mask] * h, c[..., components][mask] * h
        result[f"{name}_difference_RMS_B0"] = float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=-1))))
        result[f"{name}_difference_relative_L2"] = float(np.linalg.norm(a - b) / np.linalg.norm(a))
    result["current_feedback_run_L2_n0"] = float(np.linalg.norm(f[..., 6:9][mask] * h))
    result["current_control_L2_n0"] = float(np.linalg.norm(c[..., 6:9][mask] * h))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("feedback_native", type=Path)
    parser.add_argument("control_native", type=Path)
    parser.add_argument("--step", type=int, default=38)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(args.feedback_native, args.control_native, args.step)
    serialized = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
