#!/usr/bin/env python3
"""Read saved staggered fields and audit the numerical rotating magnetic source.

This measures the current prototype; it neither changes nor reruns it. Dipole
fits describe the seed-field family and are not an inference of Vela's field.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


OFFSETS_B = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])


def potential(points, moment_axis, radius):
    moment = np.eye(3)[moment_axis]
    r2 = np.sum(points * points, axis=-1)
    outside = r2 >= radius ** 2
    f = np.empty_like(r2)
    f[outside] = r2[outside] ** -1.5
    f[~outside] = (2.5 - 1.5 * r2[~outside] / radius ** 2) / radius ** 3
    return np.cross(moment, points) * f[..., None]


def discrete_curl_matrix(points, component, radius, dx):
    # A_j(m) is linear in each component of m. A centered derivative of the
    # edge potential is exactly the native compatible Yee initialization.
    a, b = ((2, 1), (0, 2), (1, 0))[component]
    result = np.empty((len(points), 3))
    for m in range(3):
        ep = np.zeros(3); ep[b] = dx / 2
        eq = np.zeros(3); eq[a] = dx / 2
        result[:, m] = ((potential(points + ep, m, radius)[:, a] -
                          potential(points - ep, m, radius)[:, a]) / dx -
                         (potential(points + eq, m, radius)[:, b] -
                          potential(points - eq, m, radius)[:, b]) / dx)
    return result


def sample_staggered_b(field, points, origin, dx):
    # field is (n,n,n,Bcomponent[,momentcomponent]). Never treat components as
    # collocated and never extrapolate through the computational boundary.
    extra = field.shape[4:]
    result = np.empty((len(points), 3) + extra)
    for c in range(3):
        coordinate = (points - origin) / dx - OFFSETS_B[c]
        left = np.floor(coordinate).astype(int)
        fraction = coordinate - left
        if np.any(left < 0) or np.any(left >= np.array(field.shape[:3]) - 1):
            raise ValueError("Magnetic sampling would extrapolate")
        value = np.zeros((len(points),) + extra)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    weight = ((fraction[:, 0] if i else 1 - fraction[:, 0]) *
                              (fraction[:, 1] if j else 1 - fraction[:, 1]) *
                              (fraction[:, 2] if k else 1 - fraction[:, 2]))
                    weight = weight.reshape((len(points),) + (1,) * len(extra))
                    value += weight * field[left[:, 0] + i, left[:, 1] + j,
                                             left[:, 2] + k, c]
        result[:, c] = value
    return result


def spin_phase(t, omega, spinup):
    s = np.clip(t / spinup, 0, 1)
    return omega * (spinup * (s ** 3 - .5 * s ** 4) + max(0, t - spinup))


def fit(design, values, expected, reference):
    moment = np.linalg.lstsq(design, values, rcond=None)[0]
    norm = np.linalg.norm(moment)
    predicted = design @ moment
    phase = np.arctan2(moment[1], moment[0])
    expected_phase = np.arctan2(expected[1], expected[0])
    phase_error = np.arctan2(np.sin(phase - expected_phase),
                             np.cos(phase - expected_phase))
    cosine = np.dot(moment, expected) / (norm * np.linalg.norm(expected))
    return {
        "moment_B0_RLC3": moment.tolist(), "magnitude": float(norm),
        "magnitude_relative_to_initial_seed": float(norm / reference),
        "inclination_deg": float(np.degrees(np.arctan2(np.hypot(*moment[:2]), moment[2]))),
        "phase_deg_mod_360": float(np.degrees(phase) % 360),
        "phase_error_deg_vs_rigid_rotation_at_B_time": float(np.degrees(phase_error)),
        "orientation_error_deg_vs_rigid_rotation_at_B_time": float(np.degrees(np.arccos(np.clip(cosine, -1, 1)))),
        "relative_fit_residual_L2": float(np.linalg.norm(values - predicted) /
                                           max(np.linalg.norm(values), 1e-300)),
        "design_condition_number": float(np.linalg.cond(design)),
    }


def audit(native, sphere_samples=4096):
    files = sorted(native.glob("state-*.json"))
    if not files:
        raise ValueError("No native states found")
    first = json.loads(files[0].read_text())
    n = first["active_resolution"][0]
    dx, R = first["dx"], first["stellar_radius"]
    origin = np.asarray(first["origin"])
    g = first["ghost_cells_each_side"]
    grid = np.stack(np.meshgrid(np.arange(n), np.arange(n), np.arange(n),
                                indexing="ij"), axis=-1)
    positions, bases, body_masks = [], [], []
    for c in range(3):
        points = origin + (grid + OFFSETS_B[c]) * dx
        positions.append(points)
        bases.append(discrete_curl_matrix(points.reshape(-1, 3), c, R, dx).reshape(n, n, n, 3))
        body_masks.append(np.linalg.norm(points, axis=-1) < first.get("anchor_guaranteed_B_radius", R-dx))
    # Last two axes are the staggered field component and source moment axis.
    basis_grid = np.stack(bases, axis=-2)
    body_design = np.concatenate([bases[c][body_masks[c]] for c in range(3)])
    j = np.arange(sphere_samples)
    nz = 1 - 2 * (j + .5) / sphere_samples
    phi = 2.399963229728653 * j
    directions = np.column_stack((np.sqrt(1 - nz ** 2) * np.cos(phi),
                                  np.sqrt(1 - nz ** 2) * np.sin(phi), nz))
    shell_radii = (R - dx, R, R + 2 * dx)
    shell_matrices = []
    for radius in shell_radii:
        basis = sample_staggered_b(basis_grid, directions * radius, origin, dx)
        radial = np.einsum("nc,nck->nk", directions, basis)
        tangent = basis - directions[..., None] * radial[:, None, :]
        shell_matrices.append((radius, radial, tangent.reshape(-1, 3)))

    seed_magnitude = R ** 3 / 2  # supplied polar_field=1, independently checked below
    chi = np.radians(first["inclination_deg"])
    rows = []
    for path in files:
        s = json.loads(path.read_text())
        for key in ("field_shape", "origin", "dx", "stellar_radius", "inclination_deg"):
            if s[key] != first[key]:
                raise ValueError(f"State geometry changed: {key}")
        raw = np.fromfile(native / s["fields_file"], dtype="<f8").reshape(s["field_shape"])
        B = raw[g:g+n, g:g+n, g:g+n, 3:6] * dx
        if not np.all(np.isfinite(B)):
            raise ValueError("Nonfinite magnetic state")
        tB = s["momentum_B_and_J_time"]
        phase = spin_phase(tB, s["final_omega"], s["spinup_time"])
        expected = seed_magnitude * np.array([np.sin(chi) * np.cos(phase),
                                               np.sin(chi) * np.sin(phase), np.cos(chi)])
        values = np.concatenate([B[..., c][body_masks[c]] for c in range(3)])
        body = fit(body_design, values, expected, seed_magnitude)
        shells = []
        for radius, radial_design, tangent_design in shell_matrices:
            field = sample_staggered_b(B, directions * radius, origin, dx)
            radial = np.sum(field * directions, axis=-1)
            tangent = field - radial[:, None] * directions
            shells.append({"radius_RLC": float(radius),
                            "region": "interior" if radius < R else "stellar_surface" if radius == R else "exterior_effective_moment",
                            "radial": fit(radial_design, radial, expected, seed_magnitude),
                            "tangential": fit(tangent_design, tangent.reshape(-1), expected, seed_magnitude)})
        rows.append({"completed_step": s["completed_step"],
                     "position_E_time": s["position_and_E_time"],
                     "magnetic_time": tB, "expected_phase_at_B_time_rad": float(phase),
                     "interior_fit": body, "spherical_shell_fits": shells,
                     "fields_sha256": hashlib.sha256((native / s["fields_file"]).read_bytes()).hexdigest()})
    initial = rows[0]["interior_fit"]
    if abs(initial["magnitude"] / seed_magnitude - 1) > 1e-10 or initial["relative_fit_residual_L2"] > 1e-10:
        raise ValueError("Initial state does not match the independently reconstructed finite dipole")
    body_rows = [r["interior_fit"] for r in rows]
    final = body_rows[-1]
    steady_rows = [r for r in rows if r["magnetic_time"] >= first["spinup_time"]]
    phase_rate = None
    if len(steady_rows) >= 2:
        times = np.array([r["magnetic_time"] for r in steady_rows])
        phases = np.array([r["expected_phase_at_B_time_rad"] +
                           np.radians(r["interior_fit"]["phase_error_deg_vs_rigid_rotation_at_B_time"])
                           for r in steady_rows])
        phase_rate = float(np.polyfit(times, phases, 1)[0])
    pre_boundary = [r for r in rows if r["magnetic_time"] < 2.0 - R]
    return {
        "scope": "independent numerical-source audit of saved scaled kinetic prototype; no Vela inference",
        "method": "linear least-squares finite discrete-curl dipole family at actual staggered B positions; shell radial/tangential fits use separately interpolated component grids",
        "raw_to_tetrad": "B_normalized = dx * B_U; no collocated-component assumption",
        "body_mask": "each magnetic component sampled at its own Yee location inside the guaranteed anchored B core, r < anchor_core_fraction*Rstar-dx/2",
        "body_data_rows": int(len(body_design)), "sphere_samples": sphere_samples,
        "seed_magnitude_B0_RLC3": seed_magnitude, "initial_inclination_deg": float(np.degrees(chi)),
        "summary": {
            "final_interior_magnitude_relative_to_initial": final["magnitude_relative_to_initial_seed"],
            "final_interior_inclination_deg": final["inclination_deg"],
            "final_interior_inclination_change_deg": final["inclination_deg"] - np.degrees(chi),
            "final_interior_phase_error_deg": final["phase_error_deg_vs_rigid_rotation_at_B_time"],
            "final_interior_fit_residual_L2": final["relative_fit_residual_L2"],
            "max_abs_interior_magnitude_change_fraction": float(max(abs(x["magnitude_relative_to_initial_seed"]-1) for x in body_rows)),
            "max_abs_interior_inclination_change_deg": float(max(abs(x["inclination_deg"]-np.degrees(chi)) for x in body_rows)),
            "max_abs_interior_phase_error_deg": float(max(abs(x["phase_error_deg_vs_rigid_rotation_at_B_time"]) for x in body_rows)),
            "interior_moment_phase_rate_after_spinup": phase_rate,
            "prescribed_final_rotation_rate": first["final_omega"],
            "before_nominal_boundary_light_time_state": {"magnetic_time": pre_boundary[-1]["magnetic_time"],
                                            "fit": pre_boundary[-1]["interior_fit"]} if pre_boundary else None,
        },
        "boundary_and_source_risk": {
            "star_radius_RLC": R, "star_radius_cells": R/dx,
            "box_half_extent_RLC": 2.5, "zero_match_layer_width_RLC": .5,
            "zero_match_layer_inner_face_RLC": 2.0,
            "nominal_continuum_initial_boundary_to_surface_light_time": 2.0-R,
            "nominal_continuum_source_wave_roundtrip_light_time": 2*(2.0-R),
            "final_time": rows[-1]["position_E_time"],
            "causally_isolated_for_entire_run": False,
            "boundary_or_resolution_sensitivity_measured": False,
            "interpretation": "The source-repair experiment prescribes a rotating interior magnetization through a compatible compact interior electric field constructed from a rotating edge potential and a vertex scalar potential. No B reset is applied; exterior fields remain native PIC evolved. The guaranteed core fit tests the imposed maintained magnetization. Surface/exterior behavior still includes source staircasing, the blend layer, loading and possible later boundary influence. No single-run fit separates those causes. Exterior fitted moments also include genuine induction/wave/plasma structure and are not conserved source moments.",
        },
        "limits": [
            "B and fitted magnetic phase are evaluated at t-dt/2, not the E/position timestamp",
            "The fit seed has the model's finite interior current profile, not a uniform-field or singular-dipole replacement",
            "Tangential fits near the stellar surface can have large residuals and should not be mistaken for direct magnetization measurements",
            "No changed-resolution, source-radius, spinup, loading or boundary-control run was performed in this audit",
            "Boundary times are continuum c=1 light-travel estimates, not a proof of exact finite-difference isolation; numerical stencil tails or constraint errors require explicit boundary controls",
            "A finite-grid seed-field conductor is not a validated fixed-obliquity magnetized neutron star atmosphere",
        ],
        "states": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native_directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.native_directory)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
