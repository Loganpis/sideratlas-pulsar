"""Audit actual published magnetic-H EOS tables along a calculated column.

Trilinear interpolation in logB,logT,logR; no extrapolation.  This bounds
neither bound-bound/bound-free monochromatic opacity nor the optical flux.
Tables: Potekhin & Chabrier 2003/2004, Ioffe author distribution.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def read_table(path):
    blocks = {}
    current = None
    logb = None
    for line in gzip.decompress(path.read_bytes()).decode('ascii').splitlines():
        try:
            values = [float(item) for item in line.split()]
        except ValueError:
            continue
        if len(values) == 2 and 4 < values[0] < 9 and 10 < values[1] < 16:
            current, logb = values
            blocks[current] = []
        elif len(values) == 14 and current is not None:
            blocks[current].append(values)
    logt = np.array(sorted(blocks))
    data = np.array([blocks[t] for t in logt])
    logr = data[0, :, 0]
    assert np.allclose(data[:, :, 0], logr[None, :])
    assert np.all(np.diff(logt) > 0) and np.all(np.diff(logr) > 0)
    return logb, logt, logr, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('column', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    parent = Path(__file__).resolve().parent
    table_paths = [parent / f'hmag12_{i}.dat.gz' for i in [4, 5]]
    parsed = [read_table(path) for path in table_paths]
    b = np.log10(3e12)
    weight = (b - parsed[0][0]) / (parsed[1][0] - parsed[0][0])
    column = np.load(args.column)
    t = column['temperature_k']
    rho = column['density_g_cm3'] if 'density_g_cm3' in column else 1.67262192369e-24 * 1.3313e14 * column['mass_column_g_cm2'] / (2 * 1.380649e-16 * t)
    points = np.c_[np.log10(t), np.log10(rho / (t / 1e6)**3)]
    values = []
    for _, logt, logr, data in parsed:
        interp = RegularGridInterpolator((logt, logr), data[:, :, 1:], bounds_error=False, fill_value=np.nan)
        values.append(interp(points))
    value = (1 - weight) * values[0] + weight * values[1]
    valid = np.all(np.isfinite(value), axis=1)
    # data after removal of logR has original x(H) column8 at index7.
    rows = []
    for i in range(len(t)):
        row = dict(column_g_cm2=float(column['mass_column_g_cm2'][i]), temperature_k=float(t[i]),
                   density_g_cm3=float(rho[i]), within_table_domain=bool(valid[i]))
        if valid[i]:
            row.update(atomic_fraction=float(value[i, 7]), ground_state_atomic_fraction=float(value[i, 8]),
                       molecular_fraction_per_proton=float(value[i, 9]), perturbed_fraction=float(value[i, 10]),
                       eos_pressure_over_ideal_ionized=float(value[i, 1] / 2),
                       rosseland_parallel_cgs=float(10**value[i, 11]))
        rows.append(row)
    report = dict(
        input_column=str(args.column), input_sha256=hashlib.sha256(args.column.read_bytes()).hexdigest(),
        source_urls=[f'http://www.ioffe.ru/astro/NSG/Hmagnet/{path.name}' for path in table_paths],
        tables=[dict(name=path.name, bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in table_paths],
        interpolation='linear in table quantities over logT, logR, logB; no extrapolation',
        neutral_opacity_error_bounded=False, optical_intensity_error_bounded=False,
        valid_points=int(np.sum(valid)), total_points=len(t),
        maximum_atomic_fraction=float(np.max(value[valid, 7])) if np.any(valid) else None,
        maximum_molecular_fraction_per_proton=float(np.max(value[valid, 9])) if np.any(valid) else None,
        rows=rows,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
