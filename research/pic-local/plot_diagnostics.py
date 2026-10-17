#!/usr/bin/env python3
"""Plot saved audit measurements; this figure is not a photon image."""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audit', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = json.loads(args.audit.read_text())
    states = report['snapshots']
    t = [s['position_and_E_time'] for s in states]
    plt.rcParams.update({'font.family': 'DejaVu Serif', 'mathtext.fontset': 'cm',
                         'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    ax = axes[0, 0]
    ax.semilogy(t, [s['divB_rms_over_rms_B_per_rlc'] for s in states], '.-', color='#156082')
    ax.set(title='Discrete magnetic divergence in the interior shell',
           ylabel=r'$\mathrm{rms}(\nabla\!\cdot B)/\mathrm{rms}(B)\ [R_{LC}^{-1}]$')
    ax = axes[0, 1]
    ratios = np.array([s['gauss_relative_l2_symmetric'] if
                      s['gauss_relative_l2_symmetric'] is not None else np.nan for s in states])
    ax.plot(t, 100*ratios, '.-', color='#9a3412')
    ax.set(title='Electric constraint residual in the interior shell',
           ylabel=r'$100\,\Vert\nabla\!\cdot E-\rho\Vert/(\Vert\nabla\!\cdot E\Vert+\Vert\rho\Vert)\ [\%]$')
    ax.text(.03, .97, 'Undefined near a zero reference; no convergence claim',
            va='top', transform=ax.transAxes, fontsize=8)
    ax = axes[1, 0]
    ax.plot(t, [s['species_counts']['1'] for s in states], '-', label='Electrons', color='#156082')
    ax.plot(t, [s['species_counts']['2'] for s in states], '--', label='Positrons', color='#9a3412')
    ax.set(title='Live numerical particles', ylabel='Macroparticles')
    ax.legend(frameon=False)
    ax = axes[1, 1]
    ax.plot(t, [s['gamma_max'] for s in states], '-', label='Maximum', color='#156082')
    ax.plot(t, [s['gamma_mean'] for s in states], '--', label='Mean', color='#9a3412')
    ax.set(title='Particle Lorentz factors in the reduced model', ylabel=r'$\gamma$')
    ax.legend(frameon=False)
    for ax in axes.flat:
        ax.set_xlabel(r'Completed position/electric-field time $[R_{LC}/c]$')
        ax.grid(alpha=.18)
    fig.suptitle('Local 3D PIC prototype — numerical checks, not a Vela optical image', fontsize=14)
    fig.supxlabel('48³ mesh · prescribed pair supply · enlarged numerical star · magnetic field and momenta at t − Δt/2', fontsize=9)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
