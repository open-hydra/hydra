#!/usr/bin/env python3
"""Axis and wall Mach number of a JPL nozzle run vs Chang (1980).

Usage:
    jpl_compare.py [FIELD ...] [--ref JPL_chang.dat] [--xmin 0.2] [--plot]

FIELD defaults to OUTPUT/gas-field.tec next to this script; with several fields, each one is
reported. For every field, the axis (first cell row) is compared with "Chang MPL Axis" and the wall
(last cell row) with "Chang MPL Wall" (MPL = the 20 micron particles), over X > xmin:

    M = |u| / sqrt(g R T)          X = (x - x_throat) / D*

Each cell row has its own x (the mesh is curved), taken at the cell centre of the k = 0 face.
"""
import argparse
import functools
import re
import sys
from pathlib import Path

import numpy as np

# Find ORION's Python reader even when it is not installed in the active environment
HERE = Path(__file__).resolve().parent
for parent in [HERE, *HERE.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break

try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "ORION package not found. Run 'conda activate base' or install ORION in the current Python environment."
    )

X_THROAT = 0.064207              # throat abscissa [m]
D_STAR   = 2.0 * 0.02033740      # throat diameter [m]
ROWS = {'axis': (0, 'Chang MPL Axis'), 'wall': (-1, 'Chang MPL Wall')}


@functools.lru_cache(maxsize=None)
def load_field(path):
    """Block 1 of a MOSE field: nodal (x, y) and the cell-centred variables, by name."""
    xb, yb, _, vb, names = read_TEC(str(path))
    var = dict(zip(names[3:], vb[0]))
    return xb[0], yb[0], var


def mach_row(path, jcell):
    """X and Mach number along cell row jcell (0 = axis, -1 = wall)."""
    xn, _, var = load_field(path)
    j = jcell % var['T'].shape[1]
    xc = 0.25 * (xn[:-1, j, 0] + xn[1:, j, 0] + xn[:-1, j + 1, 0] + xn[1:, j + 1, 0])
    u, v, w = (var[k][:, j, 0] for k in ('u', 'v', 'w'))
    M = np.sqrt(u * u + v * v + w * w) / np.sqrt(var['g'][:, j, 0] * var['R'][:, j, 0] * var['T'][:, j, 0])
    return (xc - X_THROAT) / D_STAR, M


@functools.lru_cache(maxsize=None)
def chang_zone(path, name):
    """(X, M) of one zone of the reference file (1-D zones, BLOCK packing: all X, then all M)."""
    lines = Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        if f'ZONE T="{name}"' in line:
            n = int(re.search(r'I\s*=\s*(\d+)', lines[i + 1]).group(1))
            vals = []
            for line2 in lines[i + 2:]:
                try:
                    vals += [float(t) for t in line2.split()]
                except ValueError:
                    continue
                if len(vals) >= 2 * n:
                    break
            return np.array(vals[:n]), np.array(vals[n:2 * n])
    raise KeyError(f'zone "{name}" not in {path}')


def compare(field, ref, zone, jcell, xmin=0.2):
    """Exit Mach, RMS and mean of (run - reference) over X > xmin."""
    X, M = mach_row(field, jcell)
    xr, mr = chang_zone(ref, zone)
    sel = X > xmin
    d = M[sel] - np.interp(X[sel], xr, mr)
    return float(M[-1]), float(np.sqrt((d ** 2).mean())), float(d.mean()), float(mr[-1])


def report(run, ref, label, zone='Chang MPL Axis', jcell=0, xmin=0.2):
    """Print and return (exit Mach, RMS, mean) for one row (kept for compare_3solver.py)."""
    exit_m, rms, mean, ref_exit = compare(run, ref, zone, jcell, xmin)
    print(f'{label}')
    print(f'   exit-plane Mach : {exit_m:.4f}    (ref exit {ref_exit:.3f})')
    print(f'   vs {zone}, X>{xmin:g} : RMS {rms:.4f}   mean {mean:+.4f}')
    return exit_m, rms, mean


def plot(fields, ref):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5))
    for row, (jcell, zone) in ROWS.items():
        xr, mr = chang_zone(ref, zone)
        ax.plot(xr, mr, 'o', markersize=4, alpha=0.7, label=zone)
        for f in fields:
            X, M = mach_row(f, jcell)
            ax.plot(X, M, linewidth=2, label=f'{f} ({row})' if len(fields) > 1 else row)
    ax.set_xlabel('(x - x*) / D*')
    ax.set_ylabel('Mach')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('fields', nargs='*', default=[str(HERE / 'OUTPUT' / 'gas-field.tec')],
                        help='MOSE gas field(s) (default: OUTPUT/gas-field.tec next to this script)')
    parser.add_argument('--ref', default=str(HERE / 'JPL_chang.dat'), help='Chang (1980) reference curves')
    parser.add_argument('--xmin', type=float, default=0.2, help='compare over X > xmin (default 0.2)')
    parser.add_argument('--plot', action='store_true', help='plot Mach vs X with the reference')
    args = parser.parse_args()

    for f in args.fields:
        if not Path(f).is_file():
            sys.exit(f'{f}: not found (is sol-diter set in [MOSE-IO]?)')

    for f in args.fields:
        print(f'{f}  (X > {args.xmin:g})')
        print(f'   {"row":<5} {"exit M":>7} {"ref":>6} {"RMS":>7} {"mean":>8}')
        for row, (jcell, zone) in ROWS.items():
            exit_m, rms, mean, ref_exit = compare(f, args.ref, zone, jcell, args.xmin)
            print(f'   {row:<5} {exit_m:7.4f} {ref_exit:6.3f} {rms:7.4f} {mean:+8.4f}')

    if args.plot:
        plot(args.fields, args.ref)


if __name__ == '__main__':
    main()
