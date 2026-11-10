#!/usr/bin/env python3
"""The JPL nozzle runs against the curves of Chang, Hourng & Chien (1996), Fig. 11: gas Mach number on the axis
and on the wall. REPORT MODE: every row is printed, none is judged, the exit status is 0.

usage: jpl_chang_gate.py [--ref JPL_chang.dat] [--xmin 0.2] [<label>=<run dir or gas-field.tec> ...]
       without a run: OUTPUT/gas-field.tec of the three JPL cases of this tree (IGLOO, ICE, both), those that exist.

Definitions (those of JPL-Lagrangian-20micron/jpl_compare.py; the reader is tec_cc.read_cc):
  X = (x - x_throat)/D*, x of a cell = mean of the four nodes of its k = 0 face; axis = first cell row, wall = last;
  M = |u|/sqrt(g R T). A reference zone is interpolated onto the run's cells; a row uses the cells with
  xmin < X <= the last station of the reference zone. The exit Mach number is taken at that last station.
Reference zones: OG = gas alone, MPS = 1 micron radius, MPL = 10 micron radius (the case the runs reproduce)."""
import hashlib
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tec_cc import read_cc  # noqa: E402

X_THROAT, D_STAR = 0.064207, 2.0 * 0.02033740
ROWS = (('axis', 0, 'Axis'), ('wall', -1, 'Wall'))
DEFAULT = (('IGLOO', 'JPL-Lagrangian-20micron'), ('ICE', 'JPL-Eulerian-20micron'), ('both', 'JPL-3solver-20micron'))


def chang_zone(ref, name):
    """(X, M) of one zone of the reference file."""
    lines = open(ref).read().splitlines()
    for i, line in enumerate(lines):
        if 'ZONE T="%s"' % name in line:
            n = int(re.search(r'I\s*=\s*(\d+)', lines[i + 1]).group(1))
            vals = []
            for l2 in lines[i + 2:]:
                try:
                    vals += [float(t) for t in l2.split()]
                except ValueError:
                    continue
                if len(vals) >= 2 * n:
                    break
            return np.array(vals[:n]), np.array(vals[n:2 * n])
    raise KeyError('%s: no zone "%s"' % (ref, name))


def mach_row(field, jcell):
    """(X, M) along one cell row of a gas field read by read_cc."""
    x = field['x']
    j = jcell % field['T'].shape[1]
    xc = 0.25 * (x[0, j, :-1] + x[0, j, 1:] + x[0, j + 1, :-1] + x[0, j + 1, 1:])
    u, v, w = (field[k][0, j, :] for k in ('u', 'v', 'w'))
    mach = np.sqrt(u * u + v * v + w * w) / np.sqrt(field['g'][0, j, :] * field['R'][0, j, :] * field['T'][0, j, :])
    return (xc - X_THROAT) / D_STAR, mach


def against(X, M, zone, xmin):
    """Run minus reference zone on the cells the zone covers downstream of xmin."""
    xr, mr = zone
    s = (X > xmin) & (X <= xr[-1])
    d = M[s] - np.interp(X[s], xr, mr)
    return dict(n=int(s.sum()), x0=X[s].min(), x1=X[s].max(), rms=np.sqrt((d ** 2).mean()), mean=d.mean(), std=d.std(),
                exit=float(np.interp(xr[-1], X, M)), refexit=mr[-1])


def main():
    args, ref, xmin, runs = sys.argv[1:], None, 0.2, []
    while args:
        a = args.pop(0)
        if a == '--ref':
            ref = args.pop(0)
        elif a == '--xmin':
            xmin = float(args.pop(0))
        elif '=' in a:
            runs.append(tuple(a.split('=', 1)))
        else:
            sys.exit(__doc__)
    test = os.path.dirname(HERE)
    if ref is None:
        ref = os.path.join(test, 'JPL-Lagrangian-20micron', 'JPL_chang.dat')
    if not runs:
        runs = [(lab, os.path.join(test, case)) for lab, case in DEFAULT]
        runs = [(lab, p) for lab, p in runs if os.path.isfile(os.path.join(p, 'OUTPUT', 'gas-field.tec'))]
    print('== JPL nozzle against Chang, Hourng & Chien (1996), Fig. 11 -- report mode: nothing is judged')
    print('reference %s (md5 %s), rows on %g < X <= the last reference station'
          % (ref, hashlib.md5(open(ref, 'rb').read()).hexdigest()[:8], xmin))
    zones = {(kind, Z): chang_zone(ref, 'Chang %s %s' % (kind, Z)) for kind in ('MPL', 'MPS', 'OG') for _, _, Z in ROWS}
    rows = {}
    for lab, path in runs:
        if os.path.isdir(path):
            path = os.path.join(path, 'OUTPUT', 'gas-field.tec') if os.path.isdir(os.path.join(path, 'OUTPUT')) \
                else os.path.join(path, 'gas-field.tec')
        field = read_cc(path)
        print('%s: %s' % (lab, path))
        for row, j, Z in ROWS:
            X, M = mach_row(field, j)
            rows[lab, row] = (X, M)
            a, b, o = (against(X, M, zones[kind, Z], xmin) for kind in ('MPL', 'MPS', 'OG'))
            print('  %s: %d cells, X %.3f..%.3f | run - MPL: rms %.4f mean %+.4f std %.4f | exit M %.4f, MPL %.3f (%+.4f)'
                  ' | rms against MPS %.4f, OG %.4f | margin %.2f'
                  % (row, a['n'], a['x0'], a['x1'], a['rms'], a['mean'], a['std'], a['exit'], a['refexit'],
                     a['exit'] - a['refexit'], b['rms'], o['rms'], min(b['rms'], o['rms']) / a['rms']))
    labs = [lab for lab, _ in runs]
    for i, A in enumerate(labs):
        for B in labs[i + 1:]:
            out = []
            for row, _, _ in ROWS:
                (XA, MA), (XB, MB) = rows[A, row], rows[B, row]
                if XA.shape != XB.shape or not np.allclose(XA, XB, rtol=0.0, atol=1e-12):
                    out.append('%s: different cells' % row)
                    continue
                d = (MA - MB)[XA > xmin]
                out.append('%s rms %.4f max %.4f' % (row, np.sqrt((d ** 2).mean()), np.abs(d).max()))
            print('gas Mach, %s - %s (X > %g): %s' % (A, B, xmin, ' | '.join(out)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
