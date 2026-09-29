#!/usr/bin/env python3
"""Centreline / wall Mach of a hydra-MI2 JPL run vs Chang (1980).

Kept next to the session report because $CLAUDE_JOB_DIR/tmp is cleaned between sessions and this
script was lost once already.

Conventions re-used from the Phase-8 analysis (recorded in auto-memory), NOT re-derived:
  * gas-field.tec vars: x y z rho(1) u v w p T g R mil kl   (g = gamma, R = 287.04)
    x,y,z are NODAL (I x J x K); the other 10 are CELLCENTERED ((I-1) x (J-1) x (K-1)).
    A uniform-shape parser silently produces garbage -- this bit me once (exit Mach 142738).
  * M = |v| / sqrt(g*R*T)
  * abscissa X = (x - x_throat)/D*, x_throat = 0.064207 m, r* = 0.02033740 m, D* = 2*r*
  * for dp = 20 micron the comparand is "Chang MPL Axis" (MPL = 20um, MPS = 2um)
  * j = 0 is the axis cell row, j = J-2 the wall row; each row has its OWN x (curved mesh) --
    never reuse the axis row's x for another row.

Validation (reproduce before trusting a new number):
  pre-dual-fix  OUTPUT.axisfix-dfce11c -> exit M 2.9085, axis RMS 0.4350, axis mean +0.3596
  post-dual-fix OUTPUT                 -> exit M 2.0890, axis RMS 0.0944, axis mean -0.0650
"""
import sys
import re
import numpy as np

X_THROAT = 0.064207
D_STAR = 2.0 * 0.02033740


def read_gas_field(path):
    lines = open(path).read().splitlines()
    hdr = lines[1]
    I = int(re.search(r'I\s*=\s*(\d+)', hdr).group(1))
    J = int(re.search(r'J\s*=\s*(\d+)', hdr).group(1))
    K = int(re.search(r'K\s*=\s*(\d+)', hdr).group(1))
    vals = np.array([float(t) for t in lines[2:] if t.strip()])

    n_nod = I * J * K
    n_cel = (I - 1) * (J - 1) * max(K - 1, 1)
    off = 0
    nodal = []
    for _ in range(3):                                   # x, y, z
        nodal.append(vals[off:off + n_nod].reshape(K, J, I))
        off += n_nod
    cell = []
    for _ in range(10):                                  # rho u v w p T g R mil kl
        cell.append(vals[off:off + n_cel].reshape(max(K - 1, 1), J - 1, I - 1))
        off += n_cel
    assert off == vals.size, f'leftover {vals.size - off} values'
    return nodal, cell, (I, J, K)


def centreline(path, jcell=0):
    (xn, yn, zn), c, (I, J, K) = read_gas_field(path)
    xc = 0.25 * (xn[0, jcell, :-1] + xn[0, jcell, 1:] +
                 xn[0, jcell + 1, :-1] + xn[0, jcell + 1, 1:])
    u, v, w = c[1][0, jcell], c[2][0, jcell], c[3][0, jcell]
    T, g, R = c[5][0, jcell], c[6][0, jcell], c[7][0, jcell]
    M = np.sqrt(u * u + v * v + w * w) / np.sqrt(g * R * T)
    return (xc - X_THROAT) / D_STAR, M


def chang_zone(path, name):
    lines = open(path).read().splitlines()
    for i, ln in enumerate(lines):
        if f'ZONE T="{name}"' in ln:
            n = int(re.search(r'I\s*=\s*(\d+)', lines[i + 1]).group(1))
            vals = []
            for ln2 in lines[i + 2:]:
                for t in ln2.split():
                    try:
                        vals.append(float(t))
                    except ValueError:
                        pass
                if len(vals) >= 2 * n:
                    break
            vals = np.array(vals[:2 * n])
            return vals[:n], vals[n:]                     # BLOCK packing: X then M
    raise KeyError(name)


def report(run, ref, label, zone='Chang MPL Axis', jcell=0):
    X, M = centreline(run, jcell)
    xr, mr = chang_zone(ref, zone)
    sel = X > 0.2
    d = M[sel] - np.interp(X[sel], xr, mr)
    rms, mean = float(np.sqrt((d ** 2).mean())), float(d.mean())
    print(f'{label}')
    print(f'   exit-plane Mach : {M[-1]:.4f}    (ref exit {mr[-1]:.3f})')
    print(f'   vs {zone}, X>0.2 : RMS {rms:.4f}   mean {mean:+.4f}')
    return float(M[-1]), rms, mean


if __name__ == '__main__':
    report(sys.argv[1], sys.argv[2],
           sys.argv[3] if len(sys.argv) > 3 else sys.argv[1],
           zone=sys.argv[4] if len(sys.argv) > 4 else 'Chang MPL Axis',
           jcell=int(sys.argv[5]) if len(sys.argv) > 5 else 0)
