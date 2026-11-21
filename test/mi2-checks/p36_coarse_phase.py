#!/usr/bin/env python3
"""The gas Mach number of a coarse grid level against the finest grid's, on the JPL nozzle (gate p36, long leg).

usage: p36_coarse_phase.py <coarse field> <fine field>
         Tecplot BLOCK fields of MOSE (nodal x y z, cell-centred u v w T g R ...), the coarse one on cells of 2 x 2 fine
         cells. Prints, over the coarse cells: the largest |M_coarse - M_fine|, M_fine the mean of the four sub-cells'
         Mach numbers, on the exit column (the last coarse column) and over the cells downstream of the throat
         (X = (x - x_throat)/D* > 0.2, x of a cell the mean of its k = 0 face's nodes), and the exit Mach number of the
         two fields on the axis (first row) and at the wall (last row).
       p36_coarse_phase.py --same <field A> <field B>
         Prints, per cell-centred variable, max |A - B| over the variable's scale in A, the larger of its range and of
         its largest magnitude (the largest over the variables last)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tec_cc import read_cc  # noqa: E402

X_THROAT, D_STAR = 0.064207, 2.0 * 0.02033740


def mach(d):
    return np.sqrt(d['u'] ** 2 + d['v'] ** 2 + d['w'] ** 2) / np.sqrt(d['g'] * d['R'] * d['T'])


def cell_x(d):
    x = d['x'][0]                                     # (J, I) nodes of the k = 0 plane
    return 0.25 * (x[:-1, :-1] + x[1:, :-1] + x[:-1, 1:] + x[1:, 1:])


def compare(coarse, fine):
    c, f = read_cc(coarse), read_cc(fine)
    mc, mf = mach(c)[0], mach(f)[0]                   # (J-1, I-1) cells
    jc, ic = mc.shape
    assert mf.shape == (2 * jc, 2 * ic), 'the fine field is not 2 x 2 cells per coarse cell'
    mr = 0.25 * (mf[0::2, 0::2] + mf[1::2, 0::2] + mf[0::2, 1::2] + mf[1::2, 1::2])
    X = (cell_x(c) - X_THROAT) / D_STAR
    d = np.abs(mc - mr)
    down = X > 0.2
    print('exit column: max |M_coarse - M_fine| %.6e' % d[:, -1].max())
    print('X > 0.2: max |M_coarse - M_fine| %.6e over %d cells' % (d[down].max(), down.sum()))
    print('exit Mach, axis: coarse %.6f fine %.6f | wall: coarse %.6f fine %.6f' % (mc[0, -1], mf[0, -1], mc[-1, -1],
                                                                                   mf[-1, -1]))


def same(a, b):
    A, B = read_cc(a), read_cc(b)
    worst = 0.0
    for n in A:
        if n in ('x', 'y', 'z') or n not in B:
            continue
        scale = max(float(A[n].max() - A[n].min()), float(np.abs(A[n]).max()))
        dd = float(np.abs(A[n] - B[n]).max())
        r = dd / scale if scale > 0 else dd
        worst = max(worst, r)
        print('%-8s max|d| %.6e scale %.6e max|d|/scale %.6e' % (n, dd, scale, r))
    print('largest max|d|/scale %.6e' % worst)


if __name__ == '__main__':
    if sys.argv[1] == '--same':
        same(sys.argv[2], sys.argv[3])
    else:
        compare(sys.argv[1], sys.argv[2])
