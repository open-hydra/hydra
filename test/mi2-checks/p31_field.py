#!/usr/bin/env python3
"""p31_field.py peak <field.tec> <var> | small <field.tec> <family> <rho_mat> <d>: of an OUTPUT field, the largest
value of a cell variable (peak), or the number of cells family <family> occupies (bulk density above ICE's empty
1e-6 kg/m3) with drops of diameter (6 rho / (pi n rho_mat))^(1/3) at or below <d>, and their smallest diameter
(small). Prints one line; reads only."""
import math
import re
import sys


def read(path):
    lines = open(path).read().split('\n')
    names = lines[0].split('=', 1)[1].split()
    I, J, K = (int(x) for x in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', ' '.join(lines[:2])).groups())
    vals = ' '.join(lines[2:]).split()
    cols, off = {}, 0
    for i, n in enumerate(names):
        k = I * J * K if i < 3 else (I - 1) * (J - 1) * max(K - 1, 1)
        cols[n] = [float(x) for x in vals[off:off + k]]
        off += k
    return cols


if sys.argv[1] == 'peak':
    print('%.6e' % max(read(sys.argv[2])[sys.argv[3]]))
elif sys.argv[1] == 'small':
    c, f = read(sys.argv[2]), sys.argv[3]
    rho_mat, dmax = float(sys.argv[4]), float(sys.argv[5])
    ds = [(6.0 * r / (math.pi * n * rho_mat)) ** (1.0 / 3.0)
          for r, n in zip(c['rho_p' + f], c['n_p' + f]) if r > 1e-6 and n > 0.0]
    small = [d for d in ds if d <= dmax]
    print('%d %.6e' % (len(small), min(small) if small else float('nan')))
else:
    sys.exit('p31_field.py: peak or small')
