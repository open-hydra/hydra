#!/usr/bin/env python3
"""Rewrite the two columns of INPUT/part-ic.tec that ATLAS ICB cannot set for this case.

ICB writes the condensed velocity equal to the gas velocity and the number density with pi
rounded to 3.14 (ATLAS src/ICB/builder_dp.f90). This case wants the velocity of [ice-box] up
and an exact number density, so that ICE's particle radius is exactly [ICB-Block1] dp / 2 and
the oracle of check.py follows from input.ini alone. Run it after every regeneration:

    make_case.py [case-dir]        default: the directory of this script

It reads <case-dir>/input.ini ([GPB-Phase2] rho, [ICB-Block1] dp, [ice-box] up), rewrites
the up_11 and np_11 columns of <case-dir>/INPUT/part-ic.tec in ATLAS's own number format
(every other line is kept byte for byte), and prints what it wrote. Idempotent."""
import configparser
import math
import os
import re
import sys


def read_ini(path):
    cfg = configparser.ConfigParser(inline_comment_prefixes=(';',), interpolation=None, strict=False)
    cfg.optionxform = str
    cfg.read(path)
    return cfg


def fmt(v):
    """ATLAS ICB's number format: ' 0.ddddddddddddddd' + 'E' + signed 3-digit exponent."""
    if v == 0.0:
        return ' 0.000000000000000E+000'
    s = '%.14E' % abs(v)                       # d.dddddddddddddd E+ee  (15 significant digits)
    mant, ex = s.split('E')
    digits = mant[0] + mant[2:]
    return ('-' if v < 0 else ' ') + '0.' + digits + 'E%+04d' % (int(ex) + 1)


def main():
    case = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__))
    ini = os.path.join(case, 'input.ini')
    tec = os.path.join(case, 'INPUT', 'part-ic.tec')
    cfg = read_ini(ini)
    rho_mat = float(cfg['GPB-Phase2']['rho'])
    dp = float(cfg['ICB-Block1']['dp'])
    up = float(cfg['ice-box']['up'])

    lines = open(tec).read().split('\n')
    names = [t.strip('"') for t in re.findall(r'"[^"]+"|\S+', lines[0].split('=', 1)[1])]
    I, J, K = (int(x) for x in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', lines[1]).groups())
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    start, first = {}, 2
    for n in names:
        cnt = nn if n in ('x', 'y', 'z') else nc
        start[n] = first
        first += cnt
    for key in ('rp_11', 'up_11', 'np_11'):
        if key not in start:
            sys.exit('[make_case] %s: no column %s (columns: %s)' % (tec, key, ' '.join(names)))
    if len(lines) < first:
        sys.exit('[make_case] %s: %d lines, %d expected' % (tec, len(lines), first))

    rho_p = [float(lines[start['rp_11'] + c]) for c in range(nc)]
    rp = 0.5 * dp
    vol = 4.0 / 3.0 * math.pi * rp ** 3
    n_p = [r / (rho_mat * vol) for r in rho_p]
    for c in range(nc):
        lines[start['up_11'] + c] = fmt(up)
        lines[start['np_11'] + c] = fmt(n_p[c])
    open(tec, 'w').write('\n'.join(lines))

    # what ICE will see: its radius from (rho_p, n_p, rho_mat) against dp/2
    rp_ice = [(0.75 * r / (n * math.pi * rho_mat)) ** (1.0 / 3.0) for r, n in zip(rho_p, n_p)]
    err = max(abs(x - rp) / rp for x in rp_ice)
    print('[make_case] %s: %d cells' % (os.path.relpath(tec, case), nc))
    print('[make_case]   rho_p %.15e .. %.15e kg/m3 (ATLAS, kept)' % (min(rho_p), max(rho_p)))
    print('[make_case]   up_11 set to %s m/s' % fmt(up).strip())
    print('[make_case]   np_11 set to rho_p / (rho_mat 4/3 pi rp^3): %s .. %s 1/m3  (rho_mat %g, rp %g)'
          % (fmt(min(n_p)).strip(), fmt(max(n_p)).strip(), rho_mat, rp))
    print('[make_case]   ICE radius from (rho_p, n_p, rho_mat): worst |R_p - dp/2| / (dp/2) = %.2e' % err)


if __name__ == '__main__':
    main()
