#!/usr/bin/env python3
"""Rewrite the columns of INPUT/part-ic.tec that ATLAS ICB cannot set for the ICE aluminium boxes.

Family 1, the burning aluminium, takes the velocity of [ice-box] up and an exact number density, so that ICE's
particle radius is exactly [ICB-Block1] dp / 2: ICB writes the gas velocity and rounds pi to 3.14. Family 2, the
smoke (krho 0), takes ICE's empty state, 1e-25 in each of its six columns: ICB writes zeros for a family of krho 0,
and ICE's drag on a family of zero number density hands a non-finite source. Run it after every regeneration
(test/mi2-checks/regen_case.sh <case-dir> <scratch-dir> GPB ICB BCB --install):

    make_case.py [case-dir]        default: the directory of this script

Used by test/ice-al-box, test/ice-al-box-direct and test/ice-al-freeze-box. It reads <case-dir>/input.ini
([GPB-Phase2] rho, its first value, the aluminium's; [ICB-Block1] dp; [ice-box] up) and rewrites up_11, np_11 and
the six columns of family 2 of <case-dir>/INPUT/part-ic.tec in the number format of the file (ATLAS's, through
test/ice-box/make_case.py's fmt and digits_of); every other line is kept byte for byte. Idempotent."""
import importlib.util
import math
import os
import re
import sys

EMPTY = 1.0e-25
FAMILY2 = ('rp_21', 'up_21', 'vp_21', 'wp_21', 'Tp_21', 'np_21')


def ice_box_helpers():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ice-box', 'make_case.py')
    spec = importlib.util.spec_from_file_location('ice_box_make_case', path)
    mc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mc)
    return mc


def main():
    mc = ice_box_helpers()
    case = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
    tec = os.path.join(case, 'INPUT', 'part-ic.tec')
    cfg = mc.read_ini(os.path.join(case, 'input.ini'))
    rho_mat = float(cfg['GPB-Phase2']['rho'].split()[0])
    dp = float(cfg['ICB-Block1']['dp'])
    up = float(cfg['ice-box']['up'])

    lines = open(tec).read().split('\n')
    names = [t.strip('"') for t in re.findall(r'"[^"]+"|\S+', lines[0].split('=', 1)[1])]
    I, J, K = (int(x) for x in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', lines[1]).groups())
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    start, first = {}, 2
    for n in names:
        start[n] = first
        first += nn if n in ('x', 'y', 'z') else nc
    for key in ('rp_11', 'up_11', 'np_11') + FAMILY2:
        if key not in start:
            sys.exit('[make_case] %s: no column %s (columns: %s)' % (tec, key, ' '.join(names)))
    if len(lines) < first:
        sys.exit('[make_case] %s: %d lines, %d expected' % (tec, len(lines), first))

    nd = mc.digits_of(lines[2])
    vol = 4.0 / 3.0 * math.pi * (0.5 * dp) ** 3
    n_p = []
    for c in range(nc):
        rho_p = float(lines[start['rp_11'] + c])
        n_p.append(rho_p / (rho_mat * vol))
        lines[start['up_11'] + c] = mc.fmt(up, nd)
        lines[start['np_11'] + c] = mc.fmt(n_p[-1], nd)
        for v in FAMILY2:
            lines[start[v] + c] = mc.fmt(EMPTY, nd)
    open(tec, 'w').write('\n'.join(lines))

    print('[make_case] %s: %d cells' % (os.path.relpath(tec, case), nc))
    print('[make_case]   up_11 set to %s m/s' % mc.fmt(up, nd).strip())
    print('[make_case]   np_11 set to rho_p / (rho_mat 4/3 pi (dp/2)^3): %s .. %s 1/m3  (rho_mat %g, dp %g)'
          % (mc.fmt(min(n_p), nd).strip(), mc.fmt(max(n_p), nd).strip(), rho_mat, dp))
    print('[make_case]   family 2 set to %s in its six columns' % mc.fmt(EMPTY, nd).strip())


if __name__ == '__main__':
    main()
