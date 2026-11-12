#!/usr/bin/env python3
"""JPL-3solver: gas vs Chang (axis + wall), and the two condensed phases vs each other on the same gas.
Usage: compare_3solver.py <run-dir-with-OUTPUT> [igloo-euler-file]"""
import sys, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'mi2-checks'))
sys.path.insert(0, os.path.join(HERE, '..', 'JPL-Lagrangian-20micron'))
from tec_cc import read_cc
import jpl_compare
import numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))


def first_401_field1(path):
    """Field 1 of the first dispersed 401 record: the phase's share of the inlet flow, krho/(1 - krho_tot)."""
    lines = open(path).read().splitlines()
    for i, line in enumerate(lines):
        t = line.split()
        if len(t) == 6 and all(x.isdigit() for x in t) and t[5] == '401':
            return float(lines[i + 1].split()[0])
    sys.exit('no dispersed 401 record in ' + path)


K_PHASE = first_401_field1(os.path.join(d, 'INPUT', 'partL-bc.txt'))   # each phase's share of the inlet flow
gas = d + '/OUTPUT/gas-field.tec'
ref = os.path.join(HERE, '..', 'JPL-Lagrangian-20micron', 'JPL_chang.dat')
jpl_compare.report(gas, ref, '3-solver gas (half load each phase)')
jpl_compare.report(gas, ref, '3-solver gas WALL', zone='Chang MPL Wall', jcell=39)

g = read_cc(gas)
rhou = (g['rho(1)'][0, :, 0] * g['u'][0, :, 0]).mean()
print(f'inlet gas rho*u mean {rhou:.4f}; half-load prescription {K_PHASE:.5f}*rho*u = {K_PHASE*rhou:.4f} per phase')

ice = read_cc(d + '/OUTPUT/partE-field.tec')
eul_files = sorted(glob.glob(d + '/OUTPUT/partL-euler1-sweep*.tec'), key=lambda f: int(f.rsplit('sweep', 1)[1].split('.')[0]))
igl = read_cc(sys.argv[2] if len(sys.argv) > 2 else eul_files[-1])
print('IGLOO euler file:', os.path.basename(eul_files[-1]), 'vars:', [k for k in igl if k not in ('x', 'y', 'z')])
print('ICE   field vars:', [k for k in ice if k not in ('x', 'y', 'z')])
ri, ui = ice['rho_p1'][0], ice['u_p1'][0]
kr = [k for k in igl if 'rho' in k.lower()][0]; ku = [k for k in igl if k.lower().startswith('u')][0]
rl, ul = igl[kr][0], igl[ku][0]
print(f'inlet-column condensed flux  ICE {(ri[:, 0]*ui[:, 0]).mean():.3f}   IGLOO {(rl[:, 0]*ul[:, 0]).mean():.3f}   (expected ~ {K_PHASE*rhou:.3f} each)')
# exit-plane radial profile of condensed density, both phases (i = last column)
print('exit column rho_p, j = 0(axis) .. 40(wall) every 5 cells:')
print('   ICE  ', ' '.join(f'{v:8.4f}' for v in ri[::5, -1]))
print('   IGLOO', ' '.join(f'{v:8.4f}' for v in rl[::5, -1]))
tot_i, tot_l = ri.sum(), rl.sum()
print(f'sum of exit... domain-sum rho_p: ICE {tot_i:.3f}  IGLOO {tot_l:.3f}  ratio IGLOO/ICE {tot_l/tot_i:.3f}')
