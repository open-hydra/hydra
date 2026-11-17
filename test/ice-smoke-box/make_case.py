#!/usr/bin/env python3
"""Write INPUT/ and MESH/ of ice-smoke-box from those of test/ice-burnout-box.

ATLAS GPB cannot write this case's products record: A and its smoke S are both constant-property
materials, so the reaction of A into N2 and S gives neither a molar mass. The script edits
ice-burnout-box's files instead:

    make_case.py [case-dir]        default: the directory of this script

  part-phase.txt       S, one population (ICE family 2)
  part-models.txt      d-smoke 3 um at the end of A's row
  part-products.txt    A's record: 0.5 kg of N2 and 0.5 kg of S per kg burnt, S born in the phase
                       part- with h298 -1.0e6 J/kg; S's own record (no process)
  part-properties.dat  S's zone: A's temperatures with S's own constants (rho 2500 kg/m3,
                       c 1200 J/(kg K)) and enthalpy c T + c0, c0 such that h(298.15 K) is the
                       record's h298 (-1.0e6 J/kg)
  part-ic.tec          family 2 empty: ICE's empty state, 1e-25 in each of its six columns

Every other file is ice-burnout-box's, byte for byte (MESH/ is test/ice-box's). The sibling
directory ice-burnout-box is read; input.ini, a-n2.yaml and check.py of the case are left alone.
Idempotent."""
import os
import re
import shutil
import sys

NU_S, H_S, D_S, EMPTY = 0.5, -1.0e6, 3.0e-6, 1.0e-25
RHO_S, CP_S, T_REF = 2500.0, 1200.0, 298.15


def edit(path, old, new):
    s = open(path).read()
    if s.count(old) != 1:
        sys.exit('[make_case] %s: the text to replace is found %d times: %s' % (path, s.count(old), old[:60]))
    open(path, 'w').write(s.replace(old, new))


def main():
    case = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__))
    src = os.path.join(os.path.dirname(case), 'ice-burnout-box')
    for name in ('INPUT', 'MESH'):
        if os.path.exists(os.path.join(case, name)):
            shutil.rmtree(os.path.join(case, name))
        shutil.copytree(os.path.join(src, name), os.path.join(case, name))
    I = os.path.join(case, 'INPUT') + '/'
    with open(I + 'part-phase.txt', 'a') as f:
        f.write('S 1\n')
    edit(I + 'part-models.txt', ' d-burnout=3e-06\n', ' d-burnout=3e-06 d-smoke=%g\n' % D_S)
    edit(I + 'part-products.txt', """material A combustion 1 formation 1.0000000000e+06
  N2 1.0000000000e+00 1.0000000000e+00
""", """material A combustion 1 formation 1.0000000000e+06 S %.10e %.10e part-
  N2 %.10e %.10e
material S none 0 none 0.0000000000e+00
""" % (NU_S, H_S, 1.0 - NU_S, 1.0 - NU_S))
    lines = open(I + 'part-properties.dat').read().splitlines(True)
    z = [i for i, l in enumerate(lines) if l.startswith('ZONE')]
    if len(z) != 1:
        sys.exit('[make_case] part-properties.dat: %d zones, one expected' % len(z))
    zone = lines[z[0]:]
    zone[0] = zone[0].replace('T="A"', 'T="S"')
    c0 = H_S - CP_S * T_REF
    zone[2:] = ['%.1f %.6f %.6f %.6f\n' % (T, CP_S, RHO_S, CP_S * T + c0)
                for T in (float(l.split()[0]) for l in zone[2:])]
    with open(I + 'part-properties.dat', 'a') as f:
        if not lines[-1].endswith('\n'):
            f.write('\n')
        f.writelines(zone)
    ic = open(I + 'part-ic.tec').read().splitlines(True)
    hdr = ic[0].rstrip('\n') + '   rp_21   up_21   vp_21   wp_21   Tp_21   np_21\n'
    zl = ic[1].replace('[4-9]=CELLCENTERED', '[4-15]=CELLCENTERED')
    if zl == ic[1]:
        sys.exit('[make_case] part-ic.tec: the zone line has no [4-9]')
    ni, nj, nk = (int(x) for x in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', ic[1]).groups())
    ncell = (ni - 1) * (nj - 1) * (nk - 1)
    open(I + 'part-ic.tec', 'w').writelines([hdr, zl] + ic[2:] + [' %.17E\n' % EMPTY] * (6 * ncell))
    print('[make_case] %s: INPUT/ and MESH/ written from %s' % (case, src))


if __name__ == '__main__':
    main()
