#!/usr/bin/env python3
"""prep.py mesh <case> | empty <case>: the two steps of a steady case (ice-burnout-steady) that ATLAS does
not do.
  mesh   writes <case>/MESH/mesh.tec: test/MI2/validation/JPL-Eulerian-20micron's mesh with every 5th node line in i and in j kept
         (190 x 40 cells -> 38 x 8), before ATLAS GPB ICB BCB (test/MI2/regen_case.sh <case> <scratch> GPB ICB
         BCB --install).
  empty  after the regeneration, sets family 2 of <case>/INPUT/part-ic.tec to ICE's empty state, 1e-25 in every slot
         (ICB writes zeros for krho 0, on which ICE's drag is NaN): test/steady-ice-smoke's smoke family.
Gate: `test/MI2/run.sh ice-burnout-steady`."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'JPL-Eulerian-20micron', 'MESH', 'mesh.tec')
STEP = 5
EMPTY = 1.0e-25


def mesh(case):
    lines = open(SRC).read().split('\n')
    I, J, K = (int(x) for x in re.search(r'I=(\d+), J=(\d+), K=(\d+)', lines[1]).groups())
    vals = ' '.join(lines[2:]).split()
    n = I * J * K
    if len(vals) != 3 * n or (I - 1) % STEP or (J - 1) % STEP:
        sys.exit('steady_case_prep: %s is not a %d-coarsenable one-block mesh' % (SRC, STEP))
    ii, jj = range(0, I, STEP), range(0, J, STEP)
    out = [vals[c * n + i + j * I + k * I * J] for c in range(3) for k in range(K) for j in jj for i in ii]
    head = lines[1].replace('I=%d, J=%d' % (I, J), 'I=%d, J=%d' % (len(ii), len(jj)))
    os.makedirs(os.path.join(case, 'MESH'), exist_ok=True)
    with open(os.path.join(case, 'MESH', 'mesh.tec'), 'w') as f:
        f.write(lines[0] + '\n' + head + '\n' + '\n'.join(' ' + v for v in out) + '\n')
    print('%s/MESH/mesh.tec: %d x %d cells' % (case, len(ii) - 1, len(jj) - 1))


def empty(case):
    tec = os.path.join(case, 'INPUT', 'part-ic.tec')
    lines = open(tec).read().split('\n')
    names = lines[0].split('=', 1)[1].split()
    I, J, K = (int(x) for x in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', lines[1]).groups())
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(K - 1, 1)
    vals = ' '.join(lines[2:]).split()
    cols, off = {}, 0
    for n in names:
        k = nn if n in ('x', 'y', 'z') else nc
        cols[n] = vals[off:off + k]
        off += k
    if off != len(vals) or 'rp_21' not in cols:
        sys.exit('steady_case_prep: %s has no family 2 in the layout expected' % tec)
    for v in ('rp_21', 'up_21', 'vp_21', 'wp_21', 'Tp_21', 'np_21'):
        cols[v] = ['%.17E' % EMPTY] * nc
    with open(tec, 'w') as f:
        f.write(lines[0] + '\n' + lines[1] + '\n')
        for n in names:
            for v in cols[n]:
                f.write(' %s\n' % v)
    print('%s: family 2 empty (%g, %d cells)' % (tec, EMPTY, nc))


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] not in ('mesh', 'empty'):
        sys.exit(__doc__)
    {'mesh': mesh, 'empty': empty}[sys.argv[1]](sys.argv[2].rstrip('/'))
