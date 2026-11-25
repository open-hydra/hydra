#!/usr/bin/env python3
"""ch-steady: igloo-al in pseudo-time (MOSE, ICE and IGLOO steady) with three law channels, one per route that
needs a steady run:
    x1 = IGLOO:AL2O3(L) ICE:AL2O3(L)  constant   rate=R1 birth-d=3e-6
    x2 = ICE:AL2O3(L)   IGLOO:AL2O3(L) relaxation k=K
    x3 = gas:H2O        ICE:AL2O3(L)  constant   rate=R3 birth-d=3e-6

  check.py <run-dir>

A steady run keeps no books: the claims are on the rates each block reports (V the box's volume).
  S0  the coupler reports the three channels; the run ends (GAME OVER)
  S1  x3: R3 V every block, and it is all the gas gives up (channels-to-gas)
  S2  x1: IGLOO gives up mass from the block after its first sweep, never more than R1 V (only its own mass)
  S3  x2: ICE's alumina reaches IGLOO's parcels, from IGLOO's sweep after x1's first (IGLOO reports what it took)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import mi2check as m

run = sys.argv[1]
ini = m.read_ini(run + '/input.ini')
txt = m.read_log(run + '/run.log')
ch = ini['MI2-Channels']
R1 = float(ch['x1'].split('rate=')[1].split()[0])
R3 = float(ch['x3'].split('rate=')[1].split()[0])
for line in ['x1: IGLOO AL2O3(L) -> ICE AL2O3(L), law constant (born as new drops)',
             'x2: ICE AL2O3(L) -> IGLOO AL2O3(L), law relaxation (absorbed)',
             'x3: gas H2O -> ICE AL2O3(L), law constant (born as new drops)']:
    m.claim('[MI2] channel ' + line in txt, 'S0 the coupler reports channel ' + line)
m.claim('GAME OVER' in txt, 'S0 the run ends')

V = sum(m.field(run + '/OUTPUT/partE-field.tec')[1])
steps = m.report(txt)
moved = [s['moved'] for s in steps]
worst = max(abs(x['x3'] - R3 * V) / (R3 * V) for x in moved)
m.claim(worst <= 1e-12, 'S1 x3 moves R3 V = %.6e kg/s in %d blocks, worst rel.err %.2e' % (R3 * V, len(steps), worst))
worst = max(abs(s['chM'] + R3 * V) / (R3 * V) for s in steps)
m.claim(worst <= 1e-12, 'S1 the gas gives up R3 V and nothing else, worst rel.err %.2e' % worst)

x1 = [x['x1'] for x in moved]
m.claim(x1[0] == 0.0 and all(0.0 < a <= R1 * V for a in x1[2:]),
        'S2 x1: nothing before IGLOO\'s first sweep, then %.3e to %.3e kg/s (R1 V = %.3e)'
        % (min(x1[2:]), max(x1[2:]), R1 * V))
x2 = [x['x2'] for x in moved]
k = next(i for i, a in enumerate(x1) if a > 0.0)
m.claim(x2[k] == 0.0 and all(a > 0.0 for a in x2[k + 1:]), 'S3 x2: from the sweep after x1\'s first, %.3e to %.3e '
        'kg/s' % (min(x2[k + 1:]), max(x2[k + 1:])))

m.done('ch-steady')
