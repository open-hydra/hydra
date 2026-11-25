#!/usr/bin/env python3
"""ch-gas-ice: ice-al with one law channel, steam of the gas born as alumina drops in ICE family 2:
    cond = gas:H2O ICE:AL2O3(L) constant rate=R birth-d=3e-6

  check.py <run-dir>      (ice-al's run beside it: <run-dir>/../ice-al)

MOSE steps in pseudo-time here (a channel with an ICE end needs a steady MOSE or ICE), so the gas keeps no
books: its side is judged on the rates. The channel moves R V per second (V the box's volume); everything else is
ice-al's, so the ICE claim is the run's number minus ice-al's:
  C0  the coupler reports the channel
  C1  ICE: family 2 gains R V dt over its one (time-accurate) step
  C2  the gas: every block hands it -R V of H2O (channels-to-gas), and the channel reports moving R V
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import mi2check as m

run = sys.argv[1]
base = os.path.join(os.path.dirname(os.path.normpath(run)), 'ice-al')
inp = run + '/INPUT/'
ini = m.read_ini(run + '/input.ini')
txt = m.read_log(run + '/run.log')
R = float(ini['MI2-Channels']['cond'].split('rate=')[1].split()[0])
m.claim('[MI2] channel cond: gas H2O -> ICE AL2O3(L), law constant (born as new drops)' in txt,
        'C0 the coupler reports the channel')

H = m.zones(inp + 'part-properties.dat', 3)['AL2O3(L)']
b, vol = m.field(run + '/OUTPUT/part-field.tec')
b0, _ = m.field(base + '/OUTPUT/part-field.tec')
V = sum(vol)
gain = m.family(b, vol, 2, H, 0.0)['M'] - m.family(b0, vol, 2, H, 0.0)['M']
m.judge('C1 ICE', gain, R * V * m.ice_dt(txt), 1e-9, 'family 2 gained (with - without) vs R V dt')

steps = m.report(txt)
worst = max(abs(s['chM'] + R * V) / (R * V) for s in steps if 'chM' in s)
m.claim(worst <= 1e-12, 'C2 gas: %d blocks hand it -R V = %.6e kg/s, worst rel.err %.2e'
        % (len(steps), -R * V, worst))
worst = max(abs(s['moved']['cond'] - R * V) / (R * V) for s in steps if 'cond' in s['moved'])
m.claim(worst <= 1e-12, 'C2 the channel reports moving R V, worst rel.err %.2e' % worst)

m.done('ch-gas-ice')
