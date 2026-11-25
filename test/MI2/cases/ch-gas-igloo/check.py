#!/usr/bin/env python3
"""ch-gas-igloo: igloo-al with one law channel, steam of the gas absorbed by IGLOO's parcels as alumina:
    cond = gas:H2O IGLOO:AL2O3(L) constant rate=R

  check.py <run-dir>      (igloo-al's run beside it: <run-dir>/../igloo-al)

IGLOO takes the offer in its next sweep and reports what its parcels absorbed; the gas gives up what IGLOO reports
(lagged books), so the gas's side follows IGLOO's, block by block. With V the box's volume:
  A0  the coupler reports the channel
  A1  the lag: the blocks before IGLOO's first report hand the gas nothing; the first report is one block later
  A2  steady absorption: the parcels take the offer R of the cells they cross, the same every block (1e-5), and
      no more than R V
  A3  the gas's H2O: gained (with - without) = SUM channels-to-gas dt (seven figures of the H2O gained)
  A4  the books close: mass of every species (closed-domain balance); energy handed vs gained to MOSE's floor (4e-10
      of the gas's energy, igloo-al's), and the channel adds nothing to that gap
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import mi2check as m

run = sys.argv[1]
base = os.path.join(os.path.dirname(os.path.normpath(run)), 'igloo-al')
inp = run + '/INPUT/'
ini = m.read_ini(run + '/input.ini')
txt, txt0 = m.read_log(run + '/run.log'), m.read_log(base + '/run.log')
R = float(ini['MI2-Channels']['cond'].split('rate=')[1].split()[0])
m.claim('[MI2] channel cond: gas H2O -> IGLOO AL2O3(L), law constant (absorbed)' in txt,
        'A0 the coupler reports the channel')

V = sum(m.field(run + '/OUTPUT/partE-field.tec')[1])
steps = [s for s in m.report(txt) if 'chM' in s and s['t1'] > s['t0']]
k = next(i for i, s in enumerate(steps) if s['chM'] != 0.0)
first = [s for s in m.report(txt) if 'F' in s and s['t1'] > s['t0']][0]
m.claim(k >= 1 and steps[k]['t0'] > first['t0'], 'A1 lag: %d block(s) hand nothing, the first absorption at t0 %.6e s '
        '(IGLOO\'s first sources at %.6e s)' % (k, steps[k]['t0'], first['t0']))
full = steps[k + 1:]
a = -full[-1]['chM']
worst = max(abs(-s['chM'] - a) / a for s in full)
m.claim(full and worst <= 1e-5 and 0.0 < a <= R * V, 'A2 %d blocks absorb %.9e kg/s (%.3e of R V), worst spread %.2e'
        % (len(full), a, a / (R * V), worst))

names = m.species(inp)
lost = sum(s['chM'] * (s['t1'] - s['t0']) for s in steps)
dH2O = m.balance(txt, names)['H2O'] - m.balance(txt0, names)['H2O']
m.judge('A3 gas H2O', dH2O, lost, 1e-6, 'gained (with - without) vs SUM channels-to-gas dt',
        scale=m.balance(txt, names)['H2O'])

for l in txt.split('closed-domain mass balance')[1].split('\n')[2:]:
    f = l.split()
    if len(f) != 4 or not f[0].isdigit():
        break
    m.judge('A4 mass ' + names[int(f[0]) - 1], float(f[2]), float(f[1]), 1e-6, 'gas gained vs handed')
handed, gained = m.energy_books(txt)
handed0, gained0 = m.energy_books(txt0)
E_gas = abs(m.report(txt)[0]['E0'])
m.judge('A4 energy', gained, handed, 1e-8, 'gas gained vs handed, of the gas\'s %.6e J' % E_gas, scale=E_gas)
m.judge('A4 energy gap', gained - handed, gained0 - handed0, 1e-4, 'with vs without, of the channel\'s %.3e J'
        % (handed - handed0), scale=handed - handed0)

m.done('ch-gas-igloo')
