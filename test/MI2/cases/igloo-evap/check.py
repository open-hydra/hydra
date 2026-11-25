#!/usr/bin/env python3
"""igloo-evap: water drops evaporate in IGLOO into the H2O of a closed box of N2 (the evaporation record).

  check.py <run-dir>

Four streams of 20 um drops, each evaporating whole within the box. IGLOO's sources of a sweep act on the next MOSE
block unchanged; t_src is the time of the blocks with sources.
  E0  the record: A evaporates into H2O, 1 kg per kg
  E1  mass the drops give up: SUM M dt = SUM mdot t_src (each drop burns out: IGLOO leaves 1e-6 of it)
  E2  the gas's H2O (closed-domain balance) = SUM M dt; its N2 untouched
  E3  gas energy: the books (handed vs gained, eight figures) and SUM dE = SUM Q dt (MOSE's floor: 1e-5)
(walls bound the box: the gas's momentum is not judged)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import mi2check as m

run = sys.argv[1]
inp = run + '/INPUT/'
ini = m.read_ini(run + '/input.ini')
txt = m.read_log(run + '/run.log')
n = len(ini['IGLOO-BC']['y'].split())
mdot = m.values(ini, 'IGLOO-BC', 'mdot', n)

rec = m.records(inp + 'part-products.txt', 'A')
m.claim(rec.get('evaporation', {}).get('gas') == {'H2O': 1.0}, 'E0 the record: A -> H2O, 1 kg per kg (%s)' % rec)

steps = [s for s in m.report(txt) if 'M' in s and s['t1'] > s['t0']]
m.claim(len(steps) >= 3, 'E0 %d blocks with sources' % len(steps))
t_src = sum(s['t1'] - s['t0'] for s in steps)
given = sum(s['M'] * (s['t1'] - s['t0']) for s in steps)
m.judge('E1 drops', given, sum(mdot) * t_src, 2e-6, 'SUM M dt vs SUM mdot t_src')

bal = m.balance(txt, m.species(inp))
m.judge('E2 H2O', bal.get('H2O', 0.0), given, 1e-6, 'the gas gained vs SUM M dt (seven figures)')
m.claim('N2' not in bal, 'E2 N2: untouched (no row in the closed-domain balance)')

handed, gained = m.energy_books(txt)
m.judge('E3 energy books', gained, handed, 1e-7, 'the gas gained vs handed (closed-domain bookkeeping)')
dE = sum(s['E1'] - s['E0'] for s in steps)
Qdt = sum(s['Q'] * (s['t1'] - s['t0']) for s in steps)
m.judge('E3 energy', dE, Qdt, 1e-5, 'SUM dE vs SUM Q dt')

m.done('igloo-evap')
