#!/usr/bin/env python3
"""ice-al: liquid aluminium (ICE family 1) burns with steam, 2 AL(L) + 3 H2O => AL2O3(L) + 3 H2; the alumina it
generates is born in family 2 of the same phase.

  check.py <run-dir>

Cycle 1 is one forward-Euler ICE step dt at the initial state; its source is reported as rates (ice-to-gas F Q, the
gas's share M through the record's species) and is held through cycle 2 (t2 = SUM t1 - t0 of its blocks).
  B1  mass: generated = nu_s x the aluminium lost; lost = SUM g x lost + generated; M = SUM g x lost / dt
  B2  momentum: aluminium's lost = F dt + generated's; energy (formation scale): aluminium's lost = Q dt + generated's
  B3  birth state: forward Euler takes the source at the initial state: the drops' u_p and T_p
  B4  number: n = rho / (pi/6 rho_s d_gen^3) in every cell
  G1  the gas's species over cycle 2 (closed-domain balance) = g(s) (lost / dt) t2
"""
import math
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import mi2check as m

run = sys.argv[1]
inp = run + '/INPUT/'
ini = m.read_ini(run + '/input.ini')
txt = m.read_log(run + '/run.log')
comb = m.records(inp + 'part-products.txt', 'AL(L)')['combustion']
nu_s, h298_s, h298 = float(comb['tail'][1]), float(comb['tail'][2]), comb['h298']
m.claim(comb['tail'][0] == 'AL2O3(L)' and comb['tail'][3] == 'part-', 'B0 AL(L) burns into %.10e kg of %s per kg, '
        'born in %s' % (nu_s, comb['tail'][0], comb['tail'][3]))
m.claim('ICE material AL(L) generates AL2O3(L) (family 2)' in txt, 'B0 the coupler reports the product and its family')

dt = m.ice_dt(txt)
H = m.zones(inp + 'part-properties.dat', 3)
a, vol = m.ic_field(inp + 'part-ic.tec')
b, _ = m.field(run + '/OUTPUT/part-field.tec')
al0, al1 = m.family(a, vol, 1, H['AL(L)'], h298), m.family(b, vol, 1, H['AL(L)'], h298)
gen = m.family(b, vol, 2, H['AL2O3(L)'], h298_s)
lost = al0['M'] - al1['M']
first = [s for s in m.report(txt) if 'iceF' in s][0]
F, Q, M = first['iceF'], first['iceQ'], first['iceM']
g_sum = sum(comb['gas'].values())
print('      the step: dt %.6e s, aluminium lost %.10e kg, generated %.10e kg' % (dt, lost, gen['M']))

m.judge('B1 generated', gen['M'], nu_s * lost, 1e-12, 'family 2 vs nu_s x the aluminium lost')
# the record's gas row prints eleven figures (5e-11 measured)
m.judge('B1 mass', lost, g_sum * lost + gen['M'], 2e-10, 'lost vs SUM g x lost + generated')
m.judge('B1 gas share', M, g_sum * lost / dt, 2e-10, 'ice-to-gas-mass vs SUM g x lost / dt')
m.judge('B2 momentum', al0['P'][0] - al1['P'][0], F[0] * dt + gen['P'][0], 1e-9, 'lost vs F dt + generated\'s')
m.judge('B2 energy', al0['E'] - al1['E'], Q * dt + gen['E'], 1e-7, 'lost vs Q dt + generated\'s')

up = float(ini['ice-exchange']['up'])
Tp = float(ini['ICB-Block1']['T']) * float(ini['ICB-Block1']['kT'].split()[0])
Tm = sum(b['rho_p2'][c] * vol[c] * b['T_p2'][c] for c in gen['cells']) / gen['M']
m.judge('B3 birth velocity', gen['P'][0] / gen['M'], up, 1e-9, 'the generated mean u vs the drops\' initial u_p')
m.judge('B3 birth temperature', Tm, Tp, 1e-9, 'the generated mean T vs the drops\' initial T_p')

rho = m.zones(inp + 'part-properties.dat', 2)['AL2O3(L)']
models = [l for l in open(inp + 'part-models.txt') if l.startswith('AL(L) combustion')][0]
d_gen = float(models.split('d-gen=')[1].split()[0])
worst = max(abs(b['rho_p2'][c] / b['n_p2'][c] / (math.pi / 6.0 * m.at(rho, b['T_p2'][c]) * d_gen ** 3) - 1.0)
            for c in gen['cells'])
m.claim(worst <= 1e-12, 'B4 number: rho / n = pi/6 rho_s d_gen^3 in %d cells, worst rel.err %.2e'
        % (len(gen['cells']), worst))

t2 = sum(s['t1'] - s['t0'] for s in m.report(txt) if s['cycle'] == 2)
bal = m.balance(txt, m.species(inp))
for sp, nu in comb['gas'].items():
    m.judge('G1 ' + sp, bal.get(sp, 0.0), nu * lost / dt * t2, 1e-6, 'the gas gained over cycle 2 vs g (lost/dt) t2')

m.done('ice-al')
