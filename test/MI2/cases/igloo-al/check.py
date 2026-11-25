#!/usr/bin/env python3
"""igloo-al, igloo-al-release: liquid aluminium burns with steam in IGLOO, 2 AL(L) + 3 H2O => AL2O3(L) + 3 H2; the
alumina it generates is born in ICE (phase partE).

  check.py <run-dir>

The carry record keeps a share xi of the metal's oxide on the drop. igloo-al leaves it as a residue parcel;
igloo-al-release (residue = off) hands it to ICE at burnout through the oxide's release record. With t_src the time
of the MOSE blocks with IGLOO's sources and SUM mdot the injected flow, f the metal burnt into the combustion record
per kg injected (1 - xi), and the residue r = xi nu_ox per kg injected (none with the release):
  K0  the records, the coupler's report of the product's ICE family, the residue's printed size
  K1  generated mass (ICE's field) = nu_s f SUM mdot t_src (+ xi nu_ox nu_rel SUM mdot t_src released)
  K2  the gas's species (closed-domain balance) = g(s) SUM mdot t_src
  K3  momentum: gas gained + ICE's = SUM mdot u_inj t_src (the residue's, bounded by its |u_end|: the tolerance)
  K4  energy (formation scale): gas gained + ICE's + residue's = SUM mdot (H(T0) + u_inj^2/2) t_src
  K5  number: SUM n vol = generated mass / (pi/6 rho d_gen^3)
  K6  hand-off: ICE's field holds what the igloo-gen-delivered lines report
"""
import math
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import mi2check as m

run = sys.argv[1]
case = os.path.basename(os.path.normpath(run))
inp = run + '/INPUT/'
ini = m.read_ini(run + '/input.ini')
txt = m.read_log(run + '/run.log')
n = len(ini['IGLOO-BC']['y'].split())
mdot, up = m.values(ini, 'IGLOO-BC', 'mdot', n), m.values(ini, 'IGLOO-BC', 'up', n)
T0, d0 = m.values(ini, 'IGLOO-BC', 'temp0', n), m.values(ini, 'IGLOO-BC', 'diam', n)
g2 = ini['GPB-Phase2']
xi = float(g2['share'])
released = g2.get('residue', 'on') == 'off'

# --- K0: the records ---
al = m.records(inp + 'partL-products.txt', 'AL(L)')
comb, carry = al['combustion'], al['carry']
prod, nu_s, h298_s, phase = comb['tail'][0], float(comb['tail'][1]), float(comb['tail'][2]), comb['tail'][3]
nu_ox = float(carry['tail'][1])
m.claim(prod == 'AL2O3(L)' and phase == 'partE-', 'K0 AL(L) burns into %.10e kg of %s per kg, born in %s'
        % (nu_s, prod, phase))
m.claim('IGLOO material AL(L) generates AL2O3(L)' in txt, 'K0 the coupler reports the product and its ICE family')
steps = [s for s in m.report(txt) if 'F' in s and s['t1'] > s['t0']]
t_src = sum(s['t1'] - s['t0'] for s in steps)
m.claim(len(steps) >= 3, 'K0 %d blocks with sources, t_src %.10e s' % (len(steps), t_src))
M_in = sum(mdot) * t_src

HC = m.zones(inp + 'partL-properties.dat', 3)['AL(L)']
HO = m.zones(inp + 'partE-properties.dat', 3)['AL2O3(L)']
rho_c, rho_o = [float(x) for x in g2['rho'].split()[:2]]

# --- the residue of each block: the oxide parcels' end records of the sweep that made the block's sources ---
res = {'M': 0.0, 'P': 0.0, 'E': 0.0}
if released:
    rel = m.records(inp + 'partL-products.txt', 'AL2O3(L)')['release']
    nu_rel = float(rel['tail'][1])
    m.claim(rel['tail'][3] == 'partE-', 'K0 residue = off: the oxide is released as %.10e kg per kg' % nu_rel)
    gen_per_kg = nu_s * (1.0 - xi) + xi * nu_ox * nu_rel
else:
    worst = 0.0
    for k, s in enumerate(steps):
        dt = s['t1'] - s['t0']
        end = m.outloc(run + '/OUTPUT/partL-outloc-AL2O3(L)%s.dat' % ('' if k == 0 else '-sweep%d' % k))
        for i in range(n):
            Te, ue, de = end[i + 1]
            r = xi * nu_ox
            worst = max(worst, abs(rho_o * de ** 3 / (rho_c * d0[i] ** 3) / r - 1.0))
            res['M'] += mdot[i] * r * dt
            res['P'] += mdot[i] * r * dt * ue
            res['E'] += mdot[i] * r * dt * (m.at(HO, Te) + 0.5 * ue * ue)
    m.claim(worst <= 2e-5, 'K0 residue xi nu_ox = %.9f per kg; the printed d_end gives it within %.2e (six figures)'
            % (xi * nu_ox, worst))
    gen_per_kg = nu_s * (1.0 - xi)

# --- K1-K6 ---
cols, vol = m.field(run + '/OUTPUT/partE-field.tec')
g = m.family(cols, vol, 1, HO, h298_s)
# the coupler splits the metal through the record's nu_ox, IGLOO's carry law holds its own: they agree to 1e-7 (two
# tables of atomic weights), which moves the generated mass by xi/(1 - xi) of it
m.judge('K1 generated mass', g['M'], gen_per_kg * M_in, 3e-8, 'ICE field vs SUM mdot t_src x %.9f' % gen_per_kg)

bal = m.balance(txt, m.species(inp))
for sp, nu in comb['gas'].items():
    m.judge('K2 ' + sp, bal.get(sp, 0.0), nu * M_in, 1e-6, 'the gas gained vs g SUM mdot t_src (seven figures)')

dpx = sum(s['p1'][0] - s['p0'][0] for s in steps)
F_in = sum(a * u for a, u in zip(mdot, up)) * t_src
# the outloc files print the residue's speed, not its velocity: its momentum is bounded by it, the tolerance
m.judge('K3 momentum x', dpx + g['P'][0], F_in, 1e-6 + res['P'] / F_in,
        'gas %.6e + ICE %.6e vs SUM mdot u_inj t_src (the residue: below %.3e)' % (dpx, g['P'][0], res['P']))

dE = sum(s['E1'] - s['E0'] for s in steps)
E_in = sum(a * (m.at(HC, t) + 0.5 * u * u) for a, u, t in zip(mdot, up, T0)) * t_src
E_gas = abs(m.report(txt)[0]['E0'])
# judged against the gas's total energy (MOSE sums it over the cells: its round-off is the floor, 4e-10 measured)
m.judge('K4 energy', dE + g['E'] + res['E'], E_in, 1e-8,
        'gas %.9e + ICE %.9e + residue %.9e J vs SUM mdot (H(T0) + u^2/2) t_src, of the gas\'s %.6e J'
        % (dE, g['E'], res['E'], E_gas), scale=E_gas)
m.claim(abs(g['E']) >= 1e-5 * E_gas, 'K4 non-vacuity: ICE holds %.3e J, %.1e of the gas\'s energy'
        % (abs(g['E']), abs(g['E']) / E_gas))

Tm = sum(cols['rho_p1'][c] * vol[c] * cols['T_p1'][c] for c in g['cells']) / g['M']
rho_s = m.at(m.zones(inp + 'partE-properties.dat', 2)['AL2O3(L)'], Tm)
m.judge('K5 number', g['N'], g['M'] / (math.pi / 6.0 * rho_s * float(g2['d-gen']) ** 3), 1e-12,
        'SUM n vol vs mass / (pi/6 rho d_gen^3) at %.4f K' % Tm)

delivered = sum(float(l.split()[2]) for l in txt.split('\n') if 'igloo-gen-delivered' in l)
m.judge('K6 hand-off', g['M'], delivered, 1e-12, 'ICE field mass vs the delivered lines')

m.done(case)
