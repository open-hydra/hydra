#!/usr/bin/env python3
"""The checks of gate p30: liquid aluminium burns with steam, 2 AL(L) + 3 H2O => AL2O3(L) + 3 H2, and the liquid
alumina is the smoke, born in an ICE family (the tail of the combustion record that ATLAS GPB writes).

  p30_al_smoke.py igloo|ice <run-dir> <log> [evap-box-mom-dir]

run-dir holds input.ini, INPUT/ and OUTPUT/ of the run; log its log with the [MI2-report] lines; evap-box-mom-dir is
test/evap-box-mom, whose check.py gives the readers (default: ../evap-box-mom beside this script's folder).
g(s) is the gas row of the combustion record (kg of species s per kg of aluminium burnt; H2O negative), nu_s its smoke
(kg of alumina per kg burnt), h298 and h298_s their formation enthalpies at T_ref = 298.15 K; H_c, H_o the enthalpy
tables of AL(L) and AL2O3(L) (absolute, the gas's formation scale). The smoke's energy on that scale: per kg
h298_s + H_o(T) - H_o(T_ref) + |u|^2/2, from ICE's field.

mode igloo (test/igloo-al-box: the oxide cap, xi-cap, its residue left; test/igloo-al-box-direct: direct down to
d-stop). Over the MOSE blocks with IGLOO's sources, t_src = SUM (t1 - t0), each stream of mass flow mdot injected at
u_inj, T0 ends as a residue parcel of mass fraction r (per kg injected) at T_end, u_end (OUTPUT/partL-outloc-*.dat of
the sweep whose sources the block carries): the cap's oxide, r = xi nu_ox (nu_ox = nu_s, the cap record's kg of oxide
per kg), with the metal burnt whole; or the aluminium left at d-stop, r = (d_end/d0)^3, the burnt fraction 1 - r.
With f the aluminium burnt into the combustion record per kg injected (1 - xi with the cap, 1 - r at d-stop) and f_g
the aluminium whose gas terms the gas receives (1 with the cap, whose record has the same gas row; 1 - r at d-stop):
  K1  smoke mass (ICE's field)                 = nu_s f SUM mdot t_src
  K2  the gas's species (closed-domain balance) = g(s) f_g SUM mdot t_src
  K3  gas momentum gained + smoke's + residue's = SUM mdot u_inj t_src  (the residue's |u_end|: the tolerance)
  K4  gas energy gained + smoke's + residue's (r SUM mdot t_src (H_res(T_end) + u_end^2/2))
      = SUM mdot (H_c(T0) + u_inj^2/2) t_src
  K5  (direct only) the smoke's temperature is T-burn: the drops burn at it, and the birth state is burn-weighted
  K6  number: SUM n vol = smoke mass / (pi/6 rho_s d_smoke^3)
  K7  hand-off: the field holds what the igloo-smoke-delivered lines report
mode ice (test/ice-al-box, -direct): cycle 1 is one forward-Euler ICE step dt at the initial state; family 1 is the
aluminium, family 2 the smoke (empty at the start: 1e-25). The step's source is reported as rates (ice-to-gas F Q,
ice-to-gas-mass M through the record's species, the smoke taken out) and frozen through cycle 2 (t2 = SUM t1 - t0):
  B1  mass of the step: aluminium lost (initial field - final) = M dt + smoke;  smoke = nu_s x aluminium lost
  B2  momentum: aluminium's lost = F dt + smoke's;  energy: aluminium's lost (formation scale) = Q dt + smoke's
  B3  birth state: forward Euler takes the source at the initial state: the smoke at the drops' u_p and T_p
  B4  number: n_s = rho_s / (pi/6 rho_s d_smoke^3) in every cell with smoke
  G1  the gas's species over cycle 2 (closed-domain balance) = g(s) (aluminium lost / dt) t2
One PASS/FAIL line per check; exit 1 on any FAIL."""
import math
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
mom = os.path.abspath(sys.argv[4]) if len(sys.argv) > 4 else os.path.join(here, '..', 'evap-box-mom')
sys.path.insert(0, mom)
import check as eb          # test/evap-box-mom/check.py: read_ini, values, read_report

T_REF = 298.15
TOL_M = 1e-12        # K1 K6 K7 B1 B4: sums of the field against sums of the same terms
TOL_BOOK = 1e-6      # K2 G1: the closed-domain balance prints seven figures
TOL_E = 1e-7         # B2: the energy of one ICE step against its reported rate (3e-10 measured)
TOL_E_GAS = 1e-8     # K4: of the gas's total energy (its round-off over the blocks: 4e-10 measured)
TOL_T = 1e-9         # K5 B3: a state the smoke copies
EMPTY = 1.0e-25      # ICE's empty state of a family
npass = nfail = 0


def ok(msg):
    global npass
    npass += 1
    print('PASS  ' + msg)


def ko(msg):
    global nfail
    nfail += 1
    print('FAIL  ' + msg)


def judge(name, a, b, tol, what, scale=None):
    den = abs(scale if scale is not None else b)
    e = abs(a - b) / den if den > 0.0 else (0.0 if a == b else float('inf'))
    (ok if e <= tol else ko)('%s: %s %.16e vs %.16e, rel.err %.2e (tol %.0e)' % (name, what, a, b, e, tol))
    return e


def zones(path, col):
    """{zone: {T: column col}} of an ATLAS property table (integer-kelvin rows)."""
    out, cur = {}, None
    for line in open(path):
        if line.startswith('ZONE'):
            cur = out.setdefault(re.search(r'T="([^"]*)"', line).group(1), {})
            continue
        if cur is None:
            continue
        try:
            v = [float(x) for x in line.split()]
        except ValueError:
            continue
        if len(v) > col:
            cur[v[0]] = v[col]
    return out


def at(rows, T):
    lo = float(int(T))
    return rows[lo] if T == lo else rows[lo] + (rows[lo + 1.0] - rows[lo]) * (T - lo)


def field(path):
    """Cell columns and cell volumes of an ICE Tecplot ASCII field (BLOCK: nodes, then the cell-centred variables)."""
    lines = open(path).read().split('\n')
    names = lines[0].split('=')[1].split()
    I, J, K = [int(re.search(r'\b%s=(\d+)' % c, lines[1]).group(1)) for c in 'IJK']
    vals = []
    for l in lines[2:]:
        if re.match(r'\s*[A-Za-z]', l):
            continue
        vals += [float(x) for x in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    X = [vals[q * nn:(q + 1) * nn] for q in range(3)]
    cols, off = {}, 3 * nn
    for n in names[3:]:
        cols[n], off = vals[off:off + nc], off + nc
    vol = []
    for k in range(K - 1):
        for j in range(J - 1):
            for i in range(I - 1):
                a = i + j * I + k * I * J
                vol.append((X[0][a + 1] - X[0][a]) * (X[1][a + I] - X[1][a]) * (X[2][a + I * J] - X[2][a]))
    return cols, vol


def record(path):
    """The combustion record of AL(L): h298, {species: nu}, nu_s, h298_s, the smoke's phase; the cap record's nu_ox."""
    lines = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    out = {}
    for i, f in enumerate(lines):
        if f[:3] == ['material', 'AL(L)', 'combustion']:
            out.update(h298=float(f[5]), g={t[0]: float(t[1]) for t in lines[i + 1:i + 1 + int(f[3])]},
                       smoke=f[6], nu_s=float(f[7]), h298_s=float(f[8]), phase=f[9])
        if f[:3] == ['material', 'AL(L)', 'cap']:
            out['nu_ox'] = float(f[7])
    return out


def gas_species(inp):
    """The species of INPUT/gas-phase.txt in its order (the order of the closed-domain balance's rows)."""
    return [l.split()[0] for l in open(inp + 'gas-phase.txt').read().split('\n')[1:] if l.strip()]


def balance(txt, species):
    """{species: gas gained} of the closed-domain mass balance (rows in the gas phase's order)."""
    m = re.search(r'closed-domain mass balance.*?\n.*?\n((?:\s+\d+\s+.*\n)+)', txt)
    out = {}
    if m:
        for l in m.group(1).split('\n'):
            f = l.split()
            if len(f) == 4 and f[0].isdigit():
                out[species[int(f[0]) - 1]] = float(f[2])
    return out


def smoke_totals(cols, vol, sfx, HO, h298_s):
    r, T = cols['rho_p' + sfx], cols['T_p' + sfx]
    U = [cols['u_p' + sfx], cols['v_p' + sfx], cols['w_p' + sfx]]
    # cells the smoke never reached hold a floor (ICE's empty state, or ICB's for a phase of krho 0: 1e-20 at 0 K)
    full = [c for c in range(len(r)) if r[c] > 1.0e-18 and T[c] >= 1.0]
    floor = sum(r[c] * vol[c] for c in range(len(r)) if c not in set(full))
    print('      %d cells hold no smoke (%.3e kg of floor, left out)' % (len(r) - len(full), floor))
    M = sum(r[c] * vol[c] for c in full)
    P = [sum(r[c] * q[c] * vol[c] for c in full) for q in U]
    E = sum(r[c] * vol[c] * (h298_s + at(HO, T[c]) - at(HO, T_REF) + 0.5 * sum(q[c] ** 2 for q in U)) for c in full)
    N = sum(cols['n_p' + sfx][c] * vol[c] for c in full)
    return M, P, E, N, full


def outloc(path):
    """{ID: (T, |u|, d)} of the end records of an outloc file."""
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) == 12:
            try:
                out[int(f[8])] = (float(f[3]), float(f[4]), float(f[11]))
            except ValueError:
                continue
    return out


def igloo(run, log):
    inp = run + '/INPUT/'
    ini = eb.read_ini(run + '/input.ini')
    n = len(ini['IGLOO-BC']['y'].split())
    mdot = eb.values(ini, 'IGLOO-BC', 'mdot', n)
    up = eb.values(ini, 'IGLOO-BC', 'up', n)
    T0 = eb.values(ini, 'IGLOO-BC', 'temp0', n)
    d0 = eb.values(ini, 'IGLOO-BC', 'diam', n)
    g2 = ini['GPB-Phase2']
    capped = 'xi-cap' in g2
    rec = record(inp + 'partL-products.txt')
    species = gas_species(inp)
    (ok if rec.get('smoke') == 'AL2O3(L)' and rec.get('phase') == 'partE-' else ko)(
        'K0 the record of AL(L) makes %.10e kg of %s per kg burnt, born in %s' % (rec.get('nu_s', 0), rec.get('smoke'),
                                                                           rec.get('phase')))
    txt = open(log, errors='replace').read()
    (ok if re.search(r'IGLOO material AL\(L\) burns into\s+the smoke AL2O3\(L\) \(ICE family 1\)', txt) else ko)(
        'K0 the coupler reports the smoke and its ICE family')
    steps = [s for s in eb.read_report(log) if 'F' in s and s['t1'] > s['t0']]
    t_src = sum(s['t1'] - s['t0'] for s in steps)
    (ok if len(steps) >= 3 else ko)('K0 %d blocks with sources, t_src %.10e s' % (len(steps), t_src))
    if not steps:
        return
    HC = zones(inp + 'partL-properties.dat', 3)
    HO = zones(inp + 'partE-properties.dat', 3)['AL2O3(L)']
    rho_c = float(g2['rho'].split()[0])
    # the residue of each block: the end records (at time-end) of the residue parcels of the sweep that made its
    # sources (sweep k for the k-th block); its mass exact (xi nu_ox with the cap, (d-stop/d0)^3 at d-stop), the
    # printed d_end (six figures) checked against it
    res_mass = res_E = res_P = 0.0
    xi = float(g2['xi-cap']) if capped else 0.0
    # residue = off (the gate's variant "release"): the cap handed at burnout is born as smoke through the release
    # record of the oxide (nu_rel kg of smoke per kg of oxide): no residue, every kg of oxide a smoke
    released = capped and any(l.startswith('AL(L) cap ') and 'residue=off' in l for l in open(inp + 'partL-models.txt'))
    T_end = []
    worst_r = 0.0
    for k, s in enumerate(steps if not released else []):
        dt = s['t1'] - s['t0']
        if capped:
            name = 'partL-outloc-AL2O3(L)%s.dat' % ('' if k == 0 else '-sweep%d' % k)
            rho_r, Hr = float(g2['rho'].split()[1]), HO
        else:
            name = 'partL-outloc-AL(L)%s.dat' % ('' if k == 0 else '-sweep%d' % k)
            rho_r, Hr = rho_c, HC['AL(L)']
        end = outloc(run + '/OUTPUT/' + name)
        for i in range(n):
            Te, ue, de = end[i + 1]
            r = xi * rec['nu_ox'] if capped else (float(g2['d-stop']) / d0[i]) ** 3
            worst_r = max(worst_r, abs(rho_r * de ** 3 / (rho_c * d0[i] ** 3) / r - 1.0))
            T_end.append((dt, Te))
            res_mass += mdot[i] * r * dt
            res_E += mdot[i] * r * dt * (at(Hr, Te) + 0.5 * ue * ue)
            res_P += mdot[i] * r * dt * ue
    rbar = res_mass / (sum(mdot) * t_src)
    if released:
        nu_rel = [float(l.split()[7]) for l in open(inp + 'partL-products.txt')
                  if l.split()[:3] == ['material', 'AL2O3(L)', 'release']][0]
        ok('K0 residue = off: the cap handed at burnout is born as %.10e kg of smoke per kg of oxide' % nu_rel)
        # smoke per kg injected: nu_s (1 - xi) through the combustion record, xi nu_ox nu_rel through the release
        f = (1.0 - xi) + xi * rec['nu_ox'] * nu_rel / rec['nu_s']
        fg = 1.0
    else:
        (ok if worst_r <= 2e-5 else ko)('K0 residue: %.9f per kg injected (%s); the printed d_end gives it within %.2e '
                                        '(six figures: tol 2e-5)' % (rbar, 'xi nu_ox' if capped else '(d-stop/d0)^3',
                                                                     worst_r))
        f, fg = (1.0 - xi, 1.0) if capped else (1.0 - rbar, 1.0 - rbar)
    cols, vol = field(run + '/OUTPUT/partE-field.tec')
    Ms, Ps, Es, Ns, full = smoke_totals(cols, vol, '1', HO, rec['h298_s'])
    M_in = sum(mdot) * t_src
    # with the cap the coupler splits the metal through the record's nu_ox, IGLOO's cap law holds its own: they agree
    # within the coupler's check, 1e-7 (two tables of atomic weights), which moves the smoke by xi/(1 - xi) of it
    # (6.7e-10 measured); at d-stop IGLOO's step control ends the burn at d-stop to its own tolerance (8.6e-11 measured)
    judge('K1 smoke mass', Ms, rec['nu_s'] * f * M_in, 3e-8 if capped else 1e-9,
          'ICE field vs nu_s f SUM mdot t_src (f %.6f)' % f)
    bal = balance(txt, species)
    for sp in rec['g']:
        if sp in bal:
            judge('K2 %s' % sp, bal[sp], rec['g'][sp] * fg * M_in, TOL_BOOK if capped else 1e-6,
                  'the gas gained (closed-domain balance) vs g f_g SUM mdot t_src')
        else:
            ko('K2 %s: no row in the closed-domain balance' % sp)
    F_an = sum(a * b for a, b in zip(mdot, up)) * t_src
    dpg = [sum(s['p1'][c] - s['p0'][c] for s in eb.read_report(log) if 'p1' in s) for c in range(3)]
    judge('K3 momentum x', dpg[0] + Ps[0], F_an, res_P / F_an + 1e-9,
          'gas %.6e + smoke %.6e vs SUM mdot u_inj t_src (tol: the residue\'s |u_end| %.3e of it)' % (dpg[0], Ps[0],
                                                                                                    res_P / F_an))
    Q_an = sum(a * (at(HC['AL(L)'], t) + 0.5 * u * u) for a, u, t in zip(mdot, up, T0)) * t_src
    dEg = sum(s['E1'] - s['E0'] for s in eb.read_report(log) if 'E1' in s)
    # the floor: the round-off of the gas's total energy (MOSE sums it over the cells), against which the terms are
    # small differences: the residual is judged against the gas's total energy at the start (measured 4e-10 of it)
    E_gas = abs(eb.read_report(log)[0]['E0'])
    judge('K4 energy', dEg + Es + res_E, Q_an, TOL_E_GAS,
          'gas %.9e + smoke %.9e + residue %.9e J vs SUM mdot (H_c(T0) + u^2/2) t_src, of the gas\'s %.6e J'
          % (dEg, Es, res_E, E_gas), scale=E_gas)
    (ok if abs(Es) >= 1e3 * TOL_E_GAS * E_gas else ko)(
        'K4 non-vacuity: the smoke holds %.3e J, %.1e of the gas\'s energy (a coupler that left it in the gas fails K4)'
        % (abs(Es), abs(Es) / E_gas))
    Tm = sum(cols['rho_p1'][c] * vol[c] * cols['T_p1'][c] for c in full) / Ms
    if not capped:
        # the rule of C-IGLOO9: the state weighted by the metal burnt over the residence; a direct drop burns at
        # T-burn, so its smoke is born at T-burn; the residue parcels' time-end temperature (the residence's end,
        # where this box's drops stop at d-stop and stay) is what a birth at the residence's end would give
        Tb = float(g2['T-burn'])
        Tr = sum(dt * t for dt, t in T_end) / sum(dt for dt, t in T_end)
        judge('K5 birth temperature', Tm, Tb, 1e-6,
              'the smoke\'s mass-weighted temperature vs T-burn (burn-weighted); the residue parcels\' at time-end '
              '%.4f K' % Tr)
        (ok if abs(Tr - Tb) >= 1e3 * 1e-6 * Tb else ko)(
            'K5 non-vacuity: the residence\'s end (the residue at time-end, %.4f K) is %.3f K from T-burn' % (Tr, abs(Tr - Tb)))
    elif T_end:
        print('      the smoke\'s mass-weighted temperature %.4f K (the residue parcels end at %.4f K on average)' % (
            Tm, sum(dt * t for dt, t in T_end) / sum(dt for dt, t in T_end)))
    else:
        print('      the smoke\'s mass-weighted temperature %.4f K (no residue)' % Tm)
    rho_s = at(zones(inp + 'partE-properties.dat', 2)['AL2O3(L)'], Tm)
    ds = float(g2['d-smoke'])
    judge('K6 number', Ns, Ms / (math.pi / 6.0 * rho_s * ds ** 3), TOL_M, 'SUM n vol vs mass / (pi/6 rho_s d^3)')
    dl = [0.0] * 5
    for l in txt.split('\n'):
        if 'igloo-smoke-delivered' in l:
            dl = [a + float(b) for a, b in zip(dl, l.split()[2:])]
    judge('K7 hand-off', Ms, dl[0], TOL_M, 'ICE field mass vs the delivered lines')


def ice(run, log):
    inp = run + '/INPUT/'
    ini = eb.read_ini(run + '/input.ini')
    g2 = ini['GPB-Phase2']
    rec = record(inp + 'part-products.txt')
    species = gas_species(inp)
    (ok if rec.get('smoke') == 'AL2O3(L)' and rec.get('phase') == 'part-' else ko)(
        'K0 the record of AL(L) makes %.10e kg of %s per kg burnt, born in %s' % (rec.get('nu_s', 0), rec.get('smoke'),
                                                                           rec.get('phase')))
    txt = open(log, errors='replace').read()
    (ok if 'ICE material AL(L) burns into the smoke AL2O3(L) (family 2)' in txt else ko)(
        'K0 the coupler reports the smoke and its family')
    dt = float(re.search(r'ICE  \| Iter =\s+1 \| Time = \S+ \| Delta t = (\S+)', txt).group(1))
    HC = zones(inp + 'part-properties.dat', 3)
    RH = zones(inp + 'part-properties.dat', 2)
    a, vol = field(inp + 'part-ic.tec')
    b, _ = field(run + '/OUTPUT/part-field.tec')
    a = {k.replace('rp_', 'rho_p').replace('up_', 'u_p').replace('vp_', 'v_p').replace('wp_', 'w_p')
         .replace('Tp_', 'T_p').replace('np_', 'n_p').replace('11', '1').replace('21', '2'): v for k, v in a.items()}

    def tot(c, sfx, H, h298):
        r, T = c['rho_p' + sfx], c['T_p' + sfx]
        U = [c['u_p' + sfx], c['v_p' + sfx], c['w_p' + sfx]]
        M = sum(x * v for x, v in zip(r, vol))
        P = [sum(x * q * v for x, q, v in zip(r, Q, vol)) for Q in U]
        E = sum(x * v * (h298 + at(H, t) - at(H, T_REF) + 0.5 * sum(q[i] ** 2 for q in U))
                for i, (x, v, t) in enumerate(zip(r, vol, T)))
        return M, P, E
    M0, P0, E0 = tot(a, '1', HC['AL(L)'], rec['h298'])
    M1, P1, E1 = tot(b, '1', HC['AL(L)'], rec['h298'])
    Ms, Ps, Es, Ns, full = smoke_totals(b, vol, '2', HC['AL2O3(L)'], rec['h298_s'])
    lost = M0 - M1
    rep = [l.split() for l in txt.split('\n') if '[MI2-report] ice-to-gas ' in l]
    repm = [l.split() for l in txt.split('\n') if '[MI2-report] ice-to-gas-mass' in l]
    F = [float(x) for x in rep[0][2:5]]
    Q = float(rep[0][5])
    Mg = float(repm[0][2])
    print('      the step: dt %.6e s, aluminium lost %.10e kg, smoke %.10e kg; reported M %.10e kg/s, F %.6e N, Q %.6e W'
          % (dt, lost, Ms, Mg, F[0], Q))
    judge('B1 smoke mass', Ms, rec['nu_s'] * lost, TOL_M, 'family 2 vs nu_s x the aluminium family 1 lost')
    judge('B1 mass', lost, sum(rec['g'].values()) * lost + Ms, 2e-10,
          'lost vs the record\'s gas row (sum) x lost + smoke (the row prints eleven figures: 5e-11 measured)')
    judge('B1 gas share', Mg, sum(rec['g'].values()) * lost / dt, 2e-10,
          'the reported ice-to-gas-mass (the gas\'s net gain, steam consumed) vs SUM g x lost / dt')
    judge('B2 momentum', P0[0] - P1[0], F[0] * dt + Ps[0], 1e-9, 'family 1 lost vs F dt + smoke\'s')
    judge('B2 energy', E0 - E1, Q * dt + Es, TOL_E, 'family 1 lost (formation scale) vs Q dt + smoke\'s')
    up = float(ini['ice-box']['up'])
    Tp = float(ini['ICB-Block1']['T']) * float(ini['ICB-Block1']['kT'].split()[0])
    us = Ps[0] / Ms
    Tm = sum(b['rho_p2'][c] * vol[c] * b['T_p2'][c] for c in full) / Ms
    judge('B3 birth velocity', us, up, TOL_T, 'the smoke\'s mean u vs the drops\' initial u_p')
    judge('B3 birth temperature', Tm, Tp, TOL_T, 'the smoke\'s mean T vs the drops\' initial T_p')
    ds = float([t for t in open(inp + 'part-models.txt') if t.startswith('AL(L) combustion')][0].split('d-smoke=')[1].split()[0])
    worst = 0.0
    for c in full:
        rho_s = at(RH['AL2O3(L)'], b['T_p2'][c])
        worst = max(worst, abs(b['rho_p2'][c] / b['n_p2'][c] / (math.pi / 6.0 * rho_s * ds ** 3) - 1.0))
    (ok if full and worst <= 1e-12 else ko)('B4 number: rho_s / n_s = pi/6 rho_S d^3 in %d cells, worst rel.err %.2e' %
                                           (len(full), worst))
    t2 = sum(s['t1'] - s['t0'] for s in eb.read_report(log) if s.get('cycle') == 2)
    bal = balance(txt, species)
    for sp in rec['g']:
        if sp in bal:
            judge('G1 %s' % sp, bal[sp], rec['g'][sp] * lost / dt * t2, TOL_BOOK,
                  'the gas gained over cycle 2 vs g (lost/dt) t2')
        else:
            ko('G1 %s: no row in the closed-domain balance' % sp)


def main():
    mode, run, log = sys.argv[1], os.path.abspath(sys.argv[2]), sys.argv[3]
    {'igloo': igloo, 'ice': ice}[mode](run, log)
    print('== p30_al_smoke.py: %d passed, %d failed' % (npass, nfail))
    sys.exit(0 if nfail == 0 else 1)


if __name__ == '__main__':
    main()
