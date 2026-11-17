#!/usr/bin/env python3
"""igloo-smoke-box: the smoke of a burning IGLOO material is born in an ICE family, and the gas, the smoke and the
injected streams close mass, momentum and energy.

  python3 check.py RUN_DIR LOG [EVAP_BOX_MOM_DIR]

RUN_DIR holds input.ini, INPUT/ and OUTPUT/ of the run; LOG its log with the [MI2-report] lines; the third argument is
test/evap-box-mom, whose check.py gives the readers (default: ../evap-box-mom beside this file).

The case is test/igloo-burnout-box with a smoke: per kg of A burnt, nu_s kg of alumina AL2O3(a) (the tail of A's
combustion record in INPUT/partL-products.txt) born in the ICE family of partE, the rest N2. Cycle 1 runs the MOSE
blocks with IGLOO's sources; after each block the coupler adds the smoke of the block to ICE's state; cycle 2 is one ICE
step, whose field OUTPUT/partE-field.tec is the smoke at the end (drag = NoDrag, heat-transfer = NoHeat: the smoke
keeps what it was born with, and a step of a periodic box moves no mass, momentum or energy out of it). Every stream
burns out inside the box (d-burnout), so over the blocks with sources, t_src = SUM (t1 - t0):
  K1  smoke mass (ICE field)             = nu_s SUM mdot t_src
  K2  the gas's species (the coupler's closed-domain balance) = (1 - nu_s) SUM mdot t_src
  K3  gas momentum gained (MOSE, the report's gas-momentum lines) + smoke momentum (ICE field) = SUM mdot u_inj t_src
  K4  gas energy gained (MOSE, gas-energy lines) + smoke energy (ICE field, on the gas's formation scale: per kg
      h298_s + H_s(T) - H_s(T_ref) + |u|^2/2, H_s the alumina table of INPUT/partE-properties.dat, h298_s the record's)
      = SUM mdot (h_c(T0) + u_inj^2/2 + h298 - h_c(T_ref)) t_src, h_c the enthalpy of A (INPUT/partL-properties.dat)
  K5  the state the smoke is born with: C-IGLOO9 hands, per cell residence of a parcel, the mass it burnt there and its
      velocity and temperature weighted by the metal burnt over the residence; a parcel of this box burns out in the
      cell it is injected in, so the smoke is born at the burn-weighted mean of the parcel's life. Oracle: one parcel's
      life integrated here, outside IGLOO (IGLOO's model 5: Stokes drag, Ranz-Marshall Nu, the d^n law, beta-part of
      q-comb on the drop, constant cp; gas at rest at its initial temperature), against the smoke's temperature (the
      mass-weighted energy in ICE's field, inverted on the table) and its mean velocity (momentum over mass); the end
      of the life, what a birth at the residence's end would give, is the non-vacuity (90 K and 0.246 m/s away).
  K6  number: each drop born at d-smoke, SUM n vol = smoke mass / (pi/6 rho_s d_smoke^3)
  K7  hand-off: the field holds what the coupler reports it delivered (igloo-smoke-delivered lines)
With the hand-off in rates (DECISION 21: MOSE or ICE not time-accurate; the coupler says "reaches ICE as rates"), ICE's
smoke family integrates, from empty, the rates the coupler set before ICE's block over ICE's own steps; after K0:
  R0  rate lines, one ICE block (ice-time t0 t1), no amount delivered
  R1  smoke mass above ICE's floor = M_rate (t1 - t0)        R2 momentum = P_rate (t1 - t0)
  R3  energy (formation scale, as K4) = E_rate (t1 - t0)      R4 number above the floor = N_rate (t1 - t0)
  R5  non-vacuity: what the gas did not get over the blocks (igloo-smoke-withheld) is not what ICE holds
One PASS/FAIL line per check with the numbers and the tolerance; exit 1 on any FAIL."""
import importlib.util
import math
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
mom = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'evap-box-mom')
sys.path.insert(0, mom)
import check as eb          # test/evap-box-mom/check.py: read_ini, values, read_report

TOL_M = 1e-12       # K1, K6, K7 sums over the field against the closed form (round-off of the telescoped sums)
TOL_BOOK = 1e-6     # K2 the closed-domain balance prints seven figures
TOL_P = 1e-9        # K3 momentum closure, of SUM mdot u_inj t_src (measured 5.5e-15 along x, 1.8e-11 across: the
#                         gas's y momentum, 2e-13, with fifteen figures)
TOL_E = 1e-7        # K4 energy closure, of |SUM mdot (...) t_src| (measured 7e-10: MOSE's energy round trip, as
#                         test/igloo-burn-box's D2)
TOL_T = 1.0         # K5 kelvin: the gas warms by 0.1 K over the run and the parcels' end temperature by 0.3 K
#                         (measured 0.09 K against the burn-weighted mean: IGLOO's sums are a trapezoid on its steps)
TOL_U = 1e-3        # K5 velocity, of the burn-weighted mean (measured 5.8e-5: the gas the parcels push is not at rest)
TOL_R3 = 1e-12      # R3 the smoke's energy in ICE against the rate's (measured 7.2e-15: ICE's table read back here;
#                         no MOSE in it, unlike K4)
T_REF = 298.15
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


def table(path, col):
    """{T: column col} of an ATLAS table of integer-kelvin rows, and the column's name."""
    out, names = {}, []
    for line in open(path):
        if line.startswith('VARIABLES'):
            names = re.findall(r'"([^"]*)"', line)
            continue
        try:
            v = [float(x) for x in line.split()]
        except ValueError:
            continue
        if len(v) > col:
            out[v[0]] = v[col]
    return out, names


def at(rows, T):
    lo = float(int(T))
    return rows[lo] if T == lo else rows[lo] + (rows[lo + 1.0] - rows[lo]) * (T - lo)


def field(path):
    """Cell columns and cell volumes of an ICE Tecplot ASCII field (BLOCK, nodes then cell-centred variables)."""
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


def smoke_record(path, material):
    for line in open(path):
        f = line.split()
        if f[:3] == ['material', material, 'combustion'] and len(f) == 10:
            return f[6], float(f[7]), float(f[8]), f[9], float(f[5])
    return None


def life(ini, models, gas_mu, gas_k, gas_cp, Tg, d_end):
    """One parcel's life outside IGLOO (RK4, a fixed step): its end temperature and velocity, and the burn-weighted
    mean temperature and velocity."""
    g = ini['GPB-Phase2']
    rho, cp = float(g['rho']), float(g['cp'])
    K, X, N = models['K-burn'], models['X-eff'], models['n-burn']
    beta, q, Tign = models['beta-part'], models['q-comb'], models['T-ign']
    d0 = float(ini['IGLOO-BC']['diam'].split()[0])
    u0 = float(ini['IGLOO-BC']['up'].split()[0])
    T0 = float(ini['IGLOO-BC']['temp0'].split()[0])
    rhog = float(ini['ICB-Block1']['p']) / (8314.462618 / 28.014 * Tg)
    pr = gas_mu * gas_cp / gas_k

    def rhs(s):
        m, u, T = s
        d = (6.0 * m / (math.pi * rho)) ** (1.0 / 3.0)
        re_ = rhog * abs(u) * d / gas_mu
        nu = 2.0 + 0.6 * math.sqrt(re_) * pr ** (1.0 / 3.0)
        mdot = -(rho * math.pi * K * X / (2.0 * N)) * d ** (3.0 - N) if T >= Tign else 0.0
        return [mdot, -3.0 * math.pi * gas_mu * d * u / m,
                (nu * gas_k * math.pi * d * (Tg - T) - beta * q * mdot) / (m * cp)]

    m0, mend = math.pi / 6.0 * rho * d0 ** 3, math.pi / 6.0 * rho * d_end ** 3
    s, h, wT, wu = [m0, u0, T0], 1.0e-7, 0.0, 0.0
    while s[0] > mend:
        k1 = rhs(s)
        k2 = rhs([a + 0.5 * h * b for a, b in zip(s, k1)])
        k3 = rhs([a + 0.5 * h * b for a, b in zip(s, k2)])
        k4 = rhs([a + h * b for a, b in zip(s, k3)])
        n = [a + h / 6.0 * (b + 2 * c + 2 * e + f) for a, b, c, e, f in zip(s, k1, k2, k3, k4)]
        if n[0] < mend:
            fr = (s[0] - mend) / (s[0] - n[0])
            n = [a + fr * (b - a) for a, b in zip(s, n)]
        wT += (s[0] - n[0]) * 0.5 * (s[2] + n[2])
        wu += (s[0] - n[0]) * 0.5 * (s[1] + n[1])
        s = n
    wT += s[0] * s[2]
    wu += s[0] * s[1]
    return s[2], s[1], wT / m0, wu / m0


def rates(run, txt, inp, pE, h298s):
    """R0-R5: the hand-off in rates."""
    last, rate, tice, wsum, ndel = None, None, None, [0.0] * 5, 0
    for l in txt.split('\n'):
        if '[MI2-report] igloo-smoke-rate' in l:
            last = [float(x) for x in l.split()[2:]]
        elif '[MI2-report] ice-time' in l and last is not None:
            rate, tice = last, [float(x) for x in l.split()[2:4]]
        elif '[MI2-report] igloo-smoke-withheld' in l:
            wsum = [a + float(b) for a, b in zip(wsum, l.split()[2:])]
        elif '[MI2-report] igloo-smoke-delivered' in l:
            ndel += 1
    nice = txt.count('[MI2-report] ice-time')
    nrate = txt.count('[MI2-report] igloo-smoke-rate')
    good = rate is not None and nice == 1 and ndel == 0 and tice[1] > tice[0]
    (ok if good else ko)('R0 the hand-off is in rates: %d rate lines, %d ICE block, %d amounts delivered'
                         % (nrate, nice, ndel))
    if rate is None or tice is None:
        return
    dt = tice[1] - tice[0]
    cols, vol = field(os.path.join(run, 'OUTPUT', pE + '-field.tec'))
    r, T, nn = cols['rho_p1'], cols['T_p1'], cols['n_p1']
    # ICE's floor: the cells the smoke never reached (T = 0); none (a moving smoke reaches every cell): not
    # subtracted, its whole mass is 1e-20 kg/m3 in the box, 8e-15 of the smoke's
    rf = min([a for a, t in zip(r, T) if t < 1.0] or [0.0])
    nf = min([a for a, t in zip(nn, T) if t < 1.0] or [0.0])
    U = [cols['u_p1'], cols['v_p1'], cols['w_p1']]
    Ms = sum((a - rf) * v for a, v in zip(r, vol))
    Ps = [sum(a * c * v for a, c, v in zip(r, q, vol)) for q in U]
    HS, _ = table(os.path.join(inp, pE + '-properties.dat'), 3)
    Es = sum(a * v * (h298s + at(HS, t) - at(HS, T_REF) + 0.5 * sum(q[c] ** 2 for q in U))
             for c, (a, v, t) in enumerate(zip(r, vol, T)) if t >= 1.0)
    Ns = sum((a - nf) * v for a, v in zip(nn, vol))
    print('      ICE block %.10e s; floor %.3e kg/m3, %.3e 1/m3 (subtracted)' % (dt, rf, nf))
    judge('R1 smoke mass', Ms, rate[0] * dt, TOL_M, 'ICE field above the floor vs M_rate (t1 - t0)')
    judge('R2 momentum x', Ps[0], rate[1] * dt, TOL_M, 'ICE field vs P_rate,x (t1 - t0)')
    judge('R2 momentum y', Ps[1], rate[2] * dt, TOL_M, 'ICE field vs P_rate,y (t1 - t0)', scale=rate[1] * dt)
    judge('R3 energy', Es, rate[4] * dt, TOL_R3, 'ICE field (formation scale, as K4) vs E_rate (t1 - t0)')
    judge('R4 number', Ns, rate[5] * dt, TOL_M, 'ICE field above the floor vs N_rate (t1 - t0)')
    w = wsum[0]
    (ok if Ms > 0.0 and abs(Ms - w) > 1e3 * TOL_M * abs(Ms) else ko)(
        'R5 non-vacuity: the gas did not get %.6e kg over the blocks, %.4e times what ICE holds (an amount hand-off '
        'would hold it)' % (w, w / Ms if Ms > 0.0 else float('inf')))


def main():
    run, log = os.path.abspath(sys.argv[1]), sys.argv[2]
    inp = os.path.join(run, 'INPUT')
    ini = eb.read_ini(os.path.join(run, 'input.ini'))
    pL, pE = ini['IGLOO-General']['phase'], ini['ICE-Parameters']['phase']
    n = len(ini['IGLOO-BC']['y'].split())
    mdot = eb.values(ini, 'IGLOO-BC', 'mdot', n)
    up = eb.values(ini, 'IGLOO-BC', 'up', n)
    T0 = eb.values(ini, 'IGLOO-BC', 'temp0', n)
    rec = smoke_record(os.path.join(inp, pL + '-products.txt'), 'A')
    if rec is None:
        ko('K0 the combustion record of A in INPUT/%s-products.txt has no smoke tail' % pL)
        print('== check.py: %d passed, %d failed' % (npass, nfail))
        sys.exit(1)
    sname, nus, h298s, sphase, h298 = rec
    (ok if sphase == pE + '-' else ko)('K0 the record of A makes %.10e kg of %s per kg burnt, born in the phase %s'
                                        % (nus, sname, sphase))
    row = [l.split() for l in open(os.path.join(inp, pL + '-models.txt')) if l.startswith('A combustion')][0]
    models = {k: float(v) for k, v in (t.split('=') for t in row[3:])}
    (ok if models.get('d-smoke', 0.0) > 0.0 else ko)('K0 the row of A gives d-smoke %.3e m' % models.get('d-smoke', 0))
    txt = open(log, errors='replace').read()
    (ok if re.search(r'IGLOO material A burns into\s+the smoke %s \(ICE family 1\)' % re.escape(sname), txt) else ko)(
        'K0 the coupler reports the smoke and its ICE family')
    if 'reaches ICE as rates' in txt:
        rates(run, txt, inp, pE, h298s)
        print('== check.py: %d passed, %d failed' % (npass, nfail))
        sys.exit(0 if nfail == 0 else 1)

    steps = [s for s in eb.read_report(log) if 'F' in s and s['t1'] > s['t0']]
    t_src = sum(s['t1'] - s['t0'] for s in steps)
    M_in = sum(mdot) * t_src
    (ok if len(steps) >= 3 else ko)('K0 %d blocks with sources, t_src %.10e s' % (len(steps), t_src))

    cols, vol = field(os.path.join(run, 'OUTPUT', pE + '-field.tec'))
    r, T = cols['rho_p1'], cols['T_p1']
    # cells the smoke never reached hold ICE's floor density at T = 0, off the table: their mass is reported
    empty = [c for c, t in enumerate(T) if t < 1.0]
    m_empty = sum(r[c] * vol[c] for c in empty)
    T = [t if t >= 1.0 else T_REF for t in T]
    U = [cols['u_p1'], cols['v_p1'], cols['w_p1']]
    Ms = sum(a * v for a, v in zip(r, vol))
    Ps = [sum(a * c * v for a, c, v in zip(r, q, vol)) for q in U]
    HS, _ = table(os.path.join(inp, pE + '-properties.dat'), 3)
    Es = sum(a * v * (h298s + at(HS, t) - at(HS, T_REF) + 0.5 * sum(q[c] ** 2 for q in U))
             for c, (a, v, t) in enumerate(zip(r, vol, T)))
    Ns = sum(a * v for a, v in zip(cols['n_p1'], vol))

    judge('K1 smoke mass', Ms, nus * M_in, TOL_M, 'ICE field vs nu_s SUM mdot t_src')
    print('      %d cells hold no smoke (ICE\'s floor, %.3e kg in all, at T_REF in the sums)' % (len(empty), m_empty))

    m = re.search(r'closed-domain mass balance.*?\n.*?\n\s+1\s+(\S+)\s+(\S+)\s+(\S+)', txt)
    if not m:
        ko('K2 the log carries no closed-domain mass balance')
    else:
        judge('K2 gas species', float(m.group(2)), (1.0 - nus) * M_in, TOL_BOOK,
              'N2 the gas gained (closed-domain balance) vs (1 - nu_s) SUM mdot t_src')

    F_an = sum(a * b for a, b in zip(mdot, up)) * t_src
    dpg = [sum(s['p1'][c] - s['p0'][c] for s in eb.read_report(log) if 'p1' in s) for c in range(3)]
    judge('K3 momentum x', dpg[0] + Ps[0], F_an, TOL_P,
          'gas gained %.6e + smoke %.6e vs SUM mdot u_inj t_src' % (dpg[0], Ps[0]))
    judge('K3 momentum y', dpg[1] + Ps[1], 0.0, TOL_P, 'gas %.6e + smoke %.6e vs 0' % (dpg[1], Ps[1]), scale=F_an)

    hc, _ = table(os.path.join(inp, pL + '-properties.dat'), 3)
    Q_an = sum(a * (at(hc, t) + 0.5 * u * u + h298 - at(hc, T_REF)) for a, u, t in zip(mdot, up, T0)) * t_src
    dEg = sum(s['E1'] - s['E0'] for s in eb.read_report(log) if 'E1' in s)
    judge('K4 energy', dEg + Es, Q_an, TOL_E,
          'gas gained %.9e J + smoke %.9e J (formation scale) vs SUM mdot (h_c(T0) + u^2/2 + dh) t_src' % (dEg, Es))
    (ok if abs(Es) >= 0.5 * abs(Q_an) else ko)(
        'K4 non-vacuity: the smoke holds %.3f of the energy the streams bring (a coupler that left it in the gas '
        'would miss K4 by that much)' % (abs(Es) / abs(Q_an)))

    # K5: the smoke's mean energy per kg, inverted on the table, against one parcel's life
    e_mean = sum(a * v * at(HS, t) for a, v, t in zip(r, vol, T)) / Ms
    lo, hi = 200.0, 2000.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if at(HS, mid) < e_mean else (lo, mid)
    T_smoke = 0.5 * (lo + hi)
    MU, _ = table(os.path.join(inp, 'gas-transport.dat'), 1)
    KG, _ = table(os.path.join(inp, 'gas-transport.dat'), 2)
    CPG, _ = table(os.path.join(inp, 'gas-thermo.dat'), 1)
    Tg = float(ini['ICB-Block1']['T'])
    T_end, u_end, T_burn, u_burn = life(ini, models, at(MU, Tg), at(KG, Tg), at(CPG, Tg), Tg, models['d-burnout'])
    judge('K5 birth temperature', T_smoke, T_burn, TOL_T / T_burn,
          'the smoke in ICE (K) vs the burn-weighted mean of a parcel\'s life (integrated here; its end %.3f K)' % T_end)
    us = math.sqrt(sum(p * p for p in Ps)) / Ms
    judge('K5 birth velocity', us, u_burn, TOL_U,
          'the smoke\'s mean velocity in ICE (m/s) vs the burn-weighted mean of the life (its end %.3e m/s)' % u_end)
    (ok if abs(T_burn - T_end) >= 10.0 * TOL_T and abs(u_burn - u_end) >= 100.0 * TOL_U * u_burn else ko)(
        'K5 non-vacuity: the end of the life and its burn-weighted mean are %.1f K and %.3e m/s apart, the initial '
        'temperature %.1f K away' % (abs(T_end - T_burn), abs(u_burn - u_end), abs(T_end - T0[0])))

    rho_s = at(table(os.path.join(inp, pE + '-properties.dat'), 2)[0], T_smoke)
    judge('K6 number', Ns, Ms / (math.pi / 6.0 * rho_s * models['d-smoke'] ** 3), TOL_M,
          'SUM n vol vs smoke mass / (pi/6 rho_s d_smoke^3)')

    dl = [0.0] * 5
    for l in txt.split('\n'):
        if 'igloo-smoke-delivered' in l:
            dl = [a + float(b) for a, b in zip(dl, l.split()[2:])]
    judge('K7 hand-off', Ms, dl[0], TOL_M, 'ICE field mass vs the delivered lines (the coupler\'s own sums)')

    print('== check.py: %d passed, %d failed' % (npass, nfail))
    sys.exit(0 if nfail == 0 else 1)


if __name__ == '__main__':
    main()
