#!/usr/bin/env python3
"""ice-igloo-burn: one run, two burning phases on two models, the books closed per phase.

  python3 check.py RUN_DIR LOG [ICE_BOX_DIR]

RUN_DIR holds input.ini, INPUT/, MESH/ and OUTPUT/ of the run; LOG is its log with the [MI2-report] lines; the third
argument is test/MI2/cases/ice-exchange, whose check.py gives the readers of the Tecplot files and of the transport table (default:
../ice-exchange beside this file).

The Eulerian phase (ICE, [GPB-Phase3], combustion = Beckstead) is a uniform cloud at rest with a frozen source; the
Lagrangian phase (IGLOO, [GPB-Phase2], combustion = direct) is four steady streams that burn down to their residue.
Each phase burns into its own species of the gas (its record of INPUT/<phase>-products.txt, "formation"): the
Lagrangian one into the species of a-ar.yaml, the Eulerian one into that of a-he.yaml, so each species gains what
its phase burns and the carrier gains nothing. Closed forms, from input.ini and the INPUT files (V the box volume,
T_ref = 298.15 K, dh = h298 - h_c(T_ref) of each material, h298 from its record):
  ICE     mdot = rho pi K-burn X-eff/(2 n) d^(3-n)      M_E = n_p mdot V
          Qc   = n_p Nu pi d k_g (T_g - T_p) V          Q_E = M_E (c T_p + dh - beta q) - Qc
  IGLOO   M_L = SUM mdot (1 - r)                        Q_L = SUM mdot [h_c(T0) + u^2/2 - r h_c(T_end) + (1 - r) dh]
          r = (d_end/d0)^3, the end records of the sweep whose sources the report carries
and over a block [t0, t1] the gas energy rises by (Q_E + Q_L) (t1 - t0). W, the warming, does not read the report:
the first law of the rigid box with the tables' internal energies e_s = h_s - R_s T of the three species,
SUM_s m_s e_s(T1) = m0 e_carrier(T0) + E_in, solved for T1 and compared with the mass-weighted temperature the log
prints (the gas is not uniform in composition near the streams: the tolerance).
One PASS/FAIL line per check with the numbers and the tolerance; exit 1 on any FAIL."""
import importlib.util
import math
import os
import re
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
icebox = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'ice-exchange')
spec = importlib.util.spec_from_file_location('icebox_check', os.path.join(icebox, 'check.py'))
ib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ib)     # read_ini, read_tec, read_table, flint_interp, uniform, NU

TOL_ICE = 1e-8      # T1  ICE's reported mass and energy rates against the closed forms (measured 1e-14)
TOL_IGLOO = 1e-9    # T2  IGLOO's (the residue's end state is read with seven figures)
TOL_EN = 1e-5       # T3  E1 - E0 against (Q_E + Q_L) (t1 - t0) per block, as test/MI2/cases/igloo-burn
TOL_BOOK_M = 1e-11  # T4  the coupler's closed-domain mass balance
TOL_BOOK_E = 1e-8   # T4  the coupler's closed-domain energy bookkeeping
TOL_SUM = 1e-6      # T4  mass injected over the run (seven printed figures) against M_E t + M_L t
TOL_WARM = 1e-5     # W   warming against the first law (the log prints the mass-weighted temperature with six
#                         decimals, 2e-8 of the rise, of a gas whose composition is not uniform near the streams)
TOL_SPEC = 1e-6     # T5  mass a species gained (seven printed figures) against what its phase burnt
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


def judge(name, a, b, tol, what):
    e = abs(a - b) / abs(b)
    (ok if e <= tol else ko)('%s: %s %.16e vs %.16e, rel.err %.2e (tol %.0e)' % (name, what, a, b, e, tol))


def zones(path):
    """{zone name: {T: row}} of a Tecplot table of ATLAS (integer-kelvin rows), and its variable names."""
    out, cur, names = {}, None, []
    for line in open(path):
        if line.startswith('VARIABLES'):
            names = re.findall(r'"([^"]*)"', line)
        elif line.startswith('ZONE'):
            cur = out.setdefault(re.search(r'T="([^"]*)"', line).group(1), {})
        elif cur is not None:
            try:
                v = [float(x) for x in line.split()]
            except ValueError:
                continue
            if v:
                cur[v[0]] = v
    return out, names


def at(rows, T, col):
    lo = float(int(T))
    if T == lo:
        return rows[lo][col]
    return rows[lo][col] + (rows[lo + 1.0][col] - rows[lo][col]) * (T - lo)


def record(path, material):
    """The combustion record of a material in <phase>products.txt: (energy word, h298, {species: nu_mass})."""
    lines = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    for i, f in enumerate(lines):
        if f[0] == 'material' and f[1] == material and f[2] == 'combustion':
            return f[4], float(f[5]), {g[0]: float(g[1]) for g in lines[i + 1:i + 1 + int(f[3])]}
    return 'none', 0.0, {}


def values(sec, key, n):
    v = [float(x) for x in sec[key].split()]
    return v * n if len(v) == 1 else v


def main():
    run = os.path.abspath(sys.argv[1])
    log = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'run.log')
    cfg = ib.read_ini(os.path.join(run, 'input.ini'))
    gL, gE, bc = cfg['GPB-Phase2'], cfg['GPB-Phase3'], cfg['IGLOO-BC']
    pL, pE = gL['name'], gE['name']
    lineL = [l for l in open(os.path.join(run, 'INPUT', pL + '-phase.txt')) if l.startswith('A ')][0]
    lineE = [l for l in open(os.path.join(run, 'INPUT', pE + '-phase.txt')) if l.startswith('A ')][0]
    (ok if 'combustion=direct' in lineL and 'combustion=Beckstead' in lineE else ko)(
        'P  two phases, two models: %s (IGLOO) "%s", %s (ICE) "%s"' % (pL, lineL.split()[2], pE, lineE.split()[2]))

    # ---- the gas and the box
    gic = ib.read_tec(os.path.join(run, 'INPUT', 'gas-ic.tec'))
    rho0 = ib.uniform(sum(gic[k] for k in gic if k.startswith('rho')), 'gas IC rho')
    p0 = ib.uniform(gic['p'], 'gas IC p')
    mesh = ib.read_tec(os.path.join(run, 'MESH', 'mesh.tec'))
    V = float(np.prod([mesh[k].max() - mesh[k].min() for k in ('x', 'y', 'z')]))
    gtab, _ = zones(os.path.join(run, 'INPUT', 'gas-thermo.dat'))
    ttab, _ = zones(os.path.join(run, 'INPUT', 'gas-transport.dat'))
    sp = cfg['GPB-Phase1']['species'].split()[0]      # the carrier: the one species of the gas at the start
    order = [l.split()[0] for l in open(os.path.join(run, 'INPUT', 'gas-phase.txt')).read().split('\n')[1:] if l.strip()]
    Wm = {l.split()[0]: float(l.split()[1]) for l in open(os.path.join(run, 'INPUT', 'gas-phase.txt')).read().split('\n')[1:] if l.strip()}
    enL, h298L, nuL = record(os.path.join(run, 'INPUT', pL + '-products.txt'), 'A')
    enE, h298E, nuE = record(os.path.join(run, 'INPUT', pE + '-products.txt'), 'A')
    sL, sE = (list(nuL) + ['?'])[0], (list(nuE) + ['?'])[0]
    (ok if enL == enE == 'formation' and len(nuL) == len(nuE) == 1 and len({sp, sL, sE}) == 3 and sL in order and sE in order else ko)(
        'P  each phase burns into its own species of the gas: %s into %s, %s into %s, the carrier %s' % (pL, sL, pE, sE, sp))
    txt = open(log, errors='replace').read()
    g = re.search(r'mass-weighted gas T:\s+(\S+)\s+->\s+(\S+)', txt)
    h = re.search(r'mass-weighted frozen cv \[J/kg/K\] :\s+(\S+)\s+gas mass \[kg\] :\s+(\S+)', txt)
    if not g or not h:
        ko('the log carries no gas temperature diagnostic')
        print('== check.py: %d passed, %d failed' % (npass, nfail)); sys.exit(1)
    Tg0, Tg1 = float(g.group(1)), float(g.group(2))
    R = p0 / (rho0 * Tg0)
    k0 = at(ttab[sp], Tg0, 2)

    # ---- ICE: the frozen source of the uniform cloud
    pic = ib.read_tec(os.path.join(run, 'INPUT', pE + '-ic.tec'))
    rho_p = ib.uniform(pic['rp_11'], 'cloud rho_p')
    Tp = ib.uniform(pic['Tp_11'], 'cloud T_p')
    n_p = ib.uniform(pic['np_11'], 'cloud n_p')
    rest = all(float(np.abs(pic[k]).max()) == 0.0 for k in ('up_11', 'vp_11', 'wp_11'))
    rhoE, cE = float(gE['rho']), float(gE['cp'])
    Kb, nb, Xe = float(gE['K-burn']), float(gE.get('n-burn', 1.8)), float(gE.get('X-eff', 1.0))
    beta, qE, dresE, Tign = float(gE.get('beta-part', 0.0)), float(gE['q-comb']), float(gE['d-stop']), float(gE.get('T-ign', 2350.0))
    d = (6.0 * rho_p / (math.pi * rhoE * n_p)) ** (1.0 / 3.0)
    (ok if rest and Tp >= Tign and d > dresE else ko)('P  ICE\'s cloud is at rest, above T-ign (%g K >= %g K) and larger than its residue' % (Tp, Tign))
    mdot = rhoE * math.pi * Kb * Xe / (2.0 * nb) * d ** (3.0 - nb)
    M_E = n_p * mdot * V
    Qc = n_p * ib.NU * math.pi * d * k0 * (Tg0 - Tp) * V
    Q_E = M_E * (cE * Tp + (h298E - cE * T_REF) - beta * qE) - Qc

    # ---- the report
    steps, cur = [], None
    for line in txt.split('\n'):
        f = line.split()
        if line.startswith(' [MI2-report] cycle'):
            cur = {'cycle': int(f[2]), 't0': float(f[6]), 't1': float(f[8])}
            steps.append(cur)
        elif line.startswith(' [MI2-report] gas-energy'):
            cur['E0'], cur['E1'] = float(f[2]), float(f[3])
        elif line.startswith(' [MI2-report] ice-to-gas-mass'):
            cur['ME'] = float(f[2])
        elif line.startswith(' [MI2-report] ice-to-gas'):
            cur['QE'] = float(f[5])
        elif line.startswith(' [MI2-report] igloo-to-gas'):
            cur['FL'], cur['QL'], cur['ML'] = float(f[2]), float(f[5]), float(f[6])
    src = [s for s in steps if 'ME' in s]
    lag = [s for s in steps if 'ML' in s]
    if not src or not lag:
        ko('the report has no ice-to-gas-mass or no igloo-to-gas line: %d and %d steps' % (len(src), len(lag)))
        print('== check.py: %d passed, %d failed' % (npass, nfail)); sys.exit(1)

    # ---- IGLOO: the streams, with the end records of the sweep whose sources the last report carries
    n = len(bc['y'].split())
    md, up, T0, d0 = (values(bc, k, n) for k in ('mdot', 'up', 'temp0', 'diam'))
    ptab, pn = zones(os.path.join(run, 'INPUT', pL + '-properties.dat'))
    hcol = [i for i, v in enumerate(pn) if v.startswith('Enthalpy')][0]
    hc = ptab['A']
    name = '%s-outloc-A.dat' % pL if len(lag) <= 1 else '%s-outloc-A-sweep%d.dat' % (pL, len(lag) - 1)
    end = {}
    for line in open(os.path.join(run, 'OUTPUT', name)):
        f = line.split()
        if len(f) == 12:
            try:
                end[int(f[8])] = (float(f[3]), float(f[11]))
            except ValueError:
                continue
    r = [(end[i + 1][1] / d0[i]) ** 3 for i in range(n)]
    dhL = h298L - at(hc, T_REF, hcol)
    M_L = sum(m * (1.0 - x) for m, x in zip(md, r))
    Q_L = sum(m * (at(hc, t, hcol) + 0.5 * u * u - x * at(hc, end[i + 1][0], hcol) + (1.0 - x) * dhL)
              for i, (m, u, t, x) in enumerate(zip(md, up, T0, r)))
    print('      closed forms: ICE M %.15e kg/s, Q %.15e W (Qc %.9e W); IGLOO M %.15e kg/s, Q %.15e W'
          % (M_E, Q_E, Qc, M_L, Q_L))

    # ---- T1, T2: each solver against its closed form
    judge('T1 ICE mass', src[-1]['ME'], M_E, TOL_ICE, 'reported ice-to-gas-mass vs n mdot V')
    judge('T1 ICE energy', src[-1]['QE'], Q_E, TOL_ICE, 'reported ice-to-gas Q vs M (c T_p + dh - beta q) - Qc')
    (ok if len(set(s['ME'] for s in src)) == 1 and len(set(s['QE'] for s in src)) == 1 else ko)(
        'T1 ICE\'s source is frozen: one value of M and of Q in the %d reports' % len(src))
    judge('T2 IGLOO mass', lag[-1]['ML'], M_L, TOL_IGLOO, 'reported igloo-to-gas M vs SUM mdot (1 - r)')
    judge('T2 IGLOO energy', lag[-1]['QL'], Q_L, TOL_IGLOO, 'reported igloo-to-gas Q vs its closed form')
    judge('T2 IGLOO momentum', lag[-1]['FL'], sum(m * u for m, u in zip(md, up)), max(r), 'reported Fx vs SUM mdot u_inj (to within r)')

    # ---- T3: the gas gains the sum in every block
    worst = 0.0
    for s in src:
        dt = s['t1'] - s['t0']
        q = s['QE'] + s.get('QL', 0.0)
        worst = max(worst, abs((s['E1'] - s['E0']) - q * dt) / (abs(q) * dt))
    (ok if worst <= TOL_EN else ko)('T3 the gas energy rises by (Q_E + Q_L) (t1 - t0) in each of the %d blocks: worst %.2e (tol %.0e)'
                                    % (len(src), worst, TOL_EN))
    wl = max(abs((s['E1'] - s['E0']) - s['QE'] * (s['t1'] - s['t0'])) / (s['QL'] * (s['t1'] - s['t0'])) for s in lag)
    (ok if abs(wl - 1.0) <= 1e-3 else ko)('T3 non-vacuity: without IGLOO\'s share the residual is %.6f of Q_L (t1 - t0) in the worst block' % wl)

    # ---- T4: the coupler's books and the mass of the two phases together
    t_src = sum(s['t1'] - s['t0'] for s in src)
    t_lag = sum(s['t1'] - s['t0'] for s in lag)
    rows = {int(a): (float(b), float(c), float(d_)) for a, b, c, d_ in re.findall(
        r'^\s+(\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+)\s*$', txt, re.M)}
    idle = [int(a) for a in re.findall(r'^\s+(\d+)\s+\(untouched', txt, re.M)]
    e = re.search(r'rel\.err \(bookkeeping only, datum-blind\) :\s+(\S+)', txt)
    iL, iE, iC = (order.index(s) + 1 if s in order else 0 for s in (sL, sE, sp))
    if iL not in rows or iE not in rows or not e:
        ko('T4 the log carries no closed-domain mass balance of %s and %s, or no energy bookkeeping' % (sL, sE))
    else:
        werr = max(rows[iL][2], rows[iE][2])
        (ok if werr <= TOL_BOOK_M else ko)('T4 closed-domain mass balance of %s and %s: worst rel.err %.3e (tol %.0e)' % (sL, sE, werr, TOL_BOOK_M))
        (ok if float(e.group(1)) <= TOL_BOOK_E else ko)('T4 closed-domain energy bookkeeping: rel.err %.3e (tol %.0e)' % (float(e.group(1)), TOL_BOOK_E))
        judge('T5 species of the Lagrangian phase', rows[iL][1], M_L * t_lag, TOL_SPEC, 'mass %s gained vs M_L t' % sL)
        judge('T5 species of the Eulerian phase', rows[iE][1], M_E * t_src, TOL_SPEC, 'mass %s gained vs M_E t' % sE)
        (ok if idle == [iC] else ko)('T5 the carrier %s is untouched (rows untouched: %s)' % (sp, idle))

    # ---- W: the warming, by the first law of the rigid box
    m0 = rho0 * V
    dM = M_E * t_src + M_L * t_lag

    mass = {sp: m0, sL: M_L * t_lag, sE: M_E * t_src}

    def rise(internal):
        def e_gas(s, Tk):
            return at(gtab[s], Tk, 2) - (8314.46 / Wm[s] * Tk if internal else 0.0)
        target = m0 * e_gas(sp, Tg0) + Q_E * t_src + Q_L * t_lag
        lo, hi = Tg0 - 50.0, Tg0 + 1500.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if sum(mass[s] * e_gas(s, mid) for s in mass) < target else (lo, mid)
        return 0.5 * (lo + hi) - Tg0
    dT_v, dT_p = rise(True), rise(False)
    judge('W  warming', Tg1 - Tg0, dT_v, TOL_WARM, 'gas temperature rise of the log vs the first law of the rigid box (e = h - R T)')
    (ok if abs(dT_p - dT_v) >= 100.0 * TOL_WARM * dT_v else ko)(
        'W  non-vacuity: the same oracle at constant pressure (h in place of e) gives %.6f K, %.2f %% away' % (dT_p, 100.0 * abs(dT_p - dT_v) / dT_v))

    print('== check.py: %d passed, %d failed' % (npass, nfail))
    sys.exit(0 if nfail == 0 else 1)


if __name__ == '__main__':
    main()
