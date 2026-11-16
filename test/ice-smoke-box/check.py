#!/usr/bin/env python3
"""ice-smoke-box: a burning ICE material whose reaction makes a condensed product, the smoke, born in an ICE family of
the MK closure: the books of the three bodies (the gas's share the coupler hands, the burning cloud, the smoke) and
the smoke's birth against closed forms.

usage: check.py <run-dir> <log> [ice-box-dir] [birth|gas]
  run-dir      the run as run_case.sh lays it out: input.ini, INPUT/, MESH/, OUTPUT/
  log          its log, with the [MI2-report] lines of [MI2-Coupling] report = true and ICE's step lines
  ice-box-dir  test/ice-box, whose check.py gives the readers (default: ../ice-box beside this file)
  birth        the variant of one ICE step: the smoke against its birth alone
  gas          the variant of one ICE step whose sources stay frozen through a few MOSE iterations: the gas's gains

The material A (family 1) burns on the d^n law and burns out (d-burnout); per kg burnt the record of
INPUT/part-products.txt hands the gas 1 - nu_s kg (its species) and makes nu_s kg of the smoke S (its tail: S, nu_s,
h298_S, the phase), born in family 2 at d-smoke with the burning drop's velocity and temperature. ICE advances alone
(cycle 1); the report prints after ICE step j what the sources of that step hand the gas: M_j (mass, through the
record's species), F_j (force), Q_j (energy, on the gas's scale), the smoke already taken out; the removal stores its
gas share (ice-removed-stored, -delivered). With m0 = rho_p V, E0 = m0 (c T_p + |u_p|^2/2 + h298_A - c T_ref) and the
smoke at the end M_s = SUM rho_s V, P_s = SUM rho_s u_s V, E_s = SUM rho_s (c T_s + |u_s|^2/2) V + delta_s M_s with
delta_s = h298_S - c T_ref (ICE's energy of S is c T, A's table):
  B1 mass      SUM M dt + removal + M_s = m0,   and M_s = nu_s m0            (every kg of A burnt makes nu_s kg of S)
  B2 momentum  SUM F dt + removal + P_s = m0 u_p
  B3 energy    SUM Q dt + removal + E_s = E0
  B4 size      rho_s / n_s = pi/6 rho_S d_s^3 in every cell                  (each drop born at d-smoke keeps its mass)
"birth" (N-coupling 1: one forward-Euler step of A at its initial state, d-burnout below the drops):
  S1 the smoke's mass  M_s = nu_s M_b dt, M_b = n_p V rho pi K X d0^(3-n)/(2n)   (the d^n law at d0)
  S2 its velocity u_s = u_p, its temperature T_s = T_p                        (born with the burning drop's)
  S3 its drop mass rho_s / n_s = pi/6 rho_S d_s^3
"gas" (the same step, then MOSE through the frozen sources for a time t):
  G1 the gas gains (1 - nu_s) M_b t of mass: the record's species, without the smoke
  G2 its momentum and energy gain the reported F t and Q t (the smoke's taken out); the closed-domain mass balance
Prints PASS/FAIL lines with the numbers compared and the tolerance; exits 1 on any FAIL."""
import math
import os
import re
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
icebox = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'ice-box')
sys.path.insert(0, icebox)
import check as ib          # test/ice-box/check.py: read_ini, read_tec, uniform, judge, ok, ko, info

TOL_BOOK = 1e-11     # B1-B3: sums of some hundred terms printed with 17 figures
TOL_SIZE = 1e-12     # B4, S3: a ratio of two moments the births add in proportion
TOL_BIRTH = 1e-12    # S1, S2
TOL_GAS = 1e-8       # G1, G2: the gas's gains against the frozen rates (test/ice-burn-box's B1 and B2: 1e-10 measured)
T_REF = 298.15


def fnum(s):
    return float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', s))


def record(path):
    """The combustion record of A: h298, nu_s, h298 of the smoke."""
    for l in open(path):
        f = l.split()
        if f[:3] == ['material', 'A', 'combustion']:
            if len(f) < 10:
                sys.exit('check.py: the combustion record of A has no smoke tail')
            return float(f[5]), float(f[7]), float(f[8])
    sys.exit('check.py: no combustion record of A')


def main():
    run, log = os.path.abspath(sys.argv[1]), sys.argv[2]
    mode = sys.argv[4] if len(sys.argv) > 4 else 'books'
    birth = mode == 'birth'
    cfg = ib.read_ini(os.path.join(run, 'input.ini'))
    gpb, phys = cfg['GPB-Phase2'], cfg['ICE-Physics']
    rho_mat, c_mat = (float(v) for v in (phys['density'].split()[0], phys['specific-heat'].split()[0]))
    rho_S, c_S = (float(v) for v in (phys['density'].split()[1], phys['specific-heat'].split()[1]))
    KX, nb = float(gpb['K-burn']) * float(gpb.get('X-eff', 1.0)), float(gpb.get('n-burn', 1.8))
    row = [l.split() for l in open(os.path.join(run, 'INPUT', 'part-models.txt')) if l.startswith('A combustion ')][0]
    d_s = float([t for t in row if t.startswith('d-smoke=')][0].split('=')[1])
    h298, nu_s, h298_S = record(os.path.join(run, 'INPUT', 'part-products.txt'))
    m_s = math.pi / 6.0 * rho_S * d_s ** 3
    pic = ib.read_tec(os.path.join(run, 'INPUT', 'part-ic.tec'))
    rho_p = ib.uniform(pic['rp_11'], 'part IC rho_p')
    up = np.array([ib.uniform(pic[k], 'part IC ' + k) for k in ('up_11', 'vp_11', 'wp_11')])
    Tp = ib.uniform(pic['Tp_11'], 'part IC T_p')
    n_p = ib.uniform(pic['np_11'], 'part IC n_p')
    mesh = ib.read_tec(os.path.join(run, 'MESH', 'mesh.tec'))
    V = float(np.prod([mesh[k].max() - mesh[k].min() for k in ('x', 'y', 'z')]))
    m0 = rho_p * V
    E0 = m0 * (c_mat * Tp + 0.5 * float(np.dot(up, up)) + h298 - c_mat * T_REF)
    delta_s = h298_S - c_S * T_REF
    ib.info('A: m0 %.9e kg, nu_s %.6f, smoke S at d-smoke %.3e m (%.6e kg per drop), datum offset %.6e J/kg'
            % (m0, nu_s, d_s, m_s, delta_s))

    field = ib.read_tec(os.path.join(run, 'OUTPUT', 'part-field.tec'))
    rho_s, n_s = field['rho_p2'], field['n_p2']
    u_s = np.array([field[k] for k in ('u_p2', 'v_p2', 'w_p2')])
    T_s = field['T_p2']
    ncell = rho_s.size
    Vc = V / ncell
    M_s = float(rho_s.sum()) * Vc
    P_s = (rho_s * u_s).sum(axis=1) * Vc
    E_s = float((rho_s * (c_S * T_s + 0.5 * (u_s ** 2).sum(axis=0))).sum()) * Vc + delta_s * M_s
    size = rho_s / n_s
    worst = float(np.max(np.abs(size - m_s)) / m_s)

    d0 = (6.0 * m0 / (math.pi * rho_mat * n_p * V)) ** (1.0 / 3.0)
    Mb = n_p * V * rho_mat * math.pi * KX / (2.0 * nb) * d0 ** (3.0 - nb)
    if mode == 'gas':
        steps, _ = ib.parse_report(log)
        txt = open(log, errors='replace').read()
        last = [s for s in steps if 'F' in s]
        masses = sorted(set(re.findall(r'\[MI2-report\] ice-to-gas-mass\s+(\S+)', txt)))
        (ib.ok if len(masses) == 1 and len(last) >= 2 else ib.ko)('G  frozen source: %d distinct ice-to-gas-mass values in %d blocks' % (len(masses), len(last)))
        if len(masses) != 1 or len(last) < 2:
            print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail)); sys.exit(1)
        M, F, Q = fnum(masses[0]), last[0]['F'], last[0]['Q']
        t = last[-1]['t1'] - last[0]['t0']
        gic = ib.read_tec(os.path.join(run, 'INPUT', 'gas-ic.tec'))
        rho_gas0 = ib.uniform(sum(gic[k] for k in gic if k.startswith('rho')), 'gas IC rho')
        gf = ib.read_tec(os.path.join(run, 'OUTPUT', 'gas-field.tec'))
        rho1 = ib.uniform(sum(gf[k] for k in gf if k.startswith('rho')), 'final gas rho', tol=1e-10)
        gain = (rho1 - rho_gas0) * V
        ib.judge('G1 the report', M, (1.0 - nu_s) * Mb, TOL_BIRTH, 'reported M vs (1 - nu_s) M_b, the d^n law at d0')
        ib.judge('G1 gas mass', gain, (1.0 - nu_s) * Mb * t, TOL_GAS, 'gas mass gain (final field) vs (1 - nu_s) M_b t')
        ib.judge('G2 gas momentum', (last[-1]['p1'] - last[0]['p0'])[0], F[0] * t, TOL_GAS, 'gas x-momentum gain vs Fx t')
        ib.judge('G2 gas energy', last[-1]['E1'] - last[0]['E0'], Q * t, TOL_GAS, 'gas energy gain vs Q t')
        em = re.search(r'closed-domain mass balance.*?\n.*?\n\s+1\s+\S+\s+\S+\s+(\S+)', txt)
        try:
            eb = float(em.group(1)) if em else None
        except ValueError:
            eb = None
        (ib.ok if eb is not None and eb <= 1e-11 else ib.ko)('G2 closed-domain mass balance of species 1: rel.err %s (tol 1e-11)' % (em.group(1) if em else 'absent'))
        print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
        sys.exit(0 if ib.nfail == 0 else 1)
    if birth:
        dt = float(cfg['ICE-Numerics']['dt-max'])
        ib.judge('S1 the smoke\'s mass', M_s, nu_s * Mb * dt, TOL_BIRTH, 'SUM rho_s V vs nu_s M_b dt, the d^n law at d0')
        du = float(np.max(np.abs(u_s - up[:, None])))
        dT = float(np.max(np.abs(T_s - Tp)))
        (ib.ok if du <= TOL_BIRTH * max(abs(up[0]), 1.0) and dT <= TOL_BIRTH * Tp else ib.ko)(
            'S2 born with the burning drop\'s velocity and temperature: worst |u_s - u_p| %.1e m/s, |T_s - T_p| %.1e K '
            '(T_p %.3f K; tol %.0e relative)' % (du, dT, Tp, TOL_BIRTH))
        (ib.ok if worst <= TOL_SIZE else ib.ko)('S3 each smoke drop at d-smoke: rho_s/n_s vs pi/6 rho_S d_s^3 %.6e kg, worst %.1e (tol %.0e)'
                                                % (m_s, worst, TOL_SIZE))
        print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
        sys.exit(0 if ib.nfail == 0 else 1)

    dts, Ms, Fs, Qs, Ss = [], [], [], [], []
    dsum = np.zeros(5)
    for line in open(log, errors='replace'):
        m = re.match(r'ICE\s+\| Iter =\s+(\d+) \| Time = (\S+) \| Delta t = (\S+)', line)
        if m:
            if int(m.group(1)) == len(dts) + 1:
                dts.append(float(m.group(3)))
        elif line.startswith(' [MI2-report] ice-to-gas-mass'):
            Ms.append(fnum(line.split()[2]))
        elif line.startswith(' [MI2-report] ice-to-gas '):
            f = line.split()
            Fs.append(np.array([fnum(v) for v in f[2:5]])); Qs.append(fnum(f[5]))
        elif line.startswith(' [MI2-report] ice-removed-delivered '):
            dsum = dsum + np.array([fnum(v) for v in line.split()[2:7]])
        elif line.startswith(' [MI2-report] ice-removed-stored '):
            Ss.append(np.array([fnum(v) for v in line.split()[2:7]]) + dsum)
    n = len(Ms)
    (ib.ok if n > 0 and len(dts) == n and len(Fs) == n and len(Ss) == n else ib.ko)(
        'P  %d ICE steps, %d mass lines, %d source lines and %d store lines in the report' % (len(dts), n, len(Fs), len(Ss)))
    if n == 0 or len(dts) != n or len(Fs) != n or len(Ss) != n:
        print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail)); sys.exit(1)
    dt_max = float(cfg['ICE-Numerics']['dt-max'])
    worst_dt = max(abs(d - dt_max) / dt_max for d in dts)
    (ib.ok if worst_dt <= 1e-6 else ib.ko)('P  every one of the %d ICE steps is dt-max %.6e s: worst %.1e' % (n, dt_max, worst_dt))
    R = Ss[-1]
    sM = sum(M * dt_max for M in Ms) + R[0]
    sF = sum(F * dt_max for F in Fs) + R[1:4]
    sQ = sum(Q * dt_max for Q in Qs) + R[4]
    ib.info('the gas\'s share: mass %.9e kg (the removal %.3e), the smoke: %.9e kg' % (sM, R[0], M_s))
    ib.judge('B1 mass', sM + M_s, m0, TOL_BOOK, 'SUM M dt + removal + smoke vs rho_p V')
    ib.judge('B1 the smoke', M_s, nu_s * m0, TOL_BOOK, 'SUM rho_s V vs nu_s rho_p V')
    ib.judge('B2 momentum', sF[0] + P_s[0], m0 * up[0], TOL_BOOK, 'SUM Fx dt + removal + smoke vs rho_p V u_p')
    ib.judge('B3 energy', sQ + E_s, E0, TOL_BOOK, 'SUM Q dt + removal + smoke (with its offset) vs E0')
    (ib.ok if worst <= TOL_SIZE else ib.ko)('B4 each smoke drop keeps the mass it is born with: rho_s/n_s vs %.6e kg, worst %.1e (tol %.0e)'
                                            % (m_s, worst, TOL_SIZE))
    rmax = float(np.abs(field['rho_p1']).max())
    (ib.ok if rmax <= 1e-20 else ib.ko)('B5 the burning family ends empty: largest bulk density %.1e kg/m3' % rmax)
    print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
    sys.exit(0 if ib.nfail == 0 else 1)


if __name__ == '__main__':
    main()
