#!/usr/bin/env python3
"""ice-burn: the mass, the heat of reaction and the energy of a burning ICE material reach the gas.

usage: check.py <run-dir> [log] [ice-exchange-dir]
  run-dir      the directory of the run as run_case.sh lays it out: input.ini, INPUT/, MESH/, OUTPUT/
  log          the run's log, with the [MI2-report] lines of [MI2-Coupling] report = true
  ice-exchange-dir  the directory of ice-exchange, whose check.py gives the readers (default: ../ice-exchange)

The case is test/MI2/cases/ice-exchange with its one material burning from the first source, with the model its material line
names: combustion = direct (T_p equals boiling-temperature) or combustion = Beckstead (T_p above T-ign). Cycle 1 runs one forward-Euler ICE step, so the source ICE stores is the one of the initial
state; cycle 2 keeps it frozen through 200 coupling steps of one time-accurate MOSE iteration. Closed forms, from
input.ini and the INPUT files (V the box volume, d = 2 R_p ICE's diameter, c_pg = gamma R/(gamma - 1) the gas heat
capacity ICE forms, Nu = 2 the Stokes heat law):
  B      = (c_pg (T_g - T_p) + q) / L_v                          transfer number
  mdot   = pi d (k_g / c_pg) Nu ln(1 + B)                        burn rate of one drop [kg/s]
  M      = n_p mdot V                                            mass to the gas [kg/s]
  F      = (rho_p (u_p - u_g)/tau + n_p mdot u_p) V              force on the gas, with the momentum of the burnt mass [N]
  W      = rho_p u_p (u_p - u_g)/tau V                           drag work [W]
  dh     = h298 - c T_ref, T_ref = 298.15 K                      datum offset: h298 the material's enthalpy on the gas's
                                                                 scale, from the record of INPUT/part-products.txt [J/kg]
  Q      = W + M (c T_p + |u_p|^2/2 + dh)                        energy to the gas [W]
The heat of reaction Q_r = h298 - SUM nu h_s(T_ref) is inside dh: the gas receives the species of the record.
For combustion = Beckstead (the drop keeps beta q-comb per kg burnt and exchanges heat as an inert one):
  mdot   = rho pi K-burn X-eff / (2 n) d^(3-n)                   burn rate of one drop [kg/s]
  Qc     = n_p Nu pi d k_g (T_g - T_p) V                         convective heat the gas gives the drops [W]
  Q      = W + M (c T_p + |u_p|^2/2 + dh - beta q) - Qc          energy to the gas [W]
and over cycle 2, of duration t: the gas gains M t of mass, F t of momentum and Q t of energy.
B5, the warming, does not read the coupler's report: the box is rigid, so the gas's internal energy e = h - R T (its
table; c_v, not c_p) obeys (m0 + M t) (e(T1) + |u1|^2/2) = m0 e(T0) + Q t, solved for T1 and compared with the
temperature of the final field. Per kg burnt at T_p = T_ref the gas warms by Q_r + R T, not by Q_r.
Prints PASS/FAIL lines with the numbers compared and the tolerance; exits 1 on any FAIL."""
import math
import os
import re
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
icebox = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'ice-exchange')
sys.path.insert(0, icebox)
import check as ib          # test/MI2/cases/ice-exchange/check.py: read_ini, read_tec, read_table, flint_interp, uniform, parse_report, judge

TOL_M = 1e-8        # B1  reported mass rate against the closed form, and the gas mass gain against M t
TOL_Q = 1e-8        # B2  reported energy rate against the closed form
TOL_E = 1e-8        # B2  gas energy gain against Q t (total)
TOL_BOOK_M = 1e-11  # B3  the coupler's closed-domain mass balance (round-off 3e-15 with ifx 2023.1.0, 2e-13 with gfortran 12.2.1)
TOL_BOOK_E = 1e-8   # B3  the coupler's closed-domain energy bookkeeping
T_REF = 298.15
TOL_WARM = 1e-8     # B5  final gas temperature rise against the first law of the rigid box (measured 1e-10)
NON_VAC = 0.02      # B4  the gas gains at least this fraction of its mass, and the heat of reaction is at least
#                         this fraction of its internal energy


def thermo_at(path, T):
    """cp and h of the one species of <gas>-thermo.dat at T: linear between the integer-kelvin rows."""
    rows = {}
    for line in open(path):
        f = line.split()
        if len(f) == 5:
            try:
                v = [float(x) for x in f]
            except ValueError:
                continue
            rows[v[0]] = v
    lo = float(int(T))
    a, b = rows[lo], rows[lo + 1.0]
    return a[1] + (b[1] - a[1]) * (T - lo), a[2] + (b[2] - a[2]) * (T - lo)


def record(path, material):
    """The combustion record of a material in <phase>products.txt: (energy word, h298, {species: nu_mass})."""
    lines = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    for i, f in enumerate(lines):
        if f[0] == 'material' and f[1] == material and f[2] == 'combustion':
            return f[4], float(f[5]), {g[0]: float(g[1]) for g in lines[i + 1:i + 1 + int(f[3])]}
    ib.ko('P  no combustion record of %s in %s' % (material, path))
    return 'none', 0.0, {}


def log_number(txt, pattern):
    m = re.search(pattern, txt)
    return float(m.group(1)) if m else None


def main():
    run = os.path.abspath(sys.argv[1])
    log = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'run.log')
    cfg = ib.read_ini(os.path.join(run, 'input.ini'))
    rho_mat = float(cfg['GPB-Phase2']['rho'])
    c_mat = float(cfg['GPB-Phase2']['cp'])
    phys, gpb = cfg['ICE-Physics'], cfg['GPB-Phase2']

    def key(name, default=None):
        """A burning key: the material's ([GPB-Phase2], which GPB writes to the phase and models files), else [ICE-Physics]."""
        return gpb.get(name, phys.get(name, default))
    model = key('combustion')
    if model not in ('direct', 'Beckstead'):
        ib.ko('P  the case names no burning model (combustion = %s)' % model)
    q, d_res = float(key('q-comb')), float(key('d-stop'))
    beta = 0.0
    if model == 'direct':
        Lv, Tb = float(phys.get('latent-heat')), float(key('T-burn', phys.get('boiling-temperature')))
    else:
        Kb, nb, Xe = float(key('K-burn')), float(key('n-burn', 1.8)), float(key('X-eff', 1.0))
        Tb, beta = float(key('T-ign', 2350.0)), float(key('beta-part', 0.0))
    ncyc = len(cfg['HYDRA-Parameters']['MOSE-iter'].split())

    gic = ib.read_tec(os.path.join(run, 'INPUT', 'gas-ic.tec'))
    rho0 = ib.uniform(sum(gic[k] for k in gic if k.startswith('rho')), 'gas IC rho')
    p0 = ib.uniform(gic['p'], 'gas IC p')
    ug0 = np.array([ib.uniform(gic[k], 'gas IC ' + k) for k in ('u', 'v', 'w')])
    pic = ib.read_tec(os.path.join(run, 'INPUT', 'part-ic.tec'))
    rho_p = ib.uniform(pic['rp_11'], 'part IC rho_p')
    up = np.array([ib.uniform(pic[k], 'part IC ' + k) for k in ('up_11', 'vp_11', 'wp_11')])
    Tp = ib.uniform(pic['Tp_11'], 'part IC T_p')
    n_p = ib.uniform(pic['np_11'], 'part IC n_p')
    mesh = ib.read_tec(os.path.join(run, 'MESH', 'mesh.tec'))
    Vbox = float(np.prod([mesh[k].max() - mesh[k].min() for k in ('x', 'y', 'z')]))
    field = ib.read_tec(os.path.join(run, 'OUTPUT', 'gas-field.tec'))
    R = ib.uniform(field['R'], 'gas field R', tol=1e-13)
    Tg0 = p0 / (rho0 * R)
    Tmin, mu_tab, k_tab = ib.read_table(os.path.join(run, 'INPUT', 'gas-transport.dat'))
    mu0, k0 = ib.flint_interp(Tmin, mu_tab, Tg0), ib.flint_interp(Tmin, k_tab, Tg0)
    cp_g, _ = thermo_at(os.path.join(run, 'INPUT', 'gas-thermo.dat'), Tg0)
    _, h_ref = thermo_at(os.path.join(run, 'INPUT', 'gas-thermo.dat'), T_REF)
    gam = cp_g / (cp_g - R)
    cpg = gam * R / (gam - 1.0)
    (ib.ok if Tp >= Tb * (1.0 - 1e-10) else ib.ko)('P  the drops are at the burning temperature: T_p %g K, %s %g K'
                                                   % (Tp, 'T-burn' if model == 'direct' else 'T-ign', Tb))

    Rp = (0.75 * rho_p / (n_p * math.pi * rho_mat)) ** (1.0 / 3.0)
    d = 2.0 * Rp
    tau = rho_mat * d ** 2 / (18.0 * mu0)
    slip = up - ug0
    (ib.ok if d > d_res else ib.ko)('P  the drops are larger than the residue: d %.6e m, d-stop %.6e m' % (d, d_res))
    if model == 'direct':
        B = (cpg * (Tg0 - Tp) + q) / Lv
        mdot = math.pi * d * k0 / cpg * ib.NU * math.log1p(B)
        Qc = 0.0
    else:
        B = float('nan')
        mdot = rho_mat * math.pi * Kb * Xe / (2.0 * nb) * d ** (3.0 - nb)
        Qc = n_p * ib.NU * math.pi * d * k0 * (Tg0 - Tp) * Vbox
    M_an = n_p * mdot * Vbox
    F_an = (rho_p * slip / tau + n_p * mdot * up) * Vbox
    W_an = rho_p * float(np.dot(up, slip)) / tau * Vbox
    energy, h298, nu = record(os.path.join(run, 'INPUT', 'part-products.txt'), gpb.get('material', 'A'))
    (ib.ok if energy == 'formation' and list(nu) == [cfg['GPB-Phase1']['species'].split()[0]] else ib.ko)(
        'P  the record of the material is "%s" and gives the one species of the gas (%s)' % (energy, ' '.join(nu)))
    Qr = h298 - h_ref
    dh = h298 - c_mat * T_REF
    Q_an = W_an + M_an * (c_mat * Tp + 0.5 * float(np.dot(up, up)) + dh - beta * q) - Qc
    ib.info('gas: rho %.15e, T %.15e K, cp %.9e, k %.9e; drops: d %.9e m, n %.9e 1/m3, tau %.9e s' % (rho0, Tg0, cpg, k0, d, n_p, tau))
    ib.info('closed forms (%s): B %.12f, mdot %.12e kg/s, M %.15e kg/s, Fx %.15e N, W %.9e W, Qc %.9e W, dh %.9e J/kg, Q_r %.9e J/kg, Q %.15e W'
            % (model, B, mdot, M_an, F_an[0], W_an, Qc, dh, Qr, Q_an))

    txt = open(log, errors='replace').read()
    steps, _ = ib.parse_report(log)
    last = [s for s in steps if s['cycle'] == ncyc]
    masses = sorted(set(re.findall(r'\[MI2-report\] ice-to-gas-mass\s+(\S+)', txt)))
    if not masses:
        ib.ko('B1 the report carries no ice-to-gas-mass line: the coupler hands no ICE mass')
        M = 0.0
    else:
        (ib.ok if len(masses) == 1 else ib.ko)('B1 frozen source: %d distinct ice-to-gas-mass values' % len(masses))
        M = float(masses[0])
    F, Q = last[0]['F'], last[0]['Q']
    t = last[-1]['t1'] - last[0]['t0']
    dE = last[-1]['E1'] - last[0]['E0']
    dp_ = last[-1]['p1'] - last[0]['p0']
    rho1 = ib.uniform(sum(field[k] for k in field if k.startswith('rho')), 'final gas rho', tol=1e-10)
    dM = (rho1 - rho0) * Vbox
    ib.info('report: t %.15e s over %d steps, M %.15e kg/s, Fx %.15e N, Q %.15e W' % (t, len(last), M, F[0], Q))
    ib.info('gas: mass gain %.15e kg (%.3f of its mass), energy gain %.15e J' % (dM, dM / (rho0 * Vbox), dE))

    if M != 0.0:
        ib.judge('B1 mass', M, M_an, TOL_M, 'reported M vs n mdot V')
    ib.judge('B1 mass', dM, M_an * t, TOL_M, 'gas mass gain (final field) vs n mdot V t')
    ib.judge('B1 momentum', F[0], F_an[0], TOL_M, 'reported Fx vs (rho_p u_p/tau + n mdot u_p) V')
    ib.judge('B1 momentum', dp_[0], F_an[0] * t, TOL_M, 'gas x-momentum gain vs the closed form')
    ib.judge('B2 energy', Q, Q_an, TOL_Q, 'reported Q vs W + M (c T_p + u_p^2/2 + dh - beta q) - Qc')
    if M != 0.0:
        m1 = M / (n_p * Vbox)       # the measured rate of one drop, and the constant of its law
        law = 4.0 * m1 / (math.pi * rho_mat * d) if model == 'direct' else 2.0 * nb * m1 / (rho_mat * math.pi * d ** (3.0 - nb))
        ib.info('twin: law constant %.9e (Beckstead: K-burn X-eff [m^n/s]; direct: K [m2/s]) from the rate of one drop %.12e kg/s'
                % (law, m1))
        ib.info('twin: energy per kg burnt (Q - W + Qc)/M + beta q %.9f J/kg' % ((Q - W_an + Qc) / M + beta * q))
    ib.judge('B2 energy', dE, Q_an * t, TOL_E, 'gas energy gain vs the closed form')

    em = log_number(txt, r'closed-domain mass balance.*?\n.*?\n\s+1\s+\S+\s+\S+\s+(\S+)')
    ee = log_number(txt, r'rel\.err \(bookkeeping only, datum-blind\) :\s+(\S+)')
    if em is None or ee is None:
        ib.ko('B3 the log carries no closed-domain mass balance or energy bookkeeping of the run')
    else:
        (ib.ok if em <= TOL_BOOK_M else ib.ko)('B3 closed-domain mass balance of species 1: rel.err %.3e (tol %.0e)' % (em, TOL_BOOK_M))
        (ib.ok if ee <= TOL_BOOK_E else ib.ko)('B3 closed-domain energy bookkeeping: rel.err %.3e (tol %.0e)' % (ee, TOL_BOOK_E))

    cv = cp_g - R
    (ib.ok if M_an * t >= NON_VAC * rho0 * Vbox else ib.ko)('B4 non-vacuity: mass handed %.3e kg >= %g of the gas mass %.3e kg'
                                                            % (M_an * t, NON_VAC, rho0 * Vbox))
    (ib.ok if M_an * (Qr - beta * q) * t >= NON_VAC * rho0 * Vbox * cv * Tg0 else ib.ko)(
        'B4 non-vacuity: heat of reaction handed to the gas %.3e J >= %g of the gas internal energy %.3e J'
        % (M_an * (Qr - beta * q) * t, NON_VAC, rho0 * Vbox * cv * Tg0))

    # B5: the warming, written from the record and not from the coupler's offset: 1 kg of the material at T_ref
    # carries h298 on the gas's scale (the gas at T_ref plus the heat of reaction), so at T_p h298 + c (T_p - T_ref);
    # the frozen drops keep beta q and take Qc from the gas. The box is rigid: the gas's e = h - R T receives it
    tpath = os.path.join(run, 'INPUT', 'gas-thermo.dat')
    m0, m1 = rho0 * Vbox, rho0 * Vbox + M_an * t
    ke0 = 0.5 * float(np.dot(ug0, ug0))
    ke1 = 0.5 * float(np.dot(m0 * ug0 + F_an * t, m0 * ug0 + F_an * t)) / m1 ** 2
    e_in = h298 - beta * q + c_mat * (Tp - T_REF) + 0.5 * float(np.dot(up, up))

    def rise(internal):
        def e_gas(Tk):
            return thermo_at(tpath, Tk)[1] - (R * Tk if internal else 0.0)
        target = m0 * (e_gas(Tg0) + ke0) + (M_an * e_in + W_an - Qc) * t - m1 * ke1
        lo, hi = Tg0 - 100.0, Tg0 + 1500.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if m1 * e_gas(mid) < target else (lo, mid)
        return 0.5 * (lo + hi) - Tg0
    T1 = ib.uniform(field['T'], 'final gas T', tol=1e-10)
    dT_v, dT_p = rise(True), rise(False)
    ib.judge('B5 warming', T1 - Tg0, dT_v, TOL_WARM,
             'gas temperature rise of the final field vs the first law of the rigid box (e = h - R T; per kg burnt at T_ref: Q_r - beta q + R T = %.6e J/kg)'
             % (Qr - beta * q + R * T_REF))
    (ib.ok if abs(dT_p - dT_v) >= 100.0 * TOL_WARM * dT_v else ib.ko)(
        'B5 non-vacuity: the same oracle at constant pressure (h in place of e) gives %.6f K, %.2f %% away' % (dT_p, 100.0 * abs(dT_p - dT_v) / dT_v))

    print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
    sys.exit(0 if ib.nfail == 0 else 1)


if __name__ == '__main__':
    main()
