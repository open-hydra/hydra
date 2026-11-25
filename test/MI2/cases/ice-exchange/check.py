#!/usr/bin/env python3
"""ice-exchange: the ICE source stored in cycle 1 must reach the gas exactly during cycle 2.

usage: check.py <run-dir> [log]   (log: <run-dir>/run.log)
  run-dir  the directory of the run as run_case.sh lays it out: input.ini, INPUT/, MESH/, OUTPUT/
  log      the run's log, with the [MI2-report] lines of [MI2-Coupling] report = true

Oracle, closed forms from input.ini and the INPUT files, independent of both integrators and of
the time step (u_g, T_g: the gas state of cycle 1, at which ICE evaluated its source; V: box volume):
  R_p    = (3 rho_p / (4 pi n_p rho_mat))^(1/3)             ICE's radius from its own state
  tau    = rho_mat (2 R_p)^2 / (18 mu_g)                     Stokes drag, C_d = 24/Re
  F      = rho_p (u_p - u_g) / tau * V                       force on the gas [N]
  Q_heat = 2 Nu k_g pi R_p (T_g - T_p) n_p * V, Nu = 2       heat into the particles [W]
  W      = rho_p u_p (u_p - u_g) / tau * V                   drag work on the particles [W]
  over a MOSE block of duration t:  d(sum rho u V) = F t,  d(sum rho E V) = (W - Q_heat) t
mu_g and k_g are interpolated from INPUT/gas-transport.dat as FLINT does (linear between the
integer-kelvin rows) at T_g = p / (rho R) of INPUT/gas-ic.tec, R from the gas field.
Prints PASS/FAIL lines with the numbers compared and the tolerance; exits 1 on any FAIL."""
import configparser
import math
import os
import re
import sys

import numpy as np

TOL_TAU = 1e-10    # A0  reported tau against the closed form
TOL_MOM = 1e-8     # A1  momentum: gas gain against F t (report) and against the closed form
TOL_EN = 1e-8      # A2  energy: gas gain against Q t (report) and against the closed form
TOL_EN_STEP = 1e-7  # A2  the same per step: the block's gain is 1/200 of the total while the energy
#                     snapshots carry the same bounded noise (MOSE's primitives come out of its e -> T
#                     inversion; +-1e-11 J per snapshot measured at 1 thread, 2.2e-9 of a step's gain)
TOL_YZ = 1e-10     # A1  |dp_y|, |dp_z| relative to |dp_x|
TOL_UNI = 1e-10    # A4  spread of u over the cells relative to its mean
TOL_TAB = 1e-12    # S1  the table interpolation reproduces the field's mil and kl at its T
TOL_TAUSPREAD = 1e-12   # A0  taumax - taumin relative to tau (uniform state)
TOL_RP = 1e-12     # P   ICE's radius from the IC against dp/2 (make_case.py was applied)
NON_VAC_U = 0.05   # A3  gas velocity gain, fraction of u_p
NON_VAC_Q = 1e-6   # A3  heat exchanged over the block [J]
NU = 2.0           # Stokes heat: Nu = 2 (ICE Lib_Heat.f90)
RUNIV = 8314.51    # FLINT's universal gas constant [J/(kmol K)], cross-check of the field's R

npass = nfail = 0


def ok(msg):
    global npass
    npass += 1
    print('PASS  ' + msg)


def ko(msg):
    global nfail
    nfail += 1
    print('FAIL  ' + msg)


def info(msg):
    print('      ' + msg)


def rel(a, b):
    return abs(a - b) / abs(b)


def judge(name, a, b, tol, what):
    """PASS/FAIL on |a - b| / |b| <= tol, printing both numbers."""
    e = rel(a, b)
    (ok if e <= tol else ko)('%s: %s %.16e vs %.16e, rel.err %.2e (tol %.0e)' % (name, what, a, b, e, tol))
    return e


def read_ini(path):
    cfg = configparser.ConfigParser(inline_comment_prefixes=(';',), interpolation=None, strict=False)
    cfg.optionxform = str
    cfg.read(path)
    return cfg


def read_tec(path):
    """Tecplot BLOCK file, nodal x y z then cell-centred variables (names quoted or bare) -> dict name -> 1-D array."""
    txt = open(path).read()
    m = re.search(r'VARIABLES\s*=\s*(.*)', txt)
    names = [t.strip('"') for t in re.findall(r'"[^"]+"|\S+', m.group(1))]
    z = re.search(r'ZONE.*?I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', txt, re.S)
    I, J, K = (int(v) for v in z.groups())
    body = txt[z.end():]
    body = body[body.index('\n') + 1:]
    body = body.split('DATASETAUXDATA', 1)[0]     # the auxiliary data a solver appends after the zone are not field values
    vals = np.array([float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?', body.replace('D', 'E'))])
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    out, p = {}, 0
    for n in names:
        cnt = nn if n in ('x', 'y', 'z') else nc
        out[n] = vals[p:p + cnt]
        p += cnt
    if p != len(vals):
        raise ValueError('%s: %d values, %d expected from the header' % (path, len(vals), p))
    return out


def read_table(path):
    """FLINT transport table (one zone): rows T mu k on consecutive integer kelvins."""
    rows = []
    for line in open(path):
        t = line.split()
        if len(t) != 3:
            continue
        try:
            rows.append((float(t[0]), float(t[1]), float(t[2])))
        except ValueError:
            continue
    tab = np.array(rows)
    T = tab[:, 0]
    if not np.all(T == np.arange(T[0], T[0] + len(T))):
        raise ValueError('%s: rows are not consecutive integer kelvins' % path)
    return int(T[0]), tab[:, 1], tab[:, 2]


def flint_interp(Tmin, col, T):
    """FLINT co_k_mi_lam_Wilke, single species: linear between the rows idint(T) and idint(T)+1."""
    Ti = int(T)
    Tdiff = T - Ti
    Tmax = Tmin + len(col) - 1
    i1 = max(Tmin, min(Ti, Tmax))
    i2 = max(Tmin, min(Ti + 1, Tmax))
    return col[i1 - Tmin] + (col[i2 - Tmin] - col[i1 - Tmin]) * Tdiff


def uniform(a, name, tol=0.0):
    """The value of a uniform field; FAIL if it spreads more than tol (relative)."""
    lo, hi = float(a.min()), float(a.max())
    ref = max(abs(lo), abs(hi))
    if ref == 0.0 or (hi - lo) <= tol * ref:
        return 0.5 * (lo + hi)
    ko('%s is not uniform: %.16e .. %.16e' % (name, lo, hi))
    return 0.5 * (lo + hi)


def parse_report(log):
    """[MI2-report] lines -> list of dicts (one per coupling step) and the cycles the frozen notice names."""
    steps, frozen, cur = [], [], None
    for line in open(log, errors='replace'):
        m = re.match(r'\s*\[MI2-report\] (\S+)\s*(.*)$', line)
        if m:
            key, rest = m.group(1), m.group(2).split()
            if key == 'cycle':
                cur = {'cycle': int(rest[0]), 'step': int(rest[2]), 't0': float(rest[4]), 't1': float(rest[6])}
                steps.append(cur)
            elif cur is not None and key == 'gas-momentum':
                cur['p0'] = np.array([float(v) for v in rest[0:3]])
                cur['p1'] = np.array([float(v) for v in rest[3:6]])
            elif cur is not None and key == 'gas-energy':
                cur['E0'], cur['E1'] = float(rest[0]), float(rest[1])
            elif cur is not None and key == 'ice-to-gas':
                cur['F'] = np.array([float(v) for v in rest[0:3]])
                cur['Q'], cur['taumin'], cur['taumax'] = (float(v) for v in rest[3:6])
                cur['ice_line'] = m.group(2).strip()
            continue
        m = re.search(r'\[MI2\] ICE source FROZEN this cycle .*\[cycle (\d+), step (\d+),', line)
        if m:
            frozen.append((int(m.group(1)), int(m.group(2))))
    return steps, frozen


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    run = os.path.abspath(sys.argv[1])
    log = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'run.log')
    cfg = read_ini(os.path.join(run, 'input.ini'))
    rho_mat = float(cfg['GPB-Phase2']['rho'])
    dp = float(cfg['ICB-Block1']['dp'])
    hp = cfg['HYDRA-Parameters']
    ncyc = len(hp['MOSE-iter'].split())
    nsteps_ini = int(hp['N-coupling'].split()[-1])
    ice_iter_last = int(hp['ICE-iter'].split()[-1])

    # --- the initial state ICE saw (cycle 1), from the INPUT files ------------------------------------
    gic = read_tec(os.path.join(run, 'INPUT', 'gas-ic.tec'))
    rho_cols = [k for k in gic if k.startswith('rho')]
    rho0 = uniform(sum(gic[k] for k in rho_cols), 'gas IC rho')
    p0 = uniform(gic['p'], 'gas IC p')
    ug0 = np.array([uniform(gic[k], 'gas IC ' + k) for k in ('u', 'v', 'w')])
    pic = read_tec(os.path.join(run, 'INPUT', 'part-ic.tec'))
    rho_p = uniform(pic['rp_11'], 'part IC rho_p')
    up = np.array([uniform(pic[k], 'part IC ' + k) for k in ('up_11', 'vp_11', 'wp_11')])
    Tp = uniform(pic['Tp_11'], 'part IC T_p')
    n_p = uniform(pic['np_11'], 'part IC n_p')
    mesh = read_tec(os.path.join(run, 'MESH', 'mesh.tec'))
    Vbox = float(np.prod([mesh[k].max() - mesh[k].min() for k in ('x', 'y', 'z')]))
    ncell = pic['rp_11'].size

    field = read_tec(os.path.join(run, 'OUTPUT', 'gas-field.tec'))
    # R is computed in each cell from its composition: uniform to round-off, which the 17 digits of the file show
    R = uniform(field['R'], 'gas field R', tol=1e-13)
    W = float(open(os.path.join(run, 'INPUT', 'gas-phase.txt')).read().split('\n')[1].split()[1])
    info('gas constant R %.15e J/(kg K) from the field; Runiv/W with W %.6f: %.15e (rel.diff %.1e)'
         % (R, W, RUNIV / W, rel(R, RUNIV / W)))
    Tg0 = p0 / (rho0 * R)
    Tmin, mu_tab, k_tab = read_table(os.path.join(run, 'INPUT', 'gas-transport.dat'))
    mu0, k0 = flint_interp(Tmin, mu_tab, Tg0), flint_interp(Tmin, k_tab, Tg0)
    info('initial gas: rho %.15e kg/m3, p %.15e Pa, T %.15e K, u (%g %g %g) m/s' % (rho0, p0, Tg0, *ug0))
    info('mu(T) %.15e Pa s, k(T) %.15e W/(m K) from the transport table rows %d..%d'
         % (mu0, k0, Tmin, Tmin + len(mu_tab) - 1))
    info('condensed: rho_p %.15e kg/m3 (rho_p/rho %.6f), u_p (%g %g %g) m/s, T_p %g K, n_p %.15e 1/m3'
         % (rho_p, rho_p / rho0, *up, Tp, n_p))
    info('box: %d cells, V %.15e m3' % (ncell, Vbox))

    # S1: the interpolation rule reproduces MOSE's mil and kl from the final field's own T
    Tf = float(field['T'][0])
    judge('S1 transport table', flint_interp(Tmin, mu_tab, Tf), float(field['mil'][0]), TOL_TAB,
          'mu interpolated at the final T %.12f K vs the field mil' % Tf)
    judge('S1 transport table', flint_interp(Tmin, k_tab, Tf), float(field['kl'][0]), TOL_TAB,
          'k interpolated at the final T vs the field kl')

    # P: the IC make_case.py wrote (ICE's radius is dp/2)
    Rp = (0.75 * rho_p / (n_p * math.pi * rho_mat)) ** (1.0 / 3.0)
    judge('P  particle radius', 2.0 * Rp, dp, TOL_RP, 'ICE radius x2 from (rho_p, n_p, rho_mat) vs [ICB-Block1] dp')

    # --- the closed forms --------------------------------------------------------------------------
    tau_an = rho_mat * (2.0 * Rp) ** 2 / (18.0 * mu0)
    slip = up - ug0
    F_an = rho_p * slip / tau_an * Vbox
    Qheat_an = 2.0 * NU * k0 * math.pi * Rp * (Tg0 - Tp) * n_p * Vbox
    W_an = rho_p * float(np.dot(up, slip)) / tau_an * Vbox
    Q_an = W_an - Qheat_an
    info('closed forms: tau %.15e s, F (%.15e %g %g) N, Q_heat %.15e W, W_drag %.15e W, W - Q_heat %.15e W'
         % (tau_an, F_an[0], F_an[1], F_an[2], Qheat_an, W_an, Q_an))

    # --- the report of the last cycle ----------------------------------------------------------------
    steps, frozen = parse_report(log)
    last = [s for s in steps if s['cycle'] == ncyc]
    if ice_iter_last != 0:
        ko('the last cycle runs ICE (ICE-iter %d): the source is not frozen' % ice_iter_last)
    if len(last) != nsteps_ini or any(('F' not in s or 'E1' not in s or 'p1' not in s) for s in last):
        ko('report: cycle %d has %d complete steps, %d expected' % (ncyc, len(last), nsteps_ini))
        sys.exit(1)
    ok('report: cycle %d has all %d coupling steps' % (ncyc, len(last)))
    if any(c == ncyc and st == 1 for c, st in frozen):
        ok('ICE source FROZEN notice printed at cycle %d step 1' % ncyc)
    else:
        ko('no ICE source FROZEN notice at cycle %d step 1' % ncyc)
    lines = set(s['ice_line'] for s in last)
    (ok if len(lines) == 1 else ko)('S2 frozen source: %d distinct ice-to-gas lines over %d steps'
                                    % (len(lines), len(last)))
    if any(s['t1'] <= s['t0'] for s in last) or any(last[i + 1]['t0'] != last[i]['t1'] for i in range(len(last) - 1)):
        ko('the MOSE clock is not advancing continuously through the steps')
    F, Q, taumin, taumax = last[0]['F'], last[0]['Q'], last[0]['taumin'], last[0]['taumax']
    t = last[-1]['t1'] - last[0]['t0']
    dp_ = last[-1]['p1'] - last[0]['p0']
    dE = last[-1]['E1'] - last[0]['E0']
    info('report: t %.15e s over %d steps, F (%.15e %.3e %.3e) N, Q %.15e W, tau %.15e .. %.15e s'
         % (t, len(last), F[0], F[1], F[2], Q, taumin, taumax))
    info('report: gas momentum gain (%.15e %.3e %.3e) kg m/s, energy gain %.15e J (E0 %.15e J)'
         % (dp_[0], dp_[1], dp_[2], dE, last[0]['E0']))

    # A0 tau
    (ok if abs(taumax - taumin) <= TOL_TAUSPREAD * abs(taumax) else ko)(
        'A0 tau uniform: taumin %.16e, taumax %.16e (tol %.0e)' % (taumin, taumax, TOL_TAUSPREAD))
    judge('A0 tau', taumax, tau_an, TOL_TAU, 'reported tau vs rho_mat d^2 / (18 mu)')

    # A1 momentum
    judge('A1 momentum', dp_[0], F[0] * t, TOL_MOM, 'gas x-momentum gain vs reported Fx t')
    worst = max(rel(s['p1'][0] - s['p0'][0], s['F'][0] * (s['t1'] - s['t0'])) for s in last)
    (ok if worst <= TOL_MOM else ko)('A1 momentum: every step gains Fx dt, worst rel.err %.2e (tol %.0e)'
                                     % (worst, TOL_MOM))
    judge('A1 momentum', F[0], F_an[0], TOL_MOM, 'reported Fx vs rho_p u_p / tau V')
    judge('A1 momentum', dp_[0], F_an[0] * t, TOL_MOM, 'gas x-momentum gain vs rho_p u_p / tau V t')
    yz = max(abs(dp_[1]), abs(dp_[2])) / abs(dp_[0])
    (ok if yz <= TOL_YZ else ko)('A1 momentum: |dp_y|, |dp_z| <= %.0e |dp_x| (measured %.2e)' % (TOL_YZ, yz))
    (ok if F[1] == 0.0 and F[2] == 0.0 else ko)('A1 momentum: reported Fy, Fz are zero (%g, %g)' % (F[1], F[2]))

    # A2 energy
    judge('A2 energy', dE, Q * t, TOL_EN, 'gas energy gain vs reported Q t')
    worst = max(rel(s['E1'] - s['E0'], s['Q'] * (s['t1'] - s['t0'])) for s in last)
    (ok if worst <= TOL_EN_STEP else ko)('A2 energy: every step gains Q dt, worst rel.err %.2e (tol %.0e)'
                                         % (worst, TOL_EN_STEP))
    judge('A2 energy', Q, Q_an, TOL_EN, 'reported Q vs (W_drag - Q_heat)')
    judge('A2 energy', dE, Q_an * t, TOL_EN, 'gas energy gain vs (W_drag - Q_heat) t')

    # A3 non-vacuity
    u_final = field['u']
    du = float(u_final.mean()) - ug0[0]
    (ok if du >= NON_VAC_U * up[0] else ko)('A3 non-vacuity: gas velocity gain %.6e m/s >= %g u_p (%.3e)'
                                             % (du, NON_VAC_U, NON_VAC_U * up[0]))
    (ok if Qheat_an * t >= NON_VAC_Q else ko)('A3 non-vacuity: heat exchanged Q_heat t %.6e J >= %.0e J'
                                               % (Qheat_an * t, NON_VAC_Q))
    info('final gas: u %.12e m/s (momentum gain / (rho V) %.12e), T %.9f K (%+.4f K)'
         % (float(u_final.mean()), dp_[0] / (rho0 * Vbox), Tf, Tf - Tg0))

    # A4 uniformity of the final gas
    spread = float(u_final.max() - u_final.min()) / abs(float(u_final.mean()))
    (ok if spread <= TOL_UNI else ko)('A4 uniform: u spread over the cells %.2e relative (tol %.0e)' % (spread, TOL_UNI))
    for k in ('rho(1)', 'T', 'p'):
        if k in field:
            a = field[k]
            info('A4 spread of %s: %.2e relative' % (k, float(a.max() - a.min()) / abs(float(a.mean()))))

    print('== check.py: %d passed, %d failed' % (npass, nfail))
    sys.exit(0 if nfail == 0 else 1)


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        ko('check.py could not evaluate the run: %s: %s' % (type(e).__name__, e))
        print('== check.py: %d passed, %d failed' % (npass, nfail))
        sys.exit(1)
