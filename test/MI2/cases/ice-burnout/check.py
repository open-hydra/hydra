#!/usr/bin/env python3
"""ice-burnout: a burning ICE material without a residue burns out, and the books of the whole burn close.

usage: check.py <run-dir> [log] [ice-exchange-dir]
  run-dir      the directory of the run as run_case.sh lays it out: input.ini, INPUT/, MESH/, OUTPUT/
  log          the run's log, with the [MI2-report] lines of [MI2-Coupling] report = true and ICE's step lines
  ice-exchange-dir  the directory of ice-exchange, whose check.py gives the readers (default: ../ice-exchange)

Cycle 1 advances ICE alone, one forward-Euler step per coupling step; the gas does not move. After ICE step j the
report prints what the source of that step hands the gas: the mass rate M_j, the force F_j, the energy rate Q_j.
With dt_j the step (dt-max when the case gives it: every printed step must equal it; else the six printed figures),
V the box volume, m0 = rho_p V the mass of the cloud, c its heat capacity, T_ref = 298.15 K and h298 the material's
enthalpy on the gas's scale (the record of INPUT/part-products.txt):
  SUM_j M_j dt_j = m0                                      the mass of the cloud                    [kg]
  SUM_j F_j dt_j = m0 u_p                                  its momentum                             [kg m/s]
  SUM_j Q_j dt_j = m0 (c T_p + |u_p|^2/2 + h298 - c T_ref) its energy on the gas's scale            [J]
whatever the burning law, the drag and the heat exchange do on the way: the cloud ends empty, so what it held is
what its sources handed, with what the removal took out: the amounts ICE stores for the next MOSE block (the line
ice-removed-stored) and what a MOSE block has handed of them (ice-removed-delivered). The removal: the step that
first stores an amount hands, beside what it burns, all that the cloud still held. It is the first step that leaves a drop at or below the mass of a drop of
d-burnout (the size rule), or that leaves less than it burnt (the rule of a step too long for the end of the burn);
after it the source is zero, the amount removed (stored or delivered) does not change and ICE's field is empty.
The word "frozen" as fourth argument checks the gas side of the removal on a variant of the case: d-burnout above
the diameter of the drops and one ICE step, so that the removal fires in that step: its source is the burn alone,
M = n_p V rho pi K X d0^(3-n)/(2n) (the d^n law at d0), and the rest of the cloud is stored, m0 - M dt; cycle 2 keeps
the source frozen through a few MOSE iterations of total duration t and hands the store once, in its first block:
the gas gains M t + (m0 - M dt), Fx t + the stored momentum, Q t + the stored energy (the inter-phase books).
"frozen2" is the same with two ICE steps in the coupling step (ICE-iter 2): the second finds the cloud empty, its
source is zero; the gas takes the block's rate, the steps' sources weighted by their dt, M dt/(2 dt) = M/2, over its
time t, and what the first step stored: it gains m0 + M t/2 - M dt.
"cold" checks a variant whose drops never burn (T-ign above any temperature they reach) with d-burnout above them:
no step hands mass, nothing is stored, the cloud keeps its mass.
Prints PASS/FAIL lines with the numbers compared and the tolerance; exits 1 on any FAIL."""
import math
import os
import re
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
icebox = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'ice-exchange')
sys.path.insert(0, icebox)
import check as ib          # test/MI2/cases/ice-exchange/check.py: read_ini, read_tec, uniform, judge, ok, ko, info

TOL_BOOK = 1e-11     # R1-R3 with dt-max: sums of some hundred terms printed with 17 figures (measured 1e-15 to 1e-13)
TOL_BOOK_DT = 2e-6   # R1-R3 without dt-max: ICE prints its step with six figures (5e-7)
TOL_STEP = 1e-6      # every printed step against dt-max (six printed figures)
NON_VAC = 1e-4       # R5  the removal hands at least this fraction of the cloud's mass
TOL_GAS = 1e-8       # F2  the gas's gains against the frozen rates (test/MI2/cases/ice-burn's B1 and B2: 1e-10 measured)
TOL_RATE = 1e-12     # F1  the burn rate of the step against the d^n law at d0 (one forward-Euler step)
TOL_COLD = 1e-13     # C1  the cold cloud's bulk density against the injected one
T_REF = 298.15


def fnum(s):
    """A number of the report; Fortran prints a three-digit exponent without its E."""
    return float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', s))


def record(path, material):
    lines = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    for i, f in enumerate(lines):
        if f[0] == 'material' and f[1] == material and f[2] == 'combustion':
            return f[4], float(f[5])
    return 'none', 0.0


def lines5(txt, key):
    """The five numbers (mass, momentum x y z, energy) of each report line of the key, in the order of the log."""
    return [np.array([fnum(v) for v in l.split()[2:7]]) for l in txt.splitlines() if l.startswith(' [MI2-report] ' + key + ' ')]


def frozen(run, log, m0, up, E0, dt, rho_gas0, V, Mb, iters):
    """The removal of a frozen cycle: the step's source against the burn, the store against the rest of the cloud,
    and the gas's gains against the frozen source over the MOSE time plus the store handed once."""
    steps, _ = ib.parse_report(log)
    txt = open(log, errors='replace').read()
    last = [s for s in steps if 'F' in s]
    masses = sorted(set(re.findall(r'\[MI2-report\] ice-to-gas-mass\s+(\S+)', txt)))
    (ib.ok if len(masses) == 1 and len(last) >= 2 else ib.ko)('F1 frozen source: %d distinct ice-to-gas-mass values in %d blocks' % (len(masses), len(last)))
    if len(masses) != 1 or not last:
        return
    M, F, Q = fnum(masses[0]), last[0]['F'], last[0]['Q']
    t = last[-1]['t1'] - last[0]['t0']
    deliv = lines5(txt, 'ice-removed-delivered')
    nz = [d for d in deliv if d[0] != 0.0]
    (ib.ok if len(nz) == 1 and deliv and deliv[0][0] != 0.0 else ib.ko)(
        'F1 the store reaches the gas once, in the first MOSE block after ICE: %d of %d blocks deliver' % (len(nz), len(deliv)))
    D = nz[0] if nz else np.zeros(5)
    if iters == 1:
        ib.judge('F1 burn rate', M, Mb, TOL_RATE, 'reported M vs the d^n law at d0, n_p V rho pi K X d0^(3-n)/(2n)')
        ib.judge('F1 mass', M * dt + D[0], m0, TOL_BOOK, 'the step\'s burn M dt and the store vs rho_p V')
        ib.judge('F1 momentum', F[0] * dt + D[1], m0 * up[0], TOL_BOOK, 'Fx dt and the stored momentum vs rho_p V u_p')
        ib.judge('F1 energy', Q * dt + D[4], E0, TOL_BOOK, 'Q dt and the stored energy vs rho_p V (c T_p + u_p^2/2 + h298 - c T_ref)')
        whole = m0 + Mb * (t - dt)
    else:
        # the block's rate: the first step burns M_b over dt, the second (the cloud empty) nothing, over the same dt
        ib.judge('F1 the block\'s rate', M, 0.5 * Mb, TOL_RATE,
                 'reported M vs M_b dt / (2 dt), the mean of the two ICE steps (the second finds the cloud empty)')
        ib.judge('F1 mass', D[0], m0 - Mb * dt, TOL_RATE, 'the store vs rho_p V less the first step\'s burn, the d^n law at d0')
        whole = m0 + 0.5 * Mb * t - Mb * dt
    field = ib.read_tec(os.path.join(run, 'OUTPUT', 'gas-field.tec'))
    rho1 = ib.uniform(sum(field[k] for k in field if k.startswith('rho')), 'final gas rho', tol=1e-10)
    gain = (rho1 - rho_gas0) * V
    ib.judge('F2 gas mass', gain, M * t + D[0], TOL_GAS, 'gas mass gain (final field) vs M t + the store')
    ib.judge('F2 gas momentum', (last[-1]['p1'] - last[0]['p0'])[0], F[0] * t + D[1], TOL_GAS, 'gas x-momentum gain vs Fx t + the store')
    ib.judge('F2 gas energy', last[-1]['E1'] - last[0]['E0'], Q * t + D[4], TOL_GAS, 'gas energy gain vs Q t + the store')
    ib.judge('F2 the whole cloud', gain, whole, TOL_GAS, 'gas mass gain vs %s' % (
        'rho_p V + M_b (t - dt), M_b the d^n law at d0' if iters == 1 else 'rho_p V + M_b t/2 - M_b dt'))
    em = re.search(r'closed-domain mass balance.*?\n.*?\n\s+1\s+\S+\s+\S+\s+(\S+)', txt)
    try:
        eb = float(em.group(1)) if em else None
    except ValueError:
        eb = None
    (ib.ok if eb is not None and eb <= 1e-11 else ib.ko)('F2 closed-domain mass balance of species 1: rel.err %s (tol 1e-11)' % (em.group(1) if em else 'absent'))
    (ib.ok if D[0] >= 0.5 * m0 and t != dt else ib.ko)('F3 non-vacuity: the store is %.3f of the cloud, the MOSE time %.3f of the ICE step'
                                                       % (D[0] / m0, t / dt))


def cold(run, log, m0, V, rho_p):
    """No drop burns: no mass handed, nothing stored, the cloud keeps its mass."""
    txt = open(log, errors='replace').read()
    Ms = [fnum(v) for v in re.findall(r'\[MI2-report\] ice-to-gas-mass\s+(\S+)', txt)]
    Ss = lines5(txt, 'ice-removed-stored')
    (ib.ok if Ms and all(M == 0.0 for M in Ms) else ib.ko)('C1 no step hands mass: %d steps, largest %.1e kg/s' % (len(Ms), max([abs(M) for M in Ms] or [0.0])))
    (ib.ok if Ss and all(s[0] == 0.0 for s in Ss) else ib.ko)('C1 nothing is stored: %d steps, largest %.1e kg' % (len(Ss), max([abs(s[0]) for s in Ss] or [0.0])))
    field = ib.read_tec(os.path.join(run, 'OUTPUT', 'part-field.tec'))
    rho = sum(field[k] for k in field if k.startswith('rho_p'))
    ib.judge('C1 the cloud keeps its mass', float(np.mean(rho)), rho_p, TOL_COLD, 'mean bulk density at the end vs the injected one')


def main():
    run = os.path.abspath(sys.argv[1])
    log = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'run.log')
    cfg = ib.read_ini(os.path.join(run, 'input.ini'))
    gpb, num = cfg['GPB-Phase2'], cfg['ICE-Numerics']
    rho_mat, c_mat = float(gpb['rho']), float(gpb['cp'])
    d_b = float(gpb.get('d-burnout', 0.0))
    m_b = math.pi / 6.0 * rho_mat * d_b ** 3 if d_b > 0.0 else 1e-20
    KX, nb = float(gpb['K-burn']) * float(gpb.get('X-eff', 1.0)), float(gpb.get('n-burn', 1.8))
    dt_max = float(num['dt-max']) if 'dt-max' in num else None
    (ib.ok if gpb.get('combustion') == 'Beckstead' and 'd-stop' not in gpb else ib.ko)(
        'P  the material burns on the d^n law and leaves no residue (no d-stop); removal at a drop mass of %.6e kg' % m_b)

    pic = ib.read_tec(os.path.join(run, 'INPUT', 'part-ic.tec'))
    rho_p = ib.uniform(pic['rp_11'], 'part IC rho_p')
    up = np.array([ib.uniform(pic[k], 'part IC ' + k) for k in ('up_11', 'vp_11', 'wp_11')])
    Tp = ib.uniform(pic['Tp_11'], 'part IC T_p')
    n_p = ib.uniform(pic['np_11'], 'part IC n_p')
    mesh = ib.read_tec(os.path.join(run, 'MESH', 'mesh.tec'))
    V = float(np.prod([mesh[k].max() - mesh[k].min() for k in ('x', 'y', 'z')]))
    energy, h298 = record(os.path.join(run, 'INPUT', 'part-products.txt'), gpb.get('material', 'A'))
    m0 = rho_p * V
    E0 = m0 * (c_mat * Tp + 0.5 * float(np.dot(up, up)) + h298 - c_mat * T_REF)

    mode = sys.argv[4] if len(sys.argv) > 4 else 'books'
    if mode in ('frozen', 'frozen2'):
        gic = ib.read_tec(os.path.join(run, 'INPUT', 'gas-ic.tec'))
        rho_gas0 = ib.uniform(sum(gic[k] for k in gic if k.startswith('rho')), 'gas IC rho')
        (ib.ok if dt_max is not None and m_b >= m0 / (n_p * V) else ib.ko)(
            'P  %s: d-burnout is above the drops (%.6e kg against %.6e kg per drop) and the step is dt-max' % (mode, m_b, m0 / (n_p * V)))
        d0 = (6.0 * m0 / (math.pi * rho_mat * n_p * V)) ** (1.0 / 3.0)
        Mb = n_p * V * rho_mat * math.pi * KX / (2.0 * nb) * d0 ** (3.0 - nb)
        frozen(run, log, m0, up, E0, dt_max, rho_gas0, V, Mb, 1 if mode == 'frozen' else 2)
        print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
        sys.exit(0 if ib.nfail == 0 else 1)
    if mode == 'cold':
        T_ign = float(gpb['T-ign'])
        (ib.ok if m_b >= m0 / (n_p * V) and T_ign > Tp else ib.ko)(
            'P  cold: d-burnout above the drops (%.6e kg against %.6e kg per drop), T-ign %.1f K above their %.1f K'
            % (m_b, m0 / (n_p * V), T_ign, Tp))
        cold(run, log, m0, V, rho_p)
        print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
        sys.exit(0 if ib.nfail == 0 else 1)

    # ICE's steps and the report lines that follow them, in the order of the log
    dts, Ms, Fs, Qs, Ss = [], [], [], [], []
    dsum = np.zeros(5)   # what the MOSE blocks have delivered of the store so far
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
            Ss.append(np.array([fnum(v) for v in line.split()[2:7]]) + dsum)      # the amount removed so far
    n = len(Ms)
    (ib.ok if n > 0 and len(dts) == n and len(Fs) == n and len(Ss) == n else ib.ko)(
        'P  %d ICE steps, %d mass lines, %d source lines and %d store lines in the report' % (len(dts), n, len(Fs), len(Ss)))
    if n == 0 or len(dts) != n or len(Fs) != n or len(Ss) != n:
        print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail)); sys.exit(1)
    tol = TOL_BOOK
    if dt_max is not None:
        worst = max(abs(d - dt_max) / dt_max for d in dts)
        (ib.ok if worst <= TOL_STEP else ib.ko)('P  every one of the %d ICE steps is dt-max %.6e s: worst %.1e (tol %.0e)' % (n, dt_max, worst, TOL_STEP))
        dts = [dt_max] * n
    else:
        tol = TOL_BOOK_DT
        ib.info('no dt-max: the steps are the six printed figures of ICE (%.6e to %.6e s)' % (min(dts), max(dts)))

    S = Ss[-1]           # what the removal took out, stored or delivered
    sM = sum(M * d for M, d in zip(Ms, dts)) + S[0]
    sF = sum(F * d for F, d in zip(Fs, dts)) + S[1:4]
    sQ = sum(Q * d for Q, d in zip(Qs, dts)) + S[4]
    ib.judge('R1 mass', sM, m0, tol, 'SUM M dt + the removal vs rho_p V')
    ib.judge('R2 momentum', sF[0], m0 * up[0], tol, 'SUM Fx dt + the removal vs rho_p V u_p')
    ib.judge('R3 energy', sQ, E0, tol, 'SUM Q dt + the removal vs rho_p V (c T_p + u_p^2/2 + h298 - c T_ref)')
    ib.info('twin: energy per kg handed (SUM Q dt)/(SUM M dt) %.9f J/kg (the removal included)' % (sQ / sM))

    # the removal: the step that first stores an amount
    js = next((j for j, s in enumerate(Ss) if s[0] != 0.0), None)
    if js is None:
        ib.ko('R4 no step stores a removal')
        print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail)); sys.exit(1)
    held = m0 - sum(M * d for M, d in zip(Ms[:js], dts[:js]))           # before step js + 1
    hand = Ms[js] * dts[js] + Ss[js][0]
    dj = (6.0 * max(held, 0.0) / (math.pi * rho_mat * n_p * V)) ** (1.0 / 3.0)
    burn = n_p * V * rho_mat * math.pi * KX / (2.0 * nb) * dj ** (3.0 - nb) * dts[js]
    left = held - burn
    ib.info('removal at ICE step %d of %d: the cloud held %.6e of its mass (d %.6e m), the step burns %.6e and leaves %.6e (%.6e kg per drop)'
            % (js + 1, n, held / m0, dj, burn / m0, left / m0, left / (n_p * V)))
    (ib.ok if abs(hand - held) <= max(tol, 1e-9) * m0 else ib.ko)(
        'R4 the step of the removal hands all the cloud held: %.9e vs %.9e of its mass' % (hand / m0, held / m0))
    by_size, by_step = left <= m_b * n_p * V * (1.0 + 1e-9), left <= burn * (1.0 + 1e-9)
    (ib.ok if by_size or by_step else ib.ko)('R4 the removal fires by %s' % (
        'the size rule: a drop is left with %.6e kg, at or below %.6e kg' % (left / (n_p * V), m_b) if by_size
        else 'the rule of the long step: the step leaves %.6e and burnt %.6e of the mass' % (left / m0, burn / m0) if by_step
        else 'neither rule (left %.6e kg per drop, bound %.6e kg; burnt %.6e of the mass)' % (left / (n_p * V), m_b, burn / m0)))
    if js > 0:
        dprev = (6.0 * max(held + Ms[js - 1] * dts[js - 1], 0.0) / (math.pi * rho_mat * n_p * V)) ** (1.0 / 3.0)
        burn_prev = Ms[js - 1] * dts[js - 1]
        early = held > m_b * n_p * V and held > burn_prev
        (ib.ok if early else ib.ko)('R4 the step before left the drops above both bounds (%.6e kg per drop; held %.6e, burnt %.6e of the mass; d before it %.6e m)'
                                    % (held / (n_p * V), held / m0, burn_prev / m0, dprev))
    after = Ms[js + 1:]
    (ib.ok if len(after) >= 1 and all(M == 0.0 for M in after) and all(abs(Q) < 1e-20 for Q in Qs[js + 1:]) else ib.ko)(
        'R5 after the removal the source is zero: %d steps, largest mass rate %.1e kg/s, largest energy rate %.1e W'
        % (len(after), max([abs(M) for M in after] or [0.0]), max([abs(Q) for Q in Qs[js + 1:]] or [0.0])))
    drift = max([float(np.max(np.abs(s - Ss[js]) / np.maximum(np.abs(Ss[js]), 1e-300))) for s in Ss[js + 1:]] or [0.0])
    # the store and the delivery are two sums over the 400 cells: 400 x 1.1e-16 = 4.4e-14 at worst (1.5e-14 measured)
    (ib.ok if drift <= 1e-13 else ib.ko)(
        'R5 the amount removed, stored or delivered, stays what the removal took over the %d later steps: %.1e (tol 1e-13)'
        % (len(Ss[js + 1:]), drift))
    (ib.ok if left >= NON_VAC * m0 else ib.ko)('R5 non-vacuity: the removal hands %.3e of the cloud\'s mass beside the burn of its step (at least %g)'
                                               % (left / m0, NON_VAC))
    field = ib.read_tec(os.path.join(run, 'OUTPUT', 'part-field.tec'))
    rmax = float(max(np.abs(field[k]).max() for k in field if k.startswith('rho_p')))
    (ib.ok if rmax <= 1e-20 else ib.ko)('R6 ICE\'s final field is empty: largest bulk density %.1e kg/m3' % rmax)
    txt = open(log, errors='replace').read()
    nfl = len(re.findall(r'negative density after the update|negative densities floored', txt))
    (ib.ok if nfl == 0 else ib.ko)('R6 no density was floored (%d lines)' % nfl)

    print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
    sys.exit(0 if ib.nfail == 0 else 1)


if __name__ == '__main__':
    main()
