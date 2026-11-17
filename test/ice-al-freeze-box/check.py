#!/usr/bin/env python3
"""ice-al-freeze-box: the smoke of burning aluminium freezes, and its heat of fusion reaches the gas.

usage: check.py <run-dir> <log> [ice-box-dir]
  run-dir      the run as run_case.sh lays it out: input.ini, INPUT/, MESH/, OUTPUT/
  log          its log, with the [MI2-report] lines of [MI2-Coupling] report = true and ICE's step lines
  ice-box-dir  test/ice-box, whose check.py gives the readers (default: ../ice-box beside this file)

ICE alone (cycle 1, steps of dt-max, the gas held): AL(L) (family 1) burns direct down to d-stop; per kg burnt the
record of INPUT/part-products.txt hands the gas its species and makes nu_s kg of the smoke AL2O3(L) (family 2), born
at the drop's temperature; the smoke cools to the gas's. The report prints after ICE step j what the step hands the
gas: M_j (mass), F_j (force), Q_j (energy, the gas's formation scale). Energies on that scale, per kg of a material
of table H (INPUT/part-properties.dat, absolute): h298 + H(T) - H(T_ref) + |u|^2/2 (gate p30's). With E0 the
aluminium's at the start (INPUT/part-ic.tec), E_Al its at the end, G = SUM Q_j dt_j the gas's share:
  B1 mass      SUM M dt + M_Al + M_s = M_Al0,   M_s = nu_s (M_Al0 - M_Al)
  B3 energy    G + E_Al + E_s(ICE) = E0 (of the largest term), E_s(ICE) the smoke's energy as ICE holds it: the liquid
               table's, less the
               latent part f L(T) of its solidification model (L(T) = h-fus - (H(T-melt) - H(T)) + cp-solid
               (T-melt - T), Lib_Solidification.f90 solid_de_tab), with the constants of its line of part-phase.txt
  B4 size      rho_s / n_s = pi/6 rho_S d_s^3 in every cell
  F0 the box crosses the melting point: the smoke is born above T-melt 2327 K (the drops' temperature) and ends
     below T-nuc
  F1 the smoke is frozen at the end: f = 1 and chi = 1 in every cell
  F2 the books with the database's solid: E_s(true) = SUM rho_s V (h298_s + H_a(T) - H(T_ref) + |u|^2/2), H_a the
     table of AL2O3(a) GPB writes (test/igloo-smoke-box/INPUT/partE-properties.dat, the same database and scale).
     R = G + E_Al + E_s(true) - E0, per kg of smoke and h-fus: the error of ICE's constant cp-solid, eps(T) =
     H_a(T) - H_a(T-melt) + cp-solid (T-melt - T) (F2a: R equals SUM rho_s V eps(T_s) to 1e-6 of M_s h-fus), and
     within TOL_TRUTH of zero (F2b): the gas received the heat of fusion. A smoke on the liquid's table alone (no
     solidification model) gives R = SUM rho_s V (H_a - H)(T_s): -0.87 h-fus at 1750 K.
Prints PASS/FAIL lines with the numbers compared and the tolerance; exits 1 on any FAIL."""
import importlib.util
import math
import os
import re
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
icebox = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'ice-box')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ib = load('icebox_check', os.path.join(icebox, 'check.py'))          # read_ini, read_tec, uniform, judge, ok, ko, info
p30 = load('p30_al_smoke', os.path.join(here, '..', 'mi2-checks', 'p30_al_smoke.py'))   # zones, at, record

TOL_BOOK = 1e-10     # B1, B3: sums of 1500 terms printed with 17 figures; B3 of its largest term (the gas's share and
#                      the smoke's energy, about 21 J each on the formation scale, cancel to E0, 2.5 J)
TOL_SIZE = 1e-12     # B4
TOL_EPS = 1e-6       # F2a, of M_s h-fus: the books' round-off and the tables' rows
TOL_TRUTH = 0.03     # F2b, of M_s h-fus: eps(1750 K) / h-fus = 2.5e-2 with cp-solid = cp of AL2O3(a) at T-melt
T_REF = 298.15


def fnum(s):
    return float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', s))


def tokens(run, name):
    for l in open(os.path.join(run, 'INPUT', 'part-phase.txt')):
        f = l.split()
        if f and f[0] == name:
            return dict(t.split('=', 1) for t in f[2:])
    return {}


def main():
    run, log = os.path.abspath(sys.argv[1]), sys.argv[2]
    inp = os.path.join(run, 'INPUT') + '/'
    rec = p30.record(inp + 'part-products.txt')
    H = p30.zones(inp + 'part-properties.dat', 3)
    RH = p30.zones(inp + 'part-properties.dat', 2)
    HA = p30.zones(os.path.join(here, '..', 'igloo-smoke-box', 'INPUT', 'partE-properties.dat'), 3)['AL2O3(a)']
    HC, HO = H['AL(L)'], H['AL2O3(L)']
    tk = tokens(run, 'AL2O3(L)')
    solid = tk.get('solidification') == 'on'
    Tm = float(tk.get('T-melt', 2327.0))
    hf = float(tk.get('h-fus', 0.0)) if solid else HO[Tm] - HA[Tm]
    cps = float(tk.get('cp-solid', 0.0))
    Tnuc = float(tk.get('T-nuc', 0.8 * Tm)) or 0.8 * Tm
    row = [l.split() for l in open(inp + 'part-models.txt') if l.startswith('AL(L) combustion ')][0]
    d_s = float([t for t in row if t.startswith('d-smoke=')][0].split('=')[1])
    ib.info('the smoke AL2O3(L): solidification %s, T-melt %g K, h-fus %.6e J/kg, cp-solid %g, T-nuc %g K; nu_s %.10f'
            % ('on' if solid else 'off', Tm, hf, cps, Tnuc, rec['nu_s']))

    pic = ib.read_tec(inp + 'part-ic.tec')
    mesh = ib.read_tec(os.path.join(run, 'MESH', 'mesh.tec'))
    V = float(np.prod([mesh[k].max() - mesh[k].min() for k in ('x', 'y', 'z')]))
    r0 = pic['rp_11']
    Vc = V / r0.size
    u0 = np.array([pic['up_11'], pic['vp_11'], pic['wp_11']])
    T0 = pic['Tp_11']
    M0 = float(r0.sum()) * Vc
    E0 = float(sum(r0[c] * (rec['h298'] + p30.at(HC, T0[c]) - p30.at(HC, T_REF) + 0.5 * float(np.dot(u0[:, c], u0[:, c])))
                   for c in range(r0.size))) * Vc

    fld = ib.read_tec(os.path.join(run, 'OUTPUT', 'part-field.tec'))
    r1, T1 = fld['rho_p1'], fld['T_p1']
    u1 = np.array([fld['u_p1'], fld['v_p1'], fld['w_p1']])
    M1 = float(r1.sum()) * Vc
    E1 = float(sum(r1[c] * (rec['h298'] + p30.at(HC, T1[c]) - p30.at(HC, T_REF) + 0.5 * float(np.dot(u1[:, c], u1[:, c])))
                   for c in range(r1.size))) * Vc
    rs, Ts, ns = fld['rho_p2'], fld['T_p2'], fld['n_p2']
    us = np.array([fld['u_p2'], fld['v_p2'], fld['w_p2']])
    f = fld.get('f_p2')
    chi = fld.get('chi_p2')
    Ms = float(rs.sum()) * Vc

    def lat(T):
        return hf - (p30.at(HO, Tm) - p30.at(HO, T)) + cps * (Tm - T)

    kin = [0.5 * float(np.dot(us[:, c], us[:, c])) for c in range(rs.size)]
    E_liq = float(sum(rs[c] * (rec['h298_s'] + p30.at(HO, Ts[c]) - p30.at(HO, T_REF) + kin[c]) for c in range(rs.size))) * Vc
    E_ice = E_liq - (float(sum(rs[c] * f[c] * lat(Ts[c]) for c in range(rs.size))) * Vc if f is not None else 0.0)
    E_true = float(sum(rs[c] * (rec['h298_s'] + p30.at(HA, Ts[c]) - p30.at(HO, T_REF) + kin[c]) for c in range(rs.size))) * Vc
    eps = float(sum(rs[c] * (p30.at(HA, Ts[c]) - p30.at(HA, Tm) + cps * (Tm - Ts[c])) for c in range(rs.size))) * Vc

    dts, Ms_, Qs, Rs = [], [], [], np.zeros(5)
    for line in open(log, errors='replace'):
        m = re.match(r'ICE\s+\| Iter =\s+(\d+) \| Time = (\S+) \| Delta t = (\S+)', line)
        if m:
            if int(m.group(1)) == len(dts) + 1:
                dts.append(float(m.group(3)))
        elif line.startswith(' [MI2-report] ice-to-gas-mass'):
            Ms_.append(fnum(line.split()[2]))
        elif line.startswith(' [MI2-report] ice-to-gas '):
            Qs.append(fnum(line.split()[5]))
        elif line.startswith(' [MI2-report] ice-removed-delivered '):
            Rs = Rs + np.array([fnum(v) for v in line.split()[2:7]])
    n = min(len(dts), len(Ms_), len(Qs))
    (ib.ok if n > 0 and len(dts) == len(Ms_) == len(Qs) else ib.ko)(
        'P  %d ICE steps, %d mass lines, %d source lines in the report' % (len(dts), len(Ms_), len(Qs)))
    G_m = sum(a * b for a, b in zip(Ms_[:n], dts[:n])) + Rs[0]
    G = sum(a * b for a, b in zip(Qs[:n], dts[:n])) + Rs[4]
    ib.info('aluminium %.9e -> %.9e kg, smoke %.9e kg at %.3f..%.3f K; the gas\'s share %.9e kg, %.9e J'
            % (M0, M1, Ms, float(Ts.min()), float(Ts.max()), G_m, G))
    ib.judge('B1 mass', G_m + M1 + Ms, M0, TOL_BOOK, 'SUM M dt + aluminium + smoke vs the aluminium at the start')
    ib.judge('B1 the smoke', Ms, rec['nu_s'] * (M0 - M1), TOL_BOOK, 'SUM rho_s V vs nu_s x the aluminium burnt')
    big = max(abs(G), abs(E1), abs(E_ice), abs(E0))
    e3 = abs(G + E1 + E_ice - E0) / big
    (ib.ok if e3 <= TOL_BOOK else ib.ko)(
        'B3 energy: SUM Q dt + aluminium + smoke (as ICE holds it) %.16e vs E0 %.16e J, %.2e of the largest term, %.6e J '
        '(tol %.0e)' % (G + E1 + E_ice, E0, e3, big, TOL_BOOK))
    rho_S = p30.at(RH['AL2O3(L)'], float(Ts[0]))
    worst = float(np.max(np.abs(rs / ns / (math.pi / 6.0 * rho_S * d_s ** 3) - 1.0)))
    (ib.ok if worst <= TOL_SIZE else ib.ko)('B4 each smoke drop keeps the mass it is born with: worst %.1e (tol %.0e)'
                                            % (worst, TOL_SIZE))
    Tb = float(T0.max())
    (ib.ok if Tb > Tm and float(Ts.max()) < Tnuc else ib.ko)(
        'F0 the box crosses the melting point: born at the drops\' %.3f K (> T-melt %g K), ends at %.3f K (< T-nuc %g K)'
        % (Tb, Tm, float(Ts.max()), Tnuc))
    if f is None:
        ib.ko('F1 the smoke carries no frozen fraction: no column f_p2 in the field (no solidification model on it)')
    else:
        (ib.ok if float(f.min()) == 1.0 and float(chi.min()) == 1.0 else ib.ko)(
            'F1 the smoke is frozen: f %.6f..%.6f, chi %.6f..%.6f' % (float(f.min()), float(f.max()), float(chi.min()),
                                                                     float(chi.max())))
    R = G + E1 + E_true - E0
    scale = Ms * hf
    ib.info('R = G + E_Al + E_s(true) - E0 = %.9e J = %.6e of M_s h-fus; the constant cp-solid\'s error SUM rho_s V eps(T) '
            '= %.9e J = %.6e' % (R, R / scale, eps, eps / scale))
    if solid:
        (ib.ok if abs(R - eps) / scale <= TOL_EPS else ib.ko)(
            'F2a R is the constant cp-solid\'s error: R %.9e J vs SUM rho_s V eps %.9e J, %.1e of M_s h-fus (tol %.0e)'
            % (R, eps, abs(R - eps) / scale, TOL_EPS))
    else:
        ib.ko('F2a no solidification model on the smoke: R is not the constant cp-solid\'s error')
    (ib.ok if abs(R) / scale <= TOL_TRUTH else ib.ko)(
        'F2b the gas received the heat of fusion: |R| / (M_s h-fus) = %.4e (tol %.2f, the solid cp\'s error at the end)'
        % (abs(R) / scale, TOL_TRUTH))
    print('== check.py: %d passed, %d failed' % (ib.npass, ib.nfail))
    sys.exit(0 if ib.nfail == 0 else 1)


if __name__ == '__main__':
    main()
