#!/usr/bin/env python3
"""Solidification coherence between ICE and IGLOO on one stream (driver of p8_solid_coherence_check.sh).

  p8_solid_coherence.py make  <dir> <leg>    write <dir>/iglo (IGLOO) and <dir>/ice-o1, <dir>/ice-o2 (ICE)
  p8_solid_coherence.py check <dir> <leg>    read the three runs, print PASS/FAIL lines, exit 1 on a FAIL

One generator writes IGLOO's solid-box gas and material for both solvers on one x-spacing, nx = 240 over
0.15 m: IGLOO's box 240 x 5 x 5 over 0.15 x 0.05 x 0.05 m, ICE's strip 240 x 1 on the same x-nodes. The
particles move with the gas (IGLOO injects at the gas velocity, ICE runs NoDrag), Nu = 2 in both (IGLOO
Ranz-Marshall at Re = 0, ICE Stokes), so the heat per unit mass flow of a cell at T is Q(T) dx/(u m).

  leg A  a liquid at 2400 K into a 600 K gas, h-fus = 1.07e6 (IGLOO's solid-box): cooling, nucleation,
         the freezing plateau, the solid.
  leg B  a solid at 1500 K into a 3000 K gas, h-fus = 4e5: heating to T-melt, the melting plateau, the liquid.

ICE runs steady twice: o1 first order, o2 MUSCL with van Leer. Checks, per ICE run:
  (i)   leg A: ICE's first cell with f > 0 and first cell with f = 1 against the cells of IGLOO's first
        plateau row and first row after the plateau, +-1;
        leg B: ICE's first cell with f < 1 and first cell with chi = 0 against the cells of IGLOO's
        first plateau row and first liquid row, +-1
  (ii)  every cell with both faces inside the plateau of the closed form, one cell in from each end:
        T = T-melt in ICE and in IGLOO's Eulerian T_p (1e-9 relative)
  (iii) the heat to the gas per unit mass flow in those cells, ICE n Q(T) dx/(rho u) against IGLOO's
        source.tec E over the 25 cells of the slab per unit mass flow, to 1e-6
  (iv)  the total heat per unit mass flow, e_in - e_out (ICE) against h(T_0) - h(T_exit) (IGLOO), within
        the first-order front error c (T_exit - T_g)/L 2 dx (c, L those of the exit phase, T_exit ICE's o1)
"""
import math
import os
import re
import sys

LX, LY, LZ = 0.15, 0.05, 0.05
NX, NY, NZ = 240, 5, 5
DX = LX / NX
RHO_G, U_G, MU_G, K_G, GAM, R_G, P_G = 1.2, 10.0, 1.8e-5, 0.026, 1.4, 287.0, 101325.0
RHO_P, C_L, C_S, T_M = 2950.0, 1250.0, 600.0, 2327.0
T_N = 0.8 * T_M
D_P = 30.0e-6
NU = 2.0
M_P = RHO_P * math.pi / 6.0 * D_P ** 3
A_P = math.pi * D_P * K_G * NU                 # W/K per particle
KRHO, GP = 0.34, 29.5                          # IGLOO inlet loading; ICE inlet mass flux kg/(s m^2)
LEGS = {'A': dict(T_g=600.0, T_0=2400.0, h_fus=1.07e6),
        'B': dict(T_g=3000.0, T_0=1500.0, h_fus=4.0e5)}
PLAN_TOL = {'A': 4561.8, 'B': 5726.4}          # the plan's total-heat bounds at nx = 240
TOL_T, TOL_Q = 1.0e-9, 1.0e-6


# ---------------------------------------------------------------------------------------------------
#  closed forms (x = u t)
# ---------------------------------------------------------------------------------------------------

def lengths():
    return U_G * M_P * C_L / A_P, U_G * M_P * C_S / A_P           # L_l, L_s


def events(leg):
    """Leg A: (x_n, x_s), nucleation and end of freezing; leg B: (x_m, x_l), melt start and melt end."""
    p = LEGS[leg]
    L_l, L_s = lengths()
    Tg, T0, h = p['T_g'], p['T_0'], p['h_fus']
    if leg == 'A':
        xn = L_l * math.log((T0 - Tg) / (T_N - Tg))
        f0 = C_L * (T_M - T_N) / h
        return xn, xn + U_G * (1.0 - f0) * M_P * h / (A_P * (T_M - Tg))
    xm = L_s * math.log((Tg - T0) / (Tg - T_M))
    return xm, xm + U_G * M_P * h / (A_P * (Tg - T_M))


def energy(T, f, h):
    """ICE's specific energy and IGLOO's hSolid without its datum: c_l T - f (h - (c_l - c_s)(T_m - T))."""
    return C_L * T - f * (h - (C_L - C_S) * (T_M - T))


# ---------------------------------------------------------------------------------------------------
#  writers
# ---------------------------------------------------------------------------------------------------

def linspace(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def write_properties(path):
    rows = ['%.1f %.6f %.6f %.6f' % (T, C_L, RHO_P, C_L * T) for T in range(1, 5001)]
    with open(path, 'w') as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\n'
                'VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"\nZONE T="A"\nI=5000, F=POINT\n')
        f.write('\n'.join(rows) + '\n')


def write_iglo(d, leg):
    p = LEGS[leg]
    os.makedirs(os.path.join(d, 'INPUT'))
    os.makedirs(os.path.join(d, 'OUTPUT'))
    gas = [('Roi( 1 )', RHO_G), ('U', U_G), ('V', 0.0), ('W', 0.0), ('P', P_G), ('T', p['T_g']),
           ('MIL', MU_G), ('KL', K_G), ('GAM', GAM), ('R', R_G), ('MIT', MU_G),
           ('r_p', 0.0), ('u_p', 0.0), ('v_p', 0.0), ('T_p', 0.0), ('n_p', 0.0), ('R_p', 0.0), ('odot', 0.0)]
    xs, ys, zs = linspace(0.0, LX, NX + 1), linspace(0.0, LY, NY + 1), linspace(0.0, LZ, NZ + 1)
    ncell = NX * NY * NZ
    with open(os.path.join(d, 'INPUT/solfile.tec'), 'w') as f:
        f.write(' TITLE="p8 coherence: uniform-gas box"\n')
        f.write(' VARIABLES = "X", "Y", "Z",' + ','.join('"%s"' % n for n, _ in gas) + '\n')
        f.write(' Zone T="block: 1",I=%d,J=%d,K=%d\n' % (NX + 1, NY + 1, NZ + 1))
        f.write(' ,DATAPACKING=BLOCK,  VARLOCATION=([4-21]=CELLCENTERED)\n')
        for comp in range(3):
            f.writelines('  %.15E\n' % (x, y, z)[comp] for z in zs for y in ys for x in xs)
        for _, v in gas:
            f.writelines(['  %.15E\n' % v] * ncell)
    ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (NX, NZ), 4: (NX, NZ), 5: (NX, NY), 6: (NX, NY)}
    inlet = '   %.5E   %.5E   normal,   normal,   %.5E   %.5E   %.5E   Dirac   %.5E\n' \
            % (KRHO, 1.0, p['T_0'] / p['T_g'], 0.5 * D_P, 0.0, 0.0)
    with open(os.path.join(d, 'INPUT/bc.txt'), 'w') as f:
        for face in range(1, 7):
            mend, nend = ranges[face]
            code = 401 if face == 1 else 100
            for n in range(1, nend + 1):
                for m in range(1, mend + 1):
                    f.write('   %6d%6d%6d%6d%6d%6d\n' % (1, face, m, n, 1, code))
                    if code == 401:
                        f.write(inlet)
    with open(os.path.join(d, 'INPUT/phase.txt'), 'w') as f:
        f.write('condensed-dispersed\nA 1 solidification=on h-fus=%g cp-solid=%g\n' % (p['h_fus'], C_S))
    write_properties(os.path.join(d, 'INPUT/properties.dat'))
    with open(os.path.join(d, 'input.ini'), 'w') as f:
        f.write("""; p8 coherence stream, IGLOO side: one parcel per inlet cell at the gas velocity (Re 0, Nu 2).

[IGLOO-General]
gas-file    = INPUT/solfile.tec
gas-order   = 1
print-dcell = 1
out-file    = ALL
mollify     = off

[IGLOO-ODE]
ode-solver   = H-sdirk4
relative-tol = 1e-11
absolute-tol = 1e-11

[IGLOO-Models]
drag = Stokes
heat = Ranz-Marshall

[IGLOO-BC]
ds = 10.0

[GPB-Phase1]
type = condensed-dispersed
rho  = %g
cp   = %g

[BCB-Block1]
face1 = in
face2 = out
face3 = wall
face4 = wall
face5 = wall
face6 = wall

[out]
type = outlet

[in]
type = inlet
krho = %g
dp   = %g

[wall]
type = wall
q    = 0.0
""" % (RHO_P, C_L, KRHO, D_P))


def tec_block(path, names, xn, yn, zn, cellvars, zone):
    nvar = 3 + len(names)
    lines = [' VARIABLES ="x" "y" "z" ' + ' '.join('"%s"' % v for v in names),
             ' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED)'
             % (zone, len(xn), len(yn), len(zn), nvar)]
    for comp in range(3):
        lines += [' %.15E' % (x, y, z)[comp] for z in zn for y in yn for x in xn]
    for v in cellvars:
        lines += [' %.15E' % c for c in v]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def write_ice(d, leg, order):
    p = LEGS[leg]
    os.makedirs(os.path.join(d, 'INPUT'))
    os.makedirs(os.path.join(d, 'OUTPUT'))
    xn, yn, zn = linspace(0.0, LX, NX + 1), [0.0, DX], [0.0, DX]
    rho_in = GP / U_G
    # leg A starts from the inlet state everywhere, leg B from an empty domain (ICE's O1 and O4)
    rho0, T0 = (rho_in, p['T_0']) if leg == 'A' else (1.0e-20, p['T_0'])
    n0 = rho0 / (RHO_P * 4.0 / 3.0 * math.pi * (0.5 * D_P) ** 3)
    tec_block(os.path.join(d, 'INPUT/part-ic.tec'), ['rp1', 'up1', 'vp1', 'wp1', 'Tp1', 'np1'], xn, yn, zn,
              [[v] * NX for v in (rho0, U_G, 0.0, 0.0, T0, n0)], 'B1-CD')
    tec_block(os.path.join(d, 'INPUT/gas.tec'), ['rho', 'u', 'v', 'w', 'T', 'R', 'gam', 'k', 'mu', 'dt'],
              xn, yn, zn, [[v] * NX for v in (RHO_G, U_G, 0.0, 0.0, p['T_g'], R_G, GAM, K_G, MU_G, 0.0)], 'B1-GAS')
    rec = lambda i, f, code: '%8d%8d%8d%8d%8d%8d' % (1, i, 1, 1, f, code)       # noqa: E731
    rows = [rec(1, 1, 402), '   %.10e %.10e normal normal %.10e %.10e 0.0 Dirac 0.0' % (GP, U_G, p['T_0'], 0.5 * D_P),
            rec(NX, 2, 400)]
    rows += [rec(i, f, 400) for f in (3, 4) for i in range(1, NX + 1)]
    rows += [rec(i, f, 0) for f in (5, 6) for i in range(1, NX + 1)]
    with open(os.path.join(d, 'INPUT/part-bc.txt'), 'w') as f:
        f.write('\n'.join(rows) + '\n')
    with open(os.path.join(d, 'INPUT/part-phase.txt'), 'w') as f:
        f.write('condensed-dispersed phase\nA 1 solidification=on h-fus=%g cp-solid=%g\n' % (p['h_fus'], C_S))
    muscl = order == 'o2'
    with open(os.path.join(d, 'input.ini'), 'w') as f:
        f.write('\n'.join([
            '[ICE-Parameters]', 'iter-threshold = %d' % (20 * NX), 'time-threshold = 1e+30', 'res-threshold = 0.0', '',
            '[ICE-Numerics]', 'time-scheme = RK2', 'cfl = 0.8', 'time-accurate = false',
            'space-reconstruction = %s' % ('MUSCL' if muscl else 'first-order'),
            'flux-limiter = %s' % ('vanleer' if muscl else 'none'), '',
            '[ICE-Family1]', 'closure = MK', '',
            '[ICE-Physics]', 'drag = NoDrag', 'heat-transfer = Stokes', 'density = %g' % RHO_P,
            'specific-heat = %g' % C_L, 'emissivity = 0.0', '',
            '[ICE-IO]', 'ic-format = tecplot ascii', 'sol-format = tecplot ascii', 'shell-diter = 1000000000',
            'sol-diter = 1000000000', 'res-diter = 1000000000', '']))


# ---------------------------------------------------------------------------------------------------
#  readers
# ---------------------------------------------------------------------------------------------------

def read_block(path):
    """First zone of a Tecplot BLOCK file: (names without x y z, [cell values per variable])."""
    text = open(path).read()
    head, rest = re.split(r'^\s*ZONE', text, maxsplit=1, flags=re.IGNORECASE | re.MULTILINE)
    header, body = rest.split('\n', 1)
    I, J, K = (int(re.search(r'\b%s\s*=\s*(\d+)' % a, header, re.IGNORECASE).group(1)) for a in 'IJK')
    vals = [float(t) for t in body.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(1, K - 1)
    names = re.findall(r'"([^"]*)"', head)[3:]
    return names, [vals[3 * nn + v * nc: 3 * nn + (v + 1) * nc] for v in range((len(vals) - 3 * nn) // nc)]


def load_rows(path, ncol):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) != ncol:
            continue
        try:
            r = [float(v) for v in c[:-1]] + [int(c[-1])]
        except ValueError:
            continue
        parts.setdefault(r[-1], []).append(r)
    return parts


def iglo_results(d):
    """Per x-slab Euler T_p and heat per unit mass flow, the cells of the first plateau and first liquid rows
    (each parcel's), the exit temperatures and the log."""
    _, (rho, u, v, w, T, n) = read_block(os.path.join(d, 'OUTPUT/euler1.tec'))
    _, src = read_block(os.path.join(d, 'OUTPUT/source.tec'))
    E = src[4]
    traj = load_rows(os.path.join(d, 'OUTPUT/trajectories-A.dat'), 10)
    exits = load_rows(os.path.join(d, 'OUTPUT/outloc-A.dat'), 9)
    mdot = sum(r[0][6] for r in exits.values())
    Tp = [sum(T[i + NX * c] for c in range(NY * NZ)) / (NY * NZ) for i in range(NX)]
    q = [sum(E[i + NX * c] for c in range(NY * NZ)) / mdot if mdot > 0 else float('nan') for i in range(NX)]
    first_plat, first_liq = set(), set()
    for rows in traj.values():
        rows = sorted(rows, key=lambda r: r[0])
        plat = [r for r in rows if abs(r[6] - T_M) <= 0.5e-6]
        if plat:
            first_plat.add(int(round(plat[0][0] / DX)))
            after = [r for r in rows if r[0] > plat[-1][0]]
            if after:
                first_liq.add(int(round(after[0][0] / DX)))
    log = open(os.path.join(d, 'run_out.txt'), errors='replace').read()
    return dict(Tp=Tp, q=q, first_plat=first_plat, first_after=first_liq, nparc=len(traj), nexit=len(exits),
                Texit=[r[0][3] for r in exits.values()], log=log)


def ice_results(d):
    names, var = read_block(os.path.join(d, 'OUTPUT/part-field.tec'))
    col = {n: var[i] for i, n in enumerate(names)}
    return col


# ---------------------------------------------------------------------------------------------------
#  check
# ---------------------------------------------------------------------------------------------------

class Rep(object):
    def __init__(self):
        self.npass = self.nfail = 0

    def check(self, ok, text):
        print(('PASS  ' if ok else 'FAIL  ') + text)
        if ok:
            self.npass += 1
        else:
            self.nfail += 1


def check(root, leg):
    p = LEGS[leg]
    Tg, T0, h = p['T_g'], p['T_0'], p['h_fus']
    L_l, L_s = lengths()
    x1, x2 = events(leg)
    rep = Rep()
    print('  leg %s: T_g %g K, T_0 %g K, h-fus %g; closed form %s = %.6f m, %s = %.6f m, dx = %.6f m'
          % (leg, Tg, T0, h, 'x_n' if leg == 'A' else 'x_m', x1, 'x_s' if leg == 'A' else 'x_l', x2, DX))
    home1, home2 = int(x1 // DX) + 1, int(x2 // DX) + 1
    plateau = [i for i in range(NX) if i * DX > x1 + DX and (i + 1) * DX < x2 - DX]
    q_ref = A_P * abs(T_M - Tg) * DX / (U_G * M_P)
    try:
        ig = iglo_results(os.path.join(root, 'iglo'))
    except (IOError, OSError, ValueError, IndexError, ZeroDivisionError) as exc:
        rep.check(False, 'leg %s IGLOO: no readable output (%s)' % (leg, exc))
        ig = None
    if ig is not None:
        bad = [k for k in ('stuck in cell', 'no net progress', 'Inner loop', 'above T-melt') if k in ig['log']]
        rep.check(ig['nparc'] == NY * NZ and ig['nexit'] == NY * NZ and not bad,
                  'leg %s IGLOO: %d parcels, %d exits, no stuck, no-progress, inner-loop or above-T-melt line%s'
                  % (leg, ig['nparc'], ig['nexit'], (' (found: %s)' % ', '.join(bad)) if bad else ''))
        print('       IGLOO: first plateau row in cell(s) %s, first row after the plateau in cell(s) %s, exit %.4f K'
              % (sorted(ig['first_plat']), sorted(ig['first_after']),
                 sum(ig['Texit']) / max(1, len(ig['Texit']))))
    Texit_o1 = None
    for order in ('o1', 'o2'):
        tag = 'leg %s ICE %s' % (leg, 'first order' if order == 'o1' else 'MUSCL')
        d = os.path.join(root, 'ice-' + order)
        try:
            col = ice_results(d)
            T, f, chi = col['T_p1'], col['f_p1'], col['chi_p1']
            rho, u, n = col['rho_p1'], col['u_p1'], col['n_p1']
        except (IOError, OSError, KeyError, ValueError, IndexError) as exc:
            why = ''
            logf = os.path.join(d, 'log')
            if os.path.exists(logf):
                hits = [s.strip() for s in open(logf, errors='replace') if re.search(r'ERROR|error stop|Wrong', s)]
                why = (': ' + hits[0]) if hits else ''
            rep.check(False, '%s: no solution with f and chi (%s)%s' % (tag, exc.__class__.__name__, why))
            continue
        if order == 'o1':
            Texit_o1 = T[-1]
        if leg == 'A':
            fr = next((i + 1 for i in range(NX) if f[i] > 0.0), None)
            fz = next((i + 1 for i in range(NX) if f[i] == 1.0), None)
            print('       %s: first cell with f > 0: %s, first with f = 1: %s (closed form %d, %d), exit %.4f K'
                  % (tag, fr, fz, home1, home2, T[-1]))
            if ig is not None:
                ok = (fr is not None and fz is not None and len(ig['first_plat']) == 1 and len(ig['first_after']) == 1
                      and abs(fr - min(ig['first_plat'])) <= 1 and abs(fz - min(ig['first_after'])) <= 1)
                rep.check(ok, '%s (i): freezing cells (%s, %s) against IGLOO (%s, %s), +-1'
                          % (tag, fr, fz, sorted(ig['first_plat']), sorted(ig['first_after'])))
        else:
            mush = next((i + 1 for i in range(NX) if f[i] < 1.0), None)
            liq = next((i + 1 for i in range(NX) if chi[i] == 0.0), None)
            print('       %s: first cell with f < 1: %s, first with chi = 0: %s (closed form %d, %d), exit %.4f K'
                  % (tag, mush, liq, home1, home2, T[-1]))
            if ig is not None:
                ok = (mush is not None and liq is not None and len(ig['first_plat']) == 1 and len(ig['first_after']) == 1
                      and abs(mush - min(ig['first_plat'])) <= 1 and abs(liq - min(ig['first_after'])) <= 1)
                rep.check(ok, '%s (i): melt cells (%s, %s) against IGLOO (%s, %s), +-1'
                          % (tag, mush, liq, sorted(ig['first_plat']), sorted(ig['first_after'])))
        wT = max(abs(T[i] - T_M) / T_M for i in plateau)
        wTi = max(abs(ig['Tp'][i] - T_M) / T_M for i in plateau) if ig is not None else float('inf')
        rep.check(wT <= TOL_T and wTi <= TOL_T, '%s (ii): T = T_m in the %d plateau cells, ICE %.1e, IGLOO %.1e'
                  % (tag, len(plateau), wT, wTi))
        qi = [n[i] * A_P * (T[i] - Tg) * DX / (rho[i] * u[i]) for i in range(NX)]
        if ig is not None:
            wq = max(abs(qi[i] - ig['q'][i]) / abs(ig['q'][i]) for i in plateau)
            wr = max(abs(abs(qi[i]) - q_ref) / q_ref for i in plateau)
            other = [i for i in range(NX) if i not in plateau]
            wo = max(abs(qi[i] - ig['q'][i]) for i in other)
            rep.check(wq <= TOL_Q, '%s (iii): heat per unit mass flow in the plateau cells, ICE against IGLOO %.1e '
                      '(ICE against Q(T_m) dx/(u m) = %.2f J/kg: %.1e)' % (tag, wq, q_ref, wr))
            print('       %s: elsewhere the per-cell heat differs by at most %.1f J/kg' % (tag, wo))
            ph_in = 1.0 if leg == 'B' else 0.0
            e_in = energy(T0, ph_in, h)
            e_out = energy(T[-1], f[-1], h)
            Tx = sum(ig['Texit']) / len(ig['Texit'])
            ig_heat = e_in - (energy(Tx, 1.0, h) if leg == 'A' else C_L * Tx)
            c_exit, L_exit = (C_S, L_s) if leg == 'A' else (C_L, L_l)
            tol = c_exit * abs((Texit_o1 if Texit_o1 is not None else T[-1]) - Tg) / L_exit * 2.0 * DX
            dh = abs((e_in - e_out) - ig_heat)
            rep.check(dh <= tol, '%s (iv): total heat per unit mass flow %.1f (ICE) against %.1f J/kg (IGLOO), '
                      'off by %.1f within %.1f (the plan: %.1f)' % (tag, e_in - e_out, ig_heat, dh, tol, PLAN_TOL[leg]))
    return rep


def main():
    if len(sys.argv) != 4 or sys.argv[1] not in ('make', 'check') or sys.argv[3] not in LEGS:
        sys.exit(__doc__)
    what, root, leg = sys.argv[1:]
    if what == 'make':
        write_iglo(os.path.join(root, 'iglo'), leg)
        for order in ('o1', 'o2'):
            write_ice(os.path.join(root, 'ice-' + order), leg, order)
        return 0
    rep = check(root, leg)
    print('  leg %s: %d passed, %d failed' % (leg, rep.npass, rep.nfail))
    return 1 if rep.nfail else 0


if __name__ == '__main__':
    sys.exit(main())
