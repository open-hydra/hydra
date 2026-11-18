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
  leg C  leg A on a liquid heat capacity that varies, c_l = a + b T (a = 50, b = 0.5 J/(kg K), the plan's table of
         IGLOO's solid-box-varcp): the Cp column and its integral a T + b T^2/2 as Enthalpy in both solvers' tables,
         ICE without its specific-heat key. Closed form (x = u t): the liquid x(T) = (u m/(Nu k pi d)) [b (T_0 - T)
         + (a + b T_g) ln((T_0 - T_g)/(T - T_g))], x_n = 0.033493 m, f0 = (h_l(T_m) - h_l(T_n))/h-fus = 0.477209,
         x_s = 0.061056 m, the solid from x_s with c_s, T_exit 902.5113 K, total heat 2014610.9 J/kg.

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
Leg C only (the energies on the table's chords, as both solvers interpolate it between its integer kelvins):
  (v)   the liquid upstream of nucleation, the cells up to 3 (o1) or 4 (o2) before nucleation's cell: o1 equals its
        backward-Euler march of the table (ICE's steady first-order cell i is that march at the cell's outflow face)
        to 1e-9 K, and differs from IGLOO's trajectory there by no more than the march's own error against the
        closed form plus IGLOO's row budget; o2 against IGLOO at the cell centres within IGLOO's budget plus twice
        the cell average's dx^2 |T''|/24. IGLOO's budget is that of its solid-box-varcp check: 4 x (5e-7 (2 + 2
        |T - T_g|/L(T)) + 1e-8), L(T) = u m c_l(T)/(Nu k pi d) (the F12.6 prints of x and T; 4 for the integrator's
        one-signed drift over the per-kelvin cp seams)
  and IGLOO against the closed form: every liquid row upstream of nucleation within that budget, the first plateau
  and first solid rows in the closed form's cells, the exit temperature within the budget at L_s.
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
        'B': dict(T_g=3000.0, T_0=1500.0, h_fus=4.0e5),
        'C': dict(T_g=600.0, T_0=2400.0, h_fus=1.07e6, cp=(50.0, 0.5))}
PLAN_TOL = {'A': 4561.8, 'B': 5726.4}          # the plan's total-heat bounds at nx = 240
TOL_T, TOL_Q = 1.0e-9, 1.0e-6
HALF_ULP, INT_FLOOR, SEAM_IG, TOL_MARCH = 0.5e-6, 1.0e-8, 4.0, 1.0e-9     # leg C: IGLOO's row budget, ICE's march


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
    if leg == 'C':
        xn = x_liquid(T_N, leg)
        f0 = (h_cont(T_M, leg) - h_cont(T_N, leg)) / h
        return xn, xn + U_G * (1.0 - f0) * M_P * h / (A_P * (T_M - Tg))
    if leg == 'A':
        xn = L_l * math.log((T0 - Tg) / (T_N - Tg))
        f0 = C_L * (T_M - T_N) / h
        return xn, xn + U_G * (1.0 - f0) * M_P * h / (A_P * (T_M - Tg))
    xm = L_s * math.log((Tg - T0) / (Tg - T_M))
    return xm, xm + U_G * M_P * h / (A_P * (Tg - T_M))


def energy(T, f, h):
    """ICE's specific energy and IGLOO's hSolid without its datum: c_l T - f (h - (c_l - c_s)(T_m - T))."""
    return C_L * T - f * (h - (C_L - C_S) * (T_M - T))


# leg C: the liquid's heat capacity c_l(T) = a + b T, the table's Enthalpy its exact integral a T + b T^2/2 on the
# integer kelvins, both solvers linear between the nodes (the chord model; it differs from the continuous form by
# b/8 K^2 at most, which moves the liquid's temperature by 5.1e-5 K on this stream)

def h_cont(T, leg):
    a, b = LEGS[leg]['cp']
    return a * T + 0.5 * b * T * T


def h_chord(T, leg):
    i = min(max(int(T), 1), 4999)
    return h_cont(i, leg) + (h_cont(i + 1, leg) - h_cont(i, leg)) * (T - i)


def x_liquid(T, leg):
    """Distance at which the liquid reaches T: x = (u m/(Nu k pi d)) [b (T_0 - T) + (a + b T_g) ln((T_0 - T_g)/(T - T_g))]."""
    a, b = LEGS[leg]['cp']
    Tg, T0 = LEGS[leg]['T_g'], LEGS[leg]['T_0']
    return U_G * M_P / A_P * (b * (T0 - T) + (a + b * Tg) * math.log((T0 - Tg) / (T - Tg)))


def T_liquid(x, leg):
    """The inverse of x_liquid, by bisection (x_liquid falls monotonically from T_0 to T_g)."""
    lo, hi = LEGS[leg]['T_g'] + 1e-9, LEGS[leg]['T_0']
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if x_liquid(mid, leg) > x:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def energy_tab(T, f, leg):
    """ICE's energy on the table without its datum: e_l(T) - f L(T), L(T) = h_fus - (e_l(T_m) - e_l(T)) + c_s (T_m - T);
    at f = 1 it is IGLOO's solid enthalpy h_l(T_m) - h_fus - c_s (T_m - T)."""
    h = LEGS[leg]['h_fus']
    return h_chord(T, leg) - f * (h - (h_chord(T_M, leg) - h_chord(T, leg)) + C_S * (T_M - T))


def be_liquid(leg, n):
    """ICE's first-order steady march of the liquid, cell by cell on the table's chords (backward Euler in x):
    e_l(T_i) - e_l(T_i-1) = Nu k pi d (T_g - T_i) dx/(u m), T_0 the inlet's. Returns T_1..T_n."""
    Tg = LEGS[leg]['T_g']
    k = A_P * DX / (U_G * M_P)
    out, Tp = [], LEGS[leg]['T_0']
    for _ in range(n):
        lo, hi = Tg, Tp
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if h_chord(mid, leg) - h_chord(Tp, leg) + k * (mid - Tg) > 0.0:
                hi = mid
            else:
                lo = mid
        Tp = 0.5 * (lo + hi)
        out.append(Tp)
    return out


def interp4(xs, ys, x):
    """Cubic Lagrange interpolation on the four rows around x (rows sorted by x)."""
    import bisect
    k = min(max(bisect.bisect_left(xs, x) - 2, 0), len(xs) - 4)
    tot = 0.0
    for j in range(k, k + 4):
        w = 1.0
        for m in range(k, k + 4):
            if m != j:
                w *= (x - xs[m]) / (xs[j] - xs[m])
        tot += w * ys[j]
    return tot


# ---------------------------------------------------------------------------------------------------
#  writers
# ---------------------------------------------------------------------------------------------------

def linspace(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def write_properties(path, leg='A'):
    if 'cp' in LEGS[leg]:
        a, b = LEGS[leg]['cp']
        rows = ['%.1f %.6f %.6f %.6f' % (T, a + b * T, RHO_P, h_cont(T, leg)) for T in range(1, 5001)]
    else:
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
    # field 1 is krho/(1 - krho_tot) over the dispersed phases of the cell (one here), field 10 G the cell's fixed
    # mass flux (none): IGLOO injects field 1 times (mdotGas + G area)
    inlet = '   %.16E   %.5E   normal,   normal,   %.5E   %.5E   %.5E   Dirac   %.5E   %.5E\n' \
            % (KRHO / (1.0 - KRHO), 1.0, p['T_0'] / p['T_g'], 0.5 * D_P, 0.0, 0.0, 0.0)
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
    write_properties(os.path.join(d, 'INPUT/properties.dat'), leg)
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
    if 'cp' in p:
        write_properties(os.path.join(d, 'INPUT/part-properties.dat'), leg)
    muscl = order == 'o2'
    with open(os.path.join(d, 'input.ini'), 'w') as f:
        f.write('\n'.join([
            '[ICE-Parameters]', 'iter-threshold = %d' % (20 * NX), 'time-threshold = 1e+30', 'res-threshold = 0.0', '',
            '[ICE-Numerics]', 'time-scheme = RK2', 'cfl = 0.8', 'time-accurate = false',
            'space-reconstruction = %s' % ('MUSCL' if muscl else 'first-order'),
            'flux-limiter = %s' % ('vanleer' if muscl else 'none'), '',
            '[ICE-Family1]', 'closure = MK', '',
            '[ICE-Physics]', 'drag = NoDrag', 'heat-transfer = Stokes', 'density = %g' % RHO_P,
            ] + ([] if 'cp' in p else ['specific-heat = %g' % C_L]) + ['emissivity = 0.0', '',
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
    # The zone's data end at the next record: a further zone or the auxiliary data a solver appends.
    body = re.split(r'^\s*(?:ZONE|DATASETAUXDATA)\b', body, maxsplit=1, flags=re.IGNORECASE | re.MULTILINE)[0]
    I, J, K = (int(re.search(r'\b%s\s*=\s*(\d+)' % a, header, re.IGNORECASE).group(1)) for a in 'IJK')
    vals = [float(t) for t in body.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(1, K - 1)
    names = [t.strip('"') for t in re.findall(r'"[^"]*"|[^\s,]+', re.search(r'VARIABLES\s*=\s*(.*)', head).group(1))][3:]
    return names, [vals[3 * nn + v * nc: 3 * nn + (v + 1) * nc] for v in range((len(vals) - 3 * nn) // nc)]


def load_rows(path, *cols):
    """Rows of an IGLOO particle file by parcel: {ID: [[<cols...>, ID], ...]}. The columns are found by name in the
    variables= line, so their position and any further column (t, bc, d_p, f) do not matter."""
    parts, idx = {}, None
    for line in open(path):
        if line.lstrip().lower().startswith('variables'):
            names = re.findall(r'"([^"]*)"', line)
            idx = [names.index(c) for c in cols + ('ID',)]
            continue
        c = line.split()
        if idx is None or len(c) <= max(idx):
            continue
        try:
            r = [float(c[i]) for i in idx[:-1]] + [int(c[idx[-1]])]
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
    traj = load_rows(os.path.join(d, 'OUTPUT/trajectories-A.dat'), 'X', 'T')
    exits = load_rows(os.path.join(d, 'OUTPUT/outloc-A.dat'), 'T', 'mdot')
    mdot = sum(r[0][1] for r in exits.values())
    Tp = [sum(T[i + NX * c] for c in range(NY * NZ)) / (NY * NZ) for i in range(NX)]
    q = [sum(E[i + NX * c] for c in range(NY * NZ)) / mdot if mdot > 0 else float('nan') for i in range(NX)]
    first_plat, first_liq = set(), set()
    for rows in traj.values():
        rows = sorted(rows, key=lambda r: r[0])
        plat = [r for r in rows if abs(r[1] - T_M) <= 0.5e-6]
        if plat:
            first_plat.add(int(round(plat[0][0] / DX)))
            after = [r for r in rows if r[0] > plat[-1][0]]
            if after:
                first_liq.add(int(round(after[0][0] / DX)))
    log = open(os.path.join(d, 'run_out.txt'), errors='replace').read()
    rows = sorted(traj[min(traj)], key=lambda r: r[0]) if traj else []
    return dict(Tp=Tp, q=q, first_plat=first_plat, first_after=first_liq, nparc=len(traj), nexit=len(exits),
                Texit=[r[0][0] for r in exits.values()], log=log, rows=rows,
                allrows=[r for v in traj.values() for r in v])


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
          % (leg, Tg, T0, h, 'x_m' if leg == 'B' else 'x_n', x1, 'x_l' if leg == 'B' else 'x_s', x2, DX))
    home1, home2 = int(x1 // DX) + 1, int(x2 // DX) + 1
    plateau = [i for i in range(NX) if i * DX > x1 + DX and (i + 1) * DX < x2 - DX]
    q_ref = A_P * abs(T_M - Tg) * DX / (U_G * M_P)
    try:
        ig = iglo_results(os.path.join(root, 'iglo'))
    except (IOError, OSError, ValueError, IndexError, ZeroDivisionError) as exc:
        rep.check(False, 'leg %s IGLOO: no readable output (%s)' % (leg, exc))
        ig = None
    if ig is not None:
        bad = [k for k in ('stuck in cell', 'no net progress', 'Inner loop') if k in ig['log']]
        rep.check(ig['nparc'] == NY * NZ and ig['nexit'] == NY * NZ and not bad,
                  'leg %s IGLOO: %d parcels, %d exits, no stuck, no-progress or inner-loop line%s'
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
        if leg in ('A', 'C'):
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
            Tx = sum(ig['Texit']) / len(ig['Texit'])
            if leg == 'C':
                e_in, e_out = energy_tab(T0, 0.0, leg), energy_tab(T[-1], f[-1], leg)
                ig_heat = e_in - energy_tab(Tx, 1.0, leg)
            else:
                e_in = energy(T0, ph_in, h)
                e_out = energy(T[-1], f[-1], h)
                ig_heat = e_in - (energy(Tx, 1.0, h) if leg == 'A' else C_L * Tx)
            c_exit, L_exit = (C_S, L_s) if leg in ('A', 'C') else (C_L, L_l)
            tol = c_exit * abs((Texit_o1 if Texit_o1 is not None else T[-1]) - Tg) / L_exit * 2.0 * DX
            dh = abs((e_in - e_out) - ig_heat)
            ref = ('the plan: %.1f' % PLAN_TOL[leg]) if leg in PLAN_TOL else \
                ('closed form %.1f' % (h_cont(T0, leg) - h_cont(T_M, leg) + h + C_S * (T_M - T_exit_cf(leg))))
            rep.check(dh <= tol, '%s (iv): total heat per unit mass flow %.1f (ICE) against %.1f J/kg (IGLOO), '
                      'off by %.1f within %.1f (%s)' % (tag, e_in - e_out, ig_heat, dh, tol, ref))
        if leg == 'C':
            check_liquid(rep, tag, order, T, ig, leg, home1)
    if leg == 'C' and ig is not None:
        check_closed_form(rep, ig, leg, home1, home2)
    return rep


def T_exit_cf(leg):
    """Leg C's closed-form exit temperature: the solid from x_s, T - T_g = (T_m - T_g) exp(-(x - x_s)/L_s)."""
    _, x2 = events(leg)
    return LEGS[leg]['T_g'] + (T_M - LEGS[leg]['T_g']) * math.exp(-(LX - x2) / lengths()[1])


def L_liquid(T, leg):
    """Leg C's local relaxation length u m c_l(T)/(Nu k pi d)."""
    a, b = LEGS[leg]['cp']
    return U_G * M_P * (a + b * T) / A_P


def budget_iglo(T, L):
    """IGLOO's row budget on a varying-cp stream (its solid-box-varcp check.py): SEAM x (the F12.6 prints of x and T,
    5e-7 (2 + 2 |T - T_g|/L), plus its integrator floor 1e-8); SEAM = 4 covers the one-signed drift of the integrator
    stepping over the per-kelvin cp seams (IGLOO measured 1.12-1.16 of the bare budget)."""
    return SEAM_IG * (HALF_ULP * (2.0 + 2.0 * abs(T - LEGS['C']['T_g']) / L) + INT_FLOOR)


def check_closed_form(rep, ig, leg, home1, home2):
    """Leg C, IGLOO against the closed form of the table: the liquid rows upstream of nucleation, the event cells, the
    exit temperature."""
    rows = [r for r in ig['allrows'] if r[0] <= (home1 - 3) * DX + 1e-9]
    res = [(abs(r[1] - T_liquid(r[0], leg)), budget_iglo(r[1], L_liquid(r[1], leg))) for r in rows]
    w = max(res, key=lambda v: v[0] / v[1])
    rep.check(len(rows) >= 25 * (home1 - 3) and all(v[0] <= v[1] for v in res),
              'leg %s IGLOO: liquid temperature on its %d rows upstream of nucleation against the closed form of the '
              'table, worst %.2e K = %.2f of its budget' % (leg, len(rows), w[0], w[0] / w[1]))
    ok = ig['first_plat'] == {home1} and ig['first_after'] == {home2}
    rep.check(ok, 'leg %s IGLOO: first plateau and first solid rows in cells %s, %s, the closed form\'s %d, %d'
              % (leg, sorted(ig['first_plat']), sorted(ig['first_after']), home1, home2))
    Tx, Tcf = sum(ig['Texit']) / len(ig['Texit']), T_exit_cf(leg)
    b = budget_iglo(Tcf, lengths()[1])
    rep.check(abs(Tx - Tcf) <= b, 'leg %s IGLOO: exit temperature %.4f K against the closed form %.4f K, off by %.1e '
              '= %.2f of its budget' % (leg, Tx, Tcf, abs(Tx - Tcf), abs(Tx - Tcf) / b))


def check_liquid(rep, tag, order, T, ig, leg, home1):
    """Leg C, the liquid upstream of nucleation, ICE against IGLOO's trajectory (cubic through its rows). First order:
    ICE's steady march is backward Euler in x on the table's chords (be_liquid), so cell i sits at its outflow face and
    |ICE - IGLOO| <= |march - closed form| + IGLOO's budget; the march itself to 1e-6 K. MUSCL: at the cell centres, within
    IGLOO's budget plus twice the cell average's own term dx^2 |T''|/24. Excluded: the cells within 2 (first order) or 3
    (MUSCL) of nucleation's cell, where the reconstruction sees the recalescence."""
    xs = [r[0] for r in ig['rows']]
    ts = [r[1] for r in ig['rows']]
    if order == 'o1':
        n = home1 - 3
        be = be_liquid(leg, n)
        res = []
        for i in range(n):
            xf = (i + 1) * DX
            ex = abs(T[i] - interp4(xs, ts, xf)) - abs(be[i] - T_liquid(xf, leg))
            res.append((ex, budget_iglo(T[i], L_liquid(T[i], leg)), abs(T[i] - be[i])))
        w = max(range(n), key=lambda i: res[i][0] / res[i][1])
        wm = max(range(n), key=lambda i: res[i][2])
        rep.check(all(v[0] <= v[1] for v in res) and res[wm][2] <= TOL_MARCH,
                  '%s (v): liquid temperature in cells 1-%d: ICE against its backward-Euler march of the table %.1e K '
                  '(cell %d, bound %.0e); against IGLOO at the outflow faces, beyond the march\'s own error, %.2e K '
                  '(cell %d) = %.2f of IGLOO\'s budget' % (tag, n, res[wm][2], wm + 1, TOL_MARCH, res[w][0], w + 1,
                                                          res[w][0] / res[w][1]))
    else:
        n = home1 - 4
        res = []
        for i in range(n):
            xc = (i + 0.5) * DX
            Tc = T_liquid(xc, leg)
            curv = (T_liquid(xc + DX, leg) - 2.0 * Tc + T_liquid(xc - DX, leg)) / DX ** 2 if i > 0 else \
                (T_liquid(xc + 2 * DX, leg) - 2.0 * T_liquid(xc + DX, leg) + Tc) / DX ** 2
            b = budget_iglo(T[i], L_liquid(T[i], leg)) + 2.0 * DX ** 2 * abs(curv) / 24.0
            res.append((abs(T[i] - interp4(xs, ts, xc)), b))
        w = max(range(n), key=lambda i: res[i][0] / res[i][1])
        rep.check(all(v[0] <= v[1] for v in res), '%s (v): liquid temperature in cells 1-%d against IGLOO at their '
                  'centres, worst %.2e K in cell %d = %.3f of its bound' % (tag, n, res[w][0], w + 1, res[w][0] / res[w][1]))


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
