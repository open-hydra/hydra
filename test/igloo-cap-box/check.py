#!/usr/bin/env python3
"""igloo-cap-box: an IGLOO material that burns with an oxide cap, its residue left in the box: the gas species, the
mass and the energy the coupler hands the gas, against closed forms.

  python3 check.py RUN_DIR LOG [EVAP_BOX_MOM_DIR]

RUN_DIR holds input.ini, INPUT/ and OUTPUT/ of the run; LOG is its log with the [MI2-report] lines; the third
argument is test/evap-box-mom, whose check.py gives the readers (default: ../evap-box-mom beside this file).

Streams of a metal A (mass flow mdot, injected at u_inj and T0) burn out; a share xi (xi-cap) of the oxide forms on
the drop as a cap of the material OX and is left as a residue (residue = on). INPUT/part-products.txt gives three
records: "combustion" of A (gas species per kg of A, row g, sum 1), "cap" of A (nu kg of OX per kg of A and the gas
species it draws, row c, sum 1 - nu) and "release" of OX (not used here: the residue stays). IGLOO's cap law forms
xi nu kg of OX per kg of A burnt, so per kg of A the cap takes xi kg and 1 - xi kg burns into the gas:
  M      = SUM mdot (1 - xi nu)                           net mass to the gas                       [kg/s]
  m_s    = SUM mdot [(1 - xi) g(s) + xi c(s)]             each gas species s                        [kg/s]
  Q      = SUM mdot [h_A(T0) + u_inj^2/2] - SUM xi nu mdot [h_OX(T_end) + u_end^2/2]                [W]
  d_res  = d0 (xi nu rho_A / rho_OX)^(1/3)                the residue of a drop                     [m]
h_A and h_OX are the tables of INPUT/part-properties.dat (on the scale of the gas: the materials give h0); T_end and
u_end are the end records of the residues (OUTPUT/part-outloc-OX*.dat) of the sweep whose sources the report
carries. IGLOO integrates the cap beside the metal: the residue's mass flow is xi nu mdot within its step control,
so M and the species are met within TOL_CAP and not to round-off; the split of the species between the three rows
is exact for the mass the report gives (S2).
One PASS/FAIL line per check with the numbers and the tolerance; exit 1 on any FAIL."""
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
mom = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'evap-box-mom')
sys.path.insert(0, mom)
import check as eb          # test/evap-box-mom/check.py: read_ini, values, read_report

TOL_CAP = 1e-6      # S1  M, the species and Q against the closed forms (IGLOO's step control on the cap: 2e-8 measured;
#                         the species are printed with seven figures)
TOL_SPLIT = 1e-6    # S2  the species against the three rows for the reported mass (seven printed figures)
TOL_EN = 1e-5       # S3  E1 - E0 = Q (t1 - t0) per block (test/igloo-burn-box's D2)
TOL_BOOK_M = 1e-9   # S4  the coupler's closed-domain mass balance
TOL_BOOK_E = 1e-5   # S4  the coupler's closed-domain energy bookkeeping
TOL_RES = 1e-5      # S5  the residue's diameter and mass flow (six printed figures)
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


def records(path):
    """{(material, process): (h298, {species: nu_mass}, product, nu_product, h298 of the product)}."""
    lines = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    out = {}
    for i, f in enumerate(lines):
        if f[0] == 'material':
            tail = (f[6], float(f[7]), float(f[8])) if f[2] == 'cap' else ('', 0.0, 0.0)
            out[(f[1], f[2])] = (float(f[5]), {g[0]: float(g[1]) for g in lines[i + 1:i + 1 + int(f[3])]}) + tail
    return out


def outloc(path):
    """{ID: (T, |u|, mdot, d)} of the end records of a part-outloc file."""
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) == 12:
            try:
                out[int(f[8])] = (float(f[3]), float(f[4]), float(f[6]), float(f[11]))
            except ValueError:
                continue
    return out


def main():
    run, log = os.path.abspath(sys.argv[1]), sys.argv[2]
    ini = eb.read_ini(os.path.join(run, 'input.ini'))
    g = ini['GPB-Phase2']
    mat, prod = g['material'].split()[0], g['product'].split()[0]
    n = len(ini['IGLOO-BC']['y'].split())
    mdot, up, T0, d0 = (eb.values(ini, 'IGLOO-BC', k, n) for k in ('mdot', 'up', 'temp0', 'diam'))
    xi = float(g['xi-cap'])
    rhoA, rhoP = (float(v) for v in g['rho'].split()[:2])
    order = [l.split()[0] for l in open(os.path.join(run, 'INPUT', 'gas-phase.txt')).read().split('\n')[1:] if l.strip()]
    rec = records(os.path.join(run, 'INPUT', 'part-products.txt'))
    if (mat, 'combustion') not in rec or (mat, 'cap') not in rec:
        ko('C  INPUT/part-products.txt has no combustion or no cap record of %s' % mat)
        print('== check.py: %d passed, %d failed' % (npass, nfail)); sys.exit(1)
    h298, rowg = rec[(mat, 'combustion')][:2]
    _, rowc, pname, nu, h298P = rec[(mat, 'cap')]
    (ok if g.get('residue', '').split()[0] == 'on' and pname == prod and abs(sum(rowg.values()) - 1.0) < 1e-9
        and abs(sum(rowc.values()) + nu - 1.0) < 1e-9 and all(s in order for s in list(rowg) + list(rowc)) else ko)(
        'C  %s burns with a cap of %s (xi %.3f, %.10f kg per kg), residue on; the two rows close the mass and name species of the gas'
        % (mat, pname, xi, nu))

    ptab, pnames = zones(os.path.join(run, 'INPUT', 'part-properties.dat'))
    hcol = [i for i, v in enumerate(pnames) if v.startswith('Enthalpy')][0]
    hA, hP = ptab[mat], ptab[prod]
    (ok if abs(at(hA, T_REF, hcol) - h298) <= 1e-9 * abs(h298) and abs(at(hP, T_REF, hcol) - h298P) <= 1e-6 + 1e-9 * abs(h298P) else ko)(
        'C  the tables of %s and %s are on the scale of the gas: h(T_ref) %.9e and %.3e J/kg, the records %.9e and %.3e'
        % (mat, prod, at(hA, T_REF, hcol), at(hP, T_REF, hcol), h298, h298P))

    steps = [s for s in eb.read_report(log) if 'F' in s]
    if not steps:
        ko('S1 the report has no igloo-to-gas line')
        print('== check.py: %d passed, %d failed' % (npass, nfail)); sys.exit(1)
    last = steps[-1]
    nsrc = len(steps)
    name = 'part-outloc-%s.dat' % prod if nsrc <= 1 else 'part-outloc-%s-sweep%d.dat' % (prod, nsrc - 1)
    end = outloc(os.path.join(run, 'OUTPUT', name))
    (ok if len(end) == n else ko)('S5 one residue per stream: %d end records of %s for %d streams' % (len(end), prod, n))
    if len(end) != n:
        print('== check.py: %d passed, %d failed' % (npass, nfail)); sys.exit(1)

    A = sum(mdot)
    M_an = sum(m * (1.0 - xi * nu) for m in mdot)
    Q_an = sum(m * (at(hA, t, hcol) + 0.5 * u * u) for m, u, t in zip(mdot, up, T0)) \
        - sum(xi * nu * m * (at(hP, end[i + 1][0], hcol) + 0.5 * end[i + 1][1] ** 2) for i, m in enumerate(mdot))
    rate = {s: A * ((1.0 - xi) * rowg.get(s, 0.0) + xi * rowc.get(s, 0.0)) for s in set(rowg) | set(rowc)}
    print('      closed forms: M %.15e kg/s, Q %.15e W; per kg of %s: %s'
          % (M_an, Q_an, mat, ', '.join('%s %.9f' % (s, rate[s] / A) for s in sorted(rate))))

    # ---- S1: the report against the closed forms
    judge('S1 mass', last['M'], M_an, TOL_CAP, 'reported M vs SUM mdot (1 - xi nu)')
    judge('S1 energy', last['Q'], Q_an, TOL_CAP, 'reported Q vs SUM mdot (h_A(T0) + u^2/2) - SUM xi nu mdot (h_OX(T_end) + u_end^2/2)')
    Fx, bound = sum(m * u for m, u in zip(mdot, up)), sum(xi * nu * m * end[i + 1][1] for i, m in enumerate(mdot))
    (ok if abs(last['F'][0] - Fx) <= bound * (1.0 + 1e-6) else ko)(
        'S1 momentum: reported Fx %.9e N is SUM mdot u_inj %.9e N within the momentum the residues keep, %.3e N'
        % (last['F'][0], Fx, bound))

    # ---- the mass each species gained over the run (the closed-domain table of the log)
    txt = open(log, errors='replace').read()
    rows = {int(a): (float(b), float(c), float(d)) for a, b, c, d in re.findall(
        r'^\s+(\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+)\s*$', txt, re.M)}
    t_src = sum(s['t1'] - s['t0'] for s in steps)
    net = sum(s['M'] * (s['t1'] - s['t0']) for s in steps)                # what the report says of the mass
    gone = A * t_src - net                                                # the oxide the caps gained
    for s in sorted(rate):
        i = order.index(s) + 1
        if i not in rows:
            ko('S1 species %s has no row in the closed-domain mass balance' % s)
            continue
        judge('S1 species %s' % s, rows[i][0], rate[s] * t_src, TOL_CAP, 'mass handed over the run vs SUM mdot ((1 - xi) g + xi c) t')
        split = rowg.get(s, 0.0) * (A * t_src - gone / nu) + rowc.get(s, 0.0) / nu * gone
        judge('S2 split %s' % s, rows[i][0], split, TOL_SPLIT, 'mass handed vs g (A - O/nu) + c O/nu, O from the reported mass')
    one = {s: A * (1.0 - xi * nu) * rowg.get(s, 0.0) for s in rowg}
    worst = min(abs(rate[s] - one[s]) / abs(rate[s]) for s in rowg)
    (ok if worst >= 100.0 * TOL_CAP else ko)('S2 non-vacuity: one row for the net mass would miss every species by at least %.3f' % worst)

    # ---- S3: the gas gains Q (t1 - t0) in each block; S4: the coupler's books
    wabs, wrel = 0.0, 0.0
    for s in steps:
        dt = s['t1'] - s['t0']
        res = abs((s['E1'] - s['E0']) - s['Q'] * dt)
        wabs, wrel = max(wabs, res), max(wrel, res / (abs(s['Q']) * dt))
    (ok if len(steps) >= 3 and wrel <= TOL_EN else ko)(
        'S3 the gas energy rises by Q (t1 - t0) in each of the %d blocks with sources: worst residual %.2e J, %.2e of Q (t1 - t0) (tol %.0e)'
        % (len(steps), wabs, wrel, TOL_EN))
    e = re.search(r'rel\.err \(bookkeeping only, datum-blind\) :\s+(\S+)', txt)
    werr = max([r[2] for r in rows.values()] or [1.0])
    (ok if len(rows) == len(rate) and werr <= TOL_BOOK_M else ko)(
        'S4 the %d species of the records close in the mass balance: worst rel.err %.3e (tol %.0e)' % (len(rows), werr, TOL_BOOK_M))
    (ok if e and float(e.group(1)) <= TOL_BOOK_E else ko)('S4 closed-domain energy bookkeeping: rel.err %s (tol %.0e)' % (e.group(1) if e else 'absent', TOL_BOOK_E))

    # ---- S5: the residues
    for i in range(n):
        d_an = d0[i] * (xi * nu * rhoA / rhoP) ** (1.0 / 3.0)
        if abs(end[i + 1][3] - d_an) > TOL_RES * d_an or abs(end[i + 1][2] - xi * nu * mdot[i]) > TOL_RES * xi * nu * mdot[i]:
            ko('S5 residue %d: d %.6e m (closed form %.6e), mass flow %.6e kg/s (xi nu mdot %.6e)' % (i + 1, end[i + 1][3], d_an, end[i + 1][2], xi * nu * mdot[i]))
            break
    else:
        ok('S5 the %d residues have the diameter d0 (xi nu rho_A/rho_OX)^(1/3) %.6e m and the mass flow xi nu mdot (tol %.0e)'
           % (n, d0[0] * (xi * nu * rhoA / rhoP) ** (1.0 / 3.0), TOL_RES))

    print('      twin: residue per kg injected %.9f kg; net mass to the gas per kg injected %.9f kg'
          % (sum(end[i + 1][2] for i in range(n)) / A, last['M'] / A))
    print('      twin: species per kg injected %s' % ' '.join('%s %+.9f' % (s, rows[order.index(s) + 1][0] / (A * t_src)) for s in sorted(rate)
                                                             if order.index(s) + 1 in rows))
    print('== check.py: %d passed, %d failed' % (npass, nfail))
    sys.exit(0 if nfail == 0 else 1)


if __name__ == '__main__':
    main()
