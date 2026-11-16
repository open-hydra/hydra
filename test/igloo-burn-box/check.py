#!/usr/bin/env python3
"""igloo-burn-box: the mass, the momentum, the heat of reaction and the energy of a burning IGLOO material reach
the gas.

  python3 check.py RUN_DIR LOG [EVAP_BOX_MOM_DIR [law]]

RUN_DIR holds input.ini, INPUT/ and OUTPUT/ of the run; LOG is its log with the [MI2-report] lines; the third
argument is test/evap-box-mom, whose check.py gives the readers (default: ../evap-box-mom beside this file).

Two cases share this check: test/igloo-burn-box (test/evap-box-mom with a gas of one species and its drops burning,
combustion = Beckstead, into that species: a-n2.yaml) and test/igloo-burn-box-rx (graphite drops burning to CO2 in
N2, O2 and CO2: c-o2.yaml). The record of INPUT/part-products.txt is "formation" in both. Every stream burns out inside
the box and inside time-end, so everything it carries in ends in the gas, whatever the burning law does on the way
(the share beta-part of the heat that the law leaves on the drop returns with it). Closed forms, from input.ini and
the INPUT files, per stream of mass flow mdot injected at u_inj and T0 (h_c the material's enthalpy table of
INPUT/part-properties.dat, T_ref = 298.15 K):
  M  = SUM mdot                                         mass to the gas                    [kg/s]
  Fx = SUM mdot u_inj                                   force on the gas                   [N]
  dh = h298 of the record - h_c(T_ref)                  the burnt mass on the gas's enthalpy scale [J/kg]
  Q  = SUM mdot (h_c(T0) + u_inj^2/2 + dh)              energy to the gas [W]
and over a MOSE block [t0, t1] with these sources in force the gas energy rises by Q (t1 - t0). The gas gains
nu(s) M of each species s of the record, and warms by the heat of reaction Q_r = h298 - SUM nu h_s(T_ref).
With a residue (d-stop > 0: test/igloo-burn-box-beckstead, test/igloo-burn-box-direct) a stream of mass flow mdot
leaves r mdot in the box as inert drops, r = (d_end/d0)^3, which end at time-end at T_end (OUTPUT/part-outloc-*.dat):
  M  = SUM mdot (1 - r)        Q = SUM mdot [h_c(T0) + u_inj^2/2 - r h_c(T_end) + (1 - r) dh]
and Fx = SUM mdot u_inj to within r (the residue ends at the gas's speed).
The word "law" as fourth argument checks the law itself on a run whose time-end falls inside the burn and whose
drops enter almost at rest: the diameter of each parcel at the end of the first sweep (the gas still at its initial
state), against
  Beckstead   d^n = d0^n - K-burn X-eff t
  direct      d^2 = d0^2 - K t,  K = 4 Nu k_g ln(1 + B) / (rho c_pg),  B = (c_pg (T_g - T_p) + q) / L_v,  Nu = 2
              (Ranz-Marshall at the injection speed is above 2 by 0.6 Re^1/2 Pr^1/3: the tolerance), T_p held.
D6 (a gas of one species G, which the record gives) is the warming of the rigid box by the first law, with the
table's internal energy e = h - R T (c_v, not c_p): (m0 + M t) e(T1) = m0 e(T0) + E_in, solved for T1 and compared
with the mass-weighted temperature the log prints. E_in is written from the record, not from the coupler's offset:
1 kg of the material at T_ref carries h298 on the gas's scale, so at T h298 + h_c(T) - h_c(T_ref). Per kg burnt at
T_ref the gas warms by Q_r + R T, not by Q_r.
One PASS/FAIL line per check with the numbers and the tolerance; exit 1 on any FAIL."""
import glob
import math
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
mom = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(here, '..', 'evap-box-mom')
sys.path.insert(0, mom)
import check as eb          # test/evap-box-mom/check.py: read_ini, values, read_report

TOL_RATE = 1e-12    # D1  M, Fx and Q of the report against the closed forms (telescoped sums: round-off)
TOL_REC = 1e-10     # D1  the same for a capped material handed whole to the gas (residue = off): its mass is the sum
#                         of three rows of the products file, which prints eleven figures (3e-12 measured)
TOL_EN = 1e-5       # D2  E1 - E0 = Q (t1 - t0) per block, against the larger of Q and the heat-of-reaction rate
#                         Q_r M: evap-box-mom's M3a tolerance (MOSE's energy round trip loses about 1e-10 J per block)
TOL_BOOK_M = 1e-9   # D3  the coupler's closed-domain mass balance (the gas gains 1e-4 of its mass: round-off 1e-12)
TOL_BOOK_E = 1e-5   # D3  the coupler's closed-domain energy bookkeeping (as D2)
TOL_WARM = 1e-2     # D5  m_gas c_v dT against the heat of reaction (frozen c_v; the mass arrives 30 K colder: 1e-3)
TOL_WARM_X = 5e-5   # D6  warming of a one-species case against the first law of the rigid box: the log prints the
#                         final gas mass with five figures (1.2e-5) and two temperatures with six decimals (9e-6 of
#                         a rise of 0.11 K)
TOL_RES = 1e-9      # D1  with a residue: d_end and T_end are read with seven figures (3 r 5e-7 = 2e-10 of M)
TOL_LAW = 2e-6      # L   d at time-end against the law: seven printed figures of d, IGLOO's step control
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
    """Column col of a table at T: linear between the integer-kelvin rows."""
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
    eb.die('no combustion record of %s in %s' % (material, path))


def outloc(path):
    """{ID: (T, t, d)} of the end records of a part-outloc file."""
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) == 12:
            try:
                out[int(f[8])] = (float(f[3]), float(f[9]), float(f[11]))
            except ValueError:
                continue
    return out


def law(run, ini, material, model, d0, rho):
    """The diameter each parcel has at time-end of the first sweep, against the burning law."""
    end = outloc(os.path.join(run, 'OUTPUT', 'part-outloc-%s.dat' % material))   # the first sweep's file has no number
    g, b = ini['GPB-Phase2'], ini['IGLOO-BC']
    Tg, pg = float(ini['ICB-Block1']['T']), float(ini['ICB-Block1']['p'])
    u0 = float(b['up'].split()[0])
    (ok if len(end) == len(d0) else ko)('L  %d end records for %d streams' % (len(end), len(d0)))
    worst, kmean = 0.0, 0.0
    for i, d in enumerate(d0):
        Te, te, de = end[i + 1]
        if model == 'Beckstead':
            n, KX = float(g.get('n-burn', 1.8)), float(g['K-burn']) * float(g.get('X-eff', 1.0))
            d_an = (d ** n - KX * te) ** (1.0 / n)
            k_meas, tol = (d ** n - de ** n) / te, TOL_LAW
        else:
            p = ini.get('IGLOO-Properties', {})
            Lv, Tb, q = float(p.get('Lv')), float(g.get('T-burn', p.get('boiling-temperature'))), float(g['q-comb'])
            tt, _ = zones(os.path.join(run, 'INPUT', 'gas-transport.dat'))
            th, _ = zones(os.path.join(run, 'INPUT', 'gas-thermo.dat'))
            sp = ini['GPB-Phase1']['species'].split()[0]
            mu, kg, cpg = at(tt[sp], Tg, 1), at(tt[sp], Tg, 2), at(th[sp], Tg, 1)
            Rg = 8314.46 / 28.0134      # R of N2, for the Reynolds number of the bound only
            Re = pg / (Rg * Tg) * u0 * d / mu
            dNu = 0.6 * math.sqrt(Re) * (mu * cpg / kg) ** (1.0 / 3.0)
            B = (cpg * (Tg - Tb) + q) / Lv
            K = 8.0 * kg * math.log1p(B) / (rho * cpg)
            d_an = math.sqrt(d * d - K * te)
            k_meas, tol = (d * d - de * de) / te, 0.5 * dNu + TOL_LAW
            (ok if abs(Te - Tb) <= 1e-6 else ko)('L  parcel %d keeps its burning temperature: %.6f K at the end, %.6f K' % (i + 1, Te, Tb))
        e = abs(de - d_an) / d_an
        worst, kmean = max(worst, e), kmean + k_meas / len(d0)
        (ok if e <= tol and d_an < 0.95 * d else ko)(
            'L  parcel %d at t %.9e s: d %.6e m vs the law %.9e m, rel.err %.2e (tol %.1e)' % (i + 1, te, de, d_an, e, tol))
    print('      twin: law constant %.9e (Beckstead: K-burn X-eff [m^n/s]; direct: K [m2/s]) from the diameters at time-end' % kmean)
    print('== check.py: %d passed, %d failed' % (npass, nfail))
    sys.exit(0 if nfail == 0 else 1)


def main():
    run, log = os.path.abspath(sys.argv[1]), sys.argv[2]
    ini = eb.read_ini(os.path.join(run, 'input.ini'))
    n = len(ini['IGLOO-BC']['y'].split())
    mdot = eb.values(ini, 'IGLOO-BC', 'mdot', n)
    up = eb.values(ini, 'IGLOO-BC', 'up', n)
    T0 = eb.values(ini, 'IGLOO-BC', 'temp0', n)
    material = ini['GPB-Phase2'].get('material', 'A').split()[0]
    species = ini['GPB-Phase1']['species'].split()
    energy, h298, nu = record(os.path.join(run, 'INPUT', 'part-products.txt'), material)
    model = ini['GPB-Phase2'].get('combustion', 'none').split()[0]
    (ok if model in ('Beckstead', 'direct') else ko)('C  the material burns with combustion = %s' % model)
    d0 = eb.values(ini, 'IGLOO-BC', 'diam', n)
    if len(sys.argv) > 4 and sys.argv[4] == 'law':
        law(run, ini, material, model, d0, float(ini['GPB-Phase2']['rho']))
    d_res = float(ini['GPB-Phase2'].get('d-stop', 0.0))
    (ok if energy == 'formation' and all(s in species for s in nu) else ko)(
        'C  the record of %s is "%s" and its species (%s) are species of the gas (%s)'
        % (material, energy, ' '.join(nu), ' '.join(species)))

    ptab, pnames = zones(os.path.join(run, 'INPUT', 'part-properties.dat'))
    hcol = [i for i, v in enumerate(pnames) if v.startswith('Enthalpy')][0]
    hc = ptab[material]
    gtab, gnames = zones(os.path.join(run, 'INPUT', 'gas-thermo.dat'))
    gcol = [i for i, v in enumerate(gnames) if v.startswith('Enthalpy') or v == 'h'][0] if any(
        v.startswith('Enthalpy') or v == 'h' for v in gnames) else 2
    hgas = {s: at(gtab[s], T_REF, gcol) for s in nu}
    q = 0.0
    dh = h298 - at(hc, T_REF, hcol)
    Qr = h298 - sum(nu[s] * hgas[s] for s in nu)
    r, Tend, tol_m, tol_f = [0.0] * n, list(T0), TOL_RATE, TOL_RATE
    if any(l.split()[:3] == ['material', material, 'cap'] for l in open(os.path.join(run, 'INPUT', 'part-products.txt'))) \
            and 'xi-cap' in ini['GPB-Phase2']:
        tol_m = TOL_REC
    if d_res > 0.0:
        # the end records of the sweep whose sources the last report carries: the k-th block with sources runs on
        # sweep k (part-outloc-<m>.dat is sweep 1, part-outloc-<m>-sweep<i>.dat sweep i + 1)
        nsrc = len([s for s in eb.read_report(log) if 'F' in s])
        name = 'part-outloc-%s.dat' % material if nsrc <= 1 else 'part-outloc-%s-sweep%d.dat' % (material, nsrc - 1)
        end = outloc(os.path.join(run, 'OUTPUT', name))
        r = [(end[i + 1][2] / d0[i]) ** 3 for i in range(n)]
        Tend = [end[i + 1][0] for i in range(n)]
        tol_m, tol_f = TOL_RES, max(r)
        (ok if all(abs(end[i + 1][2] - d_res) <= 1e-6 * d_res for i in range(n)) else ko)(
            'C  every parcel ends at the residue diameter d-stop %.6e m (r = %.6e of its mass)' % (d_res, max(r)))
    M_an = sum(m * (1.0 - x) for m, x in zip(mdot, r))
    F_an = sum(m * u for m, u in zip(mdot, up))
    Q_an = sum(m * (at(hc, t, hcol) + 0.5 * u * u - x * at(hc, te, hcol) + (1.0 - x) * (dh + q))
               for m, u, t, x, te in zip(mdot, up, T0, r, Tend))
    Q_res = sum(m * x * (at(hc, t, hcol) - at(hc, te, hcol) + 0.5 * u * u) for m, u, t, x, te in zip(mdot, up, T0, r, Tend))
    Q_dry = sum(m * (at(hc, t, hcol) + 0.5 * u * u) for m, u, t in zip(mdot, up, T0))
    Q_dry = Q_dry - M_an * Qr
    print('      closed forms: M %.15e kg/s, Fx %.15e N, dh %.9e J/kg, Q_r %.9e J/kg, Q %.15e W (without the heat of '
          'reaction %.6e W)' % (M_an, F_an, dh, Qr, Q_an, Q_dry))

    steps = [s for s in eb.read_report(log) if 'F' in s]
    if not steps:
        ko('D1 the report has no igloo-to-gas line')
        print('== check.py: %d passed, %d failed' % (npass, nfail))
        sys.exit(1)
    worst = {'M': 0.0, 'F': 0.0, 'Q': 0.0, 'E': 0.0}
    for s in steps:
        worst['M'] = max(worst['M'], abs(s['M'] - M_an) / M_an)
        worst['F'] = max(worst['F'], abs(s['F'][0] - F_an) / F_an)
        worst['Q'] = max(worst['Q'], abs(s['Q'] - Q_an) / abs(Q_an))
    last = steps[-1]
    judge('D1 mass', last['M'], M_an, tol_m, 'reported M vs SUM mdot (1 - r)')
    judge('D1 momentum', last['F'][0], F_an, tol_f, 'reported Fx vs SUM mdot u_inj')
    judge('D1 energy', last['Q'], Q_an, tol_m, 'reported Q vs SUM mdot (h_c(T0) + u^2/2 - r h_c(T_end) + (1 - r) dh)')
    if d_res == 0.0:
        (ok if max(worst['M'], worst['F'], worst['Q']) <= tol_m else ko)(
            'D1 every one of the %d reports: worst rel.err M %.2e, Fx %.2e, Q %.2e (tol %.0e)'
            % (len(steps), worst['M'], worst['F'], worst['Q'], tol_m))
    print('      twin: energy per kg burnt (Q - residue)/M %.9f J/kg' % ((last['Q'] - Q_res) / last['M']))

    # the report of a step gives the rates in force during its block (IGLOO solves between blocks, relaxation off)
    nblocks = len(steps)
    wabs = 0.0
    for s_ in steps:
        dt = s_['t1'] - s_['t0']
        res = abs((s_['E1'] - s_['E0']) - s_['Q'] * dt)
        wabs = max(wabs, res)
        worst['E'] = max(worst['E'], res / (max(abs(s_['Q']), abs(Qr) * M_an) * dt))
    (ok if nblocks >= 3 and worst['E'] <= TOL_EN else ko)(
        'D2 the gas energy rises by Q (t1 - t0) in each of the %d blocks with sources: worst residual %.2e J, %.2e of '
        'max(Q, Q_r M) (t1 - t0) (tol %.0e)' % (nblocks, wabs, worst['E'], TOL_EN))

    txt = open(log, errors='replace').read()
    m = re.search(r'closed-domain mass balance.*?\n.*?\n\s+1\s+\S+\s+\S+\s+(\S+)', txt)
    e = re.search(r'rel\.err \(bookkeeping only, datum-blind\) :\s+(\S+)', txt)
    if not m or not e:
        ko('D3 the log carries no closed-domain mass balance or energy bookkeeping')
    else:
        (ok if float(m.group(1)) <= TOL_BOOK_M else ko)('D3 closed-domain mass balance of species 1: rel.err %.3e (tol %.0e)'
                                                        % (float(m.group(1)), TOL_BOOK_M))
        (ok if float(e.group(1)) <= TOL_BOOK_E else ko)('D3 closed-domain energy bookkeeping: rel.err %.3e (tol %.0e)'
                                                        % (float(e.group(1)), TOL_BOOK_E))
    (ok if abs(M_an * Qr) >= 0.1 * abs(Q_dry) else ko)(
        'D4 non-vacuity: the heat of reaction is %.3f of the energy the gas receives without it' % (abs(M_an * Qr) / abs(Q_dry)))
    if True:
        # the species of the record: every row of the closed-domain mass balance, and the ratio of the masses handed
        rows = re.findall(r'^\s+(\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+)\s*$', txt, re.M)
        inj = {int(r[0]): float(r[1]) for r in rows}
        werr = max([float(r[3]) for r in rows] or [1.0])
        (ok if len(rows) == len(nu) and werr <= TOL_BOOK_M else ko)(
            'D4 the %d species of the record close in the mass balance: worst rel.err %.3e (tol %.0e)' % (len(rows), werr, TOL_BOOK_M))
        order = sorted(nu, key=lambda s: -nu[s])
        if len(inj) == 2:
            a, b = max(inj.values()), min(inj.values())
            judge('D4 species', a / b, nu[order[0]] / nu[order[-1]], 1e-6, 'mass produced / mass consumed (seven printed figures each) vs the row (%s / %s)'
                  % (order[0], order[-1]))
        # the heat of reaction warms the gas: the log's mass-weighted temperatures, frozen c_v and gas mass
        g = re.search(r'mass-weighted gas T:\s+(\S+)\s+->\s+(\S+)', txt)
        h = re.search(r'mass-weighted frozen cv \[J/kg/K\] :\s+(\S+)\s+gas mass \[kg\] :\s+(\S+)', txt)
        t_src = sum(s_['t1'] - s_['t0'] for s_ in steps)
        if len(species) == 1:
            pass        # a gas of one species: D6, the first law, in place of this estimate with a frozen c_v
        elif not g or not h:
            ko('D5 the log carries no gas temperature diagnostic')
        else:
            warm = float(h.group(2)) * float(h.group(1)) * (float(g.group(2)) - float(g.group(1)))
            judge('D5 heat of reaction', warm, Qr * M_an * t_src, TOL_WARM,
                  'm_gas c_v dT vs Q_r M t (Q_r %.6e J/kg, dT %.3f K)' % (Qr, float(g.group(2)) - float(g.group(1))))

    if len(species) == 1:
        g = re.search(r'mass-weighted gas T:\s+(\S+)\s+->\s+(\S+)', txt)
        h = re.search(r'mass-weighted frozen cv \[J/kg/K\] :\s+(\S+)\s+gas mass \[kg\] :\s+(\S+)', txt)
        if not g or not h:
            ko('D6 the log carries no gas temperature diagnostic')
        else:
            # the gas mass and the frozen c_v of the log are those of the final state
            Tg0, Tg1, cv, m1 = float(g.group(1)), float(g.group(2)), float(h.group(1)), float(h.group(2))
            sp = species[0]
            Rg = at(gtab[sp], Tg1, 1) - cv
            t_src = sum(s_['t1'] - s_['t0'] for s_ in steps)
            m0 = m1 - M_an * t_src
            E_in = t_src * sum(m * (at(hc, t, hcol) + 0.5 * u * u - x * at(hc, te, hcol)
                                    + (1.0 - x) * (h298 - at(hc, T_REF, hcol)))
                               for m, u, t, x, te in zip(mdot, up, T0, r, Tend))

            def rise(internal):
                def e_gas(Tk):
                    return at(gtab[sp], Tk, gcol) - (Rg * Tk if internal else 0.0)
                target = m0 * e_gas(Tg0) + E_in
                lo, hi = Tg0 - 50.0, Tg0 + 500.0
                for _ in range(200):
                    mid = 0.5 * (lo + hi)
                    lo, hi = (mid, hi) if m1 * e_gas(mid) < target else (lo, mid)
                return 0.5 * (lo + hi) - Tg0
            dT_v, dT_p = rise(True), rise(False)
            judge('D6 warming', Tg1 - Tg0, dT_v, TOL_WARM_X,
                  'gas temperature rise of the log vs the first law of the rigid box (e = h - R T)')
            (ok if abs(dT_p - dT_v) >= 100.0 * TOL_WARM_X * dT_v else ko)(
                'D6 non-vacuity: the same oracle at constant pressure (h in place of e) gives %.6f K, %.2f %% away'
                % (dT_p, 100.0 * abs(dT_p - dT_v) / dT_v))

    print('== check.py: %d passed, %d failed' % (npass, nfail))
    sys.exit(0 if nfail == 0 else 1)


if __name__ == '__main__':
    main()
