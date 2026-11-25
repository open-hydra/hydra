#!/usr/bin/env python3
"""igloo-two-materials: two injected IGLOO materials, each handing the gas its own species
(test/MI2/cases/igloo-two-materials: water drops evaporate into H2O, graphite particles burn to CO2 with the gas's O2).

  check.py <run-dir> [log]

run-dir holds input.ini, INPUT/ and OUTPUT/ of the run; log its log with the [MI2-report] lines (default
run-dir/run.log).
Every material injects the DB streams of [IGLOO-BC] (mass flow mdot_i, diameter d0_i, temperature T0_i). The k-th
MOSE block with IGLOO's sources (duration dt_k = t1 - t0 of its report line) carries the sources of sweep k (sweep 0
writes the untagged files): the mass rate material m loses there is IGLOO's own record,
  r_m(k) = SUM_i mdot_i (1 - rho_m(T_end) d_end^3 / (rho_m(T0_i) d0_i^3))
over the end records of OUTPUT/part-outloc-<m>[-sweep<k>].dat (burnt out, or stopped at time-end), rho_m the Density
column of the material's zone of part-properties.dat, and L_m = SUM_k r_m(k) dt_k what it lost over the run. nu(m,s) is
the row of its record in part-products.txt (kg of species s per kg of m converted).
  M1  setup: two materials, each injected (groups >= 1), each with a record; the records' species do not overlap
  M2  every sweep applied: one end record per stream and material, and both materials lose mass
  M3  per block: SUM_m (SUM_s nu(m,s)) r_m(k) = the coupler's M of the block's igloo-to-gas line
  M4  per species, the gas gained (the closed-domain balance) = SUM_m nu(m,s) L_m: each material's mass reaches its
      own species, and only those
  M5  the coupler's own audit: injected = gas gained per species
  M6  energy: the closed-domain bookkeeping (sourceEn + datum terms of both materials) against the gas's total energy
  M7  the case discriminates: the two materials lose clearly different amounts, so a swap of the material index moves
      every species total of M4 by more than its tolerance
One PASS/FAIL line per check; exit 1 on any FAIL."""
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..'))
import mi2check as eb      # test/MI2/mi2check.py: read_ini, values, read_report

# Measured on hydra 220487a (ifx), 1 thread: the run is byte-identical at 1 and 2 threads, so the floor of M3/M4 is
# the six figures IGLOO prints d_end with and the seven of the closed-domain balance (time-end 0.02 as shipped, 0.04
# and none: M3 3.5e-6, 6.7e-7, 1.7e-9; M4 1.3e-6, 1.2e-7, 2.1e-7); M6 7e-10..1.4e-9 (igloo-evap 2.9e-9).
# The material index swapped in IGLOOxMOSE (sourceMass(nm+1-imat)): M4 and M5 5.2e-1 on all three species, M6 9.6e-6.
TOL_BLOCK = 2e-5     # M3
TOL_SPECIES = 1e-5   # M4
TOL_AUDIT = 1e-10    # M5: the coupler's own audit (its printed rel.err; measured <= 1.3e-13)
TOL_E = 1e-8         # M6
MIN_SPLIT = 0.25     # M7: measured 0.47 as shipped; 0.16 at time-end 0.04, 0.00 without time-end (both burn out)
npass = nfail = 0


def ok(msg):
    global npass
    npass += 1
    print('PASS  ' + msg)


def ko(msg):
    global nfail
    nfail += 1
    print('FAIL  ' + msg)


def relerr(a, b):
    den = max(abs(a), abs(b))
    return abs(a - b) / den if den > 0.0 else 0.0


def materials(inp):
    """[(name, groups)] of INPUT/part-phase.txt."""
    rows = [l.split() for l in open(inp + 'part-phase.txt') if l.strip()]
    return [(f[0], int(f[1])) for f in rows[1:]]


def records(inp):
    """{material: (process, {species: nu})} of INPUT/part-products.txt."""
    lines = [l.split() for l in open(inp + 'part-products.txt') if l.strip() and not l.startswith('#')]
    out = {}
    for i, f in enumerate(lines):
        if f[0] == 'material':
            out[f[1]] = (f[2], {t[0]: float(t[1]) for t in lines[i + 1:i + 1 + int(f[3])]})
    return out


def density(inp):
    """{material: {T: rho}} of INPUT/part-properties.dat, the column named Density."""
    out, cur, col = {}, None, None
    for line in open(inp + 'part-properties.dat'):
        if line.startswith('VARIABLES'):
            col = [n.strip(' "') for n in line.split('=', 1)[1].split(',')].index('Density')
            continue
        if line.startswith('ZONE'):
            cur = out.setdefault(re.search(r'T="([^"]*)"', line).group(1), {})
            continue
        if cur is None:
            continue
        try:
            v = [float(x) for x in line.split()]
        except ValueError:
            continue
        if len(v) > col:
            cur[v[0]] = v[col]
    return out


def at(rows, T):
    lo = float(int(T))
    return rows[lo] if T == lo else rows[lo] + (rows[lo + 1.0] - rows[lo]) * (T - lo)


def gas_species(inp):
    return [l.split()[0] for l in open(inp + 'gas-phase.txt').read().split('\n')[1:] if l.strip()]


def balance(txt, species):
    """{species: (injected, gained, rel.err)} of the closed-domain mass balance; an untouched species: None."""
    m = re.search(r'closed-domain mass balance.*?\n.*?\n((?:\s+\d+\s+.*\n)+)', txt)
    out = {}
    if m:
        for l in m.group(1).split('\n'):
            f = l.split()
            if len(f) == 4 and f[0].isdigit():
                out[species[int(f[0]) - 1]] = (float(f[1]), float(f[2]), float(f[3]))
            elif len(f) > 1 and f[0].isdigit() and 'untouched' in l:
                out[species[int(f[0]) - 1]] = None
    return out


def outloc(path):
    """{ID: (T, d)} of the end records of an outloc file (12 columns, out-time on)."""
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) == 12:
            try:
                out[int(f[8])] = (float(f[3]), float(f[11]))
            except ValueError:
                continue
    return out


def main():
    run = sys.argv[1]
    log = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'run.log')
    inp = run + '/INPUT/'
    ini = eb.read_ini(run + '/input.ini')
    n = len(ini['IGLOO-BC']['y'].split())
    mdot = eb.values(ini, 'IGLOO-BC', 'mdot', n)
    d0 = eb.values(ini, 'IGLOO-BC', 'diam', n)
    T0 = eb.values(ini, 'IGLOO-BC', 'temp0', n)
    txt = open(log, errors='replace').read()

    mats = materials(inp)
    rec = records(inp)
    rho = density(inp)
    names = [m for m, _ in mats]
    sets = [set(rec[m][1]) if m in rec else set() for m in names]
    disjoint = len(names) == 2 and not (sets[0] & sets[1])
    good = (len(mats) == 2 and all(g >= 1 for _, g in mats) and all(m in rec for m in names) and disjoint
            and re.search(r'IGLOO materials\s+-->\s+2\b', txt) is not None)
    (ok if good else ko)('M1 two injected materials, each with its record, their species apart: %s' % '; '.join(
        '%s (%d group(s), %s: %s)' % (m, g, rec[m][0] if m in rec else 'NO RECORD',
                                      ' '.join('%s %+.6f' % kv for kv in sorted(rec[m][1].items())) if m in rec else '-')
        for m, g in mats))
    if not good:
        return

    steps = [s for s in eb.read_report(log) if 'M' in s and s['t1'] > s['t0']]
    (ok if len(steps) >= 3 else ko)('M2 %d MOSE blocks with IGLOO\'s sources (t_src %.10e s)'
                                    % (len(steps), sum(s['t1'] - s['t0'] for s in steps)))
    if not steps:
        return
    L = {m: 0.0 for m in names}
    worst_block, nrec_bad, idle = 0.0, [], []
    for k, s in enumerate(steps):
        dt = s['t1'] - s['t0']
        rk = {}
        for m in names:
            path = '%s/OUTPUT/part-outloc-%s%s.dat' % (run, m, '' if k == 0 else '-sweep%d' % k)
            end = outloc(path) if os.path.exists(path) else {}
            if sorted(end) != list(range(1, n + 1)):
                nrec_bad.append('%s sweep %d: %d end records' % (m, k, len(end)))
                rk[m] = 0.0
                continue
            r = 0.0
            for i in range(n):
                Te, de = end[i + 1]
                r += mdot[i] * (1.0 - at(rho[m], Te) * de ** 3 / (at(rho[m], T0[i]) * d0[i] ** 3))
            rk[m] = r
            if not r > 0.0:
                idle.append('%s sweep %d' % (m, k))
            L[m] += r * dt
        Mk = sum(sum(rec[m][1].values()) * rk[m] for m in names)
        worst_block = max(worst_block, relerr(Mk, s['M']))
    (ok if not nrec_bad and not idle else ko)(
        'M2 every applied sweep: %d end records per material, each material losing mass%s' % (
            n, '' if not (nrec_bad or idle) else ': ' + '; '.join(nrec_bad + ['%s loses nothing' % x for x in idle])))
    (ok if worst_block <= TOL_BLOCK else ko)(
        'M3 per block, IGLOO\'s records against the coupler\'s igloo-to-gas mass: worst rel.err %.2e (tol %.0e, '
        '%d blocks)' % (worst_block, TOL_BLOCK, len(steps)))
    for m in names:
        print('      %s lost %.10e kg over the blocks (IGLOO\'s records)' % (m, L[m]))

    species = gas_species(inp)
    bal = balance(txt, species)
    expect = {sp: sum(rec[m][1].get(sp, 0.0) * L[m] for m in names) for sp in species}
    for sp in species:
        if expect[sp] == 0.0:
            (ok if bal.get(sp, 'missing') is None else ko)(
                'M4 %s: no material converts to it, the gas holds it unchanged (%s)' % (
                    sp, 'untouched' if bal.get(sp, 'missing') is None else bal.get(sp)))
            continue
        if not bal.get(sp):
            ko('M4 %s: no row in the closed-domain mass balance' % sp)
            continue
        e = relerr(bal[sp][1], expect[sp])
        src = ' + '.join('%+.6f x %s' % (rec[m][1][sp], m) for m in names if sp in rec[m][1])
        (ok if e <= TOL_SPECIES else ko)('M4 %s gained %.6e kg = %s lost %.10e kg: rel.err %.2e (tol %.0e)'
                                         % (sp, bal[sp][1], src, expect[sp], e, TOL_SPECIES))
        (ok if bal[sp][2] <= TOL_AUDIT else ko)('M5 %s: the coupler\'s audit, injected %.6e vs gained %.6e, rel.err '
                                                '%.3e (tol %.0e)' % (sp, bal[sp][0], bal[sp][1], bal[sp][2], TOL_AUDIT))
    m6 = re.search(r'rel\.err \(bookkeeping only, datum-blind\) :\s+(\S+)', txt)
    e6 = float(m6.group(1)) if m6 else float('inf')
    (ok if e6 <= TOL_E else ko)('M6 energy handed (sourceEn and both materials\' datum terms) against the gas\'s '
                                'total energy gained: rel.err %.3e (tol %.0e)' % (e6, TOL_E))
    a, b = (L[m] for m in names)
    split = abs(a - b) / max(abs(a), abs(b))
    (ok if split >= MIN_SPLIT else ko)('M7 the materials lose %.4e and %.4e kg: they differ by %.2f of the larger '
                                       '(>= %.2f), far above M4\'s tolerance' % (a, b, split, MIN_SPLIT))


if __name__ == '__main__':
    main()
    print('== p17: %d passed, %d failed' % (npass, nfail))
    sys.exit(1 if nfail else 0)
