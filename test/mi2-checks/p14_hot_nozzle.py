#!/usr/bin/env python3
"""p14: the inputs (G1) and the run witnesses (G2) of the hot JPL nozzle cases (JPL-3solver-hot-2micron,
JPL-3solver-hot-20micron); the comparison of ICE and IGLOO is JPL-3solver-hot-2micron/compare_freezing.py.
usage: p14_hot_nozzle.py g1 <case-dir> 2um|20um
       p14_hot_nozzle.py g2 <run-dir> <log> 2um|20um
       p14_hot_nozzle.py inject <input.ini>   (the exit probes of W1b, appended to a scratch copy)
G1 (the case's INPUT/ and input.ini): (a) line 1 of part{E,L}-phase.txt 'condensed-dispersed modeling=eulerian|
   lagrangian', line 2 identical in both = 'A 1' + {solidification=on, T-melt=2327, h-fus=1090000, T-nuc=1861.6 (2um) |
   2227 (20um), cp-solid=1360}; (b) both properties tables Cp 1380 and Density 4004.62 on all 5000 rows, the files
   identical; (c) the bc files identical, every 401 record share 0.15/0.70, kV 1, kT 1, radius 1e-6 | 1e-5; (d)
   input.ini: [ICE-Physics] emissivity 0, density 4004.62, specific-heat 1380, T0 2500 | 2400 in [ICB-Block1] and [inj],
   IGLOO-iter 2, [ICE-Parameters] res-threshold 0, [ICE-IO] sol-diter 125, sol-overwrite false, no [IGLOO-Properties]
   emissivity; (e) both products files 'material A none 0 none 0.0000000000e+00' (no stoichiometry).
G2 (a run): hygiene (exit 0, GAME OVER, no failure signature, one 'All particles out of domain' per coupling step);
   solidification active (ICE f_p1, chi_p1; IGLOO f; a parcel of the last sweep with 3 plateau rows or more);
   W1 the MOSE density residual's maximum in each of the last 10 blocks <= 1e-8 of its first value;
   W1b the axis exit probe's Mach range over the last MOSE block <= 1e-7 (p18's R2; the wall's printed), not judged
       without probe files;
   W2 ICE's two last fields (ICE_ITERATION last and last - sol-diter) on the G3 cells (rows 1-10, both solvers' bulk
       density above 1e-3 of its peak): |dT_p1| <= 1e-3 K, |df_p1|, |dchi_p1| <= 1e-6, rho_p1, u_p1, v_p1, n_p1 <= 1e-4 of
       the field's maximum, no chi_p1 across 1/2; over every cell printed; ICE's residual history printed;
   W3 IGLOO's exit-plane records (bc 400) of the last sweep N against N-1 and N-5, per parcel: T within 0.01 K, f within
       1e-6 (2um) | 1e-5 (20um), the nucleation and plateau-end cells identical; 20um: the wall records' (bc 301) T too.
Each check prints PASS / FAIL / SKIP; '== G<n> n passed, m failed'; exit 1 on a failure."""
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'JPL-3solver-hot-2micron'))
import compare_freezing as cf               # noqa: E402
import jpl_chang_gate as jg                 # noqa: E402
from tec_cc import read_cc                  # noqa: E402

DESIGN = {'2um': dict(TN=1861.6, T0=2500.0, RAD=1e-6, W3F=1e-6, wall=False),
          '20um': dict(TN=2227.0, T0=2400.0, RAD=1e-5, W3F=1e-5, wall=True)}
W1_TOL, W1B_TOL, W2_T, W2_F, W2_REL, W3_T = 1e-8, jg.R2_AXIS, 1e-3, 1e-6, 1e-4, 0.01
EXTRA_SIGNS = ('no net progress', 'stuck in cell', 'Inner loop', 'marking gone')


class Count:
    def __init__(self, tag):
        self.tag, self.p, self.f = tag, 0, 0

    def __call__(self, ok, text):
        if ok is None:
            print('SKIP  %s %s' % (self.tag, text))
            return
        self.p, self.f = self.p + bool(ok), self.f + (not ok)
        print('%s  %s %s' % ('PASS' if ok else 'FAIL', self.tag, text))

    def done(self):
        print('== %s %d passed, %d failed' % (self.tag, self.p, self.f))
        return 1 if self.f else 0


# ------------------------------------------------------------------------------------------------ G1
def g1(case, design):
    D = DESIGN[design]
    want = {'solidification': 'on', 'T-melt': 2327.0, 'h-fus': 1090000.0, 'T-nuc': D['TN'], 'cp-solid': 1360.0}
    I = os.path.join(case, 'INPUT')
    ok = Count('G1')

    def lines(name):
        return open(os.path.join(I, name)).read().splitlines()

    pe, pl = lines('partE-phase.txt'), lines('partL-phase.txt')
    ok(pe[0] == 'condensed-dispersed modeling=eulerian' and pl[0] == 'condensed-dispersed modeling=lagrangian',
       '(a) line 1: %r / %r' % (pe[0], pl[0]))
    toks = pl[1].split() if len(pl) > 1 else []
    got = {}
    for t in toks[2:]:
        k, _, v = t.partition('=')
        try:
            got[k] = v if k == 'solidification' else float(v)
        except ValueError:
            got[k] = v
    ok(len(pe) > 1 and pe[1] == pl[1] and toks[:2] == ['A', '1'] and got == want,
       '(a) line 2 equal in both and = A 1 + %s: %r / %r' % (want, pe[1] if len(pe) > 1 else None,
                                                             pl[1] if len(pl) > 1 else None))
    props = []
    for n in ('partE-properties.dat', 'partL-properties.dat'):
        rows = [l.split() for l in lines(n)[4:] if l.strip()]
        props.append(len(rows) == 5000 and all(float(r[1]) == 1380.0 and float(r[2]) == 4004.62 for r in rows))
    ok(all(props) and lines('partE-properties.dat') == lines('partL-properties.dat'),
       '(b) properties: Cp 1380, Density 4004.62 on every row of both, the files identical')
    be, bl = lines('partE-bc.txt'), lines('partL-bc.txt')
    recs = [bl[i + 1].replace(',', ' ').split() for i in range(0, len(bl) - 1, 2) if bl[i].split()[-1] == '401']
    rad = sorted(set(float(r[5]) for r in recs))
    share = sorted(set(float(r[0]) for r in recs))
    k = 0.15 / 0.70
    ok(be == bl and bool(recs) and all(abs(float(r[0]) - k) <= 1e-15 * k and float(r[1]) == 1.0 and float(r[4]) == 1.0
                                      and abs(float(r[5]) - D['RAD']) <= 1e-15 * D['RAD'] for r in recs),
       '(c) bc files identical, %d records 401: share %s, radius %s (want %.3g)' % (len(recs), share, rad, D['RAD']))
    ini = open(os.path.join(case, 'input.ini')).read()

    def key(sec, kk):
        m = re.search(r'^\[%s\][^\[]*?^\s*%s\s*=\s*([^;\n]*)' % (re.escape(sec), re.escape(kk)), ini, re.M | re.S)
        return m.group(1).strip() if m else None

    d = {('ICE-Physics', 'emissivity'): '0', ('ICE-Physics', 'density'): '4004.62', ('ICE-Physics', 'specific-heat'): '1380',
         ('ICB-Block1', 'T0'): '%g' % D['T0'], ('inj', 'T0'): '%g' % D['T0'], ('HYDRA-Parameters', 'IGLOO-iter'): '2',
         ('ICE-Parameters', 'res-threshold'): '0', ('ICE-IO', 'sol-diter'): '125', ('ICE-IO', 'sol-overwrite'): 'false'}

    def same(a, b):
        if a == b:
            return True
        try:
            return float(a) == float(b)
        except (TypeError, ValueError):
            return False

    bad = {'%s %s' % sk: key(*sk) for sk, v in d.items() if not same(key(*sk), v)}
    ok(not bad and key('IGLOO-Properties', 'emissivity') is None, '(d) input.ini keys; wrong or missing: %s' % bad)
    pr = [l for n in ('partE-products.txt', 'partL-products.txt') for l in lines(n) if l.startswith('material ')]
    ok(pr == ['material A none 0 none 0.0000000000e+00'] * 2, '(e) products records %s' % pr)
    return ok.done()


# ------------------------------------------------------------------------------------------------ G2
def ini_int(run, sec, key):
    v = cf.ini_key(run, sec, key)
    return int(v) if v else None


def keyed(path, col):
    d = {}
    for line in open(path):
        t = line.split()
        if len(t) > col:
            try:
                d[int(float(t[0]))] = float(t[col])
            except ValueError:
                pass
    return d


def events(path, mesh, Tm, has_f):
    """{ID: (nucleation cell, plateau-end cell)} of a trajectory file (cells (i, j) 0-based, None when absent)."""
    cols = ('X', 'Y', 'T') + (('f',) if has_f else ())
    out = {}
    pts, owner = [], []
    for pid, rows in cf.load_rows(path, *cols).items():
        a = np.array([r[:len(cols)] for r in rows], float)
        f = a[:, 3] if has_f else np.zeros(len(a))
        k = np.flatnonzero((np.abs(a[:, 2] - Tm) <= cf.HALF_ULP) & (f > 0))
        ev = en = None
        if len(k) and k[0] > 0:
            ev = len(pts)
            pts.append(0.5 * (a[k[0] - 1, :2] + a[k[0], :2]))
            if k[-1] + 1 < len(a):
                en = len(pts)
                pts.append(0.5 * (a[k[-1], :2] + a[k[-1] + 1, :2]))
        owner.append((pid, ev, en))
    if pts:
        p = np.array(pts)
        ii, jj = mesh.locate(p[:, 0], p[:, 1])
    for pid, ev, en in owner:
        out[pid] = (None if ev is None else (int(ii[ev]), int(jj[ev])), None if en is None else (int(ii[en]), int(jj[en])))
    return out


def g2(run, log, design):
    D = DESIGN[design]
    ok = Count('G2')
    out = os.path.join(run, 'OUTPUT')
    if not os.path.isdir(out):
        print('SKIP  G2 no OUTPUT/ in %s: the run checks are not judged' % run)
        return 0
    # hygiene
    txt = open(log, errors='replace').read() if os.path.isfile(log) else ''
    m = re.findall(r'# exit=(\d+)', txt)
    bad = [t for t in jg.HYGIENE + EXTRA_SIGNS if t in txt]
    n_out, ncpl = txt.count('All particles out of domain'), ini_int(run, 'HYDRA-Parameters', 'N-coupling')
    ok(bool(txt) and (not m or m[-1] == '0') and 'GAME OVER' in txt and not bad and n_out == ncpl,
       'hygiene (%s): exit %s, GAME OVER %d, failure signatures %s, \'All particles out of domain\' %d (N-coupling %s)'
       % (log, m[-1] if m else 'not recorded', txt.count('GAME OVER'), bad or 'none', n_out, ncpl))
    # solidification active
    it, efile = cf.ice_field(run)
    E = read_cc(efile)
    sweep, F = cf.igloo_files(run)
    has_f = cf.has_column(F['traj'], 'f')
    mat = cf.material(run)
    Tm = mat['Tm']
    nplat = 0
    if has_f:
        for rows in cf.load_rows(F['traj'], 'T', 'f').values():
            if sum(1 for r in rows if abs(r[0] - Tm) <= cf.HALF_ULP and r[1] > 0) >= 3:
                nplat += 1
    ok('f_p1' in E and 'chi_p1' in E and has_f and nplat >= 1,
       'solidification active: ICE f_p1 %s, chi_p1 %s (%s); IGLOO f %s, %d parcels of sweep %d with 3 plateau rows or more'
       % ('f_p1' in E, 'chi_p1' in E, os.path.basename(efile), has_f, nplat, sweep))
    # W1 and the per-block series of the gas witnesses
    blk = ini_int(run, 'HYDRA-Parameters', 'MOSE-iter')
    res = keyed(os.path.join(out, 'gas-residual-history.dat'), 2)
    first = res[min(res)]
    nblk = max(res) // blk
    probe = {}
    for p in ('axis', 'wall'):
        f = os.path.join(out, 'exit_%s.txt' % p)
        probe[p] = keyed(f, 1) if os.path.isfile(f) else None
    print('      per MOSE block (%d iterations): axis M range | wall M range | density residual block max / first' % blk)
    series = []
    for b in range(1, nblk + 1):
        lo, hi = (b - 1) * blk, b * blk
        rng = []
        for p in ('axis', 'wall'):
            v = [x for i, x in probe[p].items() if lo < i <= hi] if probe[p] else []
            rng.append(max(v) - min(v) if v else None)
        w1 = max(r for i, r in res.items() if lo < i <= hi) / first
        series.append((b, rng[0], rng[1], w1))
        print('      %4d | %10s | %10s | %.3e' % (b, '%.3e' % rng[0] if rng[0] is not None else '-',
                                                  '%.3e' % rng[1] if rng[1] is not None else '-', w1))
    w1w = max(s[3] for s in series[-10:])
    ok(w1w <= W1_TOL, 'W1 MOSE density residual: worst block maximum of the last 10 blocks %.3e of the first value <= %.0e'
       % (w1w, W1_TOL))
    if series[-1][1] is None:
        ok(None, 'W1b no probe file, not judged')
    else:
        tail = [(b, np.log10(a)) for b, a, _, _ in series[-10:] if a]
        sl = np.polyfit([b for b, _ in tail], [l for _, l in tail], 1)[0] if len(tail) > 1 else float('nan')
        ok(series[-1][1] <= W1B_TOL, 'W1b axis exit Mach range over the last block %.3e <= %.0e (wall %.3e; last 10 blocks '
           'axis %.3e .. %.3e, slope of log10 %.3f per block, printed)' % (series[-1][1], W1B_TOL, series[-1][2],
                                                                           min(s[1] for s in series[-10:]),
                                                                           max(s[1] for s in series[-10:]), sl))
    # W2
    sd = ini_int(run, 'ICE-IO', 'sol-diter')
    fields = {}
    for f in glob.glob(os.path.join(out, 'partE-field*.tec')):
        fields.setdefault(cf.trailer(f, 'ICE_ITERATION'), f)
    prev = it - sd
    if prev not in fields:
        ok(None, 'W2 no ICE field at ICE_ITERATION %d beside %d (fields at %s): not judged'
           % (prev, it, sorted(k for k in fields if k is not None)))
    else:
        A, B = read_cc(fields[prev]), E
        L = read_cc(F['euler'])
        rE = B['rho_p1'][0]
        rL = L[[n for n in L if n.startswith('rho')][0]][0]
        g3 = np.zeros_like(rE, bool)
        g3[:cf.INNER] = True
        g3 &= (rE > cf.MASS * rE.max()) & (rL > cf.MASS * rL.max())
        parts, okw = [], True
        for v, tol, kind in (('T_p1', W2_T, 'abs'), ('f_p1', W2_F, 'abs'), ('chi_p1', W2_F, 'abs'), ('rho_p1', W2_REL, 'rel'),
                             ('u_p1', W2_REL, 'rel'), ('v_p1', W2_REL, 'rel'), ('n_p1', W2_REL, 'rel')):
            if v not in B or v not in A:
                parts.append('%s absent' % v)
                continue
            dlt = np.abs(B[v][0] - A[v][0])
            if kind == 'rel':
                dlt = dlt / max(np.abs(B[v][0]).max(), 1e-300)
            g, a = float(dlt[g3].max()), float(dlt.max())
            okw = okw and g <= tol
            parts.append('%s %.2e (all cells %.2e)' % (v, g, a))
        if 'chi_p1' in B and 'chi_p1' in A:
            cross = int((((A['chi_p1'][0] >= 0.5) != (B['chi_p1'][0] >= 0.5)) & g3).sum())
            okw = okw and cross == 0
            parts.append('chi across 1/2: %d' % cross)
        ok(okw, 'W2 ICE ICE_ITERATION %d against %d on the %d G3 cells (T <= %.0e K, f, chi <= %.0e, others <= %.0e of the '
           'maximum): %s' % (it, prev, int(g3.sum()), W2_T, W2_F, W2_REL, '; '.join(parts)))
    eres = keyed(os.path.join(out, 'partE-residual-history.dat'), 2)
    iblk = ini_int(run, 'HYDRA-Parameters', 'ICE-iter')
    if eres and iblk:
        n = max(eres)
        last = [eres[i] for i in sorted(eres) if i > n - iblk]
        print('      ICE density residual (printed): first %.4e, last block max %.4e min %.4e, last %.4e'
              % (eres[min(eres)], max(last), min(last), eres[n]))
    # W3
    if sweep < 5:
        ok(None, 'W3 fewer than 6 sweeps: not judged')
    else:
        mesh = cf.Mesh(E['x'][0], E['y'][0])

        def tagged(kind, n):
            return os.path.join(out, 'partL-%s-A%s.dat' % (kind, '-sweep%d' % n if n else ''))

        def recs(n, bc):
            fcol = cf.has_column(tagged('outloc', n), 'f')
            r = cf.load_rows(tagged('outloc', n), 'T', 'bc', *(('f',) if fcol else ()))
            return {pid: (rows[0][0], rows[0][2] if fcol else None) for pid, rows in r.items() if int(rows[0][1]) == bc}, fcol

        cur, fcol = recs(sweep, 400)
        ev_cur = events(tagged('trajectories', sweep), mesh, Tm, has_f) if has_f else None
        for back in (1, 5):
            old, _ = recs(sweep - back, 400)
            common = sorted(set(cur) & set(old))
            dT = max(abs(cur[p][0] - old[p][0]) for p in common) if common else float('nan')
            okT = common and set(cur) == set(old) and dT <= W3_T
            txt3 = ('sweep %d against %d: %d exit-plane records (%d in both), max |dT| %.2e K <= %.2f'
                    % (sweep, sweep - back, len(cur), len(common), dT, W3_T))
            okF = okE = True
            if fcol and has_f:
                dF = max(abs(cur[p][1] - old[p][1]) for p in common) if common else float('nan')
                okF = dF <= D['W3F']
                ev_old = events(tagged('trajectories', sweep - back), mesh, Tm, has_f)
                moved = [p for p in ev_cur if p in ev_old and ev_cur[p] != ev_old[p]]
                okE = not moved
                txt3 += ', max |df| %.2e <= %.0e, parcels whose nucleation or plateau-end cell moved %d' % (dF, D['W3F'],
                                                                                                            len(moved))
            else:
                txt3 += '; f and event cells: SKIP (no f column)'
            if D['wall']:
                cw, _ = recs(sweep, 301)
                ow, _ = recs(sweep - back, 301)
                cm = sorted(set(cw) & set(ow))
                dw = max(abs(cw[p][0] - ow[p][0]) for p in cm) if cm else 0.0
                okT = okT and set(cw) == set(ow) and dw <= W3_T
                txt3 += '; wall records %d, max |dT| %.2e K' % (len(cw), dw)
            ok(bool(okT and okF and okE), 'W3 IGLOO ' + txt3)
    return ok.done()


def inject(ini):
    """p18's two exit probes (jpl_chang_gate.PROBE_CELLS: the Mach number of the last cell column, axis and wall rows)
    appended to a scratch input.ini, PROBE_PER_BLOCK rows per MOSE block, else 10 (a block of 250 iterations: every 25)."""
    txt = open(ini).read()
    if re.search(r'^\s*\[MOSE-Probes\]', txt, re.M):
        sys.exit('%s: has a [MOSE-Probes] section already' % ini)
    blk = int(re.search(r'^\s*MOSE-iter\s*=\s*(\d+)', txt, re.M).group(1))
    per = next((n for n in (jg.PROBE_PER_BLOCK, 10) if blk % n == 0), None)
    if per is None:
        sys.exit('%s: MOSE-iter %d is a multiple of neither %d nor 10' % (ini, blk, jg.PROBE_PER_BLOCK))
    out = '\n[MOSE-Probes]\n' + ''.join('probe%d = %s\n' % (n + 1, p) for n, (p, _) in enumerate(jg.PROBE_CELLS))
    for p, cell in jg.PROBE_CELLS:
        out += '\n[%s]\nvariables = M\nindex-position = %s\nditer = %d\n' % (p, cell, blk // per)
    open(ini, 'a').write(out)
    return 0


def main():
    a = sys.argv[1:]
    if len(a) == 2 and a[0] == 'inject':
        sys.exit(inject(a[1]))
    if len(a) == 3 and a[0] == 'g1':
        sys.exit(g1(a[1], a[2]))
    if len(a) == 4 and a[0] == 'g2':
        sys.exit(g2(a[1], a[2], a[3]))
    sys.exit(__doc__)


if __name__ == '__main__':
    main()
