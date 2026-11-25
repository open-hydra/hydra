#!/usr/bin/env python3
"""Refutation of compare_freezing.py: each check must fail on a run doctored against it, and only that check.
usage: refute_compare.py <run-dir> <scratch> [--design full|partial]

<run-dir> is a correct run (compare_freezing.py passes it). In <scratch>, one copy per check holds input.ini, INPUT/
(linked) and the five output files compare_freezing.py reads (linked, but the one doctored); compare_freezing.py runs on
the untouched copy (every check PASS) and on each doctored copy, where the targeted check must FAIL and no other. The
size of every doctoring is derived from what compare_freezing.py measures on the untouched run (its tolerances, N_front,
the measured offsets), never a fixed number:
  T1  IGLOO's Eulerian T_p in T1's cells moved toward T0 by (tolerance + delta_hot + 0.01) (T0 - T_p)
  T2  ICE's front moved upstream in rows 1-10 by N_front + 1 cells (plus the largest offset of the other sign): the cells
      before it on the plateau (chi_p1 1, T_p1 T_melt, f_p1 on the row's plateau line)
  T3  ICE's T_p1 x (1 + 1e-6) in its plateau cells; T3-parcels: IGLOO's T x (1 + 1e-6) on one inner parcel's plateau
      rows but its first and last (its events, and every other parcel, unchanged)
  T4  (full) each inner parcel's plateau ends N_front + 2 cells before ICE's (its later plateau rows solid at T_melt - 1)
  T5  (full) an exit-plane record outside the inner region at T_melt + 3 K; an ICE exit cell outside it at f_p1 0.9;
      a wall record at T0 - 1 K must change nothing
  T5p (partial) an exit-plane record outside the inner region at T_melt - 1 K; an ICE exit cell outside it at f_p1 1
  T6  (full) IGLOO's inner exit-plane temperatures moved by the tolerance + the measured difference + 1 K
  C1  IGLOO's plateau f scaled about f0 so that the row with the largest ratio leaves C1's 2 % by 1 %
  C2  (full) an exit-plane record at f 0.99; C2p (partial) IGLOO's inner exit f moved by the tolerance + the measured
      difference + 0.001
Prints one line per copy and '== refute_compare: n passed, m failed'; exit 1 on any failure."""
import argparse
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_freezing as cf          # noqa: E402

CF = os.path.join(HERE, 'compare_freezing.py')


# ------------------------------------------------------------------------------------------------ Tecplot BLOCK files
class Tec:
    """A Tecplot BLOCK file as written by ICE and IGLOO: the header and the trailer kept verbatim, the values as one array."""

    def __init__(self, path):
        text = open(path).read()
        z = re.search(r'^\s*ZONE.*$', text, re.M | re.I)
        self.head = text[:z.end() + 1]
        rest = text[z.end() + 1:]
        m = re.search(r'^\s*DATASETAUXDATA', rest, re.M)
        body, self.tail = (rest[:m.start()], rest[m.start():]) if m else (rest, '')
        self.vals = np.array([float(t) for t in body.split()])
        I, J, K = (int(re.search(r'\b%s\s*=\s*(\d+)' % a, z.group(0), re.I).group(1)) for a in 'IJK')
        self.nn, self.nc, self.shape = I * J * K, (I - 1) * (J - 1) * max(1, K - 1), (max(1, K - 1), J - 1, I - 1)
        names = re.findall(r'"[^"]+"|[^\s,]+', re.search(r'VARIABLES\s*=\s*(.*)', text).group(1))
        self.names = [n.strip('"') for n in names]

    def var(self, name):
        v = self.names.index(name) - 3
        s = 3 * self.nn + v * self.nc
        return self.vals[s:s + self.nc].reshape(self.shape)[0]     # a view: (J-1, I-1)

    def write(self, path):
        with open(path, 'w') as f:
            f.write(self.head)
            f.write(''.join(' %.17E\n' % v for v in self.vals))
            f.write(self.tail)


# ------------------------------------------------------------------------------------------------ IGLOO point files
def edit_points(src, dst, edit):
    """Copy an IGLOO particle file, calling edit(cols, tokens, k) on each data line (cols: name -> index; k: the
    line's ordinal within its parcel); edit returns the new tokens or None for the line unchanged."""
    cols, seen, out = None, {}, []
    for line in open(src):
        low = line.lstrip().lower()
        if low.startswith('variables'):
            names = re.findall(r'"([^"]*)"', line)
            cols = {n: i for i, n in enumerate(names)}
            out.append(line)
            continue
        t = line.split()
        if cols is None or len(t) < len(cols) or low.startswith('zone'):
            out.append(line)
            continue
        try:
            pid = int(t[cols['ID']])
            float(t[0])
        except ValueError:
            out.append(line)
            continue
        k = seen.get(pid, 0)
        seen[pid] = k + 1
        new = edit(cols, t, k)
        out.append(line if new is None else ' ' + ' '.join(new) + '\n')
    open(dst, 'w').write(''.join(out))


# ------------------------------------------------------------------------------------------------ copies
def files(run):
    _, ef = cf.ice_field(run)
    _, F = cf.igloo_files(run)
    return dict(ice=ef, traj=F['traj'], exit=F['exit'], euler=F['euler'])


def make_copy(run, scratch, name):
    d = os.path.join(scratch, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, 'OUTPUT'))
    shutil.copy2(os.path.join(run, 'input.ini'), d)
    os.symlink(os.path.abspath(os.path.join(run, 'INPUT')), os.path.join(d, 'INPUT'))
    for f in files(run).values():
        os.symlink(os.path.abspath(f), os.path.join(d, 'OUTPUT', os.path.basename(f)))
    return d


def target(d, f):
    """The path to write a doctored file to (the link replaced)."""
    p = os.path.join(d, 'OUTPUT', os.path.basename(f))
    os.remove(p)
    return p


def verdicts(d, design):
    p = subprocess.run([sys.executable, CF, d, '--design', design], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    v = {}
    for line in p.stdout.splitlines():
        m = re.match(r'^(PASS|FAIL|SKIP)\s+(\S+)\s', line)
        if m:
            v[m.group(2)] = m.group(1)
    return v, p.stdout


# ------------------------------------------------------------------------------------------------ doctorings
def doctor_T1(run, d, R):
    tf = files(run)['euler']
    t = Tec(tf)
    T = t.var([n for n in t.names if n.startswith('T')][0])
    a = cf.DELTA_TOL + R['delta_hot'] + 0.01
    m = R['t1_mask']
    T[m] += a * (R['T0'] - T[m])
    t.write(target(d, tf))
    return 'IGLOO T_p + %.4f (T0 - T_p) in %d cells' % (a, int(m.sum()))


def plateau_line(fE, xc, j):
    ic = np.flatnonzero((fE[j] > 0) & (fE[j] < 1))
    ic = ic[1:-1] if len(ic) >= 3 else ic
    return np.polyfit(xc[j, ic], fE[j, ic], 1)


def doctor_T2(run, d, R, mesh):
    tf = files(run)['ice']
    t = Tec(tf)
    chi, f, T = t.var('chi_p1'), t.var('f_p1'), t.var('T_p1')
    fE = f.copy()
    note = []
    for j in range(cf.INNER):
        iF, N = R['ice_front'][j], R['n_front'][j]
        gaps = [i - iF for _, i, jj in R['events'] if jj == j]
        s = N + 1 + max(0, -min(gaps)) if gaps else N + 1
        lo = max(0, iF - s)
        line = plateau_line(fE, mesh.xc, j)
        for i in range(lo, iF):
            chi[j, i] = 1.0
            T[j, i] = R['mat']['Tm']
            f[j, i] = min(max(np.polyval(line, mesh.xc[j, i]), 1e-6), 1 - 1e-6)
        note.append(s)
    t.write(target(d, tf))
    return 'ICE front moved upstream by %s cells in rows 1-%d' % (note, cf.INNER)


def doctor_T3(run, d, R):
    tf = files(run)['ice']
    t = Tec(tf)
    f, T = t.var('f_p1'), t.var('T_p1')
    m = (f > 0) & (f < 1)
    T[m] *= 1 + 1e-6
    t.write(target(d, tf))
    return 'ICE T_p1 x (1 + 1e-6) in %d plateau cells' % int(m.sum())


def doctor_T3_parcels(run, d, R):
    tf = files(run)['traj']
    Tm = R['mat']['Tm']
    rows = cf.load_rows(tf, 'T', 'f')
    inner = sorted({pid for pid, _, j in R['events'] if 0 <= j < cf.INNER})
    pick = next((pid for pid in inner if sum(1 for r in rows[pid] if 0 < r[1] < 1) >= 4), None)
    if pick is None:
        return None
    ks = [k for k, r in enumerate(rows[pick]) if 0 < r[1] < 1]
    mid = set(ks[1:-1])

    def edit(cols, t, k):
        if int(t[cols['ID']]) != pick or k not in mid:
            return None
        t = list(t)
        t[cols['T']] = '%.6f' % (Tm * (1 + 1e-6))
        return t
    edit_points(tf, target(d, tf), edit)
    return 'IGLOO T = T_melt x (1 + 1e-6) on %d plateau rows of parcel %d' % (len(mid), pick)


def doctor_T4(run, d, R, mesh):
    tf = files(run)['traj']
    Tm = R['mat']['Tm']
    rows = cf.load_rows(tf, 'X', 'Y', 'T', 'f')
    conv = {}
    for pid, i_end, j in R['ends']:
        if not 0 <= j < cf.INNER or R['ice_solid'][j] < 0:
            continue
        a = np.array([r[:4] for r in rows[pid]])
        plat = np.flatnonzero((np.abs(a[:, 2] - Tm) <= cf.HALF_ULP) & (a[:, 3] > 0))
        ii, _ = mesh.locate(a[plat, 0], a[plat, 1])
        cut = R['ice_solid'][j] - (R['n_front'][j] + 2)
        k = plat[(ii >= cut) & (np.arange(len(plat)) >= 3)]
        if len(k) and k[0] > plat[0] + 2:
            conv[pid] = set(k.tolist())
    if not conv:
        return None

    def edit(cols, t, k):
        pid = int(t[cols['ID']])
        if pid in conv and k in conv[pid]:
            t = list(t)
            t[cols['T']] = '%.6f' % (Tm - 1.0)
            t[cols['f']] = '%.6f' % 1.0
            return t
        return None
    edit_points(tf, target(d, tf), edit)
    return 'plateau end of %d inner parcels moved before ICE\'s first f_p1 = 1 minus (N_front + 2) cells' % len(conv)


def exit_records(run):
    ef = files(run)['exit']
    out = []
    for pid, rows in cf.load_rows(ef, 'Y', 'T', 'bc').items():
        for k, r in enumerate(rows):
            out.append((pid, k, r[0], r[1], int(r[2])))
    return out


def y10(run):
    E = cf.read_cc(cf.ice_field(run)[1])
    y = E['y'][0]
    return y[cf.INNER, y.shape[1] - 1]


def doctor_record(run, d, pick, col, value, R):
    """One exit record (pick(records) -> (ID, ordinal)) with column col set to value."""
    ef = files(run)['exit']
    sel = pick(exit_records(run))
    if sel is None:
        return None

    def edit(cols, t, k):
        if int(t[cols['ID']]) == sel[0] and k == sel[1]:
            t = list(t)
            t[cols[col]] = '%.6f' % value
            return t
        return None
    edit_points(ef, target(d, ef), edit)
    return 'exit record of parcel %d: %s = %.6f' % (sel[0], col, value)


def outer_exit(recs, yl):
    c = [(r[0], r[1]) for r in recs if r[4] == 400 and r[2] > yl]
    return c[len(c) // 2] if c else None


def doctor_ice_exit(run, d, R, value):
    tf = files(run)['ice']
    t = Tec(tf)
    rho, f, T = t.var('rho_p1'), t.var('f_p1'), t.var('T_p1')
    ie = rho.shape[1] - 1
    js = [j for j in range(cf.INNER, rho.shape[0]) if rho[j, ie] > cf.MASS * rho.max()]
    if not js:
        return None
    j = js[len(js) // 2]
    f[j, ie] = value
    T[j, ie] = R['mat']['Tm']
    t.write(target(d, tf))
    return 'ICE exit cell (190, %d): f_p1 = %g, T_p1 = T_melt' % (j + 1, value)


def doctor_T6(run, d, R):
    ef = files(run)['exit']
    yl = y10(run)
    diff = R['t6_diff']
    sg = 1.0 if diff >= 0 else -1.0
    m = R['t6_tol'] + abs(diff) + 1.0

    def edit(cols, t, k):
        if int(float(t[cols['bc']])) == 400 and float(t[cols['Y']]) <= yl:
            t = list(t)
            t[cols['T']] = '%.6f' % (float(t[cols['T']]) - sg * m)
            return t
        return None
    edit_points(ef, target(d, ef), edit)
    return 'IGLOO inner exit-plane T %+.3f K' % (-sg * m)


def doctor_C1(run, d, R):
    tf = files(run)['traj']
    Tm, f0 = R['mat']['Tm'], R['f0']
    done = [c for c in R['c1'] if c[1] is not None]
    rmax = max(c[3] for c in done)
    k = (1 + cf.C1_TOL + 0.01) / (1 + rmax)

    def edit(cols, t, kk):
        T, f = float(t[cols['T']]), float(t[cols['f']])
        if abs(T - Tm) <= cf.HALF_ULP and 0 < f < 1:
            t = list(t)
            t[cols['f']] = '%.6f' % (f0 + k * (f - f0))
            return t
        return None
    edit_points(tf, target(d, tf), edit)
    return 'IGLOO plateau f -> f0 + %.4f (f - f0)' % k


def doctor_C2p(run, d, R):
    ef = files(run)['exit']
    yl = y10(run)
    diff = R['c2p_diff']                       # ICE - IGLOO
    sg = -1.0 if diff >= 0 else 1.0            # move IGLOO away from ICE
    m = R['c2p_tol'] + abs(diff) + 0.001

    def edit(cols, t, k):
        if int(float(t[cols['bc']])) == 400 and float(t[cols['Y']]) <= yl:
            t = list(t)
            t[cols['f']] = '%.6f' % (float(t[cols['f']]) + sg * m)
            return t
        return None
    edit_points(ef, target(d, ef), edit)
    return 'IGLOO inner exit-plane f %+.4f' % (sg * m)


# ------------------------------------------------------------------------------------------------ driver
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('run')
    ap.add_argument('scratch')
    ap.add_argument('--design', choices=('full', 'partial'), default='full')
    a = ap.parse_args()
    run, sc, design = os.path.abspath(a.run), os.path.abspath(a.scratch), a.design
    os.makedirs(sc, exist_ok=True)
    R = cf.analyse(run, design)
    mesh = cf.Mesh(cf.read_cc(cf.ice_field(run)[1])['x'][0], cf.read_cc(cf.ice_field(run)[1])['y'][0])
    Tm = R['mat']['Tm']
    yl = y10(run)
    npass = nfail = 0

    def report(ok, text):
        nonlocal npass, nfail
        npass, nfail = npass + ok, nfail + (not ok)
        print('%s  %s' % ('PASS' if ok else 'FAIL', text))

    # the round trip of the writers: an unchanged rewrite reads back identical
    d = make_copy(run, sc, 'roundtrip')
    for key in ('ice', 'euler'):
        f = files(run)[key]
        Tec(f).write(target(d, f))
    for key in ('traj', 'exit'):
        f = files(run)[key]
        edit_points(f, target(d, f), lambda c, t, k: list(t))
    A, B = cf.read_cc(files(run)['ice']), cf.read_cc(files(d)['ice'])
    same = all(np.array_equal(A[n], B[n]) for n in A)
    A, B = cf.read_cc(files(run)['euler']), cf.read_cc(files(d)['euler'])
    same = same and all(np.array_equal(A[n], B[n]) for n in A)
    for key, cols in (('traj', ('X', 'Y', 'T', 'f')), ('exit', ('Y', 'T', 'mdot', 'bc', 'f'))):
        same = same and cf.load_rows(files(run)[key], *cols) == cf.load_rows(files(d)[key], *cols)
    report(same, 'writers: an unchanged rewrite of the four files reads back identical')

    v, out = verdicts(make_copy(run, sc, 'untouched'), design)
    bad = sorted(k for k, s in v.items() if s == 'FAIL')
    report(not bad, 'untouched copy: %s' % (', '.join('%s %s' % kv for kv in sorted(v.items()))))
    if bad:
        print(out)

    full = design == 'full'
    cases = [('T1', 'T1', lambda d: doctor_T1(run, d, R)),
             ('T2', 'T2', lambda d: doctor_T2(run, d, R, mesh)),
             ('T3', 'T3', lambda d: doctor_T3(run, d, R)),
             ('T3-parcels', 'T3', lambda d: doctor_T3_parcels(run, d, R))]
    if full:
        cases += [('T4', 'T4', lambda d: doctor_T4(run, d, R, mesh)),
                  ('T5-record', 'T5', lambda d: doctor_record(run, d, lambda r: outer_exit(r, yl), 'T', Tm + 3.0, R)),
                  ('T5-ice-cell', 'T5', lambda d: doctor_ice_exit(run, d, R, 0.9)),
                  ('T5-wall-record', None, lambda d: doctor_record(
                      run, d, lambda r: next(((x[0], x[1]) for x in r if x[4] == 301), None), 'T', R['T0'] - 1.0, R)),
                  ('T6', 'T6', lambda d: doctor_T6(run, d, R)),
                  ('C1', 'C1', lambda d: doctor_C1(run, d, R)),
                  ('C2', 'C2', lambda d: doctor_record(run, d, lambda r: outer_exit(r, yl), 'f', 0.99, R))]
    else:
        cases += [('T5p-record', 'T5p', lambda d: doctor_record(run, d, lambda r: outer_exit(r, yl), 'T', Tm - 1.0, R)),
                  ('T5p-ice-cell', 'T5p', lambda d: doctor_ice_exit(run, d, R, 1.0)),
                  ('C1', 'C1', lambda d: doctor_C1(run, d, R)),
                  ('C2p', 'C2p', lambda d: doctor_C2p(run, d, R))]
    for name, tgt, doc in cases:
        d = make_copy(run, sc, name)
        what = doc(d)
        if what is None:
            report(False, '%s: no record or cell to doctor' % name)
            continue
        v, out = verdicts(d, design)
        fails = sorted(k for k, s in v.items() if s == 'FAIL')
        ok = fails == ([tgt] if tgt else [])
        report(ok, '%s (%s): FAIL %s, expected %s' % (name, what, fails or 'none', tgt or 'none'))
        if not ok:
            print('\n'.join('      ' + l for l in out.splitlines() if re.match(r'^(PASS|FAIL|SKIP)', l)))
    print('== refute_compare: %d passed, %d failed' % (npass, nfail))
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
