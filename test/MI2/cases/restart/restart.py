#!/usr/bin/env python3
"""MOSE + ICE restart inside hydra-MI2: the checks of check.sh.

  restart.py exact <one-go OUTPUT> <restarted OUTPUT> <part-1 OUTPUT> <MOSE sol-diter> <ICE sol-diter>
                [<MOSE restart iteration> <ICE restart iteration>]   (a killed part 1)
                                                                            PASS/FAIL lines, exit 1 on a FAIL
  restart.py way_out <log> <run dir>    follows, in <run dir>, the way out the last stop of <log> names, item by
                                            item ("keep a copy of X (...)", "copy A to B[ and C to D]", "set K = V in
                                            [S][ and [T]]"); prints each action, exit 1 on an item it cannot follow
  restart.py perturb <tec file>         adds 1e-6 to the first coordinate of the first zone (a mesh other than
                                            the run's)

The run in one go against part 1 + the restart. gas-field.tec, part-field.tec and part-coupling.tec (what ICE's last
coupling step hands the gas and leaves for its next): every header line equal, the output numbers included, every value
of every variable of every zone equal. One quantity is compared as zero instead, where it is below 1e-20 everywhere in
both runs: the component ICE holds at zero on a 2-D mesh, w (part-field `w_p<g>`), what ICE hands of it (the
z-momentum entries `S4`, `B4`, `G4`, `GA4`, `R4` of part-coupling) and its residual columns (the 4th norm of the domain
and of each family in part-residual-history). ICE sets it to zero at every set-up, so a restart cannot carry the
round-off its run leaves there; the variables and columns so compared are named. Every line of gas-diagnostic.tec and
part-residuals.tec equal; a file absent from both runs (a time-accurate solver writes no residual field) is noted, not
passed. Each residual history equal to the one-go file row for row, part 1's rows first (of a killed part 1, its rows up
to the restart iteration: each solver drops the later ones on a restart). The end of each run is its end-of-run
output: in part 1, the run in one go and the restart, no residual history holds an iteration twice in a row, and each
solver's output number is the outputs its iterations made, ITERATION // sol-diter (ITERATION and the number from the
field file). Non-vacuous: part 1's fields differ from the end's, and the restart's counters are those of the whole run.
A missing file is a FAIL.
"""
import glob
import os
import re
import sys

FIELDS = ('gas-field.tec', 'part-field.tec', 'part-coupling.tec')
EXACT = ('gas-diagnostic.tec', 'part-residuals.tec')
HISTORIES = ('gas-residual-history.dat', 'part-residual-history.dat')
COUNTER = re.compile(r'^DATASETAUXDATA (MOSE|ICE)_OUTPUT_NUMBER = "([^"]+)"')
ZERO = 1e-20


class Rep(object):
    def __init__(self):
        self.npass = self.nfail = 0

    def check(self, ok, text):
        print(('PASS  ' if ok else 'FAIL  ') + text)
        if ok:
            self.npass += 1
        else:
            self.nfail += 1


def lines(d, f):
    """The lines of a file; [] when it is missing (every check on it then fails)."""
    if not os.path.isfile(os.path.join(d, f)):
        print('  missing: %s' % os.path.join(d, f))
        return []
    with open(os.path.join(d, f)) as fh:
        return fh.read().split('\n')


def aux(d, f, key):
    """A DATASETAUXDATA value of a field file, or None."""
    for line in lines(d, f):
        m = re.match(r'^DATASETAUXDATA %s = "([^"]+)"' % key, line)
        if m:
            return m.group(1)
    return None


def first_diff(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return i + 1, x.strip()[:60], y.strip()[:60]
    return min(len(a), len(b)) + 1, '(length %d)' % len(a), '(length %d)' % len(b)


def tec(ls):
    """A Tecplot BLOCK file as (header and auxiliary lines, [(name, [value strings of every zone, zone after zone])],
    number of zones), coordinates included; None when the values do not fill the zones exactly."""
    names = [t.strip('"') for t in re.findall(r'"[^"]+"|[^\s,]+', re.search(r'VARIABLES\s*=\s*(.*)', ls[0]).group(1))]
    num = re.compile(r'^\s*[-+]?(\d|\.\d|inf|nan)', re.I)
    head = [x for x in ls if x.strip() and not num.match(x)]
    zones = []
    for x in ls:
        if x.lstrip().startswith('ZONE'):
            I, J, K = (int(g) for g in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', x).groups())
            zones.append((I, J, K, []))
        elif zones and num.match(x):
            zones[-1][3].extend(x.split())
    per = dict((n, []) for n in names)
    for I, J, K, tok in zones:
        nn, nc = I * J * K, (I - 1) * (J - 1) * max(1, K - 1)
        p = 0
        for n in names:
            cnt = nn if n.lower() in ('x', 'y', 'z') else nc
            per[n].extend(tok[p:p + cnt])
            p += cnt
        if p != len(tok):
            return None
    return head, [(n, per[n]) for n in names], len(zones)


def small(v):
    return all(abs(float(x)) < ZERO for x in v)


# ICE's zeroed component on a 2-D mesh and what ICE hands of it: the only quantity compared as zero
HELD = re.compile(r'^(w|S4|B4|G4|GA4|R4)_p\d+$')


def field(rep, one, rst, f):
    a, b = lines(one, f), lines(rst, f)
    if len(a) < 2 or len(b) < 2:
        rep.check(False, '%s: missing' % f)
        return
    ta, tb = tec(a), tec(b)
    if ta is None or tb is None:
        rep.check(False, '%s: the values do not fill its zones (%s)' % (f, 'one go' if ta is None else 'restart'))
        return
    (ha, va, za), (hb, vb, zb) = ta, tb
    ca = [(i, COUNTER.match(x)) for i, x in enumerate(ha) if COUNTER.match(x)]
    cb = [(i, COUNTER.match(x)) for i, x in enumerate(hb) if COUNTER.match(x)]
    step = float(cb[0][1].group(2)) - float(ca[0][1].group(2)) if ca and cb else None
    head = len(ha) == len(hb) and all(x == y for x, y in zip(ha, hb) if not COUNTER.match(x))
    ok_c = (len(ca) == 1 and len(cb) == 1 and ca[0][0] == cb[0][0] and step == 0.0) or (not ca and not cb)
    held, bad = [], []
    if [n for n, _ in va] != [n for n, _ in vb] or za != zb:
        bad.append('the variables or zones differ')
    for (n, x), (_, y) in zip(va, vb):
        if x == y:
            continue
        if HELD.match(n) and len(x) == len(y) and small(x) and small(y):
            held.append(n)
            continue
        rng = [float(v) for v in x]
        scale = max(max(rng) - min(rng), 1e-12 * max(abs(v) for v in rng), 1e-300)
        d = max(abs(float(u) - float(v)) for u, v in zip(x, y)) / scale if len(x) == len(y) else float('inf')
        bad.append('%s by up to %.2e of its range' % (n, d))
    num = ('header equal, the output number %s against %s%s' % (cb[0][1].group(2), ca[0][1].group(2),
                                                                 '' if step == 0.0 else ' (off by %s)' % step)
           if ca and cb else 'header equal' if not ca and not cb else 'an output number in one file only')
    rep.check(head and ok_c and not bad, '%s: %d variables, %d zone(s), %s%s; %s%s' % (
        f, len(va), za, 'every value equal' if not bad else 'values differ',
        ' (ICE\'s zeroed component compared as zero, below %g everywhere in both runs: %s)' % (ZERO, ', '.join(held))
        if held else '', num if head else 'header differs', '' if not bad else ' -- ' + '; '.join(bad[:4])))


def rows(ls):
    return [x.split() for x in ls if x.strip()]


def held_columns(one, ncol):
    """The 1-based columns of part-residual-history.dat that hold ICE's w norm: the 4th norm of the domain and, with
    several families, of each family (iteration, time, then nres norms per block of columns)."""
    names = [n for n, _ in (tec(lines(one, 'part-field.tec')) or ((), [], 0))[1]]
    fam = len(set(int(m.group(1)) for m in (re.match(r'^rho_p(\d+)$', n) for n in names) if m))
    if fam < 1:
        return []
    nres = (ncol - 2) // (1 + fam) if fam > 1 else ncol - 2
    return [6 + nres * k for k in range(1 + fam if fam > 1 else 1)]


def history(rep, one, rst, part1, f, cut=None):
    a, b, p1 = rows(lines(one, f)), rows(lines(rst, f)), rows(lines(part1, f))
    if cut is not None:
        p1 = [r for r in p1 if int(r[0]) <= cut]          # a killed part 1: its rows up to the restart iteration
    zero = []
    if f.startswith('part') and a:
        zero = [c - 1 for c in held_columns(one, len(a[0]))
                if all(abs(float(r[c - 1])) < ZERO for r in a + b if c - 1 < len(r))]
    def same(r, s):
        return len(r) == len(s) and all(x == y or (c in zero and abs(float(x)) < ZERO and abs(float(y)) < ZERO)
                                        for c, (x, y) in enumerate(zip(r, s)))
    used = [c for c in zero if any(len(r) > c and len(s) > c and r[c] != s[c] for r, s in zip(a, b))]
    ok = (len(a) == len(b) and len(a) > 1 and all(same(r, s) for r, s in zip(a, b))
          and len(b) >= len(p1) and all(same(r, s) for r, s in zip(b, p1)))
    rep.check(ok, '%s: %d rows against the run in one go\'s %d, %s%s, part 1\'s %d rows first' % (
        f, len(b), len(a), 'row for row equal' if ok else 'rows differ',
        ', column(s) %s (ICE\'s w norm) compared as zero (below %g in both runs)'
        % (','.join(str(c + 1) for c in used), ZERO) if used else '', len(p1)))


def once(rep, runs, f):
    """The end of a run writes its last row once: no history holds a row twice in a row."""
    out, n = [], 0
    for name, d in runs:
        r = rows(lines(d, f))
        dup = [x[0] for x, y in zip(r, r[1:]) if x[0] == y[0]]
        n += len(dup) + (0 if len(r) > 1 else 1)
        out.append('%s %d of %d%s%s' % (name, len(dup), len(r), ' (iteration %s)' % ', '.join(dup[:3]) if dup else '',
                                        '' if len(r) > 1 else ' (fewer than 2 rows)'))
    rep.check(n == 0, '%s: no iteration written twice in a row (repeated: %s)' % (f, '; '.join(out)))


def counters(rep, runs, sd):
    """Each solver's output number is the outputs its iterations made: ITERATION // sol-diter."""
    for f, key, n in (('gas-field.tec', 'MOSE', sd[0]), ('part-field.tec', 'ICE', sd[1])):
        out, ok = [], True
        for name, d in runs:
            it, num = aux(d, f, key + '_ITERATION'), aux(d, f, key + '_OUTPUT_NUMBER')
            good = it is not None and num is not None and int(float(num)) == int(float(it)) // n
            ok = ok and good
            out.append('%s %s at iteration %s' % (name, '-' if num is None else int(float(num)),
                                                   '-' if it is None else int(float(it))))
        rep.check(ok, '%s: the %s output number is the outputs written, iterations // sol-diter %d (%s)'
                  % (f, key, n, '; '.join(out)))


def exact(one, rst, part1, sdm, sdk, cuts=(None, None)):
    rep = Rep()
    for f in FIELDS:
        field(rep, one, rst, f)
    for f in EXACT:
        if not os.path.isfile(os.path.join(one, f)) and not os.path.isfile(os.path.join(rst, f)):
            print('  note: %s absent from both runs (not judged)' % f)
            continue
        a, b = lines(one, f), lines(rst, f)
        rep.check(a == b and len(a) > 1, '%s: %d lines, %s%s' % (
            f, len(a), 'byte for byte equal' if a == b else 'differ',
            '' if a == b else ' -- first difference at line %d: %s | %s' % first_diff(a, b)))
    for f in HISTORIES:
        history(rep, one, rst, part1, f, cuts[0] if f.startswith('gas') else cuts[1])
    runs = (('part 1', part1), ('one go', one), ('restart', rst))
    for f in HISTORIES:
        once(rep, runs, f)
    counters(rep, runs, (int(sdm), int(sdk)))
    # non-vacuous: the restart went on from part 1 and made the whole run's iterations
    moved = all(lines(part1, f) != lines(rst, f) for f in FIELDS[:2])
    cyc = [(aux(one, f, 'CYCLE'), aux(rst, f, 'CYCLE'), aux(part1, f, 'CYCLE')) for f in FIELDS[:2]]
    ok = moved and all(c[0] == c[1] and c[0] is not None and c[2] is not None and int(c[2]) < int(c[0]) for c in cyc)
    rep.check(ok, 'the restart went on from part 1: its fields moved, and the iteration counts of the restarted '
              'files equal the run in one go\'s (MOSE %s, ICE %s; part 1 %s, %s)'
              % (cyc[0][1], cyc[1][1], cyc[0][2], cyc[1][2]))
    return rep


def way_out(log, d):
    """Follow the way out of the last stop in log, literally, in the run directory d."""
    text = [x for x in open(log, errors='replace').read().split('\n') if 'start a new run from them: ' in x]
    if not text:
        print('  way out: no stop with a way out in ' + log)
        return 1
    items = text[-1].split('start a new run from them: ', 1)[1].strip().split('; ')
    for it in items:
        m = re.match(r'keep a copy of (\S+) \(', it)
        if m:
            src = os.path.join(d, m.group(1))
            open(src + '.kept', 'wb').write(open(src, 'rb').read())
            print('  way out: kept ' + m.group(1) + ' as ' + m.group(1) + '.kept')
            continue
        if it.startswith('copy '):
            for pair in it[5:].split(' and '):
                a, b = pair.split(' to ')
                srcs = sorted(glob.glob(os.path.join(d, a)))
                if not srcs:
                    print('  way out: nothing matches ' + a)
                    return 1
                if b.endswith('/'):
                    os.makedirs(os.path.join(d, b), exist_ok=True)
                    for f in srcs:
                        open(os.path.join(d, b, os.path.basename(f)), 'wb').write(open(f, 'rb').read())
                elif len(srcs) == 1:
                    open(os.path.join(d, b), 'wb').write(open(srcs[0], 'rb').read())
                else:
                    print('  way out: ' + a + ' names several files and ' + b + ' is no directory')
                    return 1
                print('  way out: copied %s (%d file(s)) to %s' % (a, len(srcs), b))
            continue
        m = re.match(r'set (\S+) = (.+?) in ((\[[^\]]+\])( and \[[^\]]+\])*)$', it)
        if m:
            key, val = m.group(1), m.group(2)
            ini = os.path.join(d, 'input.ini')
            ls = open(ini).read().split('\n')
            for sec in re.findall(r'\[([^\]]+)\]', m.group(3)):
                h = [i for i, x in enumerate(ls) if x.strip() == '[' + sec + ']']
                if not h:
                    ls += ['', '[' + sec + ']', key + ' = ' + val]
                    continue
                i, done = h[0] + 1, False
                while i < len(ls) and not ls[i].startswith('['):
                    if re.match(re.escape(key) + r'\s*=', ls[i]):
                        ls[i], done = key + ' = ' + val, True
                    i += 1
                if not done:
                    ls.insert(h[0] + 1, key + ' = ' + val)
            open(ini, 'w').write('\n'.join(ls))
            print('  way out: set %s = %s in %s' % (key, val, m.group(3)))
            continue
        print('  way out: an item not followed: ' + it)
        return 1
    return 0


def perturb(f):
    """The first coordinate of the first zone, plus 1e-6: the hand-off of another mesh."""
    ls = open(f).read().split('\n')
    z = next(i for i, x in enumerate(ls) if x.lstrip().startswith('ZONE'))
    i = next(k for k in range(z + 1, len(ls)) if re.match(r'^\s*[-+]?(\d|\.\d)', ls[k]))
    t = ls[i].split()
    t[0] = '%.17E' % (float(t[0]) + 1e-6)
    ls[i] = ' ' + ' '.join(t)
    open(f, 'w').write('\n'.join(ls))
    return 0


def main():
    if len(sys.argv) == 4 and sys.argv[1] == 'way_out':
        return way_out(sys.argv[2], sys.argv[3])
    if len(sys.argv) == 3 and sys.argv[1] == 'perturb':
        return perturb(sys.argv[2])
    if len(sys.argv) in (7, 9) and sys.argv[1] == 'exact':
        rep = exact(*sys.argv[2:7], cuts=(int(sys.argv[7]), int(sys.argv[8])) if len(sys.argv) == 9 else (None, None))
        print('  exact: %d passed, %d failed' % (rep.npass, rep.nfail))
        return 1 if rep.nfail else 0
    sys.exit(__doc__)


if __name__ == '__main__':
    sys.exit(main())
