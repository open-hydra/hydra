#!/usr/bin/env python3
"""MOSE + ICE restart inside hydra-MI2: the checks of p33_restart_check.sh.

  p33_restart.py exact <one-go OUTPUT> <restarted OUTPUT> <part-1 OUTPUT>   PASS/FAIL lines, exit 1 on a FAIL

The run in one go against part 1 + the restart. gas-field.tec, part-field.tec and part-coupling.tec (what ICE's last
coupling step hands the gas and leaves for its next): every header line equal but each solver's output-number record,
equal or one higher after the restart (one higher while MI2's final postprocess, an intermediate output, writes one
output more than the run made), every value of every variable equal. One class of values is compared as zero instead:
a variable below 1e-20 everywhere in both runs (a component ICE holds at zero, w on a 2-D mesh, and what ICE hands of
it: ICE sets it to zero at every set-up, so a restart cannot carry the round-off its run leaves there); the variables
so compared are named. Every line of gas-diagnostic.tec and part-residuals.tec equal, or the file absent from both
runs (a time-accurate solver writes none). Each residual history equal to the one-go file once consecutive repeated rows are collapsed, with as
many rows or one more (one more while that same postprocess writes the last row twice, and the restart keeps part 1's
copy), a column below 1e-20 throughout both runs compared as zero (named). Non-vacuous: part 1's fields differ from the
end's, and the restart's counters are those of the whole run. A missing file is a FAIL.
"""
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
    """A one-zone Tecplot BLOCK file as (header and auxiliary lines, [(name, [value strings])]), coordinates included."""
    names = [t.strip('"') for t in re.findall(r'"[^"]+"|[^\s,]+', re.search(r'VARIABLES\s*=\s*(.*)', ls[0]).group(1))]
    z = next(i for i, x in enumerate(ls) if x.lstrip().startswith('ZONE'))
    I, J, K = (int(g) for g in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', ls[z]).groups())
    num = re.compile(r'^\s*[-+]?(\d|\.\d|inf|nan)', re.I)
    head = [x for x in ls if x.strip() and not num.match(x)]
    vals = ' '.join(x for x in ls[z + 1:] if num.match(x)).split()
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(1, K - 1)
    out, p = [], 0
    for n in names:
        cnt = nn if n.lower() in ('x', 'y', 'z') else nc
        out.append((n, vals[p:p + cnt]))
        p += cnt
    return head, out


def small(v):
    return all(abs(float(x)) < ZERO for x in v)


def field(rep, one, rst, f):
    a, b = lines(one, f), lines(rst, f)
    if len(a) < 2 or len(b) < 2:
        rep.check(False, '%s: missing' % f)
        return
    ha, va = tec(a)
    hb, vb = tec(b)
    ca = [(i, COUNTER.match(x)) for i, x in enumerate(ha) if COUNTER.match(x)]
    cb = [(i, COUNTER.match(x)) for i, x in enumerate(hb) if COUNTER.match(x)]
    step = float(cb[0][1].group(2)) - float(ca[0][1].group(2)) if ca and cb else None
    head = len(ha) == len(hb) and all(x == y for x, y in zip(ha, hb) if not COUNTER.match(x))
    ok_c = (len(ca) == 1 and len(cb) == 1 and ca[0][0] == cb[0][0] and step in (0.0, 1.0)) or (not ca and not cb)
    held, bad = [], []
    if [n for n, _ in va] != [n for n, _ in vb]:
        bad.append('the variables differ')
    for (n, x), (_, y) in zip(va, vb):
        if x == y:
            continue
        if len(x) == len(y) and small(x) and small(y):
            held.append(n)
            continue
        rng = [float(v) for v in x]
        scale = max(max(rng) - min(rng), 1e-12 * max(abs(v) for v in rng), 1e-300)
        d = max(abs(float(u) - float(v)) for u, v in zip(x, y)) / scale if len(x) == len(y) else float('inf')
        bad.append('%s by up to %.2e of its range' % (n, d))
    how = {0.0: 'equal', 1.0: 'one higher after the restart: the extra output of the final postprocess'}
    num = ('header equal but the output number, %s against %s (%s)' % (cb[0][1].group(2), ca[0][1].group(2),
                                                                        how.get(step, 'off by %s' % step))
           if ca and cb else 'header equal' if not ca and not cb else 'an output number in one file only')
    rep.check(head and ok_c and not bad, '%s: %d variables, %s%s; %s%s' % (
        f, len(va), 'every value equal' if not bad else 'values differ',
        ' (compared as zero, below %g everywhere in both runs: %s)' % (ZERO, ', '.join(held)) if held else '',
        num if head else 'header differs', '' if not bad else ' -- ' + '; '.join(bad[:4])))


def rows(ls):
    return [x.split() for x in ls if x.strip()]


def history(rep, one, rst, part1, f):
    a, b, p1 = rows(lines(one, f)), rows(lines(rst, f)), rows(lines(part1, f))
    ua, ub = collapse(a), collapse(b)
    zero = [c for c in range(len(a[0]) if a else 0)
            if all(abs(float(r[c])) < ZERO for r in a + b if c < len(r))]
    def same(r, s):
        return len(r) == len(s) and all(x == y or (c in zero and abs(float(x)) < ZERO and abs(float(y)) < ZERO)
                                        for c, (x, y) in enumerate(zip(r, s)))
    used = [c for c in zero if any(len(r) > c and len(s) > c and r[c] != s[c] for r, s in zip(ua, ub))]
    extra = len(b) - len(a)
    ok = (len(ua) == len(ub) and all(same(r, s) for r, s in zip(ua, ub)) and extra in (0, 1) and len(a) > 1
          and len(b) >= len(p1) and all(same(r, s) for r, s in zip(b, p1)))
    how = {0: 'as many rows', 1: 'one row more: the last row of part 1 written twice by the final postprocess'}
    rep.check(ok, '%s: %d rows against the run in one go\'s %d, equal once repeated rows are collapsed (%d against %d)%s, '
              'part 1\'s %d rows first (%s)' % (f, len(b), len(a), len(ub), len(ua),
                                               ', column(s) %s compared as zero (below %g in both runs)'
                                               % (','.join(str(c + 1) for c in used), ZERO) if used else '', len(p1), how.get(extra, '%+d rows' % extra)))


def exact(one, rst, part1):
    rep = Rep()
    for f in FIELDS:
        field(rep, one, rst, f)
    for f in EXACT:
        if not os.path.isfile(os.path.join(one, f)) and not os.path.isfile(os.path.join(rst, f)):
            rep.check(True, '%s: absent from both runs' % f)
            continue
        a, b = lines(one, f), lines(rst, f)
        rep.check(a == b and len(a) > 1, '%s: %d lines, byte for byte equal%s' % (
            f, len(a), '' if a == b else ' -- first difference at line %d: %s | %s' % first_diff(a, b)))
    for f in HISTORIES:
        history(rep, one, rst, part1, f)
    # non-vacuous: the restart went on from part 1 and made the whole run's iterations
    moved = all(lines(part1, f) != lines(rst, f) for f in FIELDS[:2])
    cyc = [(aux(one, f, 'CYCLE'), aux(rst, f, 'CYCLE'), aux(part1, f, 'CYCLE')) for f in FIELDS[:2]]
    ok = moved and all(c[0] == c[1] and c[0] is not None and c[2] is not None and int(c[2]) < int(c[0]) for c in cyc)
    rep.check(ok, 'the restart went on from part 1: its fields moved, and the iteration counts of the restarted '
              'files equal the run in one go\'s (MOSE %s, ICE %s; part 1 %s, %s)'
              % (cyc[0][1], cyc[1][1], cyc[0][2], cyc[1][2]))
    return rep


def collapse(rows_):
    """The rows with consecutive repeats removed."""
    return [r for i, r in enumerate(rows_) if i == 0 or r != rows_[i - 1]]


def main():
    if len(sys.argv) == 5 and sys.argv[1] == 'exact':
        rep = exact(*sys.argv[2:])
        print('  exact: %d passed, %d failed' % (rep.npass, rep.nfail))
        return 1 if rep.nfail else 0
    sys.exit(__doc__)


if __name__ == '__main__':
    sys.exit(main())
