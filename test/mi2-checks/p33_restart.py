#!/usr/bin/env python3
"""MOSE + ICE restart inside hydra-MI2: the checks of p33_restart_check.sh.

  p33_restart.py exact <one-go OUTPUT> <restarted OUTPUT> <part-1 OUTPUT>   PASS/FAIL lines, exit 1 on a FAIL
  p33_restart.py plain <one-go OUTPUT> <restarted OUTPUT>                   XFAIL/FAIL line, exit 1 on a FAIL

exact: the run in one go against part 1 + the restart, both on the primed schedule (a coupling step with MOSE-iter 0
at the restart point). Every line of gas-field.tec and part-field.tec equal but each solver's output-number record,
equal or one higher after the restart (today one higher: MI2's final postprocess, an intermediate output, writes one
output more than the run made); every line of gas-diagnostic.tec and part-residuals.tec equal; each residual history
equal to the one-go file once consecutive repeated rows are collapsed, with as many rows or one more (today one more:
that same postprocess writes the last row twice, and the restart keeps part 1's copy). Which of the two holds is
printed. Non-vacuous: part 1's fields differ from the end's, and the restart's counters are those of the whole run.
plain: the same without the priming step. The coupler's ICE source is not stored, so the first MOSE block after the
restart runs without it: the fields must differ (XFAIL, the largest difference over each variable's range printed);
equal fields mean the coupler now restores or primes its source, and the gate's priming step can go (FAIL).
"""
import os
import re
import sys

FIELDS = ('gas-field.tec', 'part-field.tec')
EXACT = ('gas-diagnostic.tec', 'part-residuals.tec')
HISTORIES = ('gas-residual-history.dat', 'part-residual-history.dat')
COUNTER = re.compile(r'^DATASETAUXDATA (MOSE|ICE)_OUTPUT_NUMBER = "([^"]+)"')


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


def exact(one, rst, part1):
    rep = Rep()
    for f in FIELDS:
        a, b = lines(one, f), lines(rst, f)
        ca = [(i, COUNTER.match(x)) for i, x in enumerate(a) if COUNTER.match(x)]
        cb = [(i, COUNTER.match(x)) for i, x in enumerate(b) if COUNTER.match(x)]
        same = len(a) == len(b) and all(x == y for x, y in zip(a, b) if not COUNTER.match(x))
        step = float(cb[0][1].group(2)) - float(ca[0][1].group(2)) if ca and cb else None
        ok_c = len(ca) == 1 and len(cb) == 1 and ca[0][0] == cb[0][0] and step in (0.0, 1.0)
        nd = sum(1 for x, y in zip(a, b) if x != y)
        where = '' if same else ', first difference at line %d: %s | %s' % first_diff(
            [x for x in a if not COUNTER.match(x)], [y for y in b if not COUNTER.match(y)])
        how = {0.0: 'equal', 1.0: 'one higher after the restart: the extra output of the final postprocess'}
        rep.check(same and ok_c, '%s: %d lines, every one equal but the output number, %s against %s (%s): %d '
                  'line(s) differ%s' % (f, len(a), cb[0][1].group(2) if cb else '-', ca[0][1].group(2) if ca else '-',
                                        how.get(step, 'off by %s' % step), nd, where))
    for f in EXACT:
        a, b = lines(one, f), lines(rst, f)
        rep.check(a == b, '%s: %d lines, byte for byte equal%s' % (f, len(a), '' if a == b else
                                                                   ' -- first difference at line %d: %s | %s'
                                                                   % first_diff(a, b)))
    for f in HISTORIES:
        a, b, p1 = lines(one, f), lines(rst, f), lines(part1, f)
        a, b, p1 = [x for x in a if x.strip()], [x for x in b if x.strip()], [x for x in p1 if x.strip()]
        ua, ub = collapse(a), collapse(b)
        extra = len(b) - len(a)
        ok = ua == ub and extra in (0, 1) and len(a) > 1 and b[:len(p1)] == p1
        how = {0: 'as many rows', 1: 'one row more: the last row of part 1 written twice by the final postprocess'}
        rep.check(ok, '%s: %d rows against the run in one go\'s %d, equal once repeated rows are collapsed (%d against '
                  '%d), part 1\'s %d rows first (%s)' % (f, len(b), len(a), len(ub), len(ua), len(p1),
                                                       how.get(extra, '%+d rows' % extra)))
    # non-vacuous: the restart went on from part 1 and made the whole run's iterations
    moved = all(lines(part1, f) != lines(rst, f) for f in FIELDS)
    cyc = [(aux(one, f, 'CYCLE'), aux(rst, f, 'CYCLE'), aux(part1, f, 'CYCLE')) for f in FIELDS]
    ok = moved and all(c[0] == c[1] and c[0] is not None and c[2] is not None and int(c[2]) < int(c[0]) for c in cyc)
    rep.check(ok, 'the restart went on from part 1: its fields moved, and the iteration counts of the restarted '
              'files equal the run in one go\'s (MOSE %s, ICE %s; part 1 %s, %s)'
              % (cyc[0][1], cyc[1][1], cyc[0][2], cyc[1][2]))
    return rep


def collapse(rows):
    """The rows with consecutive repeats removed."""
    return [r for i, r in enumerate(rows) if i == 0 or r != rows[i - 1]]


def read_tec(path):
    """Cell variables of a one-zone Tecplot BLOCK file: {name: [values]} (x, y, z dropped)."""
    txt = open(path).read().split('DATASETAUXDATA')[0]
    names = [t.strip('"') for t in re.findall(r'"[^"]+"|[^\s,]+', re.search(r'VARIABLES\s*=\s*(.*)', txt).group(1))]
    z = re.search(r'ZONE.*?I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', txt, re.S)
    I, J, K = (int(g) for g in z.groups())
    body = txt[z.end():]
    vals = [float(v) for v in body[body.index('\n') + 1:].split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(1, K - 1)
    out, p = {}, 0
    for n in names:
        cnt = nn if n.lower() in ('x', 'y', 'z') else nc
        if cnt == nc:
            out[n] = vals[p:p + cnt]
        p += cnt
    return out


def plain(one, rst):
    worst = (0.0, '')
    for f in FIELDS:
        a, b = read_tec(os.path.join(one, f)), read_tec(os.path.join(rst, f))
        for n in a:
            rng = max(a[n]) - min(a[n])
            scale = max(rng, 1e-12 * max(abs(v) for v in a[n]), 1e-300)
            d = max(abs(x - y) for x, y in zip(a[n], b[n])) / scale
            if d > worst[0]:
                worst = (d, '%s %s' % (f.split('-')[0], n))
    if worst[0] > 0.0:
        print('XFAIL plain restart (no priming step): the fields differ from the run in one go by up to %.2e of a '
              'variable\'s range (%s): the coupler\'s ICE source is not stored, the first MOSE block after the restart '
              'runs without it' % worst)
        return 0
    print('FAIL  plain restart (no priming step): the fields equal the run in one go: the coupler now restores its '
          'ICE source; drop the priming step from p33 and assert the plain restart')
    return 1


def main():
    if len(sys.argv) == 5 and sys.argv[1] == 'exact':
        rep = exact(*sys.argv[2:])
        print('  exact: %d passed, %d failed' % (rep.npass, rep.nfail))
        return 1 if rep.nfail else 0
    if len(sys.argv) == 4 and sys.argv[1] == 'plain':
        return plain(*sys.argv[2:])
    sys.exit(__doc__)


if __name__ == '__main__':
    sys.exit(main())
