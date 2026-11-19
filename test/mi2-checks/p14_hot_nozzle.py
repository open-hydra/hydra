#!/usr/bin/env python3
"""The inputs of the hot JPL nozzle cases (JPL-3solver-hot-2micron, JPL-3solver-hot-20micron), check G1 of p14; it runs
alone on a case directory (or a mutated copy) and needs nothing but its INPUT/ and input.ini.
usage: p14_hot_nozzle.py g1 <case-dir> 2um|20um
G1: (a) line 1 of part{E,L}-phase.txt 'condensed-dispersed modeling=eulerian|lagrangian', line 2 identical in both =
   'A 1' + {solidification=on, T-melt=2327, h-fus=1090000, T-nuc=1861.6 (2um) | 2227 (20um), cp-solid=1360}; (b) both
   properties tables Cp 1380 and Density 4004.62 on all 5000 rows, the files identical; (c) the bc files identical, every
   401 record share 0.15/0.70, kV 1, kT 1, radius 1e-6 | 1e-5; (d) input.ini: [ICE-Physics] emissivity 0, density
   4004.62, specific-heat 1380, T0 2500 | 2400 in [ICB-Block1] and [inj], IGLOO-iter 2, [ICE-Parameters] res-threshold 0,
   [ICE-IO] sol-diter 125, sol-overwrite false, no [IGLOO-Properties] emissivity; (e) both products files
   'material A none 0 none 0.0000000000e+00' (no stoichiometry).
One PASS / FAIL line per check, '== G1 n passed, m failed', exit 1 on a failure."""
import os
import re
import sys

DESIGN = {'2um': dict(TN=1861.6, T0=2500.0, RAD=1e-6, W3F=1e-6, wall=False),
          '20um': dict(TN=2227.0, T0=2400.0, RAD=1e-5, W3F=1e-5, wall=True)}


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


def main():
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == 'g1' and a[2] in DESIGN:
        sys.exit(g1(a[1], a[2]))
    sys.exit(__doc__)


if __name__ == '__main__':
    main()
