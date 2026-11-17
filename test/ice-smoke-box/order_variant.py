#!/usr/bin/env python3
"""order_variant.py <case-dir> <out-dir> <igfirst|iglast>: test/ice-smoke-box with the burning material A on the IG
closure, in one of the two family orders (after the audit S7's make_order_variants.py):
  igfirst  A = family 1 (IG), S = family 2 (MK)
  iglast   S = family 1 (MK), A = family 2 (IG): ICE steps A last, so the smoke's birth after the family loop must
           convert S with S's own closure, not with the IG one the loop left
The phase file, the zones of the property table, the [ICE-Physics] per-material keys and the families' closures follow
the order; A's IG state is the case's MK state with a granular pressure of 1e-25 (a cold cloud), S empty (1e-25 in each
slot). Reads <case-dir>, writes only <out-dir>."""
import os
import shutil
import sys

EMPTY = 1.0e-25


def main():
    src, out, tag = sys.argv[1].rstrip('/'), sys.argv[2].rstrip('/'), sys.argv[3]
    if tag not in ('igfirst', 'iglast'):
        sys.exit('order_variant.py: igfirst or iglast')
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    shutil.copytree(src + '/INPUT', out + '/INPUT')
    shutil.copytree(src + '/MESH', out + '/MESH')
    shutil.copy2(src + '/input.ini', out + '/input.ini')
    # the initial condition: two families of MK, the case's A then S
    lines = open(src + '/INPUT/part-ic.tec').read().splitlines()
    names = lines[0].split('=')[1].split()
    vals = [float(v) for l in lines[2:] for v in l.split()]
    nnode = 21 * 21 * 2
    ncell = (len(vals) - 3 * nnode) // (len(names) - 3)
    assert 3 * nnode + ncell * (len(names) - 3) == len(vals)
    col = {n: vals[3 * nnode + (i - 3) * ncell:3 * nnode + (i - 2) * ncell] for i, n in enumerate(names) if i >= 3}
    A = [('rp', col['rp_11']), ('up', col['up_11']), ('vp', col['vp_11']), ('wp', col['wp_11']),
         ('Pp', [EMPTY] * ncell), ('Tp', col['Tp_11']), ('np', col['np_11'])]
    S = [(b, [EMPTY] * ncell) for b in ('rp', 'up', 'vp', 'wp', 'Tp', 'np')]
    fams = [('A', A), ('S', S)] if tag == 'igfirst' else [('S', S), ('A', A)]
    vn, data = [], []
    for k, (_, cols) in enumerate(fams, start=1):
        for b, c in cols:
            vn.append('%s_%d1' % (b, k))
            data += list(c)
    z = lines[1].replace('[4-15]=CELLCENTERED', '[4-%d]=CELLCENTERED' % (3 + len(vn)))
    assert z != lines[1]
    with open(out + '/INPUT/part-ic.tec', 'w') as f:
        f.write(' VARIABLES = x y z   ' + '   '.join(vn) + '\n' + z + '\n')
        for v in vals[:3 * nnode] + data:
            f.write(' %.17E\n' % v)
    # the phase file and the property table in the families' order
    pl = open(src + '/INPUT/part-phase.txt').read().splitlines()
    assert pl[1:] == ['A 1 combustion=Beckstead', 'S 1'], pl
    line = {'A': pl[1], 'S': pl[2]}
    open(out + '/INPUT/part-phase.txt', 'w').write('\n'.join([pl[0]] + [line[m] for m, _ in fams]) + '\n')
    t = open(src + '/INPUT/part-properties.dat').read()
    ia, is_ = t.index('ZONE T="A"'), t.index('ZONE T="S"')
    assert ia < is_
    zone = {'A': t[ia:is_], 'S': t[is_:]}
    zone = {k: v if v.endswith('\n') else v + '\n' for k, v in zone.items()}
    open(out + '/INPUT/part-properties.dat', 'w').write(t[:ia] + ''.join(zone[m] for m, _ in fams))
    # the deck: per-material keys and closures in the families' order
    ini = open(src + '/input.ini').read()
    for key in ('density       = ', 'specific-heat = ', 'emissivity    = '):
        i = ini.index('\n' + key) + 1
        j = ini.index('\n', i)
        v = ini[i + len(key):j].split()
        assert len(v) == 2, (key, v)
        if tag == 'iglast':
            ini = ini[:i] + key + v[1] + ' ' + v[0] + ini[j:]
    old = '[ICE-Family1]\nclosure = MK\n\n[ICE-Family2]\nclosure = MK\n'
    assert ini.count(old) == 1
    clo = ('IG', 'MK') if tag == 'igfirst' else ('MK', 'IG')
    ini = ini.replace(old, '[ICE-Family1]\nclosure = %s\n\n[ICE-Family2]\nclosure = %s\n' % clo)
    open(out + '/input.ini', 'w').write(ini)
    print('order_variant.py: %s written (families %s, closures %s)' % (out, ' '.join(m for m, _ in fams), ' '.join(clo)))


if __name__ == '__main__':
    main()
