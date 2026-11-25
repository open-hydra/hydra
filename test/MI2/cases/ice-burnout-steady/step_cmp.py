#!/usr/bin/env python3
"""step_cmp.py <run A> <run B> <tol> [file ...]: two runs of one steady case compared field by field
(the second with a local time step changed): for each OUTPUT Tecplot ASCII
field (default: every *-field.tec of A's OUTPUT/) and each cell-centred variable, the largest difference over the cells
divided by the largest magnitude of the variable in A (its peak; a velocity component's peak is the peak speed of
its phase, so a component that is round-off, the swirl of a wedge, is judged against the flow); prints one line per
file with the worst variable, and a last line `worst <value> <file> <variable>`; exit 1 when the worst exceeds tol.
A family's variables (<name>_p<k>) other than its density are compared, and their peaks taken, only in the cells the
family occupies (rho_p<k> above ICE's rho_empty, 1e-6 kg/m3, in either run): an empty cell's temperature and
velocity are ratios of floors (the cells beyond a steady burn-out front). Reads only."""
import glob
import os
import re
import sys

RHO_EMPTY = 1e-6


def read(path):
    """Cell-centred variables of a one-zone Tecplot ASCII BLOCK file (the first three nodal)."""
    lines = open(path, errors='replace').read().split('\n')
    head = ' '.join(lines[:2])
    names = re.findall(r'"([^"]+)"', lines[0]) or lines[0].split('=', 1)[1].split()
    m = re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', head)
    I, J, K = (int(x) for x in m.groups())
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(K - 1, 1)
    vals = []
    for l in lines[2:]:
        if re.match(r'\s*[A-Za-z]', l) and not re.match(r'\s*[-+]?[0-9.]', l):
            if l.strip().startswith(('ZONE', 'DATASETAUXDATA', 'AUXDATA')):
                break
            continue
        vals += l.split()
    out, off = {}, 0
    for i, n in enumerate(names):
        k = nn if i < 3 else nc
        if i >= 3:
            out[n] = [float(v.replace('D', 'E')) for v in vals[off:off + k]]
        off += k
    return out


def main():
    a, b, tol = sys.argv[1].rstrip('/'), sys.argv[2].rstrip('/'), float(sys.argv[3])
    files = sys.argv[4:] or sorted(os.path.basename(f) for f in glob.glob(a + '/OUTPUT/*-field.tec'))
    worst = (0.0, '', '')
    for f in files:
        A, B = read(a + '/OUTPUT/' + f), read(b + '/OUTPUT/' + f)
        w = (0.0, '')
        # the cells each family occupies, in either run: its variables other than its density are judged there only
        occ = {}
        for n in A:
            m = re.match(r'rho(_p\d+)$', n)
            if m:
                occ[m.group(1)] = [x > RHO_EMPTY or y > RHO_EMPTY for x, y in zip(A[n], B[n])]

        def cells(n):
            m = re.match(r'.+?(_p\d+)$', n)
            if m and not n.startswith('rho_') and m.group(1) in occ:
                return [c for c, o in enumerate(occ[m.group(1)]) if o]
            return range(len(A[n]))
        speed = {}
        for n in A:
            m = re.match(r'([uvw])(_p\d+)?$', n)
            if m:
                k = m.group(2) or ''
                speed[k] = max([speed.get(k, 0.0)] + [abs(A[n][c]) for c in cells(n)])
        for n in A:
            cs = cells(n)
            if not cs:
                continue
            peak = max(abs(A[n][c]) for c in cs)
            m = re.match(r'([uvw])(_p\d+)?$', n)
            if m:
                peak = speed[m.group(2) or '']
            if peak == 0.0:
                d = max(abs(B[n][c]) for c in cs)
            else:
                d = max(abs(A[n][c] - B[n][c]) for c in cs) / peak
            if d > w[0]:
                w = (d, n)
        tot = len(next(iter(A.values()))) if A else 0
        held = ', '.join('%s %d of %d cells' % (k[1:], sum(v), tot) for k, v in sorted(occ.items()))
        print('%-28s worst %-10s %.3e%s' % (f, w[1], w[0], ('   (occupied: ' + held + ')') if held else ''))
        if w[0] > worst[0]:
            worst = (w[0], f, w[1])
    print('worst %.6e %s %s' % worst)
    sys.exit(0 if worst[0] <= tol else 1)


if __name__ == '__main__':
    main()
