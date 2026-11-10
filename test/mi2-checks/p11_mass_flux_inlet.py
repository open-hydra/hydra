#!/usr/bin/env python3
"""Numeric checks of p11_mass_flux_inlet_check.sh: what the coupler hands IGLOO at a gas inlet, and what IGLOO
injects from it. Prints PASS/FAIL lines with the numbers compared; exits 1 on any FAIL.
  usage: p11_mass_flux_inlet.py <run-dir> <log> 403 <g>      gas inlet 403 with the prescribed mass flux g [kg/(s m2)]
         p11_mass_flux_inlet.py <run-dir> <log> 401 <g>      gas inlet 401 (T0, p0); g is the 403 twin's flux, a bound only
<run-dir> holds INPUT/{gas-bc.txt,part-bc.txt}, MESH/mesh.tec and OUTPUT/{gas-field.tec,part-euler1-sweep*.tec}.
  report : every ' [MI2-report] igloo-inlet-mdot <sum> <n>' line of the log has n = the inlet records of gas-bc.txt
           (codes 401-407, 420) and, for 403, sum = g * SUM(cell%area) to 1e-8 relative, cell%area recomputed from the
           mesh exactly as IGLOO does (mean of the two i-faces of the first-column hexahedron, each quad split along
           its 1-3 diagonal); for 401, 0.5 <= sum / (g * SUM(area)) <= 2 (MOSE's computed inflow, positive).
  flux   : the inlet-column condensed flux mean_j(rho_p u_p)[i=1] of the last IGLOO euler file equals
           krho/(1-krho) * mean_j(rho u)[i=1] of the gas field within 5 % (krho from part-bc.txt), both above a
           10 kg/(s m2) floor; for 403 it also equals krho/(1-krho) * g within 5 %. The column mean is not weighted
           by the cell areas and follows the gas of the start-up: 0.990 to 1.032 of krho/(1-krho) * g over one to
           four sweeps of 50 to 800 MOSE iterations; a loading applied twice or not at all is 30 % off or more."""
import glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tec_cc import read_cc

INLET_CODES = {401, 402, 403, 404, 405, 406, 407, 420}   # the records MOSExIGLOO hands over
GUIDE = {1: (1, 2, 3, 4), 2: (7, 6, 5, 8)}               # IGLOO's face-vertex table, faces 1 and 2
FLUX_TOL = 0.05                                          # relative band of the inlet-column condensed flux
fails = 0


def verdict(ok, text):
    global fails
    if not ok:
        fails += 1
    print(('PASS  ' if ok else 'FAIL  ') + text)


def header_lines(path):
    """(code list) of the record headers of a bc file: six integers on a line."""
    codes = []
    for line in open(path, encoding='utf-8', errors='replace'):
        t = line.split()
        if len(t) == 6 and all(x.isdigit() for x in t):
            codes.append(int(t[5]))
    return codes


def krho_of(part_bc):
    """krho of the first dispersed inlet record: column 1 of the payload line after a 401 header."""
    lines = open(part_bc, encoding='utf-8', errors='replace').read().splitlines()
    for k, line in enumerate(lines):
        t = line.split()
        if len(t) == 6 and all(x.isdigit() for x in t) and int(t[5]) == 401:
            return float(lines[k + 1].split()[0])
    raise SystemExit('FAIL  no dispersed 401 record in ' + part_bc)


def igloo_face1_areas(mesh):
    """cell%area of the face-1 cells (i = 1, j = 1..Ny, k = 1) as IO.f90 computes it."""
    m = read_cc(mesh)
    x, y, z = m['x'], m['y'], m['z']                       # (K, J, I) nodal
    node = lambda i, j, k: np.array([x[k, j, i], y[k, j, i], z[k, j, i]])
    ny = x.shape[1] - 1

    def quad_area(V, f):                                   # computeArea: triangles (1,2,3) + (1,4,3) of face f
        g = GUIDE[f]; area = 0.0
        for i in (2, 4):
            v1 = V[g[i - 1]] - V[g[0]]; v2 = V[g[2]] - V[g[i - 1]]
            area = np.linalg.norm(np.cross(v1, v2)) + area
        return area / 2

    areas = []
    for j in range(1, ny + 1):                             # hexahedron (i=1, j, k=1): getVertices order
        i, k = 1, 1
        V = {1: node(i - 1, j - 1, k - 1), 2: node(i - 1, j - 1, k), 3: node(i - 1, j, k), 4: node(i - 1, j, k - 1),
             5: node(i, j - 1, k - 1), 6: node(i, j - 1, k), 7: node(i, j, k), 8: node(i, j, k - 1)}
        areas.append(0.5 * (quad_area(V, 1) + quad_area(V, 2)))
    return areas


def report_lines(log):
    out = []
    for line in open(log, encoding='utf-8', errors='replace'):    # the log carries UTF-8 whatever the locale
        t = line.split()
        if len(t) == 4 and t[0] == '[MI2-report]' and t[1] == 'igloo-inlet-mdot':
            out.append((float(t[2]), int(t[3])))
    return out


def column_flux(path, rho_key, u_key):
    d = read_cc(path)
    kr = [k for k in d if rho_key(k)][0]; ku = [k for k in d if u_key(k)][0]
    return float((d[kr][0, :, 0] * d[ku][0, :, 0]).mean())


def main():
    if len(sys.argv) != 5 or sys.argv[3] not in ('401', '403'):
        sys.exit(__doc__)
    run, log, mode, g = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    n_inlet = sum(1 for c in header_lines(run + '/INPUT/gas-bc.txt') if c in INLET_CODES)
    krho = krho_of(run + '/INPUT/part-bc.txt')
    areas = igloo_face1_areas(run + '/MESH/mesh.tec')
    gA = 0.0
    for a in areas:                                        # the coupler's record-order sum
        gA += g * a
    print(f'  inlet records {n_inlet}, krho {krho}, SUM(cell%area) {sum(areas):.10e} m2, g*SUM(area) {gA:.16e} kg/s')

    # --- report lines ---
    rep = report_lines(log)
    if not rep:
        verdict(False, f'{mode}: no " [MI2-report] igloo-inlet-mdot" line in the log (binary without the coupling report)')
    for s, n in rep:
        print(f'  report: sum {s:.16e} kg/s  n {n}  sum/(g*SUM(area)) {s / gA:.6e}')
    verdict(bool(rep) and all(n == n_inlet for _, n in rep),
            f'{mode}: every report line counts n = {n_inlet} inlet records' + ('' if rep else ' (no line)'))
    if mode == '403':
        verdict(bool(rep) and all(abs(s - gA) <= 1e-8 * gA for s, _ in rep),
                f'403: every report sum = g * SUM(cell%area) = {gA:.10e} kg/s to 1e-8 relative' +
                (f' (worst rel. diff {max(abs(s - gA) / gA for s, _ in rep):.3e})' if rep else ' (no line)'))
    else:
        verdict(bool(rep) and all(0.5 * gA <= s <= 2.0 * gA for s, _ in rep),
                f'401: every report sum is MOSE\'s positive computed inflow within [0.5, 2] x g*SUM(area)')

    # --- inlet-column condensed flux ---
    # IGLOO's final write: part-euler1.tec after one sweep, part-euler1-sweep<k>.tec after several
    eul = sorted(glob.glob(run + '/OUTPUT/part-euler1*.tec'),
                 key=lambda f: int(f.rsplit('sweep', 1)[1].split('.')[0]) if 'sweep' in f else 0)
    gas = run + '/OUTPUT/gas-field.tec'
    if not eul or not os.path.exists(gas):
        verdict(False, f'{mode}: inlet-column flux not measurable (missing ' +
                ('IGLOO euler file' if not eul else 'gas-field.tec') + ': the run did not reach its final write)')
        return
    cond = column_flux(eul[-1], lambda k: 'rho' in k.lower(), lambda k: k.lower().startswith('u'))
    gflx = column_flux(gas, lambda k: k == 'rho(1)', lambda k: k == 'u')
    load = krho / (1.0 - krho)
    print(f'  {os.path.basename(eul[-1])}: inlet-column rho_p*u_p {cond:.4f}; gas rho*u {gflx:.4f}; '
          f'krho/(1-krho)*rho*u {load * gflx:.4f}; krho/(1-krho)*g {load * g:.4f} kg/(s m2)')
    verdict(cond > 10.0 and gflx > 10.0 and abs(cond - load * gflx) <= FLUX_TOL * load * gflx,
            f'{mode}: inlet-column condensed flux = krho/(1-krho) x gas rho*u within {100 * FLUX_TOL:.0f} % '
            f'(ratio {cond / (load * gflx):.4f})')
    if mode == '403':
        verdict(abs(cond - load * g) <= FLUX_TOL * load * g,
                f'403: inlet-column condensed flux = krho/(1-krho) x g within {100 * FLUX_TOL:.0f} % '
                f'(ratio {cond / (load * g):.4f})')


if __name__ == '__main__':
    main()
    sys.exit(1 if fails else 0)
