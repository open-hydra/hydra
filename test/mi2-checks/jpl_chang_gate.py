#!/usr/bin/env python3
"""The JPL nozzle runs against the curves of Chang, Hourng & Chien (1996), Fig. 11: gas Mach number on the axis
and on the wall.

usage: jpl_chang_gate.py [--ref JPL_chang.dat] [--xmin 0.2] [<label>=<run dir or gas-field.tec> ...]
         REPORT MODE: every row is printed, none is judged, the exit status is 0. Without a run: OUTPUT/gas-field.tec
         of the three JPL cases of this tree (IGLOO, ICE, both), those that exist.
       jpl_chang_gate.py --live [--require-probes] IGLOO=<run dir> ICE=<run dir> both=<run dir>
         LIVE MODE (gate p18): the rows R1-R6 of each run judged against the bounds of MEASURED (below); a label names
         the case whose bounds apply (IGLOO = JPL-Lagrangian-20micron, ICE = JPL-Eulerian-20micron, both =
         JPL-3solver-20micron); the rows that compare cases run when the three labels are given. The run's log is
         <run dir>.log. Exit status 1 when a row fails.
       jpl_chang_gate.py --measure IGLOO=<run dir> ICE=<run dir> both=<run dir>
         prints the MEASURED block for reference runs of the three cases (the rows R2 come from the probe files).
       jpl_chang_gate.py --inject <input.ini>
         appends the two exit probes of R2 to a case's input.ini (a scratch copy: p18_jpl_chang_check.sh does it).

Definitions (those of JPL-Lagrangian-20micron/jpl_compare.py; the reader is tec_cc.read_cc):
  X = (x - x_throat)/D*, x of a cell = mean of the four nodes of its k = 0 face; axis = first cell row, wall = last;
  M = |u|/sqrt(g R T). A reference zone is interpolated onto the run's cells; a row uses the cells with
  xmin < X <= the last station of the reference zone. The exit Mach number is taken at that last station.
  Where the reference curve of the runs' case has a Mach drop (its maximum, then a minimum more than DROP_MIN below
  it: the axis), two more lines: the row on the cells outside the stretch between those two stations, and the drop
  itself -- maximum, minimum and the place of half height, of the run and of the reference. The reference has six
  stations over that stretch: a pointwise difference there measures the place of the front more than the curve.
Reference zones: OG = gas alone, MPS = 1 micron radius, MPL = 10 micron radius (the case the runs reproduce).
Rows of the live mode:
  R1 health: '# exit=0', GAME OVER, none of the failure signatures; with IGLOO one 'All particles out of domain' per
     coupling step (N-coupling of the run's input.ini) and no 'marking gone'.
  R2 steady: the two exit probes written (OUTPUT/exit_axis.txt, exit_wall.txt; with --require-probes); the range of
     the axis probe's M over the last MOSE block below an absolute bound: on a converged run the range is the floor,
     which moves by up to a factor 10 from block to block and with the thread count (1e-10 to 2e-8 with IGLOO, 1e-14
     to 1e-12 without), against 1e-6 to 2e-5 on the axis of a run that is not steady. The wall probe's range (no such
     gap: 2e-9 on a half-length run) and the density residual's last-block maximum over its first-block maximum are
     printed, not judged.
  R3 absolute: |run - MPL| rows: axis rms outside the reference's drop, axis exit M, wall rms, wall exit M, the drop's
     place (half height X) and depth (maximum - minimum), each within its bound.
  R4 class: rms against MPL below half the smaller rms against MPS and OG, both rows.
  R5 offset: the mean of run - MPL (axis outside the drop, wall) within a tolerance of the measured mean; the std
     within its bound.
  R6 solvers: the shares k of the inflow (field 1 of each phase's first 401 record) sum to the cases' loading; the
     inlet column (i = 1), cell by cell, condensed flux rho_p u_p over k times the LOCAL gas flux rho u: ICE's 40 cells;
     IGLOO's projection on its 38 inner cells and its two end cells apart; with two phases SUM rho_p of IGLOO's
     projection over ICE's field; across the cases, the three-solver run's mean offsets between those of ICE alone and
     IGLOO alone, and the gas Mach of ICE alone against IGLOO alone."""
import glob
import hashlib
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tec_cc import read_cc  # noqa: E402

X_THROAT, D_STAR = 0.064207, 2.0 * 0.02033740
DROP_MIN = 0.1            # a reference curve has a drop when it falls by more than this after its maximum
ROWS = (('axis', 0, 'Axis'), ('wall', -1, 'Wall'))
DEFAULT = (('IGLOO', 'JPL-Lagrangian-20micron'), ('ICE', 'JPL-Eulerian-20micron'), ('both', 'JPL-3solver-20micron'))
HYGIENE = ('[ERROR]', ' FAIL', 'on-physical state', 'Integration failed at', 'p<0 detected', 'ABORT:')  # p10, p11, p19

# ---- live mode: the rules and the measured values (the only numbers p18 judges with) -------------------------------
R3_FACTOR = 1.2   # R3: bound = 1.2 x the |value| of the reference run (round plan 7.3, R3)
R3_FLOOR = 0.01   # R3: a bound is never below 0.01 in Mach (round plan 7.3, R3)
R4_MARGIN = 2.0   # R4: rms against MPL < 1/2 min(rms against MPS, against OG) (round plan 7.3, R4)
R5_AXIS = 0.010   # R5: plan 7.3, axis outside the drop; smallest RED move 0.0886 (27 % loading), 2x and threads 1.1e-9
R5_WALL = 0.006   # R5: plan 7.3, wall; smallest RED move 0.0218 (27 % loading), 2x run and threads 3.6e-10
R5_STD = 1.2      # R5: std of run - MPL within 1.2 x the measured std (plan 7.3: 'std bounded', no number; R3's factor)
R6_FACTOR = 1.4   # R6: |ratio - 1| within 1.4 x the measured |ratio - 1| (round plan 7.3, R6)
R2_AXIS = 1e-7    # R2: 6 x the converged floor 1.7e-8 (1, 4 threads, 2x runs), 1/12 of the half-length ICE case 1.2e-6
R2_FACTOR = 1.2   # R2, wall: printed beside the range, not judged: the range is the floor (with IGLOO 1e-10 to 2e-8)
K_PRESCRIBED = 0.30 / 0.70  # R6: the cases' loading krho 0.30 (0.15 per phase with two) as share k = krho/(1 - krho)
K_TOL = 1e-12     # R6: relative; the 401 records print k with 17 digits
# MEASURED: the values of the reference runs, printed by --measure; re-taken there when a pin moves a JPL row
# MEASURED-BEGIN
# one-thread runs of the three cases (with the probes) on hydra main 35b36fb, ifx 2023.1.0, 2026-11-13; their
# fields equal byte for byte those of c4df4a0 (the paper gas): bounds from these values by the rules above
MEASURED = {
    'IGLOO': {  # JPL-Lagrangian-20micron
        'axis_out_rms': 0.0418134,  # R3: rms of run - MPL, axis cells outside the reference drop
        'axis_exit': -0.053923,  # R3: axis exit M - MPL at the last reference station
        'wall_rms': 0.0191081,  # R3: rms of run - MPL, wall
        'wall_exit': -0.0188546,  # R3: wall exit M - MPL
        'drop_half': 0.0133912,  # R3: half-height X of the drop - MPL (interpolated; a minimum is on a cell)
        'drop_depth': 0.00926615,  # R3: depth of the drop (max - min) - MPL
        'axis_out_mean': -0.0403269,  # R5: mean of run - MPL, axis outside the drop
        'axis_out_std': 0.0110498,  # R5: its std
        'wall_mean': -0.0167405,  # R5: mean of run - MPL, wall
        'wall_std': 0.00921278,  # R5: its std
        'r2_axis': 4.51484e-09,  # R2, reference only (R2_AXIS bounds it): axis probe M range, last block
        'r2_wall': 6.6024e-10,  # R2 (printed only): range of the wall exit probe M over the last block
        'r6_igloo_lo': 0.986009,  # R6: lowest ratio of IGLOO's 38 inner inlet cells
        'r6_igloo_hi': 1.02792,  # R6: highest ratio of IGLOO's 38 inner inlet cells
        'r6_igloo_axis': 1.05206,  # R6: IGLOO's inlet cell on the axis
        'r6_igloo_wall': 1.05551,  # R6: IGLOO's inlet cell on the wall
    },
    'ICE': {  # JPL-Eulerian-20micron
        'axis_out_rms': 0.0176352,  # R3: rms of run - MPL, axis cells outside the reference drop
        'axis_exit': -0.00909924,  # R3: axis exit M - MPL at the last reference station
        'wall_rms': 0.0131392,  # R3: rms of run - MPL, wall
        'wall_exit': -0.0108249,  # R3: wall exit M - MPL
        'drop_half': 0.0304751,  # R3: half-height X of the drop - MPL (interpolated; a minimum is on a cell)
        'drop_depth': 0.0258507,  # R3: depth of the drop (max - min) - MPL
        'axis_out_mean': -0.013818,  # R5: mean of run - MPL, axis outside the drop
        'axis_out_std': 0.0109574,  # R5: its std
        'wall_mean': -0.00957433,  # R5: mean of run - MPL, wall
        'wall_std': 0.00899832,  # R5: its std
        'r2_axis': 3.8991e-13,  # R2, reference only (R2_AXIS bounds it): axis probe M range, last block
        'r2_wall': 3.9968e-14,  # R2 (printed only): range of the wall exit probe M over the last block
        'r6_ice': 0.000761028,  # R6: max |ratio - 1| of ICE's inlet cells (part, 1.0000..1.0008)
    },
    'both': {  # JPL-3solver-20micron
        'axis_out_rms': 0.0287272,  # R3: rms of run - MPL, axis cells outside the reference drop
        'axis_exit': -0.0310335,  # R3: axis exit M - MPL at the last reference station
        'wall_rms': 0.0162315,  # R3: rms of run - MPL, wall
        'wall_exit': -0.0152066,  # R3: wall exit M - MPL
        'drop_half': 0.0216245,  # R3: half-height X of the drop - MPL (interpolated; a minimum is on a cell)
        'drop_depth': 0.014665,  # R3: depth of the drop (max - min) - MPL
        'axis_out_mean': -0.0271995,  # R5: mean of run - MPL, axis outside the drop
        'axis_out_std': 0.00924327,  # R5: its std
        'wall_mean': -0.0134274,  # R5: mean of run - MPL, wall
        'wall_std': 0.00911957,  # R5: its std
        'r2_axis': 1.23669e-09,  # R2, reference only (R2_AXIS bounds it): axis probe M range, last block
        'r2_wall': 3.3231e-10,  # R2 (printed only): range of the wall exit probe M over the last block
        'r6_ice': 0.000743324,  # R6: max |ratio - 1| of ICE's inlet cells (partE, 1.0000..1.0007)
        'r6_igloo_lo': 0.985568,  # R6: lowest ratio of IGLOO's 38 inner inlet cells
        'r6_igloo_hi': 1.02977,  # R6: highest ratio of IGLOO's 38 inner inlet cells
        'r6_igloo_axis': 1.05206,  # R6: IGLOO's inlet cell on the axis
        'r6_igloo_wall': 1.05761,  # R6: IGLOO's inlet cell on the wall
        'r6_rho': 1.04699,  # R6: SUM rho_p of IGLOO's projection over ICE's field
    },
    'cross': {  # ICE alone against IGLOO alone, gas Mach on X > 0.2
        'ice_igloo_axis': 0.0337456,  # R6: rms of the difference, axis
        'ice_igloo_wall': 0.00719239,  # R6: rms of the difference, wall
    },
}
# MEASURED-END


def chang_zone(ref, name):
    """(X, M) of one zone of the reference file."""
    lines = open(ref).read().splitlines()
    for i, line in enumerate(lines):
        if 'ZONE T="%s"' % name in line:
            n = int(re.search(r'I\s*=\s*(\d+)', lines[i + 1]).group(1))
            vals = []
            for l2 in lines[i + 2:]:
                try:
                    vals += [float(t) for t in l2.split()]
                except ValueError:
                    continue
                if len(vals) >= 2 * n:
                    break
            return np.array(vals[:n]), np.array(vals[n:2 * n])
    raise KeyError('%s: no zone "%s"' % (ref, name))


def mach_row(field, jcell):
    """(X, M) along one cell row of a gas field read by read_cc."""
    x = field['x']
    j = jcell % field['T'].shape[1]
    xc = 0.25 * (x[0, j, :-1] + x[0, j, 1:] + x[0, j + 1, :-1] + x[0, j + 1, 1:])
    u, v, w = (field[k][0, j, :] for k in ('u', 'v', 'w'))
    mach = np.sqrt(u * u + v * v + w * w) / np.sqrt(field['g'][0, j, :] * field['R'][0, j, :] * field['T'][0, j, :])
    return (xc - X_THROAT) / D_STAR, mach


def against(X, M, zone, xmin):
    """Run minus reference zone on the cells the zone covers downstream of xmin."""
    xr, mr = zone
    s = (X > xmin) & (X <= xr[-1])
    d = M[s] - np.interp(X[s], xr, mr)
    return dict(n=int(s.sum()), x0=X[s].min(), x1=X[s].max(), rms=np.sqrt((d ** 2).mean()), mean=d.mean(), std=d.std(),
                exit=float(np.interp(xr[-1], X, M)), refexit=mr[-1])


def drop(X, M, xlo, xhi=None):
    """The Mach drop of a curve: its maximum on xlo < X (<= xhi), the minimum after it, the place of half height.
    None when the maximum is the last point or the fall after it is below DROP_MIN."""
    s = np.where((X > xlo) & (X <= (X[-1] if xhi is None else xhi)))[0]
    i = s[np.argmax(M[s])]
    if i == len(M) - 1:
        return None
    j = i + int(np.argmin(M[i:]))
    if M[i] - M[j] < DROP_MIN:
        return None
    half = 0.5 * (M[i] + M[j])
    k = i + int(np.argmax(M[i:j + 1] <= half))                   # the first point at or below half height
    xh = X[k - 1] + (M[k - 1] - half) / (M[k - 1] - M[k]) * (X[k] - X[k - 1])
    return dict(xmax=X[i], mmax=M[i], xmin=X[j], mmin=M[j], xhalf=xh)


def gas_field_path(path):
    if os.path.isdir(path):
        return os.path.join(path, 'OUTPUT', 'gas-field.tec') if os.path.isdir(os.path.join(path, 'OUTPUT')) \
            else os.path.join(path, 'gas-field.tec')
    return path


def report(runs, ref, xmin):
    print('== JPL nozzle against Chang, Hourng & Chien (1996), Fig. 11 -- report mode: nothing is judged')
    print('reference %s (md5 %s), rows on %g < X <= the last reference station'
          % (ref, hashlib.md5(open(ref, 'rb').read()).hexdigest()[:8], xmin))
    zones = {(kind, Z): chang_zone(ref, 'Chang %s %s' % (kind, Z)) for kind in ('MPL', 'MPS', 'OG') for _, _, Z in ROWS}
    rows = {}
    for lab, path in runs:
        path = gas_field_path(path)
        field = read_cc(path)
        print('%s: %s' % (lab, path))
        for row, j, Z in ROWS:
            X, M = mach_row(field, j)
            rows[lab, row] = (X, M)
            a, b, o = (against(X, M, zones[kind, Z], xmin) for kind in ('MPL', 'MPS', 'OG'))
            print('  %s: %d cells, X %.3f..%.3f | run - MPL: rms %.4f mean %+.4f std %.4f | exit M %.4f, MPL %.3f (%+.4f)'
                  ' | rms against MPS %.4f, OG %.4f | margin %.2f'
                  % (row, a['n'], a['x0'], a['x1'], a['rms'], a['mean'], a['std'], a['exit'], a['refexit'],
                     a['exit'] - a['refexit'], b['rms'], o['rms'], min(b['rms'], o['rms']) / a['rms']))
            xr, mr = zones['MPL', Z]
            dr = drop(xr, mr, xmin)
            if dr is None:
                continue
            out = (X > xmin) & (X <= xr[-1]) & ((X <= dr['xmax']) | (X >= dr['xmin']))
            d = M[out] - np.interp(X[out], xr, mr)
            print('  %s outside the drop (X <= %.3f or X >= %.3f): %d cells | run - MPL: rms %.4f mean %+.4f std %.4f'
                  % (row, dr['xmax'], dr['xmin'], out.sum(), np.sqrt((d ** 2).mean()), d.mean(), d.std()))
            du = drop(X, M, xmin, dr['xmin'])
            if du is None:
                print('  %s drop: none in the run | MPL: max %.3f at X %.3f, min %.3f at X %.3f, half height at X %.3f'
                      % (row, dr['mmax'], dr['xmax'], dr['mmin'], dr['xmin'], dr['xhalf']))
                continue
            print('  %s drop: max %.4f at X %.3f, min %.4f at X %.3f, half height at X %.3f | MPL: max %.3f at X %.3f, '
                  'min %.3f at X %.3f, half height at X %.3f | run - MPL: max %+.4f, min %+.4f, depth %+.4f, '
                  'half height X %+.3f'
                  % (row, du['mmax'], du['xmax'], du['mmin'], du['xmin'], du['xhalf'], dr['mmax'], dr['xmax'],
                     dr['mmin'], dr['xmin'], dr['xhalf'], du['mmax'] - dr['mmax'], du['mmin'] - dr['mmin'],
                     (du['mmax'] - du['mmin']) - (dr['mmax'] - dr['mmin']), du['xhalf'] - dr['xhalf']))
    labs = [lab for lab, _ in runs]
    for i, A in enumerate(labs):
        for B in labs[i + 1:]:
            out = []
            for row, _, _ in ROWS:
                (XA, MA), (XB, MB) = rows[A, row], rows[B, row]
                if XA.shape != XB.shape or not np.allclose(XA, XB, rtol=0.0, atol=1e-12):
                    out.append('%s: different cells' % row)
                    continue
                d = (MA - MB)[XA > xmin]
                out.append('%s rms %.4f max %.4f' % (row, np.sqrt((d ** 2).mean()), np.abs(d).max()))
            print('gas Mach, %s - %s (X > %g): %s' % (A, B, xmin, ' | '.join(out)))
    return 0


# ---- the quantities of one run ------------------------------------------------------------------------------------
def ini_int(run, key):
    m = re.search(r'^\s*%s\s*=\s*(\d+)' % re.escape(key), open(os.path.join(run, 'input.ini')).read(), re.M)
    return int(m.group(1)) if m else None


def chang_rows(run, zones, xmin):
    """The R3, R4, R5 quantities of a run, and its two Mach rows (for the comparison across cases)."""
    field = read_cc(gas_field_path(run))
    q = {}
    for row, j, Z in ROWS:
        X, M = mach_row(field, j)
        a, b, o = (against(X, M, zones[kind, Z], xmin) for kind in ('MPL', 'MPS', 'OG'))
        q[row] = dict(X=X, M=M, rms=a['rms'], mean=a['mean'], std=a['std'], exit=a['exit'] - a['refexit'],
                      margin=min(b['rms'], o['rms']) / a['rms'])
        xr, mr = zones['MPL', Z]
        dr = drop(xr, mr, xmin)
        if dr is None:
            continue
        out = (X > xmin) & (X <= xr[-1]) & ((X <= dr['xmax']) | (X >= dr['xmin']))
        d = M[out] - np.interp(X[out], xr, mr)
        q[row].update(out_rms=np.sqrt((d ** 2).mean()), out_mean=d.mean(), out_std=d.std())
        du = drop(X, M, xmin, dr['xmin'])
        q[row].update(half=None if du is None else du['xhalf'] - dr['xhalf'],
                      depth=None if du is None else (du['mmax'] - du['mmin']) - (dr['mmax'] - dr['mmin']),
                      xmin_run=None if du is None else du['xmin'], xmin_ref=dr['xmin'])
    return q


def inlet_k(bc):
    """Field 1 of the first 401 record of a bc file: the phase's share k of the inflow."""
    lines = open(bc).read().splitlines()
    for i, line in enumerate(lines):
        t = line.split()
        if len(t) == 6 and all(x.isdigit() for x in t) and t[5] == '401':
            return float(lines[i + 1].split()[0].rstrip(',').replace('D', 'E'))
    raise ValueError('%s: no 401 record' % bc)


def column(d, pick):
    k = [n for n in d if pick(n)]
    if len(k) != 1:
        raise ValueError('ambiguous variable: %s' % k)
    return np.asarray(d[k[0]])


def inlet(run):
    """R6 quantities: per condensed phase the ratio of its inlet-column flux to k x the local gas flux; SUM rho_p."""
    g = read_cc(gas_field_path(run))
    gf = np.asarray(g['rho(1)'])[0, :, 0] * np.asarray(g['u'])[0, :, 0]
    out = {}
    for f in glob.glob(os.path.join(run, 'OUTPUT', '*-field.tec')):
        ph = os.path.basename(f)[:-len('-field.tec')]
        if ph == 'gas':
            continue
        c = read_cc(f)
        rho = column(c, lambda n: n.lower().startswith('rho_p'))
        u = column(c, lambda n: n.lower().startswith('u_p'))
        k = inlet_k(os.path.join(run, 'INPUT', ph + '-bc.txt'))
        out['ice'] = dict(r=rho[0, :, 0] * u[0, :, 0] / (k * gf), rho=rho.sum(), k=k, phase=ph)
    eul = glob.glob(os.path.join(run, 'OUTPUT', '*-euler1*.tec'))
    if eul:
        last = sorted(eul, key=lambda p: int(p.rsplit('sweep', 1)[1].split('.')[0]) if 'sweep' in p else 0)[-1]
        ph = os.path.basename(last).split('-euler1')[0]
        c = read_cc(last)
        rho = column(c, lambda n: 'rho' in n.lower())
        u = column(c, lambda n: n.lower().startswith('u'))
        k = inlet_k(os.path.join(run, 'INPUT', ph + '-bc.txt'))
        out['igloo'] = dict(r=rho[0, :, 0] * u[0, :, 0] / (k * gf), rho=rho.sum(), k=k, phase=ph,
                            file=os.path.basename(last))
    return out


PROBE_CELLS = (('exit_axis', '1 190 1 1'), ('exit_wall', '1 190 40 1'))  # block i j k: the last cell column, j 1 and 40
PROBE_PER_BLOCK = 20   # R2: probe rows per MOSE block (round plan 7.3, R2: diter = block/20)


def inject(ini):
    """Append the R2 probes to an input.ini: M at the exit cell of the axis and of the wall row, every block/20."""
    txt = open(ini).read()
    if re.search(r'^\s*\[MOSE-Probes\]', txt, re.M):
        sys.exit('%s: has a [MOSE-Probes] section already' % ini)
    blk = int(re.search(r'^\s*MOSE-iter\s*=\s*(\d+)', txt, re.M).group(1))
    if blk % PROBE_PER_BLOCK:
        sys.exit('%s: MOSE-iter %d is not a multiple of %d' % (ini, blk, PROBE_PER_BLOCK))
    out = '\n[MOSE-Probes]\n' + ''.join('probe%d = %s\n' % (n + 1, p) for n, (p, _) in enumerate(PROBE_CELLS))
    for p, cell in PROBE_CELLS:
        out += '\n[%s]\nvariables = M\nindex-position = %s\nditer = %d\n' % (p, cell, blk // PROBE_PER_BLOCK)
    open(ini, 'a').write(out)
    return 0


def probes(run):
    """R2 quantities: last-block range of M at the exit probes (None without probe files), the residual ratio."""
    blk = ini_int(run, 'MOSE-iter')
    out = {}
    for p in ('axis', 'wall'):
        f = os.path.join(run, 'OUTPUT', 'exit_%s.txt' % p)
        if not os.path.isfile(f):
            out[p] = None
            continue
        rows = {}
        for line in open(f):
            t = line.split()
            if len(t) >= 2:
                rows[int(float(t[0]))] = float(t[1])          # the last row is written twice: keyed on the iteration
        it = np.array(sorted(rows))
        m = np.array([rows[i] for i in it])
        last = m[it > it[-1] - blk]
        out[p] = float(last.max() - last.min())
    h = np.loadtxt(os.path.join(run, 'OUTPUT', 'gas-residual-history.dat'))
    out['res_ratio'] = float(h[-blk:, 2].max() / h[:blk, 2].max())
    out['res_last'] = float(h[-1, 2])
    return out


def health(run, label):
    """R1: (ok, text) from the run's log."""
    log = run.rstrip('/') + '.log'
    if not os.path.isfile(log):
        return False, 'no log %s' % log
    txt = open(log, errors='replace').read()
    m = re.findall(r'# exit=(\d+)', txt)
    bad = [t for t in HYGIENE if t in txt]
    ok = bool(m) and m[-1] == '0' and 'GAME OVER' in txt and not bad
    text = 'exit %s, GAME OVER %d, failure signatures %s' % (m[-1] if m else '?', txt.count('GAME OVER'), bad or 'none')
    if label in ('IGLOO', 'both'):
        n, out, gone = ini_int(run, 'N-coupling'), txt.count('All particles out of domain'), txt.count('marking gone')
        ok = ok and out == n and gone == 0
        text += ", 'All particles out of domain' %d (N-coupling %s), 'marking gone' %d" % (out, n, gone)
    return ok, text


# ---- live mode ----------------------------------------------------------------------------------------------------
class Verdict:
    def __init__(self):
        self.passed = self.failed = 0

    def __call__(self, ok, text):
        if ok:
            self.passed += 1
        else:
            self.failed += 1
        print(('PASS  ' if ok else 'FAIL  ') + text)


def r3_bound(v):
    return max(R3_FACTOR * abs(v), R3_FLOOR)


def live(runs, ref, xmin, require_probes):
    labels = dict(runs)
    if not set(labels) <= {'IGLOO', 'ICE', 'both'} or not MEASURED:
        sys.exit('--live takes IGLOO=, ICE=, both= and a filled MEASURED block')
    print('== p18: the JPL nozzle against Chang, Hourng & Chien (1996), Fig. 11 -- live rows R1-R6')
    print('reference %s (md5 %s)' % (ref, hashlib.md5(open(ref, 'rb').read()).hexdigest()[:8]))
    zones = {(kind, Z): chang_zone(ref, 'Chang %s %s' % (kind, Z)) for kind in ('MPL', 'MPS', 'OG') for _, _, Z in ROWS}
    v, q = Verdict(), {}
    for lab, run in runs:
        m = MEASURED[lab]
        print('-- %s: %s' % (lab, run))
        ok, text = health(run, lab)
        v(ok, '%s R1 health: %s' % (lab, text))
        p = probes(run)
        for row in ('axis', 'wall'):
            if p[row] is None:
                if require_probes:
                    v(False, '%s R2 %s exit probe: no probe file' % (lab, row))
                else:
                    print('  %s R2 %s exit probe: no probe file, not judged' % (lab, row))
                continue
            if row == 'axis':
                v(p[row] <= R2_AXIS, '%s R2 axis exit Mach, range over the last block %.3e <= %.0e'
                  % (lab, p[row], R2_AXIS))
                continue
            if require_probes:
                v(True, '%s R2 %s exit probe: written' % (lab, row))
            print('  %s R2 %s exit Mach, range over the last block %.3e (1.2 x the reference run %.3e): reported, not '
                  'judged' % (lab, row, p[row], R2_FACTOR * m['r2_' + row]))
        print('  %s R2 density residual, last-block max over first-block max %.2e (last %.2e): reported, not judged'
              % (lab, p['res_ratio'], p['res_last']))
        c = chang_rows(run, zones, xmin)
        q[lab] = c
        a, w = c['axis'], c['wall']
        for name, val, key in (('axis rms outside the drop', a['out_rms'], 'axis_out_rms'),
                               ('axis exit M - MPL', a['exit'], 'axis_exit'),
                               ('wall rms', w['rms'], 'wall_rms'), ('wall exit M - MPL', w['exit'], 'wall_exit'),
                               ('drop: half-height X - MPL', a['half'], 'drop_half'),
                               ('drop: depth - MPL', a['depth'], 'drop_depth')):
            b = r3_bound(m[key])
            if val is None:
                v(False, '%s R3 %s: no drop in the run (bound %.4f)' % (lab, name, b))
                continue
            v(abs(val) <= b, '%s R3 %s: |%+.4f| <= %.4f' % (lab, name, val, b))
        if a['xmin_run'] is not None:
            print('  %s R3 drop minimum at X %.3f (MPL %.3f): reported; the place is judged at half height' %
                  (lab, a['xmin_run'], a['xmin_ref']))
        for row in ('axis', 'wall'):
            v(c[row]['margin'] > R4_MARGIN, '%s R4 %s margin (rms against MPS, OG over against MPL) %.2f > %g'
              % (lab, row, c[row]['margin'], R4_MARGIN))
        for name, val, std, key, tol in (('axis outside the drop', a['out_mean'], a['out_std'], 'axis_out', R5_AXIS),
                                         ('wall', w['mean'], w['std'], 'wall', R5_WALL)):
            mm, ms = m[key + '_mean'], m[key + '_std']
            v(abs(val - mm) <= tol, '%s R5 %s mean %+.4f within %+.4f +- %.3f' % (lab, name, val, mm, tol))
            v(std <= R5_STD * ms, '%s R5 %s std %.4f <= %.4f' % (lab, name, std, R5_STD * ms))
        i = inlet(run)
        need = {'IGLOO': ('igloo',), 'ICE': ('ice',), 'both': ('ice', 'igloo')}[lab]
        miss = [s for s in need if s not in i]
        if miss:
            v(False, '%s R6 condensed phases: no output of %s' % (lab, ' and '.join(s.upper() for s in miss)))
        if i:
            ks = sum(s['k'] for s in i.values())
            v(abs(ks / K_PRESCRIBED - 1) <= K_TOL, '%s R6 prescribed loading: inflow share k of the 401 records, '
              'summed over the phases, %.6f = %.6f' % (lab, ks, K_PRESCRIBED))
        if 'ice' in i:
            r, b = i['ice']['r'], R6_FACTOR * m['r6_ice']
            v(np.abs(r - 1).max() <= b, '%s R6 ICE inlet (%s, k %.6f), 40 cells: %.4f..%.4f, |r - 1| <= %.4f'
              % (lab, i['ice']['phase'], i['ice']['k'], r.min(), r.max(), b))
        if 'igloo' in i:
            r = i['igloo']['r']
            lo, hi = 1 - R6_FACTOR * (1 - m['r6_igloo_lo']), 1 + R6_FACTOR * (m['r6_igloo_hi'] - 1)
            v(lo <= r[1:-1].min() and r[1:-1].max() <= hi, '%s R6 IGLOO inlet (%s, k %.6f), 38 inner cells: %.4f..%.4f '
              'within %.4f..%.4f' % (lab, i['igloo']['file'], i['igloo']['k'], r[1:-1].min(), r[1:-1].max(), lo, hi))
            for name, val, key in (('axis cell', r[0], 'r6_igloo_axis'), ('wall cell', r[-1], 'r6_igloo_wall')):
                b = R6_FACTOR * abs(m[key] - 1)
                v(abs(val - 1) <= b, '%s R6 IGLOO inlet, %s: %.4f, |r - 1| <= %.4f' % (lab, name, val, b))
        if 'ice' in i and 'igloo' in i:
            s, b = i['igloo']['rho'] / i['ice']['rho'], R6_FACTOR * abs(m['r6_rho'] - 1)
            v(abs(s - 1) <= b, '%s R6 SUM rho_p IGLOO / ICE %.4f, |s - 1| <= %.4f' % (lab, s, b))
    if {'IGLOO', 'ICE', 'both'} <= set(q):
        for row, key in (('axis', 'out_mean'), ('wall', 'mean')):
            lo, hi = sorted((q['IGLOO'][row][key], q['ICE'][row][key]))
            val = q['both'][row][key]
            v(lo <= val <= hi, 'cross R6 three solvers between the two: %s mean %+.4f in [%+.4f, %+.4f]'
              % ('axis (outside the drop)' if row == 'axis' else row, val, lo, hi))
        for row in ('axis', 'wall'):
            A, B = q['ICE'][row], q['IGLOO'][row]
            if A['X'].shape != B['X'].shape:
                v(False, 'cross R6 ICE - IGLOO gas Mach, %s: different cells' % row)
                continue
            d = (A['M'] - B['M'])[A['X'] > xmin]
            rms, b = np.sqrt((d ** 2).mean()), R6_FACTOR * MEASURED['cross']['ice_igloo_' + row]
            v(rms <= b, 'cross R6 ICE - IGLOO gas Mach, %s rms %.4f <= %.4f' % (row, rms, b))
    print('== %d passed, %d failed' % (v.passed, v.failed))
    return 1 if v.failed else 0


def measure(runs, ref, xmin):
    """The MEASURED block of reference runs of the three cases, to paste between the markers."""
    labels = dict(runs)
    zones = {(kind, Z): chang_zone(ref, 'Chang %s %s' % (kind, Z)) for kind in ('MPL', 'MPS', 'OG') for _, _, Z in ROWS}
    q = {}
    print('MEASURED = {')
    for lab, case in DEFAULT:
        run = labels[lab]
        c, p, i = chang_rows(run, zones, xmin), probes(run), inlet(run)
        q[lab] = c
        a, w = c['axis'], c['wall']
        print("    '%s': {  # %s" % (lab, case))
        for key, val, what in (
                ('axis_out_rms', a['out_rms'], 'R3: rms of run - MPL, axis cells outside the reference drop'),
                ('axis_exit', a['exit'], 'R3: axis exit M - MPL at the last reference station'),
                ('wall_rms', w['rms'], 'R3: rms of run - MPL, wall'),
                ('wall_exit', w['exit'], 'R3: wall exit M - MPL'),
                ('drop_half', a['half'], 'R3: half-height X of the drop - MPL (interpolated; a minimum is on a cell)'),
                ('drop_depth', a['depth'], 'R3: depth of the drop (max - min) - MPL'),
                ('axis_out_mean', a['out_mean'], 'R5: mean of run - MPL, axis outside the drop'),
                ('axis_out_std', a['out_std'], 'R5: its std'),
                ('wall_mean', w['mean'], 'R5: mean of run - MPL, wall'),
                ('wall_std', w['std'], 'R5: its std'),
                ('r2_axis', p['axis'], 'R2, reference only (R2_AXIS bounds it): axis probe M range, last block'),
                ('r2_wall', p['wall'], 'R2 (printed only): range of the wall exit probe M over the last block')):
            print("        '%s': %s,  # %s" % (key, 'None' if val is None else '%.6g' % val, what))
        if 'ice' in i:
            print("        'r6_ice': %.6g,  # R6: max |ratio - 1| of ICE's inlet cells (%s, %.4f..%.4f)"
                  % (np.abs(i['ice']['r'] - 1).max(), i['ice']['phase'], i['ice']['r'].min(), i['ice']['r'].max()))
        if 'igloo' in i:
            r = i['igloo']['r']
            print("        'r6_igloo_lo': %.6g,  # R6: lowest ratio of IGLOO's 38 inner inlet cells" % r[1:-1].min())
            print("        'r6_igloo_hi': %.6g,  # R6: highest ratio of IGLOO's 38 inner inlet cells" % r[1:-1].max())
            print("        'r6_igloo_axis': %.6g,  # R6: IGLOO's inlet cell on the axis" % r[0])
            print("        'r6_igloo_wall': %.6g,  # R6: IGLOO's inlet cell on the wall" % r[-1])
        if 'ice' in i and 'igloo' in i:
            print("        'r6_rho': %.6g,  # R6: SUM rho_p of IGLOO's projection over ICE's field"
                  % (i['igloo']['rho'] / i['ice']['rho']))
        print('    },')
    print("    'cross': {  # ICE alone against IGLOO alone, gas Mach on X > %g" % xmin)
    for row in ('axis', 'wall'):
        d = (q['ICE'][row]['M'] - q['IGLOO'][row]['M'])[q['ICE'][row]['X'] > xmin]
        print("        'ice_igloo_%s': %.6g,  # R6: rms of the difference, %s" % (row, np.sqrt((d ** 2).mean()), row))
    print('    },')
    print('}')
    return 0


def main():
    args, ref, xmin, runs, mode, req = sys.argv[1:], None, 0.2, [], 'report', False
    while args:
        a = args.pop(0)
        if a == '--ref':
            ref = args.pop(0)
        elif a == '--xmin':
            xmin = float(args.pop(0))
        elif a in ('--live', '--measure'):
            mode = a[2:]
        elif a == '--require-probes':
            req = True
        elif a == '--inject':
            return inject(args.pop(0))
        elif '=' in a:
            runs.append(tuple(a.split('=', 1)))
        else:
            sys.exit(__doc__)
    test = os.path.dirname(HERE)
    if ref is None:
        ref = os.path.join(test, 'JPL-Lagrangian-20micron', 'JPL_chang.dat')
    if mode == 'live':
        return live(runs, ref, xmin, req)
    if mode == 'measure':
        return measure(runs, ref, xmin)
    if not runs:
        runs = [(lab, os.path.join(test, case)) for lab, case in DEFAULT]
        runs = [(lab, p) for lab, p in runs if os.path.isfile(os.path.join(p, 'OUTPUT', 'gas-field.tec'))]
    return report(runs, ref, xmin)


if __name__ == '__main__':
    sys.exit(main())
