#!/usr/bin/env python3
"""ICE against IGLOO on the hot JPL nozzle: the freezing of the same alumina in both condensed solvers, on one gas.
usage: compare_freezing.py <run-dir> [--design full|partial] [--log <file>]

<run-dir> holds input.ini, INPUT/ and OUTPUT/ (a run_case.sh scratch directory or a case with its reference).
ICE: of OUTPUT/partE-field*.tec the one whose DATASETAUXDATA ICE_ITERATION is largest. IGLOO: the last sweep N of
OUTPUT/partL-*-sweep<N>.* (the untagged files are sweep 0): its trajectories, exit records and Eulerian field.
Cells (i, j): i = 1..190 along the nozzle, j = 1..40 from the axis; the inner region is j <= 10. --design is the
caller's: full = the drops freeze completely before the exit (JPL-3solver-hot-2micron), partial = they leave on the
plateau (JPL-3solver-hot-20micron). One PASS / FAIL / SKIP line per check, exit 1 on any FAIL:
  T0  the frozen-fraction columns exist (ICE f_p1, chi_p1; IGLOO f): without them the run is read as all liquid
  T1  liquid stage: delta_hot = max |T_ICE - T_IGLOO| / (T0 - mean T_p) <= 2 x 0.0155 over the inner cells where both
      carry mass, both are liquid, max(3, mollify-passes) or more cells upstream of IGLOO's first plateau cell of the
      row (IGLOO's Eulerian field is mollified over +-passes cells), cooled > 5 K
  T2  nucleation front: each IGLOO parcel's event cell (the midpoint of its last liquid and first plateau rows) in
      rows 1-10 within N_front of ICE's front (first chi_p1 >= 1/2 of the row);
      N_front = 3 + ceil(delta_hot (T0 - T_nuc) / (s dx)), s ICE's liquid slope over the 3 cells before its front
  T3  plateau temperature: ICE cells with 0 < f_p1 < 1 and IGLOO's Eulerian T_p in the inner cells whose +-passes
      box is crossed by plateau rows only (or none) at T_melt to 1e-9 relative; ICE plateau cells in rows 1-10 >= 50
      (full) / 260 (partial)
  T4  (full) end of the plateau: each inner parcel's end cell (midpoint of its last plateau row and the next row)
      within N_front + 1 of ICE's first f_p1 = 1 of the row
  T5  (full) every exit: IGLOO's exit-plane records (bc 400) T <= T_melt - 1 K; ICE's exit cells (bulk density above
      1e-3 of its peak) f_p1 = 1, T_p1 < T_melt; T5p (partial): those records and cells on the plateau, 0 < f < 1
  T6  (full) inner exit temperature, mass-flow weighted: |T_ICE - T_IGLOO| <= delta_hot (T0 - T_exit)
  C1  plateau slope df/dx per inner row, to 2 %: IGLOO's least-squares slope through each parcel's plateau rows in the
      row (its first and last left out) against ICE's through its f_p1 interpolated along the row at the same X, inside
      ICE's plateau cells but its first and last; the mean over the parcels of each
  C2  (full) every exit-plane record f = 1; C2p (partial) inner exit f, mass-flow weighted:
      |f_ICE - f_IGLOO| <= 0.02 (f_IGLOO - f0) + s_f N_front dx
Printed, not judged: the heat released per unit mass flow by each phase and delta over every row."""
import argparse
import glob
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mi2-checks'))
from tec_cc import read_cc                     # noqa: E402
from p8_solid_coherence import load_rows       # noqa: E402

INNER = 10              # rows j <= 10
MASS = 1e-3             # a cell carries mass: bulk density above 1e-3 of the solver's own peak
DROP = 5.0              # T1: cooled by more than 5 K
DELTA_TOL = 2 * 0.0155  # T1: twice the cold reference's inner maximum
HALF_ULP = 0.5e-6       # an IGLOO row on the plateau: |T - T_melt| <= half the last printed digit (p8)
TOL_T = 1e-9            # T3: relative
C1_TOL = 0.02
T3_COUNT = {'full': 50, 'partial': 260}


# ------------------------------------------------------------------------------------------------ inputs
def ini_key(run, sec, key):
    txt = open(os.path.join(run, 'input.ini')).read()
    m = re.search(r'^\[%s\][^\[]*?^\s*%s\s*=\s*([^;\n]*)' % (re.escape(sec), re.escape(key)), txt, re.M | re.S)
    return m.group(1).strip() if m else None


def material(run):
    """T_melt, T_nuc, h_fus, c_s from line 2 of INPUT/partL-phase.txt; c_l from the properties table."""
    tok = dict(t.split('=', 1) for t in open(os.path.join(run, 'INPUT', 'partL-phase.txt')).read().splitlines()[1].split()[2:]
               if '=' in t)
    rows = [l.split() for l in open(os.path.join(run, 'INPUT', 'partL-properties.dat')).read().splitlines()[4:] if l.strip()]
    return dict(Tm=float(tok.get('T-melt', 2327)), Tn=float(tok.get('T-nuc', 0.8 * float(tok.get('T-melt', 2327)))),
                hf=float(tok.get('h-fus', 'nan')), cs=float(tok.get('cp-solid', 'nan')), cl=float(rows[0][1]))


def mollify_passes(run):
    """IGLOO's binomial passes on its Eulerian output: [IGLOO-General] mollify (off -> 0), mollify-passes (default 8,
    IGLOO Lib_Mollify.f90 DEFAULT_MOLLIFY_PASSES)."""
    w = (ini_key(run, 'IGLOO-General', 'mollify') or 'on').lower()
    if w in ('off', 'false', 'no', '0', '.false.'):
        return 0
    n = ini_key(run, 'IGLOO-General', 'mollify-passes')
    return int(n) if n else 8


def trailer(path, key):
    with open(path, 'rb') as f:
        f.seek(0, 2)
        f.seek(max(0, f.tell() - 4096))
        tail = f.read().decode(errors='replace')
    m = re.search(r'DATASETAUXDATA\s+%s\s*=\s*"([^"]+)"' % key, tail)
    return float(m.group(1)) if m else None


def ice_field(run):
    files = sorted(glob.glob(os.path.join(run, 'OUTPUT', 'partE-field*.tec')))
    its = [(trailer(f, 'ICE_ITERATION'), f) for f in files]
    its = [t for t in its if t[0] is not None]
    if not its:
        sys.exit('no OUTPUT/partE-field*.tec with an ICE_ITERATION trailer in %s' % run)
    best = max(t[0] for t in its)
    return best, [f for it, f in its if it == best][0]


def igloo_files(run):
    sweeps = [int(m.group(1)) for f in glob.glob(os.path.join(run, 'OUTPUT', 'partL-*-sweep*.*'))
              for m in [re.search(r'-sweep(\d+)\.', f)] if m]
    n = max(sweeps) if sweeps else 0
    tag = '-sweep%d' % n if n else ''
    o = os.path.join(run, 'OUTPUT')
    return n, dict(traj=os.path.join(o, 'partL-trajectories-A%s.dat' % tag), exit=os.path.join(o, 'partL-outloc-A%s.dat' % tag),
                   euler=os.path.join(o, 'partL-euler1%s.tec' % tag))


def has_column(path, name):
    for line in open(path):
        if line.lstrip().lower().startswith('variables'):
            return name in re.findall(r'"([^"]*)"', line)
    return False


# ------------------------------------------------------------------------------------------------ geometry
class Mesh:
    """The block's k = 0 node plane (x axial, y radial): cell centroids, lengths along the row, exit face areas and a
    point locator (convex quads, binned)."""

    def __init__(self, x, y):
        self.x, self.y = x, y                      # (J + 1, I + 1)
        self.nj, self.ni = x.shape[0] - 1, x.shape[1] - 1
        q = [(x[:-1, :-1], y[:-1, :-1]), (x[:-1, 1:], y[:-1, 1:]), (x[1:, 1:], y[1:, 1:]), (x[1:, :-1], y[1:, :-1])]
        self.q = q
        self.xc = sum(a for a, _ in q) / 4
        self.yc = sum(b for _, b in q) / 4
        wx, wy = 0.5 * (x[:-1, :-1] + x[1:, :-1]), 0.5 * (y[:-1, :-1] + y[1:, :-1])
        ex, ey = 0.5 * (x[:-1, 1:] + x[1:, 1:]), 0.5 * (y[:-1, 1:] + y[1:, 1:])
        self.dx = np.hypot(ex - wx, ey - wy)       # length along the row (west to east face midpoints)
        cx = np.stack([a for a, _ in q], -1).reshape(-1, 4)
        cy = np.stack([b for _, b in q], -1).reshape(-1, 4)
        self.cx, self.cy = cx, cy
        self.x0, self.y0 = x.min(), y.min()
        self.hx = (x.max() - self.x0) / 400 + 1e-15
        self.hy = (y.max() - self.y0) / 150 + 1e-15
        self.bins = {}
        bx0 = np.floor((cx.min(1) - self.x0) / self.hx).astype(int)
        bx1 = np.floor((cx.max(1) - self.x0) / self.hx).astype(int)
        by0 = np.floor((cy.min(1) - self.y0) / self.hy).astype(int)
        by1 = np.floor((cy.max(1) - self.y0) / self.hy).astype(int)
        for c in range(cx.shape[0]):
            for bx in range(bx0[c], bx1[c] + 1):
                for by in range(by0[c], by1[c] + 1):
                    self.bins.setdefault((bx, by), []).append(c)

    def locate(self, px, py):
        """(i, j) 0-based of the cell holding each point; -1 outside (the nearest centroid within half a cell is taken
        for a point on the boundary)."""
        px, py = np.asarray(px, float), np.asarray(py, float)
        out = np.full(px.shape, -1, int)
        bx = np.floor((px - self.x0) / self.hx).astype(int)
        by = np.floor((py - self.y0) / self.hy).astype(int)
        key = bx * 1000003 + by
        order = np.argsort(key, kind='stable')
        ks = key[order]
        cuts = np.flatnonzero(np.diff(ks)) + 1
        for grp in np.split(order, cuts):
            if not len(grp):
                continue
            cand = self.bins.get((bx[grp[0]], by[grp[0]]))
            if not cand:
                continue
            cand = np.array(cand)
            X, Y = px[grp][:, None, None], py[grp][:, None, None]
            ax, ay = self.cx[cand][None], self.cy[cand][None]
            bxn, byn = np.roll(ax, -1, -1), np.roll(ay, -1, -1)
            cr = (bxn - ax) * (Y - ay) - (byn - ay) * (X - ax)
            scale = np.hypot(bxn - ax, byn - ay) * 1e-9
            inside = (cr >= -scale).all(-1)
            hit = inside.argmax(1)
            ok = inside[np.arange(len(grp)), hit]
            out[grp[ok]] = cand[hit[ok]]
        miss = np.flatnonzero(out < 0)
        if len(miss):
            xc, yc = self.xc.ravel(), self.yc.ravel()
            for m in miss:
                d = np.hypot(xc - px[m], yc - py[m])
                c = int(d.argmin())
                if d[c] <= 0.75 * max(self.dx.ravel()[c], 1e-12):
                    out[m] = c
        i = np.where(out >= 0, out % self.ni, -1)
        j = np.where(out >= 0, out // self.ni, -1)
        return i, j


# ------------------------------------------------------------------------------------------------ analysis
def analyse(run, design):
    """Every measured number and the verdicts: {'checks': [(id, verdict, text)], ...}."""
    R = dict(run=run, design=design, checks=[])

    def check(cid, ok, text):
        v = 'SKIP' if ok is None else ('PASS' if ok else 'FAIL')
        R['checks'].append((cid, v, text))

    mat = material(run)
    passes = mollify_passes(run)
    R['mollify_passes'] = passes
    Tm, Tn, hf, cs, cl = mat['Tm'], mat['Tn'], mat['hf'], mat['cs'], mat['cl']
    T0 = float(ini_key(run, 'inj', 'T0'))
    f0 = cl * (Tm - Tn) / hf
    R.update(mat=mat, T0=T0, f0=f0)
    it, efile = ice_field(run)
    E = read_cc(efile)
    sweep, F = igloo_files(run)
    L = read_cc(F['euler'])
    R.update(ice_file=os.path.basename(efile), ice_iteration=it, sweep=sweep)
    mesh = Mesh(E['x'][0], E['y'][0])
    nj, ni = mesh.nj, mesh.ni
    g = lambda d, n: d[n][0]                                       # noqa: E731  (J, I) cell arrays
    lname = lambda p: [n for n in L if n.startswith(p)][0]           # noqa: E731
    TE, rE, uE = g(E, 'T_p1'), g(E, 'rho_p1'), g(E, 'u_p1')
    TL, rL = L[lname('T')][0], L[lname('rho')][0]
    traj_f = has_column(F['traj'], 'f')
    exit_f = has_column(F['exit'], 'f')
    witness = 'f_p1' in E and 'chi_p1' in E and traj_f and exit_f
    fE = g(E, 'f_p1') if 'f_p1' in E else np.zeros_like(TE)
    chiE = g(E, 'chi_p1') if 'chi_p1' in E else np.zeros_like(TE)
    check('T0', witness, 'frozen-fraction columns: ICE f_p1 %s, chi_p1 %s; IGLOO trajectories f %s, exit records f %s '
          '(ICE %s, ICE_ITERATION %d; IGLOO sweep %d)' % ('f_p1' in E, 'chi_p1' in E, traj_f, exit_f,
                                                         R['ice_file'], it, sweep))
    mE, mL = rE > MASS * rE.max(), rL > MASS * rL.max()
    jj = np.arange(nj)[:, None] * np.ones((1, ni), int)
    ii = np.ones((nj, 1), int) * np.arange(ni)[None, :]
    inner = jj < INNER
    g3 = inner & mE & mL
    R['g3_cells'] = int(g3.sum())

    # ---- IGLOO parcels: rows, plateau rows, segment cells
    cols = ('X', 'Y', 'T') + (('f',) if traj_f else ())
    parts = load_rows(F['traj'], *cols)
    P = {}
    seg_x, seg_y, seg_pid, seg_plat, seg_endplat = [], [], [], [], []
    for pid, rows in parts.items():
        a = np.array([r[:len(cols)] for r in rows], float)
        X, Y, T = a[:, 0], a[:, 1], a[:, 2]
        f = a[:, 3] if traj_f else np.zeros(len(a))
        plat = (np.abs(T - Tm) <= HALF_ULP) & (f > 0)
        P[pid] = dict(X=X, Y=Y, T=T, f=f, plat=plat)
        if len(a) < 2:
            continue
        both = plat[:-1] & plat[1:]
        for s in (0.02, 0.5, 0.98):
            seg_x.append(X[:-1] + s * (X[1:] - X[:-1]))
            seg_y.append(Y[:-1] + s * (Y[1:] - Y[:-1]))
            seg_pid.append(np.full(len(a) - 1, pid))
            seg_plat.append(both)
            seg_endplat.append(plat[1:])
    seg_x, seg_y = np.concatenate(seg_x), np.concatenate(seg_y)
    seg_plat, seg_endplat = np.concatenate(seg_plat), np.concatenate(seg_endplat)
    si, sj = mesh.locate(seg_x, seg_y)
    okc = si >= 0
    crossed = np.zeros((nj, ni), bool)
    notplat = np.zeros((nj, ni), bool)
    crossed[sj[okc], si[okc]] = True
    notplat[sj[okc & ~seg_plat], si[okc & ~seg_plat]] = True
    iglo_plat_cells = crossed & ~notplat                      # crossed by plateau segments only
    # first plateau cell of each row: the most upstream cell where a segment ends on the plateau
    iP = np.full(nj, 10 ** 6, int)                            # no plateau row in a row: every cell is upstream
    sel = okc & seg_endplat
    for j in range(nj):
        s = sel & (sj == j)
        if s.any():
            iP[j] = si[s].min()
    R['iglo_first_plateau_cell'] = [int(v) if v < ni else None for v in iP[:INNER]]
    R['unlocated'] = int((~okc).sum())

    # events per parcel
    ev, en = [], []
    for pid, p in P.items():
        k = np.flatnonzero(p['plat'])
        if not len(k) or k[0] == 0:
            continue
        a, b = k[0] - 1, k[0]
        ev.append((pid, 0.5 * (p['X'][a] + p['X'][b]), 0.5 * (p['Y'][a] + p['Y'][b])))
        last = k[-1]
        if last + 1 < len(p['X']):
            en.append((pid, 0.5 * (p['X'][last] + p['X'][last + 1]), 0.5 * (p['Y'][last] + p['Y'][last + 1])))
    evi, evj = mesh.locate([e[1] for e in ev], [e[2] for e in ev]) if ev else (np.array([], int), np.array([], int))
    eni, enj = mesh.locate([e[1] for e in en], [e[2] for e in en]) if en else (np.array([], int), np.array([], int))
    R['parcels'] = len(P)
    R['events'] = [(ev[k][0], int(evi[k]), int(evj[k])) for k in range(len(ev))]
    R['ends'] = [(en[k][0], int(eni[k]), int(enj[k])) for k in range(len(en))]

    # ---- T1 liquid stage. IGLOO's Eulerian field is mollified over +-passes cells in i and j: its cells nearer than
    # that to IGLOO's first plateau cell carry the recalescence of the plateau, so they are left out with the 3 of the
    # front itself; the plan's 3-cell set is printed beside it.
    excl = max(3, passes)
    upstream = ii <= (iP[:, None] - excl)
    liqE = (chiE < 0.5) & (fE == 0)
    Tbar = 0.5 * (TE + TL)
    drop = T0 - Tbar
    t1_3 = g3 & liqE & (ii <= (iP[:, None] - 3)) & (drop > DROP)
    d3 = np.where(t1_3, np.abs(TE - TL) / np.where(drop > 0, drop, 1), 0.0)
    R.update(delta_3cells=float(d3.max()) if t1_3.any() else 0.0, t1_3cells=int(t1_3.sum()), t1_exclusion=excl)
    t1 = g3 & liqE & upstream & (drop > DROP)
    delta = np.where(t1, np.abs(TE - TL) / np.where(drop > 0, drop, 1), 0.0)
    dhot = float(delta.max()) if t1.any() else 0.0
    R.update(t1_cells=int(t1.sum()), delta_hot=dhot, t1_mask=t1)
    if t1.any():
        jm, im = np.unravel_index(int(np.argmax(delta)), delta.shape)
        where = 'at cell (%d, %d), x %.4f m: T_ICE %.3f K, T_IGLOO %.3f K' % (im + 1, jm + 1, mesh.xc[jm, im], TE[jm, im],
                                                                            TL[jm, im])
    else:
        where = 'no cell'
    check('T1', t1.any() and dhot <= DELTA_TOL, 'liquid stage: delta_hot %.5f <= %.4f over %d inner liquid cells %d or '
          'more cells upstream of IGLOO\'s first plateau cell (%s)' % (dhot, DELTA_TOL, int(t1.sum()), excl, where))
    allrows = mE & mL & liqE & upstream & (drop > DROP)
    dall = np.where(allrows, np.abs(TE - TL) / np.where(drop > 0, drop, 1), 0.0)
    R['delta_all_rows'] = float(dall.max()) if allrows.any() else 0.0

    # ---- ICE fronts, N_front per inner row
    iF = np.array([int(np.flatnonzero(chiE[j] >= 0.5)[0]) if (chiE[j] >= 0.5).any() else -1 for j in range(nj)])
    iE1 = np.array([int(np.flatnonzero(fE[j] == 1.0)[0]) if (fE[j] == 1.0).any() else -1 for j in range(nj)])
    Nf, dxF, sF = np.zeros(nj, int), np.zeros(nj), np.zeros(nj)
    for j in range(INNER):
        i = iF[j]
        if i < 3:
            Nf[j] = 3
            continue
        xs, ts = mesh.xc[j, i - 3:i], TE[j, i - 3:i]
        s = abs(np.polyfit(xs, ts, 1)[0])
        sF[j], dxF[j] = s, mesh.dx[j, i]
        Nf[j] = 3 + int(math.ceil(dhot * (T0 - Tn) / (s * dxF[j]))) if s > 0 else 3
    R.update(ice_front=iF[:INNER].tolist(), ice_solid=iE1[:INNER].tolist(), n_front=Nf[:INNER].tolist(),
             slope_front=sF[:INNER].tolist(), dx_front=dxF[:INNER].tolist())

    # ---- T2 nucleation front
    rows2 = [(pid, i, j) for pid, i, j in R['events'] if 0 <= j < INNER]
    bad2, gaps = [], []
    for pid, i, j in rows2:
        if iF[j] < 0:
            bad2.append('parcel %d row %d: no ICE front' % (pid, j + 1))
            continue
        gaps.append((i - iF[j], j))
        if abs(i - iF[j]) > Nf[j]:
            bad2.append('parcel %d row %d: event cell %d, ICE front %d, gap %d > N_front %d' % (pid, j + 1, i + 1, iF[j] + 1,
                                                                                             i - iF[j], Nf[j]))
    nofront = [j + 1 for j in range(INNER) if iF[j] < 0]
    R['t2_gaps'] = gaps
    per_row = [sum(1 for _, _, j in rows2 if j == r) for r in range(INNER)]
    txt2 = ('nucleation front: %d IGLOO events in rows 1-%d (per row %s); ICE front cells %s; N_front %s; gap event - '
            'front min %s max %s' % (len(rows2), INNER, per_row, [int(v) + 1 if v >= 0 else None for v in iF[:INNER]],
                                     Nf[:INNER].tolist(), min(g_ for g_, _ in gaps) if gaps else None,
                                     max(g_ for g_, _ in gaps) if gaps else None))
    if nofront:
        txt2 += '; rows without an ICE front %s' % nofront
    if bad2:
        txt2 += '; ' + '; '.join(bad2[:5]) + (' (+%d more)' % (len(bad2) - 5) if len(bad2) > 5 else '')
    check('T2', bool(rows2) and not nofront and not bad2, txt2)

    # ---- T3 plateau temperature
    platE = (fE > 0) & (fE < 1)
    errE = float((np.abs(TE[platE] - Tm) / Tm).max()) if platE.any() else 0.0
    nE = int((platE & inner).sum())
    # IGLOO's Eulerian T_p is judged where the whole mollifier support (+-passes cells) is crossed by plateau rows only
    # or by none
    okL = iglo_plat_cells | ~crossed
    supp = np.ones_like(okL)
    for dj in range(-passes, passes + 1):
        for di in range(-passes, passes + 1):
            sh = np.ones_like(okL)
            src = okL[max(0, dj):nj + min(0, dj), max(0, di):ni + min(0, di)]
            sh[max(0, -dj):nj - max(0, dj), max(0, -di):ni - max(0, di)] = src
            supp &= sh
    cellL = iglo_plat_cells & supp & inner
    errL = float((np.abs(TL[cellL] - Tm) / Tm).max()) if cellL.any() else 0.0
    R.update(t3_ice_cells=nE, t3_iglo_cells=int(cellL.sum()))
    need = T3_COUNT[design]
    check('T3', nE >= need and errE <= TOL_T and errL <= TOL_T,
          'plateau temperature: ICE %d cells with 0 < f_p1 < 1 (%d in rows 1-%d, >= %d), max |T_p1 - T_melt|/T_melt '
          '%.3e; IGLOO %d inner cells crossed by plateau rows only (%d of them with the mollifier\'s +-%d box so), max '
          '|T_p - T_melt|/T_melt %s; <= %.0e' % (int(platE.sum()), nE, INNER, need, errE, int((iglo_plat_cells & inner).sum()),
                                               int(cellL.sum()), passes, '%.3e' % errL if cellL.any() else 'not judged (no cell)',
                                               TOL_T))

    # ---- T4 end of the plateau
    if design == 'full':
        rows4 = [(pid, i, j) for pid, i, j in R['ends'] if 0 <= j < INNER]
        bad4, gaps4 = [], []
        for pid, i, j in rows4:
            if iE1[j] < 0:
                bad4.append('parcel %d row %d: no ICE f_p1 = 1' % (pid, j + 1))
                continue
            gaps4.append((i - iE1[j], j))
            if abs(i - iE1[j]) > Nf[j] + 1:
                bad4.append('parcel %d row %d: end cell %d, ICE %d, gap %d > %d' % (pid, j + 1, i + 1, iE1[j] + 1,
                                                                                 i - iE1[j], Nf[j] + 1))
        R['t4_gaps'] = gaps4
        txt4 = ('end of the plateau: %d inner parcels; ICE first f_p1 = 1 cells %s; tolerance N_front + 1 %s; gap min %s '
                'max %s' % (len(rows4), [int(v) + 1 if v >= 0 else None for v in iE1[:INNER]], (Nf[:INNER] + 1).tolist(),
                            min(g_ for g_, _ in gaps4) if gaps4 else None, max(g_ for g_, _ in gaps4) if gaps4 else None))
        if bad4:
            txt4 += '; ' + '; '.join(bad4[:5]) + (' (+%d more)' % (len(bad4) - 5) if len(bad4) > 5 else '')
        check('T4', bool(rows4) and not bad4, txt4)
    else:
        check('T4', None, 'end of the plateau (design partial)')

    # ---- exit records and exit cells
    ecols = ('Y', 'T', 'mdot', 'bc') + (('f',) if exit_f else ())
    recs = [r for rows in load_rows(F['exit'], *ecols).values() for r in rows]
    rec = np.array([r[:len(ecols)] for r in recs], float) if recs else np.zeros((0, len(ecols)))
    ry, rT, rm, rbc = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 3].astype(int)
    rf = rec[:, 4] if exit_f else np.zeros(len(rec))
    ex = rbc == 400
    wall = rbc == 301
    ie = ni - 1
    massx = mE[:, ie]
    y = mesh.y
    area = math.pi * (y[1:, ie + 1] ** 2 - y[:-1, ie + 1] ** 2)
    wE = rE[:, ie] * uE[:, ie] * area
    y10 = y[INNER, ie + 1]
    exin = ex & (ry <= y10)
    R.update(exit_records=int(ex.sum()), wall_records=int(wall.sum()), exit_inner_records=int(exin.sum()))

    def wmean(v, w):
        return float((v * w).sum() / w.sum()) if w.sum() > 0 else float('nan')

    TiE = wmean(TE[:INNER, ie], wE[:INNER])
    TiL = wmean(rT[exin], rm[exin])
    fiE = wmean(fE[:INNER, ie], wE[:INNER])
    fiL = wmean(rf[exin], rm[exin])
    R.update(exit_T_ice=TiE, exit_T_iglo=TiL, exit_f_ice=fiE, exit_f_iglo=fiL)
    if design == 'full':
        hot = ex & (rT > Tm - 1.0)
        badE = massx & ~((fE[:, ie] == 1.0) & (TE[:, ie] < Tm))
        check('T5', ex.any() and not hot.any() and not badE.any(),
              'every exit solid: %d exit-plane records, %d above T_melt - 1 K (max T %.3f K); %d wall records (bc 301, '
              'not judged; max T %.3f K); %d ICE exit cells with mass, %d not solid (f_p1 = 1, T_p1 < T_melt)'
              % (int(ex.sum()), int(hot.sum()), rT[ex].max() if ex.any() else float('nan'), int(wall.sum()),
                 rT[wall].max() if wall.any() else float('nan'), int(massx.sum()), int(badE.sum())))
    else:
        offL = ex & ~((np.abs(rT - Tm) <= HALF_ULP) & (rf > 0) & (rf < 1))
        badE = massx & ~((fE[:, ie] > 0) & (fE[:, ie] < 1))
        check('T5p', ex.any() and not offL.any() and not badE.any(),
              'no exit solid: %d exit-plane records, %d off the plateau (|T - T_melt| > %.1e K or f not in (0, 1)); %d '
              'ICE exit cells with mass, %d with f_p1 not in (0, 1)' % (int(ex.sum()), int(offL.sum()), HALF_ULP,
                                                                        int(massx.sum()), int(badE.sum())))

    # ---- T6 inner exit temperature
    if design == 'full':
        Tex = 0.5 * (TiE + TiL)
        tol6 = dhot * (T0 - Tex)
        R.update(t6_diff=TiE - TiL, t6_tol=tol6)
        check('T6', exin.any() and abs(TiE - TiL) <= tol6,
              'inner exit temperature: ICE %.4f K, IGLOO %.4f K (%d records, Y <= %.5f m), |difference| %.4f K <= '
              'delta_hot (T0 - T_exit) = %.5f x %.2f = %.4f K' % (TiE, TiL, int(exin.sum()), y10, abs(TiE - TiL), dhot,
                                                                 T0 - Tex, tol6))
    else:
        check('T6', None, 'inner exit temperature (design partial)')

    # ---- C1 plateau slope
    c1 = []
    if not (traj_f and 'f_p1' in E):
        check('C1', None, 'plateau slope: no frozen-fraction column')
    else:
        plat_rows = {}                                     # row -> [(X, f) of one parcel's plateau rows in it]
        for pid, p in P.items():
            k = np.flatnonzero(p['plat'] & (p['f'] < 1))
            if len(k) < 3:
                continue
            pi_, pj_ = mesh.locate(p['X'][k], p['Y'][k])
            for j in range(INNER):
                kk = k[pj_ == j]
                if len(kk) >= 3:
                    kk = kk[1:-1]
                if len(kk) >= 2:
                    plat_rows.setdefault(j, []).append((p['X'][kk], p['f'][kk]))
        for j in range(INNER):
            ic = np.flatnonzero(platE[j])
            ic = ic[1:-1] if len(ic) >= 3 else ic[:0]
            if len(ic) < 2 or j not in plat_rows:
                c1.append((j, None, None, None))
                continue
            # df/dx varies along the plateau and a parcel crosses a row over part of it: each parcel's slope is set
            # beside ICE's over the same points (ICE's f_p1 interpolated along the row at the parcel's X), inside ICE's
            # plateau cells but its first and last
            lo, hi = mesh.xc[j, ic].min(), mesh.xc[j, ic].max()
            sl, se = [], []
            for x, f in plat_rows[j]:
                m = (x >= lo) & (x <= hi)
                if m.sum() < 2 or np.ptp(x[m]) <= 0:
                    continue
                sl.append(np.polyfit(x[m], f[m], 1)[0])
                se.append(np.polyfit(x[m], np.interp(x[m], mesh.xc[j, ic], fE[j, ic]), 1)[0])
            if not sl:
                c1.append((j, None, None, None))
                continue
            sE, sL = float(np.mean(se)), float(np.mean(sl))
            c1.append((j, sE, sL, sL / sE - 1 if sE else float('inf')))
        done = [c for c in c1 if c[1] is not None]
        R['c1'] = c1
        worst = max(done, key=lambda c: abs(c[3])) if done else None
        check('C1', bool(done) and all(abs(c[3]) <= C1_TOL for c in done),
              'plateau slope df/dx: %d inner rows compared%s; worst row %s: ICE %s /m, IGLOO %s /m, %s (<= %.0f %%)'
              % (len(done), '' if len(done) == INNER else ' (rows without enough points: %s)'
                 % [c[0] + 1 for c in c1 if c[1] is None],
                 worst[0] + 1 if worst else None, '%.4f' % worst[1] if worst else '-', '%.4f' % worst[2] if worst else '-',
                 '%+.3f %%' % (100 * worst[3]) if worst else '-', 100 * C1_TOL))

    # ---- C2 / C2p
    if not exit_f:
        check('C2' if design == 'full' else 'C2p', None, 'exit f: no f column in the exit records')
    elif design == 'full':
        nb = int((ex & (rf != 1.0)).sum())
        check('C2', ex.any() and nb == 0, 'exit f = 1: %d exit-plane records, %d with f != 1 (min f %.6f)'
              % (int(ex.sum()), nb, rf[ex].min() if ex.any() else float('nan')))
    else:
        sf = float(np.mean([c[1] for c in c1 if c[1] is not None])) if any(c[1] is not None for c in c1) else 0.0
        ndx = float(max(Nf[j] * dxF[j] for j in range(INNER)))
        tol = 0.02 * (fiL - f0) + abs(sf) * ndx
        R.update(c2p_diff=fiE - fiL, c2p_tol=tol)
        check('C2p', exin.any() and abs(fiE - fiL) <= tol,
              'inner exit f: ICE %.5f, IGLOO %.5f, |difference| %.5f <= 0.02 (f_IGLOO - f0) + s_f N_front dx = 0.02 x '
              '%.4f + %.4f /m x %.3e m = %.5f (f0 %.4f)' % (fiE, fiL, abs(fiE - fiL), fiL - f0, abs(sf), ndx, tol, f0))

    # ---- printed, not judged: heat released per unit mass flow, delta over every row
    def e(T, f):
        return np.where(f <= 0, hf + cl * (T - Tm), np.where(f >= 1, cs * (T - Tm), (1 - f) * hf))

    wIn = rE[:, 0] * uE[:, 0] * math.pi * (y[1:, 0] ** 2 - y[:-1, 0] ** 2)
    qE = wmean(e(TE[:, 0], fE[:, 0]), wIn) - wmean(e(TE[:, ie], fE[:, ie]), wE)
    first = {pid: (p['T'][0], p['f'][0]) for pid, p in P.items()}
    ein = []
    exrows = load_rows(F['exit'], 'T', 'mdot', 'bc', *(('f',) if exit_f else ()))
    for pid, rows in exrows.items():
        for r in rows:
            if int(r[2]) == 400 and pid in first:
                fo = r[3] if exit_f else 0.0
                ein.append((float(e(np.array(first[pid][0]), np.array(first[pid][1]))),
                            float(e(np.array(r[0]), np.array(fo))), r[1]))
    qL = (sum(a * w for a, _, w in ein) - sum(b * w for _, b, w in ein)) / sum(w for _, _, w in ein) if ein else float('nan')
    R.update(q_ice=qE, q_iglo=qL)
    R['info'] = ['heat released per unit mass flow, inlet to exit plane: ICE %.4e J/kg, IGLOO %.4e J/kg (%+.3f %%)'
                 % (qE, qL, 100 * (qL / qE - 1) if qE else float('nan')),
                 'delta over every row (both carry mass, both liquid, upstream of IGLOO\'s plateau, cooled > 5 K): max %.5f'
                 % R['delta_all_rows'],
                 'T1 with the plan\'s set (3 or more cells upstream of IGLOO\'s first plateau cell; not judged: IGLOO mollifies '
                 'its Eulerian field over +-%d cells): delta %.5f over %d cells' % (passes, R['delta_3cells'], R['t1_3cells']),
                 'inner exit: T ICE %.3f / IGLOO %.3f K; f ICE %.5f / IGLOO %.5f; G3 cells %d; IGLOO parcels %d, events %d, '
                 'plateau ends %d; trajectory points outside every cell %d of %d' % (TiE, TiL, fiE, fiL, R['g3_cells'], len(P),
                                                                                   len(ev), len(en), R['unlocated'], len(okc))]
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('run')
    ap.add_argument('--design', choices=('full', 'partial'), default='full')
    ap.add_argument('--log', help='the run\'s log: its build stamp and end are printed')
    a = ap.parse_args()
    if a.log:
        txt = open(a.log, errors='replace').read()
        print('log %s: %s; GAME OVER %d' % (a.log, '; '.join(l.strip() for l in txt.splitlines()
                                                            if l.startswith('# [MI2] built from hydra')) or 'no stamp',
                                            txt.count('GAME OVER')))
    R = analyse(a.run, a.design)
    print('run %s, design %s: ICE %s (ICE_ITERATION %d), IGLOO sweep %d; T0 %g K, T_melt %g K, T_nuc %g K, f0 %.4f'
          % (a.run, a.design, R['ice_file'], R['ice_iteration'], R['sweep'], R['T0'], R['mat']['Tm'], R['mat']['Tn'],
             R['f0']))
    for line in R['info']:
        print('  ' + line)
    for cid, v, text in R['checks']:
        print('%-4s  %-3s %s' % (v, cid, text))
    nf = sum(1 for _, v, _ in R['checks'] if v == 'FAIL')
    print('== compare_freezing: %d passed, %d failed, %d skipped' % (sum(1 for _, v, _ in R['checks'] if v == 'PASS'), nf,
                                                                    sum(1 for _, v, _ in R['checks'] if v == 'SKIP')))
    sys.exit(1 if nf else 0)


if __name__ == '__main__':
    main()
