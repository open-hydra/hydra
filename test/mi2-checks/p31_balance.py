#!/usr/bin/env python3
"""p31_balance.py: gate p31's balance of a converged steady field, which step-independence alone does not see (a wrong
model that does not depend on the step passes p31's pair of runs). One line, exit 1 when the relative difference
exceeds tol. Reads only: the run's MESH/mesh.tec, OUTPUT/<field>, INPUT/ and the run's log.
  smoke   <run> <log> <field> <k> <tol>: ICE's smoke births. Family k's mass leaving the domain (its boundary cells'
          flux rho u.A through each boundary face whose ghost copies the cell, Lib_Ghost's ghost_rule: an outflow
          record, 400; a wall, 300/301 other than face 3, or an inlet, 401-403, that the cell's velocity points out of;
          the wedge faces, 200, and the axis mirror) against nu_s SUM burn = nu_s/(1 - nu_s) M, M the last [MI2-report] ice-to-gas-mass (the gas's
          share 1 - nu_s of every kg ICE burns), nu_s the burning material's smoke yield in kg per kg (its cap record's
          product in INPUT/part-products.txt).
  igloo   <run> <log> <field> <k> <tol>: IGLOO's smoke born in ICE (ext). Family k's mass leaving, as above, against the
          last [MI2-report] igloo-smoke-rate's mass.
  burnout <run> <log> <tol>: the steady burn-out's books. M (ice-to-gas-mass; the deck's material burns to gas only)
          against nu_g (SUM burn vol + SUM f out), from the field of family 1: burn the Beckstead rate n rho_m pi K X/(2 n)
          d^(3 - n) of every cell at or above T-ign (INPUT/part-models.txt), out the flux of each flagged cell's own
          state through each of its faces that it points out of (interior faces, and boundary faces whose ghost copies
          the cell, as above), f the fraction of it that burns at the cell's exit: a cell flagged when it holds drops
          (rho above 1e-25, n > 0) at or above T-ign within ICE's band of the removal size (d-burnout, else m_burnout's
          drop), f as ICE's burnout_fraction, its band read from the run's log (the setup's line "burnout_band =",
          ICE's one copy of it). The cell volumes from the mesh (divergence theorem over the faces' bilinear-quad area
          vectors).
  agree   <run> <log> <TA log>: the steady burn-out against the time-accurate one, which keeps the hard switch: the
          steady run's metal leaving through the outflow records (type 400), within 1e-6 of the mass to the gas (exit
          1 otherwise: the time-accurate switch lets nothing out; the steady band lets 1 - f of a band drop's outflow
          leave at the outlet face, 2.97e-8 of the mass to the gas at [0.90, 1] on the shipped deck at first order, and
          1e-6 is ten times that rounded up, as the gate's tolerances are), and, stated and not judged, the steady run's mass and energy to the gas (its last ice-to-gas-mass;
          ice-to-gas's energy) against the time-accurate cycle's mean and range, the smallest and the largest of its
          ICE blocks over the second half (ice-to-gas-mass + ice-removed-rate's mass; ice-to-gas's energy +
          ice-removed-rate's), each block once: the coupler reports a block's hand-off in its cycle and again in the
          frozen cycle after it, and the repeat is dropped.
  agree-in <run> <log> <TA log>: the same where the two models coincide: the outlet, and the steady mass and energy
          to the gas each inside the time-accurate run's range, judged (exit 1 outside).
  band    <log> <Mod_BC_Fluxes.f90> <text>...: the band in one form everywhere: the log's "burnout_band =" line, the
          source's `burnout_band(2) = [a_R8, b_R8]` and each text's "`burnout_band` = [a, b]"; exit 1 if one is missing
          or they differ.
A copying ghost makes the face state the boundary cell's, at first order and under MUSCL alike (the slope toward a
copy is zero)."""
OUT_TYPES = (400, 300, 301, 401, 402, 403)
BAND = None   # the band, read from the run's log (ICE's burnout_band); set only by a planted test
import math
import os
import re
import sys

import numpy as np


def read(path):
    """A one-zone Tecplot ASCII BLOCK file: (I, J, K) nodes and every variable (nodal x y z, the rest cell-centred)."""
    lines = open(path, errors='replace').read().split('\n')
    names = re.findall(r'"([^"]+)"', lines[0]) or lines[0].split('=', 1)[1].split()
    I, J, K = (int(x) for x in re.search(r'I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', ' '.join(lines[:2])).groups())
    vals = []
    for l in lines[2:]:
        if l.strip() and re.match(r'\s*[A-Za-z]', l):
            continue
        vals += l.split()
    out, off = {}, 0
    for i, n in enumerate(names):
        k = I * J * K if i < 3 else (I - 1) * (J - 1) * max(K - 1, 1)
        out[n] = np.array([float(v.replace('D', 'E')) for v in vals[off:off + k]])
        off += k
    return (I, J, K), out


def last(log, key, col=2):
    v = None
    for l in open(log, errors='replace'):
        if l.startswith(' [MI2-report] ' + key + ' '):
            v = float(l.split()[col])
    if v is None:
        sys.exit('p31_balance: no [MI2-report] %s line in %s' % (key, log))
    return v


class Mesh:
    def __init__(self, run):
        (I, J, K), m = read(run + '/MESH/mesh.tec')
        self.I, self.J, self.K = I, J, K
        self.P = np.stack([m['x'], m['y'], m['z']], axis=1).reshape(K, J, I, 3)

    def quad(self, ci, cj, ck, f):
        """The four nodes of face f (1 i-, 2 i+, 3 j-, 4 j+, 5 k-, 6 k+) of cell (ci, cj, ck), zero-based."""
        P = self.P
        if f in (1, 2):
            i = ci + (f - 1)
            return [P[ck, cj, i], P[ck, cj + 1, i], P[ck + 1, cj + 1, i], P[ck + 1, cj, i]]
        if f in (3, 4):
            j = cj + (f - 3)
            return [P[ck, j, ci], P[ck, j, ci + 1], P[ck + 1, j, ci + 1], P[ck + 1, j, ci]]
        k = ck + (f - 5)
        return [P[k, cj, ci], P[k, cj, ci + 1], P[k, cj + 1, ci + 1], P[k, cj + 1, ci]]

    def centre(self, ci, cj, ck):
        return self.P[ck:ck + 2, cj:cj + 2, ci:ci + 2].reshape(8, 3).mean(axis=0)

    def area_out(self, ci, cj, ck, f):
        """The face's area vector (bilinear quad: half the cross product of its diagonals), pointing out of the cell."""
        a, b, c, d = self.quad(ci, cj, ck, f)
        A = 0.5 * np.cross(c - a, d - b)
        if np.dot(A, (a + b + c + d) / 4.0 - self.centre(ci, cj, ck)) < 0:
            A = -A
        return A, (a + b + c + d) / 4.0

    def vol(self, ci, cj, ck):
        """Divergence theorem: V = 1/3 SUM_f x_f . A_f over the six faces, outward."""
        return sum(np.dot(x, A) for A, x in (self.area_out(ci, cj, ck, f) for f in range(1, 7))) / 3.0


def records(path, k=1, nfam=1):
    """(i, j, k, face, type) of family k's boundary records in a phase's bc file (1-based cells): the file lists every
    family's records in turn, family 1 first, the same faces for each."""
    out = []
    for l in open(path):
        w = l.split()
        if len(w) == 6 and all(re.match(r'-?\d+$', x) for x in w):
            out.append((int(w[1]), int(w[2]), int(w[3]), int(w[4]), int(w[5])))
    if len(out) % nfam:
        sys.exit('p31_balance: %d records in %s do not split over %d families' % (len(out), path, nfam))
    m = len(out) // nfam
    return out[(k - 1) * m:k * m]


def families(run, name):
    """The number of families of a field file (its rho_p<k> variables)."""
    head = open(run + '/OUTPUT/' + name, errors='replace').readline()
    return len(re.findall(r'rho_p\d+', head))


def field(run, name, k):
    (I, J, K), c = read(run + '/OUTPUT/' + name)
    shape = (K - 1, J - 1, I - 1)
    g = lambda v: c['%s_p%d' % (v, k)].reshape(shape)
    return g('rho'), np.stack([g('u'), g('v'), g('w')], axis=-1), g('T'), g('n')


def copies(t, f):
    """A boundary face whose ghost copies its cell when the cell's velocity points out (ICE's ghost_rule)."""
    return t in OUT_TYPES and not (t in (300, 301) and f == 3)


def outflow(run, name, k, bcfile):
    """Family k's mass leaving through the boundary faces whose ghost copies the cell, and entering through the
    outflow records (diagnostic)."""
    M = Mesh(run)
    rho, U, _, _ = field(run, name, k)
    out = inn = 0.0
    for i, j, kk, f, t in records(run + '/INPUT/' + bcfile, k, families(run, name)):
        if not copies(t, f):
            continue
        A, _ = M.area_out(i - 1, j - 1, kk - 1, f)
        q = rho[kk - 1, j - 1, i - 1] * np.dot(U[kk - 1, j - 1, i - 1], A)
        if q > 0:
            out += q
        elif t == 400:
            inn -= q
    return out, inn


def products_yield(path):
    """The burning material's smoke yield, kg of its cap product per kg (the material record's product field)."""
    for l in open(path):
        w = l.split()
        if len(w) >= 8 and w[0] == 'material' and w[2] == 'combustion':
            return float(w[7]), w[6]
    sys.exit('p31_balance: no burning material with a product in ' + path)


def judge(what, a, b, tol):
    rel = abs(a - b) / max(abs(b), 1e-300)
    print('%s: %.10e against %.10e, relative %.3e (tol %.0e)' % (what, a, b, rel, tol))
    sys.exit(0 if rel <= tol else 1)


def models(path):
    for l in open(path):
        w = l.split()
        if len(w) >= 3 and not l.startswith('#') and w[1] == 'combustion':
            kv = dict(x.split('=', 1) for x in w[3:] if '=' in x)
            return w[2], kv
    sys.exit('p31_balance: no combustion record in ' + path)


def rho_mat(path, T):
    """The material's density at T from its properties table (Temperature, Cp, Density, ...)."""
    tab = []
    for l in open(path):
        w = l.split()
        if len(w) >= 3 and re.match(r'-?[0-9.]', w[0]):
            tab.append((float(w[0]), float(w[2])))
    tab = np.array(tab)
    return np.interp(T, tab[:, 0], tab[:, 1])


def band_of(log):
    """The band ICE ran with: its setup line ' - Burnout --> steady, burnout_band = lo hi x ...'."""
    for l in open(log, errors='replace'):
        m = re.search(r'burnout_band =\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)', l)
        if m:
            return float(m.group(1)), float(m.group(2))
    sys.exit('p31_balance: no "burnout_band =" line in %s' % log)


def fraction(m, mb, band):
    """ICE's burnout_fraction: 1 at or below band[0]**3 mb, 0 at or above band[1]**3 mb, a cubic step in the
    diameter ratio between."""
    lo, hi = band
    if m <= lo ** 3 * mb:
        return 1.0
    if m >= hi ** 3 * mb:
        return 0.0
    t = (hi - (m / mb) ** (1.0 / 3.0)) / (hi - lo)
    return t * t * (3.0 - 2.0 * t)


def burnout(run, log, tol):
    model, kv = models(run + '/INPUT/part-models.txt')
    if model != 'Beckstead':
        sys.exit('p31_balance burnout: the books are written for Beckstead only (%s)' % model)
    Kb, nb, X, Tign = (float(kv[x]) for x in ('K-burn', 'n-burn', 'X-eff', 'T-ign'))
    db = float(kv.get('d-burnout', 0.0))
    band = BAND or band_of(log)
    M = Mesh(run)
    rho, U, T, n = field(run, 'part-field.tec', 1)
    rm = rho_mat(run + '/INPUT/part-properties.dat', T)
    nk, nj, ni = rho.shape
    frac = np.zeros(rho.shape)
    burn = out = 0.0
    for ck in range(nk):
        for cj in range(nj):
            for ci in range(ni):
                r, nn, t = rho[ck, cj, ci], n[ck, cj, ci], T[ck, cj, ci]
                if not (r > 1.0000001e-25 and nn > 0.0):
                    continue
                d = (6.0 * r / (math.pi * rm[ck, cj, ci] * nn)) ** (1.0 / 3.0)
                if t >= Tign:
                    burn += nn * rm[ck, cj, ci] * math.pi * Kb * X / (2.0 * nb) * d ** (3.0 - nb) * M.vol(ci, cj, ck)
                    mb = math.pi / 6.0 * rm[ck, cj, ci] * db ** 3 if db > 0.0 else 1e-20
                    frac[ck, cj, ci] = fraction(r, mb * nn, band)
    bc = {}
    for i, j, kk, f, t in records(run + '/INPUT/part-bc.txt', 1, families(run, 'part-field.tec')):
        bc[(i - 1, j - 1, kk - 1, f)] = t
    for ck, cj, ci in zip(*np.nonzero(frac > 0.0)):
        for f in range(1, 7):
            inside = {1: ci > 0, 2: ci < ni - 1, 3: cj > 0, 4: cj < nj - 1, 5: ck > 0, 6: ck < nk - 1}[f]
            t = bc.get((ci, cj, ck, f))
            if not inside and not copies(t, f):
                continue          # the wedge faces (200) and the axis mirror: nothing leaves through them
            A, _ = M.area_out(ci, cj, ck, f)
            q = rho[ck, cj, ci] * np.dot(U[ck, cj, ci], A)
            if q > 0:
                out += frac[ck, cj, ci] * q
    Mg = last(log, 'ice-to-gas-mass')
    print('cells flagged %d (%d of them at f = 1); SUM burn vol %.10e, SUM f out of the flagged cells %.10e kg/s' % (
        (frac > 0.0).sum(), (frac >= 1.0).sum(), burn, out))
    judge('ice-to-gas-mass against SUM burn + SUM f out', Mg, burn + out, tol)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'smoke':
        run, log, name, k, tol = sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6])
        nu, prod = products_yield(run + '/INPUT/part-products.txt')
        o, i = outflow(run, name, k, 'part-bc.txt')
        Mg = last(log, 'ice-to-gas-mass')
        print('smoke %s (nu_s %.10e): leaving %.10e kg/s, entering through the outflow records %.3e' % (prod, nu, o, i))
        judge('smoke leaving against nu_s/(1 - nu_s) ice-to-gas-mass', o, nu / (1.0 - nu) * Mg, tol)
    elif mode == 'igloo':
        run, log, name, k, tol = sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6])
        o, i = outflow(run, name, k, name.split('-')[0] + '-bc.txt')
        print('ICE smoke from IGLOO: leaving %.10e kg/s, entering through the outflow records %.3e' % (o, i))
        judge('smoke leaving against igloo-smoke-rate', o, last(log, 'igloo-smoke-rate'), tol)
    elif mode == 'burnout':
        burnout(sys.argv[2], sys.argv[3], float(sys.argv[4]))
    elif mode == 'band':
        log, src, texts = sys.argv[2], sys.argv[3], sys.argv[4:]
        got = [('the log', band_of(log))]
        m = re.search(r'burnout_band\(2\)\s*=\s*\[\s*([0-9.]+)_R8\s*,\s*([0-9.]+)_R8\s*\]', open(src).read())
        got.append((os.path.basename(src), (float(m.group(1)), float(m.group(2))) if m else None))
        for t in texts:
            m = re.findall(r'`burnout_band` = \[([0-9.]+), ([0-9.]+)\]', open(t).read())
            got.append((os.path.basename(t), (float(m[0][0]), float(m[0][1])) if m and len(set(m)) == 1 else None))
        for w, b in got:
            print('band in %s: %s' % (w, '%g, %g' % b if b else 'MISSING (or two values)'))
        ok = all(b for _, b in got) and len(set(b for _, b in got)) == 1
        print('the band is %s in all %d places' % ('the same' if ok else 'NOT the same', len(got)))
        sys.exit(0 if ok else 1)
    elif mode in ('agree', 'agree-in'):
        run, log, ta = sys.argv[2], sys.argv[3], sys.argv[4]
        Mg, Eg = last(log, 'ice-to-gas-mass'), last(log, 'ice-to-gas', 5)
        series, energy, rm = [], [], None
        for l in open(ta, errors='replace'):
            if l.startswith(' [MI2-report] ice-removed-rate'):
                rm = [float(x) for x in l.split()[2:7]]
            elif l.startswith(' [MI2-report] ice-to-gas ') and rm is not None:
                energy.append(float(l.split()[5]) + rm[4])
            elif l.startswith(' [MI2-report] ice-to-gas-mass') and rm is not None:
                series.append(float(l.split()[2]) + rm[0])
        if len(series) < 4 or len(energy) != len(series):
            sys.exit('p31_balance agree: %d time-accurate blocks (%d energies) in %s' % (len(series), len(energy), ta))
        # one value per ICE block: the coupler reports a block's hand-off in its cycle and again in the next, where ICE
        # is frozen (ICE-iter = 0) and the pair repeats bit for bit; the repeat is dropped
        nrep = len(series)
        keep = [k for k in range(nrep) if k == 0 or (series[k], energy[k]) != (series[k - 1], energy[k - 1])]
        series, energy = [series[k] for k in keep], [energy[k] for k in keep]
        print('time-accurate: %d hand-off reports, %d ICE blocks once the repeat of each in the frozen cycle after it '
              'is dropped' % (nrep, len(series)))
        if len(series) < 4:
            sys.exit('p31_balance agree: %d time-accurate ICE blocks in %s' % (len(series), ta))
        half, ehalf = series[len(series) // 2:], energy[len(energy) // 2:]
        mean, emean = sum(half) / len(half), sum(ehalf) / len(ehalf)
        print('time-accurate: %d blocks, the second half %d: mass mean %.10e in [%.10e, %.10e]; energy mean %.6f in '
              '[%.6f, %.6f] W' % (len(series), len(half), mean, min(half), max(half), emean, min(ehalf), max(ehalf)))
        M = Mesh(run)
        rho, U, _, _ = field(run, 'part-field.tec', 1)
        out = 0.0
        for i, j, kk, f, t in records(run + '/INPUT/part-bc.txt', 1, families(run, 'part-field.tec')):
            if t == 400:
                A, _ = M.area_out(i - 1, j - 1, kk - 1, f)
                out += max(rho[kk - 1, j - 1, i - 1] * np.dot(U[kk - 1, j - 1, i - 1], A), 0.0)
        print('steady: metal leaving through the outflow records %.3e kg/s (%.3e of the mass to the gas)' % (out, out / Mg))
        if not out <= 1e-6 * abs(Mg):
            print('the steady run lets metal out of the outlet')
            sys.exit(1)
        print('steady: mass to the gas %.10e (%+.3e of the time-accurate mean), energy %.6f W (%+.3e)' % (
            Mg, Mg / mean - 1, Eg, Eg / emean - 1))
        inside = min(half) <= Mg <= max(half) and min(ehalf) <= Eg <= max(ehalf)
        word = 'judged' if mode == 'agree-in' else 'stated'
        print('%s: the steady mass to the gas %+.3e of the time-accurate mean (its range %+.3e to %+.3e), the '
              'energy %+.3e (%+.3e to %+.3e): %s' % (word, Mg / mean - 1, min(half) / mean - 1, max(half) / mean - 1,
                                                     Eg / emean - 1, min(ehalf) / emean - 1, max(ehalf) / emean - 1,
                                                     'inside both ranges' if inside else 'outside'))
        sys.exit(0 if mode == 'agree' or inside else 1)
    else:
        sys.exit('p31_balance.py: smoke, igloo, burnout, agree, agree-in or band')
