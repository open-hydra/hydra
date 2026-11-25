"""Readers and verdicts shared by the check.py of the hydra-MI2 test cases (test/mi2/*/check.py).

A check prints one PASS or FAIL line per claim and ends with done(), which exits 1 on any FAIL. Text only, no numpy.
"""
import re
import sys

T_REF = 298.15      # the datum of the formation enthalpies [K]
EMPTY = 1.0e-18     # below this density an ICE cell holds a floor, not mass
_count = [0, 0]


# --- verdicts --------------------------------------------------------------------------------------------------------

def claim(passed, msg):
    _count[0 if passed else 1] += 1
    print(('PASS  ' if passed else 'FAIL  ') + msg)


def judge(name, a, b, tol, what, scale=None):
    """PASS when |a - b| <= tol |scale| (scale defaults to b)."""
    den = abs(b if scale is None else scale)
    e = abs(a - b) / den if den > 0.0 else (0.0 if a == b else float('inf'))
    claim(e <= tol, '%s: %s %.12e vs %.12e, rel.err %.2e (tol %.0e)' % (name, what, a, b, e, tol))
    return e


def done(case):
    print('== %s: %d passed, %d failed' % (case, _count[0], _count[1]))
    sys.exit(1 if _count[1] or not _count[0] else 0)


# --- input.ini -------------------------------------------------------------------------------------------------------

def read_ini(path):
    """{section: {key: value}}; ';' starts a comment."""
    ini, sec = {}, None
    for raw in open(path):
        line = raw.split(';', 1)[0].strip()
        if line.startswith('['):
            sec = ini.setdefault(line.strip('[]').strip(), {})
        elif '=' in line and sec is not None:
            k, v = line.split('=', 1)
            sec[k.strip()] = v.strip()
    return ini


def values(ini, sec, key, n=None, default=None):
    """The numbers of a key, a single one broadcast to n."""
    s = ini.get(sec, {}).get(key, default)
    if s is None:
        die('input.ini: [%s] %s is missing' % (sec, key))
    v = [float(x) for x in str(s).split()]
    return v * n if n and len(v) == 1 else v


def die(msg):
    print('FAIL  ' + msg)
    sys.exit(1)


# --- the run's log ---------------------------------------------------------------------------------------------------

def read_log(path):
    return open(path, errors='replace').read()


def report(txt):
    """The [MI2-report] lines grouped per coupling step: cycle, step, t0, t1, p0, p1, E0, E1 and the rates a step
    applies (F, Q, M of IGLOO; iceF, iceQ, iceM of ICE; chF, chQ, chM of the channels), channel moved {name: kg/s}."""
    steps, cur = [], None
    for line in txt.split('\n'):
        if '[MI2-report]' not in line:
            continue
        f = line.split()
        tag = f[1]
        if tag == 'cycle':
            cur = {'cycle': int(f[2]), 'step': int(f[4]), 't0': float(f[6]), 't1': float(f[8]), 'moved': {}}
            steps.append(cur)
        elif cur is None:
            continue
        elif tag == 'gas-momentum':
            cur['p0'], cur['p1'] = [float(x) for x in f[2:5]], [float(x) for x in f[5:8]]
        elif tag == 'gas-energy':
            cur['E0'], cur['E1'] = float(f[2]), float(f[3])
        elif tag == 'igloo-to-gas':
            cur['F'], cur['Q'], cur['M'] = [float(x) for x in f[2:5]], float(f[5]), float(f[6])
        elif tag == 'ice-to-gas':
            cur['iceF'], cur['iceQ'] = [float(x) for x in f[2:5]], float(f[5])
        elif tag == 'ice-to-gas-mass':
            cur['iceM'] = float(f[2])
        elif tag == 'channels-to-gas':
            cur['chF'], cur['chQ'], cur['chM'] = [float(x) for x in f[2:5]], float(f[5]), float(f[6])
        elif tag == 'channel':
            cur['moved'][f[2]] = float(f[-1])
    return steps


def read_report(path):
    """report() of a log file."""
    return report(read_log(path))


def species(inp):
    """The gas species of INPUT/gas-phase.txt in order (the rows of the closed-domain balance)."""
    return [l.split()[0] for l in open(inp + 'gas-phase.txt').read().split('\n')[1:] if l.strip()]


def balance(txt, names):
    """{species: kg the gas gained} of the closed-domain mass balance."""
    m = re.search(r'closed-domain mass balance.*?\n.*?\n((?:\s+\d+\s+.*\n)+)', txt)
    out = {}
    for l in (m.group(1).split('\n') if m else []):
        f = l.split()
        if len(f) == 4 and f[0].isdigit():
            out[names[int(f[0]) - 1]] = float(f[2])
    return out


def energy_books(txt):
    """(handed to the gas, the gas's total energy gained) [J] of the closed-domain energy bookkeeping."""
    a = re.search(r'handed to gas .*?:\s+(\S+)', txt)
    b = re.search(r'gas total energy gained .*?:\s+(\S+)', txt)
    return (float(a.group(1)), float(b.group(1))) if a and b else (None, None)


def ice_dt(txt):
    """The time step of ICE's first iteration."""
    return float(re.search(r'ICE  \| Iter =\s+1 \| Time = \S+ \| Delta t = (\S+)', txt).group(1))


# --- ATLAS tables ----------------------------------------------------------------------------------------------------

def zones(path, col):
    """{zone: {T: column col}} of an ATLAS property table (integer-kelvin rows)."""
    out, cur = {}, None
    for line in open(path):
        if line.startswith('ZONE'):
            cur = out.setdefault(re.search(r'T="([^"]*)"', line).group(1), {})
            continue
        try:
            v = [float(x) for x in line.split()]
        except ValueError:
            continue
        if cur is not None and len(v) > col:
            cur[v[0]] = v[col]
    return out


def at(rows, T):
    """rows interpolated at T (linear between integer kelvins)."""
    lo = float(int(T))
    return rows[lo] if T == lo else rows[lo] + (rows[lo + 1.0] - rows[lo]) * (T - lo)


def records(path, material):
    """The records of a material in a products.txt (material <name> <kind> <terms> <datum> <h298> [tail], then the
    terms): {kind: {'gas': {species: nu}, 'h298', 'tail'}}. The tail: product, nu, h298 (and phase of a generated
    product) of a combustion, carry or release."""
    lines = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    out = {}
    for i, f in enumerate(lines):
        if f[0] == 'material' and f[1] == material:
            n = int(f[3])
            out[f[2]] = {'h298': float(f[5]), 'tail': f[6:],
                         'gas': {t[0]: float(t[1]) for t in lines[i + 1:i + 1 + n]}}
    return out


# --- solver output ---------------------------------------------------------------------------------------------------

def field(path):
    """Cell columns and cell volumes of an ICE Tecplot ASCII field (BLOCK: nodes, then cell-centred variables)."""
    lines = open(path).read().split('\n')
    names = lines[0].split('=')[1].split()
    I, J, K = [int(re.search(r'\b%s=(\d+)' % c, lines[1]).group(1)) for c in 'IJK']
    vals = []
    for l in lines[2:]:
        if not re.match(r'\s*[A-Za-z]', l):
            vals += [float(x) for x in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    X = [vals[q * nn:(q + 1) * nn] for q in range(3)]
    cols, off = {}, 3 * nn
    for n in names[3:]:
        cols[n], off = vals[off:off + nc], off + nc
    vol = []
    for k in range(K - 1):
        for j in range(J - 1):
            for i in range(I - 1):
                a = i + j * I + k * I * J
                vol.append((X[0][a + 1] - X[0][a]) * (X[1][a + I] - X[1][a]) * (X[2][a + I * J] - X[2][a]))
    return cols, vol


def ic_field(path):
    """ICB's initial field under the names of ICE's output (rp_11 -> rho_p1, ...)."""
    cols, vol = field(path)
    ren = {'rp_': 'rho_p', 'up_': 'u_p', 'vp_': 'v_p', 'wp_': 'w_p', 'Tp_': 'T_p', 'np_': 'n_p'}
    out = {}
    for k, v in cols.items():
        for a, b in ren.items():
            if k.startswith(a):
                k = b + k[len(a)]
        out[k] = v
    return out, vol


def family(cols, vol, fam, H, h298):
    """Mass, momentum (3), energy on the formation scale (h298 + H(T) - H(T_REF) + |u|^2/2 per kg), number and the
    occupied cells of an ICE family; cells below EMPTY are left out."""
    s = str(fam)
    r, T = cols['rho_p' + s], cols['T_p' + s]
    U = [cols['u_p' + s], cols['v_p' + s], cols['w_p' + s]]
    full = [c for c in range(len(r)) if r[c] > EMPTY and T[c] >= 1.0]
    M = sum(r[c] * vol[c] for c in full)
    P = [sum(r[c] * q[c] * vol[c] for c in full) for q in U]
    E = sum(r[c] * vol[c] * (h298 + at(H, T[c]) - at(H, T_REF) + 0.5 * sum(q[c] ** 2 for q in U)) for c in full)
    N = sum(cols['n_p' + s][c] * vol[c] for c in full) if 'n_p' + s in cols else 0.0
    return {'M': M, 'P': P, 'E': E, 'N': N, 'cells': full}


def outloc(path):
    """{parcel ID: (T, |u|, d)} of the end records of an IGLOO outloc file."""
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) == 12:
            try:
                out[int(f[8])] = (float(f[3]), float(f[4]), float(f[11]))
            except ValueError:
                pass
    return out


def read_snapshot(path):
    return snapshot(path)


def snapshot(path):
    """{parcel ID: X Y Z U V W T d m ID t} of an IGLOO snapshot file."""
    rows = {}
    for line in open(path):
        f = line.split()
        if len(f) == 11:
            try:
                v = [float(x) for x in f]
            except ValueError:
                continue
            rows[int(v[9])] = v
    return rows
