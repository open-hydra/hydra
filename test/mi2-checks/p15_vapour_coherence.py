#!/usr/bin/env python3
"""The checks of gate p15: ICE and IGLOO read the vapour data of a condensed phase from the same columns of the
property table ATLAS GPB writes (<phase>-properties.dat), and take the same values from them.

  p15_vapour_coherence.py static <hydra-root>
  p15_vapour_coherence.py run <igloo-run-dir> <igloo-log> <ice-run-dir> <ice-log>

static: the two readers' sources in the checkout (IGLOO src/lib/Lib_Thermodynamics.f90, ICE src/lib/io/Load_Table.f90)
  S1  the column names each classifier (classify_table_tokens) accepts, every spelling: the same in both
  S2  the grammar both print on a refused header: the same text
  S3  the shared reader routines (checks of the header, the nodes, the columns, Psat, the vapour columns, the positive
      columns, and their messages): the same code in both, comments, blanks, continuations and case aside
  S4  every column of the property tables ATLAS wrote for the cases of test/ (their VARIABLES lines) is a name both
      readers accept; the names both accept that no table carries are printed as the expected gap
run: the same table, INPUT/part-properties.dat of test/evap-box-ta (water with Psat, Lv, Mv and Cpv columns), read by
  IGLOO in one run and by ICE in another; each solver's setup lines and the coupler's datum lines against the table's
  columns, taken by name and interpolated as both solvers do (linear between the integer-kelvin rows):
  R0  both runs clean, the same bytes read, and each solver says it took the vapour data from the table
  R1  Lv(298.15 K): IGLOO's line, ICE's line, and the coupler's "Lv (already inside the source)" for each solver
  R2  Lv(boiling-temperature) in both solvers' lines
  R3  Psat(boiling-temperature)/boiling-pressure in both solvers' lines
  R4  the datum hOff = h(Tmin) - cp Tmin of the Enthalpy_abs column: the coupler's line for IGLOO, ICE's line
  R5  Mv and Cpv(boiling-temperature): ICE prints the values; IGLOO prints only that it took both columns (a limit of
      its log, said in the line)
  R6  the coupler's datum offsets of the two solvers differ by the table's hOff: ICE credits the material cp T, IGLOO
      its absolute enthalpy, and both the same latent heat
One PASS/FAIL line per check (INFO lines are not checks); exit 1 on any FAIL."""
import hashlib
import os
import re
import sys

T_REF = 298.15
PATM = 101325.0
# Measured on hydra 220487a (ifx), evap-box-ta's table, two runs identical: the deviation from the table is the print
# precision, at most 2.4e-6 relative for six figures (ES12.5), 2.6e-5 for the Psat ratio (F7.4) and 2.3e-8 for eight
# figures (ES16.7, R6 2.1e-8). RED (Lv renamed in ICE's copy, latent-heat 2.26e6 given to ICE): R1 7.5e-2, R6 1.1e-2.
TOL5 = 2e-5          # a value printed with six figures (ES12.5)
TOL_RATIO = 1e-4     # Psat ratio printed as F7.4 (5e-5 absolute)
TOL8 = 1e-6          # a value printed with eight figures (ES16.7)
npass = nfail = 0


def ok(msg):
    global npass
    npass += 1
    print('PASS  ' + msg)


def ko(msg):
    global nfail
    nfail += 1
    print('FAIL  ' + msg)


def relerr(a, b):
    den = max(abs(a), abs(b))
    return abs(a - b) / den if den > 0.0 else 0.0


# ----------------------------------------------------------------------------------------------------------- static
def strip_fortran(text):
    """Comments, continuations, blanks and case out of a Fortran text (string literals kept)."""
    out = []
    for line in text.split('\n'):
        q, o = None, ''
        for c in line:
            if q:
                o += c
                if c == q:
                    q = None
            elif c in '"\'':
                q, o = c, o + c
            elif c == '!':
                break
            else:
                o += c
        out.append(o)
    s = ' '.join(out).replace('&', ' ')
    return re.sub(r'\s+', ' ', s).strip().lower()


def routines(text):
    """{name: source} of the module procedures of a Fortran file."""
    pat = r'^[ \t]*(?:pure[ \t]+|elemental[ \t]+)*(?:logical[ \t]+|integer[ \t]+)?(?:subroutine|function)[ \t]+(\w+)' \
          r'.*?^[ \t]*end[ \t]+(?:subroutine|function)[ \t]+\1\b'
    return {m.group(1).lower(): m.group(0) for m in re.finditer(pat, text, re.S | re.M | re.I)}


def labels(src):
    """The case labels of classify_table_tokens: [tuple of spellings] in source order."""
    found = [tuple(re.findall(r"'([^']*)'", m)) for m in re.findall(r"^\s*case\s*\(([^)]*)\)", src, re.M)]
    return [t for t in found if t]


def grammar(text):
    m = re.search(r"character\(len=\*\), parameter :: grammar = &(.*?)\n\s*\n", text, re.S)
    if not m:
        return None
    return ''.join(re.findall(r"'([^']*)'", m.group(1)))


def table_names(path):
    for line in open(path, errors='replace'):
        if line.startswith('VARIABLES'):
            return [n.strip(' "\n') for n in line.split('=', 1)[1].split(',')]
    return []


def static(root):
    fi = os.path.join(root, 'src/lib/IGLOO/src/lib/Lib_Thermodynamics.f90')
    fc = os.path.join(root, 'src/lib/ICE/src/lib/io/Load_Table.f90')
    for f in (fi, fc):
        if not os.path.exists(f):
            ko('S0 reader source not checked out: %s' % f)
            return
    ti, tc = open(fi).read(), open(fc).read()
    ri, rc = routines(ti), routines(tc)
    li = labels(ri.get('classify_table_tokens', ''))
    lc = labels(rc.get('classify_table_tokens', ''))
    si, sc = {s for t in li for s in t}, {s for t in lc for s in t}
    same = bool(li) and {frozenset(t) for t in li} == {frozenset(t) for t in lc}   # the branches in any order
    (ok if same else ko)('S1 the column names both classifiers accept, every spelling: %s%s' % (
        ' '.join(t[0] for t in li), '' if same else ' | IGLOO only: %s; ICE only: %s; grouped differently: %s'
        % (' '.join(sorted(si - sc)) or '-', ' '.join(sorted(sc - si)) or '-', 'yes' if si == sc else '-')))
    gi, gc = grammar(ti), grammar(tc)
    (ok if gi and gi == gc else ko)('S2 the header grammar both print on a refusal is the same text (%d characters)'
                                    % len(gi or ''))
    shared = sorted((set(ri) & set(rc)) - {'classify_table_tokens'})
    differ = [r for r in shared if strip_fortran(ri[r]) != strip_fortran(rc[r])]
    (ok if len(shared) >= 10 and not differ else ko)(
        'S3 %d shared reader routines, the same code in both (comments, blanks and case aside): %s%s' % (
            len(shared), ' '.join(shared), '' if not differ else ' | DIFFER: ' + ' '.join(differ)))
    accepted = {s for t in li for s in t} & {s for t in lc for s in t}
    canon = [t[0] for t in li]
    written = {}
    tdir = os.path.join(root, 'test')
    for case in sorted(os.listdir(tdir)):
        inp = os.path.join(tdir, case, 'INPUT')
        if not os.path.isdir(inp) or os.path.islink(os.path.join(tdir, case, 'INPUT')):
            continue
        for f in sorted(os.listdir(inp)):
            if f.endswith('properties.dat'):
                for n in table_names(os.path.join(inp, f))[1:]:
                    written.setdefault(n, []).append(case)
    unread = sorted(n for n in written if n not in accepted)
    (ok if written and not unread else ko)(
        'S4 every column of the property tables ATLAS wrote for the %d cases of test/ is read by both: %s%s' % (
            len({c for v in written.values() for c in v}), ' '.join(sorted(written)),
            '' if not unread else ' | READ BY NEITHER OR ONE: ' + ' '.join(unread)))
    gap = [n for n in canon if n not in written and n not in ('Enthalpy',)]
    gap = [n for n in gap if not any(s in written for t in li if t[0] == n for s in t)]
    print('INFO  expected gap: both readers accept %s, which no ATLAS table carries (GPB writes none of them for a '
          'dispersed phase; register 3.42): a property they give comes from the solvers\' own keys' % ', '.join(gap))
    io_py = os.path.join(root, 'utils/ATLAS/src/GPB/condensed/io.py')
    if os.path.exists(io_py):
        src = open(io_py).read()
        i = src.find("'dispersed' in type")
        j = src.find('\n        else:', i)
        if i < 0 or j < 0:
            print('INFO  utils/ATLAS GPB writer (io.py): its dispersed branch not found')
        else:
            seg = src[i:j]
            names = sorted(set(re.findall(r'\\"(\w+)\\"', seg)) | ({'Enthalpy', 'Enthalpy_abs'} if 'h_label' in seg else set()))
            print('INFO  utils/ATLAS GPB writer (io.py), dispersed layout: %s' % ' '.join(names))
    else:
        print('INFO  utils/ATLAS not checked out: the names ATLAS writes are taken from the stored tables only')


# -------------------------------------------------------------------------------------------------------------- run
def table(path):
    """{name: {T: value}} of the first zone of a property table, columns by name."""
    names, rows = None, {}
    for line in open(path):
        if line.startswith('VARIABLES'):
            names = [n.strip(' "\n') for n in line.split('=', 1)[1].split(',')]
            continue
        if line.startswith('ZONE') and rows:
            break
        try:
            v = [float(x) for x in line.split()]
        except ValueError:
            continue
        if names and len(v) == len(names):
            rows[v[0]] = v
    return {n: {T: r[i] for T, r in rows.items()} for i, n in enumerate(names)}


def at(col, T):
    lo = float(int(T))
    return col[lo] if T == lo else col[lo] + (col[lo + 1.0] - col[lo]) * (T - lo)


def ini_value(path, sec, key):
    cur = None
    for raw in open(path):
        line = raw.split(';', 1)[0].strip()
        if line.startswith('['):
            cur = line.strip('[]').strip()
        elif cur == sec and '=' in line and line.split('=', 1)[0].strip() == key:
            return float(line.split('=', 1)[1].split()[0])
    return None


def num(pat, txt):
    m = re.search(pat, txt)
    return float(m.group(1)) if m else None


def judge(name, what, got, want, tol):
    if got is None:
        ko('%s %s: not printed' % (name, what))
        return
    e = relerr(got, want)
    (ok if e <= tol else ko)('%s %s %.8e, table %.8e: rel.err %.1e (tol %.0e)' % (name, what, got, want, e, tol))


def reconciliation(txt, solver):
    """(Lv, credited, offset) of the coupler's datum block of the solver's material 1, None when absent."""
    i = txt.find('enthalpy-datum reconciliation')
    if i < 0:
        return None
    block = txt[i:txt.find('HYDRA-MI2 solvers', i)]
    lv = num(r'Lv \(already inside the source\)\s*:\s*(\S+)', block)
    cr = num(r'credited by %s\s*\(h \+ Lv\)\s*:\s*(\S+)' % solver, block)
    off = num(r'=> datum offset applied by MI2\s*:\s*(\S+)', block)
    return None if None in (lv, cr, off) else (lv, cr, off)


def run(ig_run, ig_log, ic_run, ic_log):
    ti = open(ig_log, errors='replace').read()
    tc = open(ic_log, errors='replace').read()
    fi, fc = ig_run + '/INPUT/part-properties.dat', ic_run + '/INPUT/part-properties.dat'
    clean = all('GAME OVER' in t and '[ERROR]' not in t for t in (ti, tc))
    md5 = [hashlib.md5(open(f, 'rb').read()).hexdigest() for f in (fi, fc)]
    took_i = re.search(r'\[material 1\] vapour data tabulated from properties\.dat: Lv, Mv, Cpv', ti) is not None
    took_c = re.search(r'vapour data from the table: Lv\(298\.15 K\)', tc) is not None
    (ok if clean else ko)('R0 both runs clean (GAME OVER, no [ERROR] line)')
    (ok if md5[0] == md5[1] else ko)('R0 the same table read by both: md5 %s / %s' % (md5[0][:8], md5[1][:8]))
    (ok if took_i and took_c else ko)('R0 each solver takes its vapour data from the table: IGLOO %s, ICE %s' % (
        'yes' if took_i else 'NO', 'yes' if took_c else 'NO'))
    tab = table(fi)
    for n in ('Enthalpy_abs', 'Cp', 'Psat', 'Lv', 'Mv', 'Cpv'):
        if n not in tab:
            ko('R0 the table has no %s column' % n)
            return
    Tb_i = ini_value(ig_run + '/input.ini', 'IGLOO-Properties', 'Tboil')
    Tb_c = ini_value(ic_run + '/input.ini', 'ICE-Physics', 'boiling-temperature')
    Tmin = min(tab['Lv'])
    lv_ref, lv_bi, lv_bc = at(tab['Lv'], T_REF), at(tab['Lv'], Tb_i), at(tab['Lv'], Tb_c)

    judge('R1', 'IGLOO Lv(298.15 K)', num(r'\[material 1\] vapour data tabulated.*?Lv\(298\.15 K\) =\s*(\S+)', ti), lv_ref,
          TOL5)
    judge('R1', 'ICE   Lv(298.15 K)', num(r'vapour data from the table: Lv\(298\.15 K\) = (\S+) J/kg', tc), lv_ref, TOL5)
    ri, rc = reconciliation(ti, 'IGLOO'), reconciliation(tc, 'ICE')
    judge('R1', 'the coupler\'s Lv of IGLOO', ri[0] if ri else None, lv_ref, TOL8)
    judge('R1', 'the coupler\'s Lv of ICE  ', rc[0] if rc else None, lv_ref, TOL8)
    judge('R2', 'IGLOO Lv(boiling-temperature %.2f K)' % Tb_i,
          num(r'\[material 1\] vapour data tabulated.*?Lv\(boiling-temperature\) =\s*(\S+)', ti), lv_bi, TOL5)
    judge('R2', 'ICE   Lv(boiling-temperature %.2f K)' % Tb_c, num(r'Lv\(boiling-temperature\) = (\S+) J/kg, Mv', tc),
          lv_bc, TOL5)
    judge('R3', 'IGLOO Psat(boiling-temperature)/boiling-pressure',
          num(r'\[material 1\] p_sat tabulated.*?boiling-pressure =\s*(\S+)', ti), at(tab['Psat'], Tb_i) / PATM,
          TOL_RATIO)
    judge('R3', 'ICE   Psat(boiling-temperature)/1 atm', num(r'Psat from the table, psat\(boiling-temperature\) = (\S+) atm',
                                                         tc), at(tab['Psat'], Tb_c) / PATM, TOL5)
    hoff = tab['Enthalpy_abs'][Tmin] - tab['Cp'][Tmin] * Tmin
    judge('R4', 'IGLOO hOff (the coupler\'s line)', num(r'\[MI2\] IGLOO material 1: table datum absolute, hOff =\s*(\S+)',
                                                        ti), hoff, TOL8)
    judge('R4', 'ICE   hOff', num(r'material 1: density .*?enthalpy absolute \(hOff = (\S+) J/kg\)', tc), hoff, TOL5)
    judge('R5', 'ICE   Mv', num(r'Mv = (\S+) kg/kmol', tc), tab['Mv'][Tmin], TOL5)
    judge('R5', 'ICE   Cpv(boiling-temperature)', num(r'Cpv\(boiling-temperature\) = (\S+) J/\(kg K\)', tc),
          at(tab['Cpv'], Tb_c), TOL5)
    (ok if took_i else ko)('R5 IGLOO took the Mv and Cpv columns (its log names them, without their values)')
    if ri and rc:
        d = rc[2] - ri[2]
        e = abs(d - hoff) / abs(hoff)
        (ok if e <= TOL8 else ko)('R6 the coupler\'s datum offsets, ICE %.8e and IGLOO %.8e, differ by the table\'s hOff '
                                  '(ICE credits cp T, IGLOO the absolute enthalpy, both the same Lv): %.8e vs %.8e, '
                                  'rel.err %.1e (tol %.0e)' % (rc[2], ri[2], d, hoff, e, TOL8))
    else:
        ko('R6 the coupler\'s datum block is missing in %s' % ('the IGLOO log' if not ri else 'the ICE log'))

if __name__ == '__main__':
    if len(sys.argv) >= 3 and sys.argv[1] == 'static':
        static(sys.argv[2])
    elif len(sys.argv) >= 6 and sys.argv[1] == 'run':
        run(*sys.argv[2:6])
    else:
        print(__doc__)
        sys.exit(2)
    print('== p15: %d passed, %d failed' % (npass, nfail))
    sys.exit(1 if nfail else 0)
