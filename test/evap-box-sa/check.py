#!/usr/bin/env python3
"""Checks of the evap-box-sa gate on one run: the closed-domain mass balance (R1), the energy
bookkeeping (R2), the Spalart-Allmaras field (R3) and the log hygiene with the RANS setup (R4).
usage: check.py <run.log> <OUTPUT dir> <INPUT dir>   -- prints PASS/FAIL lines, exits 1 on any failure"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mi2-checks'))
from tec_cc import read_cc  # noqa: E402

# tolerances and non-vacuity bounds
MASS_RELERR = 1e-12     # species-1 row of the mass balance (evap-box-ta closes to 2.2e-16)
MASS_MIN    = 5e-9      # kg injected over the run (evap-box-ta injects 1.495e-8)
EN_RELERR   = 1e-7      # energy bookkeeping (evap-box-ta: 2.8e-9)
EN_MIN      = 1e-2      # J handed to the gas (evap-box-ta: 0.236)
SA_DEV      = 5e-2      # |mi_tilde/mi_tilde0 - 1| everywhere: production and the injection-cell density change
                        # move it by at most 3.5 % here (measured on the reference run: +1.2 % / -0.16 %)
MUT_RATIO   = (0.1, 0.4)  # mit/mil: the far-field SA state chi about 3 gives 0.21 (reference run 0.199-0.221)
BAD = ['[ERROR]', ' FAIL', 'Non-physical state', 'Integration failed at', 'NaN or p<0 detected', 'ABORT:']

npass = nfail = 0


def ok(msg):
    global npass
    npass += 1
    print('PASS  ' + msg)


def ko(msg):
    global nfail
    nfail += 1
    print('FAIL  ' + msg)


def check(cond, msg):
    (ok if cond else ko)(msg)


def read_ic_last_var(path):
    """Last cell-centred variable of an ATLAS ICB Tecplot BLOCK file (its solver names are unquoted)."""
    txt = open(path).read()
    head, zone, body = txt.split('\n', 2)
    names = re.sub(r'"', ' ', head.split('=', 1)[1]).split()
    I, J, K = (int(re.search(r'\b%s\s*=\s*(\d+)' % a, zone).group(1)) for a in 'IJK')
    vals = np.array(body.split(), dtype=float)
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    nvar = len(names) - 3
    if vals.size != 3 * nn + nvar * nc:
        raise ValueError('%s: %d values for %d variables' % (path, vals.size, len(names)))
    return names[-1], vals[3 * nn + (nvar - 1) * nc:]


log, out, inp = sys.argv[1:4]
txt = open(log, errors='replace').read()

# R4a: hygiene
hits = [b for b in BAD if b in txt]
exit0 = re.search(r'^# exit=0\b', txt, re.M) is not None
check(exit0 and 'GAME OVER' in txt and not hits,
      'log hygiene: exit 0, GAME OVER, none of %s (exit0=%s, found=%s)' % (BAD, exit0, hits))
# R4b: the run is RANS
ns = re.search(r'Equations\s+Navier-Stokes', txt) is not None
sa = re.search(r'Viscous model\s+Turbulent, SA\b', txt) is not None
check(ns and sa, 'MOSE runs Navier-Stokes with the SA closure (Navier-Stokes=%s, Turbulent SA=%s)' % (ns, sa))

# R1: closed-domain mass balance, species 1
m = re.search(r'closed-domain mass balance.*?\n.*?\n\s+1\s+(\S+)\s+(\S+)\s+(\S+)', txt, re.S)
if m is None:
    ko('closed-domain mass balance row of species 1 not found in the log')
else:
    inj, gained, rel = (float(v) for v in m.groups())
    check(rel <= MASS_RELERR, 'mass balance species 1: injected %.6e kg, gained %.6e kg, rel.err %.3e <= %.0e'
          % (inj, gained, rel, MASS_RELERR))
    check(inj >= MASS_MIN, 'mass injected %.6e kg >= %.0e kg (non-vacuity)' % (inj, MASS_MIN))

# R2: energy bookkeeping
h = re.search(r'handed to gas\s+\(sourceEn \+ datum\) \[J\] :\s+(\S+)', txt)
e = re.search(r'rel\.err \(bookkeeping only, datum-blind\) :\s+(\S+)', txt)
if h is None or e is None:
    ko('energy bookkeeping lines not found in the log')
else:
    handed, rel = float(h.group(1)), float(e.group(1))
    check(rel <= EN_RELERR, 'energy bookkeeping: handed %.8e J, rel.err %.3e <= %.0e' % (handed, rel, EN_RELERR))
    check(abs(handed) >= EN_MIN, 'energy handed to the gas |%.4e| J >= %.0e J (non-vacuity)' % (handed, EN_MIN))

# R3: the SA field of the last written gas field
field = os.path.join(out, 'gas-field.tec')
if not os.path.isfile(field):
    ko('no gas field written: %s missing' % field)
else:
    head = open(field).readline()
    t = re.search(r'SOLUTIONTIME=([^,\s]+)', open(field).read())
    tend = re.findall(r'^MOSE \| Iter =\s*(\d+) \| Time = (\S+)', txt, re.M)
    print('  field %s at t=%s (last MOSE line: iter %s, t=%s)'
          % (field, t.group(1) if t else '?', tend[-1][0] if tend else '?', tend[-1][1] if tend else '?'))
    check('mi_tilde' in re.sub(r'[",]', ' ', head.split('=', 1)[1]).split(), 'the gas field carries the SA variable mi_tilde (nprim > np)')
    d = read_cc(field)
    cellvars = [k for k in d if k not in ('x', 'y', 'z')]
    nonfinite = [k for k in cellvars if not np.isfinite(d[k]).all()]
    check(not nonfinite, 'every field variable finite everywhere (non-finite: %s)' % nonfinite)
    if 'mi_tilde' in d and not nonfinite:
        name0, mt0 = read_ic_last_var(os.path.join(inp, 'gas-ic.tec'))
        mt0 = mt0.reshape(d['mi_tilde'].shape)
        dev = d['mi_tilde'] / mt0 - 1.0
        k, j, i = np.unravel_index(np.argmax(np.abs(dev)), dev.shape)
        check(name0 == 'mi_t' and np.abs(dev).max() <= SA_DEV,
              'SA variable within %.0f %% of its initial value everywhere: min %.3e, max %.3e (largest at cell i=%d j=%d k=%d)'
              % (100 * SA_DEV, dev.min(), dev.max(), i + 1, j + 1, k + 1))
        if 'mit' in d and 'mil' in d:
            r = d['mit'] / d['mil']
            check(MUT_RATIO[0] <= r.min() and r.max() <= MUT_RATIO[1],
                  'eddy viscosity ratio mit/mil in [%.2f, %.2f]: min %.4f, max %.4f' % (MUT_RATIO + (r.min(), r.max())))
        else:
            ko('mit/mil not in the field (sol-variables must hold transport)')

print('== %d passed, %d failed' % (npass, nfail))
sys.exit(0 if nfail == 0 else 1)
