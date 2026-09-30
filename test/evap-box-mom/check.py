#!/usr/bin/env python3
"""evap-box-mom: momentum and energy bookkeeping of IGLOO's source coupling, read from the MI2 coupling report.

  python3 check.py [RUN_DIR] [--log PATH] [--tol-mom 1e-8] [--tol-en 1e-5]

RUN_DIR holds input.ini, INPUT/ and OUTPUT/ of a run (default: the current directory); --log is the run's log
(default RUN_DIR/run.log). One PASS/FAIL line per check, with the worst step, its numbers and the tolerance;
exit 1 on any FAIL. Text only (no numpy).

The report ([MI2-Coupling] report = true) gives, per coupling step, the MOSE block [t0, t1], the gas momentum
p = SUM rho u vol and total energy E = SUM rhoE vol at both ends, and, once IGLOO has produced sources, the totals
the coupler applies during the block: F = SUM sourceMom [N] and Q = SUM sourceEn [W]. IGLOOxMOSE adds these
absolute rates to MOSE's residual, unchanged over the block (IGLOO solves between blocks, relaxation off), and the
box is periodic in x and y with symmetry planes in z, so no boundary force acts on the gas:

  M1   p1 - p0 = F (t1 - t0), component by component: x to tol-mom relative, y and z to tol-mom of |Fx (t1 - t0)|
       (Fy and Fz are 1e-4 of Fx or zero); also summed over the run, and zero change over the sourceless block
  M3a  E1 - E0 = Q (t1 - t0) to tol-en relative. MOSE's energy update loses about 7e-12 J per iteration here
       (measured: per iteration, absent without the heat exchange, absent in momentum; consistent with the
       conservative energy being rebuilt at every RK stage from the pressure of FLINT's cons2prim, a Newton
       iteration on T that exits at 1e-6): 3 to 4e-6 of |Q (t1 - t0)| over a 50-iteration block, the floor
  M2   Fx is the x-momentum the streams give up, SUM_streams mdot (u_inj - u_end):
       (a) against SUM mdot u_inj alone, to 1e-2 (u_end is the gas velocity, below 1e-3 m/s: a missing
           volume factor is 1e7, a sign is -1, a component swap is 0);
       (b) against IGLOO's own end state: u_end, v_end from the snapshot file of the sweep that produced the
           block's sources, to 1e-5 of SUM mdot u_inj (the snapshot prints velocities to 1e-6 m/s); Fy likewise
  M3b  Q = SUM_streams mdot [cp (T0 - T_end) + (u_inj^2 - |v_end|^2)/2], IGLOO's energy deposit of a parcel
       that keeps its mass (no latent heat, no datum), with T_end from the same snapshot, to 1e-6 relative
  NV   non-vacuity: |px1 - px0| at least half the momentum the streams inject over the block, Fx at least half
       SUM mdot u_inj, |Q| at least half SUM mdot cp (T_gas - T0), one snapshot row per stream at time-end,
       the blocks contiguous in time (t0 of a step is t1 of the previous one)
"""
import os
import re
import sys


def die(msg):
    print('FAIL  ' + msg)
    sys.exit(1)


def read_ini(path):
    """Minimal INI reader: {section: {key: value}}; ';' starts a comment; values keep their spacing."""
    ini, sec = {}, None
    for raw in open(path):
        line = raw.split(';', 1)[0].strip()
        if not line:
            continue
        if line.startswith('['):
            sec = line.strip('[]').strip()
            ini.setdefault(sec, {})
        elif '=' in line and sec is not None:
            k, v = line.split('=', 1)
            ini[sec][k.strip()] = v.strip()
    return ini


def values(ini, sec, key, n=None, default=None):
    """A key's numbers as a list; a single value is broadcast to n; default (broadcast) when absent."""
    s = ini.get(sec, {}).get(key)
    if s is None:
        if default is None:
            die('input.ini: [%s] %s is missing' % (sec, key))
        v = [float(default)]
    else:
        v = [float(x) for x in s.split()]
    if n is not None and len(v) == 1:
        v = v * n
    if n is not None and len(v) != n:
        die('input.ini: [%s] %s has %d values, %d streams' % (sec, key, len(v), n))
    return v


def read_report(path):
    """The [MI2-report] lines grouped per coupling step."""
    steps, cur = [], None
    for line in open(path, errors='replace'):
        if '[MI2-report]' not in line:
            continue
        f = line.split()
        tag = f[1]
        if tag == 'cycle':
            cur = {'cycle': int(f[2]), 'step': int(f[4]), 't0': float(f[6]), 't1': float(f[8])}
            steps.append(cur)
        elif cur is None:
            continue
        elif tag == 'gas-momentum':
            cur['p0'] = [float(x) for x in f[2:5]]
            cur['p1'] = [float(x) for x in f[5:8]]
        elif tag == 'gas-energy':
            cur['E0'], cur['E1'] = float(f[2]), float(f[3])
        elif tag == 'igloo-to-gas':
            cur['F'] = [float(x) for x in f[2:5]]
            cur['Q'], cur['M'] = float(f[5]), float(f[6])
    return steps


def read_snapshot(path):
    """Rows X Y Z U V W T d m ID t of a snapshot file, keyed by parcel ID."""
    rows = {}
    for line in open(path):
        f = line.split()
        if len(f) != 11:
            continue
        try:
            v = [float(x) for x in f]
        except ValueError:
            continue
        rows[int(v[9])] = v
    return rows


def main():
    args = sys.argv[1:]
    run, log, tol_mom, tol_en = '.', None, 1e-8, 1e-5
    i = 0
    while i < len(args):
        a = args[i]
        if a == '--log':
            log = args[i + 1]; i += 2
        elif a == '--tol-mom':
            tol_mom = float(args[i + 1]); i += 2
        elif a == '--tol-en':
            tol_en = float(args[i + 1]); i += 2
        else:
            run = a; i += 1
    log = log or os.path.join(run, 'run.log')
    for p in (os.path.join(run, 'input.ini'), log):
        if not os.path.isfile(p):
            die('missing file ' + p)

    # --- the case: streams, material, gas ---------------------------------------------------------------
    ini = read_ini(os.path.join(run, 'input.ini'))
    xs = values(ini, 'IGLOO-BC', 'x')
    ns = len(xs)
    mdot = values(ini, 'IGLOO-BC', 'mdot', ns)
    uinj = values(ini, 'IGLOO-BC', 'up', ns, 0.0)
    vinj = values(ini, 'IGLOO-BC', 'vp', ns, 0.0)
    winj = values(ini, 'IGLOO-BC', 'wp', ns, 0.0)
    T0 = values(ini, 'IGLOO-BC', 'temp0', ns)
    tend = values(ini, 'IGLOO-General', 'time-end')[0]
    phase = ini.get('IGLOO-General', {}).get('phase', '')
    Tgas = values(ini, 'ICB-Block1', 'T')[0]
    cp = None
    for sec, kv in ini.items():
        if sec.startswith('GPB-Phase') and kv.get('name') == phase:
            cp = float(kv['cp'])
    if cp is None:
        die('input.ini: no [GPB-Phase*] section named %r with a cp' % phase)
    prefix = phase + '-' if phase else ''
    with open(os.path.join(run, 'INPUT', prefix + 'phase.txt')) as fh:
        mat = fh.read().split('\n')[1].split()[0]
    P_in = sum(mdot[s] * uinj[s] for s in range(ns))             # x-momentum injected per second [N]
    Q_scale = sum(mdot[s] * cp * abs(Tgas - T0[s]) for s in range(ns))   # heat to thermalise the streams [W]
    print('  case: %d streams, SUM mdot %.6e kg/s, SUM mdot u_inj %.6e N, cp %.1f, T0 %.1f K, gas %.1f K, '
          'time-end %.4e s, material %s' % (ns, sum(mdot), P_in, cp, T0[0], Tgas, tend, mat))

    # --- the report --------------------------------------------------------------------------------------
    steps = read_report(log)
    if not steps:
        die('no [MI2-report] lines in ' + log + ' (binary without the report, or report = false)')
    gated = [s for s in steps if 'F' in s]
    if not gated:
        die('no coupling step with IGLOO sources in the report (%d steps)' % len(steps))
    for s in steps:
        for k in ('p0', 'p1', 'E0', 'E1'):
            if k not in s:
                die('step %d: incomplete report lines' % s['step'])

    # --- per-step table (information) ------------------------------------------------------------------
    print('  step      dt [s]        dpx [kg m/s]     Fx dt [kg m/s]   dpx/(Fx dt)-1    dpy-Fy dt /|Fx dt|'
          '   dE [J]          Q dt [J]        dE/(Q dt)-1')
    for s in steps:
        dt = s['t1'] - s['t0']
        dp = [s['p1'][c] - s['p0'][c] for c in range(3)]
        dE = s['E1'] - s['E0']
        if 'F' not in s:
            print('  %4d  %.6e  %+.6e  (no sources)     -                -                 %+.6e  -' %
                  (s['step'], dt, dp[0], dE))
            continue
        Fdt = [s['F'][c] * dt for c in range(3)]
        Qdt = s['Q'] * dt
        print('  %4d  %.6e  %+.10e  %+.10e  %+.3e  %+.3e  %+.10e  %+.10e  %+.3e' %
              (s['step'], dt, dp[0], Fdt[0], dp[0] / Fdt[0] - 1 if Fdt[0] else float('nan'),
               (dp[1] - Fdt[1]) / abs(Fdt[0]) if Fdt[0] else float('nan'), dE, Qdt,
               dE / Qdt - 1 if Qdt else float('nan')))

    results = []

    def check(ok, what):
        results.append(ok)
        print(('PASS  ' if ok else 'FAIL  ') + what)

    def worst(items):
        """(value, step) with the largest value."""
        return max(items, key=lambda t: t[0])

    # --- M1: momentum bookkeeping per block --------------------------------------------------------------
    for c, name in enumerate('xyz'):
        errs = []
        for s in gated:
            dt = s['t1'] - s['t0']
            Fxdt = abs(s['F'][0] * dt)
            err = abs((s['p1'][c] - s['p0'][c]) - s['F'][c] * dt)
            errs.append((err / Fxdt if Fxdt > 0 else float('inf'), s['step'], err, s['F'][c] * dt))
        e, st, ea, ref = worst(errs)
        if c == 0:
            check(e <= tol_mom, 'M1x  px1-px0 = Fx (t1-t0) on %d blocks: worst rel %.3e at step %d '
                  '(|err| %.3e on %.6e kg m/s; tol %.0e)' % (len(gated), e, st, ea, ref, tol_mom))
        else:
            check(e <= tol_mom, 'M1%s  p%s1-p%s0 = F%s (t1-t0) to %.0e of |Fx (t1-t0)|: worst %.3e at step %d '
                  '(|err| %.3e, F%s dt %+.3e)' % (name, name, name, name, tol_mom, e, st, ea, name, ref))
    tot_dp = gated[-1]['p1'][0] - gated[0]['p0'][0]
    tot_F = sum(s['F'][0] * (s['t1'] - s['t0']) for s in gated)
    e = abs(tot_dp - tot_F) / abs(tot_F) if tot_F else float('inf')
    check(e <= tol_mom, 'M1s  px over the gated blocks: %.10e vs SUM Fx dt %.10e, rel %.3e (tol %.0e)' %
          (tot_dp, tot_F, e, tol_mom))
    ref = abs(gated[0]['F'][0] * (gated[0]['t1'] - gated[0]['t0']))
    nosrc = [s for s in steps if 'F' not in s]
    if nosrc:
        e = max(max(abs(s['p1'][c] - s['p0'][c]) for c in range(3)) / ref for s in nosrc)
        check(e <= tol_mom, 'M1o  no momentum change over the %d sourceless block(s): worst %.3e of the first '
              'gated |Fx dt| (tol %.0e)' % (len(nosrc), e, tol_mom))
    gaps = [(abs(steps[k + 1]['t0'] - steps[k]['t1']), steps[k + 1]['step']) for k in range(len(steps) - 1)]
    if gaps:
        g, st = worst(gaps)
        check(g == 0.0, 'NVt  blocks contiguous in time: worst t0(s) - t1(s-1) = %.3e at step %d (must be 0)' %
              (g, st))

    # --- M3a: energy bookkeeping per block ----------------------------------------------------------------
    errs = []
    for s in gated:
        dt = s['t1'] - s['t0']
        Qdt = s['Q'] * dt
        err = abs((s['E1'] - s['E0']) - Qdt)
        errs.append((err / abs(Qdt) if Qdt else float('inf'), s['step'], err, Qdt))
    e, st, ea, ref = worst(errs)
    check(e <= tol_en, 'M3a  E1-E0 = Q (t1-t0) on %d blocks: worst rel %.3e at step %d (|err| %.3e J on %.6e J; '
          'tol %.0e)' % (len(gated), e, st, ea, ref, tol_en))

    # --- M2/M3b: the deposit against the streams and IGLOO's own end states ------------------------------
    e, st = worst([(abs(s['F'][0] - P_in) / P_in, s['step']) for s in gated])
    check(e <= 1e-2, 'M2a  Fx = SUM mdot u_inj (%.6e N): worst rel %.3e at step %d (tol 1e-2)' % (P_in, e, st))
    nv_F = min((s['F'][0] / P_in, s['step']) for s in gated)
    check(nv_F[0] >= 0.5, 'NVF  Fx at least half SUM mdot u_inj: min Fx/(SUM mdot u_inj) %.4f at step %d' % nv_F)

    ex, ey, eq, snaps_ok, missing = [], [], [], True, []
    for s in gated:
        k = s['step'] - 2                                  # the solve of step s-1 wrote sweep tag k
        tag = '' if k == 0 else '-sweep%d' % k
        path = os.path.join(run, 'OUTPUT', '%ssnapshot-%s%s.dat' % (prefix, mat, tag))
        if not os.path.isfile(path):
            missing.append(path); snaps_ok = False
            continue
        rows = read_snapshot(path)
        if sorted(rows) != list(range(1, ns + 1)) or any(abs(r[10] - tend) > 1e-9 * tend for r in rows.values()):
            snaps_ok = False
            print('  snapshot %s: parcels %s, times %s' % (os.path.basename(path), sorted(rows),
                                                            sorted(set(r[10] for r in rows.values()))))
        Fx_ref = sum(mdot[s_ - 1] * (uinj[s_ - 1] - rows[s_][3]) for s_ in rows)
        Fy_ref = sum(mdot[s_ - 1] * (vinj[s_ - 1] - rows[s_][4]) for s_ in rows)
        Q_ref = sum(mdot[s_ - 1] * (cp * (T0[s_ - 1] - rows[s_][6])
                                    + 0.5 * (uinj[s_ - 1] ** 2 + vinj[s_ - 1] ** 2 + winj[s_ - 1] ** 2
                                             - rows[s_][3] ** 2 - rows[s_][4] ** 2 - rows[s_][5] ** 2))
                    for s_ in rows)
        ex.append((abs(s['F'][0] - Fx_ref) / P_in, s['step'], s['F'][0], Fx_ref, tag or '-sweep0'))
        ey.append((abs(s['F'][1] - Fy_ref) / P_in, s['step'], s['F'][1], Fy_ref, tag or '-sweep0'))
        eq.append((abs(s['Q'] - Q_ref) / abs(Q_ref), s['step'], s['Q'], Q_ref, tag or '-sweep0'))
    check(snaps_ok, 'NVs  one snapshot row per stream at time-end for every producing sweep (%d sweeps%s)' %
          (len(gated), '; missing ' + ', '.join(missing) if missing else ''))
    if ex:
        e, st, a, b, tg = worst(ex)
        check(e <= 1e-5, 'M2b  Fx = SUM mdot (u_inj - u_end), u_end from the snapshot: worst %.3e of SUM mdot '
              'u_inj at step %d (%.10e vs %.10e N, %s; tol 1e-5)' % (e, st, a, b, tg))
        e, st, a, b, tg = worst(ey)
        check(e <= 1e-5, 'M2c  Fy = SUM mdot (v_inj - v_end), v_end from the snapshot: worst %.3e of SUM mdot '
              'u_inj at step %d (%.6e vs %.6e N, %s; tol 1e-5)' % (e, st, a, b, tg))
        e, st, a, b, tg = worst(eq)
        check(e <= 1e-6, 'M3b  Q = SUM mdot [cp (T0 - T_end) + (u_inj^2 - |v_end|^2)/2], T_end and v_end from '
              'the snapshot: worst rel %.3e at step %d (%.10e vs %.10e W, %s; tol 1e-6)' % (e, st, a, b, tg))

    # --- NV: the effect is large enough to be measured ----------------------------------------------------
    r = min(((s['p1'][0] - s['p0'][0]) / (P_in * (s['t1'] - s['t0'])), s['step']) for s in gated)
    check(r[0] >= 0.5, 'NVp  |px1-px0| at least half the momentum the streams inject over the block: min ratio '
          '%.4e at step %d' % r)
    r = min((abs(s['Q']) / Q_scale, s['step']) for s in gated)
    check(r[0] >= 0.5, 'NVq  |Q| at least half SUM mdot cp |T_gas - T0| (%.4e W): min ratio %.4f at step %d' %
          (Q_scale, r[0], r[1]))

    n_ok = sum(1 for x in results if x)
    print('  %d checks passed, %d failed' % (n_ok, len(results) - n_ok))
    sys.exit(0 if all(results) else 1)


if __name__ == '__main__':
    main()
