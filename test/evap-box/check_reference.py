#!/usr/bin/env python3
"""Independent reference for hydra's test/evap-box: the droplet ODE IGLOO integrates, solved again.

The case (REFERENCE.md): four DB streams of 20 um water drops at 320 K and 0.5 m/s into quiescent N2 at
350 K and 101325 Pa, in a closed 0.1 x 0.1 x 0.005 m box; Stokes drag, Ranz-Marshall heat, CEM evaporation
with the Langmuir-Knudsen interface, p_sat from Clausius-Clapeyron or the material's Psat column. Each drop
stops within about half a millimetre and evaporates in place to burnout. Per stream this script integrates

    dx/dt = u,  du/dt = 3 pi mu_g d (u_g - u)/m                                          (Stokes, [CGW78])
    m c_l dT/dt = pi d k_g Nu (T_g - T) + mdot L_v,   Nu = 2 + 0.6 Re^1/2 Pr^1/3                    ([RM52])
    dm/dt = mdot = -pi d k_g/(c_pg Le) Sh ln(1 + B_M),  Sh = 2 + 0.6 Re^1/2 Sc^1/3,  B_M = (Y_s - Y_inf)/(1 - Y_s)
    X_s = p_sat(T)/p - (2 L_K/d) beta,  beta = -mdot Pr/(2 pi mu_g d)                     (CEM + LK, [MHB98])

with its own adaptive Dormand-Prince 5(4) integrator (rtol 1e-11) and an exact solve of the LK fixed point,
in the gas IGLOO sampled: the initial condition for sweep 0 (uniform), OUTPUT/gas-field.tec for the last
sweep (cell-wise constant, gas-order 1, faces crossed as IGLOO crosses them). Gas transport and
thermodynamics are rebuilt from INPUT/gas-*.dat as the gas solver builds them (integer-kelvin tables, Wilke
mixing). Then, per sweep:

  E  per stream: exit position (the stop point); exit T and, when out-time prints it, exit time at burnout
     (IGLOO consumes a drop at m <= min(1e-6 m0, 1e-20 kg)), paired with the exit mass when a final trajectory
     row (print-dtime) or the exit row's diameter carries it, else bracketed by the reference over m in
     [0.9, 1] x the burnout mass; the number of scatter samples, one at every quantum up to the exit
  S  every scatter sample, paired with the reference state at the same droplet mass: x, y, u, v, T within
     their print budget (the E13.6 mass rounding times d(.)/dm, plus the F12.6 rounding) plus what a sample on
     a solver step's chord, t from its quantum (t within the C tolerance), carries: t (|m' y''/m''| + |y'|), and
     t |y'| on x, y (the position is cubic on the step); on the plateau, for a few thermal relaxation times after
     the drop temperature crosses a node of the Psat column, plus what a solver step across the node can leave in
     T (the change of slope of the table over that step, relaxing after)
  Q  plateau samples against the closed-form wet-bulb balance at their diameter, beyond the reference's lag
  C  the sample clock: IGLOO places the k-th sample of stream n at t_k = dt_q (k - frac(n phi)),
     dt_q = L_box/(v_inj fsample-traj), inside the solver step that spans it (the columns linear on the step,
     the position cubic); the reference reaches the sample's mass within 1e-6 s of t_k
  G  last sweep: the source field against the injected mass, momentum and energy of the consumed parcels; the
     Eulerian field's sum n_p V and sum rho_p V against sum npdot t_life and sum npdot int m dt; the gas mass
     the coupler reports gained per species (run log) against the final gas field minus the initial one
  L  the run log, every sweep: no line of IGLOO's LK-interface: warning (a parcel whose Langmuir-Knudsen fixed
     point was not bracketed carries the bracket's estimate of its rate, not the model's)

The last sweep's per-stream checks and the Eulerian integrals are gated only when OUTPUT/gas-field.tec was
written at the MOSE iteration that sweep sampled (the last line of gas-residual-history.dat); otherwise they
are printed as information.

Usage: check_reference.py [RUN_DIR] [--ini PATH] [--props PATH] [--log PATH] [--sweeps 0,last]
                          [--psat auto|cc|table] [--interface auto|LK|VLE]
--ini and --props fix the modelled case independently of the run (the falsification runs). Exit 0 when every
gated check passes. Python 3, numpy only.
"""
import argparse
import copy
import math
import os
import re
import sys

import numpy as np

RU_IGLOO = 8314.46        # IGLOO's universal gas constant (Lib_Evaporation, IO)
RU_FLINT = 8314.46261815324   # the gas solver's (FLINT Lib_ThermoTransport: the exact SI value)
PATM = 101325.0
XS_CAP = 1.0e-12
EPS_BURN, M_BURN_ABS = 1.0e-6, 1.0e-20   # IGLOO consumes a drop at m <= min(EPS_BURN m0, M_BURN_ABS) [kg]
M_TAIL = 1.0e-3           # the reference integrates down to M_TAIL x the burnout mass
TRAJ_SAMPLE = 100         # IGLOO's default of [IGLOO-General] fsample-traj: scatter points per stream over L_box/v_inj
PHI = 0.6180339887498949  # scatter phase offset per ID

# ------------------------------------------------------------------ tolerances (measured floors, REFERENCE.md)
TOL = {
    "ratio": 2.0,          # S: residual / print budget, every sample and column (measured 0.99)
    "ratio_cells": 10.0,   # S on a cell-wise gas: plus the face-crossing switch time (measured 5.1)
    "T_qs": 1.0e-3,        # Q: K beyond the reference's own quasi-steady lag
    "t_lead": 1.0e-6,      # C: s a sample may precede its quantum time
    "t_lag": 1.0e-6,       # C: s a sample may trail it (IGLOO places it at the quantum, inside its step)
    "t_exit": 1.0e-6,      # E: relative exit time, when out-time prints it
    "src_rel": 1.0e-9,     # G: source-field sums vs the injected fluxes of consumed parcels
    "eul_rel": 1.0e-4,     # G: Eulerian residence and mass integrals vs the reference
}


# ------------------------------------------------------------------ input readers
def read_ini(path):
    sec, out = None, {}
    for raw in open(path):
        line = raw.split(";")[0].strip()
        if not line:
            continue
        m = re.match(r"\[(.+)\]$", line)
        if m:
            sec = m.group(1).strip()
            out.setdefault(sec, {})
            continue
        if "=" in line and sec is not None:
            k, v = line.split("=", 1)
            out[sec][k.strip()] = v.strip()
    return out


def floats(s):
    return [float(t) for t in s.split()]


def read_columns_table(path):
    """Tecplot POINT table: returns (variable names, zones) with each zone a 2-D array."""
    names, zones, cur = [], [], []
    for raw in open(path):
        line = raw.strip()
        if not line:
            continue
        if line.upper().startswith("VARIABLES"):
            names = [t.strip('"') for t in re.findall(r'"[^"]+"|[^\s,]+', line.split("=", 1)[1])]
            continue
        if line.upper().startswith("ZONE"):
            if cur:
                zones.append(np.array(cur))
            cur = []
            continue
        if line.upper().startswith(("TITLE", "I=")):
            continue
        try:
            cur.append([float(t) for t in line.split()])
        except ValueError:
            continue
    if cur:
        zones.append(np.array(cur))
    return names, zones


def read_tec_block(path):
    """BLOCK Tecplot (3 nodal coordinates + cell-centred variables): returns names and a dict of arrays."""
    lines = open(path).read().splitlines()
    names = re.findall(r'"([^"]+)"|(\S+)', lines[0].split("=", 1)[1])
    names = [a or b for a, b in names]
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    zone = lines[zi]
    ni, nj, nk = (int(re.search(r"\b%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1)
    # the zone's values end at the auxiliary data a solver appends after them
    zend = next((i for i in range(zi + 1, len(lines)) if lines[i].lstrip().startswith("DATASETAUXDATA")), len(lines))
    vals = np.array(" ".join(lines[zi + 1:zend]).split(), dtype=float)
    out, pos = {}, 0
    for n, name in enumerate(names):
        cnt = nn if n < 3 else nc
        arr = vals[pos:pos + cnt]
        pos += cnt
        shape = (nk, nj, ni) if n < 3 else (max(nk - 1, 1), nj - 1, ni - 1)
        out[name] = arr.reshape(shape)
    return names, out, (ni, nj, nk)


def read_particles(path):
    """IGLOO trajectory, scatter or outloc rows by ID, file order kept, and {name: column} of the file's variables=
    line: the ID (and the exit time t) are found by name, so further columns (bc, d_p, weight) do not matter."""
    rows, cols = {}, {}
    if not os.path.exists(path):
        return rows, cols
    for raw in open(path):
        if raw.lstrip().lower().startswith("variables"):
            cols = {n: i for i, n in enumerate(re.findall(r'"([^"]*)"', raw))}
            continue
        c = raw.split()
        if len(c) < 9 or not c[0][0] in "0123456789-":
            continue
        try:
            v = [float(t) for t in c]
        except ValueError:
            continue
        idcol = cols.get("ID", 9 if len(c) in (10, 11) and "outloc" not in os.path.basename(path) else 8)
        rows.setdefault(int(v[idcol]), []).append(v)
    return rows, cols


# ------------------------------------------------------------------ gas properties (FLINT's evaluation)
class GasTables:
    def __init__(self, inp):
        ph = [l.split() for l in open(os.path.join(inp, "gas-phase.txt")).read().splitlines()[1:] if l.strip()]
        self.species = [p[0] for p in ph]
        self.wm = np.array([float(p[1]) for p in ph])
        self.ri = RU_FLINT / self.wm
        _, tz = read_columns_table(os.path.join(inp, "gas-thermo.dat"))
        _, rz = read_columns_table(os.path.join(inp, "gas-transport.dat"))
        self.tmin = int(round(tz[0][0, 0]))
        self.cp = np.array([z[:, 1] for z in tz])
        self.mu = np.array([z[:, 1] for z in rz])
        self.k = np.array([z[:, 2] for z in rz])
        self.tmax = self.tmin + self.cp.shape[1] - 1
        ns = len(self.wm)
        self.mm025 = np.array([[(self.wm[i] / self.wm[j]) ** -0.25 for j in range(ns)] for i in range(ns)])
        self.isq8 = np.array([[1.0 / math.sqrt(8.0 * (1.0 + self.wm[i] / self.wm[j])) for j in range(ns)]
                              for i in range(ns)])

    def _interp(self, tab, T):
        ti = int(T)
        td = T - ti
        i1 = max(self.tmin, min(ti, self.tmax)) - self.tmin
        i2 = max(self.tmin, min(ti + 1, self.tmax)) - self.tmin
        return tab[:, i1] + (tab[:, i2] - tab[:, i1]) * td

    def state(self, rhoi, u, v, w, p):
        """rho, u, v, w, T, mu, k, gamma, R from the partial densities and pressure, as MI2 hands them over."""
        rhoi = np.asarray(rhoi, dtype=float)
        rho = rhoi.sum()
        R = float((self.ri * rhoi / rho).sum())
        T = p / (rho * R)
        cp = float((rhoi / rho * self._interp(self.cp, T)).sum())
        gam = cp / (cp - R)
        mui, ki = self._interp(self.mu, T), self._interp(self.k, T)
        wmt = 1.0 / float((rhoi / rho / self.wm).sum())
        xi = rhoi * wmt / (rho * self.wm)
        ns = len(xi)
        fi = np.empty((ns, ns))
        for i in range(ns):
            for j in range(ns):
                ratio = math.sqrt(mui[i]) / math.sqrt(mui[j] + 1e-20) * self.mm025[i, j]
                fi[i, j] = (1.0 + ratio) ** 2 * self.isq8[i, j]
        den = fi @ xi
        mu = float((xi * mui / den).sum())
        k = float((xi * ki / den).sum())
        return (rho, u, v, w, T, mu, k, gam, R)


# ------------------------------------------------------------------ the droplet model
class Droplet:
    def __init__(self, ini, phase_tokens, props, psat_mode, interface):
        prop = ini.get("IGLOO-Properties", {})
        names, zones = props
        tab = zones[0]
        col = {n: i for i, n in enumerate(names)}
        self.tmin = int(round(tab[0, 0]))
        # L_v and M_v as IGLOO takes them: the [IGLOO-Properties] value, else the table's Lv and Mv columns when their
        # first row is non-zero; a value given both ways is refused, as IGLOO refuses it. L_v(T) at the droplet
        # temperature, its value at 298.15 K in the Clausius-Clapeyron constant
        self.lv_tab = None
        ini_v = {key: float(prop.get(key, "0")) for key in ("Lv", "Mv")}
        for key in ("Lv", "Mv"):
            if ini_v[key] != 0.0 and key in col and np.any(tab[:, col[key]] != 0.0):
                raise SystemExit("[FAIL] %s given both in [IGLOO-Properties] and as a column of the properties table "
                                 "(IGLOO refuses it)" % key)
        if ini_v["Lv"] != 0.0:
            self.Lv = ini_v["Lv"]
        elif "Lv" in col and tab[0, col["Lv"]] != 0.0:
            self.lv_tab = tab[:, col["Lv"]]
            self.Lv = self.table_at(self.lv_tab, 298.15)
        else:
            raise SystemExit("[FAIL] no latent heat: neither [IGLOO-Properties] Lv nor an Lv column in the table")
        if ini_v["Mv"] != 0.0:
            self.Mv = ini_v["Mv"]
        elif "Mv" in col and tab[0, col["Mv"]] != 0.0:
            self.Mv = float(tab[0, col["Mv"]])
        else:
            raise SystemExit("[FAIL] no vapour molar mass: neither [IGLOO-Properties] Mv nor an Mv column in the table")
        tb = prop.get("boiling-temperature", prop.get("Tboil"))
        self.Tboil = float(tb)
        self.Le = float(prop.get("Le", "1.0"))
        self.Yinf = float(prop.get("Yinf", "0.0"))
        self.alphaE = float(phase_tokens.get("alpha-e", "1.0"))
        self.lvmv = self.Lv * self.Mv / RU_IGLOO
        self.rho_l = tab[0, col["Density"]]
        self.cp_l = tab[0, col["Cp"]]
        if np.ptp(tab[:, col["Density"]]) != 0.0 or np.ptp(tab[:, col["Cp"]]) != 0.0:
            raise SystemExit("[FAIL] the reference assumes constant liquid cp and density (the case's table)")
        hcol = col["Enthalpy_abs"] if "Enthalpy_abs" in col else col["Enthalpy"]
        self.hoff = tab[0, hcol] - self.cp_l * tab[0, 0]      # IGLOO: hOff = h(Tmin) - cp*Tmin
        self.psat_tab = tab[:, col["Psat"]] if "Psat" in col else None
        if psat_mode == "auto":
            psat_mode = "table" if self.psat_tab is not None and np.any(self.psat_tab > 0.0) else "cc"
        if psat_mode == "table" and self.psat_tab is None:
            raise SystemExit("[FAIL] --psat table but the properties table has no Psat column")
        self.psat_mode = psat_mode
        models = ini.get("IGLOO-Models", {})
        if interface == "auto":
            interface = models.get("interface", "VLE")
        self.lk = interface == "LK"
        for key, want in (("drag", "Stokes"), ("heat", "Ranz-Marshall"), ("evaporation", "CEM")):
            if models.get(key) != want:
                raise SystemExit("[FAIL] the reference implements %s = %s, the case sets %s" % (key, want, models.get(key)))
        if models.get("blowing", "none") != "none":
            raise SystemExit("[FAIL] the reference implements blowing = none")

    def table_at(self, tab, T):
        """A table column at T, linear between the integer-kelvin rows (IGLOO's tableValue)."""
        lo, hi = self.tmin, self.tmin + len(tab) - 1
        tc = min(max(T, float(lo)), float(hi))
        if tc >= hi:
            return tab[-1]
        i = int(tc)
        a = tab[i - lo]
        return a + (tab[i + 1 - lo] - a) * (tc - i)

    def lv(self, T):
        return self.Lv if self.lv_tab is None else self.table_at(self.lv_tab, T)

    def psat(self, T):
        if self.psat_mode == "cc":
            return PATM * math.exp(-self.lvmv * (1.0 / T - 1.0 / self.Tboil))
        return self.table_at(self.psat_tab, T)

    def mdot(self, gas, T, d, Re):
        rho_g, _, _, _, Tg, mu, k, gam, R = gas
        cpg = gam * R / (gam - 1.0)
        p = rho_g * R * Tg
        mg = RU_IGLOO / R
        ps = self.psat(T)
        xs_eq = 1.0 - XS_CAP if ps >= p else min(ps / p, 1.0 - XS_CAP)
        pr = mu * cpg / k
        sc = pr * self.Le
        sh = 2.0 + 0.6 * math.sqrt(max(Re, 0.0)) * sc ** (1.0 / 3.0)
        pref = -math.pi * d * (k / (cpg * self.Le)) * sh

        def cem(xs):
            ys = xs * self.Mv / (xs * self.Mv + (1.0 - xs) * mg)
            if ys <= self.Yinf:
                return 0.0
            bm = (ys - self.Yinf) / (1.0 - ys)
            return pref * math.log1p(bm) if bm > 0.0 else 0.0

        md_eq = cem(xs_eq)
        if not self.lk or md_eq == 0.0:
            return md_eq
        lkd2 = 2.0 * mu * math.sqrt(2.0 * math.pi * T * RU_IGLOO / self.Mv) / (self.alphaE * sc * p) / d
        coef = lkd2 * pr / (2.0 * math.pi * mu * d)          # X_s = X_s,eq + coef*mdot (mdot < 0)

        def g(md):
            return md - cem(max(xs_eq + coef * md, 0.0))

        lo, hi = md_eq, 0.0                                    # g(lo) <= 0 <= g(hi)
        glo, ghi = g(lo), g(hi)
        if glo >= 0.0:
            return lo
        for _ in range(200):                                   # Illinois regula falsi
            mid = hi - ghi * (hi - lo) / (ghi - glo)
            gm = g(mid)
            if abs(gm) <= 1e-15 * abs(mid) or hi - lo <= 1e-16 * abs(lo):
                return mid
            if (gm < 0.0) == (glo < 0.0):
                lo, glo = mid, gm
                ghi *= 0.5
            else:
                hi, ghi = mid, gm
                glo *= 0.5
        return 0.5 * (lo + hi)

    def continued(self, node, rising):
        """The drop with the line of the table interval it leaves at `node` continued over the next interval, in every
        column the reference reads at the drop temperature (Psat, and L_v when tabulated)."""
        other = copy.copy(self)
        i = node - self.tmin
        for name in ("psat_tab", "lv_tab"):
            tab = getattr(self, name)
            if tab is None:
                continue
            tab = tab.copy()
            if rising:
                tab[i + 1] = 2.0 * tab[i] - tab[i - 1]
            else:
                tab[i - 1] = 2.0 * tab[i] - tab[i + 1]
            setattr(other, name, tab)
        return other

    def h(self, T):
        return self.cp_l * T + self.hoff

    def diameter(self, m):
        return (6.0 * m / (math.pi * self.rho_l)) ** (1.0 / 3.0) if m > 0.0 else 0.0

    def mass(self, d):
        return math.pi / 6.0 * self.rho_l * d ** 3

    def rhs(self, s, gas):
        x, y, u, v, T, m = s
        rho_g, ug, vg, _, Tg, mu, k, gam, R = gas
        d = self.diameter(m)
        du, dv = ug - u, vg - v
        slip = math.hypot(du, dv)
        Re = rho_g * slip * d / mu
        fac = 3.0 * math.pi * mu * d * Re / (Re + 1e-20) / m    # Stokes: (pi/8) (24/(Re+toll)) d Re mu / m
        pr = mu * gam * R / ((gam - 1.0) * k)
        nu = 2.0 + 0.6 * math.sqrt(Re) * pr ** (1.0 / 3.0)
        md = self.mdot(gas, T, d, Re)
        q = nu * k * math.pi * d * (Tg - T)
        return np.array([u, v, fac * du, fac * dv, (q + md * self.lv(T)) / (m * self.cp_l), md])

    def t_quasi_steady(self, gas, d, u=0.0, v=0.0):
        """T at which pi d k Nu (Tg - T) + mdot(T) Lv(T) = 0 (the wet-bulb balance at diameter d, slip from u, v)."""
        rho_g, ug, vg, _, Tg, mu, k, gam, R = gas
        slip = math.hypot(ug - u, vg - v)
        Re = rho_g * slip * d / mu
        pr = mu * gam * R / ((gam - 1.0) * k)
        nu = 2.0 + 0.6 * math.sqrt(Re) * pr ** (1.0 / 3.0)

        def f(T):
            return nu * k * math.pi * d * (Tg - T) + self.mdot(gas, T, d, Re) * self.lv(T)

        lo, hi = 200.0, Tg
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0.0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-12:
                break
        return 0.5 * (lo + hi)


# ------------------------------------------------------------------ gas fields
class UniformGas:
    def __init__(self, state):
        self.st = state

    def cell_of(self, x, y, u=0.0, v=0.0):
        return (0, 0)

    def gas(self, cell):
        return self.st

    def bounds(self, cell):
        return (-1e300, 1e300, -1e300, 1e300)


class CellGas:
    """Piecewise-constant gas per cell (gas-order 1) on the uniform 2-D grid of gas-field.tec."""

    def __init__(self, tables, path):
        names, f, (ni, nj, nk) = read_tec_block(path)
        self.xn = f[names[0]][0, 0, :]
        self.yn = f[names[1]][0, :, 0]
        rn = [n for n in names if n.startswith("rho(")]
        self.states = {}
        for j in range(nj - 1):
            for i in range(ni - 1):
                rhoi = [f[n][0, j, i] for n in rn]
                self.states[(i, j)] = tables.state(rhoi, f["u"][0, j, i], f["v"][0, j, i], f["w"][0, j, i],
                                                   f["p"][0, j, i])
        self.ni, self.nj = ni - 1, nj - 1

    def cell_of(self, x, y, u=0.0, v=0.0):
        """Lowest-index cell containing (x, y), faces inclusive (IGLOO's full-sweep search order)."""
        i = int(np.searchsorted(self.xn, x, side="left")) - 1
        j = int(np.searchsorted(self.yn, y, side="left")) - 1
        return (min(max(i, 0), self.ni - 1), min(max(j, 0), self.nj - 1))

    def gas(self, cell):
        return self.states[cell]

    def bounds(self, cell):
        i, j = cell
        return (self.xn[i], self.xn[i + 1], self.yn[j], self.yn[j + 1])


# ------------------------------------------------------------------ Dormand-Prince 5(4)
A = [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9],
     [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
     [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
     [35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]]
C = [0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0]
E = [71 / 57600, 0.0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40]


def dp_step(f, y, h):
    ks = []
    for i in range(7):
        yi = y + h * sum((a * k for a, k in zip(A[i], ks)), np.zeros_like(y)) if i else y
        ks.append(f(yi))
    ynew = y + h * sum(a * k for a, k in zip(A[6], ks[:6]))
    err = h * sum(e * k for e, k in zip(E, ks))
    return ynew, err


def burnout_mass(m0):
    return min(EPS_BURN * m0, M_BURN_ABS)


def integrate(drop, field, s0, m_stop, rtol=1e-11, tmax=1.0):
    """Adaptive DP45 from s0 until m <= m_stop; the gas cell switches when the drop leaves it.
    Returns the list of (t, state) at every accepted step, the last one at m = m_stop (check_sweep stops at
    M_TAIL x the burnout mass: IGLOO's exit step can end below the burnout mass)."""
    atol = np.array([1e-13, 1e-13, 1e-12, 1e-12, 1e-9, 1e-26])
    t, s = 0.0, np.array(s0, dtype=float)
    cell = field.cell_of(s[0], s[1], s[2], s[3])
    out = [(t, s.copy(), cell)]
    h = 1e-7
    ncross = 0
    while t < tmax:
        gas = field.gas(cell)
        f = lambda z: drop.rhs(z, gas)
        ynew, err = dp_step(f, s, h)
        sc = atol + rtol * np.maximum(np.abs(s), np.abs(ynew))
        en = math.sqrt(float(np.mean((err / sc) ** 2)))
        if en > 1.0 or not np.all(np.isfinite(ynew)) or ynew[5] <= 0.0:
            h *= max(0.2, 0.9 * en ** -0.2) if np.isfinite(en) and en > 0 else 0.2
            continue
        x0, x1, y0, y1 = field.bounds(cell)
        if not (x0 <= ynew[0] <= x1 and y0 <= ynew[1] <= y1):
            # locate the face crossing: shrink the step onto it, then switch to the neighbour
            lo, hi = 0.0, h
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                ym, _ = dp_step(f, s, mid)
                if x0 <= ym[0] <= x1 and y0 <= ym[1] <= y1:
                    lo = mid
                else:
                    hi = mid
            ynew, _ = dp_step(f, s, hi)
            t += hi
            s = ynew
            ci, cj = cell
            ci += int(s[0] > x1) - int(s[0] < x0)
            cj += int(s[1] > y1) - int(s[1] < y0)
            cell = (ci, cj)
            ncross += 1
            if ncross > 100000:
                raise SystemExit("[FAIL] reference: more than 1e5 cell crossings (chattering on a face)")
            out.append((t, s.copy(), cell))
            continue
        if ynew[5] <= m_stop:
            lo, hi = 0.0, h                                    # land the last step on m = m_stop
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                ym, _ = dp_step(f, s, mid)
                if ym[5] > m_stop:
                    lo = mid
                else:
                    hi = mid
            ynew, _ = dp_step(f, s, hi)
            t += hi
            out.append((t, ynew.copy(), cell))
            return out
        t += h
        s = ynew
        out.append((t, s.copy(), cell))
        h *= min(5.0, 0.9 * en ** -0.2) if en > 0 else 5.0
    raise SystemExit("[FAIL] reference integration did not reach burnout")


# ------------------------------------------------------------------ pairing and comparisons
F_HALF = 0.5e-6           # half-ULP of IGLOO's F12.6 columns (x, y, z, u, v, w, T)
NODE_STEP = 1.0           # S: thermal relaxation times a solver step across a node of the Psat column can last
NODE_TIMES = 5.0          # S: thermal relaxation times after such a node over which the node allowance is followed


def e_half(v):
    """Half-ULP of an E13.6 value (0.dddddd E+ee): d, m, mdot."""
    if v == 0.0:
        return 0.0
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def state_at(drop, field, traj, i, dt):
    t0, s0, cell = traj[i]
    if dt <= 0.0:
        return s0
    gas = field.gas(cell)
    y, _ = dp_step(lambda z: drop.rhs(z, gas), s0, dt)
    return y


def second_derivative(drop, gas, s, h=1.0e-7):
    """d2(state)/dt2 at s, by a central difference of the right-hand side over +-h along the reference."""
    f = lambda z: drop.rhs(z, gas)
    sp, sm = dp_step(f, s, h)[0], dp_step(f, s, -h)[0]
    return (f(sp) - f(sm)) / (2.0 * h)


def pair_by_mass(drop, field, traj, m_target):
    """Reference time, state and gas cell at which m = m_target (m decreases along traj)."""
    ms = [p[1][5] for p in traj]
    if m_target >= ms[0]:
        return traj[0][0], traj[0][1], traj[0][2]
    if m_target <= ms[-1]:
        return traj[-1][0], traj[-1][1], traj[-1][2]
    lo, hi = 0, len(ms) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ms[mid] >= m_target:
            lo = mid
        else:
            hi = mid
    a, b = 0.0, traj[hi][0] - traj[lo][0]
    for _ in range(80):
        mid = 0.5 * (a + b)
        if state_at(drop, field, traj, lo, mid)[5] >= m_target:
            a = mid
        else:
            b = mid
    return traj[lo][0] + 0.5 * (a + b), state_at(drop, field, traj, lo, 0.5 * (a + b)), traj[lo][2]


def node_windows(drop, field, traj, m_b):
    """Plateau crossings of the nodes of the Psat and L_v columns by the reference: (m_hi, m_lo, m[], dT[]) each.

    The slopes of the tabulated p_sat and L_v change at every node, and the solver's implicit step can go across one
    whole: at the end of that step the drop sits between the reference and the reference with the line of the
    interval it left continued, and from there it relaxes onto the former at the thermal rate. dT[] bounds that along
    the mass m[]: the distance between the two references at equal mass while the crossing step can last (NODE_STEP
    relaxation times past the node), and the largest such distance, decayed at the relaxation rate, afterwards;
    the window ends NODE_TIMES relaxation times after the step that crosses. m_hi and m_lo are its ends."""
    wins = []
    if drop.psat_mode != "table" and drop.lv_tab is None:
        return wins
    lo, hi = drop.tmin, drop.tmin + len(drop.lv_tab if drop.lv_tab is not None else drop.psat_tab) - 1
    for i in range(len(traj) - 1):
        s, cell = traj[i][1], traj[i][2]
        n0, n1 = int(math.floor(s[4])), int(math.floor(traj[i + 1][1][4]))
        node = max(n0, n1)
        if abs(n1 - n0) != 1 or not lo < node < hi:
            continue
        gas = field.gas(cell)
        if abs(s[4] - drop.t_quasi_steady(gas, drop.diameter(s[5]), s[2], s[3])) >= 0.02:
            continue                                           # not on the plateau: the steps are short there
        rising = n1 > n0
        dT = -1.0e-4 if rising else 1.0e-4                     # the relaxation rate, inside the interval it leaves
        f = lambda T: drop.rhs(np.array([s[0], s[1], s[2], s[3], T, s[5]]), gas)[4]
        lam = abs((f(s[4] + dT) - f(s[4])) / dT)
        alt = drop.continued(node, rising)
        h = NODE_TIMES / lam / 200.0
        t_node = (node - s[4]) / (traj[i + 1][1][4] - s[4]) * (traj[i + 1][0] - traj[i][0])   # from this point on
        ya, yr = s.copy(), s.copy()
        ma, ta, mr, tr = [s[5]], [s[4]], [s[5]], [s[4]]
        for _ in range(200):
            ya, _ = dp_step(lambda z: alt.rhs(z, gas), ya, h)
            yr, _ = dp_step(lambda z: drop.rhs(z, gas), yr, h)
            if not (np.all(np.isfinite(ya)) and np.all(np.isfinite(yr))) or min(ya[5], yr[5]) <= m_b:
                break
            ma.append(ya[5])
            ta.append(ya[4])
            mr.append(yr[5])
            tr.append(yr[4])
        ma, ta, mr, tr = (np.array(v[::-1]) for v in (ma, ta, mr, tr))          # by rising mass
        apart = np.abs(np.interp(mr, ma, ta) - tr)[::-1]                         # at equal mass, in time
        env, e, decay = [], 0.0, math.exp(-lam * h)
        for k, a in enumerate(apart):
            e = max(e * decay, a) if k * h - t_node <= NODE_STEP / lam else e * decay
            env.append(e)
        wins.append((mr[-1], mr[0], mr, np.array(env[::-1])))
    return wins


class Report:
    def __init__(self):
        self.nfail = 0

    def check(self, name, ok, text, gated=True):
        if gated:
            self.nfail += not ok
            tag = "PASS" if ok else "FAIL"
        else:
            tag = "info"
        print("  [%s] %-3s %s" % (tag, name, text))


def box_length(run):
    _, f, _ = read_tec_block(os.path.join(run, "OUTPUT", "gas-field.tec"))
    names = list(f.keys())[:3]
    lo = [f[n].min() for n in names]
    hi = [f[n].max() for n in names]
    return math.sqrt(sum((b - a) ** 2 for a, b in zip(lo, hi)))


def check_sweep(rep, run, drop, field, sweep_tag, streams, dt_q, gated, label):
    """Per-stream comparison of one sweep's IGLOO output with the reference integration."""
    out = os.path.join(run, "OUTPUT")
    scat, _ = read_particles(os.path.join(out, "part-scatter-A%s.dat" % sweep_tag))
    exits, ecol = read_particles(os.path.join(out, "part-outloc-A%s.dat" % sweep_tag))
    trajrows, _ = read_particles(os.path.join(out, "part-trajectories-A%s.dat" % sweep_tag))
    print("\n== %s (%s)%s" % ("sweep 0" if not sweep_tag else sweep_tag.strip("-"), label,
                              "" if gated else " -- per-stream checks informational"))
    results = {}
    w = {"ratio": 0.0, "arg": "", "qs": 0.0, "lead": -1e300, "lag": 0.0, "n": 0, "node": 0}
    for sid, s0, mdot in streams:
        m_b = burnout_mass(s0[5])
        traj = integrate(drop, field, s0, M_TAIL * m_b)    # on past the burnout mass: an exit step may end below it
        wins = node_windows(drop, field, traj, m_b)
        t_life, s_end, _ = pair_by_mass(drop, field, traj, m_b)
        results[sid] = (t_life, traj)
        frac = math.modf(sid * PHI)[0]
        rows = scat.get(sid, [])
        tplat, dmin_plat, dmax_plat = [], 1e300, 0.0
        for k, r in enumerate(rows, start=1):
            meas = {"x": r[0], "y": r[1], "u": r[3], "v": r[4], "T": r[6]}
            mm = r[8]
            tt, st, cell = pair_by_mass(drop, field, traj, mm)
            gas = field.gas(cell)
            der = drop.rhs(st, gas)
            d2 = second_derivative(drop, gas, st)
            dm = e_half(mm)
            ref = {"x": st[0], "y": st[1], "u": st[2], "v": st[3], "T": st[4]}
            idx = {"x": 0, "y": 1, "u": 2, "v": 3, "T": 4}
            # just after a node of the Psat column: what a solver step across it can leave in T at this mass
            node = max([float(np.interp(mm, wm, wd)) for m_hi, m_lo, wm, wd in wins if m_lo <= mm <= m_hi], default=0.0)
            w["node"] += node > 1e-7
            tk = dt_q * (k - frac)
            off = min(abs(tt - tk), max(TOL["t_lead"], TOL["t_lag"]))    # the sample's offset from its quantum
            for key in ref:
                floor = 1e-7 if key == "T" else 1e-9
                i = idx[key]
                # the print rounding, plus what the offset carries on the step: the position is cubic on the step
                # (off |y'|), the other columns linear (off |m' y''/m'' - y'| at most)
                if key in ("x", "y"):
                    place = off * abs(der[i])
                else:
                    place = off * (abs(der[5] * d2[i] / d2[5]) + abs(der[i])) if d2[5] != 0.0 else 0.0
                budget = abs(der[i] / der[5]) * dm + F_HALF + floor + place + (node if key == "T" else 0.0)
                rr = abs(meas[key] - ref[key]) / budget
                if rr > w["ratio"]:
                    w["ratio"], w["arg"] = rr, "stream %d sample %d %s: %.6f vs %.6f (budget %.1e)" % (
                        sid, k, key, meas[key], ref[key], budget)
            w["lead"] = max(w["lead"], tk - tt)
            w["lag"] = max(w["lag"], tt - tk)
            d_meas = drop.diameter(mm)
            tqs_ref = drop.t_quasi_steady(gas, drop.diameter(st[5]), st[2], st[3])
            if abs(st[4] - tqs_ref) < 0.02:                     # the plateau: the reference sits on its QS balance
                tqs = drop.t_quasi_steady(gas, d_meas, meas["u"], meas["v"])
                w["qs"] = max(w["qs"], abs(meas["T"] - tqs) - abs(st[4] - tqs_ref))
                tplat.append(meas["T"])
                dmin_plat, dmax_plat = min(dmin_plat, d_meas), max(dmax_plat, d_meas)
            w["n"] += 1
        # exit record and burnout
        ex = exits.get(sid)
        if ex is None:
            rep.check("E", False, "stream %d: no exit record" % sid, gated)
            continue
        e = ex[-1]
        dx = max(abs(e[0] - s_end[0]), abs(e[1] - s_end[1]))
        # the exit mass: a final trajectory row's (print-dtime), else the exit row's diameter, else the burnout bracket
        last = [r for r in trajrows.get(sid, []) if r[8] <= m_b]
        dcol = ecol.get("d<sub>p</sub>")
        m_exit, dm_exit = None, 0.0
        if last and last[-1][8] > 0.0:                       # IGLOO clamps the mass slot at 0: no exit mass then
            m_exit, dm_exit = last[-1][8], e_half(last[-1][8])
        elif dcol is not None and e[dcol] > 0.0:
            m_exit = drop.mass(e[dcol])
            dm_exit = 3.0 * m_exit * e_half(e[dcol]) / e[dcol]
        if m_exit is not None:
            te, se, ce = pair_by_mass_beyond(drop, field, traj, m_exit)
            der = drop.rhs(se, field.gas(ce))
            bT = abs(der[4] / der[5]) * dm_exit + F_HALF + 1e-6
            okT = abs(e[3] - se[4]) <= bT
            tdesc = "T %.6f vs %.6f at the exit mass %.6e (budget %.1e)" % (e[3], se[4], m_exit, bT)
            t_lo = t_hi = te
        else:
            # no exit mass printed: the exit step ends somewhere in [M_TAIL, 1] x m_b, where the reference is carried
            t9 = traj[-1][0]
            tb = [st[4] for t, st, _ in traj if t >= t_life] + [s_end[4]]
            okT = min(tb) - 1e-6 <= e[3] <= max(tb) + 1e-6
            tdesc = "T %.6f in [%.6f, %.6f] (the reference over m in [%.0e, 1] x %.3e)" % (
                e[3], min(tb), max(tb), M_TAIL, m_b)
            t_lo, t_hi = t_life, t9
        # a sample at every quantum up to the exit: (k - frac) dt_q <= t_exit, on IGLOO's clock when it prints it
        et = e[ecol["t"]] if "t" in ecol else None
        n_lo, n_hi = (int(math.floor(t / dt_q + frac)) for t in ((et, et) if et is not None else (t_lo, t_hi)))
        okn = n_lo <= len(rows) <= n_hi
        text = "stream %d: exit (%.6f, %.6f) vs (%.6f, %.6f); %s; %d samples (expect %s)" % (
            sid, e[0], e[1], s_end[0], s_end[1], tdesc, len(rows),
            "%d" % n_lo if n_lo == n_hi else "%d to %d" % (n_lo, n_hi))
        okt = True
        if et is not None:                                   # out-time: the exit time is printed
            okt = t_lo * (1 - TOL["t_exit"]) <= et <= t_hi * (1 + TOL["t_exit"])
            text += "; exit time %.8e s vs %s" % (et, "%.8e" % t_lo if t_lo == t_hi else "[%.8e, %.8e]" % (t_lo, t_hi))
        rep.check("E", dx <= F_HALF + 1e-6 and okT and okn and okt, text, gated)
        if tplat:
            print("        plateau: T %.4f..%.4f K over d %.2f..%.2f um; reference lifetime %.8e s, x-stop %.6f m"
                  % (min(tplat), max(tplat), dmin_plat * 1e6, dmax_plat * 1e6, t_life, s_end[0]))
    rtol_s = TOL["ratio"] if isinstance(field, UniformGas) else TOL["ratio_cells"]
    rep.check("S", w["ratio"] <= rtol_s,
              "%d samples paired by mass: worst residual/budget %.3f (tol %.1f) -- %s; %d just after a node of the "
              "Psat column" % (w["n"], w["ratio"], rtol_s, w["arg"], w["node"]), gated)
    rep.check("Q", w["qs"] <= TOL["T_qs"],
              "plateau samples vs the closed-form wet-bulb balance at their diameter: excess over the reference's "
              "own lag %.2e K (tol %.0e)" % (w["qs"], TOL["T_qs"]), gated)
    rep.check("C", w["lead"] <= TOL["t_lead"] and w["lag"] <= TOL["t_lag"],
              "sample clock: at most %.2e s before its quantum (tol %.0e), %.2e s after (tol %.0e)"
              % (w["lead"], TOL["t_lead"], w["lag"], TOL["t_lag"]), gated)
    return results


def pair_by_mass_beyond(drop, field, traj, m_target):
    """As pair_by_mass, continuing the last step past the burnout mass when m_target is below it."""
    if m_target >= traj[-1][1][5]:
        return pair_by_mass(drop, field, traj, m_target)
    t0, s0, cell = traj[-1]
    gas = field.gas(cell)
    f = lambda z: drop.rhs(z, gas)
    a, b = 0.0, 1e-6
    while dp_step(f, s0, b)[0][5] > m_target:
        b *= 2.0
    for _ in range(80):
        mid = 0.5 * (a + b)
        if dp_step(f, s0, mid)[0][5] > m_target:
            a = mid
        else:
            b = mid
    return t0 + 0.5 * (a + b), dp_step(f, s0, 0.5 * (a + b))[0], cell


def check_global(rep, run, sweep_tag, streams, results, drop, gas_exact):
    out = os.path.join(run, "OUTPUT")
    src = os.path.join(out, "part-source%s.tec" % sweep_tag)
    eul = os.path.join(out, "part-euler1%s.tec" % sweep_tag)
    if not os.path.exists(src):
        print("  [note] %s absent: no global check for this sweep" % os.path.basename(src))
        return
    names, f, _ = read_tec_block(src)
    vol = cell_volumes(f, names)
    mdot_in = sum(s[2] for s in streams)
    p_in = sum(s[2] * s[1][2] for s in streams)
    e_in = sum(s[2] * (drop.h(s[1][4]) + 0.5 * (s[1][2] ** 2 + s[1][3] ** 2)) for s in streams)
    sw, sfx, sfy, se = (float(f[n].sum()) for n in names[3:5] + names[5:6] + names[7:8])
    ok = (abs(sw - mdot_in) <= TOL["src_rel"] * mdot_in and abs(sfx - p_in) <= TOL["src_rel"] * abs(p_in)
          and abs(sfy) <= TOL["src_rel"] * abs(p_in) and abs(se - e_in) <= TOL["src_rel"] * abs(e_in))
    rep.check("G", ok, "source field vs the injected fluxes of the consumed parcels: wdot %.9e / %.9e kg/s, "
              "Fx %.9e / %.9e N, Fy %.2e, E %.9e / %.9e W"
              % (sw, mdot_in, sfx, p_in, sfy, se, e_in))
    if os.path.exists(eul):
        en, ef, _ = read_tec_block(eul)
        npv = float((ef[en[8]] * vol).sum())
        rhov = float((ef[en[3]] * vol).sum())
        npdot = [s[2] / s[1][5] for s in streams]
        ref_np = sum(nd * results[s[0]][0] for nd, s in zip(npdot, streams))
        ref_rho = sum(nd * mass_integral(results[s[0]][1]) for nd, s in zip(npdot, streams))
        ok = abs(npv - ref_np) <= TOL["eul_rel"] * ref_np and abs(rhov - ref_rho) <= TOL["eul_rel"] * ref_rho
        rep.check("G", ok, "Eulerian field: sum n_p V %.8e vs sum npdot t_life %.8e (%.1e); sum rho_p V %.8e vs "
                  "sum npdot int m dt %.8e (%.1e) (tol %.0e)" % (npv, ref_np, abs(npv / ref_np - 1), rhov, ref_rho,
                                                              abs(rhov / ref_rho - 1), TOL["eul_rel"]), gas_exact)


def check_coupler(rep, run, inp, log, gated):
    """The coupler's closed-domain lines against the final gas field: gas gained per species."""
    if not os.path.exists(log):
        print("  [note] %s absent: the coupler's closed-domain lines are not checked" % log)
        return
    lines = open(log, errors="replace").read().splitlines()
    rows = []
    for i, l in enumerate(lines):
        if "closed-domain mass balance" in l:
            for r in lines[i + 2:]:
                t = r.split()
                if len(t) != 4 or not t[0].isdigit():
                    break
                rows.append((int(t[0]), float(t[1]), float(t[2])))
    if not rows:
        print("  [note] %s carries no closed-domain lines" % log)
        return
    names, f, _ = read_tec_block(os.path.join(run, "OUTPUT", "gas-field.tec"))
    vol = cell_volumes(f, names)
    _, g, _ = read_tec_block(os.path.join(inp, "gas-ic.tec"))
    ric = [float(g[n].flat[0]) for n in g if re.match(r"rho\d+$", n)]   # the partial densities (not rhoNuTilde)
    ok, parts = True, []
    for isp, inj, gained in rows:
        field = float(((f["rho(%d)" % isp] - ric[isp - 1]) * vol).sum())
        ok = ok and abs(gained - field) <= e_half(gained) + 1e-18
        parts.append("species %d %.6e vs %.9e" % (isp, gained, field))
    rep.check("G", ok, "gas gained as the coupler reports it (%s) vs the final gas field minus the initial: %s"
              % (os.path.basename(log), "; ".join(parts)), gated)


def check_log(rep, log):
    """IGLOO's LK-interface: warning in the run log, whatever the sweep."""
    if not os.path.exists(log):
        print("  [note] %s absent: its LK-interface: lines are not checked" % log)
        return
    hits = [l.strip() for l in open(log, errors="replace") if "LK-interface:" in l]
    rep.check("L", not hits, "%s: %d LK-interface: line(s)%s" % (os.path.basename(log), len(hits),
                                                                 " -- " + hits[0] if hits else ""))


def cell_volumes(f, names):
    x, y, z = (f[n] for n in names[:3])
    dx = np.diff(x[0, 0, :])
    dy = np.diff(y[0, :, 0])
    dz = np.diff(z[:, 0, 0]) if z.shape[0] > 1 else np.array([1.0])
    return dz[:, None, None] * dy[None, :, None] * dx[None, None, :]


def mass_integral(traj):
    ts = np.array([p[0] for p in traj])
    ms = np.array([p[1][5] for p in traj])
    return float(np.trapz(ms, ts))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", nargs="?", default=".")
    ap.add_argument("--ini", default=None, help="input.ini to model (default RUN/input.ini)")
    ap.add_argument("--sweeps", default="0,last")
    ap.add_argument("--psat", default="auto", choices=("auto", "cc", "table"))
    ap.add_argument("--interface", default="auto", choices=("auto", "LK", "VLE"))
    ap.add_argument("--props", default=None, help="properties table to model (default RUN/INPUT/<phase>-properties.dat)")
    ap.add_argument("--log", default=None, help="run log: MI2's closed-domain lines, IGLOO's warnings (default RUN/run.log)")
    a = ap.parse_args()
    run = a.run
    ini = read_ini(a.ini or os.path.join(run, "input.ini"))
    inp = os.path.join(run, "INPUT")
    phase = ini.get("IGLOO-General", {}).get("phase", "part").split()[0]
    tok = {}
    for l in open(os.path.join(inp, "%s-phase.txt" % phase)).read().splitlines()[1:]:
        for t in l.split()[2:]:
            if "=" in t:
                k, v = t.split("=", 1)
                tok[k] = v
    props = read_columns_table(a.props or os.path.join(inp, "%s-properties.dat" % phase))
    drop = Droplet(ini, tok, props, a.psat, a.interface)
    tables = GasTables(inp)
    bc = ini["IGLOO-BC"]
    xs, ys = floats(bc["x"]), floats(bc["y"])
    n = len(xs)

    def rep_list(key, default):
        v = floats(bc.get(key, default))
        return v * n if len(v) == 1 else v

    zs, ups, ds, mds, t0s = (rep_list(k, dflt) for k, dflt in
                             (("z", "0"), ("up", "0"), ("diam", "0"), ("mdot", "0"), ("temp0", "0")))
    vps, wps = rep_list("vp", "0"), rep_list("wp", "0")
    streams = []
    for i in range(n):
        m0 = math.pi / 6.0 * drop.rho_l * ds[i] ** 3
        streams.append((i + 1, np.array([xs[i], ys[i], ups[i], vps[i], t0s[i], m0]), mds[i]))
    ncoup = int(ini.get("HYDRA-Parameters", {}).get("N-coupling", "1"))
    vref = sum(math.sqrt(s[1][2] ** 2 + s[1][3] ** 2) for s in streams) / n
    # the scatter quantum: fsample-traj of the INI the run read, IGLOO's default when it is absent
    run_ini = os.path.join(run, "input.ini")
    nsample = int((read_ini(run_ini) if os.path.exists(run_ini) else ini).get("IGLOO-General", {}).get(
        "fsample-traj", TRAJ_SAMPLE))
    dt_q = box_length(run) / vref / nsample
    rep = Report()
    print("evap-box reference: %s   p_sat %s   interface %s" % (os.path.abspath(run), drop.psat_mode,
                                                                   "LK" if drop.lk else "VLE"))
    for sw in a.sweeps.split(","):
        k = ncoup - 1 if sw == "last" else int(sw)
        tag = "" if k == 0 else "-sweep%d" % k
        if k == 0:
            _, f, _ = read_tec_block(os.path.join(inp, "gas-ic.tec"))
            names = list(f.keys())
            rn = [nm for nm in names if re.match(r"rho\d+$", nm)]   # the partial densities (not rhoNuTilde)
            st = tables.state([float(f[nm].flat[0]) for nm in rn], float(f["u"].flat[0]), float(f["v"].flat[0]),
                              float(f["w"].flat[0]), float(f["p"].flat[0]))
            for nm in names[3:]:
                if np.ptp(f[nm]) != 0.0:
                    raise SystemExit("[FAIL] the sweep-0 reference assumes a uniform initial gas (%s varies)" % nm)
            field, label = UniformGas(st), "initial gas: T %.6f K, rho %.6f kg/m3, mu %.6e Pa s, k %.6e W/m/K, gamma %.6f" % (
                st[4], st[0], st[5], st[6], st[7])
            gated = True
        elif k == ncoup - 1:
            gf = os.path.join(run, "OUTPUT", "gas-field.tec")
            field = CellGas(tables, gf)
            # a steady field carries its iteration as the solution time, negative when ORION prints the sign
            it_file = abs(float(re.search(r"SOLUTIONTIME=\s*([-+0-9.EeDd]+)", open(gf).read(2000)).group(1).replace("D", "E")))
            hist = [l.split() for l in open(os.path.join(run, "OUTPUT", "gas-residual-history.dat")) if l.strip()]
            it_last = int(hist[-1][0])
            gated = int(round(it_file)) == it_last
            label = "OUTPUT/gas-field.tec at MOSE iteration %d; the sweep sampled iteration %d" % (it_file, it_last)
        else:
            print("\n  [note] sweep %d: its gas field is not in OUTPUT; skipped" % k)
            continue
        results = check_sweep(rep, run, drop, field, tag, streams, dt_q, gated, label)
        check_global(rep, run, tag, streams, results, drop, gated)
        if k == ncoup - 1 and k > 0:
            check_coupler(rep, run, inp, a.log or os.path.join(run, "run.log"), gated)
    print("\n== run log")
    check_log(rep, a.log or os.path.join(run, "run.log"))
    print("\n%s: %d check(s) failed" % ("PASS" if rep.nfail == 0 else "FAIL", rep.nfail))
    return 0 if rep.nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
