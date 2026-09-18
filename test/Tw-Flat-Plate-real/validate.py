#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation of the three real-fluid conjugate flat-plate cases:

    ks_kf_40_real/   ks_kf_400_real/   ks_kf_4000_real/

Each case is a laminar flat plate with conjugate heat transfer, real-fluid
(CoolProp) Air with Cp/mu/k forced constant so that Pr = mu*cp/k = 0.70
(see each case's set_constant_cp.py / set_constant_transport.py). The cases
differ only in the solid conductivity k_s, i.e. the ratio k_s/k_f.

The interface temperature Tw(x) extracted from each simulation
(<case>/OUTPUT/hydra-wall.tec, fluid-side wall) is compared against the
Lindstedt & Karvinen analytical correlation, in the same 3-panel style as
Tw-Flat-Plate/analytical.py (the ideal-gas reference).

Run from anywhere:  python3 validate_ks_kf_real.py
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def read_tec_wall(path):
    """
    Read the first ZONE from a Tecplot ASCII BLOCK-format wall file.

    Header format:
      VARIABLES = "x" "y" "z" "y+" "tauX" "tauY" "tauZ" "pw" "Tw" "qw"
      ZONE T = ..., I=N, J=1, K=2, DATAPACKING=BLOCK,
           VARLOCATION=([1-3]=NODAL,[4-10]=CELLCENTERED), ...

    NODAL       : x, y, z              -> I*J*K values each
    CELLCENTERED: remaining variables  -> (I-1) values each

    Returns dict with 'x' (cell-centre streamwise coord) and 'Tw'.
    """
    with open(path) as f:
        lines = f.readlines()

    var_names = re.findall(r'"([^"]+)"', lines[0])

    zone_hdr = lines[1]
    I = int(re.search(r'\bI\s*=\s*(\d+)', zone_hdr).group(1))
    J = int(re.search(r'\bJ\s*=\s*(\d+)', zone_hdr).group(1))
    K = int(re.search(r'\bK\s*=\s*(\d+)', zone_hdr).group(1))

    tokens = []
    for line in lines[2:]:
        if line.strip().upper().startswith('ZONE'):
            break
        for tok in line.split():
            try:
                tokens.append(float(tok))
            except ValueError:
                pass

    n_nodal = I * J * K
    n_nodal_vars = 3
    n_cc_vars = len(var_names) - n_nodal_vars
    n_cc = (len(tokens) - n_nodal_vars * n_nodal) // n_cc_vars

    offset = 0
    parsed = {}
    for name in var_names[:n_nodal_vars]:
        parsed[name] = np.array(tokens[offset: offset + n_nodal])
        offset += n_nodal
    for name in var_names[n_nodal_vars:]:
        parsed[name] = np.array(tokens[offset: offset + n_cc])
        offset += n_cc

    x_nodes = parsed['x'][:I]
    x_cc = 0.5 * (x_nodes[:-1] + x_nodes[1:])
    return {'x': x_cc, 'Tw': parsed['Tw']}


def tw_laminar_forced(x, T_B, T_inf, U_inf, rho_f, mu_f, k_f, cp_f, k_s, t):
    """
    Local interface temperature Tw(x) for a conjugate flat plate
    (Lindstedt & Karvinen correlation). All in SI, temperatures in K.
    """
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0.0):
        raise ValueError("All x values must be > 0.")

    nu = mu_f / rho_f
    Pr = mu_f * cp_f / k_f
    Re_x = U_inf * x / nu

    Bi_f = 0.332 / 0.5 * (Re_x ** 0.5) * (Pr ** (1.0 / 3.0)) * (k_f * t) / (k_s * x)
    theta_w = Bi_f * (0.5 + Bi_f) / (1.0 + Bi_f) ** 2
    Tw = T_B - (T_B - T_inf) * theta_w
    return {"x": x, "Re_x": Re_x, "Pr": Pr, "Bi_f": Bi_f, "theta_w": theta_w, "Tw": Tw}


if __name__ == "__main__":

    # ================================================================
    # Common parameters (shared by the three cases)
    # from atlas.ini + set_constant_*.py + mesh.ini
    # ================================================================
    T_inf = 293.0                  # free-stream static temperature   [K]
    p_inf = 1.0e5                  # ambient static pressure          [Pa]
    g_inf = 81.6055                # inflow mass flux  rho*U          [kg/(m2 s)]
    R     = 287.0                  # specific gas constant            [J/(kg K)]

    mu_inf = 1.79352e-4            # imposed constant viscosity       [Pa s]
    cp_inf = 1004.5                # imposed constant cp              [J/(kg K)]
    Pr     = 0.70                  # imposed laminar Prandtl          [-]
    k_inf  = mu_inf * cp_inf / Pr  # fluid conductivity  k = mu cp/Pr [W/(m K)]

    rho_inf = p_inf / (R * T_inf)  # free-stream density              [kg/m3]
    U_inf   = g_inf / rho_inf      # free-stream velocity             [m/s]

    delta = 0.05                   # plate thickness (Block2 height)  [m]
    L     = 1.0                    # plate length                     [m]
    T_B   = 393.0                  # bottom-wall temperature ([T_plate] BC) [K]
    ReL   = g_inf * L / mu_inf     # plate-length Reynolds number

    # -- per-case solid conductivity k_s (from each case's atlas.ini) and
    #    frozen reference wall file in reference/ (solid-side converged wall,
    #    copied from <case>/OUTPUT/hydra-wall.tec) --
    cases = {
        40:   {"sol": "reference/ks_kf_40_wall.tec",   "k_s": 10.2746,  "ylim": (310, 395)},
        400:  {"sol": "reference/ks_kf_400_wall.tec",  "k_s": 102.746,  "ylim": (365, 395)},
        4000: {"sol": "reference/ks_kf_4000_wall.tec", "k_s": 1027.46,  "ylim": (390, 395)},
    }
    ks_kf_ratios = [40, 400, 4000]

    print("=" * 56)
    print("Real-fluid conjugate flat plate  --  shared state")
    print("=" * 56)
    print(f"  T_inf = {T_inf:.1f} K   T_B = {T_B:.1f} K")
    print(f"  rho   = {rho_inf:.4f} kg/m3   U = {U_inf:.3f} m/s")
    print(f"  mu    = {mu_inf:.5e} Pa s   cp = {cp_inf:.1f} J/kg/K")
    print(f"  k_f   = {k_inf:.5f} W/m/K   Pr = {Pr}")
    print(f"  Re_L  = {ReL:.3e}   Ma = {U_inf/np.sqrt(1.4*R*T_inf):.3f}")
    print("=" * 56)

    # ================================================================
    # Analytical + numerical solution, one panel per ratio
    # ================================================================
    x = np.linspace(1e-3 * L, L, 600)          # skip leading-edge singularity

    fig, axes = plt.subplots(3, 1, figsize=(8, 13), sharex=True)
    fig.suptitle(
        "Real-fluid conjugate flat-plate interface temperature\n"
        rf"($T_\infty={T_inf:.0f}\,$K, $T_B={T_B:.0f}\,$K, "
        rf"$Re_L={ReL:.2e}$, $Pr={Pr}$)",
        fontsize=16, fontweight="bold",
    )
    fig.subplots_adjust(top=0.88)

    for ax, ratio in zip(axes, ks_kf_ratios):
        c   = cases[ratio]
        k_s = c["k_s"]

        # -- analytical (Lindstedt & Karvinen) --
        res   = tw_laminar_forced(x, T_B, T_inf, U_inf, rho_inf,
                                  mu_inf, k_inf, cp_inf, k_s, delta)
        Re_x  = res["Re_x"]
        Tw_LK = res["Tw"]

        # -- numerical (fluid-side wall, frozen reference) --
        tec_path = os.path.join(HERE, c["sol"])
        tec = read_tec_wall(tec_path)
        Re_x_num = tec['x'] * ReL / L

        # -- error metrics --
        err  = tec['Tw'] - np.interp(Re_x_num, Re_x, Tw_LK)
        rmse = np.sqrt(np.mean(err ** 2))
        i_mid = np.argmin(np.abs(Re_x - ReL / 2.0))
        print(f"\n  ks/kf = {ratio:>4d}   k_s = {k_s:.4f} W/m/K  (k_s/k_f = {k_s/k_inf:.1f})")
        print(f"    Tw analytic / numeric (x=L/2) = "
              f"{Tw_LK[i_mid]:.3f} / {np.interp(ReL/2.0, Re_x_num, tec['Tw']):.3f} K")
        print(f"    max|err| = {np.max(np.abs(err)):.3f} K   RMSE = {rmse:.3f} K")

        # -- plot --
        ax.plot(Re_x, Tw_LK, color="#2471A3", linewidth=2.5,
                label="Lindstedt correlation")
        ax.plot(Re_x_num, tec['Tw'], "o", color="#FF1A1A", markersize=5,
                markevery=5, label="Numerical (HYDRA)")
        ax.axhline(T_B, color="k", linestyle="--", linewidth=1.0,
                   label=f"$T_B = {T_B:.0f}$ K")
        ax.set_ylabel("$T_w$ [K]", fontsize=14)
        ax.set_title(f"$k_s / k_f = {ratio}$", fontsize=14)
        ax.tick_params(labelsize=12)
        ax.legend(fontsize=10, loc="lower right")
        ax.set_ylim(*c["ylim"])
        ax.grid(True, alpha=0.4)

    print("=" * 56)

    axes[-1].set_xlabel("$Re_x$", fontsize=13)
    axes[-1].set_xlim(0, ReL)
    axes[-1].xaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(1e5))
    for ax in axes:
        ax.xaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    plt.tight_layout()
    out = os.path.join(HERE, "validate.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  figure written: {os.path.relpath(out, HERE)}")
    plt.show()
