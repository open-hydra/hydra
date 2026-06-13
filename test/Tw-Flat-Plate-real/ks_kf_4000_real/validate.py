#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation of the ks_kf_4000_real conjugate flat-plate simulation.

The case is a laminar flat plate with conjugate heat transfer:
  - real-fluid (CoolProp) Air, with Cp/mu/k forced constant so that
    Pr = mu*cp/k = 0.70  (see set_constant_cp.py / set_constant_transport.py)
  - solid plate of thickness t, bottom wall held at T_B (BC [T_plate])
  - fluid free stream at T_inf, mass flux g_inf -> Re_L = 4.55e5

The interface temperature Tw(x) extracted from the numerical solution
(OUTPUT/hydra-wall.tec, fluid-side wall) is compared against the
Lindstedt & Karvinen analytical correlation for the conjugate problem.

Run:  python3 validate.py
"""

import re
import numpy as np
import matplotlib.pyplot as plt


def read_tec_wall(path):
    """
    Read the first ZONE from a Tecplot ASCII BLOCK-format wall file.

    Expected header format:
      VARIABLES = "x" "y" "z" "y+" "tauX" "tauY" "tauZ" "pw" "Tw" "qw"
      ZONE T = ..., I=N, J=1, K=2, DATAPACKING=BLOCK,
           VARLOCATION=([1-3]=NODAL,[4-10]=CELLCENTERED), ...

    NODAL       : x, y, z              -> I*J*K values each
    CELLCENTERED: remaining variables  -> (I-1) values each

    Returns
    -------
    dict with keys:
      'x'  : cell-centre streamwise coordinates (length I-1)
      'Tw' : wall temperature at cell centres   (length I-1)
    """
    with open(path) as f:
        lines = f.readlines()

    # Parse variable names from line 1
    var_names = re.findall(r'"([^"]+)"', lines[0])

    # Parse zone dimensions from line 2
    zone_hdr = lines[1]
    I = int(re.search(r'\bI\s*=\s*(\d+)', zone_hdr).group(1))
    J = int(re.search(r'\bJ\s*=\s*(\d+)', zone_hdr).group(1))
    K = int(re.search(r'\bK\s*=\s*(\d+)', zone_hdr).group(1))

    # Collect numeric tokens from Zone 1 only (stop at second ZONE header)
    tokens = []
    for line in lines[2:]:
        if line.strip().upper().startswith('ZONE'):
            break
        for tok in line.split():
            try:
                tokens.append(float(tok))
            except ValueError:
                pass

    # Layout: 3 NODAL vars x (I*J*K) + (nvar-3) CELLCENTERED vars x n_cc
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

    # Cell-centre x: midpoints of the first K-layer node x-values
    x_nodes = parsed['x'][:I]
    x_cc = 0.5 * (x_nodes[:-1] + x_nodes[1:])

    return {'x': x_cc, 'Tw': parsed['Tw']}


def tw_laminar_forced(x, T_B, T_inf, U_inf, rho_f, mu_f, k_f, cp_f, k_s, t):
    """
    Local interface temperature Tw(x) for a conjugate flat plate with:
      - bottom wall at constant temperature T_B
      - top wall cooled by laminar forced convection
      - 1D conduction through the solid thickness
      - Lindstedt & Karvinen analytical correlation

    All temperatures in K, properties in SI.
    """
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0.0):
        raise ValueError("All x values must be > 0.")

    nu = mu_f / rho_f                 # kinematic viscosity [m^2/s]
    Pr = mu_f * cp_f / k_f            # Prandtl number [-]
    Re_x = U_inf * x / nu             # local Reynolds number

    if np.any(Re_x >= 5e5):
        print("Warning: some x have Re_x >= 5e5; laminar correlation may be invalid there.")

    # Laminar mean-Biot-type parameter
    Bi_f = 0.332 / 0.5 * (Re_x ** 0.5) * (Pr ** (1.0 / 3.0)) * (k_f * t) / (k_s * x)

    # Non-dimensional wall temperature theta_w = (T_B - Tw)/(T_B - T_inf)
    theta_w = Bi_f * (0.5 + Bi_f) / (1.0 + Bi_f) ** 2
    Tw = T_B - (T_B - T_inf) * theta_w

    return {"x": x, "Re_x": Re_x, "Pr": Pr, "Bi_f": Bi_f, "theta_w": theta_w, "Tw": Tw}


if __name__ == "__main__":

    # ================================================================
    # Case parameters (from atlas.ini + set_constant_*.py + mesh.ini)
    # ================================================================
    # -- Free-stream / fluid --
    T_inf = 293.0                 # free-stream static temperature   [K]
    p_inf = 1.0e5                 # ambient static pressure          [Pa]
    g_inf = 81.6055               # inflow mass flux  rho*U          [kg/(m2 s)]
    R     = 287.0                 # specific gas constant            [J/(kg K)]

    mu_inf = 1.79352e-4           # imposed constant viscosity       [Pa s]
    cp_inf = 1004.5               # imposed constant cp              [J/(kg K)]
    Pr     = 0.70                 # imposed laminar Prandtl          [-]
    k_inf  = mu_inf * cp_inf / Pr # fluid conductivity  k = mu cp/Pr [W/(m K)]

    rho_inf = p_inf / (R * T_inf)          # free-stream density     [kg/m3]
    U_inf   = g_inf / rho_inf              # free-stream velocity     [m/s]

    # -- Solid plate (phase 'plate' in atlas.ini) --
    k_s   = 1027.46               # solid thermal conductivity       [W/(m K)]
    delta = 0.05                  # plate thickness (Block2 height)  [m]
    L     = 1.0                   # plate length                     [m]
    T_B   = 393.0                 # bottom-wall temperature ([T_plate] BC) [K]

    ReL    = g_inf * L / mu_inf            # plate-length Reynolds number
    ks_kf  = k_s / k_inf

    print("=" * 56)
    print("ks_kf_4000_real  --  case summary")
    print("=" * 56)
    print(f"  T_inf  = {T_inf:8.2f} K     T_B = {T_B:8.2f} K")
    print(f"  rho    = {rho_inf:8.4f} kg/m3  U   = {U_inf:8.3f} m/s")
    print(f"  mu     = {mu_inf:.5e} Pa s   cp  = {cp_inf:.1f} J/kg/K")
    print(f"  k_f    = {k_inf:.5f} W/m/K   Pr  = {Pr}")
    print(f"  k_s    = {k_s:.4f} W/m/K   t   = {delta} m   L = {L} m")
    print(f"  Re_L   = {ReL:.3e}     Ma  = {U_inf/np.sqrt(1.4*R*T_inf):.3f}")
    print(f"  k_s/k_f = {ks_kf:.2f}")
    print("=" * 56)

    # ================================================================
    # Analytical solution
    # ================================================================
    x = np.linspace(1e-3 * L, L, 600)          # skip leading-edge singularity
    res   = tw_laminar_forced(x, T_B, T_inf, U_inf, rho_inf,
                              mu_inf, k_inf, cp_inf, k_s, delta)
    Re_x  = res["Re_x"]
    Tw_LK = res["Tw"]

    i_mid = np.argmin(np.abs(Re_x - ReL / 2.0))
    print(f"  Bi_mean               = {np.mean(res['Bi_f']):.4f}")
    print(f"  Tw (analytic, x=L/2)  = {Tw_LK[i_mid]:.3f} K")

    # ================================================================
    # Numerical solution from the simulation (fluid-side wall)
    # ================================================================
    tec = read_tec_wall("OUTPUT/hydra-wall.tec")
    Re_x_num = tec['x'] * ReL / L              # x_cc -> Re_x

    # Interpolate analytic onto numerical x to compute the error
    Tw_LK_at_num = np.interp(Re_x_num, Re_x, Tw_LK)
    err = tec['Tw'] - Tw_LK_at_num
    rmse = np.sqrt(np.mean(err ** 2))
    print(f"  Tw (numeric, x=L/2)   = {np.interp(ReL/2.0, Re_x_num, tec['Tw']):.3f} K")
    print(f"  max |Tw_num - Tw_an|  = {np.max(np.abs(err)):.3f} K")
    print(f"  RMSE                  = {rmse:.3f} K")
    print("=" * 56)

    # ================================================================
    # Plot
    # ================================================================
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(
        "Conjugate flat plate  -  ks_kf_4000_real\n"
        rf"($T_\infty={T_inf:.0f}\,$K, $T_B={T_B:.0f}\,$K, "
        rf"$Re_L={ReL:.2e}$, $Pr={Pr}$, $k_s/k_f={ks_kf:.0f}$)",
        fontsize=14, fontweight="bold",
    )
    fig.subplots_adjust(top=0.88)

    ax.plot(Re_x, Tw_LK, color="#2471A3", linewidth=2.5,
            label="Lindstedt correlation")
    ax.plot(Re_x_num, tec['Tw'], "o", color="#FF1A1A", markersize=5,
            markevery=5, label="Numerical (HYDRA)")
    ax.axhline(T_B, color="k", linestyle="--", linewidth=1.0,
               label=f"$T_B = {T_B:.0f}$ K")
    ax.set_ylabel("$T_w$ [K]", fontsize=13)
    ax.set_xlabel("$Re_x$", fontsize=13)
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(True, alpha=0.4)
    ax.set_xlim(0, ReL)
    ax.xaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(1e5))
    ax.xaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    plt.tight_layout()
    plt.savefig("validate.png", dpi=150, bbox_inches="tight")
    print("  figure written: validate.png")
    plt.show()
