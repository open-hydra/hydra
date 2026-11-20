"""
Conjugate flat plate (ideal gas, hydra-MF): the interface temperature Tw(x) of the
three cases ks_kf_40, ks_kf_400, ks_kf_4000 against the Lindstedt & Karvinen
correlation.

    python3 analytical.py           the frozen solutions, reference/ks_kf_*_wall.tec
                                    -> Tw_analytical.png
    python3 analytical.py --runs    the cases' own runs, <case>/OUTPUT/hydra-wall.tec
                                    -> Tw_analytical_runs.png (a case not run is skipped)

The errors are measured where the correlation holds, x > 2 cm, as on the V&V page;
the leading-edge cells, where the correlation is singular, are printed apart.
Run from anywhere.
"""

import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
X_TRUST = 0.02      # [m] leading-edge cells (x below) left out of the error metrics, as on the V&V page


def read_tec_wall(path):
    """
    Read the first ZONE from a Tecplot ASCII BLOCK-format wall file.

    Expected header format:
      VARIABLES = "x" "y" "z" "y+" "tauX" "tauY" "tauZ" "pw" "Tw" "qw"
      ZONE T = ..., I=N, J=1, K=2, DATAPACKING=BLOCK,
           VARLOCATION=([1-3]=NODAL,[4-10]=CELLCENTERED), ...

    NODAL   : x, y, z  → I*J*K values each
    CELLCENTERED: remaining 7 variables  → (I-1) values each (150 for I=151)

    Returns
    -------
    dict with keys:
      'x'   : cell-centre streamwise coordinates (length I-1)
      'Tw'  : wall temperature at cell centres    (length I-1)
    """
    with open(path) as f:
        lines = f.readlines()

    # Parse variable names from line 1
    var_names = re.findall(r'"([^"]+)"', lines[0])
    if not var_names:   # names without quotes (VARIABLES = x y z y+ ...), as the current hydra writes hydra-wall.tec
        var_names = lines[0].split("=", 1)[1].split()

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

    # Layout: 3 NODAL vars × (I*J*K) + 7 CELLCENTERED vars × n_cc
    n_nodal = I * J * K          # nodes per NODAL variable  (302)
    n_nodal_vars = 3             # x, y, z
    n_cc_vars = len(var_names) - n_nodal_vars          # 7
    n_cc = (len(tokens) - n_nodal_vars * n_nodal) // n_cc_vars  # 150

    # Unpack NODAL → CELLCENTERED data in order
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


def tw_laminar_forced(
    x,
    T_B,
    T_inf,
    U_inf,
    rho_f,
    mu_f,
    k_f,
    cp_f,
    k_s,
    t
):
    """
    Compute local interface temperature Tw(x) for a conjugate flat plate
    with:
      - bottom wall at constant temperature T_B
      - top wall cooled by laminar forced convection
      - 1D conduction through the solid thickness
      - Lindstedt & Karvinen analytical correlation

    Parameters
    ----------
    x : array_like
        Streamwise locations [m], must be > 0.
    T_B : float
        Temperature at the back side of the solid [K or °C]
    T_inf : float
        Free-stream fluid temperature [K or °C]
    U_inf : float
        Free-stream velocity [m/s]
    rho_f : float
        Fluid density [kg/m^3]
    mu_f : float
        Fluid dynamic viscosity [Pa·s]
    k_f : float
        Fluid thermal conductivity [W/m/K]
    cp_f : float
        Fluid specific heat [J/kg/K]
    k_s : float
        Solid thermal conductivity [W/m/K]
    t : float
        Solid thickness [m]

    Returns
    -------
    results : dict
        Dictionary containing x, Re_x, Pr, Bi_f, theta_w, Tw
    """
    x = np.asarray(x, dtype=float)

    if np.any(x <= 0.0):
        raise ValueError("All x values must be > 0.")
    if U_inf <= 0 or rho_f <= 0 or mu_f <= 0 or k_f <= 0 or cp_f <= 0 or k_s <= 0 or t <= 0:
        raise ValueError("All transport/material properties and thickness must be positive.")

    # Fluid properties
    nu = mu_f / rho_f                     # kinematic viscosity [m^2/s]
    Pr = mu_f * cp_f / k_f                # Prandtl number [-]

    # Local Reynolds number
    Re_x = U_inf * x / nu

    # Optional laminar validity check for a flat plate
    if np.any(Re_x >= 5e5):
        print("Warning: some x locations have Re_x >= 5e5; laminar flat-plate correlation may not be valid there.")

    # Laminar mean-Biot-type parameter from the paper
    # Bi_f = 0.332/0.5 * Re_x^(1/2) * Pr^(1/3) * (k_f * t) / (k_s * x)
    Bi_f = 0.332/0.5 * (Re_x ** (1.0 / 2.0)) * (Pr ** (1.0 / 3.0)) * (k_f * t) / (k_s * x)

    # Non-dimensional wall temperature:
    # theta_w = (T_B - T_w)/(T_B - T_inf)
    theta_w = Bi_f * (0.5 + Bi_f) / (1.0 + Bi_f) ** 2

    # Local interface temperature
    Tw = T_B - (T_B - T_inf) * theta_w

    return {
        "x": x,
        "Re_x": Re_x,
        "Pr": Pr,
        "Bi_f": Bi_f,
        "theta_w": theta_w,
        "Tw": Tw,
    }


if __name__ == "__main__":
    runs = "--runs" in sys.argv[1:]

    # ================================================================
    # SECTION 1 — Boundary Conditions Analysis  (from setup.f90)
    # ================================================================
    # Gas properties
    Pr    = 0.7            # Prandtl number
    Ma    = 0.2            # inlet Mach number
    ReL   = 4.55e5         # Reynolds number based on plate length
    
    gam   = 1.4
    R     = 287.0          # specific gas constant            [J/(kg·K)]
    
    T_inf = 293.0          # free-stream static temperature   [K]
    p_inf = 1.0e5          # free-stream static pressure                  [Pa]

    # Bottom-plate boundary condition
    T_B = T_inf + 100.0    # solid bottom-wall temperature    [K]

    # Geometry
    L     = 1.0            # plate length                     [m, simulation units]
    delta = 0.05           # plate thickness                  [m, simulation units]


    # --- Derived fluid thermodynamic state ---
    a_inf   = np.sqrt(gam * R * T_inf)     # free-stream speed of sound        [m/s]
    U_inf   = Ma * a_inf                   # free-stream velocity              [m/s]
    rho_inf = p_inf / (R * T_inf)          # free-stream density               [kg/m³]
    g_inf   = rho_inf * U_inf              # free-stream mass flux             [kg/(m²·s)]
    mu_inf  = g_inf * L / ReL              # free-stream dynamic viscosity     [Pa·s]
    nu_inf  = mu_inf / rho_inf             # free-stream kinematic viscosity   [m²/s]
    cp_inf  = gam / (gam - 1.0) * R        # free-stream specific heat (cp)    [J/(kg·K)]
    k_inf   = mu_inf * cp_inf / Pr         # free-stream conductivity        [W/(m·K)]
    pt_inf  = p_inf * (1.0 + (gam - 1.0) / 2.0 * Ma**2) ** (gam / (gam - 1.0))

    print("=" * 56)
    print("=" * 56)
    print(f"GAS PROPERTIES")
    print()
    print(f"  gam  = {gam}")
    print(f"  R    = {R}  J/(kg K)")
    print(f"  mu   = {mu_inf:>8.5e}  Pa s")
    print(f"  kf   = {k_inf:>7.5f}  W/(m K)")
    print(f"  Pr   = {Pr}")
    print()
    print("=" * 56)
    print("=" * 56)
    print()
    print("=" * 56)
    print("=" * 56)
    print("BOUNDARY CONDITIONS SUMMARY")
    print()
    print(f"  -- INFLOW --")
    print(f"  g  = {g_inf:>8.4f}  kg/(m2 s)")
    print(f"  pt = {pt_inf/1e5:>8.5f}  bar")
    print(f"  T  = {T_inf:>8.1f}  K")
    print()
    print(f"  -- AMBIENT --")
    print(f"  p  = {p_inf/1e5:>8.1f}  bar")
    print(f"  pt = {pt_inf/1e5:>8.5f}  bar")
    print(f"  T  = {T_inf:>8.1f}  K")
    print()
    print(f"  -- FLOW STATE --")
    print(f"  rho  = {rho_inf:>12.5f}  kg/m3")
    print(f"  U    = {U_inf:>12.4f}  m/s")
    print(f"  ReL  = {ReL:>12.2e}")
    print()
    print(f"  -- SOLID GEOMETRY & BC --")
    print(f"  Plate length  L  = {L}  m")
    print(f"  Plate thick   d  = {delta}  m")
    print(f"  Bottom temp  T_B = {T_B}  K")
    print("=" * 56)

    # ================================================================
    # SECTION 2 — Analytical solution for ks/kf = 40, 400, 4000
    # ================================================================
    # Streamwise coordinate (skip the leading edge singularity)
    x = np.linspace(1e-3 * L, L, 600)

    ks_kf_ratios = [40, 400, 4000]
    # Solid conductivity as set in each case's INPUT/plate-properties.dat (the nominal ratio times k_f would be
    # 10.2948, 0.2 % more); the same values as test/make-vv-figures.py, so the numbers below are the V&V page's
    k_solid = {40: 10.2746, 400: 102.746, 4000: 1027.46}

    # The solution file of each ratio: the frozen one, or the case's own run with --runs
    def wall_file(ratio):
        if runs:
            return os.path.join(HERE, f"ks_kf_{ratio}", "OUTPUT", "hydra-wall.tec")
        return os.path.join(HERE, "reference", f"ks_kf_{ratio}_wall.tec")

    print(f"  solutions: {'the cases own runs, <case>/OUTPUT/hydra-wall.tec' if runs else 'the frozen ones, reference/'}")
    print(f"  errors over x > {100 * X_TRUST:.0f} cm (Re_x > {X_TRUST * ReL / L:.2g}); the leading-edge cells apart")

    fig, axes = plt.subplots(3, 1, figsize=(8, 13), sharex=True)
    fig.suptitle(
        "Conjugate flat-plate interface temperature \n"
        r"($T_\infty = 293\,\mathrm{K}$, $T_B = 393\,\mathrm{K}$, "
        r"$Re_L = 4.55\times10^5$, $Pr = 0.7$)",
        fontsize=16,
        fontweight="bold",
    )
    # Add vertical space between the suptitle and the first subplot
    fig.subplots_adjust(top=0.88)

    ylims = {40: (310, 395), 400: (365, 395), 4000: (390, 395)}

    found = 0
    for ax, ratio in zip(axes, ks_kf_ratios):
        k_s = k_solid[ratio]   # solid conductivity for this case [W/(m·K)]

        # -- Lindstedt & Karvinen correlation (accounts for 2-D solid effects) --
        res = tw_laminar_forced(
            x=x,
            T_B=T_B,
            T_inf=T_inf,
            U_inf=U_inf,
            rho_f=rho_inf,
            mu_f=mu_inf,
            k_f=k_inf,
            cp_f=cp_inf,
            k_s=k_s,
            t=delta,
        )
        Re_x  = res["Re_x"]
        Bi_f  = res["Bi_f"]
        Tw_LK = res["Tw"]

        # -- Print mid-chord values --
        tw_mid = np.interp(L / 2, x, Tw_LK)
        print(f"\n  ks/kf = {ratio:>4d}   k_s = {k_s:.4f} W/(m K)  (k_s/k_f = {k_s / k_inf:.2f})")
        print(f"    Bi_mean = {np.mean(Bi_f):.4f}")
        print(f"    Tw (Lindstedt, x=L/2)  = {tw_mid:.3f} K")

        # -- Plot the correlation --
        ax.plot(Re_x, Tw_LK, color="#2471A3", linewidth=2.5,
                label="Lindstedt Correlation")
        ax.axvspan(0, X_TRUST * ReL / L, color="0.5", alpha=0.15,
                   label=f"$x < {100 * X_TRUST:.0f}$ cm: not in the errors")

        # -- Numerical solution from .tec file (Zone 1) --
        path = wall_file(ratio)
        if not os.path.isfile(path):
            print(f"    {os.path.relpath(path, HERE)} not found: case not run, skipped")
        else:
            found += 1
            tec = read_tec_wall(path)
            xc, tw = tec['x'], tec['Tw']
            Re_x_num = xc * ReL / L     # x_cc → Re_x  (= x_cc * U_inf / nu)
            err = tw - tw_laminar_forced(xc, T_B, T_inf, U_inf, rho_inf,
                                         mu_inf, k_inf, cp_inf, k_s, delta)["Tw"]
            m = xc > X_TRUST
            print(f"    Tw (numerical, x=L/2)  = {np.interp(L / 2, xc, tw):.3f} K")
            print(f"    x > {100 * X_TRUST:.0f} cm ({m.sum()} cells): max|err| = {np.max(np.abs(err[m])):.3f} K"
                  f"   RMSE = {np.sqrt(np.mean(err[m] ** 2)):.3f} K")
            print(f"    x < {100 * X_TRUST:.0f} cm ({(~m).sum()} cells, leading edge): "
                  f"largest |err| = {np.max(np.abs(err[~m])):.3f} K (not judged)")
            ax.plot(Re_x_num, tw, "o", color="#FF1A1A", markersize=5,
                    markevery=5, label="Numerical")
        ax.axhline(T_B, color="k", linestyle="--", linewidth=1.0,
                   label=f"$T_B = {T_B:.0f}$ K")
        ax.set_ylabel("$T_w$ [K]", fontsize=14)
        ax.set_title(f"$k_s / k_f = {ratio}$",
                     fontsize=14)
        ax.tick_params(labelsize=12)
        ax.legend(fontsize=10, loc="lower right")
        ax.set_ylim(*ylims[ratio])
        ax.grid(True, alpha=0.4)

    if runs and not found:
        sys.exit("no case has OUTPUT/hydra-wall.tec: run the cases first")

    axes[-1].set_xlabel("$Re_x$", fontsize=13)
    axes[-1].set_xlim(0, 5e5)
    axes[-1].xaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(1e5))
    for ax in axes:
        ax.xaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    plt.tight_layout()
    out = os.path.join(HERE, "Tw_analytical_runs.png" if runs else "Tw_analytical.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  figure written: {os.path.relpath(out, HERE)}")
    plt.show()
