#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regenerate the figures embedded in the Hydra Verification & Validation page
(docs/vv/conjugate-flat-plate.md) from the frozen wall solutions shipped with
the test cases.

Sources
-------
  Tw-Flat-Plate/reference/ks_kf_*_wall.tec        hydra-MF  (MOSE + FUSS)
  AF-flatplate-Tw/reference/ks_kf_*_wall.tec      hydra-AF  (ARES + FUSS)

Both are `OUTPUT/hydra-wall.tec` files produced with `coupling-wall = true`.
Each holds two zones: the solid-side interface (`plate-B1F4`, written first)
and the fluid-side interface (`gas-B1F3` / `fluid-B1F3`). The solid zone
carries only Tw and qw; its y+/tau/pw columns are the placeholder 1.78625e-7.

The reference is the Lindstedt & Karvinen conjugate flat-plate correlation,
identical to the one in Tw-Flat-Plate/analytical.py.

Theming
-------
The SVGs are written to be *inlined* into the page (`{% include %}`), not
linked with <img>, so the Material theme's CSS can reach into them. They use
transparent backgrounds and `currentColor` for every structural element —
text, spines, ticks, grid — so the figures follow the light/dark theme
toggle. Data series use the theme's indigo/cyan/orange accents, which read
on both backgrounds. This matches docs/user/images/hmf-*.svg.

Usage
-----
    python3 test/make-vv-figures.py          # writes docs/vv/images/*.svg
    python3 test/make-vv-figures.py --show   # also display interactively

No solver run is required: the figures are built from the versioned solutions.
"""

import os
import re
import sys

import numpy as np
import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, os.pardir, "docs", "vv", "images"))

RATIOS = (40, 400, 4000)

# --------------------------------------------------------------------------- #
# Palette — Material indigo / cyan / deep-orange, legible on light and dark.
# THEME is a sentinel swapped for `currentColor` when the SVG is post-processed.
# --------------------------------------------------------------------------- #
THEME = "#010101"
C_MF = "#5c6bc0"        # hydra-MF  (indigo 400, theme primary family)
C_AF = "#ff7043"        # hydra-AF  (deep orange 400)
C_RATIO = {40: "#5c6bc0", 400: "#26c6da", 4000: "#ff7043"}

FONT_STACK = ("Roboto,-apple-system,BlinkMacSystemFont,Helvetica,Arial,sans-serif")

matplotlib.rcParams.update({
    "svg.fonttype": "none",        # keep real <text> nodes so CSS can recolor them
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "savefig.transparent": True,
    "text.color": THEME,
    "axes.labelcolor": THEME,
    "axes.edgecolor": THEME,
    "xtick.color": THEME,
    "ytick.color": THEME,
    "axes.titlecolor": THEME,
    "font.size": 10,
})


# --------------------------------------------------------------------------- #
# Tecplot reader
# --------------------------------------------------------------------------- #
def read_tec_zones(path):
    """Read every ZONE of a Tecplot ASCII BLOCK wall file.

    Layout per zone: 3 NODAL variables (x, y, z) of I*J*K values each,
    followed by the remaining variables CELLCENTERED. Returns a list of dicts,
    one per zone, each with the variable arrays plus 'xc' (cell-centre
    streamwise coordinate) and 'name'.
    """
    with open(path) as fh:
        lines = fh.readlines()

    var_names = re.findall(r'"([^"]+)"', lines[0])
    if not var_names:   # names without quotes (VARIABLES = x y z y+ ...), as hydra-AF writes hydra-wall.tec now
        var_names = lines[0].split("=", 1)[1].split()
    starts = [i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE")]

    zones = []
    for zi, s in enumerate(starts):
        e = starts[zi + 1] if zi + 1 < len(starts) else len(lines)
        hdr = lines[s]
        name = re.search(r"T\s*=\s*([^,]+),", hdr).group(1).strip()
        I, J, K = (int(re.search(r"\b%s\s*=\s*(\d+)" % t, hdr).group(1)) for t in "IJK")

        tokens = []
        for line in lines[s + 1:e]:
            for tok in line.split():
                try:
                    tokens.append(float(tok))
                except ValueError:
                    pass

        n_nodal, n_nodal_vars = I * J * K, 3
        n_cc_vars = len(var_names) - n_nodal_vars
        n_cc = (len(tokens) - n_nodal_vars * n_nodal) // n_cc_vars

        off, z = 0, {"name": name}
        for nm in var_names[:n_nodal_vars]:
            z[nm] = np.array(tokens[off:off + n_nodal]); off += n_nodal
        for nm in var_names[n_nodal_vars:]:
            z[nm] = np.array(tokens[off:off + n_cc]); off += n_cc

        xn = z["x"][:I]
        z["xc"] = 0.5 * (xn[:-1] + xn[1:])
        zones.append(z)
    return zones


# --------------------------------------------------------------------------- #
# Lindstedt & Karvinen conjugate flat-plate correlation
# --------------------------------------------------------------------------- #
def tw_laminar_forced(x, T_B, T_inf, U_inf, rho_f, mu_f, k_f, cp_f, k_s, t):
    """Interface temperature Tw(x) for a conjugate laminar flat plate."""
    x = np.asarray(x, dtype=float)
    Pr = mu_f * cp_f / k_f
    Re_x = U_inf * x / (mu_f / rho_f)
    Bi_f = 0.332 / 0.5 * np.sqrt(Re_x) * Pr ** (1.0 / 3.0) * (k_f * t) / (k_s * x)
    theta_w = Bi_f * (0.5 + Bi_f) / (1.0 + Bi_f) ** 2
    return {"Re_x": Re_x, "Bi_f": Bi_f, "Tw": T_B - (T_B - T_inf) * theta_w}


# --------------------------------------------------------------------------- #
# Shared case state (identical for both variants by construction)
# --------------------------------------------------------------------------- #
R, GAM, PR = 287.0, 1.4, 0.7
T_INF, P_INF, T_B = 293.0, 1.0e5, 393.0
L, DELTA, RE_L, MA = 1.0, 0.05, 4.55e5, 0.2

RHO = P_INF / (R * T_INF)
U = MA * np.sqrt(GAM * R * T_INF)
MU = RHO * U * L / RE_L
CP = GAM / (GAM - 1.0) * R
KF = MU * CP / PR

# Solid conductivity as set in each case's atlas.ini
KS = {40: 10.2746, 400: 102.746, 4000: 1027.46}

VARIANTS = {
    "MF": dict(
        path=os.path.join(HERE, "Tw-Flat-Plate", "reference", "ks_kf_%d_wall.tec"),
        label="hydra-MF (MOSE + FUSS)", color=C_MF, marker="o",
    ),
    "AF": dict(
        path=os.path.join(HERE, "AF-flatplate-Tw", "reference", "ks_kf_%d_wall.tec"),
        label="hydra-AF (ARES + FUSS)", color=C_AF, marker="x",
    ),
}

X_AN = np.linspace(1e-3 * L, L, 4000)   # analytic sampling; skips x = 0 singularity
X_TRUST = 0.02 * L                      # leading-edge cells excluded from metrics


def load(variant, ratio):
    """Return (Re_x_num, Tw_num, Re_x_an, Tw_an) for one case."""
    solid = read_tec_zones(VARIANTS[variant]["path"] % ratio)[0]
    an = tw_laminar_forced(X_AN, T_B, T_INF, U, RHO, MU, KF, CP, KS[ratio], DELTA)
    return (solid["xc"] * RE_L / L / RE_SCALE, solid["Tw"],
            an["Re_x"] / RE_SCALE, an["Tw"])


SUP = {"-": "\u207b", "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3",
       "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078",
       "9": "\u2079"}


def sup(n):
    """Unicode superscript for an integer exponent (no mathtext)."""
    return "".join(SUP[c] for c in str(int(n)))


def log_ticks(ax):
    """Label a log axis as 10<sup>n</sup> without mathtext."""
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(
            lambda v, _: "10%s" % sup(round(np.log10(v))) if v > 0 else ""))


# Figures plot Re_x in units of 1e5 so no scientific-notation offset (which
# matplotlib renders as mathtext) is needed anywhere.
RE_SCALE = 1e5


def style_axes(ax):
    ax.set_xlim(0, RE_L / RE_SCALE)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.grid(True, color=THEME, alpha=0.15, linewidth=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.8)


def legend(ax, **kw):
    lg = ax.legend(frameon=False, **kw)
    return lg


# --------------------------------------------------------------------------- #
# SVG post-processing: make the figure follow the Material light/dark theme
# --------------------------------------------------------------------------- #
def theme_svg(path):
    """Rewrite a matplotlib SVG so it inherits the page's theme colours.

    - structural strokes/fills (the THEME sentinel) become `currentColor`
    - every <text> gets fill="currentColor" unless it carries its own colour
    - figure/axes background patches become transparent
    - the DejaVu metrics font is swapped for the theme's font stack
    """
    with open(path) as fh:
        svg = fh.read()

    svg = svg.replace("fill:%s" % THEME, "fill:currentColor")
    svg = svg.replace("stroke:%s" % THEME, "stroke:currentColor")
    svg = svg.replace("fill:#ffffff", "fill:none")
    svg = svg.replace("font-family:DejaVu Sans", "font-family:%s" % FONT_STACK)

    # <text ...> nodes: colour them with currentColor unless a fill is set inline
    def tag_text(m):
        attrs = m.group(1)
        if "fill=" in attrs:
            return m.group(0)
        return '<text fill="currentColor"%s>' % attrs

    svg = re.sub(r"<text([^>]*)>", tag_text, svg)

    # activate the .matplotlib-axis rule already present in stylesheets/extra.css
    svg = svg.replace('<g id="matplotlib.axis_', '<g class="matplotlib-axis" id="matplotlib.axis_')

    with open(path, "w") as fh:
        fh.write(svg)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    fig.savefig(p, bbox_inches="tight", transparent=True)
    plt.close(fig)
    theme_svg(p)
    print("  wrote %s" % os.path.relpath(p, os.path.join(HERE, os.pardir)))


# --------------------------------------------------------------------------- #
# Figure 1 — interface temperature, both variants against the correlation
# --------------------------------------------------------------------------- #
XLABEL = "Re_x   (\u00d710\u2075)"


def fig_profiles():
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.9))
    ylim = {40: (310, 396), 400: (368, 396), 4000: (390, 393.6)}

    for ax, ratio in zip(axes, RATIOS):
        _, _, rex_a, tw_a = load("MF", ratio)
        ax.plot(rex_a, tw_a, color=THEME, lw=2.0, alpha=0.75,
                label="Lindstedt correlation", zorder=2)
        for key in ("MF", "AF"):
            v = VARIANTS[key]
            rex_n, tw_n, _, _ = load(key, ratio)
            ax.plot(rex_n, tw_n, v["marker"], color=v["color"], ms=5,
                    mew=1.4, markevery=8, label=v["label"], zorder=3)
        ax.axhline(T_B, color=THEME, ls=":", lw=1.0, alpha=0.5, zorder=1)
        ax.annotate("TB", xy=(0.06, T_B), xytext=(0, 3), textcoords="offset points",
                    fontsize=9, alpha=0.7, color=THEME)
        ax.set_title("ks / kf = %d" % ratio, fontsize=11, pad=8)
        ax.set_xlabel(XLABEL, fontsize=10)
        ax.set_ylim(*ylim[ratio])
        style_axes(ax)

    axes[0].set_ylabel("Tw   [K]", fontsize=10)
    legend(axes[0], fontsize=8.5, loc="lower right")
    fig.tight_layout()
    save(fig, "conjugate-flat-plate-profiles.svg")


# --------------------------------------------------------------------------- #
# Figure 2 — deviation from the correlation
# --------------------------------------------------------------------------- #
def fig_error():
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.6), sharey=True)

    for ax, ratio in zip(axes, RATIOS):
        for key in ("MF", "AF"):
            v = VARIANTS[key]
            rex_n, tw_n, rex_a, tw_a = load(key, ratio)
            err = np.abs(tw_n - np.interp(rex_n, rex_a, tw_a))
            ax.semilogy(rex_n, err, color=v["color"], lw=1.8, label=v["label"])
        ax.axvspan(0, X_TRUST * RE_L / L / RE_SCALE, color=THEME, alpha=0.08,
                   lw=0, zorder=0)
        ax.set_title("ks / kf = %d" % ratio, fontsize=11, pad=8)
        ax.set_xlabel(XLABEL, fontsize=10)
        ax.set_ylim(1e-3, 3e1)
        style_axes(ax)
        log_ticks(ax)

    axes[0].set_ylabel("| Tw \u2212 Tw,LK |   [K]", fontsize=10)
    legend(axes[0], fontsize=8.5, loc="lower left")
    axes[2].annotate("shaded: leading-edge singularity\nof the correlation",
                     xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                     fontsize=8, alpha=0.65, color=THEME)
    fig.tight_layout()
    save(fig, "conjugate-flat-plate-error.svg")


# --------------------------------------------------------------------------- #
# Figure 3 — cross-solver consistency (hydra-AF vs hydra-MF)
# --------------------------------------------------------------------------- #
def fig_cross_solver():
    fig, ax = plt.subplots(figsize=(8.4, 3.9))
    for ratio in RATIOS:
        rex, tw_mf, _, _ = load("MF", ratio)
        _, tw_af, _, _ = load("AF", ratio)
        ax.plot(rex, tw_af - tw_mf, color=C_RATIO[ratio], lw=1.9,
                label="ks / kf = %d" % ratio)
    ax.axhline(0.0, color=THEME, lw=0.9, ls=":", alpha=0.6)
    ax.axvspan(0, X_TRUST * RE_L / L / RE_SCALE, color=THEME, alpha=0.08,
               lw=0, zorder=0)
    ax.set_xlabel(XLABEL, fontsize=10)
    ax.set_ylabel("Tw (hydra-AF) \u2212 Tw (hydra-MF)   [K]", fontsize=10)
    ax.set_ylim(-0.85, 0.18)   # leading-edge cells run off-scale (shaded band)
    legend(ax, fontsize=9, loc="lower right")
    style_axes(ax)
    fig.tight_layout()
    save(fig, "conjugate-flat-plate-cross-solver.svg")


# --------------------------------------------------------------------------- #
# Metrics printed to stdout — the numbers quoted on the V&V page
# --------------------------------------------------------------------------- #
def report():
    print("\nFree stream: rho=%.4f kg/m3  U=%.3f m/s  mu=%.5e Pa s  k_f=%.5f W/m/K"
          % (RHO, U, MU, KF))
    print("             Re_L=%.3e  Ma=%.2f  Pr=%.2f  T_inf=%.0f K  T_B=%.0f K\n"
          % (RE_L, MA, PR, T_INF, T_B))

    hdr = ("%-9s %6s %9s %11s %11s %9s %9s"
           % ("variant", "ks/kf", "k_s", "Tw(L/2) an", "Tw(L/2) num", "max|e|*", "RMSE*"))
    print(hdr); print("-" * len(hdr))
    for key in ("MF", "AF"):
        for ratio in RATIOS:
            rex_n, tw_n, rex_a, tw_a = load(key, ratio)
            e = tw_n - np.interp(rex_n, rex_a, tw_a)
            m = rex_n > X_TRUST * RE_L / L / RE_SCALE
            print("%-9s %6d %9.4f %11.3f %11.3f %9.3f %9.3f"
                  % ("hydra-" + key, ratio, KS[ratio],
                     np.interp(RE_L / 2 / RE_SCALE, rex_a, tw_a), np.interp(RE_L / 2 / RE_SCALE, rex_n, tw_n),
                     np.max(np.abs(e[m])), np.sqrt(np.mean(e[m] ** 2))))
    print("* excluding Re_x < %.0e (leading-edge singularity of the correlation)"
          % (X_TRUST * RE_L / L))

    print("\nCross-solver (AF - MF) and interface continuity:")
    hdr = "%6s %13s %12s %14s %14s" % ("ks/kf", "max|dTw| [K]", "RMS dTw [K]",
                                       "cont. MF [K]", "cont. AF [K]")
    print(hdr); print("-" * len(hdr))
    for ratio in RATIOS:
        rex, tw_mf, _, _ = load("MF", ratio)
        _, tw_af, _, _ = load("AF", ratio)
        m = rex > X_TRUST * RE_L / L / RE_SCALE
        d = tw_af - tw_mf
        cont = []
        for key in ("MF", "AF"):
            z = read_tec_zones(VARIANTS[key]["path"] % ratio)
            cont.append(np.max(np.abs(z[0]["Tw"] - z[1]["Tw"])))
        print("%6d %13.3f %12.3f %14.2e %14.2e"
              % (ratio, np.max(np.abs(d[m])), np.sqrt(np.mean(d[m] ** 2)), cont[0], cont[1]))

    print("\nWall temperature range (whole plate):")
    for key in ("MF", "AF"):
        for ratio in RATIOS:
            _, tw, _, _ = load(key, ratio)
            print("  hydra-%s  ks/kf=%-5d  Tw = %.2f .. %.2f K" % (key, ratio, tw.min(), tw.max()))


if __name__ == "__main__":
    print("Regenerating Hydra V&V figures")
    fig_profiles()
    fig_error()
    fig_cross_solver()
    report()
    if "--show" in sys.argv:
        plt.show()
