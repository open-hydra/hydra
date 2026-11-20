#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set a CONSTANT viscosity and conductivity in INPUT/fluid-transport.dat (FLINT
real-fluid table), the constant properties of the analytical reference:
mu = 1.79352e-4 Pa s, k = mu cp / Pr with cp = 1004.5 J/(kg K), Pr = 0.70.
set_constant_cp.py sets the same cp in fluid-thermo.dat.

fluid-transport.dat is a Tecplot BLOCK file with 4 variables (pressure,
enthalpy, viscosity, conductivity): 3 header lines then 4 blocks of I*J*K
values, one per line. The pressure/enthalpy grid is kept verbatim.

Needed only after ATLAS has rebuilt a case's INPUT/: the versioned table, in
common/, is already flattened, and the cases reach it through links, which
this script refuses to write through.

Original backed up to INPUT/fluid-transport.dat.orig (once).

Usage, from a case folder:  python3 ../set_constant_transport.py
       or from anywhere:    python3 set_constant_transport.py <case folder>
"""

import os
import re
import shutil
import sys

CASE = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.getcwd()
F    = os.path.join(CASE, "INPUT", "fluid-transport.dat")

# --------------------------------------------------------------------------- #
# constants to impose  (edit here)
# --------------------------------------------------------------------------- #
VISCOSITY = 1.79352e-4        # [Pa s] molecular viscosity (MOSE/NASA value)
PRANDTL   = 0.70                   # [-]    laminar Prandtl number
CP        = 1004.5                 # [J/kg/K] cp = gamma*R/(gamma-1) = 1.4*287/0.4
CONDUCTIVITY = VISCOSITY * CP / PRANDTL   # [W/m/K]  k = mu*cp/Pr  (or set directly)


def main():
    if not os.path.isfile(F):
        sys.exit(f"no {F}: run from a case folder, or give the case folder as argument")
    if os.path.islink(F):
        sys.exit(f"{os.path.relpath(F, CASE)} is a link (to the shared, already flattened table "
                 "in common/): nothing to do. Rebuild INPUT/ with ATLAS in a copy first.")

    with open(F) as fh:
        lines = fh.read().splitlines()

    header, data = lines[:3], lines[3:]
    zone = header[2]

    def dim(tag):
        m = re.search(tag + r"\s*=\s*(\d+)", zone)
        if not m:
            raise ValueError(f"cannot find {tag} in ZONE line: {zone}")
        return int(m.group(1))

    I, J, K = dim("I"), dim("J"), dim("K")
    nval = I * J * K
    if len(data) != 4 * nval:
        raise ValueError(f"expected {4*nval} data lines (4 x {nval}), got {len(data)}")

    P = data[0:nval]              # Pressure  (kept verbatim)
    H = data[nval:2 * nval]       # Enthalpy  (kept verbatim)
    visc = [repr(VISCOSITY)]    * nval
    cond = [repr(CONDUCTIVITY)] * nval

    bak = F + ".orig"
    if not os.path.exists(bak):
        shutil.copy2(F, bak)
        print(f"  backup written: {os.path.relpath(bak, CASE)}")
    else:
        print(f"  backup already exists: {os.path.relpath(bak, CASE)} (kept)")

    with open(F, "w") as fh:
        fh.write("\n".join(header + P + H + visc + cond) + "\n")

    print(f"  fluid-transport.dat updated  (grid {I}x{J}x{K}, {nval} pts):")
    print(f"    Viscosity    = {VISCOSITY:.6e} Pa s   (constant)")
    print(f"    Conductivity = {CONDUCTIVITY:.6e} W/m/K (constant, Pr={PRANDTL}, cp={CP})")
    print(f"    Pressure / Enthalpy grid preserved")


if __name__ == "__main__":
    main()
