#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
F    = os.path.join(HERE, "INPUT", "fluid-transport.dat")

# --------------------------------------------------------------------------- #
# constants to impose  (edit here)
# --------------------------------------------------------------------------- #
VISCOSITY = 1.79352e-4        # [Pa s] molecular viscosity (MOSE/NASA value)
PRANDTL   = 0.70                   # [-]    laminar Prandtl number
CP        = 1004.5                 # [J/kg/K] cp = gamma*R/(gamma-1) = 1.4*287/0.4
CONDUCTIVITY = VISCOSITY * CP / PRANDTL   # [W/m/K]  k = mu*cp/Pr  (or set directly)


def main():
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
        print(f"  backup written: {os.path.relpath(bak, HERE)}")
    else:
        print(f"  backup already exists: {os.path.relpath(bak, HERE)} (kept)")

    with open(F, "w") as fh:
        fh.write("\n".join(header + P + H + visc + cond) + "\n")

    print(f"  transport.dat updated  (grid {I}x{J}x{K}, {nval} pts):")
    print(f"    Viscosity    = {VISCOSITY:.6e} Pa s   (constant)")
    print(f"    Conductivity = {CONDUCTIVITY:.6e} W/m/K (constant, Pr={PRANDTL}, cp={CP})")
    print(f"    Pressure / Enthalpy grid preserved")
    print(f"  -> with g=59.29: Re/m = g/mu = {59.29/VISCOSITY:.3e}  (target 5.0e6, M=0.2)")


if __name__ == "__main__":
    main()
