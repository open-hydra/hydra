# ice-igloo-burn-box — one run, two burning phases, two models

One run of `hydra-MI2` with two condensed phases that burn on different models, each solved by its own solver: the
Lagrangian phase `partL` (IGLOO, `combustion = direct`, the burning section of `test/igloo-burn-box-direct`) and the
Eulerian phase `partE` (ICE, `combustion = Beckstead`, the burning section of `test/ice-burn-box-beckstead`). The
model is a word on each material's line of its phase file, so nothing in the run is common to the two but the gas.

## The case

A box of quiescent N2 at 400 K, periodic in x and y (`MESH/` is `test/evap-box-mom`'s). The gas has three species,
N2, Ar and He, the last two absent at the start: each phase burns into its own, `partL` into Ar
(`particle-mechanism = a-ar`, `a-ar.yaml`: `A => Ar`) and `partE` into He (`a-he.yaml`: `A => He`), each material
with `h0 = 1.0e6`. ICE's cloud is uniform and at rest (loading 0.3 of the gas density, 20 µm, 300 K); IGLOO injects
the four streams of `test/igloo-burn-box-direct`. Cycle 1 runs one forward-Euler ICE step, so the source ICE stores
is that of the initial state; cycle 2 keeps it frozen (`ICE-iter = 0`) while IGLOO is solved before each block of 50
time-accurate MOSE iterations. `INPUT/` is what ATLAS writes for `input.ini` and the two reaction files
(`GPB BCB ICB`): two phase files and two products files, `partL-*` and `partE-*`.

## The oracle

Each solver against its own closed form, with `V` the box volume, `T_ref` = 298.15 K and `dh = h298 - h_c(T_ref)`
of each material (`h298` from its record):

```
ICE     mdot = rho pi K-burn X-eff/(2 n) d^(3-n)      M_E = n_p mdot V
        Qc   = n_p Nu pi d k_g (T_g - T_p) V          Q_E = M_E (c T_p + dh - beta q) - Qc
IGLOO   M_L  = SUM mdot (1 - r)                       Q_L = SUM mdot [h_c(T0) + u^2/2 - r h_c(T_end) + (1 - r) dh]
        r = (d_end/d0)^3, the end records of the sweep whose sources the report carries
```

Over a block `[t0, t1]` the gas energy rises by `(Q_E + Q_L)(t1 - t0)`. Over the run Ar gains `M_L t`, He gains
`M_E t` and N2 nothing: the species balance closes per phase. The gas warms as the first law of the rigid box says,
with the internal energies of the three species.

## Running it

```
test/mi2-checks/p24_two_solver_burn_check.sh [scratch-root] [exe]
```

runs the case, then `check.py <run-dir> <log>`; it also checks that ICE's model given by `[ICE-Physics] combustion`
instead of its material line writes the same bytes and makes the coupler print one warning: the material `A` is in
both solvers, burns on two models, and one phase file does not say so.
