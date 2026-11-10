# ice-box — the frozen ICE source on a uniform gas

A gate for the MOSE–ICE coupling of `hydra-MI2`: what `ICExMOSE` hands the gas (momentum
`-source(2:4)·vol`, energy `-source(nbase-1)·vol`, sign, factor and slot) and what `MOSExICE`
hands ICE (`rho, u, T, mu, k` in the layout its closures read) are checked against closed forms
that follow from `input.ini` alone.

## The case

`MESH/mesh.tec` is the 20×20×1 box of `test/evap-box-ta` (0.1 × 0.1 × 0.005 m). Faces 1 and 2
are periodic in x (code 201 for every phase), faces 3–6 symmetry planes. The gas is quiescent
N2 at 101325 Pa and 400 K. The ICE phase (closure MK, one family) is uniform: bulk density 0.3
of the gas density, particle diameter 20 µm, velocity 1 m/s along x, temperature 300 K,
material density 1000 kg/m³, cp 4184 J/(kg K), Stokes drag (`C_d = 24/Re`) and Stokes heat
(`Nu = 2`), no evaporation, no radiation.

Schedule (`[HYDRA-Parameters]`): cycle 1 runs no MOSE and one forward-Euler ICE step, so the
source ICE stores is the one evaluated at the initial state; cycle 2 runs 200 coupling steps of
one time-accurate MOSE iteration each with that source frozen (`ICE-iter = 0`). One iteration
per step, not one block of 200, so that the coupling report covers every time step.

## The oracle

With `V` the box volume, `u_g`, `T_g` the gas state ICE saw in cycle 1 and `rho_mat` the
material density,

```
R_p    = (3 rho_p / (4 pi n_p rho_mat))^(1/3)          ICE's radius from its own state (= d/2)
tau    = rho_mat (2 R_p)^2 / (18 mu_g)                  Stokes relaxation time
F      = rho_p (u_p - u_g) / tau · V                    force on the gas            [N]
Q_heat = 2 Nu k_g pi R_p (T_g - T_p) n_p · V, Nu = 2    heat into the particles     [W]
W      = rho_p u_p (u_p - u_g) / tau · V                drag work on the particles  [W]
```

and over cycle 2, of duration `t`: `Δ(Σ rho u V) = F t`, `Δ(Σ rho E V) = (W - Q_heat) t`.
`mu_g` and `k_g` are interpolated from `INPUT/gas-transport.dat` as FLINT does (linear between
the integer-kelvin rows) at `T_g = p/(rho R)` of `INPUT/gas-ic.tec`. Everything is independent
of both integrators and of the time step. Numbers (2026-11-10, MOSE `b3e2ba8`): `tau = 9.98867e-4 s`,
`F = 1.28168e-2 N`, `Q_heat = 1234.96 W`, `t = 3.7159e-4 s`, gas velocity gain 0.1116 m/s,
gas energy loss 0.4589 J (the gas cools by 14.4 K; the source is frozen, so it does not relax).

`check.py` (run by the gate) compares, from the `[MI2-report]` lines of cycle 2 and the output
fields: A0 the reported `tau` (uniform) with the closed form to 1e-10; A1 the gas momentum
gain with `F t` of the report and with the closed form to 1e-8 (total and every step; y, z
components zero to 1e-10 of x); A2 the gas energy gain with `Q t` of the report and with
`(W - Q_heat) t` to 1e-8 (every step to 1e-7: the energy snapshots carry ±1e-11 J from MOSE's
e→T inversion, 2e-9 of a step's gain); A3 non-vacuity (velocity gain ≥ 0.05 u_p, heat
exchanged ≥ 1e-6 J); A4 the final gas uniform to 1e-10; S1 the table interpolation reproduces
the field's `mil`, `kl`; S2 the 200 `ice-to-gas` lines identical (frozen); P the IC radius is
`dp/2`. The gate script adds A5 log hygiene, A6 the absence of the ICE mass-channel warning, and a
refusal leg: the same case with one cycle that advances MOSE and ICE together (`MOSE-iter = 200`,
`ICE-iter = 1`, `N-coupling = 1`) while both are time-accurate must stop with `advances MOSE and
ICE together` and a non-zero exit, while the case's own schedule is accepted.

## Run the gate

```bash
test/mi2-checks/p10_ice_box_check.sh [scratch-root] [exe]      # 1 thread, about 1 s
python3 -B test/ice-box/check.py <run-dir> <log>              # the numeric checks alone
```

## Regenerate the inputs

```bash
export ATLASDIR=<ATLAS ebeaf4f or later>  ATLAS_PYTHON=<python with cantera>
test/mi2-checks/regen_case.sh test/ice-box <scratch> GPB BCB ICB --install
python3 -B test/ice-box/make_case.py
```

ATLAS ICB writes the condensed velocity equal to the gas velocity and the number density with
pi rounded to 3.14 (`src/ICB/builder_dp.f90`); `make_case.py` rewrites the `up_11` column with
`[ice-box] up` and the `np_11` column with `rho_p / (rho_mat 4/3 pi (dp/2)^3)` in ATLAS's own
number format, every other line kept byte for byte. A regeneration without `--install` then
reports every file identical except `part-ic.tec` (those two columns); `make_case.py` on the
regenerated file gives the committed one byte for byte.

## Files

- `input.ini` — the case (solver sections, then the ATLAS sections); comments say what each block does
- `MESH/mesh.tec`, `INPUT/` — the mesh and the ATLAS-generated inputs (`part-ic.tec` patched as above)
- `make_case.py` — the IC patch; `check.py` — the numeric checks
- `../mi2-checks/p10_ice_box_check.sh` — the gate
