# evap-box-mom — momentum and energy bookkeeping of IGLOO's source coupling

## The case

A variant of `test/evap-box-ta`: the same 20x20x1 box of quiescent N2 (1 atm, 350 K, `MESH/mesh.tec`, 5 mm cells)
and the same four streams of 20 µm water drops injected at x = 0.02 m, y = 0.03/0.045/0.055/0.07 m with 0.5 m/s in
+x and 1e-6 kg/s each (`[IGLOO-BC]`), but

- the drops do **not** evaporate: `input.ini` sets no `evaporation` key, so IGLOO's default `NoEvaporation`
  applies (GPB has no "off" value for the material token) and the drops exchange only drag (`Stokes`) and heat
  (`Ranz-Marshall`) with the gas — no mass, no latent heat, no enthalpy datum;
- the box is periodic in x **and y** (faces 1–4, code 201 for the gas and for the dispersed phase, written by ATLAS
  from `[perx]`/`[pery]`), symmetry in z: no boundary force acts on the gas in x or y. With y symmetry walls the
  wall pressure force moves the gas y-momentum by 12–23 times |Fx (t1 − t0)| per block, which is not a coupler
  quantity; the z-faces exert no net force because the flow is uniform in z;
- a drop that has stopped (0.55 mm after injection, τ_p = 1.09 ms) never leaves the box, so `time-end = 0.05`
  ends every parcel at 50 ms of parcel time (46 τ_p, 10.6 thermal times): the final state goes to
  `OUTPUT/part-snapshot-A*.dat` (velocity 1e-5 m/s, T within 1e-3 K of the gas), and every sweep ends;
- time accurate, six coupling steps of 50 MOSE iterations, `IGLOO-iter = 1` (source mode), `relaxation = false`,
  `[MI2-Coupling] report = true`.

In source mode `IGLOOxMOSE` (`src/app/MI2/MI2_Coupling.f90`) adds IGLOO's absolute rates `sourceMom` [N] and
`sourceEn` [W] to MOSE's residual, and holds them over a MOSE block (IGLOO solves between blocks). The gas must
therefore gain exactly Fx (t1 − t0) of x-momentum and Q (t1 − t0) of energy over each block, Fx and Q being the
totals the report prints for that block. Not a physics validation: it validates the coupler's momentum and energy
bookkeeping (the mass channel of the same branch is gated by `evap-box-ta`).

## The checks (`check.py`)

From the `[MI2-report]` lines of every block that had IGLOO sources in force (steps 2–6; step 1 has none):

| id | check | tolerance | measured on G (2026-09-30) |
|---|---|---|---|
| M1x | px1 − px0 = Fx (t1 − t0) | 1e-8 relative | ≤ 1.2e-14 |
| M1y, M1z | py, pz likewise, relative to \|Fx (t1 − t0)\| (Fy is 1e-4 of Fx, Fz = 0) | 1e-8 | ≤ 3.2e-10, 0 |
| M1s, M1o | the sum over the gated blocks; no change over the sourceless block | 1e-8 | 4e-15, 0 |
| M3a | E1 − E0 = Q (t1 − t0) | 1e-5 relative | 3.1–3.8e-6 (see below) |
| M2a | Fx = Σ ṁ u_inj = 2e-6 N | 1e-2 (u_end is the gas velocity, < 1e-3 m/s) | ≤ 1.7e-4 |
| M2b, M2c | Fx = Σ ṁ (u_inj − u_end), Fy = −Σ ṁ v_end, with the end velocities of the snapshot of the producing sweep | 1e-5 of Σ ṁ u_inj (the snapshot prints 1e-6 m/s) | ≤ 3e-7, ≤ 2e-6 |
| M3b | Q = Σ ṁ [cp (T0 − T_end) + (u_inj² − \|v_end\|²)/2], T_end from the same snapshot | 1e-6 relative | ≤ 3e-8 |
| NV | \|px1 − px0\| ≥ ½ Σ ṁ u_inj (t1 − t0); Fx ≥ ½ Σ ṁ u_inj; \|Q\| ≥ ½ Σ ṁ cp (T_gas − T0) = 0.25 W; one snapshot row per stream at time-end; blocks contiguous in time | — | ratios 1.0, 1.0, 2.0 |

A log without report lines, or without a block with sources, is a FAIL (no vacuous pass).

**The energy floor.** E1 − E0 misses Q (t1 − t0) by +3.1 to +3.8e-10 J per block (3.1–3.8e-6 relative), always
with the same sign. Measured: 2× that with 100 MOSE iterations per block, the same amount per iteration at
CFL 0.1 (blocks 3× shorter), 1700× less with `heat = NoHeat` (same momentum source, no thermal change), and
nothing of the kind in the momentum channel (1e-14). That is a per-iteration loss of about 7e-12 J inside MOSE's
energy update, consistent with its round trip: `integration-variables = cons` stores primitives, and at every RK
stage `prim2cons` rebuilds ρE from the pressure that FLINT's `cons2prim` obtained by a Newton iteration on T
which exits at |Δe/e| < 1e-6 (`lib/FLINT/src/lib/Lib_ThermoTransport.f90`); ρu ↔ u has no such step. M3a's
tolerance is set above that floor; the defects it must catch (a volume factor, a sign, a swapped channel, a
missing stage) are 1e7, 2, 1 and 0.33.

## Falsification

`test/mi2-checks/bin/hydra_MI2.Gred-iglvol-2026-09-30` (the source branch of `IGLOOxMOSE` multiplying
`sourceMom` by the cell volume): M1x fails with px1 − px0 = 1.2500000e-7 × Fx (t1 − t0) on every block
(V_cell = 1.25e-7 m³), M1y (Fy is scaled the same way: 1.5e-4 of |Fx (t1 − t0)| lost), M1s and NVp fail, while
Fx, M2 and the energy checks are unchanged. The F binary (`hydra_MI2.F-2026-09-29`, no coupling report) fails at
the first check of check.py (no report lines).

## Run

    test/mi2-checks/p12_igloo_momentum_check.sh [scratch-root] [exe]      # the gate: PASS/FAIL lines, exit 1 on a FAIL
    test/mi2-checks/run_case.sh test/evap-box-mom <scratch> 1 <log> <exe>  # a run by hand (1.3 s at 1 thread)
    python3 check.py <scratch> --log <log>                                 # the checks on that run

The gate also checks the case files (faces 1–4 code 201, 5–6 code 300, in both bc files), the log (GAME OVER,
exit 0, no error or parcel-abort line, every sweep ended at time-end) and prints the run time.

## Regenerate

    export ATLASDIR=<ATLAS tree at ebeaf4f> ATLAS_PYTHON=<python with cantera>
    test/mi2-checks/regen_case.sh test/evap-box-mom <scratch> GPB BCB ICB [--install]

`INPUT/` is reproducible byte for byte from `input.ini` (a second regeneration reports every file identical). The
gas tables are those of `evap-box-ta`; the dispersed phase has a relative enthalpy table (no `h0`) and no `Psat`
column, and its material line carries no model token (`A 1`).
