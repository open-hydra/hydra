# evap-box — reference solution

## Literature

None. The case is hydra's own construction: the closed-domain mass and energy gate of the MI2 coupler
(`plan-bucket/MI2-coupler-plan.md`, Phase 6). Its `input.ini` header says so — "deliberately NOT a physics
validation case — it validates the COUPLER's bookkeeping". Nothing in its parameters comes from a paper: 20 µm
water drops at 320 K and 0.5 m/s in quiescent N2 at 350 K and 101325 Pa, with the standard IGLOO closures. The
nearest literature case, Miller–Harstad–Bellan 1998 Fig. 2 (IGLOO's `tests/evaporation/mhb98-water`), is a
1.05 mm water drop in 298 K air, so it is not this case.

The reference is therefore the one the physics admits: **a single-droplet integration** of the model IGLOO
integrates, done independently, in the gas state IGLOO sampled. `check_reference.py` builds it and compares the
run against it.

## What the case is, physically

Each of the four DB streams (x = 0.02 m; y = 0.03, 0.045, 0.055, 0.07 m; z = 0.0025 m) injects 20 µm drops at
0.5 m/s. A drop stops within 0.53 mm (Stokes, τ_p = 1.09 ms) and then evaporates in place to burnout in about
80 ms. The whole evaporation is one IGLOO segment, since the drop does not cross a cell after it stops (sweep 0).

- **Sweep 0** sees the initial gas: uniform and quiescent. The first MOSE block converges in one iteration on
  the uniform initial condition, which is a fixed point, so the gas is unchanged.
- **Sweep k** sees the gas after 200·k + 1 pseudo-time iterations with the IGLOO sources.

A closed box fed a continuous mass source has **no steady state**. With `time-accurate = false` the gas after
N pseudo-iterations is not a physical state, and MI2's closed-domain lines print `injected 0.000000E+00` and
`rel.err 1.000E+00`; MI2 itself labels them pseudo-time integrals. Even the N2, which has no source, gains
9.5e-11 kg over the run (1.9e-6 of the 4.877e-5 kg in the box). The coupler's mass and energy balance is the
time-accurate twin's (`test/evap-box-ta`), not this case's.

What *is* well posed in steady mode is each IGLOO sweep: the droplet ODE in a given gas.

## The reference

Per stream, with d = (6m/(π ρ_l))^(1/3), ρ_l = 997 and c_l = 4184 from `INPUT/part-properties.dat`:

    dx/dt = u,   du/dt = 3π μ_g d (u_g − u)/m                              Stokes drag [CGW78]
    m c_l dT/dt = π d k_g Nu (T_g − T) + ṁ L_v(T),   Nu = 2 + 0.6 Re^½ Pr^⅓     Ranz–Marshall [RM52]
    dm/dt = ṁ = −π d k_g/(c_pg Le) Sh ln(1 + B_M),   Sh = 2 + 0.6 Re^½ Sc^⅓     CEM (Spalding) [Spa53]
    B_M = (Y_s − Y_∞)/(1 − Y_s),   Y_s = X_s M_v/(X_s M_v + (1 − X_s) M_g)
    X_s = X_s,eq − (2 L_K/d) β,   β = −ṁ Pr/(2π μ_g d),   2L_K = 2 μ_g √(2π T R_u/M_v)/(α_e Sc p)    LK [MHB98]
    X_s,eq = p_sat(T)/p,   p_sat = p_atm exp(−(L_v M_v/R_u)(1/T − 1/T_boil))   or   the Psat column

These are the equations of IGLOO's `rhsEvaporation`, `interphase`, `Lib_Drag`/`Lib_Heat`, `evaporation` and
`lkCorrection`, re-typed with IGLOO's constants: R_u 8314.46, boiling clamp X_s ≤ 1 − 1e-12, the Psat column
linear between integer-kelvin nodes as `tableValue`, and burnout at m = min(1e-6 m0, 1e-20 kg) (1e-20 kg here).
The inputs are T_boil 373.15, Le 1, Y_∞ 0 and α_e 1 from `input.ini`, and L_v, M_v where IGLOO takes them: the
`Lv` and `Mv` columns ATLAS GPB writes into the properties table for the `psat-vapour` pair (L_v at the droplet
temperature, linear between the integer-kelvin rows, 2.4383e6 J/kg at 300 K; M_v 18.015), or an
`[IGLOO-Properties]` value for a table without them (IGLOO refuses both). The p_sat source is the one IGLOO uses:
the Psat column when the properties table carries a non-zero one, else Clausius–Clapeyron with L_v(298.15 K).

The independent parts:
- an adaptive Dormand–Prince 5(4) integrator at rtol 1e-11, against IGLOO's SDIRK4;
- an exact solve of the LK fixed point by regula falsi on the bracket [ṁ_eq, 0], against IGLOO's damped Picard
  iteration, which IGLOO ends by bisection on the same bracket where 30 passes do not converge;
- face crossings located by bisection;
- the gas state rebuilt from `INPUT/gas-{phase.txt,thermo.dat,transport.dat}` as the gas solver builds it: FLINT
  R_u 8314.46261815324, integer-kelvin tables, Wilke mixing. For the initial gas (pure N2, T 350.000 K) this gives
  ρ 0.975415 kg/m³, μ 2.027840e-5 Pa s, k 0.0292834 W/m/K, γ 1.398777, Pr 0.7209.

Convergence of the reference: the lifetime moves by 6.5e-11 s from rtol 1e-10 to 1e-12.

Closed-form pieces (initial gas):

| quantity | Clausius–Clapeyron p_sat | Psat column (NASA9 `H2O(L)`/`H2O`) |
|---|---|---|
| p_sat(300 K) | 3101.5 Pa | 3533.6 Pa |
| wet bulb, VLE, Re 0 (Nu = Sh = 2): T_g − T = (L_v/c_pg) ln(1 + B_M(T)) | 300.624 K | 298.962 K |
| quasi-steady T with LK at d = 20 / 10 / 5 / 2 / 1.24 µm | 300.838 / 301.049 / 301.461 / 302.623 / 303.711 K | 299.175 / 299.385 / 299.800 / 300.972 / 302.068 K |
| K = 8 k_g (T_g − T_wb)/(ρ_l L_v); d0²/K | 4.71e-9 m²/s; 84.9 ms | 4.87e-9 m²/s; 82.1 ms |
| integrated lifetime to m = 1e-15 kg; T there | 83.447 ms; 303.687 K | 80.558 ms; 302.045 K |
| stop point | 0.0205331 m | 0.0205322 m |

The stop point with d and T frozen is x0 + u0 τ_p = 0.0205463 m: the drop evaporates while it stops.

## What `check_reference.py` checks

Two sweeps are checked by default (`--sweeps 0,last`): sweep 0 against the initial gas, and the last sweep
against `OUTPUT/gas-field.tec` (cell-wise constant, gas-order 1). Sweeps 1–3 are not checkable: their gas is not
written.

| id | check | tolerance and measured floor |
|---|---|---|
| E | per stream: the exit point (the stop point) | 1.5e-6 m (F12.6 rounding) |
| E | the exit T at burnout: paired with the exit mass when a final trajectory row carries it (`print-dtime` set) or the exit row's diameter (IGLOO writes `d_p` on its exit rows), else inside the reference's T over m ∈ [0.9, 1.0] × the burnout mass; the reference integrates on to 1e-3 of the burnout mass, since the step that crosses it can end well below it | the print budget |
| E | the exit time, when `out-time = on` prints it: at the exit mass as above, else inside the reference's times over the same mass range | 1e-6 relative |
| E | the number of scatter samples: one at each quantum up to the exit, ⌊frac(nφ) + t_exit/dt_q⌋ | exact at a paired exit mass; else between the bracket's ends |
| S | every scatter sample paired with the reference state at the same droplet mass: x, y, u, v, T | within 2× the budget: \|d(·)/dm\| × the E13.6 rounding of m, + the F12.6 rounding, + what a sample placed inside a solver step carries: IGLOO's columns are linear on the step (the position cubic), so paired by mass a column y is off by δ (m′y″/m″ − y′) and the position by δ y′, δ the sample's offset from its quantum (from t_k to the reference's time at the sample's mass, capped at C's 1e-6 s): the terms \|δ\| (\|m′y″/m″\| + \|y′\|) and \|δ y′\|. On a cell-wise gas 10× (measured with IGLOO 2133e13: 5.1 (CC) and 1.6 (Psat column), at a sample just after a face crossing, where IGLOO's 1e-8 m containment tolerance shifts the cell switch by ~1e-6 s). On the plateau, a sample within 5 thermal relaxation times after the drop temperature crosses a node of the Psat column has added to its T budget what a solver step across the node can leave, because the implicit step can go across a node whole: the distance at equal mass between the reference and the reference with the line of the interval it left continued, while the crossing step can last (taken as one relaxation time past the node), and the largest such distance decayed at the relaxation rate afterwards. Measured with IGLOO 2133e13 on sweep 0 with the gas at 350 K: a 2.5e-5 s step (half a relaxation time) across 300 K, where the slope of the column rises by 5.2 %, ends 0.24 relaxation times past the node with T 5.1e-6 K high (4.6 print budgets; allowance 7.6e-6 K there) and is back within twice the print 1.5 relaxation times past the node, six steps later; the allowance peaks at 1e-4 K one relaxation time past that node and is 2e-6 K at five. With the gas 1e-6 K colder the step before the node is halved and T stays within the print. Not covered at the default sampling (a sample every 2.83 ms): with `fsample-traj = 20000` (one every 14 µs; IGLOO cc125f0, sweep 0, first stream) the row fails at 13 samples on the approach to the plateau (d 19.2–18.4 µm) at 2.0–3.9×, as with IGLOO bdb4b91, which also failed it in the last ~60 µs of the life (45× at d 0.46 µm, 37× at 0.39 µm: its unconverged LK iteration); with cc125f0 those pass |
| Q | plateau samples against the closed-form quasi-steady (wet-bulb) balance at their diameter, beyond the reference's own lag | 1e-3 K; measured 2.7e-6 / 3.0e-6 K (sweep 0), 1.1e-4 / 8.7e-5 K (cell gas) |
| C | the scatter clock: IGLOO places sample k of stream n at t_k = dt_q (k − frac(nφ)), with dt_q = L_box/(v_inj · fsample-traj) = 2.830194 ms, inside the solver step that spans it; the reference reaches the sample's mass within 1e-6 s of t_k, before or after | 1e-6 s; the sample's mass is linear on the step, so it pairs a little early where m(t) is convex |
| G | last sweep: source field against the injected mass, momentum and energy of the consumed parcels, Σ ṁ (1, u0, c_l T0 + h_off + u0²/2) | 1e-9 relative; measured exact to the print, 4.000000000e-6 kg/s, 2.000000000e-6 N, −63.09831790 W |
| G | last sweep: Eulerian Σ n_p V = Σ ṅ_p t_life and Σ ρ_p V = Σ ṅ_p ∫m dt | 1e-4; measured 1.2e-5 / 2.7e-6 (CC), 7.2e-6 / 2.7e-6 (Psat column) |
| G | last sweep: the gas mass the coupler reports gained per species (`[MI2] closed-domain mass balance` in the run log) against Σ (ρ_i − ρ_i,0) V of the final gas field | the E13.6 print; measured exact to it (1.261670e-8 kg H2O, 9.476101e-11 kg N2). Bookkeeping, not conservation: see above |
| L | the run log (`--log`, else `RUN_DIR/run.log`), whatever the sweeps: no line of IGLOO's `LK-interface:` warning, which IGLOO prints once per material and sweep when the Langmuir–Knudsen fixed point of a parcel was not bracketed (that parcel's rate is the bracket's estimate, not the model's) | none; measured none (the four boxes) |

**Gating.** The last sweep's per-stream checks, the Eulerian integrals and the coupler's lines count only when
`gas-field.tec` was written at the MOSE iteration the sweep sampled (the last line of `gas-residual-history.dat`).
The case writes it at iteration 801, the last one (`sol-diter = 267`, see below), so the last sweep is gated. A run
that writes it earlier (`sol-diter = 200` writes it at iteration 800) gets these checks printed as information: with
the older gas the drop velocities are 1–2 % off, because the pseudo-time flow is still developing, and the gas holds
1.3e-3 less H2O than the coupler reports. The source-field identity does not depend on the gas and is always gated.

**What it cannot see.**
- The closures themselves against experiment. They are re-typed, so this verifies the integration path, how the
  pieces combine, the gas sampling and the tabulated p_sat, not whether CEM + LK is right. IGLOO's own suite
  pins the closures (`test_evaporation`, `test_interface_lk`, `mhb98-water`).
- The vapour already in the gas. IGLOO takes Y_∞ from `[IGLOO-Properties]` (0), not from the local gas, so the
  reference does too.
- The steady gas field, which has no physical meaning.
- The time history between injection and burnout, except through the droplet mass. The scatter rows carry no
  time, so the history is checked against the mass (S) and against the sample clock (C), which holds a drop that
  evaporates faster or more slowly than the reference to 1e-6 s at every sample. IGLOO writes a trajectory row per
  segment, not per step, and the evaporation is a single segment, so `print-dtime` adds no intermediate rows.
  Without `out-time` the lifetime is known to one scatter quantum (3.4 %) through the count, and on the last sweep
  to 1e-4 through the Eulerian integral.

**Output keys.** The case's `input.ini` sets three keys that change what is written, not what is computed; they make
the last sweep checkable and the lifetimes exact:
- `[IGLOO-General] out-time = on` adds the parcel time to the exit and trajectory rows;
- `[IGLOO-General] print-dtime = 1.0e-3` writes a trajectory row at the end of every segment past the next print
  time; here that is one final row at the exit mass, which pairs the exit T and time exactly;
- `[MOSE-IO] sol-diter = 267` writes `gas-field.tec` at iteration 801, the one the last sweep sampled: 267 divides
  the run's 801 MOSE iterations (`sol-diter = 1` gives the same file).

## Measured

All runs at 1 thread. The figures of this section, of the closed-form table above and of the falsification table
below are those of IGLOO before the burnout mass of 1e-20 kg, the samples placed inside the step and the L_v column
(constant L_v 2.462478e6 J/kg, burnout at 1e-15 kg, a sample at the end of the step past its quantum).

**At IGLOO cc125f0** (the four boxes, one thread): every check passes, E, S, Q, C, G and L, gated on sweep 0 of the
four and on the last sweep of evap-box (the other three write `gas-field.tec` at iteration 0: their last sweep is
informational). The exit T is the reference's at the exit mass, 9.94e-21 to 9.99e-21 kg: 331.273957 K against 331.273959 K
on sweep 0 and 329.77–330.31 K within 1e-6 K on the last sweep (evap-box, exit time within 7e-10 s, 8e-9 relative);
331.267–331.276 K within 1.6e-5 K on sweep 0 of the boxes without `print-dtime`. IGLOO's damped Picard iteration of the
LK interface ends unconverged on drops below 0.9 µm (the last 0.8 ms of the life); cc125f0 then bisects the fixed point.
IGLOO bdb4b91 kept the last Picard pass, which overstates the evaporation rate (4.4 times at 100 nm, 13.6 times at
27 nm): E failed on every stream of both sweeps, the drop burning out 4.3–4.7 µs early at 298.3–299.1 K.

**The stored baseline** (`OUTPUT/`: Psat column, output keys, scatter files from `hydra_MI2` `137de966`, IGLOO with the
scatter fix): every check passes, every check gated on both sweeps.
- Sweep 0: exit points 0.020532 m; exit T 302.048731 K, equal to the reference at the printed exit mass; exit time
  8.05592007e-2 s against 8.05591999e-2 s; samples 29/28/29/28 (N − 1 … N); S 0.99; Q 3.0e-6 K; clock at least 0.9 µs
  after t_k and within 0.215 dt_q.
- Last sweep: exit points and exit T equal to the print, exit times within 1.7e-9 s; samples 30/31/31/30 (N − 1 … N);
  clock at least 1.6 µs after t_k and within 0.229 dt_q; S 1.6; Q 8.7e-5 K; G source exact; Eulerian integrals 7.2e-6
  and 2.7e-6; the coupler's lines exact.

**The baselines without the output keys** (`OUTPUT.pre-outputkeys-2026-09-27/`, Psat column;
`OUTPUT.pre-psat-db-2026-09-27/`, Clausius–Clapeyron; both from a binary without the scatter fix): sweep 0 is gated and
passes (exit T 302.048731 K and 303.687689 K inside their burnout brackets, S 0.99, Q 3.0e-6 and 2.7e-6 K, clock within
0.215 and 0.225 dt_q); the last sweep is informational (gas at iteration 800), G source exact.

## Falsification (each check shown RED)

| run / reference | failing checks |
|---|---|
| the Clausius–Clapeyron run against a Psat-column reference (`--psat table --props <table>`) | E (exit T 303.69 vs [302.05, 302.15] K; one sample too many), S and Q (1.66 K), C (samples 2.9 ms early) |
| the Psat-column run against a Clausius–Clapeyron reference (`--psat cc`) | E (exit T 302.05 vs [303.69, 303.79] K), S and Q (1.66 K), C (samples 1.02 dt_q late) |
| the LK run against a VLE reference (`--interface VLE`) | E (burnout T 303.69 vs 300.62 K; sample counts), S, Q (2.33 K), C (0.7 ms early) |
| the solver run with `blowing = LK`, checked against the case's `input.ini` (`--ini`) | E (exit T −0.12 K), S, Q (0.13 K), C; with the exact gas also the last sweep's E (exit time +0.64 ms), S, Q and the Eulerian integrals (7.2e-3) |
| the solver run with `drag = Schiller-Naumann`, checked against the case's `input.ini` | E (stop point 0.020507 vs 0.020533 m), S (u), C (13.5 µs early); Q passes (the plateau is thermal) |
| the solver run with `time-end = 0.05` (parcels stopped alive) | E (no exit records), G source (Σ wdot 2.89e-6 vs 4.00e-6 kg/s) |
| the stored run with its scatter files from a binary without the scatter fix (`OUTPUT.pre-scatter-fix-2026-09-27/` swapped in) | E last sweep (streams 1 and 3: 31 and 32 samples against at most 30 and 31), C last sweep (samples up to 0.61 ms before their quantum); every other check passes |
| the Psat-column run with the Clausius–Clapeyron run's log (`--log`) | G coupler (N2 9.475749e-11 vs 9.476101e-11 kg; the H2O agrees: both runs put 1.261670e-8 kg into the gas) |
| the run (IGLOO cc125f0) with IGLOO's `LK-interface:` warning line planted in its log (`--log`) | L (1 line, quoted) |

## Use

    python3 check_reference.py [RUN_DIR] [--sweeps 0,last] [--psat auto|cc|table] [--interface auto|LK|VLE]
                               [--ini PATH] [--props PATH] [--log PATH]

`RUN_DIR` is a case directory after a run (`INPUT/`, `input.ini`, `OUTPUT/`; default: the current directory, so
`python3 check_reference.py` here checks the stored baseline); `--log` is the run's log (default `RUN_DIR/run.log`).
A kept baseline is checked by linking its `INPUT`, `OUTPUT`, `input.ini` and log into one directory. `--ini` and
`--props` fix the modelled case independently of the run. Exit 0 when every gated check passes. Python 3 with numpy;
one to two minutes for the two sweeps.
