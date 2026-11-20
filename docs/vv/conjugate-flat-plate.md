# Conjugate Flat Plate

**Cases:** `test/MF-flatplate-Tw` (`hydra-MF`) and `test/AF-flatplate-Tw` (`hydra-AF`)

A laminar boundary layer develops over a conducting slab whose underside is held at a fixed temperature. Heat flows up through the solid, across the fluid–solid interface, and into the gas, so the interface temperature $T_w(x)$ is **not** prescribed — it is the outcome of the balance between conduction in the plate and convection in the boundary layer.

That makes it the natural verification case for Hydra: the quantity being checked is produced by the coupling itself, not by either solver alone.

The case is run in two variants that differ only in the fluid solver, and each variant sweeps the solid conductivity over three decades.

---

## Why this case

Hydra contributes no physics of its own. The solvers integrate their own equations; Hydra makes them agree at the interface. This case is built so that the coupling is the only thing under test.

| Feature | Activated by |
|---------|--------------|
| Fluid–solid interface exchange | `face = connection` on both phases, `coupling-wall = true` |
| Two-solver orchestration | `[HYDRA-Parameters]` with `<solver>-iter = 1`, `FUSS-iter = 1` |
| Conjugate regime sweep | the solid `k` in `[GPB-Phase2]`, three values |
| Laminar wall heat transfer | $Re_L = 4.55\times10^5 < 5\times10^5$ |
| Solver-agnostic coupling | the same problem through `src/app/MF` and `src/app/AF` |

Within a variant the three cases differ **only** in the solid conductivity — same mesh, same numerics, same free stream — so the sweep isolates the conjugate response from every other effect. Between variants only the fluid solver changes, so the comparison isolates the coupling instantiation.

---

## The two variants

| | `MF-flatplate-Tw` | `AF-flatplate-Tw` |
|---|---|---|
| Executable | `hydra-MF` | `hydra-AF` |
| Fluid solver | MOSE | ARES |
| Fluid model | ideal gas, `species = air` | real fluid, CoolProp `Air` table |
| Inflow specified as | total conditions $p_0,\ T_0$ | mass flux $g$, static $T$ |
| Coupling steps | 200 000 | 400 000 |
| `INPUT/` | versioned with the case | versioned with the case; the table is shared in `common/` |
| Comparison script | `analytical.py` | `validate.py` |

Everything else — geometry, mesh, boundary conditions, solid material, free-stream state, numerics — is identical by construction.

!!! info "Why the real-fluid variant uses constant properties"
    The analytical reference assumes constant fluid properties. Each real-fluid case therefore ships `set_constant_transport.py` and `set_constant_cp.py`, which overwrite the viscosity, conductivity and $c_p$ blocks of the generated table with $\mu = 1.79352\times10^{-4}$ Pa·s, $k = 0.25737$ W/(m·K) and $c_p = 1004.5$ J/(kg·K), giving $Pr = 0.70$ exactly.

    The pressure/enthalpy grid and every other EOS variable — density, temperature, sound speed — are left untouched, so $T(p,h)$ is still resolved through the real-fluid table and the real-fluid code path is fully exercised. Only the properties that enter the correlation are pinned.

---

## Problem setup

```text
                            ambient  (p = 1 bar)
             ┌────────────────────────────────────────────┐
             │                                            │
   inflow    │        fluid — MOSE / ARES   0.30 m        │   ambient
             │          150 × 90 cells                    │   (outlet)
             ├══════════  connection interface  ══════════┤  y = 0.05 m
  symmetry   │        solid — FUSS          0.05 m        │   symmetry
             │          150 × 30 cells                    │
             └────────────────────────────────────────────┘  y = 0
                       T_B = 393 K   (prescribed)
            x = 0                                       x = 1 m
```

The plate is heated from below at $T_B = 393$ K and cooled from above by a $T_\infty = 293$ K stream, a driving difference of $\Delta T = 100$ K. The solid's side faces are adiabatic, so all the heat entering the bottom leaves through the interface.

### Boundary conditions

| Phase | Face | Location | Condition |
|---|---|---|---|
| fluid (Block 1) | 1 | $x = 0$ | `inlet` |
| fluid | 2 | $x = 1$ m | `outlet`, $p = 10^5$ Pa |
| fluid | 3 | $y = 0.05$ m | `connection` — coupled to the solid |
| fluid | 4 | $y = 0.35$ m | `outlet`, $p = 10^5$ Pa |
| solid (Block 2) | 1, 2 | $x = 0,\ 1$ m | `symmetry` (adiabatic) |
| solid | 3 | $y = 0$ | `wall`, $T_B = 393$ K |
| solid | 4 | $y = 0.05$ m | `connection` — coupled to the fluid |

The ideal-gas variant sets the inlet with the isentropic total conditions for $M = 0.2$ ($p_0 = 1.02828\times10^5$ Pa, $T_0 = 295.34$ K); the real-fluid variant sets the equivalent mass flux directly ($g = 81.6055$ kg/(m²·s), $T = 293$ K). Both give the same free stream.

### Free-stream and material state

| Quantity | Symbol | Value | Unit |
|---|---|---|---|
| Static temperature | $T_\infty$ | 293.0 | K |
| Static pressure | $p_\infty$ | $1.0\times10^{5}$ | Pa |
| Density | $\rho_\infty$ | 1.18919 | kg/m³ |
| Velocity | $U_\infty$ | 68.623 | m/s |
| Mach number | $M$ | 0.20 | – |
| Dynamic viscosity | $\mu$ | $1.79353\times10^{-4}$ | Pa·s |
| Specific heat | $c_p$ | 1004.5 | J/(kg·K) |
| Prandtl number | $Pr$ | 0.70 | – |
| Fluid conductivity | $k_f = \mu c_p / Pr$ | 0.25737 | W/(m·K) |
| Plate-length Reynolds | $Re_L$ | $4.55\times10^{5}$ | – |
| Backside temperature | $T_B$ | 393.0 | K |
| Plate thickness | $t$ | 0.05 | m |
| Plate length | $L$ | 1.0 | m |

!!! note "Why the viscosity looks unusual"
    $\mu = 1.79353\times10^{-4}$ Pa·s is about ten times the physical value for air. It is chosen so that $Re_L = \rho_\infty U_\infty L / \mu = 4.55\times10^5$ on a 1 m plate — just below transition — while keeping the geometry, the Mach number and the temperatures at round values. The correlation depends only on $Re_x$, $Pr$ and $k_f t/(k_s x)$, so this is a consistent fully laminar setup rather than a physical air property.

### The conductivity sweep

The solid phase is the same FUSS material in both variants; only `[GPB-Phase2] k` changes.

| Case suffix | `k` [W/(m·K)] | $k_s/k_f$ | $Bi_f(L/2)$ | Regime |
|---|---:|---:|---:|---|
| `ks_kf_40` | 10.2746 | 39.9 | 0.704 | conduction-limited — plate cools strongly |
| `ks_kf_400` | 102.746 | 399.2 | 0.070 | intermediate |
| `ks_kf_4000` | 1027.46 | 3992.1 | 0.007 | near-isothermal plate |

### Grid

A 2-D plane mesh (`type = 2Dplane`, one cell of width 0.001 m in $z$) from the `[GRIB-*]` blocks:

| Block | Phase | Cells | Extent | Stretching |
|---|---|---:|---|---|
| Block 1 | fluid | $150 \times 90$ = 13 500 | 1.0 × 0.30 m | uniform in $x$; ratio 1.05 growing away from the interface |
| Block 2 | solid | $150 \times 30$ = 4 500 | 1.0 × 0.05 m | uniform in $x$; ratio 1.04 clustering towards the interface |

Both blocks share the same 150 streamwise cells, so the interface is point-matched and no interpolation enters the exchange. The first fluid cell sits at $y^+ = 0.8$–$4.4$ along the plate.

---

## Analytical reference

Lindstedt & Karvinen give a closed-form interface temperature for a laminar forced-convection boundary layer over a plate of finite conductivity and thickness with a prescribed backside temperature. With the local conjugate Biot-type parameter

$$
Bi_f(x) \;=\; \frac{0.332}{0.5}\, Re_x^{1/2}\, Pr^{1/3}\, \frac{k_f\, t}{k_s\, x},
$$

the non-dimensional interface temperature is

$$
\theta_w \;=\; \frac{T_B - T_w}{T_B - T_\infty} \;=\; \frac{Bi_f\,\bigl(0.5 + Bi_f\bigr)}{\bigl(1 + Bi_f\bigr)^{2}},
\qquad
T_w(x) \;=\; T_B - (T_B - T_\infty)\,\theta_w .
$$

The factor $0.332\,Re_x^{1/2}Pr^{1/3}$ is the local laminar Nusselt number for an isothermal plate, so $Bi_f$ compares the convective resistance of the boundary layer with the conductive resistance of the slab. The two limits are the physical sanity check:

- $k_s \to \infty$ ⟹ $Bi_f \to 0$ ⟹ $T_w \to T_B$ — an infinitely conductive plate stays isothermal at its backside temperature.
- $k_s \to 0$ ⟹ $Bi_f \to \infty$ ⟹ $T_w \to T_\infty$ — a non-conducting plate lets the surface float to the gas temperature.

The sweep walks between these limits, and $T_w$ rising towards $T_B$ as $k_s$ grows is what the coupling must reproduce.

!!! info "Source material"
    Both reference papers are versioned with the case in `test/MF-flatplate-Tw/reference/`: `Lindstedt Correlation.pdf` and `Steady-State-Coupled-Fluid-Solid.pdf`. The correlation is implemented in `analytical.py::tw_laminar_forced`.

---

## Configuration

The coupling is configured identically in both variants, with `MOSE-iter` or `ARES-iter` naming the fluid side:

```ini
[HYDRA-Parameters]
N-coupling    = 200000   ; 400000 in the real-fluid variant
FUSS-iter     = 1        ; solid sub-iterations per coupling step
MOSE-iter     = 1        ; fluid sub-iterations per coupling step
time-accurate = false    ; march to the steady conjugate state

[HYDRA-IO]
coupling-diter     = 20000
coupling-wall      = true    ; write OUTPUT/hydra-wall.tec
coupling-overwrite = false
```

With `FUSS-iter = MOSE-iter = 1`, one coupling step is a single solid RK3 step followed by a single fluid RK3 step, with an interface exchange before each. This is the tightest available coupling and the configuration the case is verified in.

Both solvers use the same numerics — RK3, `cfl = vnn = 0.55`, HLLC, MUSCL with a minmod limiter, and implicit residual smoothing at $\beta = 0.5$ — and the solid phase is defined as:

```ini
[GPB-Phase2]
type     = solid
name     = plate
material = mat
cp  = 700
rho = 10000
k   = 10.2746            ; 102.746 / 1027.46 in the other two cases
```

---

## Running the cases

=== "Ideal gas (`hydra-MF`)"

    The ideal-gas cases ship a complete `INPUT/` and `MESH/`, so no pre-processing is needed:

    ```bash
    cd test/MF-flatplate-Tw/ks_kf_40
    figaro hydra-MF solve -p 8       # 8 OpenMP threads
    ```

    Repeat for `ks_kf_400` and `ks_kf_4000`, then compare:

    ```bash
    cd test/MF-flatplate-Tw
    python3 analytical.py
    ```

=== "Real fluid (`hydra-AF`)"

    The real-fluid cases ship a complete `INPUT/` too. The real-fluid table, already flattened to the constant properties, is shared by the three cases through links into `common/`:

    ```bash
    cd test/AF-flatplate-Tw/ks_kf_40
    figaro hydra-AF solve -p 8       # 8 OpenMP threads
    ```

    Repeat for `ks_kf_400` and `ks_kf_4000`, then compare:

    ```bash
    cd test/AF-flatplate-Tw
    python3 validate.py                      # all three cases
    cd ks_kf_40 && python3 validate.py       # one case, from its own OUTPUT/
    ```

    To rebuild `INPUT/` (after a change to ATLAS or to the deck), work in a copy of the suite and remove the links first, since ATLAS would otherwise write through them into `common/`:

    ```bash
    cp -r test/AF-flatplate-Tw /tmp/AF-flatplate-Tw && cd /tmp/AF-flatplate-Tw/ks_kf_40
    find INPUT -type l -delete
    figaro ATLAS GPB ICB BCB             # 1. real-fluid table, IC, BC from input.ini
    python3 set_constant_transport.py    # 2. flatten the table
    python3 set_constant_cp.py
    ```

    Step 2 must follow step 1: ATLAS regenerates the table from scratch, overwriting any flattening already applied.

Each run writes `OUTPUT/hydra-wall.tec`. The suite-level scripts read the frozen solutions in `reference/`, so they work without a prior solve; point them at your own `OUTPUT/hydra-wall.tec` to check a fresh run.

!!! note "One mesh"
    Both variants run on the same mesh: `MESH/mesh.tec` is one file, identical in the six cases. The mesh table above gives the blocks it was generated with.

---

## Results

<figure>

  --8<-- "docs/vv/images/conjugate-flat-plate-profiles.svg"

</figure>

Both variants track the correlation over the whole plate and across all three decades of conductivity ratio, and the conjugate trend is captured correctly. At $k_s/k_f = 40$ the surface is dragged far below the backside temperature — 314.7 K at the leading edge, still some 20 K below $T_B$ at the trailing edge. At $k_s/k_f = 4000$ the plate is effectively isothermal, spanning only 390.4–392.8 K over its whole length.

| Variant | $k_s/k_f$ | $T_w(L/2)$ correlation | $T_w(L/2)$ Hydra | max $\lvert e \rvert$ | RMSE | max $\lvert e \rvert / \Delta T$ |
|---|---:|---:|---:|---:|---:|---:|
| `hydra-MF` | 40 | 363.796 K | 365.780 K | 2.009 K | 1.830 K | 2.0 % |
| `hydra-MF` | 400 | 389.493 K | 389.852 K | 0.678 K | 0.362 K | 0.68 % |
| `hydra-MF` | 4000 | 392.648 K | 392.684 K | 0.176 K | 0.047 K | 0.18 % |
| `hydra-AF` | 40 | 363.796 K | 365.427 K | 1.651 K | 1.497 K | 1.7 % |
| `hydra-AF` | 400 | 389.493 K | 389.804 K | 1.126 K | 0.357 K | 1.1 % |
| `hydra-AF` | 4000 | 392.648 K | 392.679 K | 0.259 K | 0.053 K | 0.26 % |

Errors are measured against $\Delta T = T_B - T_\infty = 100$ K and exclude $Re_x < 9\times10^{3}$ — the first few cells, where the correlation itself is singular.

<figure>

  --8<-- "docs/vv/images/conjugate-flat-plate-error.svg"

</figure>

The deviation is largest in the strongly conjugate case and falls by roughly an order of magnitude per decade of $k_s/k_f$, which is the expected ordering: as the plate approaches isothermal, both the correlation and the simulation converge on the trivial answer $T_w = T_B$, so there is less for either to get wrong. The residual $\approx 2$ K offset at $k_s/k_f = 40$ is the regime where the correlation's own assumptions — a 1-D conduction path through the slab, no streamwise conduction — are weakest. Routing the fluid state through a tabulated equation of state does not degrade the result: the two variants sit on top of one another.

### Interface continuity

`OUTPUT/hydra-wall.tec` contains two zones: the solid-side interface (`plate-B1F4`, written first) and the fluid-side interface (`gas-B1F3` or `fluid-B1F3`). Each side's wall temperature is computed independently by its own solver's coupled-wall routine, so their difference measures how tightly the coupling has converged:

| $k_s/k_f$ | `hydra-MF` | `hydra-AF` |
|---:|---:|---:|
| 40 | $2.8\times10^{-8}$ K | $1.1\times10^{-13}$ K |
| 400 | $6.9\times10^{-10}$ K | $1.1\times10^{-13}$ K |
| 4000 | $3.5\times10^{-11}$ K | $5.7\times10^{-14}$ K |

Temperature continuity across the interface holds to round-off. This is the coupling's core invariant, and it is satisfied far below any physical or discretization error.

!!! warning "Reading the solid-side zone"
    Only `Tw` and `qw` are meaningful in the solid zone. Its `y+`, `tauX`, `tauY`, `tauZ` and `pw` columns are filled with the placeholder constant $1.78625\times10^{-7}$, since those quantities have no solid-side meaning. The interface heat flux `qw` is the shared coupled flux and is written identically to both zones, so it is not an independent continuity check.

### Cross-solver consistency

This is the check only Hydra can perform. The two variants solve the same physical problem with the same solid solver and the same constant properties, differing only in which fluid solver and which coupled executable is used. Any difference is attributable to the fluid solver and its coupling instantiation, not to the physics.

<figure>

  --8<-- "docs/vv/images/conjugate-flat-plate-cross-solver.svg"

</figure>

| $k_s/k_f$ | max $\lvert \Delta T_w \rvert$ | RMS $\Delta T_w$ | max $\lvert \Delta T_w \rvert / \Delta T$ |
|---:|---:|---:|---:|
| 40 | 0.668 K | 0.342 K | 0.67 % |
| 400 | 0.693 K | 0.099 K | 0.69 % |
| 4000 | 0.102 K | 0.013 K | 0.10 % |

ARES and MOSE agree on the interface temperature to better than 0.7 K out of a 100 K driving difference, everywhere outside the leading-edge cells (shaded; the $k_s/k_f = 40$ curve runs off-scale there). The difference is a smooth, small, negative offset — the real-fluid path runs marginally cooler — rather than a divergence that grows along the plate, which is what a coupling inconsistency would look like.

That residual offset is expected: ARES inverts a $400\times400$ tabulated EOS where MOSE evaluates closed-form ideal-gas relations, so the two fluid states differ at the level of the table's interpolation error even with the transport properties pinned identically.

---

## What this verifies

- The **fluid–solid interface exchange** transports the correct temperature and flux: the conjugate response follows an independent analytical solution across three decades of conductivity ratio.
- The **coupling loop converges to a consistent steady state** — the two sides of the interface agree on the wall temperature to $10^{-8}$ K or better.
- The **coupling is solver-agnostic**: `src/app/MF` and `src/app/AF` reproduce the same conjugate solution to within 0.7 % of the driving temperature difference.
- The **real-fluid path survives coupling** — tabulated EOS inversion at the interface does not degrade the solution relative to the correlation.
- The **limiting behaviour** is right: increasing $k_s$ drives $T_w \to T_B$, decreasing it drives $T_w \to T_\infty$.

## Acceptance criteria

A framework change should leave all of the following intact:

| Check | Threshold |
|---|---|
| Interface temperature continuity | $< 10^{-6}$ K |
| Deviation from the correlation, $k_s/k_f = 40$ | $< 3$ % of $\Delta T$ |
| Deviation from the correlation, $k_s/k_f \ge 400$ | $< 1.5$ % of $\Delta T$ |
| `hydra-AF` vs. `hydra-MF` interface temperature | $< 1$ % of $\Delta T$ |
| Limiting behaviour | $T_w \to T_B$ as $k_s/k_f$ grows |

These are the observed margins rounded up, not tolerances derived from a formal uncertainty analysis. See [Testing](../development/testing.md) for which checks to run against which kind of change.

## Limitations

- The correlation is **singular at the leading edge** ($Bi_f \to \infty$ as $x \to 0$), so the first few cells are excluded from the error metrics. This is a property of the reference, not of the solution.
- The comparison is **steady-state only**. Transient conjugate coupling (`time-accurate = true`) is not covered.
- The case is **laminar by construction** ($Re_L < 5\times10^5$). Turbulent conjugate heat transfer is exercised in the solver-level suites.
- The real-fluid table is deliberately **flattened to constant properties**, so this case does not test conjugate heat transfer with strong real-fluid property variation. For that regime see the [ARES HTD case](https://open-hydra.github.io/ARES/vv/htd/), which is uncoupled but uses the full table.
- The runs are **serial-equivalent**: the wall output is assembled on rank 0, and the case is small enough that the MPI decomposition of the interface is not stressed. Parallel coupling behaviour is covered in [Parallel Execution](../user/parallel.md).

## Reproducing the figures

Every figure and number on this page is derived from the frozen wall solutions in each suite's `reference/` directory, so no solver run is required:

```bash
python3 test/make-vv-figures.py     # regenerates docs/vv/images/*.svg and prints the tables
```
