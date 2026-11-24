# Multiphase Coupling (`hydra-MI2`)

`hydra-MI2` couples a gas with one or two dispersed condensed phases (drops or particles) in one executable:

| Solver | Solves | Phase it takes | Documentation |
|--------|--------|----------------|---------------|
| **MOSE** | the gas (Euler/Navier-Stokes, ideal-gas mixture) | the `ideal-gas` phase | [MOSE User Guide](https://open-hydra.github.io/MOSE/user/) |
| **ICE** | an Eulerian condensed phase (moment closures per family) | the condensed phase named by `[ICE-Parameters] phase` | [ICE User Guide](https://open-hydra.github.io/ICE/user/) |
| **IGLOO** | a Lagrangian condensed phase (parcels tracked through the gas) | the condensed phase named by `[IGLOO-General] phase` | [IGLOO User Guide](https://open-hydra.github.io/IGLOO/user/) |

Either condensed solver can be switched off: a case runs MOSE + ICE, MOSE + IGLOO, or all three. The coupler hands
each condensed solver the gas state, takes back what the condensed phases do to the gas (momentum, energy and, when a
material evaporates or burns, mass and species), and moves the smoke of a burning material into an ICE family. Each
solver keeps its own models, numerics and `[<solver>-*]` sections; this page documents what `hydra-MI2` itself reads,
checks and prints.

---

## Build

```bash
figaro build hydra-MI2 --compilers=intel --use-openmp
figaro compile hydra-MI2          # incremental rebuild of an existing build/MI2
```

The executable is written to `bin/hydra_MI2`. `figaro build hydra-MI2` also checks out the two condensed-phase
solvers, `src/lib/ICE` and `src/lib/IGLOO`. ICE is fetched over SSH and is marked `update = none` in `.gitmodules`,
so a recursive clone skips it: the build checks it out on request, and needs access to its repository.

!!! warning "OpenMP only"
    `hydra-MI2` refuses to configure with `--use-mpi`. The configure step stops with:

    ```text
    hydra-MI2 does not support MPI. MOSE partitions over blocks but IGLOO requires the full gas field on every rank. Rebuild without --use-mpi.
    ```

    MOSE distributes its blocks over the MPI ranks, while IGLOO tracks each parcel through the whole gas field and
    needs that field, identical, on every rank. A block-partitioned coupling would hand IGLOO a field with holes. Run
    `hydra-MI2` on one node with OpenMP threads (`figaro hydra-MI2 -p <n> solve`).

---

## Run

```bash
cd my_case
figaro hydra-MI2 -p 8 solve
```

`figaro` lists `INPUT/*phase.txt` into `filelist.txt`, copies `bin/hydra_MI2` to the case's `bin/hydra-MI2` when the
master is newer, and runs it in the case directory (see [Running Cases](using.md)). To run the executable directly,
write `filelist.txt` first (`ls INPUT/*phase.txt > filelist.txt`): without it the run prints
`Error opening filelist.txt` and ends at once.

### Case directory

A case is prepared with ATLAS from the `[GRIB-*]`, `[GPB-*]`, `[BCB-*]` and `[ICB-*]` sections of its `input.ini`.
Each phase has its own file set, named after the phase (`<name>-`):

```text
my_case/
├── input.ini                 # [HYDRA-*], [MI2-Coupling], [MOSE-*], [ICE-*], [IGLOO-*] and the ATLAS sections
├── MESH/
└── INPUT/
    ├── gas-phase.txt         # the ideal-gas phase: MOSE
    ├── gas-thermo.dat  gas-transport.dat  gas-bc.txt  gas-ic.tec  ...
    ├── partE-phase.txt       # a condensed phase solved by ICE
    ├── partE-bc.txt  partE-ic.tec  partE-properties.dat
    ├── partE-products.txt    # what a converting material becomes (ATLAS GPB), read by hydra-MI2
    ├── partL-phase.txt       # a condensed phase solved by IGLOO
    ├── partL-bc.txt  partL-ic.tec  partL-properties.dat
    ├── partL-products.txt
    └── partL-models.txt      # the constants of the burning and breakup models (ATLAS GPB), read by IGLOO
```

The names (`gas`, `partE`, `partL`) are those of the `[GPB-Phase*] name` keys. The run writes `OUTPUT/`, each
solver's files carrying its phase prefix.

---

## Phase routing

Before any solver is set up, `hydra-MI2` assigns every `INPUT/*phase.txt` listed in `filelist.txt` to a solver,
from line 1 of the file (the type word ATLAS writes):

- `ideal-gas` → MOSE. Exactly one such phase is required.
- a word containing `dispersed` (`condensed-dispersed`, `liquid-dispersed`, `solid-dispersed`) → ICE or IGLOO, by
  the phase key of each solver:

    ```ini
    [ICE-Parameters]
    phase = partE     ; the [GPB-Phase*] name of the phase ICE solves

    [IGLOO-General]
    phase = partL     ; the [GPB-Phase*] name of the phase IGLOO solves
    ```

- any other word containing `solid` (a bulk solid) stops the run: `hydra-MI2` has no solid solver (use `hydra-MF`).

A solver takes a phase only when it is enabled (its `-iter` key is present, see [Coupling schedule](#coupling-schedule)).
Line 1 may also carry ATLAS's `modeling=eulerian` or `modeling=lagrangian` token; when present it must agree with the
solver that claims the phase (`eulerian` for ICE, `lagrangian` for IGLOO).

With a phase key absent, the solver takes the only condensed phase whose token does not contradict it, and prints a
`DEPRECATED` line asking for the key; with more or fewer than one candidate the run stops. When both condensed solvers
run, each needs its own phase: one phase cannot feed both (its mass would reach the gas twice).

### Files each solver reads

| Solver | Files of `INPUT/` (prefix `<name>-`) | Checked by `hydra-MI2` before setup |
|--------|--------------------------------------|-------------------------------------|
| MOSE | the gas phase's files | — |
| ICE | `phase.txt`, `bc.txt`, `ic.tec` (or `ic.szplt`, `ic.plt`), `properties.dat`, `models.txt` | `bc.txt` required; no `properties.dat`: a note (ICE then uses the constant `[ICE-Physics]` `density` and `specific-heat`, with one material only); no initial condition: a warning (needed unless `[ICE-Parameters] newrun = false`) |
| IGLOO | `phase.txt`, `bc.txt`, `properties.dat`, `models.txt` | `bc.txt` and `properties.dat` required |
| `hydra-MI2` | `products.txt` of each condensed phase | see [Evaporation and combustion](#evaporation-and-combustion) |

`models.txt` exists only when a material has constants of a burning or a breakup model; the solvers read it,
`hydra-MI2` does not. Its format and the keys it carries are in the solvers' input pages
([ICE](https://open-hydra.github.io/ICE/user/input/), [IGLOO](https://open-hydra.github.io/IGLOO/user/input/)) and
in [ATLAS GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/condensed-solid/#the-constants-of-the-models-namemodelstxt).
Breakup is IGLOO's: ICE has no breakup model and ignores the `breakup=` word of a material line and a `breakup` row,
one warning each, and the run goes on; ATLAS writes both for any phase whose `[GPB-Phase*]` gives breakup keys, so
give them to the phase IGLOO solves. The oxide cap, IGLOO's too, is refused by ICE (see [Oxide cap](#oxide-cap-igloo)).

### Naming rules

A phase file is `INPUT/<name>-phase.txt`, or `INPUT/phase.txt` for an unnamed phase. No name may contain `-` (the
prefix is split at the first dash), and a name must fit the solver's prefix buffer. When a case has more than one
phase file, every phase must be named, names must be distinct and no name may be a substring of another (ATLAS
associates phases by substring).

### The routing table

The run prints the assignment it made:

```text
 HYDRA-MI2 phase routing
 -  INPUT/gas-phase.txt  (ideal-gas)  --> MOSE   prefix "gas-"
 -  INPUT/partE-phase.txt  (condensed-dispersed, modeling=eulerian)  --> ICE  prefix "partE-"  (from its phase key)
 -  INPUT/partL-phase.txt  (condensed-dispersed, modeling=lagrangian)  --> IGLOO  prefix "partL-"  (from its phase key)
```

A phase no enabled solver names is listed as `UNUSED (no enabled solver names it)`; a phase taken without its key
as `LEGACY fallback`.

---

## Coupling schedule

`hydra-MI2` reads the common `[HYDRA-Parameters]` and `[HYDRA-IO]` keys (see [Hydra Input](input.md)) and four keys of
its own in `[HYDRA-Parameters]`. Lists hold one value per **cycle**; the number of cycles is the number of values of
`MOSE-iter`.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `N-coupling` | `1000000000` | Coupling steps of each cycle (one value per cycle) |
| `MOSE-iter` | `1` | MOSE iterations per coupling step (one value per cycle; `0` = MOSE not advanced) |
| `ICE-iter` | `0` | ICE iterations per coupling step (one value per cycle). **Absent: ICE is not run.** `0` in a cycle: ICE is not advanced |
| `IGLOO-iter` | `0` | IGLOO mode per cycle: `0` off, `1` source coupling, `2` Eulerian coupling. **Absent: IGLOO is not run.** Values below 1 are taken as `0`, above 1 as `2` |
| `relaxation` | `false` | Under-relax IGLOO's source on the gas |
| `relaxation-coeff` | `0.5` | Under-relaxation coefficient: one value, or one per cycle |
| `time-accurate` | `false` | Labels the end-of-run balances only (see [Steady and time-accurate runs](#steady-and-time-accurate-runs)) |

`ICE-iter` and `IGLOO-iter` must hold as many values as `MOSE-iter`; with `relaxation = true`, `relaxation-coeff` holds
one value or one per cycle. The keys must be in `[HYDRA-Parameters]`, whose name is case-sensitive.

### One coupling step

Each coupling step runs, in this order:

1. **MOSE**: `MOSE-iter` iterations, with the sources the condensed solvers produced last applied to the gas at every
   Runge-Kutta stage.
2. **IGLOO** (when `IGLOO-iter` ≠ 0): the gas state is handed to IGLOO, which tracks its parcels through the whole
   field once and produces new sources.
3. **ICE** (when `ICE-iter` ≠ 0): the gas state is handed to ICE, which makes `ICE-iter` iterations and produces new
   sources.

A cycle ends after `N-coupling` coupling steps. A solver's own stopping condition (`res-threshold`, ...) ends only
the current block of its iterations, never the cycle.

The run ends as a run of each solver alone ends. The last iteration the schedule gives a solver (the last of its
cycle's last coupling step, no later cycle advancing it) is its last: MI2 sets the solver's own iteration limit there.
After the last cycle each solver writes its end-of-run output once: the residual and probe rows of its last iteration,
the solution due there (with its output number), then its final solution, `OUTPUT/<prefix>field.tec`, whatever
`sol-diter`. Each solver's output number is then the number of outputs it wrote, and a run whose iterations are no
multiple of `sol-diter` can be restarted from its end. A solver that restarted and makes no iteration writes no
end-of-run output: its files are those it restarted from. A re-read of `input.ini` that takes a solver's later
iterations away after its last one ran (as an ordinary iteration) leaves that iteration's row written twice, and its
solution counted twice when one was due there: the end-of-run output repeats it (not detected).

### Frozen sources

A solver that has produced sources keeps them in force when a later cycle does not advance it (`ICE-iter = 0` or
`IGLOO-iter = 0`): the gas goes on receiving them. The run says so at the first coupling step of such a cycle, then
each time the cumulative MOSE iteration count crosses a multiple of `[MOSE-IO] shell-diter`:

```text
 [MI2] IGLOO source FROZEN this cycle (IGLOO-iter = 0): the gas keeps the sources produced in cycle 1, coupling step 3  [cycle 2, step 1, MOSE iter 150]
```

After an ICE [restart](#restart), ICE's sources not yet renewed are those read at the restart, and the notice says so.

### Under-relaxation

With `relaxation = true` and a coefficient \(r\) strictly between 0 and 1, IGLOO's source \(S\) on the gas is blended
with the one applied at the previous Runge-Kutta stage, \(S \leftarrow c\,S + (1-c)\,S_{\text{prev}}\). Over the
first \(1-r\) of a MOSE block of \(N\) iterations, at iteration \(i\),

$$
c = x\left(\frac{1}{1-r} - e^{-2x/r}\right),\qquad x = \frac{i}{N},
$$

clipped to \([0, 1]\); \(c = r\) when \(N = 1\), and \(c = 1\) otherwise. ICE's sources are never relaxed.

### Steady and time-accurate runs

`hydra-MI2` has no time mode of its own: MOSE and ICE each take theirs from `time-accurate` in `[MOSE-Numerics]` and
`[ICE-Numerics]`. `[HYDRA-Parameters] time-accurate` only labels the end-of-run balances (with it off they are
pseudo-time integrals, not physical amounts).

- A cycle that advances MOSE and ICE together (`MOSE-iter > 0` and `ICE-iter > 0`) is refused while either is
  time-accurate: the two solvers do not share a time step, so their physical times would part. Advance them in
  separate cycles, as `test/ice-burn-box` does:

    ```ini
    [HYDRA-Parameters]
    N-coupling    = 1 200
    MOSE-iter     = 0 1
    ICE-iter      = 1 0
    relaxation    = false
    time-accurate = true
    ```

- When MOSE or ICE steps in pseudo-time, what the coupler moves between ICE and the gas or IGLOO (ICE's burn-out
  removals, IGLOO's smoke) is handed as **rates**, so that a converged steady field does not depend on the local time
  steps. When both are time-accurate it is handed as **amounts**.
- A time-accurate ICE block of several steps hands the gas the mean of its steps' sources, each weighted by its time
  step.

A steady three-solver schedule (`test/JPL-3solver-20micron`):

```ini
[HYDRA-Parameters]
N-coupling    = 50
MOSE-iter     = 500
IGLOO-iter    = 2
ICE-iter      = 250
relaxation    = false
time-accurate = false
```

### MOSE's grid levels

`[MOSE-Multigrid] levels` above 1 needs ATLAS's coarse boundary files (`[ATLAS-Parameters] MG-levels`) and
`level1-iter`, MOSE's iteration limit of the finest grid when `levels` is above 1 (the value `iter-threshold` would
have): without it the run stops after MOSE's set-up. IGLOO and ICE always see the finest grid.

- Start on the coarse grids (`level<m>-iter` > 0 for an m > 1; steady runs): MOSE steps its coarsest level first and
  works down to the finest, as it does alone (a level until its `level<m>-iter`, or until its density residual has
  fallen by `res-threshold`). `MOSE-iter` counts the iterations of every level. Before each hand-off, report and the
  end-of-run output the coupler fills the finest grid from the level MOSE steps: MOSE's prolongation level by level,
  then the finest grid's derived variables, local time step, ghost cells, fluxes and boundary records (the inflow mass
  flow IGLOO takes), with no time step and without the condensed phases' sources. On the level MOSE steps the coupler
  adds the finest grid's sources summed over each coarse cell's sub-cells: IGLOO's rates summed; ICE's per-volume
  rates times the sub-cells' volumes (on a curved mesh a coarse cell's own volume is not its sub-cells'); IGLOO's
  Eulerian field (`IGLOO-iter = 2`) as density and number per volume, velocity and temperature per mass, evaluated with
  the coarse gas. The level's totals are the finest grid's (for a material whose cp varies with T, those of
  density x temperature, not of its enthalpy). A run that ends on a coarse level writes the finest grid, filled, under
  the name MOSE gives that level's output (`gas-field-level<m>.tec`; its diagnostic holds the filled grid's residual
  without the condensed sources). Refused at set-up: a time-accurate run (MOSE does not carry the physical time from
  one level to the next), and `relaxation = true` (it blends two steps' sources cell by cell of the stepped level).
- Multigrid cycle (`cycle = v` or `w`, every `level<m>-iter` = 0 for m > 1; MOSE has it for steady runs without
  chemistry): MOSE steps its finest grid from the first iteration. A coarse level of the cycle takes IGLOO's and ICE's
  sources summed over its sub-cells in the residual MOSE restricts to it, and the converged gas is the one of a run
  without the cycle. With a coarse start, the cycle's levels below the one MOSE steps get nothing more either.

### Other refusals

- Multigrid: `[HYDRA-Multigrid] levels` above 1 stops the run (IGLOO has no grid levels).
- `[MOSE-Multigrid]`: a start on the coarse grids in a time-accurate run or with `relaxation = true`, and
  `levels` above 1 without `level1-iter` ([MOSE's grid levels](#moses-grid-levels)).
- `[ICE-Multigrid] levels` above 1 stops the run: the coupler hands the gas to ICE's finest grid and takes ICE's source
  from it, while ICE steps a coarse level on that level's own grid (and needs ATLAS's `<name>-bc<level>.txt` for it).
- `IGLOO-iter = 2` in any cycle while an IGLOO material evaporates or burns: the Eulerian coupling hands the gas
  momentum and energy only. Use `IGLOO-iter = 1` in every cycle that runs IGLOO.

### Reading `input.ini` again

Every `[HYDRA-IO] ini-diter` coupling steps (default 1000) the values of the `[HYDRA-*]` keys are read again, so a
long run can be steered without a restart; the number of values of a list is not. The re-read checks the registry
(required keys and the allowed range of each scalar key): `N-coupling` and `time-accurate` must then be present (`[ERROR] Required parameter not set: <key>`); `FUSS-iter` and
the multigrid levels are not required. `[MI2-Coupling]` and the phase keys are read once, at start.

### Restart

`newrun = false` in `[MOSE-Parameters]` and `[ICE-Parameters]` starts each solver from the files a run left in
`OUTPUT/`: MOSE from `gas-field.tec`, ICE from `<prefix>field.tec`, each with its iteration and output counters (with
the default `sol-overwrite = true` these are the solver's last output, also of a run that was killed). The
schedule of the restarted run is what remains of the run: a run of four coupling steps, and a run of two restarted for
two more, end with the same files, the output numbers and the residual histories included. ICE sets the components a mesh holds at zero (`w` on a 2-D mesh) to zero at every
set-up, so the round-off a run leaves there is not carried over.

What ICE's last coupling step hands the gas is not in its field file: it is the source of the stages of ICE's last step
(averaged over its last block when ICE is time-accurate), not a function of the state ICE writes; nor is what ICE's next
step starts from in a coupled run, the relaxation times that bound its steps and a smoke family's birth rate. MI2
writes them to `OUTPUT/<prefix>coupling.tec` (Tecplot ASCII, 17 digits; per family the source, the heat, the smoke's
birth, the removals' store, the relaxation time and, where they exist, a smoke family's birth rate and the removals'
rate; ICE's iteration as `CYCLE`, MOSE's as `MI2_MOSE_ITERATION`) at the end of every coupling step after which both
solvers' fields on disk are of their current iteration, and at the end of the run. Each write goes to
`<prefix>coupling.tec.part`, then is renamed: a run killed during the write keeps the hand-off it had. An ICE restart
reads it back before the first MOSE block:

```text
 [MI2] ICE restarts with what its last coupling step handed the gas, from OUTPUT/part-coupling.tec (ICE iteration 100, MOSE iteration 200)
```

A run stopped after such a coupling step (killed, out of wall time, an `error stop`) is so continued exactly: the
restart equals the run in one go, the residual histories included (each solver drops the rows the stopped run wrote
past the iteration it restarts from). Hand-offs are written only at step ends where both solvers write their fields:
make each solver's `sol-diter` divide its `X-iter`, or a stopped run may leave none. A new ICE run (`newrun = true`)
removes the hand-off an earlier run left, once every set-up check has passed.

With `sol-overwrite = false` in `[MOSE-IO]` or `[ICE-IO]`, a restart reads the unnumbered `<prefix>field.tec` that
the end of the last completed run wrote, not the numbered outputs that follow it. The hand-off of that end is the one
that matches, so the numbered outputs in between write none. A run that never ended (its first run stopped) has
numbered fields and no hand-off: the rule below applies to it.

#### A restart that cannot be continued exactly

Two things make a restart other than the run in one go:

- **IGLOO iterates** (`IGLOO-iter` ≠ 0 in some cycle, and `newrun = false` in MOSE or ICE). IGLOO keeps neither its
  parcels nor the sources it hands the gas from one run to the next, so it starts empty on the restored fields.
- **ICE's hand-off cannot be used.** This happens when the file
  - is missing (a run stopped before its first step end with both outputs, or fields of an older binary);
  - is of another coupling step than the fields (the run stopped after a solver wrote its field in the middle of a
    step: ICE's field of another iteration than `CYCLE`, or MOSE's of another than `MI2_MOSE_ITERATION`);
  - is not a file MI2 wrote in full;
  - does not match ICE's set-up (blocks, dimensions, mesh coordinates, variables);
  - or a solver restarts without the counters of the run its field comes from (a time-accurate run from a steady
    field, the other way round, or a field without counters).

  The first MOSE block then runs without ICE's source.

Such a restart stops at set-up, before any iteration, when it **continues a time-accurate run**. That is, a cycle in
which IGLOO (or ICE) iterates is time-accurate and the solver that keeps that time resumes its counters:

- the gas the sources go to steps in physical time: `[MOSE-Numerics] time-accurate = true`, MOSE iterates in the
  schedule (sources stay in force until the next solve, whichever cycle MOSE iterates in) and MOSE resumes its
  counters; or
- ICE iterates in that cycle with `[ICE-Numerics] time-accurate = true`, and ICE resumes its counters.

`[HYDRA-Parameters] time-accurate` is not read for this (it only labels the end-of-run balances), nor is IGLOO's
coupling mode (`IGLOO-iter` 1 or 2): IGLOO's own solve has no time step. A time-accurate run started from a steady
field resumes no counters: it is a new run from that field, not a continuation.

- **Stop**, in a time-accurate continuation. The message names the cycle, the reason and the way out: a new run from
  the restored fields (`test/evap-box-ta` restarted):

    ```text
     [ERROR] MI2: newrun = false with IGLOO, and IGLOO iterates in cycle 1, a time-accurate cycle: the gas IGLOO hands its sources to steps in physical time ([MOSE-Numerics] time-accurate = true) and continues the restored run (MOSE resumes its counters). ... To go on from these fields, start a new run from them: keep a copy of INPUT/gas-ic.tec (the initial condition it replaces); copy OUTPUT/gas-field.tec to INPUT/gas-ic.tec; set newrun = true in [MOSE-Parameters]
    ```

  The way out is in the format each field was written in. A `.vtm` is copied with the block files it lists:
  `copy OUTPUT/gas-field.vtm to INPUT/gas-ic.vtm and OUTPUT/vtk/gas-fieldBlock*.vts to INPUT/vtk/`. When a solver's
  `ic-format` differs from its `sol-format`, the way out also says `set ic-format = <sol-format> in [MOSE-IO]` (or
  `[ICE-IO]`).
- **Stop**, too, when ICE's hand-off cannot be used and ICE iterates in no cycle (`ICE-iter = 0` everywhere). The gas
  receives ICE's source only once ICE has solved or the hand-off was read, so the run would end without it.
- **One warning** otherwise (a steady schedule, or counters not resumed); the run goes on:

    ```text
     [WARNING] MI2: ICE's hand-off is not used: OUTPUT/part-coupling.tec is missing (...). The first MOSE block runs without ICE's source. No time-accurate run is continued here (a steady schedule, or counters not resumed): a steady run converges again, a time-accurate one is a new run from the restored fields
     [WARNING] MI2: newrun = false with IGLOO: IGLOO keeps no parcels and no sources from one run to the next, so it starts empty on the restored fields: a new IGLOO start. ...
    ```

  A steady restart without the hand-off converges to the state of the run in one go (measured on
  `steady-ice-smoke` and `steady-ice-burnout`: 20 coupling steps restarted for 60 against 80 in one go, the
  fields equal to 1e-13 and 7e-11 of their range).
- **Nothing** is said when nothing is lost: no cycle runs IGLOO and the hand-off is read.

Not checked: the order of the materials (two cases with the same families per material in another order are taken for
one) and, for a time-accurate ICE, the time of the hand-off.

---

## `[MI2-Coupling]`

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `report` | `false` | Print the `[MI2-report]` lines at every coupling step (see [Coupling report](#coupling-report)) |
| `datum-correction` | `true` | Apply the enthalpy-datum offsets (see [Energy datum](#energy-datum)). `false` computes and prints them without applying them: diagnostic use only |
| `products-<material>` | — | No longer read: a non-empty value stops the run. A material's gas species come from `INPUT/<name>-products.txt` |

```ini
[MI2-Coupling]
report = true
```

---

## What is exchanged

### Gas to the condensed solvers

At each coupling step in which it runs, each condensed solver receives the gas state of every interior cell:

- **IGLOO**: the density of each gas species, the velocity, temperature, ratio of specific heats, gas constant,
  viscosity and conductivity; and, on the gas inlet faces, the gas mass flow \(\dot m\) [kg/s] of each face cell,
  which scales a particle inlet that injects in proportion to the gas (particle boundary 401). The mass flow is the
  one MOSE computes at the face for gas boundaries 401, 402, 405, 406 and 407, and the prescribed mass flux times the
  face area for 403, 404 and 420 (for a solid-propellant grain in its inflow mode, the burn law's flux at the cell
  pressure). A particle inlet 401 on a gas boundary of another kind stops the run. When IGLOO solves before MOSE's
  first iteration of the run (`MOSE-iter = 0` in the first cycle, or a restart that runs IGLOO first), the coupler
  first has MOSE evaluate the finest grid as an iteration does before its update (fluxes and boundary records, no
  update, no condensed source): the mass flow is then that of the initial field.
- **ICE**: density, velocity, temperature, gas constant, ratio of specific heats, viscosity and conductivity.

When IGLOO runs, the gas boundary records of `INPUT/<gas>-bc.txt` must be grouped by block, then by face 1 to 6.

### Condensed solvers to the gas

The sources are applied to MOSE's residual at every Runge-Kutta stage:

| Exchange | IGLOO, `IGLOO-iter = 1` (source) | IGLOO, `IGLOO-iter = 2` (Eulerian) | ICE |
|----------|----------------------------------|------------------------------------|-----|
| Momentum, energy | IGLOO's source field | from IGLOO's Eulerian fields (velocity, temperature, density, number density of each family) at the current gas state | ICE's momentum and energy sources of every family |
| Mass, species | each converting material's mass through its row of `products.txt` | none (refused with a converting material) | each converting material's mass through its row of `products.txt` |
| Enthalpy datum | offset on the converted mass | — | offset on the converted mass |
| Smoke | withheld from the gas, born in an ICE family | — | withheld from the gas (ICE keeps it in its smoke family) |
| Burn-out | — | — | what ICE removes at burn-out, through the material's row ([Burn-out in ICE](#burn-out-in-ice)) |

A material **converts** when it evaporates or burns. A material that does neither hands the gas momentum and energy
only. A row of `products.txt` gives the kg of each gas species exchanged per kg converted (negative: consumed).
When IGLOO runs in the last cycle that runs it, IGLOO computes both its source field and its Eulerian fields, so its
output carries both; the gas receives the mode of that cycle's `IGLOO-iter`.

**Consumed species.** A species a reaction consumes cannot leave a cell faster than the rate that empties the cell
over the stage's step (for ICE, that rate is shared between its families that consume the species). Each clamp is
counted; the count is printed at the end of the run, since a clamped reaction is mass-limited and no longer the
stoichiometry of the record.

**Non-finite sources.** After each solve, a non-finite value in what IGLOO or ICE hands the gas stops the run,
naming the solver, block, cell and slot.

---

## Evaporation and combustion

### `<name>-products.txt`

ATLAS GPB writes `INPUT/<name>-products.txt` for each condensed phase from its `[GPB-Phase*]` section (the vapour of
an evaporating material, the reactions of `particle-mechanism` for a burning one); the record format is in
[ATLAS GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/output/#prefixproductstxt). `hydra-MI2` reads it after each condensed
solver's setup and keeps, for each material, the record of the process the solver runs on it: `combustion` when it
burns, else `evaporation` when it evaporates. It stops when:

- a material evaporates or burns and the file is missing, or has no record of that process for it;
- the file was written for another gas phase (its `gas` line names another prefix or species count);
- a record names a material the solver does not have, or a gas species the gas phase does not have; a material of
  the solver is missing from the file;
- a record is malformed, or a row (with its smoke) does not sum to 1 within 1e-8.

Most of these messages end with `run ATLAS GPB again` or name the `[GPB-Phase*]` key to add. At setup the run prints
the rows it uses:

```text
 [MI2] IGLOO material -> gas species [kg of gas per kg of material converted], INPUT/part-products.txt
```

When no material of a solver converts, it prints `no material evaporates or burns: no mass is handed to the gas.`

A material name present in both ICE and IGLOO that burns on different models in the two, where a phase file does not
name the model on the material line (the solver's section default decided), gets a warning: write `combustion` in
both `[GPB-Phase*]` sections (`none` for a material that does not burn).

### Oxide cap (IGLOO)

An IGLOO material that burns with an oxide cap needs, besides its `combustion` record, a `cap` record (per kg of
material its cap reaction takes: the gas species exchanged and the kg of oxide formed, summing to 1) and, when its cap
is handed to the gas at burn-out (`residue = off`), a `release` record of the oxide. The run stops when one is missing, when the cap
record names another product than IGLOO's, when it does not close the mass, or when its kg of oxide per kg differs
from the ratio of IGLOO's cap law by more than 1e-7 relative.

### Burn-out in ICE

When ICE removes the last drops of a burning family that leaves no residue, what they held (mass, momentum, energy)
is stored and handed to the gas through the material's row: as a rate in a steady run, spread over the next MOSE block
in a time-accurate one.

### Smoke

A burning material whose reaction makes a condensed product, the **smoke**, has it born in an ICE family: the smoke
material is named, with its kg per kg burnt and the ICE phase, at the end of the material's `combustion` record (or of
the `release` record of an IGLOO cap handed at burn-out). The smoke takes mass, momentum and energy from what the
burning material hands over; the gas receives the rest. What ATLAS GPB requires of a smoke when it writes the files is
in [ATLAS GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/condensed-solid/#the-smoke-of-a-burning-material).
For a smoke of either solver the run stops unless:

- the smoke material is a material of ICE's phase, the phase the record names;
- it has exactly one family, of closure `MK`, and neither evaporates nor burns;
- the burning material's row of `models.txt` gives `d-smoke` > 0, the diameter each smoke drop is born with.

For the smoke of an IGLOO material, also:

- ICE runs;
- `relaxation = false`;
- `IGLOO-iter = 1` (a burning material already refuses `2`).

An IGLOO smoke is born with the velocity and temperature of the parcels that burnt in the cell; where no parcel gave
them, at rest at 298.15 K, with a warning. It reaches ICE as a source of its family over ICE's own steps (rates) or,
when MOSE and ICE are both time-accurate, as the amounts each MOSE block made, added after the block. The setup prints
one line per smoke:

```text
 [MI2] IGLOO material A burns into the smoke AL2O3(a) (ICE family 1): ...
```

### Energy datum

MOSE's species enthalpies include their heats of formation; ICE and IGLOO hold a material's enthalpy on their own
datum (the property table's). For each converting material, `hydra-MI2` adds to the energy the gas receives a constant
offset per kg converted, evaluated at \(T_\text{ref} = 298.15\) K:

- evaporation: \(\Delta h = \sum_s \nu_s\, h_s(T_\text{ref}) - \left(h_\text{solver}(T_\text{ref}) + L_v(T_\text{ref})\right)\);
- combustion: \(\Delta h = h_{298} - h_\text{solver}(T_\text{ref})\), with \(h_{298}\) the material's formation-inclusive
  enthalpy of its record;

and the same for a cap oxide and for a smoke. The setup prints the table under
`[MI2] enthalpy-datum reconciliation [J per kg of material converted]`. Two warnings size what the offset cannot fix:
an evaporating material with a constant latent heat (give it as the `Lv` column of the property table), and a
Beckstead `q-comb` that differs by more than 5 % from the heat of the reaction from the material at its `T-ign`.

---

## Run output

### Build stamp

The first lines of every run name the sources the executable was built from, one per repository (hydra, MOSE, ICE,
IGLOO, ORION, OSLO, FLINT, FiNeR):

```text
# [MI2] built from hydra 86dabae
# [MI2] built from ICE fbecd78
```

Each line holds the short commit, preceded by the release tag for ORION, OSLO and FLINT when they sit on one,
followed by `-dirty` when tracked files were modified, or `unknown` without git. The stamp is refreshed at every build.

### Setup

Before the first cycle the run prints the coupling schedule as resolved (`HYDRA-MI2 coupling schedule`: cycles,
which solvers run, each list after clamping), the [routing table](#the-routing-table), the species rows and datum
offsets of each solver, the smoke lines, and the blocks and materials each solver built (`HYDRA-MI2 solvers`).
Each cycle starts with `HYDRA-MI2: CYCLE <c> / <n>`.

### End-of-run balances

After `HYDRA-MI2: FINAL POSTPROCESS`, the run prints the mass sections when a converting material handed mass to
the gas, and the energy sections when IGLOO's source coupling or a converting ICE material handed it energy:

- `[MI2] mass injected into the gas over the run [kg]`, per species;
- `[MI2] closed-domain mass balance (valid only if NO boundary flux)`: per species, injected against the gain of the
  gas, with the relative error;
- `[MI2] closed-domain ENERGY bookkeeping (valid only if NO boundary flux)`: the energy handed to the gas, of which the
  datum correction, against the gain of the gas;
- `[MI2] enthalpy-datum diagnostic`: the change of the mass-weighted gas temperature and the change the datum
  correction accounts for.

These balances close only for a closed domain; with `[HYDRA-Parameters] time-accurate = false` they are pseudo-time
integrals. Coupling steps in IGLOO's Eulerian mode are excluded, with a warning. A run ends with `GAME OVER`.

### Coupling report

With `[MI2-Coupling] report = true`, each coupling step prints `[MI2-report]` lines (SI units; totals over the
interior cells):

| Line | Values |
|------|--------|
| `cycle <c> step <s> t0 <t> t1 <t>` | the MOSE block's start and end time |
| `gas-momentum` | gas momentum before and after the block (3 + 3) |
| `gas-energy` | gas total energy before and after the block |
| `igloo-to-gas` | IGLOO's force (3), energy rate, mass rate on the gas (source coupling) |
| `ice-to-gas` | ICE's force (3), energy rate, and the range of its relaxation time |
| `ice-to-gas-mass` | ICE's mass rate on the gas (a converting material) |
| `ice-removed-stored` | what ICE's burn-out removals still hold for the gas: mass, momentum (3), energy |
| `ice-removed-rate` | at a block's start, steady run: the removals' rate the gas receives |
| `ice-removed-delivered` | after a block, time-accurate run: the removals delivered |
| `igloo-smoke-withheld`, `igloo-smoke-rate` | steady run: the smoke kept out of the gas over the block; the rates handed to ICE (mass, momentum (3), energy, number) |
| `igloo-smoke-delivered` | time-accurate run: the smoke added to ICE after the block |
| `igloo-inlet-mdot` | the gas inflow handed to IGLOO [kg/s] and the number of face records |
| `ice-time` | steady run: ICE's time at the start and end of its block |

---

## Messages

A stop prints the message and ends the run with a non-zero exit status. Most messages begin with `[ERROR]`; the
table quotes their characteristic words.

| When | Message (start) | Cause |
|------|-----------------|-------|
| input | `"<key>" is not an option of section [HYDRA-Parameters], or it holds no value` | a list key outside `[HYDRA-Parameters]` or empty |
| input | `ICE-iter has <n> values but MOSE-iter has <m>` (same for `IGLOO-iter`) | list lengths differ |
| input | `relaxation-coeff must have either 1 value or one per IGLOO-iter cycle` | |
| input | `hydra-MI2 does not support multigrid` | `[HYDRA-Multigrid] levels` > 1 |
| routing | `no INPUT/*phase.txt listed in filelist.txt`; `cannot open phase file`; `phase file ... is empty` | |
| routing | `cannot parse phase file name`; `phase name ... contains "-"`; `is an UNNAMED phase but the case has other phase files`; `two phase files share the name`; `is a substring of`; `does not fit ...'s prefix buffer` | [naming rules](#naming-rules) |
| routing | `is a SOLID BULK phase`; `unknown phase type on line 1`; `malformed modeling token`; `modeling = "..." ... expected "lagrangian" (IGLOO) or "eulerian" (ICE)` | line 1 of a phase file |
| routing | `two ideal-gas phases`; `no ideal-gas phase file found` | |
| routing | `[...] phase = <name> names no condensed phase file`; `... is enabled but [...] phase is not set, and <n> condensed phase file(s) could be meant` | phase keys |
| routing | `is tagged "modeling=..." by ATLAS but [...] phase hands it to ...` | token and solver disagree |
| routing | `ICE and IGLOO both resolve to` | one phase for both solvers |
| routing | `<solver> owns phase "<name>" but <file> does not exist` | a required file |
| products | `[MI2-Coupling] products-<material> is no longer read` | retired key |
| products | `... material <m> and INPUT/<name>-products.txt is missing`; `... has no combustion record` / `no evaporation record` | |
| products | `was written for another gas phase`; `names the material ...`; `names the gas species ...`; `material ... is missing from`; `holds ...` | the file does not match the case: run ATLAS GPB again |
| products | `the row of ... with its smoke does not sum to 1` | |
| products | `burns with an oxide cap and ... has no cap record`; `the cap record of ... makes ...`; `does not close the mass`; `forms ... kg of oxide per kg, and the cap law of IGLOO holds ...`; `hands its cap to the gas at burnout (residue = off) and ... has no release record`; `the release record of ... does not sum to 1` | [oxide cap](#oxide-cap-igloo) |
| setup | `an IGLOO material evaporates or burns and a cycle has IGLOO-iter = 2` | |
| setup | `burns into the smoke ... and its row of INPUT/<name>-models.txt gives no d-smoke`; `which is not a material of ICE's phase`; `is born in the phase ..., and ICE runs the phase ...`; `has <n> families in ICE`; `whose closure is ...`; `evaporates or burns in ICE: a smoke does neither` | [smoke](#smoke) |
| setup | `the smoke has no family: ICE is not run`; `relaxation = true blends two coupling steps` | smoke of an IGLOO material |
| setup | `cycle <c> advances MOSE and ICE together with time-accurate on` | [time modes](#steady-and-time-accurate-runs) |
| setup | `ICE restarts ([ICE-Parameters] newrun = false) and its hand-off cannot be used: ...` (a time-accurate continuation, or ICE in no cycle) | [restart](#a-restart-that-cannot-be-continued-exactly) |
| setup | `[MOSE-Multigrid] level<m>-iter > 0 (m > 1) starts MOSE on its coarse grids, and the run is time-accurate ...` | run steady, or set `level<m>-iter = 0` for every m > 1 |
| setup | `[MOSE-Multigrid] level<m>-iter > 0 (m > 1) starts MOSE on its coarse grids, and relaxation = true ...` | set `relaxation = false`, or `level<m>-iter = 0` for every m > 1 |
| setup | `[MOSE-Multigrid] level1-iter is not set: with levels > 1 it is the iteration limit of MOSE's finest grid ...` | set `level1-iter` (the value `iter-threshold` would have) |
| setup | `[ICE-Multigrid] levels = <n>: the coupler hands the gas to ICE's finest grid ...` | ICE's coarse levels are not coupled: use `levels = 1` |
| setup | `newrun = false with IGLOO, and IGLOO iterates in cycle <c>, a time-accurate cycle` | a cold IGLOO start does not continue a time-accurate run: start a new run from the fields ([restart](#a-restart-that-cannot-be-continued-exactly)) |
| step | `MOSExIGLOO: IGLOO gasblock <b> extends beyond MOSE's ghosted range` | the gas and particle meshes disagree |
| step | `record <i> of INPUT/<gas>-bc.txt is of block ...` | gas boundary records out of order |
| step | `particles injected in proportion to the gas mass flow (particle BC 401) on a gas boundary of kind ...` | |
| step | `ICE family <g> (material <m>) loses mass and INPUT/<name>-products.txt has no record for its material` | |
| step | `ICE handed a non-finite source`; `IGLOO handed a non-finite source`; `IGLOO handed a non-finite Eulerian field` | |
| re-read | `[ERROR] Required parameter not set: <key>` | see [Reading `input.ini` again](#reading-inputini-again) |

Warnings and notes do not stop the run:

| Message (start) | Meaning |
|-----------------|---------|
| `[MI2] DEPRECATED: [...] phase not set` | a phase taken without its key |
| `[MI2] note: INPUT/<name>-properties.dat absent` | ICE uses its constant properties |
| `[MI2] WARNING: no INPUT/<name>-ic.{tec,szplt,plt}` | ICE needs an initial condition unless it restarts |
| `[WARNING] the latent heat is a constant` | a constant `Lv` against the property tables |
| `[WARNING] q-comb = ... differs from it by more than 5 %` | Beckstead's retained heat taken from another heat |
| `[WARNING] datum-correction = false` | offsets printed, not applied |
| `[WARNING] MI2: material ... burns as ... in ICE and as ... in IGLOO` | a section default decided one model |
| `[WARNING] MI2: in <n> cells ... a smoke was born with no parcel's state` | smoke born at rest, 298.15 K |
| `[MI2] ... source FROZEN this cycle` | [frozen sources](#frozen-sources) |
| `[WARNING] MI2: an overset mesh: <n> chimera (102) boundary cells on <b> block(s)` | set-up: the coupler makes no correction for a chimera overlap; every sum over the domain it prints (the gas content, what ICE and IGLOO hand the gas, the [coupling report](#coupling-report), the [end-of-run balances](#end-of-run-balances)) counts the overlap once per block that covers it |
| `[WARNING] MI2: input.ini, read again, advances <solver> after the iteration the schedule had made its last` | the rows and the solution due at that iteration were not written ([one coupling step](#one-coupling-step)) |
| `[WARNING] MI2: newrun = false with IGLOO` | IGLOO starts empty on the restored fields ([restart](#a-restart-that-cannot-be-continued-exactly)) |
| `[WARNING] MI2: ICE's hand-off is not used: ...` | the hand-off cannot be used and no time-accurate run is continued: the first MOSE block runs without ICE's source ([restart](#a-restart-that-cannot-be-continued-exactly)) |
| `[MI2] ICE starts anew: OUTPUT/<prefix>coupling.tec of an earlier run removed` | [restart](#restart) |
| `[MI2] MOSE (ICE) made no iteration in this run` | a restarted solver writes no end-of-run output |
| `[MI2] WARNING: the consumed-species floor guard clamped <n> (cell,species) extractions` | end of run: a mass-limited reaction |
| `[WARNING] <n> coupling step(s) used the EULERIAN branch` | end of run: those steps are not in the balances |
| `[WARNING] [ICE::Setup_Materials] ... is ignored: ICE has no breakup model` | a `breakup=` word or a `breakup` row in the phase ICE solves: ICE has no breakup model and goes on without it |

---

## Test cases

The cases of `test/` that run `hydra-MI2` and the gate scripts of `test/mi2-checks/` are listed in
[Testing](../development/testing.md#hydra-mi2).
