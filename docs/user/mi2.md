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
    └── partL-models.txt      # the constants of the burning models (ATLAS GPB), read by IGLOO
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

`models.txt` exists only when a material has constants of a burning model; the solvers read it, `hydra-MI2` does not.
Its format and the keys it carries are in the solvers' input pages
([ICE](https://open-hydra.github.io/ICE/user/input/), [IGLOO](https://open-hydra.github.io/IGLOO/user/input/)) and
in [ATLAS GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/condensed-solid/#the-constants-of-the-burning-models-namemodelstxt).

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
the current block of its iterations, never the cycle. After the last cycle every solver writes its final solution.

### Frozen sources

A solver that has produced sources keeps them in force when a later cycle does not advance it (`ICE-iter = 0` or
`IGLOO-iter = 0`): the gas goes on receiving them. The run says so at the first coupling step of such a cycle, then
each time the cumulative MOSE iteration count crosses a multiple of `[MOSE-IO] shell-diter`:

```text
 [MI2] IGLOO source FROZEN this cycle (IGLOO-iter = 0): the gas keeps the sources produced in cycle 1, coupling step 3 [cycle 2, step 1, MOSE iter 150]
```

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

### Other refusals

- Multigrid: `[HYDRA-Multigrid] levels` above 1 stops the run (IGLOO has no grid levels).
- `IGLOO-iter = 2` in any cycle while an IGLOO material evaporates or burns: the Eulerian coupling hands the gas
  momentum and energy only. Use `IGLOO-iter = 1` in every cycle that runs IGLOO.

### Reading `input.ini` again

Every `[HYDRA-IO] ini-diter` coupling steps (default 1000) the values of the `[HYDRA-*]` keys are read again, so a
long run can be steered without a restart; the number of values of a list is not. The re-read checks the registry
(required keys and the allowed range of each scalar key): `N-coupling` and `time-accurate` must then be present (`[ERROR] Required parameter not set: <key>`); `FUSS-iter` and
the multigrid levels are not required. `[MI2-Coupling]` and the phase keys are read once, at start.

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
  pressure). A particle inlet 401 on a gas boundary of another kind stops the run.
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
| `[MI2] WARNING: the consumed-species floor guard clamped <n> (cell,species) extractions` | end of run: a mass-limited reaction |
| `[WARNING] <n> coupling step(s) used the EULERIAN branch` | end of run: those steps are not in the balances |

---

## Test cases

The cases of `test/` that run `hydra-MI2` and the gate scripts of `test/mi2-checks/` are listed in
[Testing](../development/testing.md#hydra-mi2).
