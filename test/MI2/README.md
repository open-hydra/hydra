# hydra-MI2 tests

`cases/` holds one small case per path along which mass, momentum or energy moves among MOSE (gas), IGLOO
(Lagrangian) and ICE (Eulerian). Each case is judged by its `check.py` (or `check.sh`), which lists its claims in its
header and prints one PASS/FAIL line per claim. `validation/` holds the cases compared against references.

```bash
test/MI2/run.sh [-x exe] [-s scratch] [-t threads] [case ...]   # all of cases/ but fpe
test/MI2/run.sh fpe                                              # builds an -fpe0 binary first (FPE_EXE: use one)
```

The executable defaults to `bin/hydra_MI2`, the scratch directory to `/tmp/mi2-tests-$USER`; the scratch directory must
lie outside `test/`. Each case's PASS/FAIL lines are left in `<scratch>/<case>.check`.

## cases/

| Case | What moves | Judged |
|---|---|---|
| `igloo-evap` | IGLOO drops evaporate into the gas (evaporation record) | mass, species, energy books |
| `igloo-evap-sa` | the same, the gas RANS (Spalart-Allmaras) | the books; the turbulence variable left alone |
| `igloo-two-materials` | two IGLOO materials, each into its own species (evaporation, combustion) | each material's mass to its species |
| `igloo-burn` | IGLOO drops burn into the gas (combustion record) | mass, momentum, heat of reaction, energy against closed forms |
| `igloo-al` | IGLOO aluminium burns, carries its oxide, generates alumina in ICE | mass, species, momentum, energy across the three solvers; number; hand-off |
| `igloo-al-release` | as `igloo-al`, the carried oxide released into ICE at burnout (release record) | the same |
| `ice-exchange` | ICE drag and heat to the gas, no mass | the frozen source against closed forms |
| `ice-burn` | ICE drops burn into the gas | mass, heat of reaction, energy against closed forms |
| `ice-burnout` | ICE drops burn out; what is left at d-burnout is stored and handed whole to the gas | the books of the removal |
| `ice-burnout-steady` | `ice-burnout` in pseudo-time on a nozzle | the converged field independent of ICE's step |
| `ice-al` | ICE aluminium burns, generates alumina in another ICE family | one ICE step against closed forms; birth state; the gas's species |
| `ice-igloo-burn` | ICE and IGLOO burning in one gas | the books per phase |
| `ch-gas-ice` | `[MI2-Channels]`: gas → ICE, born as new drops | ICE's gain; the gas's rate |
| `ch-gas-igloo` | `[MI2-Channels]`: gas → IGLOO, absorbed (lagged books) | the lag; steady absorption; the gas's books against IGLOO's report |
| `ch-steady` | `[MI2-Channels]`: IGLOO → ICE, ICE → IGLOO, gas → ICE in pseudo-time | the rates of each channel |
| `refusals` | `[MI2-Channels]` lines refused at setup | the message, before any step |

The cases below test what is not a path of mass: the solvers against each other, the set-up, and the run.

| Case | What | Judged |
|---|---|---|
| `twins` | one burning section in ICE and in IGLOO (`twins/igloo`): direct, Beckstead, burnout | each its closed forms; the energy per kg and the law constant, ICE against IGLOO |
| `burnout` | `ice-burnout` with a long step, frozen (one and two ICE steps), cold; `d-stop` beside `d-burnout` | `ice-burnout/check.py` in each mode; the refusal |
| `routing` | which solver takes which `INPUT/<phase>phase.txt` | exit status and message of 20 decks |
| `ini-reread` | `[HYDRA-IO] ini-diter`: input.ini re-read during the run | OUTPUT byte-identical to the run without |
| `floor-guard` | the consumed-species guard on rates | no clamp of a legitimate sink; one warning on an empty species |
| `multigrid` | MOSE's start on its coarse grids and its V cycle, with IGLOO, ICE and the law channels | what each hands the gas on the coarse level against the finest grid; the refusals |
| `chimera` | an overset mesh | one warning with MOSE's chimera-cell count |
| `ice-bc` | ICE boundary types 301 and 403; the population count | byte identity with 300; the prescribed inflow |
| `mass-flux-inlet` | the gas inflow handed to IGLOO at a 401 or 403 inlet, also before the first gas iteration | the report against g x area; IGLOO's injected flux |
| `restart` | MOSE + ICE restarts (`newrun = false`), killed runs, IGLOO restarts | every value equal to the run in one go; the stops and their way out |
| `fpe` | `igloo-evap` at gas-order 2 under `-fpe0` | no trap |

A variant case holds its own `input.ini` and links to its base's `INPUT/` and `MESH/`. A `check.sh` builds its decks
from other cases with the helpers of `checklib.sh`. A channel case is judged
against its base's run (with the channel minus without), so the base runs first. The checks share the readers of
`mi2check.py`. ATLAS writes the inputs: `regen_case.sh` regenerates them; `run_case.sh` runs one case in a scratch
directory.

## validation/

| Case | Against | Run |
|---|---|---|
| `JPL-Lagrangian-20micron`, `JPL-Eulerian-20micron`, `JPL-3solver-20micron` | the JPL nozzle, Chang et al. (`JPL_chang.dat`) | `check_jpl.sh` |
| `JPL-3solver-hot-2micron`, `JPL-3solver-hot-20micron` | alumina freezing in the hot nozzle, ICE against IGLOO | `JPL-3solver-hot-2micron/check.sh` |
| `evap-box` | an independent integration of the drop equations (`REFERENCE.md`) | `evap-box/check.sh` |
| `vortex-cloud` | ICE and IGLOO standalone against the exact solution | `vortex-cloud/compare.py` |
