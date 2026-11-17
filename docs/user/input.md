# Input File

Hydra is configured through a single **INI-format** file named `input.ini`, located in the run directory.

In coupled runs (for example `hydra-MF`), this file contains both:

- **Hydra-level sections** (`HYDRA-*`) controlling coupling orchestration
- **Solver-level sections** (such as `MOSE-*`, `ARES-*`, `FUSS-*`, `ICE-*` and `IGLOO-*`) consumed by each solver

Hydra-level parameters define how the solvers are coordinated, while
solver-level parameters define each solver's own physics and numerics.

## File Structure

```ini
[SECTION-NAME]
parameter = value
```

Example:

```ini
[HYDRA-Parameters]
N-coupling = 5000
FUSS-iter = 1
MOSE-iter = 2
time-accurate = false

[MOSE-Numerics]
cfl = 0.9
time-scheme = RK3

[FUSS-Numerics]
vnn = 0.8
space-reconstruction = second-order
```

!!! warning "Default values"
    If a parameter is not specified, the **default value** is used silently.

!!! warning "Parameter names are case-sensitive"
    Use the exact parameter names expected by Hydra and each solver.

---

## Hydra-Level Sections

Hydra-level input as well as solver-level input may be generated automatically inside each solver repository.

Hydra-level sections are reported in the [Hydra Input Reference](registry.md) and include:

### `[HYDRA-Parameters]`

Core coupling-loop controls.

| Parameter | Default | Meaning |
|----------|---------|---------|
| `N-coupling` | `1000000000` | Number of coupling iterations (array supported) |
| `FUSS-iter` | `1` | Number of FUSS iterations per coupling step (array supported) |
| `MOSE-iter` | `1` | Number of MOSE iterations per coupling step (array supported) |
| `ARES-iter` | `1` | Number of ARES iterations per coupling step (array supported, `hydra-AF`) |
| `time-accurate` | `false` | Enables synchronized time-accurate mode |
| `time-threshold` | `1e30` | Maximum simulation time |
| `iter-threshold` | `1000000000` | Maximum total iteration threshold |

`hydra-MI2` reads four more keys in this section; there a list holds one value per cycle, and the number of cycles
is the number of values of `MOSE-iter` (see [Multiphase Coupling](mi2.md#coupling-schedule)):

| Parameter | Default | Meaning |
|----------|---------|---------|
| `ICE-iter` | `0` | Number of ICE iterations per coupling step (array, one value per cycle; absent: ICE is not run, `hydra-MI2`) |
| `IGLOO-iter` | `0` | IGLOO mode per cycle: `0` off, `1` source coupling, `2` Eulerian coupling (array; absent: IGLOO is not run, `hydra-MI2`) |
| `relaxation` | `false` | Under-relax IGLOO's source on the gas (`hydra-MI2`) |
| `relaxation-coeff` | `0.5` | Under-relaxation coefficient: one value, or one per cycle (`hydra-MI2`) |

In `hydra-MI2`, `time-accurate` only labels the end-of-run balances: MOSE and ICE take their time mode from
`[MOSE-Numerics]` and `[ICE-Numerics]`.

### `[HYDRA-IO]`

Coupling-output and runtime reload settings.

| Parameter | Default | Meaning |
|----------|---------|---------|
| `coupling-format` | `tecplot ascii` | Coupling-output file format |
| `coupling-diter` | `1000000000` | Coupling output frequency by iteration |
| `coupling-dtime` | `1e30` | Coupling output frequency by physical time |
| `coupling-overwrite` | `true` | Overwrite or append coupling wall output |
| `coupling-wall` | `false` | Enables coupled wall exchange/output |
| `ini-diter` | `1000` | Frequency for re-reading `input.ini` at runtime |

### `[HYDRA-Multigrid]`

Per-level coupling thresholds in multigrid mode.

| Parameter pattern | Meaning |
|-------------------|---------|
| `level1-iter`, `level2-iter`, ... | Iteration thresholds for each multigrid level |

`hydra-MI2` runs on one level: `levels` above 1 stops it.

### `[MI2-Coupling]`

`hydra-MI2` only.

| Parameter | Default | Meaning |
|----------|---------|---------|
| `report` | `false` | Print the `[MI2-report]` lines (what the coupler hands the gas) at every coupling step |
| `datum-correction` | `true` | Apply the enthalpy-datum offsets of the converting materials; `false` is for diagnosis only |

---

## Solver-Level Sections

`input.ini` must also include sections required by each solver you run.

- For MOSE options and accepted values, see [MOSE User Guide](https://open-hydra.github.io/MOSE/user/)
- For ARES options and accepted values, see [ARES User Guide](https://open-hydra.github.io/ARES/user/)
- For FUSS options and accepted values, see [FUSS User Guide](https://open-hydra.github.io/FUSS/user/)
- For GROOT options and accepted values, see [GROOT User Guide](https://open-hydra.github.io/GROOT/user/)
- For ICE options and accepted values, see [ICE User Guide](https://open-hydra.github.io/ICE/user/)
- For IGLOO options and accepted values, see [IGLOO User Guide](https://open-hydra.github.io/IGLOO/user/)

In other words, Hydra orchestration keys and solver keys coexist in the same INI file.

The two fluid solvers take the flow model with the same keys, so a fluid block moves from `hydra-MF` (MOSE) to `hydra-AF` (ARES) by renaming the sections:

```ini
[MOSE-Physics]              ; [ARES-Physics] in hydra-AF
equations  = navier-stokes  ; euler, navier-stokes
turbulence = SST            ; none (default) = laminar

[MOSE-Turbulence]           ; [ARES-Turbulence] in hydra-AF
Prt = 0.9
```

