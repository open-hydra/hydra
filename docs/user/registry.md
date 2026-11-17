# Input Parameters

## HYDRA-Parameter

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
| `ICE-iter` | `0` | Number of ICE iterations per coupling iteration (array supported; absent: ICE disabled, `hydra-MI2`) |
| `IGLOO-iter` | `0` | IGLOO mode per cycle: 0 = off, 1 = source coupling, 2 = Eulerian coupling (array supported; absent: IGLOO disabled, `hydra-MI2`) |
| `relaxation` | `false` | Under-relax the IGLOO → MOSE source term (`hydra-MI2`) |
| `relaxation-coeff` | `0.5` | Under-relaxation coefficient: one value, or one per coupling cycle (`hydra-MI2`) |

## HYDRA-IO

Coupling-output and runtime reload settings.

| Parameter | Default | Meaning |
|----------|---------|---------|
| `coupling-format` | `tecplot ascii` | Coupling-output file format |
| `coupling-diter` | `1000000000` | Coupling output frequency by iteration |
| `coupling-dtime` | `1e30` | Coupling output frequency by physical time |
| `coupling-overwrite` | `true` | Overwrite or append coupling wall output |
| `coupling-wall` | `false` | Enables coupled wall exchange/output |
| `ini-diter` | `1000` | Frequency for re-reading `input.ini` at runtime |

## HYDRA-Multigrid

Per-level coupling thresholds in multigrid mode.

| Parameter pattern | Meaning |
|-------------------|---------|
| `level1-iter`, `level2-iter`, ... | Iteration thresholds for each multigrid level |

## MI2-Coupling

Read by `hydra-MI2` only, outside the registry.

| Parameter | Default | Meaning |
|----------|---------|---------|
| `report` | `false` | Print the coupling totals (`[MI2-report]` lines) at every coupling step |
| `datum-correction` | `true` | Apply the enthalpy-datum offsets; `false` computes and prints them without applying them |
