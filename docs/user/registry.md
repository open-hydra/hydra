# Input Parameters

## HYDRA-Parameter

Core coupling-loop controls.

| Parameter | Default | Meaning |
|----------|---------|---------|
| `N-coupling` | `1000000000` | Number of coupling iterations (array supported) |
| `FUSS-iter` | `1` | Number of FUSS iterations per coupling step (array supported) |
| `MOSE-iter` | `1` | Number of MOSE iterations per coupling step (array supported) |
| `time-accurate` | `false` | Enables synchronized time-accurate mode |
| `time-threshold` | `1e30` | Maximum simulation time |
| `iter-threshold` | `1000000000` | Maximum total iteration threshold |

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
