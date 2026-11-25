# hydra tests

| Folder | Application | What it holds |
|---|---|---|
| `MF/flatplate-Tw/` | hydra-MF (MOSE + FUSS) | conjugate flat plate, three solid conductivities, against the Lindstedt correlation |
| `AF/flatplate-Tw/` | hydra-AF (ARES + FUSS) | the same plate with a real fluid |
| `MI2/` | hydra-MI2 (MOSE + IGLOO + ICE) | `cases/`: the coupling tests (`MI2/run.sh`); `validation/`: nozzles and references |

`make-vv-figures.py` writes the figures of `docs/vv/` from the flat plates' `reference/` solutions.
