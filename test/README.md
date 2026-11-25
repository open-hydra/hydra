# hydra tests

| Folder | Application | What it holds |
|---|---|---|
| `MF/flatplate-Tw/` | hydra-MF (MOSE + FUSS) | conjugate flat plate, three solid conductivities, against the Lindstedt correlation |
| `AF/flatplate-Tw/` | hydra-AF (ARES + FUSS) | the same plate with a real fluid |
| `MI2/` | hydra-MI2 (MOSE + IGLOO + ICE) | `cases/`: the coupling tests (`MI2/run.sh`); `validation/`: nozzles and references |

`make-vv-figures.py` writes the figures of `docs/vv/` from the flat plates' `reference/` solutions.

## ctest

A build registers the tests of its application (`test/CMakeLists.txt`); `test/tiers.txt` gives each its tier by its
time at one thread:

```bash
ctest --test-dir build/MI2 -L light               # under 5 s each
ctest --test-dir build/MI2 -L "light|medium" -j 8 # under 100 s each: about 1 min
ctest --test-dir build/MI2 -L heavy               # restart, mass-flux-inlet, fpe and the validation
```

Other labels: the kind of check (`exact`, `identity`, `contract`, `robustness`, `reference`), the application
(`MI2`) and `validation`. The runs go to `<build>/test/<application>/`.
