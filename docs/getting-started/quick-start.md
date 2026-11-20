# Quick Start

This quick start runs a Hydra executable and explains the split between:

- Hydra-level orchestration input (`HYDRA-*` sections)
- Solver-level input (`MOSE-*`, `FUSS-*`, ...)

## 1. Build `hydra-MF`

From the repository root:

```bash
source ./bootstrap.sh
figaro build hydra-MF --compilers=gnu --use-openmp
figaro compile hydra-MF
```

## 2. Enter a test case

Use one of the available test folders, for example:

```bash
cd test/MF-flatplate-Tw/ks_kf_40
```

A typical case includes:

- `input.ini` with Hydra and solver blocks
- `INPUT/` folder with solver data files
- output written under `OUTPUT/` after run

## 3. Run

```bash
figaro hydra-MF solve
```

Optional threading:

```bash
figaro hydra-MF -p 4 solve
```

## 4. Read what is happening

At runtime, Hydra coordinates solver steps according to the Hydra-level
controls in `input.ini`:

- `[HYDRA-Parameters]` for coupling iteration control
- `[HYDRA-IO]` for coupling output and runtime reload behavior
- `[HYDRA-Multigrid]` for level-dependent coupling thresholds

Solver blocks remain solver-owned and follow each solver registry.
