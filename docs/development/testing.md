# Testing

This page documents framework-level validation in the Hydra repository.

Because Hydra is an orchestration layer over solver libraries, testing
has two complementary scopes:

- framework checks (build, runtime orchestration, coupling behavior)
- solver checks (done in solver repositories)

## Current Test Assets in This Repository

```text
test/
├── hydra-MF.sh                 # legacy helper script
└── Tw-Flat-Plate/
    ├── analytical.py           # analytical comparison and plotting
    ├── ks_kf_40/
    ├── ks_kf_400/
    ├── ks_kf_4000/
    ├── SOLUTION/               # generated/collected coupling-wall solutions
    └── reference/              # reference PDFs and literature material
```

Each `ks_kf_*` case includes:

- `simulation.ini` (Hydra + solver runtime configuration)
- `atlas.ini` and `mesh.ini` (pre/post-processing setup)
- `INPUT/` and `MESH/` data used by the run

## Recommended Regression Workflow

### 1. Environment sanity

```bash
figaro doctor
```

This verifies compilers, CMake version, and expected repository layout.

### 2. Build and compile target executable

```bash
figaro build hydra-MF --compilers=gnu --use-openmp
figaro compile hydra-MF
```

### 3. Run representative benchmark case

Example case folder:

```bash
cd test/Tw-Flat-Plate/ks_kf_40
figaro hydra-MF solve
```

Repeat for `ks_kf_400` and `ks_kf_4000` when validating full coverage.

### 4. Post-run analytical check

From `test/Tw-Flat-Plate/`:

```bash
python3 analytical.py
```

This script compares numerical wall-temperature trends against the
analytical/Lindstedt correlation and generates diagnostic plots.

## What to Verify in a Framework Change

When modifying hydra-level code (`src/app/common`, `src/app/MF`, or
`scripts/`), check at least:

- build succeeds from a clean build directory
- executable starts and parses `HYDRA-*` sections correctly
- coupling loop progresses and writes expected outputs
- benchmark trend remains consistent after the change

## Framework Test Matrix

| Change type | Minimum checks |
|-------------|----------------|
| CLI/scripts (`scripts/*`) | `figaro doctor`, build, one benchmark run |
| Input/config (`src/app/common/config/*`) | parse/launch check + one benchmark run |
| Coupling internals (`src/app/MF/*`) | all three `ks_kf_*` runs + `analytical.py` |
| Build/CMake changes | clean build + run benchmark |

## Relationship to Solver-Level Tests

Hydra does not replace solver regression suites.

- Use solver repositories for detailed physics/numerics validation.
- Use Hydra tests to validate framework integration and coupling.

For solver-specific test campaigns, see the proper solver repositories.