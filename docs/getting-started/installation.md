# Installation

This page describes how to build Hydra in this repository.

Hydra uses a unified command entrypoint (`figaro`) to build both
the solvers and Hydra executables.

## Prerequisites

- Git
- CMake version 3.23 or higher
- Fortran compiler (GNU or Intel toolchain)
- Python 3
- conda
- OpenMP and MPI for parallel execution (optional)
- C and C++ compilers for building solver dependencies (optional)

## Clone

```bash
git clone --recurse-submodules https://github.com/open-hydra/hydra.git
cd hydra
```

## Environment bootstrap

Run the project bootstrap script once in your shell:

```bash
source ./bootstrap.sh
```

This runs the sanity checks, requirements installation, and sets up the environment for the current shell session. If it completes successfully, you can proceed to the next step.

## Build a Hydra executable

Example for the coupled executable `hydra-MF`:

```bash
figaro build hydra-MF --compilers=gnu --use-openmp
```

For the multiphase executable `hydra-MI2` (MOSE + ICE + IGLOO), OpenMP only:

```bash
figaro build hydra-MI2 --compilers=intel --use-openmp
```

ICE is a private repository, cloned over SSH: this build needs access to it.

!!! note
    Build options are forwarded to the underlying CMake/install scripts.
    Use `figaro` with no arguments to print the command help.

## Components and their releases

The solvers (MOSE, FUSS, ARES, ICE, IGLOO, GROOT), ATLAS and the libraries ORION, OSLO and FLINT are not part of the
clone: each lives in its own repository, and hydra requires a range of its releases (`cmake/Releases.cmake`). A range
`1.1` accepts 1.1.0 and every later 1.x, not 2.0.0; a range `0.2` accepts 0.2.x only. The build checks each component
it needs out in its folder (`src/lib/<NAME>`, `utils/ATLAS`, `lib/<NAME>`) at the release of `deps.lock`, or, when
`deps.lock` has none in the range, at the newest release on its remote, and writes it to `deps.lock`. FiNeR is the
one submodule.

```bash
figaro deps                     # every component but ICE at its release of deps.lock
figaro deps --update MOSE       # MOSE at the newest release in its range, written to deps.lock
figaro build hydra-MF --update-deps   # the same for every component of the build
```

A component folder is used as it is when it is at the release; at another release it is moved, unless it is on a
branch or has changes: the configure then stops. To build a checkout of your own (a branch under development), give
it to the configure: `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<folder>`, e.g. `-DFETCHCONTENT_SOURCE_DIR_MOSE=src/lib/MOSE`.

## Build an individual solver (optional)

You can also build solver repositories directly from Hydra:

```bash
figaro build MOSE --use-openmp
figaro build FUSS --use-mpi
```

## Next step

Proceed to [Quick Start](quick-start.md) to run a first Hydra case.
