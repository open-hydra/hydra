# Installation

This page describes how to build Hydra in this repository.

Hydra uses a unified command entrypoint (`figaro`) to build both
solver submodules and Hydra executables.

## Prerequisites

- Git with submodule support
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

It also checks out the ICE and IGLOO submodules; ICE is skipped by a recursive clone and fetched over SSH, so the
build needs access to its repository.

!!! note
    Build options are forwarded to the underlying CMake/install scripts.
    Use `figaro` with no arguments to print the command help.

## Build an individual solver (optional)

You can also build solver repositories directly from Hydra:

```bash
figaro build MOSE --use-openmp
figaro build FUSS --use-mpi
```

## Next step

Proceed to [Quick Start](quick-start.md) to run a first Hydra case.
