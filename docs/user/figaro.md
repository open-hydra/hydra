# Figaro CLI Reference

> "Pronto prontissimo son come il fulmine,  
> son il factotum della citta'."
>
> - Gioachino Rossini, *Il Barbiere di Siviglia*

Hydra employs `figaro` as the command-line entrypoint for every Hydra workflows.

It wraps common tasks for:

- environment checks
- the components at their releases
- build and compile
- simulation setup, launch, and cleanup

!!! note "Setup"
    `figaro` is defined as a shell function by hydra bootstrap scripts. After cloning the repository, run `source ./bootstrap.sh` to set up the environment and enable `figaro` in your shell session.

!!! tip "Autocompletion"
    `figaro` is provided with autocompletion support for bash and zsh. Double-tap the tab key after `figaro` to see available commands and options.

## General Syntax

```bash
figaro <command> [options] [args]
```

## Command Overview

| Command | Purpose |
|---------|---------|
| `doctor` | Run environment and toolchain checks |
| `deps ...` | Check components out at their releases (`deps.lock`) |
| `build ...` | Configure and build a hydra target or a solver |
| `compile ...` | Recompile using existing presets/build dirs |
| `ATLAS ...` | Run ATLAS preprocessing helpers and import outputs |
| `<solver> ... solve` | Run solver executable (`MOSE`, `FUSS`, `GROOT`, `hydra-*`) |
| `clean` | Remove local simulation outputs (`OUTPUT`, local `bin`, logs) |
| `kill` | Stop background run saved in `.ID` |

### `doctor`

```bash
figaro doctor
```

Checks include:

- `HYDRADIR` availability
- required tools (`cmake`, `git`, `bash`)
- compiler presence
- optional MPI/OpenMP checks
- CMake presets validity

### `deps`

```bash
figaro deps [--update] [<name>...]
```

Checks each component out in its folder at its release of `deps.lock`, or the newest in the range of
`cmake/Releases.cmake` when `deps.lock` has none; `--update` selects the newest releases in the ranges again, written
to `deps.lock`. Without a name: every component but ICE (private, cloned over SSH: name it). A build checks out what it
needs by itself (see [Installation](../getting-started/installation.md#components-and-their-releases)).

Examples:

```bash
figaro deps
figaro deps --update MOSE IGLOO
figaro deps ICE
```

### `build`

**Build hydra executable target**

```bash
figaro build hydra-MF [options]
```

Hydra build options:

- `--compilers=gnu|intel`
- `--use-openmp`
- `--use-mpi`
- `--use-tecio`
- `--use-sundials`
- `--use-cantera`
- `--type=<BuildType>`
- `--update-deps`: select again the newest releases of the build's components (the libraries and the solvers) in the
  ranges of `cmake/Releases.cmake`, written to `deps.lock` (commit it to take them)

Example:

```bash
figaro build hydra-MF --compilers=gnu --use-openmp
figaro build hydra-MI2 --compilers=intel --use-openmp
```

`hydra-MI2` is OpenMP only: it refuses to configure with `--use-mpi` (see [Multiphase Coupling](mi2.md#build)).

**Build a solver directly**

```bash
figaro build MOSE [solver install.sh options]
figaro build FUSS [solver install.sh options]
figaro build GROOT [solver install.sh options]
```

The command checks the solver out at its release (`figaro deps`), then delegates to its `install.sh build` and
injects paths to shared hydra dependencies.

### `compile`

```bash
figaro compile [--edit] <target>
```

Targets:

- hydra executable target: `hydra-MF`, `hydra-AF`, `hydra-MI2`
- solver target: `MOSE`, `FUSS`, `GROOT`, `ATLAS`

Examples:

```bash
figaro compile hydra-MF
figaro compile --edit hydra-MF
figaro compile MOSE
```

`--edit` opens the relevant `CMakePresets.json` before compilation.

## Simulation Commands

### Run executable

```bash
figaro <target> [run-options] solve
```

Supported runtime targets:

- `MOSE`
- `FUSS`
- `GROOT`
- `hydra-MF`, `hydra-AF`, `hydra-MI2`

Run options:

- `-p, --parallel <n>` set `OMP_NUM_THREADS`
- `-b, --background` run in background and store PID in `.ID`

Examples:

```bash
figaro hydra-MF solve
figaro hydra-MF -p 8 solve
figaro MOSE -b solve
figaro hydra-MI2 -p 8 solve
```

A `hydra-*` run first lists `INPUT/*phase.txt` into `filelist.txt`, which the executable reads.

## `ATLAS` helper mode

```bash
figaro ATLAS <step1> <step2> ...
figaro ATLAS import
```

This mode is a convenience wrapper around the ATLAS scripts and supports
importing generated files from `fromATLAStoSolver/` into `INPUT/`.
