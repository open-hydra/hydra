# Figaro CLI Reference

> "Pronto prontissimo son come il fulmine,  
> son il factotum della citta'."
>
> - Gioachino Rossini, *Il Barbiere di Siviglia*

Hydra employs `figaro` as the command-line entrypoint for every Hydra workflows.

It wraps common tasks for:

- environment checks
- submodule management
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
| `submodule ...` | Update or deinit git submodules |
| `build ...` | Configure and build hydra target or solver submodule |
| `compile ...` | Recompile using existing presets/build dirs |
| `ATLAS ...` | Run ATLAS preprocessing helpers and import outputs |
| `<solver> ... solve` | Run solver executable (`MOSE`, `FUSS`, `GROOT`, `PESTE`, `hydra-*`) |
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

### `submodule`

```bash
figaro submodule update [--remote] [<name>]
figaro submodule delete [--force] [<name>]
```

Examples:

```bash
figaro submodule update
figaro submodule update --remote MOSE
figaro submodule delete --force FUSS
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

Example:

```bash
figaro build hydra-MF --compilers=gnu --use-openmp
```

**Build solver submodule directly**

```bash
figaro build MOSE [solver install.sh options]
figaro build FUSS [solver install.sh options]
figaro build GROOT [solver install.sh options]
```

The command delegates to each solver's `install.sh build` and injects
paths to shared hydra dependencies.

### `compile`

```bash
figaro compile [--edit] <target>
```

Targets:

- hydra executable target: `hydra-MF`
- solver submodule target: `MOSE`, `FUSS`, `GROOT`, `ATLAS`

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
- `hydra-MF` (or another `hydra-*` executable name)

Run options:

- `-p, --parallel <n>` set `OMP_NUM_THREADS`
- `-b, --background` run in background and store PID in `.ID`

Examples:

```bash
figaro hydra-MF solve
figaro hydra-MF -p 8 solve
figaro MOSE -b solve
```

## `ATLAS` helper mode

```bash
figaro ATLAS <step1> <step2> ...
figaro ATLAS import
```

This mode is a convenience wrapper around the ATLAS scripts and supports
importing generated files from `fromATLAStoSolver/` into `INPUT/`.
