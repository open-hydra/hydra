<p align="center">
  <h1 align="center">Hydra</h1>
  <p align="center"><b>A multi-physics CFD framework for coupled compressible flow, solid heat transfer, and reactive simulations</b></p>
</p>

<p align="center">
  <a href="https://open-hydra.github.io/hydra/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran%20%7C%20C%2B%2B-734f96.svg" alt="Language: Fortran | C++">
  <a href="https://github.com/open-hydra/hydra/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
  <img src="https://img.shields.io/badge/build-CMake-064F8C.svg" alt="Build: CMake">
  <img src="https://img.shields.io/badge/parallelism-OpenMP-E84393.svg" alt="OpenMP, MPI">
</p>

---

**Hydra** is an open-source multi-physics computational framework that composes specialized CFD solvers into coupled simulation pipelines. Instead of a monolithic code, Hydra orchestrates focused solver libraries — each targeting a specific physical regime — through a unified build system and a single command-line interface.

Hydra is designed for researchers and engineers working on compressible aerodynamics, conjugate heat transfer, reactive flows, and aerothermal problems.

---

## Solver Ecosystem

Each Hydra solver is an independent library with its own repository, documentation, and verification suite. Hydra integrates them as git submodules and exposes them through the `figaro` CLI.

| Solver | Physics | Status | Documentation |
|--------|---------|--------|---------------|
| [**MOSE**](https://github.com/open-hydra/MOSE) | Compressible Euler & Navier–Stokes, turbulence (RANS), finite-rate chemistry, soot | Active | [Docs](https://open-hydra.github.io/MOSE/) |
| [**FUSS**](https://github.com/open-hydra/FUSS) | Solid heat conduction (steady & transient) | Active | [Docs](https://open-hydra.github.io/FUSS/) |
| [**GROOT**](https://github.com/open-hydra/GROOT) | Radiative heat transfer | Active | [Docs](https://open-hydra.github.io/GROOT/) |
| [**ARES**](https://github.com/open-hydra/ARES) | Real fluid dynamics | Active | [Docs](https://open-hydra.github.io/ARES/) |

### Support Libraries

| Library | Role |
|---------|------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Field I/O and structured grid management |
| [FLINT](https://github.com/MarcoGrossi92/FLINT) | Fortran library for integrated thermochemistry |
| [OSLO](https://github.com/MarcoGrossi92/OSLO) | ODE solver library |
| [FiNeR](https://github.com/szaghi/FiNeR) | Fortran INI-style input file parser |

---

## Coupled Applications

Hydra builds coupled executables where two solvers exchange interface data at run-time.

| Executable | Coupling | Use Case |
|------------|----------|----------|
| `hydra-MF` | MOSE + FUSS | Conjugate heat transfer: compressible fluid + solid heat conduction |
| `hydra-AF` | ARES + FUSS | Conjugate heat transfer: real fluid + solid heat conduction |

---

## Key Features

- **Modular multi-physics** — plug-in solver architecture; each solver evolves independently with its own CI and test suite
- **Robust numerics** — multi-block structured grids, MUSCL reconstruction, low-Mach preconditioning
- **RANS turbulence closure** — several RANS models including Spalart-Allmaras and Shear-Stress-Tensor
- **Reactive flows** — finite-rate chemistry, chemical equilibrium (NASA CEA), soot formation
- **Conjugate heat transfer** — fluid–solid coupling with consistent wall-flux and wall-temperature exchange
- **Multiphase flows** *(under development)* — Lagrangian particle tracking and Eulerian two-fluid formulations for dispersed-phase simulations
- **OpenMP + MPI parallelism** — parallel execution out of the box
- **`figaro` CLI** — a single command-line tool to configure, build, run, and manage all Hydra workflows
- **Modern Fortran + C/C++** — built with CMake ≥ 3.23, compatible with GFortran and Intel oneAPI compilers

---

## Quick Start

### Prerequisites

| Tool | Minimum version |
|------|-----------------|
| CMake | 3.23 |
| Fortran compiler | GNU `gfortran` or Intel `ifort` / `ifx` |
| C/C++ compiler | `gcc` / `g++` or Intel `icx` |
| Conda | Any recent version |
| Git | 2.x |

OpenMP support is enabled by default. MPI is optional.

### 1 — Clone the whole repo

```bash
git clone && cd hydra
```

### 2 — Bootstrap

Run the bootstrap script once to configure the environment and install the `figaro` shell function:

```bash
./bootstrap.sh
```

### 3 — Build a coupled executable

```bash
# Fluid–solid conjugate heat transfer
figaro build hydra-MF
```

Compiled binaries land in `bin/`.

### 4 — Run a case

```bash
# Prepare case with ATLAS pre-processor
figaro ATLAS GPB ICB BCB

# Launch the coupled solver
figaro hydra-MF solve
```

---

## The `figaro` CLI

`figaro` is Hydra's unified command-line interface. All common tasks — from environment checks to running simulations — go through it.

```
figaro <command> [options]
```

| Command | Description |
|---------|-------------|
| `doctor` | Run environment and toolchain sanity checks |
| `compile <solver>` | Configure and compile a Hydra solver |
| `build <solver>` | Build without reconfiguring |
| `ATLAS <type>` | Run the ATLAS pre-processor |
| `<solver> solve` | Launch a solver |
| `kill` | Terminate a running Hydra process |

---

## Documentation

Full documentation is available at **[open-hydra.github.io/hydra](https://open-hydra.github.io/hydra/)**.

| Section | Content |
|---------|---------|
| [Getting Started](https://open-hydra.github.io/hydra/getting-started/) | Installation, bootstrap, first run |
| [User Guide](https://open-hydra.github.io/hydra/user/) | Figaro CLI reference, input format, running cases |
| [Theory Guide](https://open-hydra.github.io/hydra/theory/) | Physics split and solver coupling overview |
| [Verification & Validation](https://open-hydra.github.io/hydra/vv/) | Test cases and benchmark results |
| [Developer Guide](https://open-hydra.github.io/hydra/development/) | Repository structure, contributing, testing |

Solver-specific documentation lives in the individual solver repositories (see [Solver Ecosystem](#solver-ecosystem) above).

---

## Authors & Contributors

See [AUTHORS.md](AUTHORS.md) for the full list of contributors.

**Original authors:** Marco Grossi · Vincenzo Barbato · Alessio Sereno · Alessandro Montanari · Alex Falco · Giacomo Passarani · Andrea Giacomi · Gianluca Cocirla · Marco Fabiani

---

## Citing Hydra

If you use Hydra in your research, please cite it using the metadata in [CITATION.cff](CITATION.cff) or the reference below:

```bibtex
@software{hydra,
  author  = {Grossi, Marco and Barbato, Vincenzo and Sereno, Alessio and
             Montanari, Alessandro and Falco, Alex and Passarani, Giacomo and
             Giacomi, Andrea and Cocirla, Gianluca and Fabiani, Marco},
  title   = {hydra: A multi-physics {CFD} framework for coupled compressible
             flow, solid heat transfer, and reactive simulations},
  url     = {https://github.com/open-hydra/hydra},
  license = {GPL-3.0}
}
```

---

## License

Hydra is distributed under the **GNU General Public License v3.0**. See [LICENSE](LICENSE) for the full text.
