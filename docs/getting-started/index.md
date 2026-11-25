# Getting Started

Hydra is a multi-physics CFD framework built around specialized solver libraries under `src/lib/*`. Instead of one monolithic code, Hydra composes solvers focused on different physical regimes and orchestrates them through a common build and runtime workflow.

Welcome to the Hydra getting-started guide. This section helps you
build Hydra executables and run your first workflow.

!!! info "Hydra and solver docs"
    Hydra documentation focuses on framework-level orchestration and coupling.
    Detailed solver setup remains in each solver documentation (MOSE, FUSS).

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Build Hydra and its solver dependencies from source

    [:octicons-arrow-right-24: Installation](installation.md)

-   :material-rocket-launch:{ .lg .middle } __Quick Start__

    ---

     Build and run your first Hydra case

     [:octicons-arrow-right-24: Quick start](quick-start.md)

</div>

## Scope of This Section

This getting started guide covers:

1. **Installation**
    - Hydra build workflow
    - Common build/compile commands

2. **Quick Start**
    - Run a Hydra executable
    - Understand where Hydra input stops and solver input begins

## After Quick Start

Use these pages next:

- [Core Idea](#core-idea) and [Multi-Physics Coupling](#multi-physics-coupling), below, for architecture and solver roles
- [User Guide](../user/index.md) for Hydra-level input and run workflow
- Solver docs for detailed physics and solver-local settings

---

## Core Idea

Each solver in Hydra is designed for a specific class of equations and physical assumptions. This keeps each code focused while still allowing coupled workflows.

| Solver | Physics focus | Location | Documentation |
|--------|---------------|----------|----------------|
| **MOSE** | Compressible ideal-gas flow (Euler/Navier-Stokes), turbulence, reacting flows | `src/lib/MOSE` | [MOSE docs](https://open-hydra.github.io/MOSE/) |
| **FUSS** | Solid heat conduction (steady and transient) | `src/lib/FUSS` | [FUSS docs](https://open-hydra.github.io/FUSS/) |
| **GROOT** | Radiative heat transfer | `src/lib/GROOT` | [GROOT docs](https://open-hydra.github.io/GROOT/) |
| **ARES** | Real fluid dynamics | `src/lib/ARES` | [ARES docs](https://open-hydra.github.io/ARES/) |
| **ICE** | Eulerian dispersed condensed phase (drops, particles) | `src/lib/ICE` | [ICE docs](https://open-hydra.github.io/ICE/) |
| **IGLOO** | Lagrangian dispersed condensed phase (parcels) | `src/lib/IGLOO` | [IGLOO docs](https://open-hydra.github.io/IGLOO/) |

!!! info "Hydra and solver repositories"
    Solver libraries are versioned in their own repositories; Hydra requires a range of releases of each and checks one out at build time. Hydra-level docs summarize how solvers fit together; solver-level docs provide full theory and user references.

---

## Multi-Physics Coupling

Hydra also provides coupled applications where solvers exchange data during the run. Three coupled executables exist:

| Executable | Solvers | Coupling |
|------------|---------|----------|
| `hydra-MF` (`src/app/MF`) | MOSE + FUSS | fluid-solid wall fluxes and wall temperature (conjugate heat transfer) |
| `hydra-AF` (`src/app/AF`) | ARES + FUSS | the same, with a real fluid |
| `hydra-MI2` (`src/app/MI2`) | MOSE + ICE + IGLOO | a gas and its dispersed condensed phases: the gas state to the condensed solvers, their momentum, energy, mass and species sources back to the gas ([Multiphase Coupling](../user/mi2.md)) |

---

## Documentation Map

Use the Hydra docs for framework-level orientation, then jump to solver docs for deep technical details.

| Section | Focus |
|---------|-------|
| Getting Started (this section) | Build and run workflow in this repository |
| [User Guide](../user/index.md) | Hydra-level usage and input layout, plus links to solver user docs |
| [Theory Guide](../theory/index.md) | High-level physics split by solver, plus links to solver theory docs |
| [Developer Guide](../development/index.md) | Repository structure, testing and contribution workflow |

---

## Ecosystem Components

| Component | Role |
|----------|------|
| [ATLAS](https://github.com/open-hydra/ATLAS) | Pre-processing and case preparation tools |
| Hydra (`open-hydra/hydra`) | Multi-solver integration, orchestration, coupled executables |
| Solver repos (MOSE, FUSS, GROOT, ARES, et al.) | Physics-specific implementations and detailed docs |

For licensing details, see [License](../about/license.md).
