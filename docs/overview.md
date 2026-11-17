---
title: Overview
---

# Overview

Hydra is a multi-physics CFD framework built around specialized solver libraries under `src/lib/*`. Instead of one monolithic code, Hydra composes solvers focused on different physical regimes and orchestrates them through a common build and runtime workflow.

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
    Solver libraries are versioned in their own repositories and included in Hydra as submodules. Hydra-level docs summarize how solvers fit together; solver-level docs provide full theory and user references.

---

## Multi-Physics Coupling

Hydra also provides coupled applications where solvers exchange data during the run. Three coupled executables exist:

| Executable | Solvers | Coupling |
|------------|---------|----------|
| `hydra-MF` (`src/app/MF`) | MOSE + FUSS | fluid-solid wall fluxes and wall temperature (conjugate heat transfer) |
| `hydra-AF` (`src/app/AF`) | ARES + FUSS | the same, with a real fluid |
| `hydra-MI2` (`src/app/MI2`) | MOSE + ICE + IGLOO | a gas and its dispersed condensed phases: the gas state to the condensed solvers, their momentum, energy, mass and species sources back to the gas ([Multiphase Coupling](user/mi2.md)) |

---

## Documentation Map

Use the Hydra docs for framework-level orientation, then jump to solver docs for deep technical details.

| Section | Focus |
|---------|-------|
| [Getting Started](getting-started/index.md) | Build and run workflow in this repository |
| [User Guide](user/index.md) | Hydra-level usage and input layout, plus links to solver user docs |
| [Theory Guide](theory/index.md) | High-level physics split by solver, plus links to solver theory docs |
| [Developer Guide](development/index.md) | Repository structure, testing and contribution workflow |

---

## Ecosystem Components

| Component | Role |
|----------|------|
| [ATLAS](https://github.com/open-hydra/ATLAS) | Pre-processing and case preparation tools |
| Hydra (`open-hydra/hydra`) | Multi-solver integration, orchestration, coupled executables |
| Solver repos (MOSE, FUSS, GROOT, ARES, et al.) | Physics-specific implementations and detailed docs |

For licensing details, see [License](about/license.md).
