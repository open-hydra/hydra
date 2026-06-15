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
| **PESTE** | Pressure-enthalpy compressible/reacting formulation with tabulated equilibrium thermo-chemistry | `src/lib/PESTE` | [PESTE docs](https://open-hydra.github.io/PESTE/) |
| **GROOT** | Reserved solver slot (under development) | `src/lib/GROOT` | Not published yet |

!!! info "Hydra and solver repositories"
    Solver libraries are versioned in their own repositories and included in Hydra as submodules. Hydra-level docs summarize how solvers fit together; solver-level docs provide full theory and user references.

---

## Multi-Physics Coupling

Hydra also provides coupled applications where two solvers exchange interface data during the run. The current coupled executable in `src/app/MF` coordinates:

- Fluid-side quantities from MOSE
- Solid-side quantities from FUSS
- Coupled wall fluxes and wall-temperature consistency

This is the basis for conjugate heat-transfer style simulations in the Hydra workflow.

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
| Solver repos (MOSE/FUSS/PESTE) | Physics-specific implementations and detailed docs |

For licensing details, see [License](about/license.md).
