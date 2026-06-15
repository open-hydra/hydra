# Theory Guide

This section gives a Hydra-level map of the physics solved by each
solver library. Detailed derivations, closure models, and numerical
discretization details are maintained in each solver's own documentation.

## Solver Physics at a Glance

| Solver | Core equations | Typical domain | Detailed theory |
|--------|----------------|----------------|-----------------|
| **MOSE** | Compressible Euler/Navier-Stokes for ideal-gas mixtures, with optional turbulence and reacting-flow models | External/internal aerothermodynamics, reacting gas flows | [MOSE Theory](https://open-hydra.github.io/MOSE/theory/) |
| **FUSS** | Heat conduction equation in solids (steady/transient), including source terms and thermal BC models | Solid thermal response and conjugate heat transfer solid side | [FUSS Theory](https://open-hydra.github.io/FUSS/theory/) |
| **GROOT** | Radiation transport equation | Radiative heat transfer | [GROOT Theory](https://open-hydra.github.io/GROOT/theory/) |

## Coupled Theory in Hydra

Hydra-level coupled applications (such as `hydra-MF`) combine solver
domains and enforce interface consistency through coupling boundary
conditions, source terms, and data exchange at runtime.

At a high level, this introduces:

- Cross-solver boundary coupling (fluid-side and solid-side interface states)
- Shared coupling iteration control
- Optional time-accurate synchronized stepping

For implementation details of the current coupled executable, see the
developer-side source in `src/app/MF`.

## Recommended Reading Path

1. Start from [Overview](../overview.md) for architecture and solver roles.
2. Read the solver-specific theory pages for the solver you are using.
3. Return to Hydra docs for integration and coupled-run workflow.
