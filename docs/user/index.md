# User Guide

The User Guide describes how to configure and run Hydra workflows.
Hydra can run individual solvers and coupled applications; this page
summarizes the split and points to solver-specific user documentation.

<div class="grid cards" markdown>

-   :material-console:{ .lg .middle } __Figaro CLI__

    ---

    Command reference for build, compile, run, clean, and submodule workflows.

    [:octicons-arrow-right-24: Figaro reference](figaro.md)

-   :material-file-document-outline:{ .lg .middle } __Hydra Input__

    ---

    Hydra-level control blocks (`HYDRA-*`) and how they coexist with
    solver input blocks in a single `input.ini`.

    [:octicons-arrow-right-24: Hydra input reference](input.md)

-   :material-play-circle-outline:{ .lg .middle } __Running Cases in Hydra__

    ---

    Hydra execution workflow for coupled and single-solver runs.

    [:octicons-arrow-right-24: Run guide](using.md)

-   :material-water-outline:{ .lg .middle } __Multiphase Coupling__

    ---

    `hydra-MI2`: a gas with Eulerian (ICE) and Lagrangian (IGLOO)
    condensed phases, phase routing, coupling schedule, evaporation and combustion.

    [:octicons-arrow-right-24: hydra-MI2](mi2.md)

-   :material-chip:{ .lg .middle } __Parallel Execution__

    ---

    Choosing ranks × threads, NUMA placement, measured scaling, and batch scripts.

    [:octicons-arrow-right-24: Running in parallel](parallel.md)

</div>

## Solver-Specific User Docs

Use the Hydra docs for orchestration and coupling behavior, then use
the solver docs for complete input and workflow details.

| Solver | User documentation |
|--------|---------------------|
| MOSE | [MOSE User Guide](https://open-hydra.github.io/MOSE/user/) |
| FUSS | [FUSS User Guide](https://open-hydra.github.io/FUSS/user/) |
| GROOT | [GROOT User Guide](https://open-hydra.github.io/GROOT/user/) |
| ICE | [ICE User Guide](https://open-hydra.github.io/ICE/user/) |
| IGLOO | [IGLOO User Guide](https://open-hydra.github.io/IGLOO/user/) |

## Who This Section Is For

This section is intended for:

- Users running Hydra workflows for the first time
- Users preparing coupled runs based on MOSE and FUSS, or on MOSE, ICE and IGLOO
- Users needing Hydra-level controls beyond solver-local options

If you are new to Hydra, start from [Getting Started](../getting-started/index.md).
