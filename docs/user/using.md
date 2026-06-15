# Running Solver Cases in Hydra

This page describes the Hydra-level runtime workflow: how to launch
solver executables from the framework and how to interpret the split
between framework and solver responsibilities.

For full command options and examples, see [Figaro CLI Reference](figaro.md).

**What Hydra controls**

Hydra-level behavior is configured in `input.ini` via `HYDRA-*` sections:

- Coupling iteration scheduling
- Runtime synchronization in coupled runs
- Coupling output cadence and format

See [Hydra Input](input.md) for the section-level reference.

**What each solver controls**

Physics models, numerics, and solver I/O remain solver-specific and are configured in solver-level sections. Use solver docs for complete parameter semantics.

## Typical execution modes

### A. Coupled Hydra executable

For coupled workflows (for example fluid-solid coupling):

```bash
figaro hydra-MF solve
```

This executes the Hydra application that coordinates multiple solvers
inside one coupling loop.

### B. Single solver through Hydra CLI

You can also launch an individual solver from the same command surface:

```bash
figaro MOSE solve
```

## Case layout

A run directory typically contains:

```text
my_case/
├── input.ini
├── INPUT/
└── OUTPUT/
```

The exact content of `INPUT/` depends on the solver(s) in use.

## Runtime diagnostics

A successful run should show:

- input loading status
- solver setup summary
- residual/time progression
- periodic output generation according to configured frequencies
