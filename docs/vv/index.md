# Verification & Validation

Hydra contributes no physics of its own. MOSE, FUSS, GROOT and ARES each solve their own equations and each carry their own V&V suite; what Hydra adds is the **coupling** — the exchange of interface state between two solvers inside one executable, and the iteration that drives both to a consistent steady solution.

This section therefore verifies the coupling, not the solvers.

<div class="grid cards" markdown>

-   :material-thermometer-lines:{ .lg .middle } __Conjugate Flat Plate__

    ---

    A laminar boundary layer over a conducting plate heated from below, where the interface temperature is produced by the coupling rather than prescribed. Run through both coupled executables and compared with an analytical correlation.

    [:octicons-arrow-right-24: Conjugate flat plate](conjugate-flat-plate.md)

</div>

---

## What is verified here

Three things, in increasing order of strength:

**Interface continuity.** Each side of a coupled interface computes its own wall state. Where the two sides must agree, they are checked against each other — this is the coupling's core invariant.

**Physical response.** The coupled solution is compared against an independent analytical reference over a parameter sweep wide enough to span the limiting regimes, so the check is a trend rather than a single point.

**Solver-agnosticism.** The coupling layer is instantiated once per solver pair (`src/app/MF`, `src/app/AF`). Running the same physical problem through both pairs cross-verifies the instantiations against each other. This check exists only at the Hydra level; no single-solver suite can perform it.

---

## Scope

| Level | Covered by |
|---|---|
| Interface exchange, coupling loop, coupled executables | this section |
| Fluid and solid physics, numerics, turbulence, radiation | the solver suites below |
| Parallel efficiency of the coupled run | [Parallel Execution](../user/parallel.md) |
| Which checks to run against which change | [Testing](../development/testing.md) |

Solver-level verification lives in the solver repositories:
[MOSE](https://open-hydra.github.io/MOSE/vv/) ·
[FUSS](https://open-hydra.github.io/FUSS/vv/) ·
[ARES](https://open-hydra.github.io/ARES/vv/) ·
[GROOT](https://open-hydra.github.io/GROOT/vv/)
