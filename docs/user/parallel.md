# Parallel Execution and Scalability

`hydra-MF` couples **MOSE** (compressible flow) and **FUSS** (solid conduction) in a single conjugate-heat-transfer (CHT) executable. It uses **hybrid MPI + OpenMP parallelism**:

- MPI distributes each phase across processes.
- OpenMP parallelises work within each MPI process.

## Key recommendation

> **Use one MPI rank per socket, with one OpenMP thread per physical core on that socket.**
>
> On a 4-socket × 20-core node, the default layout is therefore `4 MPI ranks × 20 OpenMP threads`.

The most important decision is **socket coverage**. Within the family of layouts that covers every socket, changing the rank/thread split typically changes performance by only 1–3%. Failing to cover all sockets can cost about 30%.

This page explains how to select a configuration, launch `hydra-MF`, verify CPU placement, interpret the scaling results, and measure performance reproducibly.

## How one coupling step works

A coupled step consists of:

```text
1 FUSS Runge-Kutta step
+ 1 MOSE Runge-Kutta step
+ 2 full fluid-solid interface exchanges
```

There are **two interface exchanges**, not one, because `BC_Coupled_Fluxes` runs once per solver sub-iteration. Both solvers are present on every MPI rank; there is no communicator split. Each solver partitions its own phase independently over the same `MPI_COMM_WORLD`.

---

## 1. Select a parallel configuration

Four quantities determine a suitable configuration:

| Symbol | Description | How to obtain it |
| --- | --- | --- |
| `S` | Sockets (NUMA domains) per node | `numactl --hardware` |
| `C` | Physical cores per socket | `lscpu` |
| `B_f` | Mesh blocks in the **fluid** phase | MDB decomposition log |
| `B_s` | Mesh blocks in the **solid** phase | MDB decomposition log |

For a core budget `N`, apply these rules in order:

1. **Cover every socket.** Use at least one MPI rank per socket. This is the dominant performance decision.
2. **Keep each OpenMP team within one socket.** Choose a thread count that divides `C`.
3. **Decompose both phases.** MOSE and FUSS are block-atomic, so each phase needs at least as many blocks as there are MPI ranks. If a phase has fewer blocks than ranks, some ranks are idle for that phase.

### Recommended layouts for a 4-socket × 20-core node

| CPU cores | Nodes | MPI ranks | Threads/rank | Configuration |
| ---: | ---: | ---: | ---: | --- |
| 20 | 1 | 4 | 5 | `4 × 5` |
| 40 | 1 | 4 | 10 | `4 × 10` |
| 80 | 1 | 4 | 20 | `4 × 20` |
| 160 | 2 | 8 | 20 | `8 × 20` |
| 240 | 3 | 12 | 20 | `12 × 20` |

### Avoid these configurations

- A single MPI rank spanning multiple sockets.
- An OpenMP team crossing a socket boundary.
- **Pure MPI at high rank counts.** Its memory footprint grows sharply; see [Section 6](#6-memory).

---

## 2. Run `hydra-MF` in parallel

### Single-node execution

```bash
export OMP_NUM_THREADS=20
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp

mpirun -n 4 -ppn 4 ./bin/hydra_MF
```

This launches 4 MPI ranks × 20 OpenMP threads = 80 CPU cores, with one rank per socket.

### Multi-node execution

Keep the **per-node** layout unchanged and add whole nodes:

```text
Node 1: 4 MPI ranks × 20 threads = 80 cores
Node 2: 4 MPI ranks × 20 threads = 80 cores
Total:  8 MPI ranks × 20 threads = 160 cores
```

```bash
export OMP_NUM_THREADS=20
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp

mpirun -n 8 -ppn 4 ./bin/hydra_MF
```

For larger runs, increase the core count by adding identically configured nodes rather than changing the layout on nodes that are already in use.

---

## 3. Strong scaling

**Strong scaling** measures how execution time changes when the **same problem** is solved on more cores.

The reference benchmark is a flat-plate CHT case containing **6,912,000 cells**:

- 5,184,000 gas cells
- 1,728,000 plate cells

The benchmark uses nodes with 4 sockets × 20 physical cores.

### Results

| Cores | Nodes | MPI × OpenMP | Wall/step | Speed-up | Parallel efficiency |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | `1 × 1` | 38.768 s | 1.0 | 100.0% |
| 4 | 1 | `4 × 1` | 9.7304 s | 4.0 | 99.6% |
| 8 | 1 | `4 × 2` | 4.9990 s | 7.8 | 97.8% |
| 16 | 1 | `4 × 4` | 2.6485 s | 15.5 | 97.2% |
| 20 | 1 | `4 × 5` | 2.1999 s | 19.4 | 96.8% |
| 40 | 1 | `4 × 10` | 1.2470 s | 37.8 | 94.5% |
| 80 | 1 | `4 × 20` | 0.83576 s | 68.9 | 86.1% |
| 160 | 2 | `8 × 20` | 0.43100 s | 132.5 | 82.8% |
| 240 | 3 | `12 × 20` | 0.28581 s | 201.0 | 83.7% |

<figure>

  {% include "user/images/hmf-strong.svg" %}

</figure>

### What the results mean

Efficiency remains above **94% through 40 cores (half a node)** and reaches **86% at 80 cores (a full node)**. Beyond one node, it is nearly flat at **83% on two nodes** and **84% on three nodes**.

The main bottleneck is **within-node memory bandwidth**, not inter-node communication. Efficiency drops by about eight percentage points from 40 to 80 cores, where the four sockets approach their memory-bandwidth limit, and then changes little as whole nodes are added.

> **Important:** the reported speed-ups and efficiencies are based on **core-cycles**, not wall-clock time alone. A single active core can run at a much higher turbo frequency than an all-core workload. Using core-cycles removes that processor frequency effect from the scaling metric. Wall-clock times are still reported as measured.

---

## 4. MPI, OpenMP, and hybrid execution

`hydra-MF` supports pure OpenMP, pure MPI, and hybrid execution. The following measurements were made in one job on one node, using the same 1-core run as the reference point.

| Cores | OpenMP `1 × N` | MPI `N × 1` | Hybrid `4 × N/4` |
| ---: | ---: | ---: | ---: |
| 8 | 98.9% | 95.3% | **97.2%** |
| 16 | 96.4% | 94.2% | **96.8%** |
| 20 | 91.4% | 92.5% | **96.2%** |
| 40 | 85.1% | 83.9% | **94.2%** |
| 80 | 79.0% | 70.1% | **86.0%** |

<figure>

  {% include "user/images/hmf-modes.svg" %}

</figure>

Up to one socket, the three execution modes are within a few percentage points of one another. Beyond one socket, the hybrid layout is clearly preferable, while pure MPI degrades the most.

### Interpreting this comparison

There are two important caveats:

1. **Below a full node, the three modes do not use the same number of sockets.** Pure OpenMP fills sockets one at a time, so at 20 cores it occupies one socket and one memory controller. MPI and hybrid runs are spread across all four sockets. Only the 80-core points are a like-for-like comparison of execution mode.
2. **Pure MPI also carries a larger memory footprint.** See [Section 6](#6-memory).

---

## 5. Choosing the rank/thread split

### Socket coverage dominates

MPI rank count and NUMA coverage are separate variables. At a fixed 40-core budget, with the rest of the node idle:

| Layout | Sockets | Wall/step | vs. best | Node memory |
| --- | ---: | ---: | ---: | ---: |
| `4 × 10` | **4** | 1.2448 s | best | 9 GB |
| `2 × 20` | 2 | 1.6346 s | +31.3% | 6 GB |
| `4 × 10` | 2 | 1.6576 s | +33.2% | 8 GB |
| `1 × 40` | 2 | 1.7246 s | +38.5% | 4 GB |
| `40 × 1` | 4 | 1.9351 s | +55.5% | 65 GB |

The two `4 × 10` rows are the controlled comparison: the MPI/OpenMP configuration and memory footprint are identical; only socket placement differs. Spreading the four ranks across four sockets is worth **33%** relative to placing them on two sockets because it exposes four memory controllers instead of two.

The `40 × 1` configuration covers all four sockets but is still slower because its memory footprint is much larger; see [Section 6](#6-memory).

### Within the socket-covered family, the split barely matters

At the same core count, different rank/thread splits that cover every socket are usually within a few percent:

| Budget | Layout | Wall/step | vs. best |
| ---: | --- | ---: | ---: |
| 80 | `4 × 20` | 0.8390 s | best |
| 80 | `8 × 10` | 0.8477 s | +1.0% |
| 80 | `16 × 5` | 0.8490 s | +1.2% |
| 80 | `80 × 1` | 1.1260 s | +34.2% |
| 160 | `16 × 10` | 0.4218 s | best |
| 160 | `8 × 20` | 0.4331 s | +2.6% |
| 240 | `12 × 20` | 0.2864 s | best |
| 240 | `24 × 10` | 0.2898 s | +1.2% |

Within the recommended family, the spread is only about 1–3%, and the best split is not consistent: fewer ranks win at 80 and 240 cores, while more ranks win at 160 cores.

> **Do not over-tune the rank/thread split. Cover all sockets first, then move on.** Pure MPI is the exception: it is a materially different configuration and is about 34% slower at 80 cores.

---

## 6. Memory

Both solvers allocate domain-wide structures on **every MPI rank**. As a result, the total node footprint grows with rank count rather than shrinking with decomposition.

| Layout | Cores | Node memory |
| --- | ---: | ---: |
| `1 × 80` | 80 | 4.09 GB |
| `4 × 20` | 80 | 9.21 GB |
| `8 × 10` | 80 | 12.68 GB |
| `16 × 5` | 80 | 24.39 GB |
| `80 × 1` | 80 | **159.37 GB** |

On a 190 GB node, `80 × 1` consumes **84% of available memory** for this 6.9-million-cell case and would leave little room for larger problems. The recommended `4 × 20` layout uses only 9.21 GB for the same work.

There are two practical consequences:

- Across nodes, the **per-node** footprint decreases as the domain is distributed: 8.57 GB/node at 160 cores and 8.14 GB/node at 240 cores.
- **Large-footprint configurations are not reproducible across nodes.** The 65 GB `40 × 1` case measured 1.34 s, 1.94 s, and 2.56 s on three nominally identical nodes, whereas the ~9 GB configurations agreed within 1.5%. Treat pure-MPI timings as node-specific.

---

## 7. Fluid/solid balance and coupling

The relative cost of the two physics remains stable over the measured range:

| Cores | Layout | MOSE (fluid) | FUSS (solid) | Coupling |
| ---: | --- | ---: | ---: | ---: |
| 1 | `1 × 1` | 94.2% | 5.6% | 0.2% |
| 40 | `4 × 10` | 94.2% | 5.6% | 0.4% |
| 80 | `4 × 20` | 94.2% | 5.4% | 0.7% |
| 160 | `8 × 20` | 94.1% | 6.7% | 0.6% |
| 240 | `12 × 20` | 92.3% | 7.9% | 0.6% |

The plate contains **25% of the cells** but accounts for only about **6% of the step time**. For this case, solid conduction is therefore much cheaper per cell than compressible Navier–Stokes. The decomposition should prioritise balancing the **fluid** phase; MDB reports the achieved balance separately for each phase.

> **The coupling is not a major scaling concern.** From 1 to 240 cores, it remains between 0.2% and 0.7% of the step. In one decomposition where 90% of the 57,600 interface entries have their two blocks on different ranks, the complete coupling still accounts for only 0.6% of the step.

The two phase decompositions should still be coordinated, but coupling is not where most of the runtime is spent.

---

## 8. Practical recommendations

### Recommended default

For a node with `S` sockets and `C` physical cores per socket, using the full node (`N = S × C`):

```bash
export OMP_NUM_THREADS=$C
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp

mpirun -n $S -ppn $S ./bin/hydra_MF
```

This gives one MPI rank per socket, with OpenMP threads filling that socket. For multi-node runs, keep this per-node layout unchanged and add whole nodes.

### Pre-run checklist

1. Decompose **both phases** with `[MDB-Phase1]` and `[MDB-Phase2]`, with at least as many blocks per phase as MPI ranks.
2. Cover **every socket**. Do not spend significant time tuning the rank/thread split after that.
3. Avoid **pure MPI** for large rank counts: it is slower, uses substantially more memory, and its timings are not portable between nodes.
4. Verify CPU placement with `I_MPI_DEBUG=4`. Do not rely on wrapping `mpirun` with `numactl` to establish the desired MPI/OpenMP placement.
5. For benchmarking, disable `res-diter`, `shell-diter`, and periodic output; keep `coupling-wall = true`; score the last timer window and use the minimum of repeated runs.

---

## 9. Measure parallel scaling

Use the solver's **internal timers**, not the wall-clock duration of the batch job.

### Enable timing

Set the following in **all three** I/O sections: `[MOSE-IO]`, `[FUSS-IO]`, and `[HYDRA-IO]`.

```text
timer-diter = 20
```

`[HYDRA-IO]` reports the coupled step using disjoint timing buckets for FUSS, MOSE, and coupling, together with the core-cycles consumed by the step:

```text
HYDRA Timing | Step 40 | 20 steps | wall/step 8.3576E-01 s | FUSS 4.5192E-02 s (  5.4 %)
             | MOSE 7.8730E-01 s ( 94.2 %) | coupling 5.6500E-03 s (  0.7 %)
             | of which MPI 0.0000E+00 s (  0.0 %)

HYDRA Cycles | core-cycles/step  1.7158E+11 | aggregate    205.30 GHz
```

### Read the timers correctly

Each timing bucket is the **maximum over ranks**, measured independently. The three buckets can therefore sum to slightly more than the reported wall-clock step. Treat them as **shares of the work**, not as an exact additive decomposition of wall time.

`core-cycles/step` is the total over all ranks and threads. Use it for efficiency measurements: under perfect scaling, core-cycles remain constant as core count increases, independent of the processor frequency chosen by power management.

### Benchmarking rules

- **Score the last timer window.** The first window includes setup transients.
- Disable periodic output and set `res-diter` and `shell-diter` beyond the run length; both default to 1.
- Keep `[HYDRA-IO] coupling-wall = true`; it gates the complete CHT coupling.
- Repeat each configuration and keep the minimum measurement.
- Measure all points of a scaling curve in the same job, and verify that the CPU placement you requested is the placement you actually received.

---

## 10. Known scalability limits

The benchmark identifies four practical limits:

1. **Filling a node costs more than crossing one.** Efficiency falls from 94.5% at 40 cores to 86.1% at 80 cores as socket memory bandwidth becomes the limit, then stays around 83–84% across two and three nodes.
2. **Pure MPI is a poor full-node configuration.** At 80 cores it reaches only 70.1% efficiency, consumes 159 GB, and its timing can vary by up to 1.9× between nominally identical nodes.
3. **Per-rank memory is not reduced by adding ranks.** Both solvers retain domain-wide allocations on every rank.
4. **Insufficient phase decomposition causes idle ranks.** If either phase has fewer mesh blocks than MPI ranks, some ranks are idle for that phase.

---

## 11. Quick reference

| Goal | Setting |
| --- | --- |
| Rank/thread layout | One rank per socket; threads = physical cores per socket |
| OpenMP placement | `OMP_PLACES=cores` |
| Intel MPI pinning | `I_MPI_PIN_DOMAIN=omp` |
| Explicit per-rank pinning | `I_MPI_PIN_DOMAIN="[0x...,0x...]"` (not `numactl`) |
| Verify placement | `I_MPI_DEBUG=4` or `ps -eLo psr,comm` |
| Coupled decomposition | `[MDB-Phase1]` + `[MDB-Phase2]` |
| Timing instrumentation | `timer-diter` in `[MOSE-IO]`, `[FUSS-IO]`, and `[HYDRA-IO]` |
| Keep coupling enabled | `[HYDRA-IO] coupling-wall = true` |
| Disable while benchmarking | `res-diter`, `shell-diter`, periodic output |

---

## Summary

For the measured `hydra-MF` case, the practical recipe is simple:

> **Cover every socket, keep OpenMP teams inside sockets, decompose both phases sufficiently, and avoid pure MPI.**

The benchmark shows that socket coverage matters far more than fine-tuning the MPI/OpenMP split, while the dominant scaling limit is memory bandwidth within a node rather than communication between nodes.
