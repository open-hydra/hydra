# Reference state (2026-09-30)

`OUTPUT/` + `run.log`: `test/mi2-checks/bin/hydra_MI2.G-2026-09-30` (md5 `880ab049cd343a7e4bb4bfd94b65e167`, the sources of 2026-09-30;
the log header stamps the checked-out trees IGLOO `e5faf0a`, ICE `494668b`), INI and inputs as shipped, **1 thread**, 7.7 s, exit 0.
Inputs from ATLAS `ebeaf4f` (gfortran build), reproducible byte for byte with `test/mi2-checks/regen_case.sh`.

| quantity | value |
|---|---|
| mass, species 1 injected = gas gained | `1.495472E-08` kg, rel.err **`2.212E-16`** |
| energy handed to gas (sourceEn + datum) | `-2.36203438E-01` J, of which datum `-2.98989125E-04` J |
| energy rel.err (bookkeeping, datum-blind) | **`2.790E-09`** |
| mass-weighted gas T | `350.000000 -> 349.050909`, delta `-9.491E-01` K |
| SA variable `mi_tilde` at iteration 950 (t 3.7387e-3 s) | `5.990428E-05 .. 6.072869E-05` (initial `6.0E-05`) |
| eddy viscosity ratio `mit/mil` | `0.1988 .. 0.2214` |

The balance lines are those of `evap-box-ta` on the same binary to the printed digits except the last two of the energy
handed (`-2.36203431E-01` there). The run ends at iteration 951 (the first coupling step converges at its first iteration),
so the stored field is the write at iteration 950. Byte identity only at `OMP_NUM_THREADS=1`.
