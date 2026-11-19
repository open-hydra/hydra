# JPL-3solver-hot-2micron

MOSE + ICE (`partE`) + IGLOO (`partL`) on the nozzle of `JPL-3solver-20micron`, at T0 = 2500 K with 2 µm alumina
(diameter) in both condensed phases and solidification on in both: the drops enter liquid (about 170 K above T-melt),
cool to T-nuc = 1861.6 K, jump to T-melt, freeze on the plateau and leave solid. Gas, mesh, loading (`krho = 0.15` per
phase), laws and schedule are the cold case's; the differences are `T0`, `dp` and the keys `solidification`, `T-melt`,
`h-fus`, `T-nuc`, `cp-solid`, identical in `[GPB-Phase2]` (ICE) and `[GPB-Phase3]` (IGLOO). No reference data: ICE and
IGLOO are compared with each other.

Inputs: `test/mi2-checks/regen_case.sh <case> <scratch> GPB BCB ICB --install` with a gfortran build of the pinned ATLAS
(`ATLASDIR` set to it); `test/mi2-checks/p14_hot_nozzle.py g1 <case> 2um` checks them. Run: `figaro hydra-MI2 -p 8
solve`, or `test/mi2-checks/run_case.sh`. Compare: `python3 compare_freezing.py <run-dir>`. Gate:
`test/mi2-checks/p14_hot_nozzle_check.sh [scratch-root] [exe]` (one thread, 45–50 min on a pcprop): it runs a copy
of the case with two exit probes of the gas Mach number added (`p14_hot_nozzle.py inject`); `P14_RUN=<run-dir>` judges
a stored run, which needs those probe files. The stored reference (`OUTPUT/`, `REFERENCE-STATE.md`) is such a run at one
thread: byte identity holds at one thread only.

The gas does not settle to the cold case's level. At one thread W1 (the MOSE density residual's maximum in each of the
last 10 blocks over its first value) is 4.9e-8 against 3.6e-9 for the cold case (blocks 41–49 fall x0.83 per block to
1.1e-8, block 50 reads 2.7e-8; in blocks of 250 iterations it stays at 1.0–1.9e-8), W1b (the axis exit Mach range over
the last block) 2.0e-8 against 2.5e-9. The user's decision of 2026-11-19: "print the steadiness figures for 2 um without
judging them, and accept the ICE-against-IGLOO comparison on a gas that still moves slightly" — every comparison number
is the same to its printed digits at one and at eight threads. The gate prints W1 and W1b with these values and the
bounds 1e-8 and 1e-7.

IGLOO's Eulerian field is mollified over +-`mollify-passes` cells (`[IGLOO-General]`, 8 by default; +-1 more with
`gas-order` 2), ICE's is not: `compare_freezing.py` leaves out max(3, passes) cells upstream of IGLOO's first plateau
cell in T1, the passes read from the case's `input.ini` (on a 3-coupling-step smoke at one thread the plan's 3-cell set
gives delta 0.045 with 8 passes and 0.006 with `mollify-passes = 0`, bound 0.031), and judges IGLOO's plateau
temperature (T3) on its parcel rows, which are not smoothed; the Eulerian T_p is printed (the 2 um plateau, about 12
cells, is shorter than the 19-cell box).

Design estimate (a quasi-1D two-phase march, not a reference): nucleation near x = 0.10 m (wall row first, from
0.088 m; axis 0.099–0.106 m), end of the plateau 0.110–0.123 m on the axis, solid at the exit near 1680 K.

The reference run (hydra 74a90a4, one thread; ICE / IGLOO): nucleation on the axis at x = 0.1062 / 0.1067 m, end of the
plateau 0.1213 / 0.1219 m, every drop solid at the exit, inner exit particle temperature 1852 / 1855 K (1860 / 1864 K
over the whole exit section); the comparison passes its 9 checks. The estimate is a single stream for the whole
section; the 2-D drops sit in the hot core of the gas (inner exit gas 1722 K against the estimate's 1525 K) and release
0.89 times the estimate's heat, the factor the cold case shows.
