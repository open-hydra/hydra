# JPL-3solver-hot-20micron

MOSE + ICE (`partE`) + IGLOO (`partL`) on the nozzle of `JPL-3solver-20micron`, at T0 = 2400 K with the cold case's
20 µm alumina (diameter) in both condensed phases and solidification on in both: the drops enter liquid (about 70 K
above T-melt) and cool to T-nuc. The supercooling is cut from both solvers' default (0.8 T-melt, 465 K) to 100 K
(`T-nuc = 2227`), a choice without data, so that the drops nucleate inside the nozzle; they leave partly frozen (no end
of the plateau, no solid exit). Gas, mesh, loading (`krho = 0.15` per phase), laws and schedule are the cold case's; the
differences are `T0` and the keys `solidification`, `T-melt`, `h-fus`, `T-nuc`, `cp-solid`, identical in
`[GPB-Phase2]` (ICE) and `[GPB-Phase3]` (IGLOO). No reference data: ICE and IGLOO are compared with each other.

Inputs: `test/MI2/regen_case.sh <case> <scratch> GPB BCB ICB --install` with a gfortran build of the pinned ATLAS
(`ATLASDIR` set to it); `test/MI2/validation/JPL-3solver-hot-2micron/hot_nozzle.py g1 <case> 20um` checks them. Run: `figaro hydra-MI2 -p 8
solve`, or `test/MI2/run_case.sh`. Compare: `python3 ../JPL-3solver-hot-2micron/compare_freezing.py <run-dir>
--design partial`.

Stopped (2026-11-19): no reference, no leg of p14 (`P14_VARIANT=20micron test/MI2/validation/JPL-3solver-hot-2micron/check.sh`
judges it on demand, red). At one thread the gas does not settle: W1 7.8e-8 and W1b 7.1e-7 in blocks of 500 MOSE
iterations, 1.7e-7 and 1.5e-6 in blocks of 250, against 1e-8 and 1e-7 (the cold case: 3.6e-9, 2.5e-9). The gas exit
Mach number on the axis switches between two states about 6.5e-7 apart at irregular coupling steps; ICE's velocity in
the axis exit cell moves with it, IGLOO's exit records do not, the cold case shows none of it, and the comparison reads
the same numbers in a block that switches and in one that does not. An A/B run of this case with ICE's source to the
gas held fixed (or with IGLOO alone) would attribute it. Record: `plan-bucket/hot-nozzle-solidification-plan.md` (the
design, the witnesses and their bounds) and the pull request that brought p14.

The comparison leaves out IGLOO's mollifier support (+-`mollify-passes` cells, 8) around the fronts: see
`../JPL-3solver-hot-2micron/README.md`.

Design estimate (a quasi-1D two-phase march, not a reference): nucleation near x = 0.116 m (axis 0.110–0.122 m, wall
row 0.103–0.105 m), exit frozen fraction about 0.28 (on the axis).

One-thread run (hydra 2afbfc8; ICE / IGLOO): nucleation on the axis at x = 0.1203 / 0.1199 m; exit frozen fraction
0.2474 in the axis cell (both), 0.2432 / 0.2434 over the inner region, below the estimate's 0.267–0.294 (the 2-D gas at
the exit is 1325 K on the axis against the estimate's 1026 K: the drops freeze less).
