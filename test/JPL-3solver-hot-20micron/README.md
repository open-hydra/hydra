# JPL-3solver-hot-20micron

MOSE + ICE (`partE`) + IGLOO (`partL`) on the nozzle of `JPL-3solver-20micron`, at T0 = 2400 K with the cold case's
20 µm alumina (diameter) in both condensed phases and solidification on in both: the drops enter liquid (about 70 K
above T-melt) and cool to T-nuc. The supercooling is cut from both solvers' default (0.8 T-melt, 465 K) to 100 K
(`T-nuc = 2227`), a choice without data, so that the drops nucleate inside the nozzle; they leave partly frozen (no end
of the plateau, no solid exit). Gas, mesh, loading (`krho = 0.15` per phase), laws and schedule are the cold case's; the
differences are `T0` and the keys `solidification`, `T-melt`, `h-fus`, `T-nuc`, `cp-solid`, identical in
`[GPB-Phase2]` (ICE) and `[GPB-Phase3]` (IGLOO). No reference data: ICE and IGLOO are compared with each other.

Inputs: `test/mi2-checks/regen_case.sh <case> <scratch> GPB BCB ICB --install` with a gfortran build of the pinned ATLAS
(`ATLASDIR` set to it). Run: `figaro hydra-MI2 -p 8 solve`, or `test/mi2-checks/run_case.sh`. Compare:
`python3 ../JPL-3solver-hot-2micron/compare_freezing.py <run-dir> --design partial`. No gate yet: the gas does not
meet the planned gate's steadiness witnesses (below).

The comparison leaves out IGLOO's mollifier support (+-`mollify-passes` cells, 8) around the fronts: see
`../JPL-3solver-hot-2micron/README.md`.

Design estimate (a quasi-1D two-phase march, not a reference): nucleation near x = 0.116 m (axis 0.110–0.122 m, wall
row 0.103–0.105 m), exit frozen fraction about 0.28.

First runs (hydra 2afbfc8, 8 threads, 2026-11-18; no reference is stored: the gas does not meet W1 and W1b, its
axis exit Mach number moving by up to 1.6e-6 from one coupling step to the next): nucleation on the axis at x = 0.120 m
(both solvers), inner exit frozen fraction 0.243 (both), below the estimate's 0.267–0.294.
