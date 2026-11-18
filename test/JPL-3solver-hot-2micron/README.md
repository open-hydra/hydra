# JPL-3solver-hot-2micron

MOSE + ICE (`partE`) + IGLOO (`partL`) on the nozzle of `JPL-3solver-20micron`, at T0 = 2500 K with 2 µm alumina
(diameter) in both condensed phases and solidification on in both: the drops enter liquid (about 170 K above T-melt),
cool to T-nuc = 1861.6 K, jump to T-melt, freeze on the plateau and leave solid. Gas, mesh, loading (`krho = 0.15` per
phase), laws and schedule are the cold case's; the differences are `T0`, `dp` and the keys `solidification`, `T-melt`,
`h-fus`, `T-nuc`, `cp-solid`, identical in `[GPB-Phase2]` (ICE) and `[GPB-Phase3]` (IGLOO). No reference data: ICE and
IGLOO are compared with each other.

Inputs: `test/mi2-checks/regen_case.sh <case> <scratch> GPB BCB ICB --install` with a gfortran build of the pinned ATLAS
(`ATLASDIR` set to it). Run: `figaro hydra-MI2 -p 8 solve`, or `test/mi2-checks/run_case.sh`. Compare:
`python3 compare_freezing.py <run-dir>`. Gate: `test/mi2-checks/p14_hot_nozzle_check.sh`.

Design estimate (a quasi-1D two-phase march, not a reference): nucleation near x = 0.10 m (wall row first, from
0.088 m; axis 0.099–0.106 m), end of the plateau 0.110–0.123 m on the axis, solid at the exit near 1680 K.
