# igloo-smoke-box — the smoke of a burning IGLOO material, born in an ICE family

`test/igloo-burnout-box` with a smoke: the four streams of A burn by `A => AL2O3(a) + 3.6396 N2` (`a-smoke.yaml`),
half of the mass into N2 and half into alumina AL2O3(a), a material of the eulerian phase `partE` (ICE, one MK
family, empty at the start), born at `d-smoke` 3 µm with the burning parcels' velocity and temperature. Cycle 1: six
MOSE blocks with IGLOO's sources, the smoke of each block joining ICE after it; cycle 2: one ICE step (NoDrag, NoHeat:
the field keeps the birth).

Closed forms and checks: `check.py` (K0-K7: the smoke's mass, the gas's species, momentum and energy of the gas, the
smoke and the streams, the smoke's temperature at birth, its number). `INPUT/` is what ATLAS writes for `input.ini`
(`GPB ICB BCB`); `MESH/` is `test/igloo-burnout-box`'s. Gate: `test/mi2-checks/p29_igloo_smoke_check.sh`.
