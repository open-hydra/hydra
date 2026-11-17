# steady-igloo-smoke

A burning IGLOO material (A => AL2O3(a) + N2) whose smoke is born in ICE's phase partE, MOSE and ICE steady on JPL's
schedule: IGLOO's smoke handed to ICE as a rate.

`MESH/mesh.tec`: `python3 test/mi2-checks/steady_case_prep.py mesh test/steady-igloo-smoke` (test/JPL-Eulerian-20micron's
nozzle, every 5th node line kept: 38 x 8 cells, 2D axisymmetric). `INPUT/`: ATLAS GPB ICB BCB of `input.ini`
(`test/mi2-checks/regen_case.sh test/steady-igloo-smoke SCRATCH GPB ICB BCB --install`). Gate:
`test/mi2-checks/p31_steady_check.sh` (the converged field does not depend on the local time steps).
