# steady-ice-smoke

A burning ICE family (AL(L), Beckstead, d-stop) makes a smoke family (AL2O3(L)) in a hot nozzle, MOSE and ICE steady
on JPL's schedule: ICE's smoke births handed as a rate (`srate`).

`MESH/mesh.tec`: `python3 test/mi2-checks/steady_case_prep.py mesh test/steady-ice-smoke` (test/JPL-Eulerian-20micron's
nozzle, every 5th node line kept: 38 x 8 cells, 2D axisymmetric). `INPUT/`: ATLAS GPB ICB BCB of `input.ini`
(`test/mi2-checks/regen_case.sh test/steady-ice-smoke SCRATCH GPB ICB BCB --install`), then
`python3 test/mi2-checks/steady_case_prep.py empty test/steady-ice-smoke` (the smoke family at ICE's empty state,
1e-25 in every slot). Gate: `test/mi2-checks/p31_steady_check.sh` (the converged field does not depend on the local
time steps).
