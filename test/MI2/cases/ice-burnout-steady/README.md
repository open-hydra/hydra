# ice-burnout-steady

A burning ICE material without a smoke (A => N2) whose drops reach d-burnout (17 um) inside the nozzle, MOSE and ICE
steady on JPL's schedule. A cell whose drops burn and are smaller than d-burnout keeps its own balance, and a
fraction f of what flows out of it burns at its exit, handed to the gas as every burn is: f rises from 0 at d-burnout
to 1 at the band's lower edge (ICE's `burnout_band` = [0.90, 1], in multiples of d-burnout; the log prints it; a
time-accurate run empties a cell at or below d-burnout at once). ICE's cfl is
0.2, not JPL's 0.4. The case converges to round-off at first order and with MUSCL-RK3 (as shipped): ICE's density
residual at the ends of the ICE blocks stays below 1e-14 after 20 coupling cycles (25 at cfl 0.1).

`MESH/mesh.tec`: `python3 test/MI2/cases/ice-burnout-steady/prep.py mesh test/MI2/cases/ice-burnout-steady` (test/MI2/validation/JPL-Eulerian-20micron's
nozzle, every 5th node line kept: 38 x 8 cells, 2D axisymmetric). `INPUT/`: ATLAS GPB ICB BCB of `input.ini`
(`test/MI2/regen_case.sh test/MI2/cases/ice-burnout-steady SCRATCH GPB ICB BCB --install`). Check:
`test/MI2/run.sh ice-burnout-steady` (`check.sh`: the case at ICE cfl 0.2 and 0.1, the converged fields equal within
1e-12 of each variable's peak).
