# igloo-burn-box-beckstead — the IGLOO twin of test/ice-burn-box-beckstead

`test/igloo-burn-box` with the burning section of `test/ice-burn-box-beckstead`, line for line (`rho = 1000`, `h0`,
`particle-mechanism = a-n2`, the law's six keys and `d-stop = 1.0e-6`), in the gas of that case (N2 at 400 K). The
four streams burn on the d^n law down to the residue diameter and stay in the box as inert drops of 1.25e-4 of their
mass until `time-end`.

Closed forms and checks: `test/igloo-burn-box/check.py` (it reads `d-stop` and the end records of
`OUTPUT/part-outloc-A-sweep*.dat`): `M = SUM mdot (1 - r)`, `Q = SUM mdot [h_c(T0) + u^2/2 - r h_c(T_end) + (1 - r)
dh]`, `r = (d_end/d0)^3`. With the word `law` as fourth argument the same script checks the law itself on a run cut
inside the burn: `d^n = d0^n - K-burn X-eff t` at `time-end`.

`INPUT/` is what ATLAS writes for `input.ini` (`GPB BCB ICB`); `MESH/` is `test/evap-box-mom`'s. Gates:
`test/mi2-checks/p22_igloo_burn_box_check.sh` and `test/mi2-checks/p23_burn_twins_check.sh`.
