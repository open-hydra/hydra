# igloo-burnout-box — the IGLOO twin of test/ice-burnout-box

`test/igloo-burn-box-beckstead` without a residue, the burning section of `test/ice-burnout-box` line for line: no
`d-stop`, `d-burnout = 3.0e-6`. The four streams burn on the d^n law; at 3 µm what is left of a drop is handed to
the gas whole, so everything a stream carries in ends in the gas.

Closed forms and checks: `test/igloo-burn-box/check.py` and `test/igloo-burn-box/README.md` (`M = SUM mdot`,
`Fx = SUM mdot u_inj`, `Q = SUM mdot (h_c(T0) + u_inj^2/2 + dh)`).

`INPUT/` is what ATLAS writes for `input.ini` (`GPB BCB ICB`); `MESH/` is `test/evap-box-mom`'s. Gate:
`test/mi2-checks/p25_burnout_check.sh`, which also compares the energy per kg handed with ICE's.
