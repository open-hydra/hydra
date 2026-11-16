# igloo-burn-box-direct — the IGLOO twin of test/ice-burn-box

`test/igloo-burn-box` with the burning section of `test/ice-burn-box`, line for line (`rho = 1000`, `h0`,
`particle-mechanism = a-n2`, `combustion = direct`, `q-comb = 1.0e6`, `d-stop = 1.0e-6`, `T-burn = 300`; the latent
heat, `[IGLOO-Properties] Lv = 2.26e6`, is the solver section's as in that case), in the gas of that case (N2 at
400 K). The four streams enter at their burning temperature, burn
on the d2 law down to the residue diameter keeping that temperature, and stay in the box as inert drops of 1.25e-4
of their mass until `time-end`.

Closed forms and checks: `test/igloo-burn-box/check.py`, as for `test/igloo-burn-box-beckstead`; in `law` mode
`d^2 = d0^2 - K t` with `K = 4 Nu k_g ln(1 + B)/(rho c_pg)`, `B = (c_pg (T_g - T_p) + q)/L_v`, `Nu = 2` (the
variant run enters the drops at 1e-6 m/s; Ranz-Marshall's excess over 2 there is the tolerance), and the parcel's
temperature unchanged.

`INPUT/` is what ATLAS writes for `input.ini` (`GPB BCB ICB`); `MESH/` is `test/evap-box-mom`'s. Gates:
`test/mi2-checks/p22_igloo_burn_box_check.sh` and `test/mi2-checks/p23_burn_twins_check.sh`.
