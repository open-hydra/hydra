# ice-burn-box-beckstead — an ICE material burning on the d^n law

`test/ice-burn-box` with the other burning model: `[GPB-Phase2]` names `combustion = Beckstead` with `K-burn = 1.39e-6`,
`n-burn = 1.8`, `X-eff = 0.5`, `T-ign = 290`, `beta-part = 0.3`, `q-comb = 1.0e6`, `d-stop = 1.0e-6` — the same lines
as `test/igloo-burn-box-beckstead`, its twin in IGLOO — with the same reaction (`particle-mechanism = a-n2`) and
`h0 = 1.0e6`. `MESH/` and `INPUT/` are `test/ice-box`'s but for `part-phase.txt`, `part-models.txt` (the constants
above, the material's row), `part-products.txt` and `part-properties.dat`.

The drops (300 K) are above `T-ign`, so the frozen source burns: each drop loses `rho pi K-burn X-eff/(2 n) d^(3-n)`,
keeps `beta-part q-comb` per kg burnt and exchanges heat with the gas as an inert drop does. Closed forms and
checks: `test/ice-burn-box/check.py` (it reads the model from the case) and `test/ice-burn-box/README.md`, with

    Q = W + M (c T_p + |u_p|^2/2 + dh - beta q) - Qc,     Qc = n_p Nu pi d k_g (T_g - T_p) V.

Gates: `test/mi2-checks/p21_ice_burn_box_check.sh` (the case against its closed forms) and
`test/mi2-checks/p23_burn_twins_check.sh` (ICE against IGLOO).
