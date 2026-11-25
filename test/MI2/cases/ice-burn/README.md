# ice-burn — a burning ICE material handed to the gas

A gate for what `hydra-MI2` hands the gas when an ICE material burns (`combustion = direct`): the mass
(`-source(1)·vol`, as the gas species of `INPUT/part-products.txt`) and the energy of the burnt mass moved to the
gas's enthalpy scale, the heat of reaction inside it, against closed forms that follow from `input.ini` and the
`INPUT/` files.

## The case

`test/MI2/cases/ice-exchange` with one change: its material burns. `MESH/` and `INPUT/` are copies of `test/MI2/cases/ice-exchange`'s but for the
four files the burning section changes or adds (`part-phase.txt`, `part-models.txt`, `part-products.txt`,
`part-properties.dat`).
`[GPB-Phase2]` of `input.ini` names the reaction, `particle-mechanism = a-n2` (`a-n2.yaml` of this folder: `A => N2`,
the one species of the gas), gives the material's enthalpy at 298.15 K on the gas's scale, `h0 = 1.0e6`, and the
burning model: `combustion = direct`, which GPB writes on the material line of `part-phase.txt`, and its constants
`q-comb = 1.0e6`, `d-stop = 1.0e-6`, `T-burn = 300`, which it writes as the material's row of `part-models.txt` (the
same lines as `test/MI2/cases/igloo-burn-box-direct`, its twin in IGLOO). The latent heat of the model is the material's: this
one has constant properties and no `Lv` column in its table, so `[ICE-Physics] latent-heat = 2.26e6` gives it. The
drops are at 300 K, the
burning temperature, so the first source ICE computes already burns. `part-products.txt` holds the record
`material A combustion 1 formation 1.0e6` with the term `N2 1.0`: the gas receives N2, and the heat of reaction is
`h0` less the enthalpy of N2 at 298.15 K, 1.0e6 J/kg. `test/MI2/cases/ice-burn-box-beckstead` is the same case with the other
model; `check.py` reads both.

Schedule as in `test/MI2/cases/ice-exchange`: one forward-Euler ICE step, then 200 coupling steps of one time-accurate MOSE
iteration on the frozen source.

## The oracle

With `V` the box volume, `d` the drop diameter, `c_pg = gamma R/(gamma - 1)`, `Nu = 2`, `q` = `q-comb`,
`L_v` = `latent-heat`, `c` the material's heat capacity, `h298` the enthalpy of the record:

```
B      = (c_pg (T_g - T_p) + q) / L_v                  transfer number
mdot   = pi d (k_g / c_pg) Nu ln(1 + B)                burn rate of one drop          [kg/s]
M      = n_p mdot V                                    mass to the gas                [kg/s]
F      = (rho_p (u_p - u_g)/tau + n_p mdot u_p) V      force on the gas               [N]
W      = rho_p u_p (u_p - u_g)/tau V                   drag work                      [W]
dh     = h298 - c 298.15                               the burnt mass on the gas's enthalpy scale [J/kg]
Q      = W + M (c T_p + u_p^2/2 + dh)                  energy to the gas              [W]
```

Over cycle 2, of duration `t`, the gas gains `M t` of mass, `F t` of momentum and `Q t` of energy. `q-comb` is in the
rate only: the heat the gas receives is in `dh`. The convective heat is absent from `Q`: the burning drop is held at
its temperature, and the heat the gas gives it returns with the vapour.

## Running it

```
test/MI2/run.sh ice-burn
```

runs the case, then `check.py <run-dir> <log>`; it also checks that the case stops when the products file has no
combustion record for the material, a record with the former word `explicit`, the former key `d-res`, or a key of
the oxide cap.
