# igloo-burn — a burning IGLOO material handed to the gas

A gate for what `hydra-MI2` hands the gas when an IGLOO material burns: the mass (`sourceMass`, as the gas species of
`INPUT/part-products.txt`), the momentum, and the energy of the burnt mass moved to the gas's enthalpy scale, the
heat of reaction inside it, against closed forms that follow from `input.ini` and the `INPUT/` files.

## The case

`a periodic variant of test/MI2/cases/igloo-evap` with two changes: the gas has one species (N2), and the material burns with IGLOO's `Beckstead`
law (`K-burn = 1.39e-6`, `n-burn = 1.8`, `X-eff = 0.5`, `T-ign = 290`, `beta-part = 0.3`, `q-comb = 1.0e6`) into that
species: `particle-mechanism = a-n2` (`a-n2.yaml` of this folder: `A => N2`) and `h0 = 1.0e6`, the material's
enthalpy at 298.15 K on the gas's scale. Four streams of 20 µm drops, 1e-6 kg/s each, enter a box of quiescent N2 at
350 K at 1 m/s and 300 K (the state of `test/MI2/cases/ice-burn`'s cloud). They ignite at once and burn out in about 5 ms
(`d0^n / (K-burn X-eff)`), well inside the box and inside `time-end = 0.05`. `part-products.txt` holds the record
`material A combustion 1 formation 1.0e6` with the term `N2 1.0`: the heat of reaction is `h0` less the enthalpy of
N2 at 298.15 K.

`INPUT/` is what ATLAS writes for `input.ini` (`GPB BCB ICB`); `MESH/` is `a periodic variant of test/MI2/cases/igloo-evap`'s.

## The oracle

Every stream burns out inside the box: what it carries in ends in the gas, whatever the law does on the way. Per
stream of mass flow `mdot`, velocity `u_inj` and temperature `T0`, with `h_c` the material's enthalpy table (on the
gas's scale: it holds `h0`) and `h298` the enthalpy of the record:

```
M  = SUM mdot                                  mass to the gas                             [kg/s]
Fx = SUM mdot u_inj                            force on the gas                            [N]
dh = h298 - h_c(298.15 K)                      the burnt mass on the gas's enthalpy scale  [J/kg] (0 here)
Q  = SUM mdot (h_c(T0) + u_inj^2/2 + dh)       energy to the gas                           [W]
```

`beta-part` is absent: the share of the heat the law leaves on the drop heats it, and returns to the gas by
convection and with the drop's mass. Over a MOSE block with these sources in force the gas energy rises by
`Q (t1 - t0)`, and the gas warms as the first law of the rigid box says.

## Running it

```
test/MI2/run.sh igloo-burn
```

runs the case, then `check.py <run-dir> <log>`; it also checks that the model given by `[IGLOO-Models] combustion`
instead of the material line writes the same bytes, and that the case stops when the products file has no
combustion record for the material.
