# ice-burnout-box — an ICE material that burns out

`test/ice-burn-box-beckstead` without a residue: the material's row of the models file carries no `d-stop`, and
`d-burnout = 3.0e-6`. A
drop burns on the d^n law until it is at or below 3 µm; what is left of the drops of a cell is then stored by ICE
(mass, momentum, energy) and handed whole to the gas over the next MOSE block, and the cell is empty. `MESH/` and
`INPUT/` are `test/ice-burn-box-beckstead`'s but for `INPUT/part-models.txt`.

## The case

Cycle 1 advances ICE alone (`N-coupling = 600`, `ICE-iter = 1`, `MOSE-iter = 0`): one forward-Euler step of
`dt-max = 1.0e-5` per coupling step, through the whole burn (5 ms on the law); the gas does not move. With
`[MI2-Coupling] report = true` the coupler prints, after each ICE step, what the source of that step would hand the
gas. Cycle 2 is one MOSE iteration, which prints the source of the last ICE step and hands the gas the stored
removal.

## The oracle

The cloud ends empty, so what it held is what its sources handed, whatever the law, the drag and the heat exchange
do on the way. With `M_j`, `F_j`, `Q_j` the rates printed after ICE step `j`, `dt` the step, `V` the box volume,
`m0 = rho_p V` the mass of the cloud, `h298` the enthalpy of the record of `INPUT/part-products.txt` and `R_M`,
`R_F`, `R_Q` what the removal took out, stored (`ice-removed-stored`) or handed by a MOSE block
(`ice-removed-delivered`):

```
SUM_j M_j dt + R_M = m0                                    mass      [kg]
SUM_j F_j dt + R_F = m0 u_p                                momentum  [kg m/s]
SUM_j Q_j dt + R_Q = m0 (c T_p + u_p^2/2 + h298 - c 298.15) energy on the gas's scale [J]
```

The removal: the step that first stores an amount stores, beside what it burns, all the cloud still held; it is the
first step that leaves a drop at or below the mass of a drop of `d-burnout`; after it the source is zero, the amount
removed does not change and ICE's field is empty.

## Running it

```
test/mi2-checks/p25_burnout_check.sh [scratch-root] [exe]
```

runs the case and `check.py <run-dir> <log>`, then its variants and the twin: a step twenty times longer (the
removal fires when a step leaves less than it burnt), `d-burnout` above the drops with one ICE step whose source is
frozen through 5 MOSE iterations (`check.py ... frozen`: the gas gains the burn over their time and the stored rest
of the cloud once), the same with two ICE steps (`frozen2`), drops that never burn (`cold`, in both boxes),
`test/igloo-burnout-box` against its closed forms and the energy per kg handed in the two solvers, and the refusal
of `d-stop` beside `d-burnout`.
