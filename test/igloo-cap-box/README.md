# igloo-cap-box — an IGLOO material that burns with an oxide cap

`test/igloo-burn-box-rx` with an oxide cap. `[GPB-Phase2]` lists two materials of constant properties: the metal `A`
(1 group, injected) and its oxide `OX` (0 groups, never injected). `A` burns with IGLOO's `Beckstead` law and
`xi-cap = 0.2`: a share of the oxide forms on the drop as a cap of `OX` (`product = OX`), which the drop leaves as a
residue (`residue = on`). `particle-mechanism = a-cap` names `a-cap.yaml` of this folder, three reactions:

```
A + O2 => CO2                  A burnt into the gas                     record "combustion" of A
A + 0.3338672 O2 => OX         A taken by the cap, with oxygen          record "cap" of A (1.8894 kg of OX per kg of A)
OX + 0.6661328 O2 => CO2       a cap handed to the gas (residue = off)  record "release" of OX
```

`h0 = 3.0e7 0` puts the two materials on the gas's enthalpy scale. IGLOO hands the coupler, beside the net mass of
`A`, the oxide its caps gained and the cap oxide handed to the gas; the coupler splits the gas species between the
three records.

## The oracle

IGLOO's cap law forms `xi nu` kg of `OX` per kg of `A` burnt (`nu` = 1.8894, the ratio of alumina to aluminium), so
per kg of `A` the cap takes `xi` kg and `1 - xi` kg burns into the gas. With `g` the row of the combustion record
and `c` that of the cap record, per stream of mass flow `mdot`:

```
M     = SUM mdot (1 - xi nu)                           net mass to the gas   [kg/s]
m_s   = SUM mdot [(1 - xi) g(s) + xi c(s)]             each gas species      [kg/s]
Q     = SUM mdot [h_A(T0) + u_inj^2/2] - SUM xi nu mdot [h_OX(T_end) + u_end^2/2]     [W]
d_res = d0 (xi nu rho_A / rho_OX)^(1/3)                the residue of a drop [m]
```

`T_end` and `u_end` are the end records of the residues (`OUTPUT/part-outloc-OX*.dat`). With `residue = off` the
cap is handed to the gas at burnout through the release record, and everything a stream carries in ends in the gas
as if no cap had formed: the closed forms of `test/igloo-burn-box`.

`INPUT/` is what ATLAS writes for `input.ini` and `a-cap.yaml` (`GPB BCB ICB`); `MESH/` is `test/evap-box-mom`'s.

## Running it

```
test/mi2-checks/p26_oxide_cap_check.sh [scratch-root] [exe]
```

runs the case and `check.py <run-dir> <log>`, the variant `residue = off` against `test/igloo-burn-box/check.py`,
and three refusals: no cap record, `residue = off` without the release record, a cap record whose oxide ratio is
not the one IGLOO's cap law holds.
