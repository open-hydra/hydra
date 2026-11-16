# igloo-burn-box-rx — graphite burning to CO2

`test/igloo-burn-box` with a reaction between species: drops of graphite (`material = C(gr)`, `thermo = NASA9`) burn
with IGLOO's `Beckstead` law in a gas of N2, O2 and CO2; `particle-mechanism = c-o2` names `c-o2.yaml` of this
folder (`C(gr) + O2 => CO2`). `INPUT/part-products.txt` holds the record `material C(gr) combustion 2 formation
8.31e-04` with the terms `O2 -2.664` and `CO2 3.664` (kg per kg of graphite): the gas loses the oxygen and gains the
dioxide, and the heat of reaction, 32.76 MJ/kg, is in the formation enthalpies of the species.

Closed forms and checks: `test/igloo-burn-box/check.py` and `test/igloo-burn-box/README.md`; beside them the two
species close in the coupler's mass balance in the ratio of the record, and the gas warms by the heat of reaction.

`INPUT/` is what ATLAS writes for `input.ini` and `c-o2.yaml` (`GPB BCB ICB`); `MESH/` is `test/evap-box-mom`'s.
Gate: `test/mi2-checks/p22_igloo_burn_box_check.sh`.
