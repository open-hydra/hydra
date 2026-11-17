# ice-al-freeze-box — the smoke of burning aluminium freezes: its heat of fusion reaches the gas

ICE alone (cycle 1, the gas held at 1750 K): liquid aluminium (family 1) burns direct down to d-stop at 2800 K (its
initial temperature, above T-burn 2791 K: the direct law holds a burning drop's temperature),
`2 AL(L) + 3 H2O => AL2O3(L) + 3 H2` (`al-h2o.yaml`); the alumina is the smoke (family 2), born at 2800 K, liquid
above T-melt 2327 K, cooling below T-nuc (1861.6 K), freezing with ICE's solidification model (the heat of fusion
and the solid's cp of nasa9.yaml at 2327 K, as GPB prints them when the model is off), down to the gas temperature.

`check.py`: the books of the three bodies (the gas's share the coupler reports, the aluminium, the smoke) on the
formation scale, and the smoke's energy against the database's solid AL2O3(a) at the end. `INPUT/` is what ATLAS
writes for `input.ini` (`test/mi2-checks/regen_case.sh test/ice-al-freeze-box SCRATCH GPB ICB BCB --install`), then
`python3 test/ice-al-box/make_case.py test/ice-al-freeze-box`: family 1's velocity and exact number density, family 2
empty (1e-25 in every slot). `MESH/` is `test/ice-box`'s. Gate: `test/mi2-checks/p32_fusion_check.sh`.
