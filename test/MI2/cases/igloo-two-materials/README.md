# igloo-two-materials — two injected IGLOO materials, each with its own gas species

`test/MI2/cases/igloo-evap`'s closed box (every face a symmetry plane, time accurate, 20 coupling steps of 50 MOSE
iterations) with two IGLOO materials, both injected on the same four DB streams of `[IGLOO-BC]` (20 µm, 1e-6 kg/s
per stream and material, 320 K):

- `water` (constant cp 4184 and rho 997, `h0` the formation enthalpy of liquid water, `p_sat`, `Lv`, `Mv` and the
  vapour cp from the pair `H2O(L)`/`H2O`) evaporates with CEM and Langmuir-Knudsen into `H2O`;
- `C(gr)` (NASA9, rho 2200) burns with IGLOO's `Beckstead` law, the constants of `test/MI2/cases/igloo-burn-box-rx`, to `CO2`
  with the gas's `O2` (`particle-mechanism = c-o2`, the file `c-o2.yaml` of this folder).

The gas is air (N2 0.77, O2 0.23) at 1 atm and 350 K. `INPUT/part-products.txt` holds one record per material:
`water evaporation latent` with `H2O 1` and `C(gr) combustion formation` with `O2 -2.664`, `CO2 3.664`, so each
gas species the drops change belongs to one material. The graphite burns out in about 5 ms; the water drops would
take about 80 ms, and `time-end = 0.02` stops them with about half their mass, so the two materials hand the gas
clearly different amounts (7.65e-9 and 1.45e-8 kg over the run): a material index swapped anywhere between IGLOO's
per-material rates and the gas species moves every species total by about half.

`[IGLOO-Models] evaporation = CEM` applies to both materials; IGLOO disables it for the burning one (its warning
"combustion set with evaporation also configured for material 2" is expected). Once a material evaporates IGLOO asks
`Tboil`, `Le`, `Yinf`, `Mv`, `Lv` and `cpv` of every material: `[IGLOO-Properties]` gives zeros where a column of
`part-properties.dat` gives the value, and the graphite's sublimation point as its unused `Tboil`.

`INPUT/` is what ATLAS writes for `input.ini` and `c-o2.yaml` (`GPB BCB ICB`, gfortran build, two regenerations
identical); `MESH/` is `test/MI2/cases/igloo-evap`'s. Gate: `test/MI2/run.sh igloo-two-materials`.
