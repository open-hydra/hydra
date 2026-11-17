# ice-al-box — liquid aluminium burns with steam, the alumina a smoke family

ICE, Beckstead down to d-stop; the alumina is the smoke, born in family 2 of ICE's phase. One reaction, `2 AL(L) + 3 H2O => AL2O3(L) + 3 H2` (`al-h2o.yaml`); the gas is N2 H2O H2.

Closed forms and checks: `test/mi2-checks/p30_al_smoke.py`. `INPUT/` is what ATLAS writes for `input.ini`
(`GPB ICB BCB`), then over ICB's field family 1's velocity and exact number density (in `test/ice-box/make_case.py`'s format) and family 2 empty (1e-25 in every slot: ICB writes zeros for krho 0, on which ICE's drag is NaN); `MESH/` is `test/ice-box`'s. Gate: `test/mi2-checks/p30_al_smoke_check.sh`.
