# Reference state (2026-09-29) — landing on hydra upstream: MOSE 95aeab6, ICE 494668b, ORION a770748

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `638977860a6b1115173d30e999cdf405` (copy `test/mi2-checks/bin/hydra_MI2.E-2026-09-29`; hydra `c3686c0`
with the MI2 port, MOSE **95aeab6**, ICE **494668b**, IGLOO **e5faf0a**, ORION `a770748`, FLINT `892c29a`, OSLO `f73983f`,
FiNeR `cfc9194`), INI and inputs as shipped, **1 thread**, 6.3 s, exit 0.
The previous reference is kept as `OUTPUT.pre-landing-2026-09-29/` and
`run.pre-landing-2026-09-29.log`.

Against the previous reference 61 of 64 artefacts differ. On the scale of each variable of `gas-field.tec`: T 5.1e-9,
rho(2) 4.1e-9, rho(1) 3.2e-7, u 1.5e-6, v 1.9e-6 (the out-of-plane velocity is round-off around zero). The log differs in one
line, the mass gate: rel.err `0.000E+00` → `2.212E-16`. Every other printed line is unchanged: energy handed to the gas
`-2.36203431E-01 J`, bookkeeping rel.err `2.790E-09`, mass-weighted gas T `350.000000 -> 349.050910`, datum offset
`-1.9992958E+04 J/kg`.

Attribution at 1 thread, one change per build, each build compared with the one before it on the four evaporation cases in
full and on `N-coupling = 3` runs of the three JPL cases (binaries `test/mi2-checks/bin/hydra_MI2.{Aprime,B,C,D,E}-2026-09-29`
and their `strict` twins):

| build | change | outputs and logs against the build before |
|---|---|---|
| A′ | the sources of the previous reference binary, built again | byte-identical |
| B | ORION `a5ca4c5` → `a770748` (the Tecplot ASCII reader rewritten) | byte-identical |
| C | MOSE `19aec48` → `95aeab6` (101 commits) | **every case moves**; with `-fp-model strict` B and C are byte-identical |
| D | hydra upstream `c3686c0`: the shared sources, FiNeR `18fa207` → `cfc9194`, FLINT `892c29a`, `lib/OSLO` | byte-identical, and so are the strict twins |
| E | ICE `9835ab1` → `494668b` (solidification), the coupler's energy slot at `nbase-1` | byte-identical, and so are the strict twins |

The move is therefore the compiler ordering the floating-point operations of the new MOSE differently: no arithmetic changed
on the paths these cases use. With a solidifying family (`solidification=on h-fus=1.09e6 cp-solid=1380` on the material line
of scratch copies of JPL-Eulerian and JPL-3solver) build E gives the gas field of the case without the tokens, to round-off
and to the ODE tolerance; the same tree with the energy read at `ncond-1` moves the gas temperature by 1.1e-2 and 6.9e-3 of
its scale.

⚠ Restart of ICE from its own output does not work at these pins: ICE `494668b` writes the variable names of its field files
without a separator (`"rho_p1""u_p1"…`), and ORION `a770748` reads that as one name (`[ICE-Parameters] newrun = false` stops
with `the initial condition holds 1 variables per cell`). The initial conditions ATLAS writes are read correctly; MOSE's
restart works. ICE's branch `merging` writes the names with a separator.

---

# Reference state (2026-09-27, later) — scatter fix: the scatter files replaced

IGLOO **d387f5a** (= `1734860` after a test-only amend; `fix(scatter)`: a step cut at a cell face and integrated again adds its scatter weight once; `Lib_Integration.f90` only, on `58b49bb`) → `bin/hydra_MI2` md5 `137de966f77b2503934c73ecb50e5305` (copy `test/mi2-checks/bin/hydra_MI2.scatter-fix-2026-09-27`).
1-thread A/B on this case's inputs, `6b899772` (`hydra_MI2.psat-db-2026-09-27`) → `137de966`: the 44 non-scatter artefacts and the log (minus the `# ` stamps) byte-identical — every balance line included; scatter
markers 2427 → 2386 in total: sweep 0 (untagged) byte-identical, sweeps 1–19 moved (sweep 1 125 → 115, sweep 12 131 → 121, the
others by 0–3; three keep their count with markers at other positions). Only the 20 scatter files were swapped in; the rest of
`OUTPUT/` and `run.log` are what both binaries write. The replaced files are in `OUTPUT.pre-scatter-fix-2026-09-27/`.

---

# Reference state (2026-09-27) — p_sat from the database (`Psat` column)

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `6b8997721313bb6c66a6eaa3d091ea3c` (copy `test/mi2-checks/bin/hydra_MI2.psat-db-2026-09-27`;
IGLOO **888e785** (= `58b49bb`, pushed; message-only amend — `run.log` stamps `IGLOO-tree=888e785`), ICE **9835ab1** (pushed), MOSE unchanged), INI as shipped, **1 thread**, 3.1 s, exit 0.

Input change (the only one): `[GPB-Phase2]` gains `psat-liquid = H2O(L)` and `psat-vapour = H2O` (the material `A` is no
species of the database, so the liquid is named; backup `input.ini.pre-psat-db`). `INPUT/` regenerated with ATLAS branch
`feat/psat-from-databases` `e666baa` (`test/mi2-checks/regen_case.sh <case> <scratch> GPB BCB ICB`, `ATLASDIR` = that
worktree; a control regeneration with the old INI reproduced all 11 files byte for byte). With the keys only
`part-properties.dat` changes: the header gains `, "Psat"` and each of its 5000 rows one `%.6e` value, nothing else
(NASA9 `H2O(L)`/`H2O`: 1.137e3 Pa at 282 K, 3.534e3 at 300 K, 4.133e4 at 350 K, 0.9854 atm at 373.15 K — x1.14 the case's
Clausius-Clapeyron line at 280–300 K, x1.05 at 350 K). Backup `INPUT.pre-psat-db-2026-09-27/`. `[IGLOO-Properties] psat = 3.0e3`
is kept: accepted and ignored. The log reports `p_sat tabulated from properties.dat (1..5000 K), psat(boiling-temperature)/Patm =
0.9854`.

Attribution, all at 1 thread (runs in scratch, `run_case.sh`):

| leg | binary | inputs | result |
|---|---|---|---|
| a | `hydra_MI2.landing-2026-09-27` (md5 `71520662…`, IGLOO `4b30e47`, ICE `3187afb`) | previous INI + INPUT | = the previous reference (`OUTPUT.pre-psat-db-2026-09-27/`), 64/64 artefacts + log |
| b | `hydra_MI2.psat-db-2026-09-27` (md5 `6b899772…`, IGLOO `888e785` = `58b49bb`, ICE `9835ab1`) | previous INI + INPUT | = a, every artefact and the log (the two commits change a refusal text only) |
| d | same | new INI + previous INPUT | = b (the two keys are read by ATLAS only) |
| c | same | new INI + INPUT | this reference |

The whole move a → c is therefore the `Psat` column: IGLOO takes `p_sat(T_p)` from the table, 14 % above the previous
Clausius-Clapeyron value where the droplets sit, so they evaporate faster, run cooler and burn out sooner.

All 64 artefacts moved. Droplets (4 parcels, a → c; IGLOO writes sweep 0 untagged and sweeps 1–19 as `-sweepN`): sweep 0 — burnout T 303.688 →
302.049 K, burnout x 0.020533 → 0.020532 m, scatter rows 118 → 114, lowest droplet T on the scatter 300.855 → 299.192 K;
final sweep 19 — burnout T 301.853–302.520 → 300.255–300.909 K, burnout x 0.019459–0.019565 → 0.019501–0.019601 m, scatter
rows 131 → 126 (the residence time), lowest droplet T 299.167 → 297.545 K. Closed-domain balances: mass 1.495472E-08 injected = gained,
rel.err **0.000E+00** (unchanged); energy handed to the gas −2.36203421E-01 → **−2.36203431E-01 J**, rel.err 2.774E-09 → **2.790E-09**;
mass-weighted gas T 350 → 349.050911 → **349.050910 K**, delta **−9.491E-01 K** (unchanged); MI2 datum offset −1.9992958E+04 J/kg and
IGLOO hOff −1.7113460E+07 J/kg unchanged. Gas field, scale-relative (max|c−a|/max|a|): H2O density `rho(1)` 6.2e-3 (the vapour cloud at
the final time — the same mass released a little earlier), N2 `rho(2)` 6.7e-5, T 8.8e-5, p 3.6e-7.

Previous reference: `OUTPUT.pre-psat-db-2026-09-27/`, `run.pre-psat-db-2026-09-27.log`.

---

# Reference state (2026-09-23)

`OUTPUT/` + `run.log`: `bin/hydra_MI2` linked to IGLOO **aba343c** (branch `work-in-progress` tip of
2026-09-23, after the five review-residual plans D/C1/C2/A/B) / ICE **7e1b3ac** + the
402/403/301 reader fix (that exact tree committed 2026-09-23 as ICE `4981fdb` (rebased 09-24 → `e05b005`), + `9bca917` (rebased 09-24 → `9bacd53`) line endings only), INI as shipped, **1 thread**, 4.66 s.

Closed-domain gates printed at the end of the run:

| quantity | value |
|---|---|
| mass, species 1 injected = gas gained | `1.495472E-08`, rel.err **`0.000E+00`** |
| energy handed to gas (sourceEn + datum) | `-2.36203421E-01 J`, of which datum `-2.98989103E-04 J` |
| energy rel.err (bookkeeping, datum-blind) | **`2.774E-09`** |
| mass-weighted gas T | `350.000000 -> 349.050911`, delta **`-9.491E-01 K`** |
| IGLOO table datum | absolute, `hOff = -1.7113460E+07 J/kg` |
| datum offset applied by MI2 | **`-1.9992958E+04 J/kg`** (route B; ~0 on an absolute table) |

**This case is `gas-order = 1`** (`input.ini:60`), so an ord2 source-reduction change (IGLOO plan A, O28)
cannot reach it — a move here after such a change is RED, not a re-baseline.

## What changed at this re-baseline, and why

Against the previous reference (`OUTPUT.igloo-cd4a3e8/`, 2026-09-18): the **run log is byte-identical** and
**63 of 64 artefacts are byte-identical**. One moved: `part-euler1-sweep19.tec`.

Cause is IGLOO **O31** (`d4b3316`, "models 2/5 deposit rho with the number rate"): `rho_p` is multiplied by a
constant **239451** in every one of the 260 non-empty cells, and `u_p`/`v_p`/`w_p`/`T_p` follow because
`finalizeEUL` divides the momentum/energy numerators by that density. `n_p` is unchanged. The old file was the
defective one; this is the fix landing, not drift. The euler field is a diagnostic — it does not feed the
coupling, which is why the gas solution and every balance number are untouched.

## Provenance of the previous reference

`OUTPUT.igloo-cd4a3e8/` + `run.igloo-cd4a3e8.log` — the 2026-09-18 run, IGLOO `cd4a3e8`. That commit is
**reachable from no branch** (it was amended away). It is pinned in IGLOO by the tag
**`reference/evap-box-ta-2026-09-18`** so `git gc` cannot drop it; without that tag the provenance of this
reference set would expire with its reflog entry. Delete the tag only when this archive is deleted.

Older archives: `OUTPUT.h0-relative-table/` (+ `run.h0-relative-table.log`) the pre-route-B relative table;
`OUTPUT.igloo-5a8eeb9-folded/` the 2026-09-02 run before the planar-slab classifier fix. Byte-identity only at
`OMP_NUM_THREADS=1`.
