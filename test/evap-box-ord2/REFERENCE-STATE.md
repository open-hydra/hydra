# Thread counts (2026-09-30)

The binary of 2026-09-30 (`test/mi2-checks/bin/hydra_MI2.G-2026-09-30`, byte-identical outputs to the reference binary)
run twice at each of 1, 2 and 8 threads on this case: every field file (`gas-field.tec`, IGLOO's Eulerian and source
fields) is byte-identical across the six runs; the parcel listings (`part-outloc-*`, `part-scatter-*`,
`part-trajectories-*`) differ only in the order of their rows (identical as sorted lists); the log differs only in the
line that prints the thread count, so every balance line is the same. The reference is a 1-thread run; any thread count
reproduces its fields.

---

# Reference state (2026-09-29) — landing on hydra upstream: MOSE 95aeab6, ICE 494668b, ORION a770748

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `638977860a6b1115173d30e999cdf405` (copy `test/mi2-checks/bin/hydra_MI2.E-2026-09-29`; hydra `c3686c0`
with the MI2 port, MOSE **95aeab6**, ICE **494668b**, IGLOO **e5faf0a**, ORION `a770748`, FLINT `892c29a`, OSLO `f73983f`,
FiNeR `cfc9194`), INI and inputs as shipped, **1 thread**, 7.4 s, exit 0.
The previous reference is kept as `OUTPUT.pre-landing-2026-09-29/` and
`run.pre-landing-2026-09-29.log`.

Against the previous reference 41 of 64 artefacts differ and the log (minus the `# ` stamps) is byte-identical. On the scale of
each variable of `gas-field.tec`: u 1.8e-11, v 1.1e-11, every other variable below 1e-14.

Why it moved: the compiler orders the floating-point operations of MOSE `95aeab6` differently from those of `19aec48`; no
arithmetic changed. The builds that isolate each change, and their `-fp-model strict` twins, are in
`test/evap-box-ta/REFERENCE-STATE.md`.

---

# Reference state (2026-09-27, later) — scatter fix: nothing moved

IGLOO **d387f5a** (= `1734860` after a test-only amend; `fix(scatter)`: a step cut at a cell face and integrated again adds its scatter weight once; `Lib_Integration.f90` only, on `58b49bb`) → `bin/hydra_MI2` md5 `137de966f77b2503934c73ecb50e5305` (copy `test/mi2-checks/bin/hydra_MI2.scatter-fix-2026-09-27`).
1-thread A/B on this case's inputs, `6b899772` (`hydra_MI2.psat-db-2026-09-27`) → `137de966`: all 64 artefacts — the 20 scatter files included — and the log (minus the `# ` stamps) byte-identical. No file was
replaced: the stored `OUTPUT/` is also the fixed binary's output.

---

# Reference state (2026-09-27) — p_sat from the database (`Psat` column)

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `6b8997721313bb6c66a6eaa3d091ea3c` (copy `test/mi2-checks/bin/hydra_MI2.psat-db-2026-09-27`;
IGLOO **888e785** (= `58b49bb`, pushed; message-only amend — `run.log` stamps `IGLOO-tree=888e785`), ICE **9835ab1** (pushed), MOSE unchanged), INI as shipped, **1 thread**, 2.9 s, exit 0.

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

41 of 64 artefacts moved: `part-outloc-A*` (20), `part-scatter-A*` (20), `part-euler1-sweep19.tec`. Droplets (4 parcels, a → c; IGLOO writes sweep 0 untagged and sweeps 1–19 as `-sweepN`): sweep 0 — burnout T 303.692 →
302.048–302.054 K, burnout x 0.020533 → 0.020532 m, scatter rows 118 → 114, lowest droplet T on the scatter 300.855 → 299.192 K;
final sweep 19 — burnout T 302.122–302.528 → 300.517–300.922 K, burnout x 0.019812–0.019836 → 0.019837–0.019860 m, scatter
rows 129 → 124 (the residence time), lowest droplet T 299.416 → 297.791 K.
Byte-identical: `gas-field.tec`, `gas-residual-history.dat`, `part-source-sweep19.tec` and the 20 trajectory files (they hold the
injection rows only). The log gains the tabulation line and nothing else: every balance number is unchanged (mass rel.err
2.212E-16, energy rel.err 2.425E-09, delta −9.492E-01 K) — the streams burn out within 1.6 mm of injection in every sweep, in the same cells
(`part-source-sweep19.tec` byte-identical), so the deposit and the gas do not see the faster rate at this resolution.

Previous reference: `OUTPUT.pre-psat-db-2026-09-27/`, `run.pre-psat-db-2026-09-27.log`.

---

# Reference state (2026-09-23) — FIRST clean baseline

`OUTPUT/` + `run.log`: `bin/hydra_MI2` linked to IGLOO **aba343c** (branch `work-in-progress` tip of
2026-09-23, after the five review-residual plans D/C1/C2/A/B) / ICE **7e1b3ac** + the
402/403/301 reader fix (that exact tree committed 2026-09-23 as ICE `4981fdb` (rebased 09-24 → `e05b005`), + `9bca917` (rebased 09-24 → `9bacd53`) line endings only), INI as shipped, **1 thread**, 4.42 s, exit 0.

| quantity | value |
|---|---|
| mass, species 1 injected = gas gained | `1.495475E-08`, rel.err **`2.212E-16`** |
| energy handed to gas (sourceEn + datum) | `-2.36203911E-01 J`, of which datum `-2.98989723E-04 J` |
| energy rel.err (bookkeeping, datum-blind) | **`2.425E-09`** |
| mass-weighted gas T | `350.000000 -> 349.050819`, delta **`-9.492E-01 K`** |
| IGLOO table datum | absolute, `hOff = -1.7113460E+07 J/kg` |
| datum offset applied by MI2 | **`-1.9992958E+04 J/kg`** |

**This case is `gas-order = 2`** (`input.ini:60`) — it is the only evap box that exercises the ord2 source
reduction, so it IS sensitive to IGLOO plan A (O28) and is the case to watch for ord2 regressions. Its sibling
`test/evap-box` and `test/evap-box-ta` are `gas-order = 1` and must stay byte-identical across ord2 changes.

## Why this is a first baseline and not a re-baseline

The previous artefacts (`OUTPUT.pre-atlas-handoff-2026-09-02/`, `run.pre-atlas-handoff-2026-09-02.log`,
`run_ord2.pre-atlas-handoff-2026-09-02.log`) were produced on 2026-09-02 against a **relative** enthalpy table:
their datum offset is `-1.7133453E+07 J/kg`, the relative-table signature. `INPUT/` was regenerated on
2026-09-18 to an **absolute** `Enthalpy_abs` table (`input.ini:104` `h0 = -15866000`), and the offset for that
table is `-1.9992958E+04`. The stored outputs therefore never corresponded to the shipped inputs and could not
be used as a comparand — they are kept as history only. The mass and energy rel.err happen to be unchanged
(`2.212E-16` / `2.425E-09`): those balances are datum-blind by construction, which is exactly why they did not
catch the mismatch.

`GAME OVER` at the end of the log is the normal terminator, not a divergence signal — it is present in every
passing MI2 run including `test/evap-box-ta`'s reference. Byte-identity only at `OMP_NUM_THREADS=1`.
