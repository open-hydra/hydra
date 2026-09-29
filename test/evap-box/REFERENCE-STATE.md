# Reference state — evap-box

## Landing on hydra upstream (2026-09-29): MOSE 95aeab6, ICE 494668b, ORION a770748

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `638977860a6b1115173d30e999cdf405` (copy `test/mi2-checks/bin/hydra_MI2.E-2026-09-29`; hydra `c3686c0`
with the MI2 port, MOSE **95aeab6**, ICE **494668b**, IGLOO **e5faf0a**, ORION `a770748`, FLINT `892c29a`, OSLO `f73983f`,
FiNeR `cfc9194`), INI and inputs as shipped, **1 thread**, 6.4 s, exit 0.
The previous reference is kept as `OUTPUT.pre-landing-2026-09-29/` and
`run.pre-landing-2026-09-29.log`.

Against the previous reference 17 of 20 artefacts differ. On the scale of each variable of `gas-field.tec`: rho(2) 8.8e-8,
T 1.1e-7, rho(1) 7.4e-6, u 5.1e-5, v 6.0e-5 (the out-of-plane velocity is round-off around zero). The log differs in five
lines: MOSE's density residual in its sixth digit from iteration 650 on, and the N2 gained over the run, `9.476101E-11` →
`9.476120E-11` kg. Reference verdict (`check_reference.py`) on the new files: **PASS**, 0 checks failed.

Why it moved: the compiler orders the floating-point operations of MOSE `95aeab6` differently from those of `19aec48`; no
arithmetic changed. The builds that isolate each change, and their `-fp-model strict` twins, are in
`test/evap-box-ta/REFERENCE-STATE.md`.

## Scatter fix (2026-09-27): `part-scatter-A-sweep{1..4}.dat` replaced

IGLOO **d387f5a** (= `1734860` after a test-only amend; `fix(scatter)`: a step cut at a cell face and integrated again adds its scatter weight once; `Lib_Integration.f90` only, on `58b49bb`) → `bin/hydra_MI2` md5 `137de966f77b2503934c73ecb50e5305` (copy `test/mi2-checks/bin/hydra_MI2.scatter-fix-2026-09-27`).
1-thread A/B on this case's inputs, `6b899772` (`hydra_MI2.psat-db-2026-09-27`) → `137de966`: the 15 non-scatter artefacts and the log (minus the `# ` stamps) byte-identical; scatter markers per file: sweep 0
(untagged) 114 → 114, byte-identical; sweeps 1–4 154 → 115, 119 → 118, 148 → 121, 124 → 122; 659 → 590 in total. Before the
swap an old-binary rerun reproduced the stored `OUTPUT/` (20 files) and log exactly. Only the scatter files were swapped in; the
rest of `OUTPUT/` and `run.log` are what both binaries write. The replaced files are in `OUTPUT.pre-scatter-fix-2026-09-27/`.

Reference verdict (`check_reference.py`, its sample count bound N − 1 … N and the scatter clock now gated on both sweeps): **PASS**,
every check gated on both sweeps; last-sweep samples 30/31/31/30, every sample at least 1.6 µs after its quantum time and within
0.229 dt_q. With the replaced files swapped back in, the last sweep fails the count on streams 1 and 3 (31 and 32 samples against
at most 30 and 31) and the clock (samples up to 0.61 ms before their quantum); the check as it was before this fix passed them.

## Output keys (2026-09-27): current `OUTPUT/` + `run.log`

INI: `[IGLOO-General] out-time = on` and `print-dtime = 1.0e-3`, `[MOSE-IO] sol-diter = 267` (was 200); backup
`input.ini.pre-outputkeys`. Run: `$HYDRADIR/bin/hydra_MI2` md5 `6b899772…` (= `test/mi2-checks/bin/hydra_MI2.psat-db-2026-09-27`),
`test/mi2-checks/run_case.sh`, **1 thread**, exit 0 (wall 24.9 s on a host at load ~260). The previous baseline is kept as
`OUTPUT.pre-outputkeys-2026-09-27/` + `run.pre-outputkeys-2026-09-27.log` (next section).

The keys change what is written, not what is computed. Against the previous baseline:
- byte-identical: the five `part-scatter-A*.dat`, `part-source-sweep4.tec`, `part-euler1-sweep4.tec`, `gas-residual-history.dat`;
- `part-outloc-A*.dat` (5 files): every row identical in its 9 columns, plus a time column `t`;
- `part-trajectories-A*.dat` (5 files): every row identical in its 10 columns, plus a time column, and one row added per parcel
  at its exit mass (4 per sweep; 3 in sweeps 1 and 3, where parcel 2 already had that row);
- `gas-field.tec`, `gas-diagnostic.tec`: now written at MOSE iteration 801, the gas the last IGLOO sweep sampled (before: 800;
  267 divides 801). Byte-identical to those of a `sol-diter = 1` run, whose 20 artefacts all equal this run's;
- `run.log`: the `run_case.sh` header (exe path; `IGLOO-tree=58b49bb`, the checkout's HEAD, which is `888e785` below after a
  message-only amend), the solution-write lines (after iterations 267, 534 and 801 and again at the final postprocess; before:
  200, 400, 600, 800) and the wall/end lines. Nothing else.

Reference verdict (`check_reference.py`; the reference is described in `REFERENCE.md`): **PASS**. Every check was gated on both
sweeps except the scatter clock (C); the clock's sweep-0-only gating and a sample-count allowance of one per face crossing existed
for the scatter double count and went with the scatter fix (section above).
- Sweep 0 (the initial gas): exit points 0.020532 m (reference 0.0205322 m); exit T 302.048731 K, equal to the reference at the
  printed exit mass 9.959410e-16 kg; exit time 8.05592007e-2 s against 8.05591999e-2 s; 114 samples within 0.99 of their print
  budget; plateau within 3.0e-6 K of the closed-form balance.
- Sweep 4 (the gas at iteration 801): exit points and exit T equal to the reference to the print, exit times within 1.7e-9 s;
  samples within 1.6 of their print budget (tolerance 10); plateau 8.7e-5 K; source field equal to the injected fluxes
  (4.000000000e-6 kg/s, 2.000000000e-6 N, -63.09831790 W); Eulerian integrals within 7.2e-6 and 2.7e-6; the coupler's
  gas-gained lines equal to the field (1.261670e-8 kg H2O, 9.476101e-11 kg N2).

To check a run: `python3 check_reference.py` from this directory (1–2 min, exit 0 = PASS), or
`python3 check_reference.py RUN_DIR [--log RUN_LOG]` for a run elsewhere (`REFERENCE.md`, Use).

Byte-identity of every file is expected on a rerun at `OMP_NUM_THREADS=1`.

## Psat-column baseline (2026-09-27): now `OUTPUT.pre-outputkeys-2026-09-27/` + `run.pre-outputkeys-2026-09-27.log`

The first baseline with p_sat from the database (`Psat` column). Until the output keys it was `OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `6b8997721313bb6c66a6eaa3d091ea3c` (copy `test/mi2-checks/bin/hydra_MI2.psat-db-2026-09-27`;
IGLOO **888e785** (= `58b49bb`, pushed; message-only amend — `run.log` stamps `IGLOO-tree=888e785`), ICE **9835ab1** (pushed), MOSE unchanged), INI as shipped (`N-coupling = 5`, `MOSE-iter = 200`, `IGLOO-iter = 1`,
steady), **1 thread**, 3.2 s, exit 0.

This case had no valid reference: the stored `OUTPUT/` and `run.log` were produced on 2026-09-02 against the inputs before
their regeneration of 2026-09-18 (the log still classified the box as an axisymmetric wedge and carried no routing table or
datum reconciliation), and no binary of this month reproduces them. They are kept as history only in `OUTPUT.2026-09-02/` and
`run.2026-09-02.log`. `OUTPUT.pre-psat-db-2026-09-27/` and `run.pre-psat-db-2026-09-27.log` hold the case on the previous inputs
with the landing binary (leg a below), reproducible; `hydra_MI2.landing-2026-09-27` and `hydra_MI2.psat-table-2026-09-26`
produce it byte for byte.

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
| a | `hydra_MI2.landing-2026-09-27` (md5 `71520662…`, IGLOO `4b30e47`, ICE `3187afb`) | previous INI + INPUT | the case on its previous inputs (kept as `OUTPUT.pre-psat-db-2026-09-27/`) |
| b | `hydra_MI2.psat-db-2026-09-27` (md5 `6b899772…`, IGLOO `888e785` = `58b49bb`, ICE `9835ab1`) | previous INI + INPUT | = a, every artefact and the log (the two commits change a refusal text only) |
| d | same | new INI + previous INPUT | = b (the two keys are read by ATLAS only) |
| c | same | new INI + INPUT | this baseline |

The whole move a → c is therefore the `Psat` column: IGLOO takes `p_sat(T_p)` from the table, 14 % above the previous
Clausius-Clapeyron value where the droplets sit, so they evaporate faster, run cooler and burn out sooner.

All 20 artefacts moved. Droplets, 4 parcels, a → c:
- sweep 0 (the untagged files; the initial gas): exit (burnout) T 303.688 → 302.049 K, exit x 0.020533 → 0.020532 m, scatter rows
  118 → 114 (29.5 → 28.5 per parcel, the residence time), lowest droplet T on the scatter 300.855 → 299.192 K;
- final sweep (`-sweep4`; the gas at MOSE iteration 801): exit T 302.657/302.220/302.242/302.742 → 301.042/300.608/300.633/301.124 K,
  exit x 0.019746–0.019877 → 0.019774–0.019896 m (the parcels drift in y, to 0.0318/0.0454/0.0547/0.0684 m), scatter rows
  128 → 124, lowest droplet T on the scatter 299.454 → 297.826 K.

The steady run's closed-domain audit is not meaningful here
(its injected column reads 0: the audit accounts injection only in the time-accurate cases); it prints H2O gained 1.261670E-08
(unchanged) and N2 9.475749E-11 → 9.476101E-11. Gas field, scale-relative: H2O density `rho(1)` 5.1e-3, N2 `rho(2)` 5.6e-5, T 7.3e-5,
p 2.7e-7.
