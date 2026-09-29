# Reference state (2026-09-29) — landing on hydra upstream: MOSE 95aeab6, ICE 494668b, ORION a770748

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `638977860a6b1115173d30e999cdf405` (copy `test/mi2-checks/bin/hydra_MI2.E-2026-09-29`; hydra `c3686c0`
with the MI2 port, MOSE **95aeab6**, ICE **494668b**, IGLOO **e5faf0a**, ORION `a770748`, FLINT `892c29a`, OSLO `f73983f`,
FiNeR `cfc9194`), INI and inputs as shipped, **1 thread**, 206 s, exit 0.
The previous reference is kept as `OUTPUT.pre-landing-2026-09-29/` and
`run.pre-landing-2026-09-29.log`.

Against the previous reference the 6 artefacts differ at round-off (on the scale of each variable of `gas-field.tec`: T 1.6e-14,
u 1.2e-14, p 9.7e-15) and the log (minus the `# ` stamps) is byte-identical. The metrics are unchanged to every printed digit:

| | axis RMS / mean / exit M | wall RMS / mean / exit M |
|---|---|---|
| **this run** | **0.0480 / −0.0068 / 2.2254** | **0.0229 / −0.0215 / 3.3784** |
| previous (`OUTPUT.pre-landing-2026-09-29/`) | 0.0480 / −0.0068 / 2.2254 | 0.0229 / −0.0215 / 3.3784 |

Why it moved: the compiler orders the floating-point operations of MOSE `95aeab6` differently from those of `19aec48`; no
arithmetic changed. The builds that isolate each change, and their `-fp-model strict` twins, are in
`test/evap-box-ta/REFERENCE-STATE.md`.
ICE moved from `9835ab1` to `494668b` (solidification) in the same landing: byte-identical on this case.

---

# Reference state (2026-09-26, later) — air transport table (ATLAS fix), axis symmetry, JAXA3

`OUTPUT/` + `run.log`: the same binary as the section below (`test/mi2-checks/bin/hydra_MI2.rebaseline-2026-09-26`,
md5 `45477f77a0fda11b97adb5b6ad2e561d`) and the same INI; the only input that changed is `INPUT/gas-transport.dat`,
regenerated with ATLAS branch `fix/mixture-transport` (`4942252`: the mixture zone is evaluated at the mixture
composition; before, it held pure O2, the last species of the per-species loop: μ +11–13 %, k +8 % at 300 K).
μ(300 K) = 1.8607e-5 Pa s now (air, Incropera Table A.4: 1.846e-5). Backup of the previous inputs: `INPUT.pre-transport-fix-2026-09-26/`.

| JPL-Eulerian | axis RMS / mean / exit M | wall RMS / exit M |
|---|---|---|
| **this run (air transport)** | **0.0480 / −0.0068 / 2.2254** | **0.0229 / 3.3784** |
| previous reference, pure-O2 transport (`OUTPUT.o2-transport-2026-09-26/`, `run.o2-transport-2026-09-26.log`) | 0.0479 / −0.0184 / 2.2109 | 0.0256 / 3.3746 |
| Chang MPL (Chang, Hourng & Chien 1996 Fig. 11) | — / — / 2.112 | — / 3.357 |

Attribution: same binary, same INI, one input file changed, so the whole move is the transport table; direction as
predicted (a less viscous gas relaxes the particles more slowly: less drag coupling, higher Mach). ICE density residual at iteration 1500: 3.12e-4 (short run, still falling).
Still different from the reference: the gas model itself (γ 1.4, cp 1070, Sutherland μ in the paper; NASA9 air
here) — the V1/V2 experiments in `src/lib/ICE/plan-bucket/wip-M/jpl-runs.md`.

---

# Reference state (2026-09-26) — FIRST recorded baseline (ICE phase mirrored at the axis, JAXA3, emissivity 0)

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `45477f77a0fda11b97adb5b6ad2e561d` (relinked 2026-09-26 17:26: IGLOO `fb75a60`
(linked, unused here), ICE **f5efc9c** (= upstream `origin/main` `53ff29f`), MOSE `d246847`), INI as shipped
(`N-coupling = 3`, `MOSE-iter = 1000`, `ICE-iter = 500`, `heat-transfer = JAXA3`, `drag = Clift-Gauvin`,
`[ICE-Physics] emissivity = 0`, `[ax] part-type = symmetry` → 300), **1 thread**, 133 s. ICE density residual at its last
iteration (1500): 3.24e-4, still falling — this is the SHORT run; convergence (coupler plan §7.3) is a separate item.
`python3 ../JPL-Lagrangian-20micron/jpl_compare.py OUTPUT/gas-field.tec ../JPL-Lagrangian-20micron/JPL_chang.dat OUTPUT/gas-field.tec`:

| | axis RMS / mean / exit M | wall RMS / exit M |
|---|---|---|
| **this run** | **0.0479 / −0.0184 / 2.2109** | **0.0256 / 3.3746** |
| R3 (2026-09-25, same INI but the axis still 200 for ICE's phase; `jpl-runs.md`) | 0.0458 / — / 2.2100 | — |
| Chang MPL | — / — / 2.112 | — / 3.357 |

## Provenance and attribution

- Before today this case had no reference record. `OUTPUT.ice7e1b3ac-kd-2026-09-02/` + `run_ice.ice7e1b3ac-kd-2026-09-02.log`
  hold the 2026-09-02 LONG run (3000 MOSE + 1500 ICE iterations per the coupler plan, Kavanau-Drake, `emiss = 1`, ICE
  `7e1b3ac`, axis 200) and `run_P0-ice7e1b3ac.log` the 2026-09-15 P0 run; neither is comparable to the present INI, so
  this is a first baseline, not a re-baseline.
- R3 → this run (+0.0021 axis RMS, exit M 2.2100 → 2.2109) = the ICE-phase axis mirror (`part-type = symmetry`, 200 → 300,
  regenerated 2026-09-25 13:18 after R3; ICE plan H decision b). The ICE library commits of 2026-09-25 (plans H, I, E,
  M P1–P7) are inert on this case up to codegen (`-fp-model strict` bit-identity, `wip-M/jpl/kd-inertness-md5.txt`).
- Still in the inputs: the ATLAS transport-table bug (pure-O2 μ) and the gas-model difference vs the reference
  (`jpl-runs.md` V1/V2). The 1-thread run is byte-reproducible; at 4 threads MOSE/ICE reruns should agree to the printed
  digits (no IGLOO ATOMIC deposition in this case — unmeasured today).

**`INPUT/` is generated** by `regen_case.sh JPL-Eulerian-20micron <scratch> GPB BCB ICB` with ATLAS `main` `76cd8a3`
(verified byte-identical on 2026-09-26; `input.ini.pre-axis` is the pre-mirror INI).
