# Reference state (2026-09-29) — landing on hydra upstream: MOSE 95aeab6, ICE 494668b, ORION a770748

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `638977860a6b1115173d30e999cdf405` (copy `test/mi2-checks/bin/hydra_MI2.E-2026-09-29`; hydra `c3686c0`
with the MI2 port, MOSE **95aeab6**, ICE **494668b**, IGLOO **e5faf0a**, ORION `a770748`, FLINT `892c29a`, OSLO `f73983f`,
FiNeR `cfc9194`), INI and inputs as shipped, **1 thread**, 3055 s on a loaded host, exit 0 (the `# exit` line
of `run.log` is missing: the harness script was edited while the run was on; the solver printed `GAME OVER` and
`/usr/bin/time` reported no failure). The previous reference is kept as `OUTPUT.pre-landing-2026-09-29/` and
`run.pre-landing-2026-09-29.log`.

Against the previous reference all 158 artefacts differ. On the scale of each variable: `gas-field.tec` v 5.4e-9, T 1.7e-9,
u 1.6e-9; `partE-field.tec` rho_p1 1.6e-9, u_p1 2.2e-10, T_p1 4.2e-11. The metrics are unchanged to every printed digit, and so
are the run witnesses (50 × "All particles out of domain", no "no net progress", no "marking gone"):

| | axis RMS / mean / exit M | wall RMS / mean / exit M |
|---|---|---|
| **this run** | **0.0484 / −0.0470 / 2.0801** | **0.0263 / −0.0244 / 3.3607** |
| previous (`OUTPUT.pre-landing-2026-09-29/`) | 0.0484 / −0.0470 / 2.0801 | 0.0263 / −0.0244 / 3.3607 |

Inlet-column condensed flux ICE 47.021 / IGLOO 47.114, domain-sum `rho_p` ICE 7017.806 / IGLOO 7214.467 (ratio 1.028): equal in both.

Why it moved: the compiler orders the floating-point operations of MOSE `95aeab6` differently from those of `19aec48`; no
arithmetic changed. The builds that isolate each change, and their `-fp-model strict` twins, are in
`test/evap-box-ta/REFERENCE-STATE.md`.
ICE moved from `9835ab1` to `494668b` (solidification) in the same landing: byte-identical on the `N-coupling = 3` run of
this case.

---

# Reference state (2026-09-27) — a 1-thread reference; scatter fix

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `137de966f77b2503934c73ecb50e5305` (copy `test/mi2-checks/bin/hydra_MI2.scatter-fix-2026-09-27`;
IGLOO **d387f5a** (= `1734860` after a test-only amend; `run.log` stamps `IGLOO-tree=d387f5a`), ICE **9835ab1**, MOSE unchanged),
INI and inputs as shipped, **1 thread**, 1807 s, exit 0. **The reference is now a 1-thread run, byte-reproducible at
`OMP_NUM_THREADS=1`**; the previous one — 4 threads, `bin/hydra_MI2` md5 `45477f77…` — is kept as `OUTPUT.pre-scatter-fix-2026-09-27/`
and `run.pre-scatter-fix-2026-09-27.log`.

The scatter fix (IGLOO: a step cut at a cell face and integrated again adds its scatter weight once). 1-thread A/B on these
inputs, `6b899772` (`hydra_MI2.psat-db-2026-09-27`) → `137de966`: the 108 non-scatter artefacts and the log (minus the `# `
stamps) byte-identical; all 50 scatter files moved, 1 098 628 → 898 945 markers in total (−18.2 %; sweep 1 21 989 → 17 993 in `partL-scatter-A-sweep1.dat`). The fix moves nothing else,
so the gas solution of this reference is `6b899772`'s at 1 thread.

Against the stored 4-thread reference the metrics agree to every printed digit (`jpl_compare.py` / `compare_3solver.py` on both
`OUTPUT/` directories):

| | axis RMS / mean / exit M | wall RMS / mean / exit M |
|---|---|---|
| **this run (1 thread)** | **0.0484 / −0.0470 / 2.0801** | **0.0263 / −0.0244 / 3.3607** |
| previous (4 threads, `OUTPUT.pre-scatter-fix-2026-09-27/`) | 0.0484 / −0.0470 / 2.0801 | 0.0263 / −0.0244 / 3.3607 |

Inlet-column condensed flux ICE 47.021 / IGLOO 47.114, domain-sum `rho_p` ICE 7017.806 / IGLOO 7214.467 (ratio 1.028): equal in both.

---

# Reference state (2026-09-26, later) — air transport table (ATLAS fix), axis symmetry, JAXA3

`OUTPUT/` + `run.log`: the same binary as the section below (`test/mi2-checks/bin/hydra_MI2.rebaseline-2026-09-26`,
md5 `45477f77a0fda11b97adb5b6ad2e561d`) and the same INI; the only input that changed is `INPUT/gas-transport.dat`,
regenerated with ATLAS branch `fix/mixture-transport` (`4942252`: the mixture zone is evaluated at the mixture
composition; before, it held pure O2, the last species of the per-species loop: μ +11–13 %, k +8 % at 300 K).
μ(300 K) = 1.8607e-5 Pa s now (air, Incropera Table A.4: 1.846e-5). Backup of the previous inputs: `INPUT.pre-transport-fix-2026-09-26/`.

| JPL-3solver | axis RMS / mean / exit M | wall RMS / exit M |
|---|---|---|
| **this run (air transport)** | **0.0484 / −0.0470 / 2.0801** | **0.0263 / 3.3607** |
| previous reference, pure-O2 transport (`OUTPUT.o2-transport-2026-09-26/`, `run.o2-transport-2026-09-26.log`) | 0.0622 / −0.0610 / 2.0597 | 0.0306 / 3.3527 |
| Chang MPL (Chang, Hourng & Chien 1996 Fig. 11) | — / — / 2.112 | — / 3.357 |

Attribution: same binary, same INI, one input file changed, so the whole move is the transport table; direction as
predicted (a less viscous gas relaxes the particles more slowly: less drag coupling, higher Mach). inlet-column condensed flux ICE 47.021, IGLOO 47.114 vs 47.012 prescribed per phase (the prescription follows the gas ρu, which moved with the viscosity); domain-sum `rho_p` IGLOO/ICE 1.028.
Still different from the reference: the gas model itself (γ 1.4, cp 1070, Sutherland μ in the paper; NASA9 air
here) — the V1/V2 experiments in `src/lib/ICE/plan-bucket/wip-M/jpl-runs.md`.

---

# Reference state (2026-09-26) — axis symmetry for both condensed phases, JAXA3, ATLAS main

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `45477f77a0fda11b97adb5b6ad2e561d` (relinked 2026-09-26 17:26 from the P7 tree:
IGLOO **fb75a60**, ICE **f5efc9c** (= upstream `origin/main` `53ff29f`), MOSE `d246847`), INI as shipped (50 × 500 MOSE,
IGLOO mode 2, 250 ICE; `heat-transfer = JAXA3` / `heat = JAXA3`, `[ICE-Physics] emissivity = 0`), 4 threads, 854 s.
`python3 compare_3solver.py <run-dir>` (always pass the directory explicitly):

| | axis RMS / mean / exit M | wall RMS / exit M |
|---|---|---|
| **this run (both phases mirrored at the axis)** | **0.0622 / −0.0610 / 2.0597** | **0.0306 / 3.3527** |
| control, same binary and inputs but `partL-type = outlet` (`OUTPUT.partL-outlet-2026-09-26/`, `run.partL-outlet-2026-09-26.log`) | 0.0569 / −0.0555 / 2.0719 | 0.0306 / 3.3528 |
| previous reference 2026-09-18 (`OUTPUT.igloo-d32a017-kd-partL-outlet-2026-09-18/`, IGLOO d32a017, Kavanau-Drake, ICE 7e1b3ac, partE axis 200) | 0.0434 / — / 2.1042 | 0.0276 / 3.3597 |
| Chang MPL (Chang, Hourng & Chien 1996 Fig. 11) | — / — / 2.112 | — / 3.357 |

Inlet-column condensed mass flux: ICE 46.957, IGLOO 47.049 vs 46.948 prescribed per phase (`0.21429·ρu`, half load each).
Exit-plane `rho_p` every 5 cells from the axis: ICE 0.3933 0.4738 0.8217 0.2535 0.0058 …, IGLOO 0.4809 0.5472 0.8340 0.2948
0.0007 … (control: IGLOO 0.4623 on the axis — the mirror keeps the parcels that used to leave); domain-sum `rho_p`
ICE 7089.9, IGLOO 7294.1, ratio 1.029. Run witnesses: 635 parcels every sweep, 50 × "All particles out of domain", zero
"no net progress" / "marking gone" (same as the 09-18 log).

## What moved since 2026-09-18, and why (attribution)

1. **Heat law** Kavanau-Drake → `JAXA3` (user, 2026-09-25, the reference's Carlson-Hoglund law): measured that day on the
   C3 coupler at 4 threads (`src/lib/ICE/plan-bucket/wip-M/jpl-runs.md` R3: 0.0544 / 2.0738, wall 0.0306 / 3.3528).
2. **ICE's phase mirrored at the axis** (`partE-type = symmetry`, 200 → 300, regenerated 2026-09-25 after R3; ICE plan H
   decision b): R3 → today's control = +0.0025 axis RMS, exit M 2.0738 → 2.0719, wall unchanged. The ICE library changes of
   2026-09-25 (plans H, I, E, M P1–P7) were proven inert on this case up to codegen (`-fp-model strict` bit-identity,
   `wip-M/jpl/kd-inertness-md5.txt`), and the JPL-Lagrangian control reproduces R3 to every printed digit, so this step
   carries the whole R3 → control move.
3. **IGLOO's phase mirrored at the axis** (`partL-type = symmetry`, 400 → 300, this re-baseline; decision in
   `plan-bucket/implemented/axis-bc-coherence-plan.md`: the reference imposes zero normal velocity and symmetric particle variables on
   the centreline, so an outlet there was wrong): control → this run = axis RMS 0.0569 → 0.0622, exit M 2.0719 → 2.0597,
   wall unchanged, IGLOO axis-row `rho_p` 0.462 → 0.481. More loading on the centreline, more drag and heat to the gas.
4. Not in this run: the ATLAS transport table still carries the pure-O2 viscosity of the `mixture` + `cantera` bug
   (`utils/ATLAS/src/GPB/ideal_gas/transport.py:91`, μ +11–13 % vs Sutherland); the reference gas itself differs (γ = 1.4,
   cp = 1070, Sutherland μ, Pr 0.72 — the V2 experiment in `jpl-runs.md`). Both remain the known causes of the axis bias
   (mean −0.061 here vs −0.029 with the reference gas on 09-25).

**`INPUT/` is generated** by `test/mi2-checks/regen_case.sh JPL-3solver-20micron <scratch> GPB BCB ICB --install` with ATLAS
`main` `76cd8a3` (the handoff branch is merged): compared with the 2026-09-25 inputs only the 190 face-3 records of
`partL-bc.txt` changed (400 → 300, no payload); backups `INPUT.pre-axis-symmetry-2026-09-26/`, `input.ini.pre-axis-symmetry`.
Same-binary 4-thread reruns differ by ~1e-9 of peak (IGLOO ATOMIC deposition); byte-identity only at `OMP_NUM_THREADS=1`.
Older archives: `OUTPUT.igloo-5a8eeb9/` (2026-09-16).

⚠ **Backup trap of `regen_case.sh --install`:** it creates `INPUT.pre-regen/` only when that folder is absent, and for
this case the folder now exists (created by the 2026-09-26 axis regeneration). A further `--install` overwrites `INPUT/`
with NO fresh backup: copy `INPUT/` aside first.
