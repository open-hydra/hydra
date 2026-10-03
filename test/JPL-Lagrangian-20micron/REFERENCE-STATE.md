# Reference state (2026-09-29) — landing on hydra upstream: MOSE 95aeab6, ICE 494668b, ORION a770748

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `638977860a6b1115173d30e999cdf405` (copy `test/mi2-checks/bin/hydra_MI2.E-2026-09-29`; hydra `c3686c0`
with the MI2 port, MOSE **95aeab6**, ICE **494668b**, IGLOO **e5faf0a**, ORION `a770748`, FLINT `892c29a`, OSLO `f73983f`,
FiNeR `cfc9194`), INI and inputs as shipped, **1 thread**, 2482 s on a loaded host, exit 0 (the `# exit` line
of `run.log` is missing: the harness script was edited while the run was on; the solver printed `GAME OVER` and
`/usr/bin/time` reported no failure). The previous reference is kept as `OUTPUT.pre-landing-2026-09-29/` and
`run.pre-landing-2026-09-29.log`.

Against the previous reference all 155 artefacts differ. On the scale of each variable of `gas-field.tec`: v 8.4e-9, u 4.1e-9,
T 3.1e-9. The metrics are unchanged to every printed digit, and so are the run witnesses (50 × "All particles out of domain",
no "no net progress", no "marking gone"):

| | axis RMS / mean / exit M | wall RMS / mean / exit M |
|---|---|---|
| **this run** | **0.0498 / −0.0479 / 2.0575** | **0.0237 / −0.0217 / 3.3654** |
| previous (`OUTPUT.pre-landing-2026-09-29/`) | 0.0498 / −0.0479 / 2.0575 | 0.0237 / −0.0217 / 3.3654 |

Why it moved: the compiler orders the floating-point operations of MOSE `95aeab6` differently from those of `19aec48`; no
arithmetic changed. The builds that isolate each change, and their `-fp-model strict` twins, are in
`test/evap-box-ta/REFERENCE-STATE.md`.

---

# Reference state (2026-09-27) — a 1-thread reference; scatter fix

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `137de966f77b2503934c73ecb50e5305` (copy `test/mi2-checks/bin/hydra_MI2.scatter-fix-2026-09-27`;
IGLOO **d387f5a** (= `1734860` after a test-only amend; `run.log` stamps `IGLOO-tree=d387f5a`), ICE **9835ab1**, MOSE unchanged),
INI and inputs as shipped, **1 thread**, 1559 s, exit 0. **The reference is now a 1-thread run, byte-reproducible at
`OMP_NUM_THREADS=1`**; the previous one — 4 threads, `bin/hydra_MI2` md5 `45477f77…` — is kept as `OUTPUT.pre-scatter-fix-2026-09-27/`
and `run.pre-scatter-fix-2026-09-27.log`.

The scatter fix (IGLOO: a step cut at a cell face and integrated again adds its scatter weight once). 1-thread A/B on these
inputs, `6b899772` (`hydra_MI2.psat-db-2026-09-27`) → `137de966`: the 105 non-scatter artefacts and the log (minus the `# `
stamps) byte-identical; all 50 scatter files moved, 1 096 850 → 896 915 markers in total (−18.2 %; sweep 1 21 943 → 17 931). The fix moves nothing else,
so the gas solution of this reference is `6b899772`'s at 1 thread.

Against the stored 4-thread reference the metrics agree to every printed digit (`jpl_compare.py` / `compare_3solver.py` on both
`OUTPUT/` directories):

| | axis RMS / mean / exit M | wall RMS / mean / exit M |
|---|---|---|
| **this run (1 thread)** | **0.0498 / −0.0479 / 2.0575** | **0.0237 / −0.0217 / 3.3654** |
| previous (4 threads, `OUTPUT.pre-scatter-fix-2026-09-27/`) | 0.0498 / −0.0479 / 2.0575 | 0.0237 / −0.0217 / 3.3654 |


---

# Reference state (2026-09-26, later) — air transport table (ATLAS fix), axis symmetry, JAXA3

`OUTPUT/` + `run.log`: the same binary as the section below (`test/mi2-checks/bin/hydra_MI2.rebaseline-2026-09-26`,
md5 `45477f77a0fda11b97adb5b6ad2e561d`) and the same INI; the only input that changed is `INPUT/gas-transport.dat`,
regenerated with ATLAS branch `fix/mixture-transport` (`4942252`: the mixture zone is evaluated at the mixture
composition; before, it held pure O2, the last species of the per-species loop: μ +11–13 %, k +8 % at 300 K).
μ(300 K) = 1.8607e-5 Pa s now (air, Incropera Table A.4: 1.846e-5). Backup of the previous inputs: `INPUT.pre-transport-fix-2026-09-26/`.

| JPL-Lagrangian | axis RMS / mean / exit M | wall RMS / exit M |
|---|---|---|
| **this run (air transport)** | **0.0498 / −0.0479 / 2.0575** | **0.0237 / 3.3654** |
| previous reference, pure-O2 transport (`OUTPUT.o2-transport-2026-09-26/`, `run.o2-transport-2026-09-26.log`) | 0.0638 / −0.0618 / 2.0393 | 0.0295 / 3.3552 |
| Chang MPL (Chang, Hourng & Chien 1996 Fig. 11) | — / — / 2.112 | — / 3.357 |

Attribution: same binary, same INI, one input file changed, so the whole move is the transport table; direction as
predicted (a less viscous gas relaxes the particles more slowly: less drag coupling, higher Mach). rho 0.955–1.021, p 0.931–1.006, T 0.974–1.000; `mil` 0.887–0.910 (the viscosity fell by 9–11 %).
Still different from the reference: the gas model itself (γ 1.4, cp 1070, Sutherland μ in the paper; NASA9 air
here) — the V1/V2 experiments in `src/lib/ICE/plan-bucket/wip-M/jpl-runs.md`.

---

# Reference state (2026-09-26) — axis symmetry for the Lagrangian phase, JAXA3, ATLAS main

`OUTPUT/` + `run.log`: `bin/hydra_MI2` md5 `45477f77a0fda11b97adb5b6ad2e561d` (relinked 2026-09-26 17:26 from the P7 tree:
IGLOO **fb75a60** (`work-in-progress`), ICE **f5efc9c** (linked, unused here), MOSE `d246847`), INI as shipped
(`IGLOO-iter = 2`, `[IGLOO-General] phase = part`, `heat = JAXA3`, `drag = Chang`), 4 threads, 731 s.
`python3 jpl_compare.py OUTPUT/gas-field.tec JPL_chang.dat OUTPUT/gas-field.tec` (axis) and
`… OUTPUT/gas-field.tec "Chang MPL Wall" 39` (wall):

| | axis RMS / mean / exit M | wall RMS / exit M |
|---|---|---|
| **this run (Lagrangian phase mirrored at the axis, `part-type = symmetry` → 300)** | **0.0638 / −0.0618 / 2.0393** | **0.0295 / 3.3552** |
| control, same binary and inputs but `part-type = outlet` → 400 (`OUTPUT.axis-outlet-2026-09-26/`, `run.axis-outlet-2026-09-26.log`) | 0.0523 / −0.0508 / 2.0629 | 0.0294 / 3.3553 |
| previous reference 2026-09-18 (`OUTPUT.igloo-d32a017-kd-axis-outlet-2026-09-18/`, `run_ghost.igloo-d32a017-kd-axis-outlet-2026-09-18.log`: IGLOO d32a017, Kavanau-Drake, axis outlet) | 0.0391 / −0.0367 / 2.0972 | 0.0261 / 3.3633 |
| Chang MPL (Chang, Hourng & Chien 1996 Fig. 11) | — / — / 2.112 | — / 3.357 |

Run witnesses: 635 parcels every sweep, 50 × "All particles out of domain", zero "no net progress" / "marking gone" —
the mirrored parcels leave through the exit plane, none is lost at the axis.

## What moved since 2026-09-18, and why (attribution)

1. **Heat law** Kavanau-Drake → `JAXA3` (user, 2026-09-25: the reference's Carlson-Hoglund law). The control run
   reproduces the 2026-09-25 measurement R3 (`src/lib/ICE/plan-bucket/wip-M/jpl-runs.md`: 2.0629 / 0.0523, wall
   3.3553 / 0.0294) to every printed digit, so the relinked binary and the library commits since (IGLOO unchanged at
   `fb75a60`; ICE not used by this case) are inert here: the whole 09-18 → control move is the heat law.
2. **Axis mirror for the Lagrangian phase** (this re-baseline; `plan-bucket/implemented/axis-bc-coherence-plan.md`): Chang, Hourng &
   Chien 1996 §4 and Chang 1980 impose zero normal velocity and symmetric particle variables on the centreline, so the
   outlet (parcels removed at the axis) was wrong. Control → this run: axis RMS 0.0523 → 0.0638, mean −0.0508 → −0.0618,
   exit M 2.0629 → 2.0393; wall unchanged (0.0294 → 0.0295). The parcels that used to leave now stay on the centreline,
   load it more and take more momentum and heat from the gas. The comparison with the reference gets WORSE on the axis
   with the current gas: the remaining bias is the gas model (the reference's γ = 1.4, cp = 1070, Sutherland μ; V2 in
   `jpl-runs.md` brings the axis RMS to 0.033 with Kavanau-Drake) and the ATLAS transport-table bug (pure-O2 μ, +11–13 %,
   `utils/ATLAS/src/GPB/ideal_gas/transport.py:91`), both still in this run's inputs.
3. Per-variable decomposition of `gas-field.tec`, this run vs control (`tec_block_diff.py`): every cell of rho, u, p, T
   differs by a small factor (rho 0.989–1.0004, p 0.9988–1.0011, T 0.9999–1.011); v and w change sign near the axis
   (ratios −9 … 1.5) — the signature of a centreline source change, not of a field-wide term.

**`INPUT/` is generated** by `test/mi2-checks/regen_case.sh JPL-Lagrangian-20micron <scratch> GPB BCB ICB --install` with
ATLAS `main` `76cd8a3` (the handoff branch is merged; no worktree needed): compared with the 2026-09-18 inputs only the
190 face-3 records of `part-bc.txt` changed (400 → 300, no payload); backups `INPUT.pre-axis-symmetry-2026-09-26/` (=
`INPUT.pre-regen/`), `input.ini.pre-axis-symmetry`. The trap wall (face 4 = 301) and the prefixed inlet keys are unchanged.
Same-binary 4-thread reruns differ by 2e-9..8e-9 of each variable's peak; byte-identity only at `OMP_NUM_THREADS=1`.

Older archives: `OUTPUT.igloo-5a8eeb9/`, `OUTPUT.igloo-a66e02f/` (0.0260 / 2.1115), `OUTPUT.pre-trapwall-axisreflect/`,
`OUTPUT.axisfix-dfce11c/`, `OUTPUT.dualfix-noclip/`, `OUTPUT.prefix-a09cf31/` (see project memory for each).

⚠ **Backup trap of `regen_case.sh --install`:** it creates `INPUT.pre-regen/` only when that folder is absent, and for
this case the folder now exists (created by the 2026-09-26 axis regeneration). A further `--install` overwrites `INPUT/`
with NO fresh backup: copy `INPUT/` aside first.
