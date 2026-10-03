# vortex-cloud — ICE, IGLOO and the exact solution on one case

ICE's verification case I (`src/lib/ICE/test/verification/I-vortex-cloud/run.py`): a particle
cloud released at the local gas velocity into a frozen solid-body vortex, `u_g = i Ω z`. With Stokes
drag and a uniform diameter every particle moves as `z(t) = C(t) z0`, so the cloud rotates by
`arg C` (lagging `Ω t`) and dilates by `|C|`, and the density is `ρ_p0(z/C)/|C|²` — exact for
both an Eulerian and a Lagrangian solver. `compare.py` solves it with both, at St = 0.01, 0.1, 1
and 10 over a quarter turn of the gas, and puts each beside the closed form:

- **ICE** (Eulerian MK moments): `run.py`'s case, built by the driver in its own work dir with
  `[ICE-IO] sol-overwrite = false` and `sol-dtime = QUARTER/8`, so one run leaves seven numbered
  intermediate solutions beside the final one; `run.py`'s `measure()` and `l1_error()` on each,
  against the closed form at that file's own `SOLUTIONTIME` (the final one is ICE's first step past
  the quarter turn, 1.571 s).
- **IGLOO** (Lagrangian parcels): the four `tests/standard/vortex-cloud-st*` cases, regenerated with
  IGLOO's `tests/tools/make_vortex_case.py` and run at `OMP_NUM_THREADS=1`: 1672 parcels stopped at
  the quarter turn (`time-end`), every trajectory row time-stamped (`out-time`). The cloud's own
  centroid and velocity slope come from `snapshot-A.dat`; the same parcels deposited on ICE's
  96 × 96 grid (triangular-shaped-cloud kernel) go through ICE's own `measure()` and `l1_error()`.
  The deposit of the exactly mapped lattice gives the kernel's own L1, printed beside it.

It is a report, not a gate: ICE's gate is its `VortexCloud` ctest, IGLOO's the four
`vortex-cloud-st*` ctest cases. IGLOO's suite also registers `vortex-cloud-3way`, which runs this
driver with `--check` against its own tree's binary and skips (exit 77) where hydra, the driver or
an ICE binary is absent.

## Running it

```bash
python3 test/vortex-cloud/compare.py                 # table and figures in test/vortex-cloud/work/
python3 test/vortex-cloud/compare.py --check         # also exit 1 when a value leaves its tolerance
python3 test/vortex-cloud/compare.py --reuse         # re-form table and figures from the runs in work/
python3 test/vortex-cloud/compare.py --reference     # write them into reference/ instead
python3 test/vortex-cloud/compare.py --no-cache      # run ICE even when a cached run matches
```

Needs `HYDRADIR` and Python 3 (matplotlib only for the figures; the table is always written). A
first run takes about 3 minutes, up to 8 on a busy host: four ICE runs (≈ 45 s each, one thread)
and four IGLOO runs. Later runs take about a minute, because ICE's runs are cached:

- ICE's output is a function of its binary and of the case the driver writes, so every run is kept in
  `test/vortex-cloud/.cache/<key>/` (untracked, listed in this directory's `.gitignore` with
  `work/`; about 11 MB for the four runs, gzip-compressed). The key is the md5 of the ICE binary,
  of every file the run reads (`input.ini` and `INPUT/`) and of the case parameters, so a rebuilt
  ICE or a changed case misses and runs again.
- A hit copies the stored solutions back and checks each against its recorded md5; an entry that
  is incomplete or fails that check is removed and ICE runs again, storing a fresh one (written in
  a temporary directory and renamed once complete). IGLOO is never cached: it is what the IGLOO
  suite tests.
- `--no-cache` runs ICE regardless and stores nothing; `rm -rf test/vortex-cloud/.cache` clears it.

| variable | default | what |
|---|---|---|
| `IGLOO_DIR` | `$HYDRADIR/src/lib/IGLOO` | the IGLOO tree whose case generator builds the IGLOO cases |
| `IGLOO_BIN` | `$IGLOO_DIR/bin/IGLOO` | the IGLOO binary |
| `ICE_BIN` | `$HYDRADIR/src/lib/ICE/bin/ICE` | ICE's standalone binary |

`bin/IGLOO` is one link target shared by every IGLOO build tree (`tests/test.sh`'s `build/verif`,
`install.sh`'s `build`, an MPI tree): whichever tree linked last provides it; `ldd bin/IGLOO | grep -c
libmpi` must print 0. hydra's own build does not write `src/lib/ICE/bin/ICE` (it
`add_subdirectory`s ICE): build ICE standalone first,

```bash
(cd src/lib/ICE && ./install.sh build --compilers=intel --use-openmp --include-orion=$HYDRADIR/lib/ORION --include-finer=$HYDRADIR/lib/third_party/FiNeR)
```

Exit codes: 0 the table is formed (with `--check`, also every value within tolerance); 1 a solver
failed, or with `--check` IGLOO's cloud departs from the closed form at the quarter turn by more than
5e-6 (in `|C|` and in lag [rad]: the per-parcel gate's 3e-6 m over X0 = 0.6 m) or ICE's dilation or
lag by more than `run.py`'s 5e-3 at its own final time; 77 (skipped) `HYDRADIR` unset, or ICE's
verification scripts, an ICE binary, IGLOO's generator or an IGLOO binary missing.

## Reference

[`reference/table.md`](reference/table.md) is the table this driver prints (first made with ICE `3187afb`
and an IGLOO with `time-end`/`out-time`; ICE `494668b` and IGLOO `e5faf0a` write the same table and the same
figures, byte for byte), with the two figures beside it: `vortex-3way-transition.svg` (lag and
`|C|` against t per St: closed form, ICE at its eight outputs, IGLOO from its time-stamped rows) and
`vortex-3way-profile.svg` (the quarter-turn density along the ray through the exact centroid). In
short, at the quarter turn:

| St | exact \|C\| / lag (t = QUARTER) | IGLOO parcels | exact \|C\| / lag (t = 1.571) | ICE | ICE L1 | IGLOO deposit L1 (kernel floor) |
|---|---|---|---|---|---|---|
| 0.01 | 1.015723 / 0.000310 | 1.015723 / 0.000310 | 1.015725 / 0.000310 | 1.015428 / 0.002476 | 0.120 | 0.0212 (0.0212) |
| 0.1 | 1.151551 / 0.025864 | 1.151551 / 0.025864 | 1.151573 / 0.025868 | 1.151011 / 0.027484 | 0.104 | 0.0170 (0.0170) |
| 1 | 1.617416 / 0.345426 | 1.617416 / 0.345426 | 1.617528 / 0.345502 | 1.616801 / 0.345841 | 0.067 | 0.0088 (0.0088) |
| 10 | 1.828830 / 0.538470 | 1.828830 / 0.538470 | 1.828992 / 0.538606 | 1.828395 / 0.538749 | 0.062 | 0.0107 (0.0107) |

IGLOO's parcels sit on the closed form to the print precision (every row within 7.1e-7 m of
`C(t) z0`; the cloud's `|C|` and lag within 2e-8); its deposited field's L1 is the deposition
kernel's alone. ICE is within 7.3e-4 in `|C|` and 2.2e-3 rad in lag, its L1 (0.06–0.12 at 96²)
coming from the numerical diffusion of its second-order transport, visible as the lower peak in the
profile figure.
