# Verification & Validation

This section documents the Verification & Validation (V&V) test suite for MOSE.  The purpose is to demonstrate that the solver produces correct results across a range of physical regimes, by comparison against analytical or well-established reference solutions.

## Test suite

| Test | Dim | Mach | Physics | Verification | Solver(s) | Ref |
|---|---|---|---|---|---|---|

## Running the tests

Each test case lives under `test/` and contains a `verify.py` script.  The script reads the solver output, compares it against the exact solution, and writes output figures.

```bash
cd test/1D/<TestName>
./MOSE.sh solve
python verify.py          # check errors and export figure
python verify.py --plot   # as above, also display figure interactively
```
