# Contributing

Thanks for contributing to Hydra.

This repository is a framework/integration layer, so contributions may
land in either:

- Hydra framework code (`src/app`, `scripts`, docs)
- solver libraries (`src/lib/MOSE`, `src/lib/FUSS`, etc) --- submit to solver repos first; once released, raise the
  range in `cmake/Releases.cmake` if hydra needs the new release (a build takes the newest in range)

Choose the repository according to ownership of the change.

## Development Workflow

1. Create a branch from `main`.
2. Keep changes focused (one topic per PR when possible).
3. Run framework checks and at least one benchmark run.
4. Open a PR with clear scope and validation notes.

## Local Setup

```bash
git clone https://github.com/open-hydra/hydra.git
cd hydra
source ./bootstrap.sh
```

## Build/Compile Loop

```bash
figaro build hydra-MF --compilers=gnu --use-openmp
figaro compile hydra-MF
```

Useful maintenance commands:

```bash
figaro doctor
figaro deps --update MOSE      # the newest MOSE release in hydra's range, into deps.lock
```

## Where to Place Changes

| Change type | Preferred location |
|-------------|--------------------|
| Coupling loop, shared runtime, HYDRA-* input | `src/app/common/*`, `src/app/MF/*` |
| Build/CLI automation | `scripts/*` |
| Solver physics/numerics internals | corresponding solver repository under `src/lib/*` |
| Framework docs | `docs/*` in this repository |

## Coding and Style Guidelines

### Fortran (framework sources)

- Use `implicit none` in every module/subroutine.
- Declare argument intent (`in`, `out`, `inout`).
- Keep naming consistent with existing hydra modules (`hydra*` prefixes).
- Preserve existing formatting style in touched files.

### Shell scripts

- Use `set -euo pipefail` for task scripts.
- Prefer reusable helpers in `scripts/utils/`.
- Keep command interfaces consistent with `figaro` conventions.

### Documentation

- Keep framework docs hydra-focused.
- Link to solver docs for solver-deep details.
- Avoid duplicating full solver manuals in this repository.

The site is built with [Zensical](https://zensical.org) from `docs/` and `mkdocs.yml`:

```bash
pip install -r requirements-docs.txt
zensical serve            # preview at http://localhost:8000
zensical build --strict   # what CI checks: fails on a broken link or anchor
```

Write links relative to the page: a link that starts with `/` would leave the published
site, and CI refuses it. A new page goes into the `nav` of `mkdocs.yml`. A figure goes inline
with a snippet, `--8<-- "docs/<path>.svg"`, inside a `<figure>`, so it follows the light and
dark themes. The home page is `docs/index.md`, plain Markdown styled by `docs/stylesheets/extra.css`; its
Hydra scheme is `docs/stylesheets/hydra.css`, the same file in every code of the suite.

## Pull Request Checklist

- [ ] Change scope is clear (framework vs solver repository)
- [ ] Build succeeds (`figaro build ...`, `figaro compile ...`)
- [ ] Relevant benchmark run completed
- [ ] Documentation updated for user-visible behavior
- [ ] No unrelated formatting or refactoring noise

## Common Contribution Scenarios

### A. hydra coupling behavior update

1. Modify `src/app/MF/*`.
2. Run `figaro compile hydra-MF`.
3. Re-run `test/MF/flatplate-Tw` benchmark cases.
4. Include before/after behavior notes in PR.

### B. HYDRA-* input extension

1. Add registration in `src/app/common/config/Register_*.f90`.
2. Wire runtime usage in framework code.
3. Update user/developer docs.
4. Validate with at least one case.

### C. Solver-level model change

If the change is purely solver internals (for example MOSE numerics),
prefer submitting it in the solver repository first. While it is unreleased, build hydra against your checkout with
`-DFETCHCONTENT_SOURCE_DIR_<NAME>=<folder>`; once released, take it (`figaro deps --update <NAME>`) and,
if hydra needs it, the lower limit in `cmake/Releases.cmake`.
