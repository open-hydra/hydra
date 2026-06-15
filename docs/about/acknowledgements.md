# Acknowledgements

Hydra is made possible by a broader open-source ecosystem of solver libraries, supporting tools, and contributor communities.

## Core Scientific and Runtime Libraries

These projects provide core capabilities used by Hydra and its solver libraries.

| Project | Repository | Role in hydra ecosystem | License |
|---------|------------|--------------------------|---------|
| FiNeR | [github.com/szaghi/FiNeR](https://github.com/szaghi/FiNeR) | INI parsing and configuration management | GPL v3.0 |
| ORION | [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format scientific I/O (Tecplot, VTK, Plot3D) | GPL v3.0 |
| FLINT | [github.com/MarcoGrossi92/FLINT](https://github.com/MarcoGrossi92/FLINT) | Thermodynamic and chemistry data/services | GPL v3.0 |
| OSLO | [github.com/MarcoGrossi92/OSLO](https://github.com/MarcoGrossi92/OSLO) | ODE integration infrastructure for stiff source terms | GPL v3.0 |

## Solver and Workflow Ecosystem

Hydra integrates or coordinates multiple solver-oriented projects:

| Project | Repository | Contribution |
|---------|------------|--------------|
| MOSE | [github.com/open-hydra/MOSE](https://github.com/open-hydra/MOSE) | Compressible fluid-flow solver library |
| FUSS | [github.com/open-hydra/FUSS](https://github.com/open-hydra/FUSS) | Solid heat-conduction solver library |
| GROOT | [github.com/open-hydra/GROOT](https://github.com/open-hydra/GROOT) | Radiative transfer equation solver |
| ATLAS | [github.com/open-hydra/ATLAS](https://github.com/open-hydra/ATLAS) | Pre-processing workflow for case preparation |

## Documentation Tooling

| Project | Website | Contribution | License |
|---------|---------|--------------|---------|
| MkDocs | [mkdocs.org](https://www.mkdocs.org) | Static documentation generation | BSD-2-Clause |
| Material for MkDocs | [squidfunk.github.io/mkdocs-material](https://squidfunk.github.io/mkdocs-material/) | Documentation theme and navigation/UI components | MIT |

## Community and Institutional Context

Hydra development benefits from:

- Research and engineering feedback from CFD users and contributors
- Open scientific software communities, especially in modern Fortran
- Academic and institutional environments that support verification and validation work

## License and Attribution Note

Hydra respects the licenses of all integrated and external dependencies.
For project-level licensing terms, see [License](license.md).

When redistributing Hydra or derivative works, include appropriate notices and preserve attribution for upstream dependencies.

---