---
title: hydra
hide:
  - navigation
  - toc
---

<div class="home-hero" markdown>

<div class="home-hero__text" markdown>

<p class="home-hero__eyebrow">Multi-Physics Framework</p>

<h1 class="home-hero__title">Hydra</h1>

<p class="home-hero__lead">Hydra composes specialized solver libraries for fluid flow, solid heat transfer, and
advanced reacting formulations, with support for coupled executions such as fluid-solid conjugate workflows.</p>

<p class="home-hero__cmd"><span>$</span> figaro hydra-MF -p 4 solve</p>

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[Verification](vv/index.md){ .md-button .home-hero__ghost }

<p class="home-hero__name"><strong>A multi-solver CFD and multi-physics framework</strong>: fluid, solid and
radiation solvers, OpenMP + MPI, a <strong>coupled</strong> strategy</p>

</div>

<div class="home-hero__shot" markdown>

<div class="home-term" markdown>

<div class="home-term__bar"><span></span><span></span><span></span></div>

```bash
# Clone and build
git clone https://github.com/open-hydra/hydra.git
cd hydra && source bootstrap.sh
figaro build hydra-MF

# Run a case
figaro hydra-MF solve
```

</div>

<p class="home-hero__caption" markdown>Start with Hydra, then choose your solver: the
[Installation guide](getting-started/installation.md) and the [Quick Start](getting-started/quick-start.md)</p>

</div>

</div>

<p class="home-kicker">The framework</p>

## Hydra at a glance { .home-heading }

<p class="home-sub" markdown>Dedicated solvers, shared workflow, and optional cross-solver coupling. Build a Hydra
executable, run a case, and then dive into solver-level docs for detailed models and input options: the
[Overview](getting-started/index.md#core-idea) names the solvers and the coupled executables.</p>

<div class="home-tools" markdown>

[:material-weather-windy: <b>Compressible Fluid Flow</b> <i>Ideal-gas Euler/Navier-Stokes, turbulence models, and reacting-flow capabilities on structured multi-block grids.</i>](theory/index.md#solver-physics-at-a-glance){ .home-tool }
[:material-thermometer: <b>Solid Heat Conduction</b> <i>Steady and transient solid conduction with multi-material support and boundary models for conjugate workflows.</i>](theory/index.md#solver-physics-at-a-glance){ .home-tool }
[:material-white-balance-sunny: <b>Radiation Transport</b> <i>Radiation transport equation for radiative heat transfer applications.</i>](theory/index.md#solver-physics-at-a-glance){ .home-tool }
[:material-swap-horizontal: <b>Coupled Execution</b> <i>Hydra applications can coordinate multiple solvers in one run, exchanging interface quantities at coupling boundaries.</i>](theory/index.md#coupled-theory-in-hydra){ .home-tool }
[:material-chip: <b>Scalable Runtime</b> <i>OpenMP and MPI support across the framework and solver libraries, suitable for workstation and HPC runs.</i>](user/parallel.md){ .home-tool }
[:material-console: <b>Single Workflow</b> <i>Common Hydra commands for build and execution, with solver-specific physics and documentation linked from one hub.</i>](user/figaro.md){ .home-tool }

</div>

<p class="home-kicker"></p>

## Hydra CFD Suite { .home-heading }

<p class="home-sub">Hydra couples solvers of gases, real fluids, solids, particles and radiation in one
simulation. It integrates specialized solver libraries and optional pre/post-processing tools; additional solvers
are under development.</p>

<div class="hydra" markdown>

<div class="hydra__case" markdown>

<p class="hydra__label">Your case</p>

<p class="hydra__file" markdown>:material-grid: structured mesh</p>

<p class="hydra__file" markdown>:material-database: [databases](https://open-hydra.github.io/ATLAS/databases/)</p>

</div>

<div class="hydra__atlas" markdown>

<p class="hydra__name">ATLAS</p>

<p class="hydra__tools" markdown>[GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/) [BCB](https://open-hydra.github.io/ATLAS/user-guide/bcb/) [ICB](https://open-hydra.github.io/ATLAS/user-guide/icb/) [STB](https://open-hydra.github.io/ATLAS/user-guide/stb/) [MDB](https://open-hydra.github.io/ATLAS/user-guide/mdb/) [KAnT](https://open-hydra.github.io/ATLAS/user-guide/kant/)</p>

</div>

<div class="hydra__link">
<span class="hydra__link-line"></span>
<span class="hydra__link-bottom">input files for each run</span>
</div>

<div class="hydra__suite is-here" markdown>

<p class="hydra__label" markdown>[Hydra](https://github.com/open-hydra) · the coupled simulation</p>

<div class="hydra__solvers" markdown>

[<b>MOSE</b> compressible flow <i>ideal-gas</i>](https://github.com/open-hydra/MOSE){ .hydra__solver .hydra__solver--gas }
[<b>GROOT</b> radiative transfer <i>no ATLAS phase</i>](https://github.com/open-hydra/GROOT){ .hydra__solver .hydra__solver--radiation }
[<b>ARES</b> real-fluid flow <i>real-fluid</i>](https://github.com/open-hydra/ARES){ .hydra__solver .hydra__solver--fluid }
[<b>FUSS</b> heat conduction in solids <i>solid</i>](https://github.com/open-hydra/FUSS){ .hydra__solver .hydra__solver--solid }
[<b>IGLOO</b> Lagrangian particles <i>dispersed</i>](https://github.com/open-hydra/IGLOO){ .hydra__solver .hydra__solver--particle }
[<b>ICE</b> Eulerian particles <i>dispersed</i>](https://github.com/open-hydra/ICE){ .hydra__solver .hydra__solver--particle }

</div>

</div>

<div class="hydra__libs" markdown>

<p class="hydra__label">Built on</p>

[<b>ORION</b> Tecplot, VTK and PLOT3D I/O](https://MarcoGrossi92.github.io/ORION/){ .hydra__lib }
[<b>FLINT</b> thermochemistry](https://MarcoGrossi92.github.io/FLINT/){ .hydra__lib }
[<b>OSLO</b> ODE integrators](https://github.com/MarcoGrossi92/OSLO){ .hydra__lib }
[<b>FiNeR</b> INI file parser](https://github.com/szaghi/FiNeR){ .hydra__lib }

</div>

</div>

<div class="home-links" markdown>

[Quick Start](getting-started/quick-start.md)
[Overview](getting-started/index.md#core-idea)
[User Guide](user/index.md)
[Input Parameters](user/registry.md)
[Theory](theory/index.md)
[Verification](vv/index.md)
[Developer Guide](development/index.md)

</div>
