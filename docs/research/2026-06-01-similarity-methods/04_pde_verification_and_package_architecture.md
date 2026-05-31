# PDE Verification And Package Architecture

Date: 2026-06-01

Bead: `scr-x9r`

Scope: best practices for verifying similarity solutions against
time-dependent PDE/free-boundary simulations and for evolving this repository
into a robust Julia package. This report separates source-backed facts from
implementation recommendations.

## Executive Findings

- The present fixed-grid method-of-lines verifier is a useful first check, but
  cone recoil needs a wider verification ladder: exact/manufactured tests,
  residual and conservation diagnostics, grid/time convergence, quantitative
  collapse metrics, and at least one independent discretization.
- Moving-coordinate or moving-mesh PDE verification is the highest-leverage
  upgrade for the current limitation that the fixed grid only reaches short
  times before tip stiffness dominates.
- Boundary-integral/free-surface computations are the natural independent
  high-fidelity comparator for inviscid potential-flow cone/wedge recoil, while
  slender-jet PDEs are the natural lower-dimensional regression family.
- The package should keep a small dependency core, expose SciML-style problem
  objects and `solve` methods, and place ModelingToolkit, MethodOfLines,
  ApproxFun, plotting, benchmark, and artifact support behind extensions or
  scripts.
- Validation benchmarks should be data products with recorded parameters,
  solver options, mesh/time settings, diagnostics, and source citations, not
  just visual figures.

## Source-Backed Facts

### Verification And Method Of Lines

The method of lines discretizes spatial variables and leaves time continuous,
turning a PDE into ODE/DAE systems that can be advanced with standard time
integrators. Scholarpedia's method-of-lines entry is citable as
DOI [`10.4249/scholarpedia.2859`](https://doi.org/10.4249/scholarpedia.2859).
This matches the current repository approach in `src/pde.jl`: nonuniform finite
differences in `z` plus an implicit ODE solve.

Oberkampf and Trucano distinguish code verification, solution verification,
model validation, errors, and uncertainties in CFD verification and validation
work. See *Progress in Aerospace Sciences* 38, 209-272,
DOI [`10.1016/S0376-0421(02)00005-2`](https://doi.org/10.1016/S0376-0421(02)00005-2),
and Sandia report SAND2002-0529,
DOI [`10.2172/793406`](https://doi.org/10.2172/793406). This distinction is
important here: matching the similarity ODE only verifies consistency of a
model/discretization, while matching Decent-King or experiments is validation
against a target physical/asymptotic model.

The method of manufactured solutions is a standard code-verification technique:
choose an analytic solution, insert it into the PDE, and add the resulting
source terms so convergence can be measured against a known answer. Roache's
"Code Verification by the Method of Manufactured Solutions" is
DOI [`10.1115/1.1436090`](https://doi.org/10.1115/1.1436090). MMS is especially
useful for this repository because the physical cone problem has singular
initial geometry and moving scales, but the implementation of derivatives,
boundary conditions, and source terms can be verified on smooth manufactured
profiles first.

### Moving Meshes, ALE, And Free Boundaries

Hirt, Amsden, and Cook introduced an arbitrary Lagrangian-Eulerian computing
method in which the grid can move independently of the material,
DOI [`10.1016/0021-9991(74)90051-5`](https://doi.org/10.1016/0021-9991(74)90051-5).
Donea, Giuliani, and Halleux applied an ALE finite-element formulation to
transient fluid-structure/free-interface dynamics and describe the grid points
as displaceable independently of fluid motion, improving interface treatment
without excessive mesh distortion,
DOI [`10.1016/0045-7825(82)90128-1`](https://doi.org/10.1016/0045-7825(82)90128-1).

Adaptive moving mesh PDEs based on equidistribution were developed by Huang,
Ren, and Russell, DOI [`10.1137/0731038`](https://doi.org/10.1137/0731038).
Huang and Russell's book *Adaptive Moving Mesh Methods* is a stable reference
for time-dependent PDEs with sharp fronts,
DOI [`10.1007/978-1-4419-7916-2`](https://doi.org/10.1007/978-1-4419-7916-2).
The source-backed implication for cone recoil is that the mesh should follow
the moving capillary-wave region and tip scale instead of resolving a rapidly
expanding physical domain with a fixed grid.

### Free-Surface Similarity And Slender-Jet Benchmarks

Keller and Miksis study time-dependent potential free-surface flows driven by
surface tension, reduce self-similar configurations to integrodifferential
systems for the free surface and potential, and compute surface waves. See
*SIAM Journal on Applied Mathematics* 43(2), 268-277,
DOI [`10.1137/0143018`](https://doi.org/10.1137/0143018).

Ting and Keller derive simplified equations for slender jets and thin sheets
with surface tension and introduce similarity solutions governed by ODEs.
See *SIAM Journal on Applied Mathematics* 50(6), 1533-1546,
DOI [`10.1137/0150090`](https://doi.org/10.1137/0150090).

Eggers and Dupont derive a one-dimensional thin, axisymmetric free-surface
Navier-Stokes approximation, compare with experiments, and study singularities
formed during pinch-off. See *Journal of Fluid Mechanics* 262, 205-221,
DOI [`10.1017/S0022112094000480`](https://doi.org/10.1017/S0022112094000480),
with arXiv record [`physics/0110081`](https://arxiv.org/abs/physics/0110081).

Papageorgiou constructs similarity solutions for inviscid, Stokes, and
Navier-Stokes jet breakup and uses slender-jet model simulations to show
selection of symmetric pinching solutions for general initial conditions.
See *Journal of Fluid Mechanics* 301, 109-132,
DOI [`10.1017/S002211209500382X`](https://doi.org/10.1017/S002211209500382X),
and the related viscous-thread model paper in *Physics of Fluids* 7, 1529-1544,
DOI [`10.1063/1.868540`](https://doi.org/10.1063/1.868540).

Eggers' review identifies one-dimensional equations and universal scaling laws
as central to free-surface breakup near singularities. See *Reviews of Modern
Physics* 69, 865-929,
DOI [`10.1103/RevModPhys.69.865`](https://doi.org/10.1103/RevModPhys.69.865).

Decent and King are the primary cone-recoil target for this package. Their
abstract states that a slender post-breakup cone is evolved by surface tension,
an asymptotic solution is found in the small cone aspect ratio, a similarity
transformation is valid for small post-bifurcation times, and a rapidly
oscillating nonlinear wave propagates away from the tip. See *IMA Journal of
Applied Mathematics* 73(1), 37-68,
DOI [`10.1093/imamat/hxm043`](https://doi.org/10.1093/imamat/hxm043).

Billingham's fat wedge/cone paper is a useful companion benchmark: its abstract
reports a long-time similarity form with `t^(2/3)` deformations and capillary
waves, with viscous damping at distances `O(t^(3/4))` from the tip. See
*Journal of Fluid Mechanics* 397, 45-71. Cambridge currently reports
DOI [`10.1017/S0022112099006047`](https://doi.org/10.1017/S0022112099006047).

Boundary-integral methods are a credible independent route for inviscid
free-surface validation. For example, Garzon, Sethian, and Gray couple a level
set representation of an inviscid pinch-off interface to a Galerkin boundary
integral solution of the axisymmetric Laplace problem, with validation against
a short-time Rayleigh-Taylor analytical solution. Stable ORNL record:
<https://www.ornl.gov/publication/numerical-simulation-non-viscous-liquid-pinch-using-coupled-level-set-boundary-integral>.
For axisymmetric liquid bridges, Volkov, Papageorgiou, and Petropoulos derive
boundary integral methods with fast series representations,
DOI [`10.1137/040604352`](https://doi.org/10.1137/040604352).

### Similarity Collapse Metrics

Bhattacharjee and Seno propose an objective quality measure for data collapse,
then minimize it to extract scaling exponents and error bars. See *Journal of
Physics A* 34, 6375-6380,
DOI [`10.1088/0305-4470/34/33/302`](https://doi.org/10.1088/0305-4470/34/33/302),
and arXiv [`cond-mat/0102515`](https://arxiv.org/abs/cond-mat/0102515).
This gives a source-backed precedent for replacing visual "snapshot collapse"
with an optimization-based scalar score.

### Julia Package Architecture And Tooling

DifferentialEquations.jl uses the `solve(prob, alg; kwargs...)` pattern for
`ODEProblem`s. Its current solver docs recommend `Rodas5P`, `Rodas4P`,
`Kvaerno5`, or `KenCarp4` for medium-tolerance stiff problems, and `QNDF` or
`FBDF` for very large or expensive semidiscretized systems; `CVODE_BDF` is also
listed as an option. Stable docs:
<https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/>.

ModelingToolkit's `PDESystem` stores equations, boundary conditions, domains,
independent/dependent variables, parameters, defaults, and optional analytic
solutions. Its PDE interface separates symbolic specification from numerical
discretization, and `discretize(sys, discretizer)` produces an `AbstractSystem`
or SciML problem. The docs list MethodOfLines.jl as the finite-difference
discretizer and note support for higher-order stencils and nonuniform grids:
<https://docs.sciml.ai/ModelingToolkit/v8.56/systems/PDESystem/>.

MethodOfLines.jl's `MOLFiniteDifference` accepts nonuniform rectilinear grids
as vectors, approximation order, advection schemes including upwind and WENO,
grid alignment choices, and optional ODAE construction:
<https://docs.sciml.ai/MethodOfLines/dev/api/discretization/>.

ApproxFun.jl represents differential and integral operators on function spaces,
including Chebyshev and Fourier spaces, and boundary functionals such as
Dirichlet operators. Its operator model is a good fit for independent spectral
collocation/residual checks of similarity ODE/BVPs:
<https://juliaapproximation.github.io/ApproxFun.jl/latest/usage/operators/>.

Julia's Pkg docs state that a `Project.toml` plus `Manifest.toml` can
instantiate the exact package environment, and that the manifest records exact
direct and indirect dependency information:
<https://pkgdocs.julialang.org/dev/toml-files/>. Pkg artifacts are immutable
content-addressed data containers declared in `Artifacts.toml`, suitable for
datasets, text, platform binaries, and other non-package data:
<https://pkgdocs.julialang.org/v1.9/artifacts/>.

PkgBenchmark expects a package benchmark suite at
`benchmark/benchmarks.jl` with a `SUITE` variable using BenchmarkTools'
dictionary-based interface:
<https://juliaci.github.io/PkgBenchmark.jl/stable/define_benchmarks/>.
BenchmarkTools supports `BenchmarkGroup` suites and recommends seeded random
inputs when randomness is used:
<https://juliaci.github.io/BenchmarkTools.jl/stable/manual/>.

## Recommendations

The following recommendations are inferences from the sources above and the
current repository state, not direct claims made by any one source.

### Verification Ladder

1. Add manufactured-solution tests for the current slender PDE operator.
   Include `R`, `u`, source terms, boundary data, and exact residuals on both
   uniform and stretched grids. These tests should measure observed spatial
   order independently of the physical cone singularity.

2. Add exact/simple similarity checks before cone recoil. Candidate checks:
   heat-equation Gaussian collapse, a smooth forced Burgers/advection-diffusion
   MMS profile, and a linear capillary-wave dispersion test on a cylinder or
   cone far field.

3. Keep the current fixed-grid MOL solver as a baseline, but compute diagnostics
   at every saved snapshot: original PDE residual, similarity ODE residual
   after rescaling, mass/flux balance over the truncated domain, positivity
   margin, stiffness/retcode statistics, and grid-quality measures.

4. Add a mapped-coordinate verifier. The simplest target is a similarity-frame
   PDE with `xi = z / ell(t)`, `ell(t) = t^(2/3)` or
   `ell(t) = a*(t - t0)^alpha`, plus optional exponent fitting. In the correct
   similarity frame, the PDE evolution should approach a steady profile; this
   avoids chasing an expanding physical length scale.

5. Add a true moving-mesh/ALE verifier after the mapped-coordinate solver is
   stable. Use a computational coordinate `eta in [0, 1]` and a map
   `z = Z(eta, t)`. Start with prescribed maps tied to `ell(t)`, then test
   equidistributed maps with monitor functions based on curvature, wave phase,
   residual, and local radius. Keep remeshing interpolation conservative for
   area `A = R^2`.

6. Add one independent discretization. Good first choices are:
   MethodOfLines-generated finite differences for a symbolic version of the
   PDE, and ApproxFun spectral residual/collocation for similarity ODE/BVP
   problems. A boundary-integral/free-surface solver is a larger project but
   would be the strongest independent inviscid check.

### PDE Formulation Choices

- Prefer conservative variables for verification runs: area `A = R^2` and
  momentum-like `A*u` where the model permits it. Primitive `R, u` is convenient
  but makes conservation diagnostics less direct and can hide mass defects.
- Treat `R <= 0`, nonfinite states, and failed retcodes as hard failures in
  benchmark runs. Silent clipping would invalidate singularity and collapse
  diagnostics.
- For stretched-grid finite differences, compare derivative operators against
  analytic functions with nonuniform-grid order tests before using them in PDE
  claims.
- For stiff capillary terms, benchmark `Rodas5P`, `TRBDF2`, `FBDF`, `QNDF`,
  and `CVODE_BDF` on the same semidiscrete problem. Record tolerances, Jacobian
  strategy, linear solver, accepted/rejected steps, and minimum radius.
- When using axial curvature or third derivatives, prefer banded/sparse
  Jacobian paths once the implementation moves beyond small grids.

### Collapse And Similarity Diagnostics

Use collapse as a quantitative regression object, not a plot-only check.

- Define a common trusted window in similarity space that excludes the
  truncated tip boundary, the far outflow boundary, and regions where the
  physical PDE has not yet resolved the similarity structure.
- For each snapshot, rescale
  `xi = z / ell(t)`, `S = R / ell(t)`, and
  `U = u / ell_dot(t)`. Interpolate all snapshots to a common `xi` grid.
- Compute weighted pairwise variance or deviation from a reference profile:
  `E_collapse = sum_i w_i ||S_i - mean(S)||^2 / sum_i w_i ||mean(S)||^2`.
  Use separate norms for profile, slope, curvature, and velocity.
- Optimize optional exponents and time offset, for example
  `ell(t) = a*(t - t0)^alpha`, and report confidence intervals or sensitivity
  bands. Bhattacharjee-Seno gives a precedent for objective collapse quality
  and exponent/error-bar extraction.
- Track capillary-wave phase and envelope separately: crest locations in `xi`,
  zero crossings of `S - epsilon*xi`, and envelope decay are often more
  discriminating than a bulk `L2` norm.
- Always report conservation and residual diagnostics beside the collapse
  score. A visually good collapse with bad mass balance is not a valid
  verification result.

### Concrete Validation Benchmarks

Recommended benchmark set, in increasing problem-specificity:

| Benchmark | Purpose | Source anchor | Candidate recorded quantities |
| --- | --- | --- | --- |
| Smooth manufactured slender PDE | Verify finite differences, BCs, source handling, solver plumbing | Roache MMS DOI `10.1115/1.1436090` | observed order, residual norms, boundary residuals |
| Linear capillary wave on a nearly cylindrical/slender far field | Verify axial curvature dispersion and wave speed/phase | Eggers review DOI `10.1103/RevModPhys.69.865`; Ting-Keller DOI `10.1137/0150090` | phase speed, amplitude error, damping if viscosity included |
| Eggers-Dupont slender jet breakup | Verify 1D free-surface PDE infrastructure and singular tracking | DOI `10.1017/S0022112094000480` | neck radius scaling, breakup time convergence, profile collapse |
| Papageorgiou viscous thread | Verify selection of self-similar pinching profiles | DOI `10.1063/1.868540`; DOI `10.1017/S002211209500382X` | similarity exponents, symmetric terminal profile, satellite/drop metrics if modeled |
| Keller-Miksis surface-tension flow | Check inviscid self-similar wave generation against a boundary-integral/asymptotic target | DOI `10.1137/0143018` | similarity profile, wave train, boundary-integral residual |
| Billingham fat wedge/cone | Test viscosity, fat-angle asymptotics, and wave damping | DOI `10.1017/S0022112099006047` | `t^(2/3)` scaling, damping distance, tip velocity behavior |
| Decent-King slender cone | Primary package target | DOI `10.1093/imamat/hxm043` | inner BVP parameters, far-field matching residual, capillary-wave phase/envelope, PDE collapse score |
| Keller-King-Ting blob formation | Bridge between local blob asymptotics and cone/jet breakup | DOI `10.1063/1.868723` | blob shape, local scaling, comparison to cone-tip limiting behavior |

Each benchmark should have a "cheap" regression configuration and a
"reference" configuration. Reference data can be stored as small text/CSV/JLD2
artifacts with hashes, while large licensed PDFs remain outside git with only
manifest/checksum metadata.

### Package Architecture

Recommended modules or namespaces:

- `Problems`: `SimilarityProblem`, `SlenderPDEProblem`,
  `PDEVerificationProblem`, `CollapseProblem`.
- `Models`: conservative and primitive slender models, curvature operators,
  nondimensionalization, source terms for MMS.
- `Solvers`: inner BVP shooting/collocation, outer/matching solvers,
  PDE/MOL backends, moving-coordinate backends.
- `Diagnostics`: residuals, conservation, collapse metrics, grid/time
  convergence, wave phase/envelope, solver retcode summaries.
- `Benchmarks`: problem constructors and reference tolerances, not large data.
- `Artifacts`: helpers for locating immutable reference datasets.

API recommendations:

- Provide constructors that return SciML problems where appropriate:
  `ODEProblem`, `NonlinearProblem`, or `BVProblem` once a BVP backend is used.
- Define package-level `solve(prob::SlenderPDEProblem, alg; kwargs...)` methods
  that delegate to SciML rather than exposing internal arrays as the main API.
- Return rich result structs: solution, physical parameters, nondimensional
  parameters, mesh, solver diagnostics, residual diagnostics, conservation
  diagnostics, and source/reference metadata.
- Keep `DifferentialEquations` or narrower OrdinaryDiffEq/SciMLBase
  dependencies in the core only if needed by the default API. Put
  ModelingToolkit, MethodOfLines, ApproxFun, BenchmarkTools, and plotting in
  extensions, test extras, or script environments.
- Use `Project.toml`/`Manifest.toml` for reproducible script and benchmark
  environments. For future Julia versions, consider the Pkg docs' warning that
  `[extras]`/`[targets]` is legacy for Julia 1.13+, but do not churn the current
  Julia 1.10 package layout just for that.
- Add `benchmark/benchmarks.jl` with a `SUITE` for solver performance:
  inner BVP solve, PDE short run, PDE reference run, collapse metric, residual
  evaluation, and figure-data generation. Benchmarks should seed any random
  inputs and record BLAS/thread settings.

### Reproducibility Design

- Every generated validation figure should have a companion metadata record:
  git commit, Julia version, project/manifest hash, solver algorithms,
  tolerances, mesh settings, physical parameters, reference citation, and
  diagnostics.
- Store reference numeric data as immutable artifacts when they are too large
  or too specialized for normal source control. Use `Artifacts.toml` with
  content hashes and download hashes; do not put licensed publisher PDFs in
  artifacts.
- Separate "paper provenance" from "computed benchmark provenance": PDFs and
  page images prove source access, while benchmark artifacts prove numerical
  reproducibility.
- For collaboration, treat benchmark updates like API changes: require a
  reason, source citation, old/new diagnostics, and an explicit tolerance
  decision.

## Access Gaps And Follow-Up Checks

- The Decent-King cone paper is marked missing in `docs/papers/README.md`; this
  report used the OUP/Crossref-facing abstract and metadata, not a local PDF.
  Quantitative Decent-King benchmark values still require the paper.
- This research pass found that `docs/papers/README.md` previously listed
  Billingham 1999 with DOI `10.1017/S0022112099006011`, but Cambridge's article
  page for "Surface-tension-driven flow in fat fluid wedges and cones" reports
  `10.1017/S0022112099006047`; `...6011` belongs to a different JFM article on
  gravity modulation. The manifest was corrected in the same research batch.
- Keller-King-Ting "Blob formation" has only a local page image per the paper
  manifest. Acquire the PDF before using quantitative blob data.
- I did not run Julia commands or inspect ignored licensed PDFs. The only
  validation command planned for this report is `git diff --check`.

## Sources Used

- Schiesser and Griffiths, "Method of lines", Scholarpedia,
  DOI [`10.4249/scholarpedia.2859`](https://doi.org/10.4249/scholarpedia.2859).
- Oberkampf and Trucano, "Verification and validation in computational fluid
  dynamics", DOI [`10.1016/S0376-0421(02)00005-2`](https://doi.org/10.1016/S0376-0421(02)00005-2);
  Sandia record DOI [`10.2172/793406`](https://doi.org/10.2172/793406).
- Roache, "Code Verification by the Method of Manufactured Solutions",
  DOI [`10.1115/1.1436090`](https://doi.org/10.1115/1.1436090).
- Hirt, Amsden, and Cook, "An Arbitrary Lagrangian-Eulerian Computing Method
  for All Flow Speeds",
  DOI [`10.1016/0021-9991(74)90051-5`](https://doi.org/10.1016/0021-9991(74)90051-5).
- Donea, Giuliani, and Halleux, "An arbitrary lagrangian-eulerian finite
  element method for transient dynamic fluid-structure interactions",
  DOI [`10.1016/0045-7825(82)90128-1`](https://doi.org/10.1016/0045-7825(82)90128-1).
- Huang, Ren, and Russell, "Moving Mesh Partial Differential Equations
  (MMPDEs) Based on the Equidistribution Principle",
  DOI [`10.1137/0731038`](https://doi.org/10.1137/0731038).
- Huang and Russell, *Adaptive Moving Mesh Methods*,
  DOI [`10.1007/978-1-4419-7916-2`](https://doi.org/10.1007/978-1-4419-7916-2).
- Keller and Miksis, "Surface tension driven flows",
  DOI [`10.1137/0143018`](https://doi.org/10.1137/0143018).
- Ting and Keller, "Slender Jets and Thin Sheets with Surface Tension",
  DOI [`10.1137/0150090`](https://doi.org/10.1137/0150090).
- Eggers and Dupont, "Drop Formation in a One-Dimensional Approximation of the
  Navier-Stokes Equation",
  DOI [`10.1017/S0022112094000480`](https://doi.org/10.1017/S0022112094000480);
  arXiv [`physics/0110081`](https://arxiv.org/abs/physics/0110081).
- Papageorgiou, "Analytical description of the breakup of liquid jets",
  DOI [`10.1017/S002211209500382X`](https://doi.org/10.1017/S002211209500382X).
- Papageorgiou, "On the breakup of viscous liquid threads",
  DOI [`10.1063/1.868540`](https://doi.org/10.1063/1.868540).
- Eggers, "Nonlinear dynamics and breakup of free-surface flows",
  DOI [`10.1103/RevModPhys.69.865`](https://doi.org/10.1103/RevModPhys.69.865).
- Decent and King, "Surface-tension-driven flow in a slender cone",
  DOI [`10.1093/imamat/hxm043`](https://doi.org/10.1093/imamat/hxm043).
- Billingham, "Surface-tension-driven flow in fat fluid wedges and cones",
  Cambridge page with DOI
  [`10.1017/S0022112099006047`](https://doi.org/10.1017/S0022112099006047).
- Garzon, Sethian, and Gray, "Numerical simulation of non-viscous liquid
  pinch-off using a coupled level set-boundary integral method", ORNL record:
  <https://www.ornl.gov/publication/numerical-simulation-non-viscous-liquid-pinch-using-coupled-level-set-boundary-integral>.
- Volkov, Papageorgiou, and Petropoulos, "Accurate and Efficient Boundary
  Integral Methods for Electrified Liquid Bridge Problems",
  DOI [`10.1137/040604352`](https://doi.org/10.1137/040604352).
- Bhattacharjee and Seno, "A measure of data collapse for scaling",
  DOI [`10.1088/0305-4470/34/33/302`](https://doi.org/10.1088/0305-4470/34/33/302);
  arXiv [`cond-mat/0102515`](https://arxiv.org/abs/cond-mat/0102515).
- DifferentialEquations.jl solver docs:
  <https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/>.
- ModelingToolkit.jl `PDESystem` docs:
  <https://docs.sciml.ai/ModelingToolkit/v8.56/systems/PDESystem/>.
- MethodOfLines.jl discretization docs:
  <https://docs.sciml.ai/MethodOfLines/dev/api/discretization/>.
- ApproxFun.jl operators docs:
  <https://juliaapproximation.github.io/ApproxFun.jl/latest/usage/operators/>.
- Julia Pkg `Project.toml`/`Manifest.toml` docs:
  <https://pkgdocs.julialang.org/dev/toml-files/>.
- Julia Pkg artifacts docs:
  <https://pkgdocs.julialang.org/v1.9/artifacts/>.
- PkgBenchmark.jl benchmark-suite docs:
  <https://juliaci.github.io/PkgBenchmark.jl/stable/define_benchmarks/>.
- BenchmarkTools.jl manual:
  <https://juliaci.github.io/BenchmarkTools.jl/stable/manual/>.
