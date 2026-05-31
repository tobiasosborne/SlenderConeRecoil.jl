# BVP and Continuation Algorithms for Similarity ODEs

Date: 2026-06-01

Bead: `scr-x9r`

Assigned scope: numerical algorithms and packages for similarity ODE
boundary-value problems, singular and semi-infinite domains, shooting,
multiple shooting, collocation, adaptive mesh/error control, homotopy and
pseudo-arclength continuation, and branch/fold handling.

## Executive Summary

The package should move its production similarity solve away from a
hand-rolled single-shooting Newton loop and toward an adaptive collocation BVP
formulation. In Julia, the most natural default is SciML
`BoundaryValueDiffEq.jl`, with `MIRK4`/`MIRK5`/`MIRK6` as the first production
path, `Ascher*` Gaussian collocation as a lineage-compatible alternative, and
`MultipleShooting` as a diagnostic fallback when event/callback-compatible IVP
integration is needed.

Continuation should be a first-class package concept. Simple parameter
homotopy can be implemented directly by reusing the previous BVP solution as
the next initial guess. Folds and branch switching should not be hand-coded
initially; they are better handled through a discretized residual passed to
`BifurcationKit.jl`, with AUTO-07p used as an external benchmark and design
reference. Chebfun-style spectral collocation is valuable for high-accuracy
reference calculations on smooth finite or mapped domains, but it is not a
drop-in replacement for adaptive mesh collocation on oscillatory, layered, or
singular similarity profiles.

For semi-infinite similarity domains, the reliable default is finite-domain
truncation with source-derived far-field or asymptotic boundary conditions,
mesh adaptation, and continuation in the truncation length. Compactification
and spectral collocation should be exposed as validation tools once the
far-field asymptotics are fixed.

## Source-Backed Facts

### Current Repository Context

`SlenderConeRecoil.jl` currently solves the inner similarity BVP by shooting
from the tip with three unknown shooting parameters, integrating a four-state
ODE with `Rodas5P`, and using a finite-difference damped Newton correction.
The README records that this solve currently reaches only about a 2 percent
far-field slope error and explicitly names continuation or homotopy as an
improvement target. See [`src/inner.jl`](../../../src/inner.jl) and
[`README.md`](../../../README.md).

The current package depends on `DifferentialEquations.jl` but not on
`BoundaryValueDiffEq.jl` or `BifurcationKit.jl`, so adding those capabilities
would be an API and dependency decision rather than just an implementation swap.
See [`Project.toml`](../../../Project.toml).

### BoundaryValueDiffEq.jl and SciML BVP Solvers

SciML's BVP solver documentation recommends MIRK methods for most cases because
they have improved stability, adaptivity, and sparsity handling. It also notes
that single shooting can be faster if it converges, but is much less robust,
while multiple shooting is generally more stable than single shooting. The same
docs list `MIRK2` through `MIRK6`, FIRK Lobatto/Radau collocation methods,
Gauss-Legendre `Ascher1` through `Ascher7` methods, and shooting/multiple
shooting solvers. Source: SciML BVP solver docs:
https://docs.sciml.ai/DiffEqDocs/dev/solvers/bvp_solve/

`BoundaryValueDiffEqShooting` provides `Shooting` and `MultipleShooting`.
`MultipleShooting` reduces the BVP to multiple IVPs, supports a configurable
number of shooting points, and has grid coarsening options. Source:
https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/shooting/

The MIRK methods in `BoundaryValueDiffEqMIRK` have defect-control adaptivity and
use sparse Jacobian handling. The docs cite Enright and Muir's MIRK defect
control work. Source:
https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/mirk/ and DOI
https://doi.org/10.1137/S1064827593251496

The FIRK methods include Radau IIA and Lobatto IIIA/B/C collocation methods,
with defect-control adaptivity for most listed methods. For large or stiff
BVPs, the docs expose `nested_nlsolve=true` to solve implicit FIRK stages with
nested nonlinear solves. Source:
https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/firk/

`BoundaryValueDiffEqAscher` implements Gauss-Legendre collocation methods with
Ascher-style error-control adaptivity and mesh refinement. The solver docs cite
Ascher-Spiteri (1994) and the Ascher-Christiansen-Russell COLSYS lineage.
Source: https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/ascher/

SciML's error-control page says adaptive collocation supports defect control
and global error control, including `DefectControl`, `GlobalErrorControl`,
`SequentialErrorControl`, `HybridErrorControl`, and `NoErrorControl`. Source:
https://docs.sciml.ai/BoundaryValueDiffEq/dev/basics/error_control/

SciML's continuation tutorial demonstrates the practical homotopy pattern:
solve an easier BVP, then reuse each intermediate solution as the initial guess
for the next, harder parameter value. Source:
https://docs.sciml.ai/BoundaryValueDiffEq/dev/tutorials/continuation/

### BifurcationKit.jl

`BifurcationKit.jl` is a Julia continuation and bifurcation package built around
Newton-Krylov solvers, pseudo-arclength continuation, Moore-Penrose and
deflated continuation, polynomial and ANM continuation, event monitoring,
bisection localization of special points, and GPU-compatible workflows. Source:
https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/capabilities/

For equilibria, BifurcationKit supports detection of branch, fold, and Hopf
points, automatic branch switching at branch points, fold/Hopf continuation
through minimally augmented formulations, codimension-2 detection, normal forms,
and automatic bifurcation diagrams. Source:
https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/capabilities/

Its PALC documentation formulates the augmented pseudo-arclength problem by
adding a hyperplane constraint to `F(x,p)=0`; the corrector solves the enlarged
system by Newton iterations, and the step size is reduced on corrector failure.
Source: https://docs.sciml.ai/BifurcationKit/v0.4/PALC/

The getting-started docs show the standard API shape: construct a bifurcation
problem, call `continuation(prob, PALC(), ContinuationPar(...))`, inspect
`br.specialpoint`, and use `bifurcationdiagram` or branch-switching calls for
connected branch exploration. Source:
https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/gettingstarted/

### AUTO-07p

AUTO is official software for continuation and bifurcation problems in ODEs.
Its home page says it handles algebraic systems and ODEs with initial
conditions, boundary conditions, and integral constraints, and includes HOMCONT
for homoclinic bifurcation analysis. Source: https://auto-07p.github.io/

The AUTO manual states that its main algorithms target continuation of ODE
systems subject to boundary conditions and integral constraints. For BVPs, it
can compute solution curves with general nonlinear boundary and integral
conditions, allow boundary conditions depending on both endpoints and
parameters, determine folds and branch points along BVP solution families,
switch branches at branch points, and compute curves of folds and branch
points. Source, manual LaTeX in official repository:
https://raw.githubusercontent.com/auto-07p/auto-07p/master/doc/auto.tex

AUTO uses pseudo-arclength continuation. Its manual defines `DS`, `DSMIN`,
`DSMAX`, and `IADS` for initial, minimum, maximum, and adaptive
pseudo-arclength step control, including retrying failed Newton/Chord steps with
half the step size. Source:
https://raw.githubusercontent.com/auto-07p/auto-07p/master/doc/auto.tex

AUTO's official repository reports release `0.9.3` as the latest GitHub
release in the current public release list, with Fortran/Python tooling and
manual updates. Source: https://github.com/auto-07p/auto-07p/releases

### Ascher-Mattheij-Russell, COLSYS, and COLNEW

Ascher, Mattheij, and Russell's book remains the standard reference for BVP ODE
conditioning, shooting, multiple shooting, finite differences, collocation,
Newton methods, mesh selection, and singular perturbation behavior. Source:
SIAM Classics DOI https://doi.org/10.1137/1.9781611971231

The COLSYS paper by Ascher, Christiansen, and Russell describes spline
collocation for mixed-order ODE BVP systems, including error estimation,
adaptive mesh selection, B-spline basis evaluation, linear algebra, and nonlinear
solves. Source DOI:
https://doi.org/10.1090/S0025-5718-1979-0521281-7

Bader and Ascher's COLNEW paper introduced a new basis implementation for
mixed-order BVP ODE solvers, replacing the B-spline basis used in COLSYS and
improving algebraic solver performance. Source DOI:
https://doi.org/10.1137/0908047

The public COLNEW Fortran source says COLNEW solves multipoint BVPs for
mixed-order ODE systems, with side conditions at ordered points, collocation at
Gaussian points, a Runge-Kutta-monomial representation, adaptive mesh
selection/error-estimation subroutines, and simple continuation by reusing
previous solution coefficients. Source:
https://www.cs.ubc.ca/~ascher/colnew.f

### MATLAB, SciPy, and Residual-Control BVP Solvers

MATLAB's `bvp4c` documentation says it implements a three-stage Lobatto IIIA
collocation formula, provides a uniformly fourth-order `C1` continuous
solution, and uses residual-based mesh selection and error control. It also
supports singular BVPs through a constant singular term matrix and smoothness
condition at the singular endpoint. Sources:
https://www.mathworks.com/help/matlab/ref/bvp4c.html and
https://www.mathworks.com/help/matlab/math/boundary-value-problems.html

Kierzenka and Shampine's bvp4c paper is "A BVP solver based on residual control
and the MATLAB PSE", ACM TOMS 27(3), 299-316. Source DOI:
https://doi.org/10.1145/502800.502801

SciPy's `solve_bvp` solves first-order BVP systems with optional unknown
parameters and an optional singular term `S*y/(x-a)`. Its docs state that it
uses a fourth-order collocation algorithm with residual control and a damped
Newton method for the collocation system. Source:
https://docs.scipy.org/doc/scipy-1.16.1/reference/generated/scipy.integrate.solve_bvp.html

### Chebfun and Spectral Collocation

Chebfun represents smooth functions adaptively by Chebyshev polynomial
interpolation and aims for about machine precision when possible. Its guide
also notes support for certain infinite intervals, with a warning that
operations involving infinities are not always as accurate or robust as finite
interval operations. Source:
https://www.chebfun.org/docs/guide/guide01.html

Chebfun `chebop` supports ODE BVPs via spectral collocation. The guide says
linear two-point BVPs can be solved by backslash, Chebfun automatically chooses
discretizations for high accuracy, and version 5 provides rectangular
collocation and ultraspherical spectral methods. Source:
https://www.chebfun.org/docs/guide/guide07.html

For nonlinear BVPs, Chebfun uses the same `chebop` syntax and refines spectral
grids until convergence; the nonlinear solve is described as a Newton, sometimes
damped Newton, iteration in function space. Source:
https://www.chebfun.org/docs/guide/guide10.html

Chebfun's boundary-layer example states that global Chebyshev collocation can
become inefficient for very rapid boundary layers and that Chebfun does not have
adaptive grid refinement in the usual sense. It recommends user-inserted
breakpoints as an a priori way to use separate Chebyshev grids on subintervals.
Source: https://www.chebfun.org/examples/ode-linear/Breakpoints.html

### Singular and Semi-Infinite Domains

For endpoint singularities of the form `y' = S*y/x + f(x,y)`, both MATLAB
`bvp4c`/`bvp5c` and SciPy `solve_bvp` expose the singular matrix formulation and
require consistency with the smoothness condition at the singular endpoint.
Sources: MATLAB singular BVP docs
https://www.mathworks.com/help/matlab/math/solve-bvp-with-singular-term.html
and SciPy `solve_bvp`
https://docs.scipy.org/doc/scipy-1.16.1/reference/generated/scipy.integrate.solve_bvp.html

For infinite or long intervals, Markowich and Ringhofer studied truncating the
infinite interval at a large finite endpoint, imposing asymptotic boundary
conditions, and using stable symmetric collocation with exponentially graded
meshes when the infinite-domain solution decays exponentially. Source DOI:
https://doi.org/10.1090/S0025-5718-1983-0679437-X

## Recommendations and Inferences for This Package

These are recommendations inferred from the sources above and the current
repository state. They are not claims made by any single cited source.

1. Make adaptive collocation the default production solver.

   The current single-shooting Newton method is useful as a diagnostic and as a
   regression baseline, but it couples convergence to a fragile IVP shooting
   map. A `BoundaryValueDiffEq.jl` collocation formulation should become the
   default for the inner similarity BVP. Start with `MIRK4` or `MIRK5`, expose
   `MIRK6` and `Ascher*` options, and preserve the current `Rodas5P` shooting
   path as `method = :shooting` or an internal comparison mode.

2. Formulate the unknown tip data as BVP parameters, not just shooting guesses.

   The current free variables `(xi0, S0, Spp0)` map naturally to BVP unknown
   parameters. The collocation residual should include tip regularity,
   far-field slope/velocity/curvature residuals, and any paper-derived matching
   constraints. This lets the nonlinear solve see the whole BVP instead of only
   the end-point shooting map.

3. Treat the semi-infinite domain explicitly.

   Add a domain policy object, for example:

   - `TruncatedDomain(xi_max, farfield = ...)`
   - `MappedDomain(map = :rational | :tanh, endpoint = Inf)`
   - `ContinuationDomain(xi_max_sequence = ...)`

   The first production path should be truncation plus asymptotic far-field
   conditions, with continuation in `xi_max`. Compactified spectral or
   Chebfun-style references should wait until the correct far-field asymptotics
   are settled.

4. Add homotopy before branch analysis.

   A simple continuation ladder should solve easier problems first and reuse the
   solution as the next initial guess. Useful homotopy axes for this repository
   are likely:

   - capillary-wave/axial-curvature term off to on;
   - smaller `xi_max` to larger `xi_max`;
   - relaxed far-field residuals to strict residuals;
   - larger cone angle or regularized problem to target small angle;
   - outer/inner matching point moved outward.

5. Use BifurcationKit for folds, branch switching, and parameter diagrams.

   Natural parameter continuation is enough for monotone homotopy paths, but it
   is not the right abstraction for folds. Once the collocation residual is
   available as `F(u, p) = 0`, add an optional BifurcationKit extension that can
   run PALC, detect folds/branch points, and continue fold curves in two
   parameters. Do not initially implement custom pseudo-arclength continuation
   inside `SlenderConeRecoil.jl`.

6. Use AUTO-07p and COLNEW as validation references, not runtime dependencies.

   AUTO is the mature external reference for continuation and BVP branch/fold
   handling. COLNEW/COLSYS are the historical reference for adaptive Gaussian
   collocation on mixed-order BVPs. For this Julia package, they should inform
   benchmarks, documentation, and validation examples, but the public API should
   first use SciML-native solvers.

7. Use Chebfun for high-accuracy smooth reference profiles.

   Chebfun is best suited here as an independent MATLAB-side reference for
   smooth finite-domain or compactified BVPs, spectral convergence checks, and
   sanity checks of asymptotic ODE reductions. It should not be the default
   production model for the capillary-wave inner BVP until layer locations and
   singular behavior are controlled.

8. Expose solver diagnostics as package data, not console output.

   Any new solve result should record:

   - solver family and algorithm;
   - nonlinear convergence status and residual norm;
   - mesh, number of intervals/nodes, and final defect or error estimate;
   - final far-field residual vector;
   - continuation path, parameter sequence, and failed steps;
   - branch/fold metadata when BifurcationKit is used.

9. Keep dependency load modular.

   `BoundaryValueDiffEq.jl` and `BifurcationKit.jl` will add compile and
   precompile cost. Prefer package extensions or optional extras if possible:
   the core package can retain the current lightweight solve path, while
   `SlenderConeRecoilBoundaryValueDiffEqExt` and
   `SlenderConeRecoilBifurcationKitExt` provide production BVP and continuation
   features.

## Julia Operational Caveats

The repository instructions and README both warn against concurrent Julia
package operations, precompilation, or test jobs in the same project
environment. This matters more if `BoundaryValueDiffEq.jl`, `NonlinearSolve.jl`,
`LinearSolve.jl`, and `BifurcationKit.jl` are added, because first-use
precompilation can be substantial.

Recommended workflow:

- run only one Julia package/test/precompile job at a time in this project;
- for parallel experiments, use separate temporary depots and avoid manifest
  writes;
- keep BVP benchmark scripts read-only with respect to `Project.toml` and
  `Manifest.toml`;
- separate package-API tests from slow continuation and branch-diagram tests.

## Candidate Upgrade Path

1. Define a source-fidelity BVP residual for the current inner similarity
   equations, including parameterized tip and far-field conditions.
2. Add `BoundaryValueDiffEq.jl` as an optional or direct dependency after a
   dependency policy decision.
3. Implement `solve_inner_bvp_collocation(; alg = MIRK4(), ...)` and compare it
   against the current shooting result on the existing regression cases.
4. Add homotopy helpers that reuse previous solutions as initial guesses.
5. Add benchmark problems before relying on the new solver: Bratu, a singular
   Emden-type problem, a boundary-layer problem, and a semi-infinite
   Falkner-Skan-like truncation test.
6. Add a BifurcationKit extension around the discretized BVP residual for PALC,
   folds, and branch points.
7. Add external validation notes or scripts for AUTO/COLNEW/Chebfun, but keep
   them outside the package runtime path.

## Access and Source Gaps

No Julia commands were needed for this report. I used public official package
documentation, official source repositories, publisher DOI landing pages, and
local repository files. I did not use TIB/VPN or acquire paywalled full-text
PDFs for algorithm papers.

The local `docs/papers/README.md` inventory is focused on capillary-flow source
papers and does not currently list local PDFs for the BVP algorithm references
surveyed here. For future source-fidelity work, algorithm PDFs worth acquiring
or recording in a separate numerical-methods manifest include:

- Ascher, Mattheij, Russell, *Numerical Solution of Boundary Value Problems for
  Ordinary Differential Equations*, DOI `10.1137/1.9781611971231`;
- Ascher, Christiansen, Russell (1979) COLSYS paper, DOI
  `10.1090/S0025-5718-1979-0521281-7`;
- Bader and Ascher (1987) COLNEW paper, DOI `10.1137/0908047`;
- Kierzenka and Shampine (2001) bvp4c paper, DOI `10.1145/502800.502801`;
- Boisvert, Muir, Spiteri (2013), DOI `10.1145/2427023.2427028`;
- Markowich and Ringhofer (1983), DOI
  `10.1090/S0025-5718-1983-0679437-X`.

## Sources

- SciML BVP solvers:
  https://docs.sciml.ai/DiffEqDocs/dev/solvers/bvp_solve/
- BoundaryValueDiffEq shooting:
  https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/shooting/
- BoundaryValueDiffEq MIRK:
  https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/mirk/
- BoundaryValueDiffEq FIRK:
  https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/firk/
- BoundaryValueDiffEq Ascher:
  https://docs.sciml.ai/BoundaryValueDiffEq/dev/solvers/ascher/
- BoundaryValueDiffEq error control:
  https://docs.sciml.ai/BoundaryValueDiffEq/dev/basics/error_control/
- BoundaryValueDiffEq continuation tutorial:
  https://docs.sciml.ai/BoundaryValueDiffEq/dev/tutorials/continuation/
- BifurcationKit capabilities:
  https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/capabilities/
- BifurcationKit getting started:
  https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/gettingstarted/
- BifurcationKit PALC:
  https://docs.sciml.ai/BifurcationKit/v0.4/PALC/
- AUTO homepage: https://auto-07p.github.io/
- AUTO manual source:
  https://raw.githubusercontent.com/auto-07p/auto-07p/master/doc/auto.tex
- AUTO releases: https://github.com/auto-07p/auto-07p/releases
- COLNEW source: https://www.cs.ubc.ca/~ascher/colnew.f
- Ascher, Mattheij, Russell SIAM book:
  https://doi.org/10.1137/1.9781611971231
- Ascher, Christiansen, Russell COLSYS paper:
  https://doi.org/10.1090/S0025-5718-1979-0521281-7
- Bader and Ascher COLNEW paper: https://doi.org/10.1137/0908047
- Enright and Muir MIRK defect control:
  https://doi.org/10.1137/S1064827593251496
- Boisvert, Muir, Spiteri global error/defect control:
  https://doi.org/10.1145/2427023.2427028
- MATLAB `bvp4c`:
  https://www.mathworks.com/help/matlab/ref/bvp4c.html
- MATLAB BVP overview:
  https://www.mathworks.com/help/matlab/math/boundary-value-problems.html
- MATLAB singular BVP example:
  https://www.mathworks.com/help/matlab/math/solve-bvp-with-singular-term.html
- Kierzenka and Shampine bvp4c paper:
  https://doi.org/10.1145/502800.502801
- SciPy `solve_bvp`:
  https://docs.scipy.org/doc/scipy-1.16.1/reference/generated/scipy.integrate.solve_bvp.html
- Chebfun guide, getting started:
  https://www.chebfun.org/docs/guide/guide01.html
- Chebfun guide, linear differential operators:
  https://www.chebfun.org/docs/guide/guide07.html
- Chebfun guide, nonlinear ODEs:
  https://www.chebfun.org/docs/guide/guide10.html
- Chebfun breakpoint/layer example:
  https://www.chebfun.org/examples/ode-linear/Breakpoints.html
- Markowich and Ringhofer long-interval collocation:
  https://doi.org/10.1090/S0025-5718-1983-0679437-X
