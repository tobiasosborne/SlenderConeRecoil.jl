# Similarity And Asymptotics Algorithms

Date: 2026-06-01

Bead: `scr-x9r`

Scope: algorithms and workflows for deriving and computing similarity
reductions and matched asymptotics, with emphasis on features that could guide
a future general `SimilarityMethods.jl` architecture.

## Local Paper And Access Notes

Checked against `docs/papers/README.md`, not against ignored PDF contents.

| Source | Local status in manifest | Access note |
| --- | --- | --- |
| Decent and King, "Surface-tension-driven flow in a slender cone", DOI [`10.1093/imamat/hxm043`](https://doi.org/10.1093/imamat/hxm043) | Missing | Primary cone-recoil source remains an acquisition gap. Metadata confirmed by the University of Birmingham page. |
| Keller and Miksis, "Surface Tension Driven Flows", DOI [`10.1137/0143018`](https://doi.org/10.1137/0143018) | Present locally, ignored | Useful source for capillary `t^(2/3)` scaling benchmarks. |
| Eggers, "Nonlinear dynamics and breakup of free-surface flows", DOI [`10.1103/RevModPhys.69.865`](https://doi.org/10.1103/RevModPhys.69.865) | Present locally, ignored | Useful review source for similarity, singularity, and free-surface validation context. |
| Lo, "Asymptotic matching by the symbolic manipulator MACSYMA", DOI [`10.1016/0021-9991(85)90059-2`](https://doi.org/10.1016/0021-9991(85)90059-2) | Not listed | Publisher access likely needed for detailed algorithm extraction. |
| Fraenkel, "On the method of matched asymptotic expansions", DOI [`10.1017/S0305004100044212`](https://doi.org/10.1017/S0305004100044212) | Not listed | Cambridge PDF preview was accessible through web search; add to manifest only if a licensed local copy is acquired. |

## Source-Backed Facts

### Scaling And Dimensional Analysis

- Dimensional analysis can be implemented as linear algebra on a dimension
  matrix: dimensionless products correspond to a nullspace basis. This is the
  organizing algorithm in the Buckingham-Pi linear-algebra formulation
  described by Curtis, Logan, and Parker, DOI
  [`10.1016/0024-3795(82)90229-4`](https://doi.org/10.1016/0024-3795(82)90229-4),
  and in package implementations such as
  [`UnitfulBuckinghamPi.jl`](https://rmsrosa.github.io/blog/2021/05/unitfulbuckinghampi/)
  and BuckinghamPy, DOI
  [`10.1016/j.softx.2021.100851`](https://doi.org/10.1016/j.softx.2021.100851).
- Practical Pi-group generation needs more than a floating nullspace. The
  `UnitfulBuckinghamPi.jl` writeup emphasizes preserving integer or rational
  exponents and using exact linear algebra where possible; BuckinghamPy
  generates multiple sets of dimensionless numbers and LaTeX/symbolic output.
- Butterfield argues that Buckingham's theorem gives necessary, not always
  sufficient, conditions for a successful dimensional analysis; the paper adds
  physical-knowledge augmentation to choose useful primitive dimensions, DOI
  [`10.1243/0954406011524748`](https://doi.org/10.1243/0954406011524748).
  This matters for package design because blindly returning any nullspace basis
  is not the same as returning a physically useful scaling.
- Keller and Miksis show surface-tension-driven potential flows where velocity
  scales like `(sigma / (rho t))^(1/3)`, implying the length scale
  `(sigma t^2 / rho)^(1/3)`, and compute free-surface shapes after reduction to
  a self-similar integral/differential problem, DOI
  [`10.1137/0143018`](https://doi.org/10.1137/0143018). This is a concrete
  regression target for any similarity-scaling module.

### Symmetry-Based Similarity Reduction

- Lie symmetry software is mature enough to use as design precedent. Maple's
  [`SimilaritySolutions`](https://www.maplesoft.com/support/help/Maple/view.aspx?path=PDEtools%2FSimilaritySolutions)
  command computes similarity solutions from a one-dimensional admitted
  symmetry group, reducing the number of independent variables by one and then
  attempting to solve the reduced problem.
- The Maple-based GeM package automatically generates determining equations for
  local, contact, higher-order, and approximate symmetry analysis, conservation
  law multipliers, and adjoint symmetries for ODE/PDE systems, DOI
  [`10.1016/j.cpc.2006.08.001`](https://doi.org/10.1016/j.cpc.2006.08.001).
- ReLie is an open-source REDUCE package for Lie group analysis of differential
  equations, covering point, conditional, contact, variational, approximate,
  and equivalence symmetries, DOI
  [`10.3390/sym13101826`](https://doi.org/10.3390/sym13101826).
- Comparative work on Lie-symmetry packages reports that systems can differ
  materially in completeness of detected symmetries; DESOLV was reported as
  especially successful for complete Lie point symmetries of PDE systems, DOI
  [`10.1016/S0010-4655(03)00348-5`](https://doi.org/10.1016/S0010-4655(03)00348-5).
  A Julia package should therefore expose verification artifacts, not just a
  final change of variables.

### Dominant-Balance And Scale Discovery

- Dominant-balance automation is best developed where terms are polynomial or
  monomial. Tropical equilibration methods formulate balance candidates as
  equality of at least two dominant monomials, typically of opposite signs, in
  each dynamic equation; this gives scaling exponents and slow-fast
  decompositions for polynomial ODE models. See the constraint-solving approach,
  DOI [`10.1186/s13015-014-0024-2`](https://doi.org/10.1186/s13015-014-0024-2),
  and the tropical-geometry review, DOI
  [`10.3389/fgene.2012.00131`](https://doi.org/10.3389/fgene.2012.00131).
- Recent tropical reduction work for chemical reaction networks explicitly
  describes symbolic algorithms to reshape and rescale networks so geometric
  singular perturbation theory can be applied, test applicability, and reduce
  models with approximate conservation laws, DOI
  [`10.1137/22M1543963`](https://doi.org/10.1137/22M1543963).
- This tropical literature is not a turnkey general PDE similarity reducer, but
  it gives a useful pattern for an implementation: enumerate candidate balances,
  reject inconsistent balances, and retain scaling metadata for later
  validation.
- Numerical renormalization group algorithms compute asymptotically
  self-similar PDE dynamics by repeatedly rescaling numerical PDE solutions and
  extracting profiles/exponents as fixed points, DOI
  [`10.1137/18M120004X`](https://doi.org/10.1137/18M120004X). This is a
  complementary route when the scaling is not obvious from symbolic balances.

### Symbolic And Formal Asymptotic Expansion

- Wolfram Language's
  [`AsymptoticDSolveValue`](https://reference.wolfram.com/language/ref/AsymptoticDSolveValue)
  computes asymptotic approximations for scalar ODEs, ODE systems, and
  parameter-dependent differential equations, with order control. The
  documentation lists Taylor, Frobenius, WKB, and boundary-layer-style
  approximations as covered categories.
- `Symbolics.jl` and `ModelingToolkit.jl` provide a Julia-native symbolic
  substrate rather than a matched-asymptotics workflow out of the box.
  `Symbolics.jl` is designed as an extensible symbolic system based on multiple
  dispatch and generic term interfaces, DOI
  [`10.1145/3511528.3511535`](https://doi.org/10.1145/3511528.3511535).
  `ModelingToolkit.jl` is a graph-transformation system for equation-based
  modeling that applies symbolic transformations and generates efficient
  numerical implementations; see DOI
  [`10.48550/arXiv.2103.05244`](https://doi.org/10.48550/arXiv.2103.05244)
  and the
  [stable documentation](https://docs.sciml.ai/ModelingToolkit/stable/).
- Lo's MACSYMA paper demonstrates that symbolic manipulation can automate the
  tedious higher-order algebra in singular perturbation matching, including
  expansions with powers and logarithms, DOI
  [`10.1016/0021-9991(85)90059-2`](https://doi.org/10.1016/0021-9991(85)90059-2).
  This is historical, but it directly implies package features: ordered
  asymptotic scales, logarithmic terms, common-part extraction, and symbolic
  bookkeeping of matching constants.

### Matched-Asymptotic Matching And Overlap Validation

- Fraenkel's review distinguishes overlap-based matching from simple use of
  Van Dyke's matching principle and gives counterexamples where the easy
  principle can be incorrect, DOI
  [`10.1017/S0305004100044212`](https://doi.org/10.1017/S0305004100044212).
  A reliable package should therefore implement explicit overlap diagnostics
  instead of assuming that equal truncated expansions at a formal limit are
  sufficient.
- The matching and composite-expansion literature treats the "common part" of
  inner and outer expansions as an object to compute, not as a plotting
  convention. This supports storing region, variable transform, expansion
  scale, and common-part metadata separately so additive composites can be
  audited term by term.

### Numerical Solvers For Similarity ODEs

- General-purpose adaptive collocation for BVPs is a proven architecture. The
  COLSYS paper implements spline collocation for mixed-order ODE BVPs with
  error estimation, adaptive mesh selection, B-spline basis evaluation, and
  nonlinear solves, DOI
  [`10.1090/S0025-5718-1979-0521281-7`](https://doi.org/10.1090/S0025-5718-1979-0521281-7).
  COLNEW is a later FORTRAN code available through Netlib and described in the
  SIAM BVP monograph appendix, DOI
  [`10.1137/1.9781611971231.appb`](https://doi.org/10.1137/1.9781611971231.appb).
- MATLAB's `bvp4c` lineage is important because it frames BVP quality in terms
  of residual control and error estimation; Kierzenka and Shampine's paper is
  available as "A BVP solver based on residual control and the MATLAB PSE", DOI
  [`10.1145/502800.502801`](https://doi.org/10.1145/502800.502801).
- Julia already has a native SciML BVP layer:
  [`BoundaryValueDiffEq.jl`](https://docs.sciml.ai/DiffEqDocs/dev/api/boundaryvaluediffeq/)
  provides shooting methods and MIRK/FIRK collocation methods under the SciML
  interface, with problem types for several boundary-value formulations.
- Spectral collocation is a strong fit for smooth similarity ODEs and far-field
  tail problems. Chebfun's rectangular spectral collocation method is the basis
  for ODE solving in Chebfun, DOI
  [`10.1093/imanum/dru062`](https://doi.org/10.1093/imanum/dru062). Chebfun
  also has literature and examples for nonlinear/singular one-dimensional BVPs,
  DOI [`10.3390/computation10070116`](https://doi.org/10.3390/computation10070116).
- `ApproxFun.jl` is the closest Julia analogue for operator-based spectral
  computation: its documentation shows boundary-value ODEs solved by composing
  differential and boundary operators, and it uses adaptive coefficient
  representations for functions. See
  [ApproxFun ODE docs](https://juliaapproximation.github.io/ApproxFun.jl/latest/generated/ODE/)
  and [ApproxFun equation docs](https://juliaapproximation.github.io/ApproxFun.jl/latest/usage/equations/).
- Taylor/series methods are useful around regular or singular endpoints where
  boundary data are better represented as local expansions than as point values.
  `TaylorSeries.jl` provides Julia Taylor expansions in one or more variables,
  DOI [`10.21105/joss.01043`](https://doi.org/10.21105/joss.01043). High-order
  Taylor ODE software with adaptive order and step size is a mature method for
  non-stiff IVPs and high precision, DOI
  [`10.1080/10586458.2005.10128904`](https://doi.org/10.1080/10586458.2005.10128904).

## Inferences And Recommendations

These are design conclusions drawn from the sources above, not direct source
claims.

### Derivation Workflow

- Treat similarity derivation as a pipeline with auditable intermediate
  products:
  1. dimensional variables and dimensions;
  2. candidate Pi groups and scale choices;
  3. scaling or Lie generators;
  4. similarity ansatz and transformed derivatives;
  5. reduced ODE/integro-ODE system;
  6. local endpoint expansions;
  7. BVP residual and boundary-condition residual.
- Prefer exact rational exponent arithmetic for scale derivations. Floating
  exponents are acceptable in numerical diagnostics, but not in canonical
  derivation output.
- Provide multiple generators for a scaling problem: nullspace Pi groups,
  user-specified repeating variables, Lie scaling generators, and
  dominant-balance candidates. Package APIs should make it clear which path
  produced a result.
- Include a "source-fidelity ledger" for derivations: every nondimensional
  variable, small parameter, far-field condition, and matching condition should
  be traceable to either a source equation, a symbolic transform, or a numerical
  fit.

### Dominant Balance

- Implement dominant-balance automation first for algebraic/polynomialized
  equation systems. Convert rational expressions to numerator/denominator form
  or explicit monomial sums where safe; reject expressions outside the
  supported grammar unless a user supplies a balance.
- Use the tropical approach as the model for enumeration:
  collect exponents of each term under a candidate scale, solve equal-order
  constraints, check sign/physics consistency, and then substitute the balance
  back into the full equations to verify that neglected terms are smaller in an
  explicit asymptotic limit.
- Store rejected balances and rejection reasons. In asymptotics, a wrong
  balance is often as useful diagnostically as the accepted one.

### Matched-Asymptotic Workflow

- A matched-asymptotics toolkit should model regions explicitly:
  `Inner`, `Outer`, and optional `Intermediate` regions with their own
  variables, validity assumptions, truncation order, and scale hierarchy.
- Matching should be more than "choose a splice point". Useful diagnostics:
  common-part extraction, overlap-window discovery, mismatch norms across a
  moving intermediate range, convergence of matching constants as the window
  moves, and sensitivity to truncation order.
- Additive composites should be represented as data structures containing the
  inner expression, outer expression, common part, and assumptions. Plotting a
  composite should be a view of that object, not the object itself.
- For numerical matching, expose at least two independent checks:
  1. asymptotic-tail fit of each region in the overlap variables;
  2. direct residual evaluation of the composite in the original governing
     equation after transforming back to physical variables.

### Numerical Similarity Solves

- Support at least three solver paths for similarity ODEs:
  shooting/continuation, adaptive collocation, and spectral/operator solve.
  The goal is cross-validation, not just solver optionality.
- Similarity BVPs commonly have unknown parameters such as front position,
  amplitude, or eigenvalue. Represent these as solve variables in a nonlinear
  residual problem rather than as ad hoc closures around an IVP solve.
- Far-field boundary conditions should support asymptotic conditions, not only
  finite endpoint values. A package should be able to encode "tail approaches
  `a*x + b + c*x^(-p)`" or "oscillatory tail with decaying envelope" and test
  finite-domain truncation sensitivity.
- For capillary/free-boundary problems, validation should include both the
  reduced similarity ODE residual and a rescaled time-dependent PDE check.
  Numerical RG/rescaling workflows can be used as independent confirmation that
  the computed profile is the attractor of the PDE, not just a solution of a
  guessed ODE.

### Julia Architecture

- Build on SciML abstractions where they match the problem:
  `NonlinearProblem` for residuals and unknown parameters,
  `BVProblem` for collocation/shooting, `ModelingToolkit` or `Symbolics` for
  transformations, and `SciMLBase` return codes/diagnostics.
- Keep asymptotic algebra decoupled from a single CAS. A small internal IR for
  scales, regions, orders, and terms can have adapters to `Symbolics.jl`,
  `SymPy`, or external CAS output.
- Reuse package-specific strengths:
  `Unitful`/`UnitfulBuckinghamPi`-style dimension metadata,
  `HomotopyContinuation.jl` for algebraic balance systems when polynomial,
  `NonlinearSolve.jl` for shooting/collocation residuals,
  `BoundaryValueDiffEq.jl` for BVPs,
  `ApproxFun.jl` for spectral BVP prototypes,
  `TaylorSeries.jl` for endpoint series,
  and `BifurcationKit.jl` for continuation in small parameters or similarity
  exponents.

## Candidate Package Capabilities

Concrete capabilities to consider for a future `SimilarityMethods.jl`:

1. Dimension-aware variable registry with exact dimension matrices, rational
   nullspace Pi groups, user-selectable repeating variables, and LaTeX/Julia
   output for the resulting scalings.
2. Scaling-generator API that can derive similarity ansatzes from dimensional
   analysis, user-specified scale groups, Lie-style infinitesimals, or
   dominant-balance enumeration.
3. Symbolic derivative transformer for common PDE ansatzes, including multiple
   dependent variables, moving fronts, unknown exponents, and small parameters.
4. Dominant-balance engine for monomial/rational equation systems with
   candidate enumeration, consistency checks, rejected-balance logs, and
   optional homotopy solving of algebraic balance equations.
5. Asymptotic-scale IR supporting powers, fractional powers, logarithms,
   exponentials, nested scales, order comparison, truncation, substitution, and
   common-part extraction.
6. Region objects for inner/outer/intermediate expansions with explicit
   variables, transformations, validity assumptions, boundary conditions, and
   source references.
7. Matched-asymptotics workspace that stores expansions, computes matching
   constants, constructs additive composites, and reports overlap diagnostics
   rather than only returning a curve.
8. BVP model builder for similarity ODEs with unknown parameters, singular
   endpoint series, far-field asymptotic boundary conditions, and automatic
   residual scaling.
9. Solver adapters for shooting, adaptive collocation, and spectral/operator
   discretization, all returning a common diagnostics object with residuals,
   mesh/domain sensitivity, and boundary-condition errors.
10. Numerical RG/rescaling workflow for PDE simulations that estimates
    similarity exponents and profiles from time-dependent data, usable as an
    independent validation path.
11. Overlap-region validator that searches intermediate windows, fits common
    parts, tracks mismatch norms as the window moves, and flags matches that
    depend strongly on arbitrary splice choices.
12. Source-fidelity and reproducibility layer: every derivation step can carry
    source equation IDs, assumptions, commands, solver tolerances, and generated
    figure metadata.
13. Benchmark suite seeded with canonical cases: Keller-Miksis capillary
    scaling, Blasius-type boundary layers, simple singular perturbation BVPs,
    thin-film/self-similar rupture examples, and the Decent-King slender-cone
    problem once the correct primary PDF is available.
