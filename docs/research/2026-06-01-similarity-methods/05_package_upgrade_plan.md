# Package Upgrade Plan

Date: 2026-06-01

Assigned bead: `scr-ts1`

Scope: synthesis plan for turning `SlenderConeRecoil.jl` into a
best-in-class repository for computational similarity methods, matched
asymptotics, and capillary free-surface benchmarks.

Inputs read:

- `README.md`
- `HANDOFF.md`
- `AGENTS.md`
- `reviews/2026-06-01-full-audit/remediation_plan.md`
- `docs/research/2026-06-01-similarity-methods/README.md`
- `docs/research/2026-06-01-similarity-methods/01_capillary_recoil_literature.md`
- `docs/research/2026-06-01-similarity-methods/02_similarity_and_asymptotics_algorithms.md`
- `docs/research/2026-06-01-similarity-methods/03_bvp_continuation_algorithms.md`
- `docs/research/2026-06-01-similarity-methods/04_pde_verification_and_package_architecture.md`

No Julia commands are required for this synthesis plan. The issue identifiers
below were converted into Beads under the `scr-8l4` execution umbrella after
review.

## Execution Principles

- Finish the review-remediation trust blockers before starting broad package
  upgrades: correct source metadata, package-loading tests, fail-loud solver
  diagnostics, validated matching, and fast/slow gates.
- Treat Decent and King (2008), DOI `10.1093/imamat/hxm043`, as the canonical
  cone-recoil source. Decent and King (2001) is useful provenance, not the
  quantitative benchmark authority.
- Keep source-fidelity prerequisites separate from algorithm upgrades. Do not
  bless current numerical values as Decent-King benchmarks until the 2008 PDF
  has been acquired and the equation ledger has been extracted.
- Keep the package core small. Add heavier SciML, spectral, benchmark, and
  artifact tooling only after a dependency/precompile policy is explicit.
- Every code-changing implementation issue should run the fast package gate.
  Solver, PDE, composite, regression, and benchmark changes should also run
  the slow gate. Do not run concurrent Julia package, precompile, or test jobs
  against the same environment.

## Milestones

| Milestone | Priority | Dependencies | Rationale | Suggested issues |
| --- | --- | --- | --- | --- |
| M0: Source-fidelity prerequisites | P0 | Review remediation chain complete | The repository cannot become authoritative while the primary cone source is missing and current equations/BCs/matching are provisional. | `UP-SF1` through `UP-SF5` |
| M1: API and provenance foundation | P1 | M0 source ledger at least drafted; P0 diagnostics from remediation complete | Stable problem/result/diagnostic objects make later solver swaps and benchmarks auditable. | `UP-API1` through `UP-API4` |
| M2: Production similarity BVP algorithms | P1 | M0 complete; M1 result/diagnostic types available | Adaptive collocation and continuation are the highest-leverage algorithm upgrades over fragile single shooting. | `UP-BVP1` through `UP-BVP5` |
| M3: Matched-asymptotic and outer-wave validation | P1 | M0 Decent-King extraction; M1 metadata; parts of M2 | Matching must be a validated asymptotic object, not only a plotting splice. | `UP-MA1` through `UP-MA4` |
| M4: PDE verification ladder | P1/P2 | Remediation PDE safety; M1 diagnostics | PDE comparison should verify operators, conservation, similarity collapse, and independent discretizations. | `UP-PDE1` through `UP-PDE5` |
| M5: Benchmarks, datasets, examples, and reproducibility | P2 | M0/M1; relevant solver/PDE pieces | Best-in-class means reproducible benchmark data, not just figures. | `UP-REP1` through `UP-REP5` |
| M6: General similarity-methods architecture | P3 | Cone benchmark stable and documented | General APIs are valuable only after the flagship cone case is trustworthy. | `UP-GEN1` through `UP-GEN5` |

## Decisions Delayed Until Decent-King 2008 PDF Is Acquired

- The exact nondimensional variables, small-aspect-ratio ordering, and
  equation numbering to use in source-fidelity tests.
- Whether the current primitive slender model in `src/slender.jl` and
  similarity reduction in `src/similarity.jl` match the paper formulation or
  are only a related reconstruction.
- The correct tip regularity conditions, unknown tip parameters, and far-field
  boundary conditions for the inner BVP.
- The correct outer expansion powers, matching constants, capillary-wave phase
  convention, envelope law, and finite-domain truncation strategy.
- Whether the current `xi0`, `S0`, wave train, outer hierarchy, and composite
  figures are source-backed enough to be retained as paper-reproduction
  outputs.
- Which Decent-King numerical values should become regression data, and what
  tolerances are justified by the source and by independent numerical solves.
- Whether collocation should target the current ODE system unchanged or a
  revised source-faithful system.
- Whether full curvature, axial-curvature truncation, or another asymptotic
  curvature representation is the paper-faithful production model.

## Must-Do Source-Fidelity Prerequisites

### `UP-SF1`: Acquire And Checksum Missing Primary Sources

Type/priority: task, P0.

Scope:

- Acquire local licensed copies of:
  `DecentKing2008_IMAJAM_73_1_37-68_hxm043.pdf`,
  `Billingham1999_JFM_397_45.pdf`, and
  `KellerKingTing1995_PoF_7_226.pdf`.
- Record access notes and SHA-256 checksums in `docs/papers/README.md`.
- Keep PDFs ignored by git; update only the tracked manifest and any fetch
  diagnostics intentionally.

Acceptance criteria:

- `docs/papers/README.md` marks the three target papers as present or records
  a precise access failure.
- Each acquired source has a checksum and DOI matching the research reports.
- Quarantined mislabeled artifacts remain quarantined or documented.
- No package source, tests, figures, or Beads state are changed in this issue.

Dependencies: review remediation source-fidelity fixes; none of the algorithm
upgrades should depend on Decent-King quantitative values before this is done.

Suggested validation gates:

- `git diff --check`
- `sha256sum docs/papers/<acquired-file>` for each local ignored PDF

Likely changed files/modules:

- `docs/papers/README.md`
- optionally `scripts/fetch_papers.mjs` if access tooling needs a documented
  fix

### `UP-SF2`: Extract A Decent-King Equation And Assumption Ledger

Type/priority: task, P0.

Scope:

- Create a source ledger mapping paper equations, assumptions, variables,
  asymptotic regions, boundary conditions, matching constants, and reported
  numerical values to package concepts.
- Separate 2008 IMA facts from 2001 IUTAM precursor facts.
- Mark any quantity that is inferred from plots or text rather than directly
  tabulated.

Acceptance criteria:

- A tracked ledger exists under `docs/research/2026-06-01-similarity-methods/`
  or `docs/roadmap/` with source equation IDs and page references.
- The ledger names all source quantities needed by inner, outer, composite, and
  PDE verification work.
- Unsupported current README/HANDOFF claims are listed for later documentation
  correction.

Dependencies: `UP-SF1`.

Suggested validation gates:

- `git diff --check`
- Manual cross-check against Decent-King 2008 PDF and local paper manifest

Likely changed files/modules:

- new ledger file under `docs/research/2026-06-01-similarity-methods/` or
  `docs/roadmap/`
- possibly `docs/method.md` if the orchestrator chooses to update narrative
  references in the same issue

### `UP-SF3`: Reconcile Slender Model And Similarity Reduction With The Ledger

Type/priority: bug, P0.

Scope:

- Compare `src/slender.jl`, `src/similarity.jl`, and `src/inner.jl` against the
  Decent-King ledger.
- Decide whether the current ODE state, curvature terms, signs,
  nondimensionalization, and regularity conditions are faithful, need
  correction, or should be labeled as a reconstructed model.
- Add source-backed unit tests for the symbolic and numerical residuals where
  the ledger gives enough detail.

Acceptance criteria:

- Each implemented governing equation has a source-ledger reference or an
  explicit "package reconstruction" label.
- Any equation/sign/condition mismatch is fixed or documented as provisional.
- Tests fail if the reconciled residual, scale, or boundary-condition formula
  regresses.
- README and `docs/method.md` no longer imply paper fidelity for unreconciled
  pieces.

Dependencies: `UP-SF2`; remediation issues for package-loading tests and
fail-loud solver diagnostics.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/slender.jl`
- `src/similarity.jl`
- `src/inner.jl`
- `src/solve_checks.jl`
- `test/test_bead3.jl`
- `test/test_bead4.jl`
- `test/test_bead5.jl`
- `test/test_numerical_regressions.jl`
- `README.md`
- `docs/method.md`

### `UP-SF4`: Reconcile Outer Expansion, Matching, And Composite Formulae

Type/priority: bug, P0.

Scope:

- Compare `src/outer.jl`, `src/outer_hierarchy.jl`, and `src/composite.jl`
  against Decent-King outer/asymptotic matching data.
- Establish the correct outer boundary data, common part, matching constants,
  and capillary-wave asymptotics.
- Retire or relabel any unsupported seeded outer solve or Laurent hierarchy.

Acceptance criteria:

- Outer and composite APIs have source-ledger references for equations and
  matching data, or are explicitly marked exploratory.
- Composite overlap subtraction matches the documented common part.
- Matching outside valid domains fails clearly.
- Regression tests cover empty overlap, bad grids, out-of-domain match points,
  and source-backed matching constants where available.

Dependencies: `UP-SF2`; remediation matching issue; `UP-SF3` for the inner
tail quantities.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/outer.jl`
- `src/outer_hierarchy.jl`
- `src/composite.jl`
- `test/test_bead6.jl`
- `test/test_bead7.jl`
- `test/test_outer_hierarchy.jl`
- `test/test_numerical_regressions.jl`
- `README.md`
- `docs/method.md`

### `UP-SF5`: Establish Source-Backed Cone Reference Data

Type/priority: task, P0.

Scope:

- Extract Decent-King numerical values and figure-derived quantities needed for
  regression tests.
- Store small reference data in a tracked text format when redistribution is
  allowed; otherwise store extraction metadata and local checksum references.
- Define cheap and reference tolerances separately.

Acceptance criteria:

- Reference data records include source DOI, page/equation/figure reference,
  extraction method, units/scaling, and tolerance rationale.
- Regression tests use source-backed values for inner parameters, far-field
  residuals, matching constants, and capillary-wave phase/envelope where the
  source permits.
- Current locally blessed values remain separate from source-backed values.

Dependencies: `UP-SF2`, `UP-SF3`, `UP-SF4`.

Suggested validation gates:

- `git diff --check`
- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `test/test_numerical_regressions.jl`
- possible new `data/reference/` or `test/reference/` files
- `docs/papers/README.md`
- `docs/method.md`

## Algorithm And Architecture Upgrades

### `UP-API1`: Define Public Problem, Parameter, And Result Types

Type/priority: feature, P1.

Scope:

- Introduce stable user-facing constructors for cone similarity, outer
  matching, composite profiles, and PDE verification problems.
- Define rich result structs with solver status, residuals, mesh/domain data,
  physical and nondimensional parameters, and source/provenance metadata.
- Preserve current functions as compatibility wrappers where practical.

Acceptance criteria:

- Public API is documented and loaded through `using SlenderConeRecoil`.
- Result objects expose diagnostics without relying on console output.
- Existing tests and scripts use exported API or clearly module-qualified
  internals.
- No broad package rename is included.

Dependencies: remediation package-load issue; remediation fail-loud diagnostics;
`UP-SF3` enough to know source metadata fields.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl` if solver
  result paths change

Likely changed files/modules:

- `src/SlenderConeRecoil.jl`
- possible new `src/problems.jl`
- possible new `src/results.jl`
- `src/inner.jl`
- `src/outer.jl`
- `src/composite.jl`
- `src/pde.jl`
- `scripts/figures.jl`
- `test/runtests.jl`
- `README.md`

### `UP-API2`: Add Source-Fidelity And Provenance Metadata Objects

Type/priority: feature, P1.

Scope:

- Add small data structures for source citations, equation IDs, assumptions,
  benchmark IDs, solver settings, and artifact hashes.
- Attach provenance to derivation outputs, solve results, benchmark data, and
  figure metadata.

Acceptance criteria:

- Every source-backed benchmark can report DOI, source equation/figure, local
  paper checksum status, and package version/commit when available.
- Figure and benchmark generation can emit provenance without duplicating
  ad-hoc metadata code.
- Metadata remains lightweight and does not force plotting or artifact
  dependencies into the core package.

Dependencies: `UP-SF2`; `UP-API1`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `julia --project=scripts scripts/figures.jl --metadata-only` when figure
  metadata paths change

Likely changed files/modules:

- possible new `src/provenance.jl`
- `src/SlenderConeRecoil.jl`
- `scripts/figures.jl`
- `figures/metadata.toml`
- `docs/papers/README.md`
- `test/runtests.jl`

### `UP-API3`: Introduce Diagnostics Utilities

Type/priority: feature, P1.

Scope:

- Centralize residual norms, boundary-condition residuals, far-field residuals,
  conservation diagnostics, collapse metrics, retcode checks, and mesh/domain
  summaries.
- Make diagnostics usable by solvers, tests, figures, and benchmarks.

Acceptance criteria:

- Inner, outer, composite, and PDE results expose comparable diagnostic fields.
- Tests assert diagnostic thresholds rather than only successful function
  return.
- Figure scripts include diagnostics in metadata and refuse failed solves.

Dependencies: remediation fail-loud diagnostics; `UP-API1`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/solve_checks.jl`
- possible new `src/diagnostics.jl`
- `src/inner.jl`
- `src/outer.jl`
- `src/composite.jl`
- `src/pde.jl`
- `scripts/figures.jl`
- `test/test_numerical_regressions.jl`

### `UP-API4`: Decide Dependency And Extension Policy

Type/priority: task, P1.

Scope:

- Decide whether `BoundaryValueDiffEq.jl`, `NonlinearSolve.jl`,
  `BifurcationKit.jl`, `ApproxFun.jl`, `MethodOfLines.jl`, `BenchmarkTools.jl`,
  and plotting tools belong in core deps, weak deps/extensions, test extras,
  script environments, or benchmark environments.
- Record precompile and test-gate implications.

Acceptance criteria:

- `Project.toml`, scripts environment, and planned benchmark environment have
  a documented dependency policy.
- Heavy optional algorithms have a path that does not make basic package load
  unnecessarily expensive.
- AGENTS/HANDOFF/README instructions remain consistent with non-concurrent
  Julia package operation guidance.

Dependencies: remediation repo hygiene; before adding any heavy dependency.

Suggested validation gates:

- `git diff --check`
- `julia --project test/runtests.jl` if `Project.toml` changes

Likely changed files/modules:

- `Project.toml`
- `README.md`
- `HANDOFF.md`
- `AGENTS.md`
- possible `benchmark/Project.toml`
- possible extension files under `ext/`

### `UP-BVP1`: Formulate A Source-Fidelity BVP Residual

Type/priority: feature, P1.

Scope:

- Represent the inner similarity BVP as a nonlinear residual with unknown tip
  parameters and source-backed boundary conditions.
- Keep the current shooting solve as a baseline implementation during the
  transition.

Acceptance criteria:

- Residual evaluation works independently of a particular solver.
- Unknown parameters such as tip location/radius/curvature are explicit solve
  variables.
- Boundary and far-field residual components are individually inspectable.
- Tests cover residual size on known solutions and failure on inconsistent
  parameters.

Dependencies: `UP-SF3`; `UP-API1`; `UP-API3`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/inner.jl`
- possible new `src/bvp_residuals.jl`
- `src/solve_checks.jl`
- `test/test_bead5.jl`
- `test/test_numerical_regressions.jl`

### `UP-BVP2`: Add Adaptive Collocation Inner Solver

Type/priority: feature, P1.

Scope:

- Implement a SciML `BVProblem` path for the inner similarity BVP, starting
  with `BoundaryValueDiffEq` MIRK methods if the dependency policy allows.
- Expose collocation as the production solve path only after diagnostics and
  source-backed residuals agree.

Acceptance criteria:

- `solve_inner(...; method = :collocation)` or equivalent returns the common
  result type with mesh/error/residual diagnostics.
- Collocation agrees with the current shooting baseline on current blessed
  cases and with Decent-King reference values once available.
- Solver options include algorithm, tolerance, mesh/domain, maximum nodes, and
  residual scaling.
- Failing collocation solves return structured failure or throw a documented
  error; no silent fallback.

Dependencies: `UP-API4`, `UP-BVP1`, `UP-SF5` for source-backed production
acceptance.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `Project.toml` or `ext/`
- `src/inner.jl`
- possible new `src/collocation.jl`
- `src/SlenderConeRecoil.jl`
- `test/test_bead5.jl`
- `test/test_numerical_regressions.jl`

### `UP-BVP3`: Add Homotopy And Domain-Continuation Helpers

Type/priority: feature, P1.

Scope:

- Implement simple continuation ladders that reuse the previous solution as
  the initial guess for the next solve.
- Candidate axes: axial-curvature term, truncation length, far-field residual
  strictness, cone angle, regularization strength, and matching point.

Acceptance criteria:

- Continuation path, parameter sequence, failed steps, and final diagnostics
  are stored in result metadata.
- Domain continuation in `xi_max` is available for semi-infinite truncation
  sensitivity.
- Tests cover at least one monotone homotopy path and one controlled failure.

Dependencies: `UP-BVP2`; `UP-API3`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/inner.jl`
- possible new `src/continuation.jl`
- `test/test_bead5.jl`
- `test/test_numerical_regressions.jl`
- `README.md`

### `UP-BVP4`: Add Multiple-Shooting Diagnostic Solver

Type/priority: feature, P2.

Scope:

- Add a multiple-shooting path as a diagnostic fallback for the inner BVP when
  collocation and single shooting disagree.
- Use SciML shooting tools if dependency policy supports them.

Acceptance criteria:

- Multiple shooting returns the common BVP result/diagnostic object.
- At least one regression compares single shooting, multiple shooting, and
  collocation on a fixed domain.
- Documentation states that multiple shooting is a diagnostic or fallback, not
  the default production solver.

Dependencies: `UP-BVP1`; `UP-BVP2`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/inner.jl`
- possible `src/collocation.jl` or `src/shooting.jl`
- `test/test_bead5.jl`
- `docs/method.md`

### `UP-BVP5`: Prototype BifurcationKit Continuation Extension

Type/priority: feature, P3.

Scope:

- Wrap the discretized BVP residual for pseudo-arclength continuation, folds,
  and branch/parameter diagrams.
- Keep this optional and extension-based unless dependency policy says
  otherwise.

Acceptance criteria:

- A small documented example can continue a simplified or source-backed BVP in
  one parameter and report fold/branch metadata if present.
- The extension is skipped cleanly when optional dependencies are unavailable.
- No custom pseudo-arclength implementation is added to core package code.

Dependencies: `UP-BVP2`; `UP-BVP3`; stable dependency policy.

Suggested validation gates:

- `julia --project test/runtests.jl`
- optional extension-specific test command documented by the implementation
  issue

Likely changed files/modules:

- `Project.toml`
- `ext/SlenderConeRecoilBifurcationKitExt.jl`
- possible `test/ext/` files
- `docs/examples/` or `docs/method.md`

### `UP-MA1`: Implement Explicit Region And Common-Part Objects

Type/priority: feature, P1.

Scope:

- Represent inner, outer, and optional intermediate regions with variables,
  transforms, validity assumptions, truncation order, source references, and
  common parts.
- Stop treating the composite profile as only an array splice.

Acceptance criteria:

- Composite objects contain inner expression/data, outer expression/data,
  common part, assumptions, and diagnostic metadata.
- Tests can inspect the common part separately from plotted output.
- Existing composite plotting continues through a compatibility layer.

Dependencies: `UP-SF4`; `UP-API1`; remediation composite matching issue.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/composite.jl`
- `src/outer.jl`
- `src/outer_hierarchy.jl`
- possible new `src/asymptotic_regions.jl`
- `test/test_bead7.jl`
- `test/test_outer_hierarchy.jl`

### `UP-MA2`: Add Overlap-Window Diagnostics

Type/priority: feature, P1.

Scope:

- Search candidate overlap windows, fit common parts, track mismatch norms as
  the window moves, and report sensitivity to truncation order.
- Replace arbitrary splice-point confidence with quantitative overlap
  diagnostics.

Acceptance criteria:

- Matching results report overlap window, mismatch norm, fit coefficients, and
  sensitivity summary.
- Tests cover no-overlap, unstable-overlap, and stable-overlap cases.
- Figures and docs label matches with diagnostic values.

Dependencies: `UP-MA1`; `UP-API3`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/composite.jl`
- possible `src/diagnostics.jl`
- `scripts/figures.jl`
- `test/test_bead7.jl`
- `test/test_numerical_regressions.jl`

### `UP-MA3`: Add Capillary-Wave Phase And Envelope Validators

Type/priority: feature, P1.

Scope:

- Implement diagnostics for zero crossings, crest locations, local wavelength,
  wave phase, and envelope decay in `S(xi) - epsilon*xi`.
- Compare against Decent-King and Keller-Miksis asymptotic expectations where
  source data permits.

Acceptance criteria:

- Wave diagnostics are robust to grid resolution and report insufficient
  resolution clearly.
- Source-backed tests exist for any claimed phase/envelope law.
- Figure metadata records wave diagnostics alongside profile data.

Dependencies: `UP-SF4`; `UP-SF5`; `UP-API3`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- possible `src/waves.jl`
- `src/outer.jl`
- `src/inner.jl`
- `scripts/figures.jl`
- `test/test_numerical_regressions.jl`

### `UP-MA4`: Repair Or Retire Higher-Order Outer Hierarchy

Type/priority: bug, P1.

Scope:

- Decide whether the current CAS-assisted outer hierarchy supports the
  source-backed epsilon ordering, including negative powers or Laurent terms.
- Implement the supported path or throw clear unsupported errors.

Acceptance criteria:

- `derive_outer_equations` and related hierarchy APIs have tests for the
  intended expansion grammar.
- Docs state exactly which expansion families are supported.
- Unsupported source-needed hierarchy claims are removed from README and
  figures.

Dependencies: remediation Laurent path; `UP-SF4`; `UP-MA1`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl` if
  numerical outer behavior changes

Likely changed files/modules:

- `src/outer_hierarchy.jl`
- `src/series.jl`
- `src/expr.jl`
- `test/test_outer_hierarchy.jl`
- `README.md`
- `docs/method.md`

### `UP-PDE1`: Add Manufactured-Solution Tests For The Slender PDE Operator

Type/priority: task, P1.

Scope:

- Add smooth manufactured profiles for `R`, `u`, source terms, and boundary
  data on uniform and stretched grids.
- Verify derivative stencils, boundary conditions, source handling, and
  observed order independently of cone singularity.

Acceptance criteria:

- MMS tests measure spatial order for the PDE operator and fail on derivative
  or boundary regressions.
- Both uniform and stretched grids are covered.
- MMS code is separate from physical cone benchmark data.

Dependencies: remediation PDE validation/cache safety; `UP-API3` helpful but
not mandatory.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/pde.jl`
- possible `src/mms.jl`
- `test/test_bead8.jl`
- `test/test_numerical_regressions.jl`

### `UP-PDE2`: Add Conservative Diagnostics And Hard Domain Failures

Type/priority: bug, P1.

Scope:

- Prefer verification diagnostics in area `A = R^2` and flux/momentum-like
  quantities where the model permits.
- Treat `R <= 0`, nonfinite states, failed retcodes, and invalid domains as
  hard benchmark failures.

Acceptance criteria:

- PDE results include mass/flux balance, positivity margin, retcode summary,
  minimum radius, step statistics, and grid-quality diagnostics.
- Tests cover invalid radius, nonmonotone grids, boundary cases, and failed
  solve handling.
- Benchmark/figure paths refuse invalid PDE results.

Dependencies: remediation PDE domain validation; `UP-API3`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/pde.jl`
- `src/solve_checks.jl`
- possible `src/diagnostics.jl`
- `test/test_bead8.jl`
- `scripts/figures.jl`

### `UP-PDE3`: Implement Quantitative Similarity-Collapse Metrics

Type/priority: feature, P1.

Scope:

- Replace visual snapshot collapse with scalar collapse scores on trusted
  similarity windows.
- Track profile, slope, curvature, velocity, capillary-wave phase, and optional
  fitted exponent/time-offset metrics.

Acceptance criteria:

- Collapse diagnostics define and record the trusted `xi` window, weights,
  interpolation grid, norms, and excluded regions.
- Tests cover perfect manufactured collapse, perturbed collapse, and invalid
  snapshot sets.
- Existing PDE verification figures report collapse score and conservation
  diagnostics.

Dependencies: `UP-PDE2`; `UP-MA3` for wave-specific metrics.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/pde.jl`
- possible `src/collapse.jl`
- `scripts/figures.jl`
- `test/test_bead8.jl`
- `test/test_numerical_regressions.jl`

### `UP-PDE4`: Add Similarity-Frame Or Mapped-Coordinate PDE Verifier

Type/priority: feature, P2.

Scope:

- Implement a PDE verification mode in `xi = z / ell(t)` or a prescribed
  mapped coordinate so the similarity profile approaches a steady state.
- Start with known `ell(t) = (sigma*t^2/rho)^(1/3)` and allow later exponent
  fitting.

Acceptance criteria:

- Mapped verifier reaches longer effective similarity times than the current
  fixed physical grid on a reference case.
- Results include residual, conservation, positivity, and collapse diagnostics.
- Tests cover coordinate-transform consistency and at least one cheap mapped
  solve.

Dependencies: `UP-PDE1`; `UP-PDE2`; `UP-PDE3`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`

Likely changed files/modules:

- `src/pde.jl`
- possible `src/mapped_pde.jl`
- `src/similarity.jl`
- `test/test_bead8.jl`
- `README.md`

### `UP-PDE5`: Add Independent Discretization Checks

Type/priority: feature, P2.

Scope:

- Add at least one independent verification path: MethodOfLines-generated
  finite differences for the PDE operator, ApproxFun spectral residuals for
  similarity ODE/BVPs, or both.
- Keep heavyweight dependencies optional or in test/script environments.

Acceptance criteria:

- Independent discretization agrees with the package implementation on a cheap
  benchmark within documented tolerances.
- Dependency loading is optional or isolated per `UP-API4`.
- Tests or scripts make failures reproducible without changing manifests
  unexpectedly.

Dependencies: `UP-API4`; `UP-PDE1`; `UP-BVP2` for spectral BVP comparison.

Suggested validation gates:

- `julia --project test/runtests.jl`
- documented optional test/script command for the independent backend

Likely changed files/modules:

- `Project.toml`, `scripts/Project.toml`, or extension environment
- possible `ext/SlenderConeRecoilMethodOfLinesExt.jl`
- possible `ext/SlenderConeRecoilApproxFunExt.jl`
- `test/` optional backend tests
- `docs/method.md`

## Validation Datasets, Documentation, Examples, Benchmarks, Reproducibility

### `UP-REP1`: Define Reference Dataset Schema And Artifact Policy

Type/priority: task, P2.

Scope:

- Define how small reference data, extracted source values, computed benchmark
  data, and large generated artifacts are stored.
- Separate paper provenance from computed benchmark provenance.

Acceptance criteria:

- A reference-data schema records source DOI, local checksum, extraction
  method, parameters, solver settings, tolerances, and data hash.
- Licensed PDFs remain untracked and are never required for normal package
  loading.
- Small reusable benchmark data has a stable tracked location or an
  `Artifacts.toml` plan with content hashes.

Dependencies: `UP-SF5`; `UP-API2`.

Suggested validation gates:

- `git diff --check`
- `julia --project test/runtests.jl` if package artifact helpers are added

Likely changed files/modules:

- possible `Artifacts.toml`
- possible `data/reference/`
- `docs/papers/README.md`
- possible `src/artifacts.jl`
- `README.md`

### `UP-REP2`: Build A Benchmark Suite With Cheap And Reference Modes

Type/priority: feature, P2.

Scope:

- Add `benchmark/benchmarks.jl` using `BenchmarkTools`/`PkgBenchmark` patterns.
- Cover inner BVP solve, residual evaluation, outer/composite matching,
  PDE short run, collapse metric, and figure-data generation.

Acceptance criteria:

- Benchmarks have seeded deterministic inputs where relevant.
- Cheap benchmarks are quick enough for routine local comparison; reference
  benchmarks are documented separately.
- Benchmark output records Julia version, package commit, thread/BLAS settings,
  solver algorithms, tolerances, and diagnostics.

Dependencies: `UP-API4`; `UP-API3`; relevant solver/PDE issues.

Suggested validation gates:

- `julia --project test/runtests.jl`
- documented benchmark command from the implementation issue

Likely changed files/modules:

- `benchmark/benchmarks.jl`
- possible `benchmark/Project.toml`
- possible `benchmark/Manifest.toml`
- `README.md`

### `UP-REP3`: Rework Figures Into Reproducible Validation Products

Type/priority: task, P2.

Scope:

- Ensure every README-visible figure is regenerated by scripts from the same
  canonical data path.
- Store figure metadata with solver diagnostics, source references, package
  version/commit, and manifest hash where practical.

Acceptance criteria:

- Figures are either explicitly provisional or source-backed.
- `scripts/figures.jl --metadata-only` updates metadata without plot binary
  churn.
- Figure generation refuses failed solver diagnostics.
- README captions do not overclaim source fidelity.

Dependencies: remediation artifact reproducibility; `UP-API2`; `UP-API3`;
source-backed figure status from `UP-SF5`.

Suggested validation gates:

- `julia --project=scripts scripts/figures.jl --metadata-only`
- `julia --project=scripts scripts/figures.jl` when figure binaries are
  intentionally regenerated
- `git diff --check`

Likely changed files/modules:

- `scripts/figures.jl`
- `figures/metadata.toml`
- tracked figure PNG/PDF files only when intentionally regenerated
- `README.md`

### `UP-REP4`: Add Source-Fidelity Documentation And User Examples

Type/priority: task, P2.

Scope:

- Create examples that show the authoritative workflows: derive scaling,
  solve cone inner BVP, compute outer/composite match, run PDE verification,
  inspect diagnostics, and reproduce a figure.
- Keep examples as runnable package workflows, not marketing prose.

Acceptance criteria:

- Examples use public API and return inspectable diagnostics.
- Documentation distinguishes source-backed Decent-King reproduction,
  package reconstruction, and exploratory extensions.
- README has a compact quick-start plus links to deeper method docs.

Dependencies: `UP-API1`; `UP-SF3`; `UP-SF4`; `UP-SF5`.

Suggested validation gates:

- `julia --project test/runtests.jl`
- example-specific run commands documented by the implementation issue
- `git diff --check`

Likely changed files/modules:

- `README.md`
- `docs/method.md`
- possible `docs/examples/`
- possible `examples/`
- `scripts/figures.jl`

### `UP-REP5`: Create A Validation Benchmark Matrix

Type/priority: task, P2.

Scope:

- Record benchmark families with cheap/reference configurations and source
  status:
  manufactured slender PDE, linear capillary waves, Keller-Miksis scaling,
  Decent-King cone, Keller-King-Ting blob, Billingham fat wedge/cone,
  Eggers-Dupont/Papageorgiou slender jet, and later modern extensions.

Acceptance criteria:

- Each benchmark row has source DOI, local access status, implemented status,
  primary recorded quantities, validation gate, and planned tolerance basis.
- Benchmarks not yet implemented are clearly marked as future data, not
  current package claims.
- Matrix drives future Beads issue creation without requiring another broad
  synthesis pass.

Dependencies: `UP-SF1`; `UP-REP1`; partial implementations as they land.

Suggested validation gates:

- `git diff --check`

Likely changed files/modules:

- `docs/roadmap/validation_benchmarks.md` or research directory equivalent
- `README.md` summary table if appropriate

## Speculative Future Extensions

These should not block the Decent-King cone upgrade unless the source ledger
shows they are needed for correctness.

### `UP-GEN1`: Add Dimension-Aware Scaling Utilities

Type/priority: feature, P3.

Scope:

- Implement exact rational dimension-matrix/nullspace utilities for Pi groups,
  Keller-Miksis scaling, and user-selected repeating variables.

Acceptance criteria:

- Keller-Miksis `L(t) = (sigma*t^2/rho)^(1/3)` and
  `U(t) = (sigma/(rho*t))^(1/3)` are derived and tested from dimensions.
- Floating exponent output is not used for canonical derivations.
- API records whether scaling came from dimensional analysis, user input, or a
  source ledger.

Dependencies: stable API policy after M1.

Suggested validation gates:

- `julia --project test/runtests.jl`

Likely changed files/modules:

- possible `src/scaling.jl`
- `src/similarity.jl`
- `test/`
- `docs/method.md`

### `UP-GEN2`: Add Asymptotic Scale And Region IR

Type/priority: feature, P3.

Scope:

- Build a small internal representation for powers, fractional powers,
  logarithms, nested scales, regions, truncation order, and common-part
  extraction.

Acceptance criteria:

- Current outer/composite workflow can use the IR without relying on ad-hoc
  string or expression manipulation.
- Unsupported expansion forms fail clearly.
- IR can adapt to the existing custom CAS and future `Symbolics.jl` adapters.

Dependencies: `UP-MA1`; `UP-MA4`.

Suggested validation gates:

- `julia --project test/runtests.jl`

Likely changed files/modules:

- `src/series.jl`
- `src/expr.jl`
- possible `src/asymptotic_ir.jl`
- `test/test_bead1.jl`
- `test/test_bead2.jl`
- `test/test_outer_hierarchy.jl`

### `UP-GEN3`: Prototype Dominant-Balance Enumeration

Type/priority: feature, P3.

Scope:

- Enumerate candidate balances for polynomial or rationalized equation systems,
  store rejected balances and reasons, and verify asymptotic dominance after
  substitution.

Acceptance criteria:

- The engine supports a deliberately narrow grammar and rejects unsupported
  systems clearly.
- At least one simple ODE/PDE scaling example and Keller-Miksis-style balance
  example are documented.
- Rejected-balance diagnostics are test-covered.

Dependencies: `UP-GEN1`; `UP-GEN2`.

Suggested validation gates:

- `julia --project test/runtests.jl`

Likely changed files/modules:

- possible `src/dominant_balance.jl`
- `src/expr.jl`
- `src/series.jl`
- `test/`
- `docs/method.md`

### `UP-GEN4`: Add A General SimilarityProblem Interface

Type/priority: feature, P3.

Scope:

- Generalize from cone-specific APIs to a package-level
  `SimilarityProblem`/`solve` interface for scaling reductions, BVP solves,
  PDE collapse, and matched asymptotics.

Acceptance criteria:

- Cone workflows remain backward compatible.
- At least two non-cone toy examples exercise the generic interface.
- Documentation is explicit that this is a generalization layer, not a package
  rename.

Dependencies: M1 through M5 complete enough to avoid premature abstraction.

Suggested validation gates:

- `julia --project test/runtests.jl`
- `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl` if cone
  workflows are affected

Likely changed files/modules:

- `src/SlenderConeRecoil.jl`
- possible `src/problems.jl`
- possible `src/solvers.jl`
- `test/`
- `README.md`
- `docs/method.md`

### `UP-GEN5`: Add Adjacent Physics Benchmark Modules

Type/priority: feature, P3/P4.

Scope:

- Add separate modules or examples for Billingham fat wedge/cone,
  Keller-King-Ting blob, Eggers-Dupont/Papageorgiou slender jet,
  Taylor-Culick/retraction, surfactant, surface viscosity, or viscoelastic
  filament benchmarks.

Acceptance criteria:

- Each physics family has explicit regime assumptions and source citations.
- No optional physics silently changes the inviscid slender-cone default.
- Benchmarks enter the validation matrix before they become user-facing claims.

Dependencies: `UP-REP5`; source acquisition for each target benchmark.

Suggested validation gates:

- `julia --project test/runtests.jl`
- relevant slow/reference benchmark gates defined by each implementation issue

Likely changed files/modules:

- possible new model files under `src/`
- `test/`
- `docs/examples/`
- `docs/roadmap/validation_benchmarks.md`
- `README.md`

## Serial Beads Execution Order

Recommended order for actual Beads:

| Order | Plan ID | Beads ID |
| --- | --- | --- |
| 1 | UP-SF1 | `scr-8l4.1` |
| 2 | UP-SF2 | `scr-8l4.2` |
| 3 | UP-SF3 | `scr-8l4.3` |
| 4 | UP-SF4 | `scr-8l4.4` |
| 5 | UP-SF5 | `scr-8l4.5` |
| 6 | UP-API4 | `scr-8l4.6` |
| 7 | UP-API1 | `scr-8l4.7` |
| 8 | UP-API2 | `scr-8l4.8` |
| 9 | UP-API3 | `scr-8l4.9` |
| 10 | UP-BVP1 | `scr-8l4.10` |
| 11 | UP-BVP2 | `scr-8l4.11` |
| 12 | UP-BVP3 | `scr-8l4.12` |
| 13 | UP-MA1 | `scr-8l4.13` |
| 14 | UP-MA2 | `scr-8l4.14` |
| 15 | UP-MA3 | `scr-8l4.15` |
| 16 | UP-MA4 | `scr-8l4.16` |
| 17 | UP-PDE1 | `scr-8l4.17` |
| 18 | UP-PDE2 | `scr-8l4.18` |
| 19 | UP-PDE3 | `scr-8l4.19` |
| 20 | UP-PDE4 | `scr-8l4.20` |
| 21 | UP-PDE5 | `scr-8l4.21` |
| 22 | UP-REP1 | `scr-8l4.22` |
| 23 | UP-REP2 | `scr-8l4.23` |
| 24 | UP-REP3 | `scr-8l4.24` |
| 25 | UP-REP4 | `scr-8l4.25` |
| 26 | UP-REP5 | `scr-8l4.26` |
| 27 | UP-BVP4 | `scr-8l4.27` |
| 28 | UP-BVP5 | `scr-8l4.28` |
| 29 | UP-GEN1 | `scr-8l4.29` |
| 30 | UP-GEN2 | `scr-8l4.30` |
| 31 | UP-GEN3 | `scr-8l4.31` |
| 32 | UP-GEN4 | `scr-8l4.32` |
| 33 | UP-GEN5 | `scr-8l4.33` |

Each child issue depends on the previous issue in the table. `scr-8l4.1`
depends on `scr-ts1`, and `scr-8l4.33` blocks `scr-8l4`, so the umbrella
cannot close while the chain is incomplete.

Rationale for this order:

- Source access and equation fidelity decide what the production equations and
  reference data are.
- API/result/diagnostic work should precede solver rewrites so new algorithms
  have stable outputs and failure semantics.
- Collocation and continuation should precede broad PDE and benchmark work so
  the flagship similarity profile is reliable.
- Matching and wave diagnostics must be in place before figures or datasets
  claim source-backed asymptotic fidelity.
- General `SimilarityMethods`-style abstractions should wait until the cone
  workflow has proven its requirements.

## Validation Gate Policy

Use these gates when converting the plan into implementation Beads:

| Change type | Minimum gate |
| --- | --- |
| Documentation or plan only | `git diff --check` |
| Core package/API changes | `julia --project test/runtests.jl` |
| Solver, BVP, outer, composite, PDE, or numerical regression changes | `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl` in addition to the fast gate |
| Broad numerical changes where both gates are needed | `SLENDER_RECOIL_TEST_GROUP=all julia --project test/runtests.jl` |
| Figure metadata only | `julia --project=scripts scripts/figures.jl --metadata-only` |
| Intentional figure regeneration | `julia --project=scripts scripts/figures.jl` |
| Benchmark changes | documented benchmark command plus fast package gate |

Do not add GitHub Actions or external CI for these gates. Keep package/test
jobs serial in the shared project environment.

## Risk Register

| Risk | Consequence | Mitigation |
| --- | --- | --- |
| Numerical risk | Collocation, continuation, outer matching, or mapped PDE changes may converge to plausible but wrong profiles, especially with oscillatory capillary tails and semi-infinite truncation. | Require source-backed residuals, independent solver comparisons, truncation/domain continuation, wave phase/envelope diagnostics, and fail-loud retcode/residual checks before updating figures or reference data. |
| Dependency/precompile risk | Adding `BoundaryValueDiffEq`, `BifurcationKit`, `ApproxFun`, `MethodOfLines`, `BenchmarkTools`, or plotting tools to core may slow package load, increase manifest churn, and worsen Julia precompile races. | Decide dependency policy before adding packages; prefer extensions, test extras, script or benchmark environments for heavy tools; run Julia package/test jobs serially. |
| Source-access risk | Decent-King 2008, Billingham 1999, or Keller-King-Ting 1995 may remain unavailable, blocking quantitative source-fidelity benchmarks. | Keep unsupported claims provisional; use 2001 IUTAM and metadata only as lower-fidelity context; record access failures and checksums; delay all Decent-King quantitative acceptance decisions until the PDF is acquired. |
| API churn risk | Premature generalization toward `SimilarityMethods.jl` could destabilize current cone workflows and tests. | Add problem/result/diagnostic types first, preserve compatibility wrappers, defer generic `SimilarityProblem` architecture until the source-backed cone benchmark is stable, and document deprecations explicitly. |
