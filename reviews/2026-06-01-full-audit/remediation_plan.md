# Review Remediation Plan

Synthesis agent: 2026-06-01. Scope: deduplicated action plan from:

- `repo_hygiene_packaging.md`
- `architecture_code_smells.md`
- `tests_edge_cases.md`
- `ground_truth_fidelity.md`
- `numerics_reproducibility.md`

No source files, tests, package operations, precompilation, or Julia commands were run for this synthesis.

## Highest-Severity Items

The review found two trust blockers that must lead remediation:

1. The cited Decent-King 2008 source is wrong or missing locally, and the README cites the wrong DOI. Until this is corrected, literature-fidelity claims cannot be trusted.
2. Numerical solve paths can silently return failed, unconverged, or non-descent results as valid solutions. Until this is fixed, figures and tests can bless invalid numerical output.

The package-load/test-harness issue is also P0 because it can hide both public API and solver-diagnostic regressions. It should be fixed before large code changes so later checks exercise the package as users load it.

## Deduplicated Finding Map

Correctness bugs:

- Wrong/mislabeled primary source and unsupported reproduction claims.
- Solver retcodes, endpoints, Newton convergence, and residuals are not checked before returning solution objects.
- Inner Newton can accept a worse line-search iterate.
- Exported outer hierarchy appears unable to handle the Laurent expansion it advertises.
- Composite overlap subtraction contradicts its own fitted-overlap documentation.
- `rescale_to_similarity` returns physical coordinates for `t <= 0` under a similarity-coordinate API.

Ground-truth and literature issues:

- Correct Decent-King 2008 paper appears to be IMA J. Appl. Math. DOI `10.1093/imamat/hxm043`, not the README QJMAM DOI.
- The implemented primitive slender-jet model is not yet derived from the cited Decent-King potential-flow formulation.
- Tip boundary conditions, outer asymptotics, matching constants, and full-curvature retention need primary-paper verification.
- PDE "verification" does not currently verify convergence to a source-backed similarity solution.

Numerical robustness:

- PDE divides by `R` without positivity or domain protection.
- PDE RHS scratch buffers are stored in `p`, creating a future reentrancy/thread-safety hazard.
- PDE stiff solve lacks sparsity/coloring/operator structure.
- Shooting Jacobian uses fixed forward finite differences and determinant-based singularity checks.
- Interpolation silently clamps out-of-domain matching queries.

Test gaps:

- Tests and scripts bypass the package module by direct source includes.
- Numerical tests are loose smoke tests, not regression, convergence, or failure tests.
- Slow solver tests are repeated and not separated from fast checks.
- PDE grid/boundary/invalid-input cases are under-tested.
- CAS edge cases, expected failures, and Laurent behavior are under-tested.

Repo tidiness and reproducibility:

- Figure script regenerates PDFs while README embeds PNGs.
- Paper downloads have no tracked manifest, checksum/access notes, or known-missing state.
- `package-lock.json` policy conflicts with `.gitignore`.
- `HANDOFF.md` and legacy agent docs are stale in places.
- Package metadata lacks an intentional Julia compat/test-target/dependency policy.
- `.gitignore` is sparse for future Julia/docs/benchmark artifacts.

Architecture and API:

- Public API is accidental and scattered through included files.
- File-order coupling leaks internal helpers across source files.
- Result types are minimal data dumps with no status/provenance metadata.
- Solver algorithms, tolerances, grids, and save policies are hard-coded.
- Generalized `SimilarityMethods.jl`-style architecture needs research before major redesign.

## Serial Execution Order

The following Beads issues were filed for remediation and chained in this order. `scr-6qi` now depends on the final item so the later research phase remains blocked until review remediation is complete.

### P0.1 Verify Decent-King Source And Correct Fidelity Claims

Beads: `scr-dre` (P0 bug)

Category: ground-truth/literature; repo documentation.

Why first: the package should not continue building tests, figures, or architecture around an unverified primary-source claim.

Acceptance criteria:

- Correct primary-source metadata is recorded for Decent and King 2008.
- The wrong or mislabeled local paper is quarantined, renamed, or clearly documented.
- README and method docs cite the correct source and do not claim unsupported reproduction fidelity.
- Claims about model derivation, tip conditions, outer matching, curvature, and PDE verification are downgraded to provisional where the source chain is not established.
- No numerical behavior changes are bundled unless needed for documentation consistency.

### P0.2 Make Tests And Scripts Load The Package Module

Beads: `scr-pjy` (P0 bug)

Category: test gaps; architecture/API; repo tidiness.

Acceptance criteria:

- `test/runtests.jl` starts with a package-load test using `SlenderConeRecoil`.
- Tests and scripts use exported API or explicit module-qualified internals instead of relying on source includes into `Main`.
- Figure generation no longer duplicates package include order.
- Public API gaps revealed by the package-load path are either exported intentionally or marked internal.
- Package-style loading is part of the local quality gate.

### P0.3 Fail Loudly On Numerical Solver Non-Success

Beads: `scr-pc5` (P0 bug)

Category: correctness bugs; numerical robustness; test gaps.

Acceptance criteria:

- Inner, outer, hierarchy, and PDE solve paths check SciML retcodes and endpoint coverage before extracting arrays.
- Inner shooting records or throws on final residual norm, iteration count, non-descent line-search failure, and near-singular Jacobian termination.
- Returned result objects or structured errors expose solver diagnostics required by tests and figures.
- Deliberately bad inputs or impossible solve settings are covered by focused tests.
- Existing figures/scripts no longer print or imply convergence merely because a function returned.

### P1.1 Fix Or De-Advertise Outer Hierarchy Laurent Path

Beads: `scr-4vi` (P1 bug)

Category: correctness bugs; ground-truth/literature; CAS/test coverage.

Acceptance criteria:

- `derive_outer_equations` has regression tests for the intended Laurent or negative-order expansion.
- Zero-leading-term negative powers either work with documented coefficients or throw a clear unsupported error.
- README/docs do not advertise hierarchy validation unless the tested path succeeds.
- Exported hierarchy APIs have an intentional public/internal status.

### P1.2 Harden PDE Domain Validation And RHS Cache Safety

Beads: `scr-y71` (P1 bug)

Category: numerical robustness; edge cases; Julia race safety.

Acceptance criteria:

- `stretched_grid`, `ddz!`, `solve_pde`, and `pde_rhs!` validate grid sizes, monotonicity, finite bounds, epsilon, and radius/domain assumptions.
- `R <= 0` is handled through a documented domain failure, transform, floor, or callback policy.
- RHS scratch-buffer ownership is made reentrant or explicitly documented and tested as serial-only.
- Boundary-condition and invalid-input tests cover the documented behavior.
- Any future parallel PDE or ensemble work is blocked until cache safety is resolved.

### P1.3 Make Outer And Composite Matching Explicit And Validated

Beads: `scr-3cn` (P1 bug)

Category: correctness bugs; numerical robustness; ground-truth/literature.

Acceptance criteria:

- Matched outer solve is either part of the intended API or deliberately internal with documented alternatives.
- Matching points outside valid domains fail clearly instead of silently clamping.
- Composite overlap subtraction matches the documented fitted-overlap behavior, or the documentation is corrected.
- Empty overlap, decreasing/duplicate grids, mismatched vector lengths, and out-of-domain match points have focused tests.
- Provisional seeded outer solves are labeled as such if they remain.

### P1.4 Add Source-Backed Numerical Regression Tests

Beads: `scr-15n` (P1 task)

Category: test gaps; ground-truth/literature; numerical robustness.

Acceptance criteria:

- Tests encode source-backed or locally blessed reference values for inner, outer/composite, and PDE behavior.
- Tolerances are scale-aware and justified.
- ODE/PDE residual checks sample multiple regions, not only one midpoint or endpoint.
- Capillary-wave and PDE similarity-collapse claims are tested if they remain documented.
- Any claim that cannot be sourced after P0.1 is removed or marked outside the current regression contract.

### P1.5 Split Fast And Slow Quality Gates

Beads: `scr-gwn` (P1 task)

Category: test gaps; repo workflow; Julia race safety.

Acceptance criteria:

- Fast deterministic tests and slow solver/convergence tests are separable by a documented command or environment variable.
- Expensive inner/PDE fixtures are reused within files where scientifically valid.
- The default local gate remains meaningful and package-loading aware.
- Slow gates include solver diagnostics from P0.3.
- AGENTS/HANDOFF document the recommended non-concurrent Julia test workflow.

### P2.1 Make Figure And Paper Artifact Reproducibility Explicit

Beads: `scr-s4e` (P2 task)

Category: repo tidiness; reproducibility; ground-truth/literature.

Acceptance criteria:

- `scripts/figures.jl` regenerates every tracked README-visible figure format, or docs consume the canonical generated format.
- Generated figure outputs include or accompany parameter, package-version, and source-commit metadata where practical.
- `docs/papers` has a tracked manifest with DOI/URL/access notes/checksum slots and known-missing state.
- Playwright/package-lock policy is explicit and reproducible.
- Licensed PDFs remain untracked unless the project deliberately changes that policy.

### P2.2 Clean Package Metadata And Repo Hygiene

Beads: `scr-i1d` (P2 task)

Category: repo tidiness; packaging; architecture/API.

Acceptance criteria:

- `Project.toml` has intentional UUID/authors/version/Julia compat/test target/dependency policy.
- Optional plotting/test tooling is separated from core dependencies where practical.
- `HANDOFF.md` and legacy planning docs no longer contradict the current test runner or artifact policy.
- `.gitignore` covers expected Julia, docs, benchmark, scratch, and coverage outputs.
- `AGENTS.md` distinguishes orchestrator/session-close duties from read-only review-agent constraints.

## P3 And Deferred Items

Do not fold the following into review remediation unless they are necessary to satisfy P0/P1 acceptance criteria:

- Renaming or repositioning the package as a general `SimilarityMethods.jl`-style library.
- Full public API redesign, submodule layout, ASCII alias strategy, and deprecation policy.
- Parameterized problem/result types beyond the minimal diagnostic metadata needed for P0.3.
- Replacing the custom symbolic AST with a broader symbolic/numeric backend.
- Choosing best-in-class BVP/collocation/continuation algorithms for production solves.
- Redesigning the PDE discretization around conservative variables, moving meshes, positivity-preserving schemes, sparse Jacobians, or matrix-free operators.
- Implementing full-curvature Decent-King equations before the primary-source research phase confirms the target formulation.

These belong to the later research and upgrade chain: `scr-x9r`, `scr-ts1`, and `scr-8l4`.

## Beads Dependency Chain

Recommended execution chain:

1. `scr-dre`
2. `scr-pjy`
3. `scr-pc5`
4. `scr-4vi`
5. `scr-y71`
6. `scr-3cn`
7. `scr-15n`
8. `scr-gwn`
9. `scr-s4e`
10. `scr-i1d`
11. Close umbrella `scr-6qi`
12. Begin research phase `scr-x9r`

Rationale: fix source truth first, force package-level loading before substantial implementation changes, make numerical failures observable, then repair advertised mathematical paths and matching, then strengthen tests, and only then clean artifact/package hygiene.

## Commands Run By Synthesis Agent

No Julia commands, package operations, precompilation, manifest updates, or tests were run.

Commands were limited to Beads/git/status/report reads and Beads issue creation/dependency updates.
