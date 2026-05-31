# Tests, Edge Cases, and Quality Gates Review

Review Agent C, 2026-06-01. Scope: static inspection of tests and implementation surfaces only. No Julia commands, package operations, precompilation, manifest updates, or full tests were run.

## Findings

### 1. High - Outer hierarchy CAS path is untested and appears unable to handle the Laurent expansion it advertises

Type: confirmed coverage defect; likely runtime defect by static inspection.

The README claims the CAS-derived outer hierarchy confirms higher-order epsilon corrections at epsilon = 0.1 (`README.md:63`). The package exports `derive_outer_equations` and related hierarchy APIs (`src/outer_hierarchy.jl:13`), and the derivation substitutes `S = epsilon*xi + epsilon^3*sigma1 + ...` before expanding `S'/S^2` (`src/outer_hierarchy.jl:51`-`src/outer_hierarchy.jl:72`). That necessarily requires a Laurent expansion because the epsilon-free coefficient of `S` is zero.

The series expander for negative powers assumes a nonzero epsilon-free part: it computes `a0 = collect_order(base, epsilon, 0)` and then forms `delta / a0` (`src/series.jl:73`-`src/series.jl:95`). For the outer hierarchy base, `a0` is `Num(0)`. The `pow` constructor eagerly evaluates numeric integer powers (`src/expr.jl:107`-`src/expr.jl:112`), so a zero base with a negative exponent has no guarded path. Tests only cover negative powers around `1 + epsilon*x` (`test/test_bead2.jl:18`-`test/test_bead2.jl:32`, `test/test_bead2.jl:55`-`test/test_bead2.jl:63`) and never call `derive_outer_equations`, `solve_outer_full`, `solve_outer_linearised`, or `eval_sexpr`.

Focused tests to add:

- A unit test for `expand_in((epsilon*x)^-2, epsilon, order)` or the intended Laurent equivalent, including negative collected orders.
- A smoke/regression test for `derive_outer_equations(order=5)` that asserts the expected nonzero orders and key first-order mass/momentum terms.
- A package-level test that imports `SlenderConeRecoil` and verifies the exported hierarchy functions are callable.

### 2. High - Solver failures and Newton non-convergence are not observable to tests

Type: confirmed testability defect and solver robustness risk.

`solve_inner_bvp` can stop because `norm(F) < newton_tol`, because the finite-difference Jacobian determinant is tiny, or because the iteration budget is exhausted, but it always returns an `InnerSolution` without a convergence flag, final residual, iteration count, or failure exception (`src/inner.jl:95`-`src/inner.jl:138`). `_shoot` calls `solve` and immediately extracts `sol.t`/`sol.u` without checking `retcode` or endpoint reach (`src/inner.jl:72`-`src/inner.jl:82`). The outer and PDE solvers follow the same pattern: `solve_outer` and `solve_outer_matched` extract arrays without checking `retcode` (`src/outer.jl:75`-`src/outer.jl:92`, `src/outer.jl:121`-`src/outer.jl:130`), and `solve_pde` does likewise (`src/pde.jl:137`-`src/pde.jl:148`).

The tests then assert weak postconditions rather than solver success. The inner test accepts a far-field slope within `0.05` of `0.1` and velocity within `0.1` (`test/test_bead5.jl:35`-`test/test_bead5.jl:48`), and its residual check uses one finite-difference midpoint with tolerance `0.1` (`test/test_bead5.jl:51`-`test/test_bead5.jl:65`). Outer tests mostly assert finite arrays, loose endpoints, or arbitrary boundedness (`test/test_bead6.jl:10`-`test/test_bead6.jl:41`). PDE tests assert construction, positivity, and finite rescaling, but not solver return status (`test/test_bead8.jl:47`-`test/test_bead8.jl:87`).

Focused tests to add:

- Unit tests for deliberately bad solver inputs once the API exposes convergence metadata or throws structured errors.
- Assertions that the ODE/PDE solvers reached the requested endpoint and returned a successful retcode.
- Regression tests for final inner BVP residual components, not only derived shape features.

### 3. High - Tests bypass the package module, so load-time, export, and public API regressions can pass

Type: confirmed quality-gate defect.

The full test entrypoint directly includes individual bead tests (`test/runtests.jl:3`-`test/runtests.jl:11`), and each bead test directly includes source files into `Main` (`test/test_bead5.jl:2`-`test/test_bead5.jl:6`, `test/test_bead8.jl:2`-`test/test_bead8.jl:9`). No test uses `using SlenderConeRecoil`. The actual package entrypoint includes all source files inside the module (`src/SlenderConeRecoil.jl:1`-`src/SlenderConeRecoil.jl:11`), which is a different load path from the tests.

This can miss missing exports, namespace coupling, duplicate-definition behavior, module-only load failures, and public API drift. It also means tests rely on unexported helpers such as `_shoot`, `inner_rhs!`, `outer_rhs!`, `_interp_scalar`, and `stretched_grid` being present in `Main`, rather than testing the package as a user would. `HANDOFF.md` is also stale about this area: it says `test/runtests.jl` only covers beads 1-4 (`HANDOFF.md:35`-`HANDOFF.md:37`), while the current entrypoint includes beads 1-8 (`test/runtests.jl:4`-`test/runtests.jl:11`).

Focused tests to add:

- A first testset that does `using SlenderConeRecoil` and checks exported API availability.
- Move internal helper tests behind explicit module qualification or a deliberate test-only include strategy.
- A documentation consistency check or review item for run instructions and stale handoff claims.

### 4. Medium-High - Numerical tolerances are often loose, dimensional, or only test "finite and bounded"

Type: confirmed regression-coverage weakness.

Several assertions would allow large scientific regressions. Inner slope tolerance is `0.05` against a target `0.1`, a 50 percent relative tolerance (`test/test_bead5.jl:39`-`test/test_bead5.jl:44`). Inner residual tolerances are absolute `0.1` at a single midpoint (`test/test_bead5.jl:51`-`test/test_bead5.jl:65`). Outer perturbations are accepted below `100.0` (`test/test_bead6.jl:37`-`test/test_bead6.jl:41`). Composite slopes are accepted below `100.0`, and overlap residuals below `10.0` (`test/test_bead7.jl:50`-`test/test_bead7.jl:71`). PDE coverage checks finite values and positivity but not accuracy against a reference solution (`test/test_bead8.jl:47`-`test/test_bead8.jl:87`).

These tests are useful smoke checks, but they are not strong regression gates for a numerical research package. They should be complemented by scale-aware tolerances, golden values, and residual norms computed across multiple points.

Focused tests to add:

- Small deterministic reference fixtures for inner, outer, composite, and PDE runs.
- Relative/absolute tolerance pairs tied to physical scales, not broad constants.
- Multiple-point residual checks for mass and momentum, including near-tip, overlap, and far-field regions.

### 5. Medium-High - Full suite has hidden slow paths and repeated expensive solves

Type: confirmed local quality-gate risk.

The test suite repeatedly recomputes stiff solver outputs. `test_bead5` calls `solve_inner_bvp` in three separate testsets (`test/test_bead5.jl:35`, `test/test_bead5.jl:52`, `test/test_bead5.jl:72`). `test_bead7` calls `solve_inner_bvp` three more times and `solve_outer_driven` three times (`test/test_bead7.jl:39`-`test/test_bead7.jl:71`). `test_bead8` calls `solve_pde` four times (`test/test_bead8.jl:47`-`test/test_bead8.jl:87`). README and handoff both identify slow solver tests (`README.md:75`, `HANDOFF.md:41`-`HANDOFF.md:44`).

This is risky in a shared Julia repo where concurrent agents may trigger precompilation, solver compilation, or long test jobs. It also encourages users to skip tests because the cheap and expensive checks are not separated.

Focused tests to add:

- Split fast unit tests from slow integration tests via an environment variable or explicit test file grouping.
- Reuse one computed inner/PDE fixture per file where acceptable.
- Add tiny RHS/unit tests that cover solver algebra without launching expensive integrations.

### 6. Medium - PDE and finite-difference boundary/invalid-input cases are mostly untested

Type: confirmed edge-case coverage gap.

`ddz!` assumes at least three grid points and indexes `z[3]`, `f[3]`, `z[N-2]`, and `f[N-2]` with no validation (`src/pde.jl:43`-`src/pde.jl:63`). `stretched_grid` accepts any `N`, `z_min`, and `z_max` without validation (`src/pde.jl:32`-`src/pde.jl:35`). Tests use only `N=50` or `N=100` and skip boundary points in derivative accuracy checks (`test/test_bead8.jl:24`-`test/test_bead8.jl:45`). The PDE RHS boundary conditions are set after derivative computation (`src/pde.jl:92`-`src/pde.jl:103`), but no test calls `pde_rhs!` directly to assert those boundary values or conservation behavior.

Focused tests to add:

- Exactness tests for constants and linear functions, including endpoints.
- Explicit error behavior for `N < 3`, non-monotone grids, `z_min >= z_max`, and nonpositive epsilon.
- Direct `pde_rhs!` tests for left and right boundary conditions and sign of the capillary-pressure term.

### 7. Medium - Ground-truth fidelity is not turned into regression tests

Type: confirmed scientific regression risk.

The project documentation records physical and literature-facing claims: Decent and King reproduction (`README.md:3`, `README.md:107`), inner blob values (`README.md:33`, `HANDOFF.md:25`-`HANDOFF.md:28`), capillary wave count and decay (`README.md:37`-`README.md:39`), PDE similarity convergence (`README.md:27`, `HANDOFF.md:27`), and known limitations (`README.md:99`-`README.md:103`). The tests do not assert published or documented constants such as tip position, tip radius, capillary-wave count, similarity collapse error, or sensitivity to the axial-curvature term.

Current physics tests mostly spot-check symbolic arithmetic or use the implementation as its own oracle (`test/test_bead3.jl:16`-`test/test_bead3.jl:90`, `test/test_bead4.jl:34`-`test/test_bead4.jl:104`). The momentum sign has one regression check through the inner solution (`test/test_bead5.jl:69`-`test/test_bead5.jl:76`), but the previous sign-error history in `HANDOFF.md:18`-`HANDOFF.md:23` argues for more direct sign and equation-residual tests.

Focused tests to add:

- Golden tests against documented inner values with tolerances appropriate to current solver limitations.
- A capillary-wave regression that verifies the dispersive `S'''` term changes the solution qualitatively.
- PDE similarity-collapse checks against a stored similarity profile for short, deterministic runs.

### 8. Medium - Symbolic CAS tests cover happy paths but not algebraic edge cases or expected failures

Type: confirmed edge-case coverage gap.

The expression tests cover basic constructors, a few differentiations, substitution, printing, equality, and hashing (`test/test_bead1.jl:8`-`test/test_bead1.jl:110`). They do not test combining like symbolic terms, repeated factors, canonical ordering stability under nested expressions, fractional or negative powers of symbolic expressions, zero/negative numeric power behavior, or expected errors for unsupported functions and nonconstant exponents. Source paths for unsupported cases exist (`src/expr_ops.jl:27`-`src/expr_ops.jl:60`), but there are no `@test_throws` assertions.

The series tests cover polynomial products and negative powers only when the leading term is nonzero (`test/test_bead2.jl:9`-`test/test_bead2.jl:70`). They do not cover Laurent terms, mixed negative and positive epsilon orders, zero leading terms, functions with epsilon-dependent arguments, or invalid expansion orders. These are directly relevant because the outer hierarchy asks `collect_order` to scan from negative orders (`src/outer_hierarchy.jl:74`-`src/outer_hierarchy.jl:84`).

Focused tests to add:

- `@test_throws` for documented unsupported differentiation and function cases.
- Canonicalization tests for repeated symbolic factors and nested additions/multiplications.
- Laurent/negative-order `collect_order` tests before relying on the outer hierarchy.

### 9. Medium-Low - Interpolation and composite helpers assume clean inputs without tests for failure modes

Type: risk.

`_interp_scalar` assumes sorted, unique, same-length vectors and clamps out-of-range queries (`src/composite.jl:90`-`src/composite.jl:112`). `inner_far_field` silently falls back to the last points when the match region is too small (`src/composite.jl:27`-`src/composite.jl:43`). `composite_solution` builds an overlap grid from `max(inner.ξ[1], outer.ξ[1])` to `min(inner.ξ[end], outer.ξ[end])` without checking that a nonempty overlap exists (`src/composite.jl:58`-`src/composite.jl:86`). Tests cover only sorted happy-path interpolation and in-range composite construction (`test/test_bead7.jl:12`-`test/test_bead7.jl:71`).

Focused tests to add:

- Explicit behavior for empty overlap, duplicate grid points, decreasing grids, mismatched vector lengths, and custom `xi_grid` outside the valid domain.
- A test that distinguishes intended clamping from accidental extrapolation in matching diagnostics.

## Open Questions

- Should slow solver and PDE checks run by default in `test/runtests.jl`, or should the default suite be fast and deterministic with a separate integration suite?
- Which published or locally trusted numerical values should become golden references for tip position, tip radius, velocity peak, wave count, and similarity-collapse error?
- Should invalid numerical inputs throw `DomainError`/`ArgumentError`, or should solvers clamp and report warnings? Tests should encode whichever policy is chosen.
- Is `outer_hierarchy.jl` intended to be public, production code? It is exported, documented, and advertised, so the test standard should currently treat it as public.
- What CI environments and Julia versions are expected? Solver tolerances may need platform-aware bounds once the package is tested through its module entrypoint.

## Commands Run

No Julia commands were run.

- `pwd && rg --files -g '!*Manifest*' -g '!*.png' -g '!*.pdf'`
- `git status --short --branch`
- `bd onboard`
- `nl -ba test/runtests.jl`
- `nl -ba test/test_bead1.jl`
- `nl -ba test/test_bead2.jl`
- `nl -ba test/test_bead3.jl`
- `nl -ba test/test_bead4.jl`
- `nl -ba test/test_bead5.jl`
- `nl -ba test/test_bead6.jl`
- `nl -ba test/test_bead7.jl`
- `nl -ba test/test_bead8.jl`
- `nl -ba README.md | sed -n '1,220p'`
- `nl -ba src/expr.jl | sed -n '1,260p'`
- `nl -ba src/expr_ops.jl | sed -n '1,260p'`
- `nl -ba src/series.jl | sed -n '1,280p'`
- `nl -ba src/SlenderConeRecoil.jl | sed -n '1,220p'`
- `nl -ba src/slender.jl | sed -n '1,260p'`
- `nl -ba src/similarity.jl | sed -n '1,300p'`
- `nl -ba src/inner.jl | sed -n '1,360p'`
- `nl -ba src/outer.jl | sed -n '1,320p'`
- `nl -ba src/composite.jl | sed -n '1,340p'`
- `nl -ba src/pde.jl | sed -n '1,420p'`
- `nl -ba src/outer_hierarchy.jl | sed -n '1,420p'`
- `nl -ba Project.toml`
- `find . -maxdepth 4 -type f \( -name '*test*' -o -name '*log*' -o -name '*out*' -o -name '*md' \) -print | sort`
- `rg -n "derive_outer_equations|solve_outer_full|solve_outer_linearised|solve_outer_matched|eval_sexpr|HierarchySolution|far_field_slope|inner_rhs!|outer_rhs!|pde_rhs!|ddz!|stretched_grid|rescale_to_similarity|verify_t_cancels|verify_slender_derivation" src test docs README.md HANDOFF.md SlenderConeRecoil_CLAUDE.md`
- `rg -n "TODO|FIXME|BUG|HACK|Known limitations|limitation|fails|slow|conver|tol|toler|maxiters|retcode|isfinite|≈|< 100|< 10|< 0\\.1|< 0\\.05|1e-" src test README.md HANDOFF.md docs`
- `find reviews -maxdepth 3 -type f -print -exec nl -ba {} \;`
- `nl -ba HANDOFF.md | sed -n '1,260p'`
- `nl -ba docs/method.md | sed -n '1,260p'`
- `rg -n "ξ₀|S₀|error|convergence|residual|Decent|King|Keller|Miksis|PDE|similarity|capillary|wave|blob|slope|velocity" README.md HANDOFF.md docs/method.md test src`
- `rg -n "include\\(|using SlenderConeRecoil|Pkg\\.test|@testset|@test |@test_throws|@test_broken|@test_skip" test src/SlenderConeRecoil.jl Project.toml README.md HANDOFF.md`
- `test -e reviews/2026-06-01-full-audit/tests_edge_cases.md; printf '%s\n' $?`
- `ls -la reviews/2026-06-01-full-audit`
- `git status --short --branch`
- `nl -ba reviews/2026-06-01-full-audit/tests_edge_cases.md | sed -n '1,260p'`
- `find reviews/2026-06-01-full-audit -maxdepth 1 -type f -printf '%f\n' | sort`
