# Findings

1. **High - confirmed defect: the inner shooting Newton can accept a worse iterate and still return a "solution".** In `src/inner.jl:106-138`, `solve_inner_bvp` exits the line search after at most 10 damping attempts, then assigns `x .= x_trial` unconditionally at `src/inner.jl:132`. If none of the trial residuals improved on `F`, the method still advances to a non-descent point. The same loop also breaks on a near-singular determinant at `src/inner.jl:117`, on residual tolerance at `src/inner.jl:108`, or after exhausting `newton_iters`, but the returned `InnerSolution` has no convergence flag, final residual, iteration count, or failure mode. Downstream figures print "Converged" solely because `solve_inner_bvp` returned (`scripts/figures.jl:317-320`). This can silently bless nonconverged shooting parameters.

2. **High - confirmed defect: ODE/PDE solver return codes are never checked before extracting results.** `_shoot` calls `solve` and immediately reads `sol.t`/`sol.u` at `src/inner.jl:74-82`. The outer solvers do the same at `src/outer.jl:82-92` and `src/outer.jl:121-130`; the full and linearized hierarchy solvers do the same at `src/outer_hierarchy.jl:119-130` and `src/outer_hierarchy.jl:148-157`; the PDE solver does the same at `src/pde.jl:137-148`. None checks `sol.retcode`, endpoint coverage, finiteness, positivity of `R`/`S`, or whether a terminal event occurred. A `MaxIters`, `DtLessThanMin`, instability, or early-stop solution would be propagated as valid numerical data.

3. **High - confirmed robustness gap: the PDE formulation has no positivity/domain protection around the `1/R` singularity.** The method-of-lines RHS forms `invR_buf .= 1.0 ./ R_orig` at `src/pde.jl:84-86` and then uses that derivative in `du` at `src/pde.jl:92-95`. `solve_pde` accepts arbitrary `epsilon`, `N`, `z_min`, `z_max`, and `t_end` at `src/pde.jl:113-115` without an `isoutofdomain`, callback, variable transform, lower-radius floor, or explicit failure if `R <= 0`. The default `z_min=0.01` gives a very small initial radius for `epsilon=0.1`, so this is not just an exotic edge case. The README already notes stiffness near the truncated tip (`README.md:101-102`), but the code does not guard the singular regime.

4. **High - confirmed numerical fragility: the shooting Jacobian uses a fixed forward finite difference and determinant test.** `solve_inner_bvp` hard-codes `delta = 1e-5` at `src/inner.jl:103-115`, regardless of the scale of `(xi0, S0, S''0)` or the local conditioning of the shooting residual. It then tests `abs(det(J)) < 1e-20` at `src/inner.jl:117`, which is a poor conditioning proxy for Newton solves. This makes the solve sensitive to parameter scaling, cancellation, and near-rank-deficient residual maps. Use a scaled central difference, automatic differentiation where possible, or `NonlinearSolve.jl` with trust-region/line-search globalization and explicit convergence diagnostics.

5. **Medium-high - risk: the PDE RHS is not reentrant or thread-safe because solver scratch buffers live in `p`.** `solve_pde` stores mutable derivative buffers in the parameter tuple at `src/pde.jl:123-131`, and `pde_rhs!` mutates them on every RHS call at `src/pde.jl:73-90`. This is allocation-friendly for a serial RHS call, but it is unsafe if SciML finite-difference Jacobian evaluation, ensembles, user-level threading, or future parallelism evaluate the same problem concurrently. The user-facing warning about Julia races is relevant here: any upgrade that parallelizes Jacobian columns or runs solves concurrently must allocate per-call cache objects or build a thread-local/cache-safe operator.

6. **Medium-high - confirmed performance gap: the PDE stiff solve uses finite-difference Jacobians without sparsity information.** `solve_pde` uses `FBDF(autodiff=false)` at `src/pde.jl:137-141`. The semi-discrete state has `2N` variables and a local finite-difference stencil, but the problem provides no sparse `jac_prototype`, coloring, matrix-free Jacobian-vector product, or mass/operator structure. For larger `N`, finite-difference Jacobian construction will scale poorly and repeatedly call a buffer-mutating RHS. This is a likely bottleneck before the package can support serious convergence studies.

7. **Medium - confirmed defect: interpolation clamps out-of-range matching points silently.** `_interp_val` returns the first or last sample when `x` is outside the tabulated interval at `src/outer.jl:133-140`, and `_interp_scalar` does the same at `src/composite.jl:94-111`. These helpers are used to seed matched outer solves (`src/outer.jl:111-122`, `src/outer_hierarchy.jl:112-120`, `src/outer_hierarchy.jl:143-149`) and to build composites. If `xi_match` falls outside the inner solve's actual interval, the code silently seeds the outer solve from an endpoint instead of failing. The helpers also assume strictly increasing, unique grids and use first-order interpolation on adaptive, oscillatory solutions.

8. **Medium - confirmed fidelity gap: the composite implementation ignores its own fitted overlap model.** `inner_far_field` fits a slope and intercept for the inner far field at `src/composite.jl:20-44`, and the docstring says the composite subtracts `slope*xi + intercept` at `src/composite.jl:51-55`. The implementation instead subtracts only `epsilon*xi` at `src/composite.jl:77-84`, and `inner_far_field` is not used by `composite_solution`. If the inner far field has a residual intercept, phase, or slope mismatch, the composite can double count overlap error while appearing smooth.

9. **Medium - risk: the outer linearized solve is not posed as a robust boundary-value or asymptotic matching problem.** `solve_outer` integrates inward from `xi_max` with `[seed, 0, 0, 0]` at `src/outer.jl:75-84`, and `solve_outer_driven` chooses `seed=epsilon^2` at `src/outer.jl:95-97`. The README describes matching to the inner state (`README.md:23-27`), while the handoff notes that the zero/small seed is not the true matching data (`HANDOFF.md:36-37`). `solve_outer_matched` exists (`src/outer.jl:99-130`) but is not exported at `src/outer.jl:18`, and tests primarily exercise the seed-based path. This is a risk to ground-truth fidelity of the outer wave field.

10. **Medium - confirmed test coverage gap: numerical tests are smoke tests with loose tolerances, not convergence or failure tests.** The inner test named "3D Newton converges" only checks `abs(slopes[end] - 0.1) < 0.05` and `abs(sol.U[end]) < 0.1` at `test/test_bead5.jl:34-48`, and the ODE residual check allows `0.1` at `test/test_bead5.jl:51-65`. PDE tests run very short integrations with `N=50`, `z_min=1.0`, and `t_end <= 0.01` at `test/test_bead8.jl:47-75`, avoiding the default near-tip stiffness. There are no tests for solver `retcode`, Newton residual monotonicity, parameter validation, mesh refinement, tolerance sensitivity, long-time failure reporting, or comparison to tabulated/literature values.

11. **Medium - confirmed performance smell: residual shooting saves and copies full trajectories for endpoint-only residuals.** Each Newton residual in `src/inner.jl:98-115` calls `_shoot`, which saves the default adaptive trajectory and constructs four arrays at `src/inner.jl:76-82`, even though the residual only uses the endpoint at `src/inner.jl:99-100`. The finite-difference Jacobian needs four full ODE solves per Newton step. Use `save_everystep=false` or an endpoint-only solve for residual evaluations, reserve full trajectory extraction for the accepted final solution, and consider `remake`/cached problems.

12. **Low-medium - confirmed edge-case defect: grid and finite-difference routines lack input validation.** `stretched_grid` and `ddz!` assume `N >= 3`, finite bounds, and strictly increasing grid values (`src/pde.jl:32-63`). `solve_pde` does not validate those assumptions before calling them (`src/pde.jl:113-116`). Invalid inputs can produce bounds errors, zero denominators, negative spacing, or NaNs rather than an actionable `ArgumentError`.

13. **Low-medium - confirmed reproducibility gap: the figure script does not regenerate all committed/readme figure assets.** `scripts/figures.jl` saves only PDF files at `scripts/figures.jl:49`, `scripts/figures.jl:74`, `scripts/figures.jl:124`, `scripts/figures.jl:201`, `scripts/figures.jl:220`, `scripts/figures.jl:290`, and `scripts/figures.jl:312`, while the README embeds PNG files at `README.md:35`, `README.md:41`, `README.md:47`, `README.md:53`, `README.md:59`, and `README.md:65`. The script also hard-codes `epsilon=0.1` and matching parameters inside plotting functions (`scripts/figures.jl:25`, `scripts/figures.jl:56`, `scripts/figures.jl:81`, `scripts/figures.jl:157`, `scripts/figures.jl:208`, `scripts/figures.jl:227-228`) and records no solver diagnostics, package versions, or source commit in the outputs.

14. **Low - confirmed API ambiguity: `rescale_to_similarity` returns physical coordinates at `t=0` under a similarity-variable API.** At `src/pde.jl:160-168`, the function returns `(pde.z, pde.R[t_idx])` when `t <= 0`, while its docstring promises `(xi, S)`. The test only checks lengths for this branch at `test/test_bead8.jl:82-86`. This can silently mix physical and similarity coordinates in downstream analysis.

# Open Questions

- What numerical target should define success for the inner BVP: final residual norm, far-field slope error, agreement with Decent-King tables/figures, or convergence under `xi_max`/tolerance continuation?
- Should the package treat the current shooting code as a demonstration, or replace it with a collocation BVP formulation (`BoundaryValueDiffEq.jl`/`BoundaryValueDiffEqShooting.jl`) with continuation in `epsilon` and `xi_max`?
- What is the intended production PDE discretization: conservative variables for `R^2`, primitive variables with positivity enforcement, a moving/adaptive mesh, spectral/compact differentiation, or a dedicated dispersive PDE integrator?
- Are figure files meant to be version-controlled golden artifacts, or should CI regenerate and compare them with numerical metadata?
- Should matching use a fitted overlap window with uncertainty diagnostics, asymptotic eigenmodes, or a global solve of the coupled inner/outer boundary conditions?

# Commands Run

- `pwd`
- `rg --files`
- `ls -la reviews reviews/2026-06-01-full-audit`
- `git status --short`
- `sed -n '1,240p' src/inner.jl`
- `sed -n '1,260p' src/outer.jl`
- `sed -n '1,280p' src/composite.jl`
- `sed -n '1,320p' src/pde.jl`
- `sed -n '1,260p' src/similarity.jl`
- `sed -n '1,260p' src/slender.jl`
- `sed -n '1,260p' src/outer_hierarchy.jl`
- `sed -n '1,260p' scripts/figures.jl`
- `sed -n '260,560p' scripts/figures.jl`
- `sed -n '1,260p' test/runtests.jl`
- `sed -n '1,260p' test/test_bead1.jl`
- `sed -n '1,280p' test/test_bead2.jl`
- `sed -n '1,280p' test/test_bead3.jl`
- `sed -n '1,260p' test/test_bead4.jl`
- `sed -n '1,260p' test/test_bead5.jl`
- `sed -n '1,260p' test/test_bead6.jl`
- `sed -n '1,260p' test/test_bead7.jl`
- `sed -n '1,320p' test/test_bead8.jl`
- `rg -n "solve\\(|Rodas|FBDF|interp|ddz|Random|savefig|reltol|abstol|maxiters|retcode|conver|det\\(|zeros\\(|range\\(" src scripts test README.md docs HANDOFF.md`
- `nl -ba src/inner.jl | sed -n '1,190p'`
- `nl -ba src/outer.jl | sed -n '1,180p'`
- `nl -ba src/composite.jl | sed -n '1,170p'`
- `nl -ba src/pde.jl | sed -n '1,190p'`
- `nl -ba src/outer_hierarchy.jl | sed -n '1,190p'`
- `nl -ba scripts/figures.jl | sed -n '1,340p'`
- `nl -ba test/test_bead5.jl | sed -n '1,150p'`
- `nl -ba test/test_bead8.jl | sed -n '1,150p'`
- `nl -ba README.md | sed -n '1,140p'`
- `nl -ba HANDOFF.md | sed -n '1,90p'`
- `nl -ba docs/method.md | sed -n '1,120p'`
- `nl -ba src/SlenderConeRecoil.jl | sed -n '1,180p'`
- `rg -n "Threads|Ensemble|parallel|BLAS|retcode|saveat|callback|isoutofdomain|Positive|remake|Nonlinear|BVP|Boundary|Shooting|DifferentialEquations" src scripts test Project.toml`
- `find docs -maxdepth 3 -type f -print`
- `git ls-files -s figures/*.png figures/*.pdf scripts/figures.jl Manifest.toml Project.toml | sed -n '1,120p'`

No Julia code, package operations, precompilation, manifest updates, or tests were run.
