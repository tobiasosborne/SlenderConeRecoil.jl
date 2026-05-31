# Architecture and Code Smells Review

Review Agent B, 2026-06-01. Scope: code architecture, API shape, maintainability, and code-smell review. I did not edit source files, other reports, or run Julia package operations, precompilation, manifest updates, or tests.

## Findings

### High Severity

1. **Confirmed defect: numerical solvers can return failed or unconverged results as if they were valid.**
   - `src/inner.jl:72-82` solves a shooting ODE and immediately extracts `sol.t`/`sol.u`; `src/inner.jl:137-138` returns an `InnerSolution` after the final shot. There is no `sol.retcode` check, integration endpoint validation, residual metadata, or exception on solver failure.
   - `src/inner.jl:106-135` exits Newton by `break` on tolerance (`src/inner.jl:108`), near-singular Jacobian (`src/inner.jl:117`), or exhausted iterations, but `InnerSolution` has no convergence flag, final residual, iteration count, or reason for termination (`src/inner.jl:17-26`).
   - `src/outer.jl:82-92`, `src/outer.jl:121-130`, `src/outer_hierarchy.jl:119-130`, `src/outer_hierarchy.jl:148-157`, and `src/pde.jl:137-148` follow the same pattern: solve, extract arrays, return a solution type without checking the solver status.
   - Impact: downstream figures/tests can silently consume partial integrations or failed nonlinear solves. For a future best-in-class similarity-methods package, this is a correctness and reproducibility blocker, not just ergonomics.

2. **Confirmed defect/risk: the package API is not actually exercised as a package; tests and scripts bypass the module namespace.**
   - The package entry point is `src/SlenderConeRecoil.jl:1-13`, but every test includes source files directly, e.g. `test/test_bead1.jl:1-2`, `test/test_bead5.jl:1-6`, `test/test_bead7.jl:1-8`, and `test/test_bead8.jl:1-9`.
   - `scripts/figures.jl:5-17` also includes the individual source files instead of using `SlenderConeRecoil`.
   - This means export behavior, module load order, package precompilation, missing imports, accidental globals, and public API shape are not covered by the visible test harness. It also encourages code outside `src/` to call internals such as `_interp_scalar` (`scripts/figures.jl:247-249`).

3. **Confirmed defect/risk: hidden file-order coupling makes `outer_hierarchy.jl` non-local and fragile.**
   - `src/outer_hierarchy.jl:109-120` calls `_interp_val` and `inner_rhs!`; `src/outer_hierarchy.jl:140-149` calls `_interp_val` and `outer_rhs!`.
   - `_interp_val` is defined in `src/outer.jl:133-140`, `outer_rhs!` in `src/outer.jl:39-65`, and `inner_rhs!` in `src/inner.jl:35-54`; none are imported through an explicit internal API.
   - The only reason this works is the top-level include order in `src/SlenderConeRecoil.jl:7-10` and the duplicated include order in `scripts/figures.jl:13-17`. This will not scale to optional components, submodules, or users loading only a subset of the methods.

4. **Confirmed defect/risk: shared mutable RHS buffers create a Julia race hazard if the PDE problem is ever evaluated concurrently or reused unsafely.**
   - `src/pde.jl:73-90` unpacks buffers from `p` and mutates them in-place on every RHS call.
   - Those buffers are allocated once and stored in the problem parameter tuple at `src/pde.jl:123-131`.
   - Current usage appears serial, but the API does not document that `pde_rhs!` is not thread-safe and does not guard against solver options, ensemble solves, callbacks, or user code that evaluates the same RHS/problem concurrently. This is exactly the kind of hidden state that becomes a race in Julia numerical code.

### Medium Severity

5. **Risk: public API is accidental, scattered, and too domain-specific for a general similarity-methods package.**
   - Exports are distributed inside included files (`src/expr.jl:5-6`, `src/series.jl:5`, `src/slender.jl:11`, `src/similarity.jl:14`, `src/inner.jl:15`, `src/outer.jl:18`, `src/composite.jl:10`, `src/outer_hierarchy.jl:13`, `src/pde.jl:15`) rather than curated at the module boundary (`src/SlenderConeRecoil.jl:1-13`).
   - The exported names mix symbolic CAS building blocks (`SExpr`, `Sym`, `Num`, `Add`, `Mul`, `Pow`, `Func`) with one-off physics functions (`solve_inner_bvp`, `solve_pde`) and diagnostics (`eval_sexpr`) without a layered API.
   - Some documented or practically used functions are not exported, including `solve_outer_matched` (`src/outer.jl:99-108`), `solve_outer_linearised` (`src/outer_hierarchy.jl:133-140`), `inner_far_field` (`src/composite.jl:20-27`), `far_field_slope` (`src/inner.jl:141-146`), `curvature_full` (`src/slender.jl:29-38`), and `verify_t_cancels` (`src/similarity.jl:157-167`). Direct includes hide this inconsistency.

6. **Risk: result types are concrete data dumps rather than extensible numerical result objects.**
   - `InnerSolution`, `OuterSolution`, `CompositeSolution`, `HierarchySolution`, and `PDESolution` store only `Vector{Float64}` fields and a few scalar parameters (`src/inner.jl:17-26`, `src/outer.jl:20-27`, `src/composite.jl:12-17`, `src/outer_hierarchy.jl:89-96`, `src/pde.jl:17-23`).
   - They omit solver algorithm, tolerances, retcode, residuals, match point, boundary conditions, scaling convention, original `ODESolution`, and provenance. They also cannot represent `BigFloat`, unitful quantities, AD dual numbers, static arrays, or alternative grids.
   - For a general similarity-methods library, these should likely become parameterized result/problem types with small, explicit metadata fields.

7. **Risk: solver APIs hard-code algorithms, tolerances, domains, and physical defaults.**
   - Inner solve defaults and tolerances are fixed in the signature and body (`src/inner.jl:95-104`, `src/inner.jl:121-135`), and the ODE method is hard-coded to `Rodas5P()` (`src/inner.jl:74-76`).
   - Outer and hierarchy solves hard-code `Rodas5P()` and tolerances (`src/outer.jl:82-84`, `src/outer.jl:121-122`, `src/outer_hierarchy.jl:119-120`, `src/outer_hierarchy.jl:148-149`).
   - PDE solve hard-codes the grid strategy, `FBDF(autodiff=false)`, tolerances, and snapshot plan (`src/pde.jl:113-141`).
   - Users cannot pass a solver, callback, save policy, tolerances, stopping criterion, linear solver, Jacobian, or validation policy without editing source.

8. **Confirmed code smell: duplicate interpolation helpers and repeated solution extraction create maintenance drift.**
   - `src/outer.jl:133-140` defines `_interp_val`; `src/composite.jl:89-112` defines `_interp` and `_interp_scalar` with overlapping behavior.
   - `src/outer_hierarchy.jl:113-116` and `src/outer_hierarchy.jl:143-146` depend on `_interp_val`, while `src/composite.jl:72-75` and `src/composite.jl:132-134` use `_interp_scalar`.
   - Array extraction from `sol.u` is repeated in `src/inner.jl:77-82`, `src/outer.jl:86-92`, `src/outer.jl:124-130`, `src/outer_hierarchy.jl:122-130`, and `src/outer_hierarchy.jl:151-157`.
   - This is small today, but it is a clear boundary for a shared interpolation/result-extraction utility or a structured solution wrapper.

9. **Risk: interpolation silently clamps outside domain and assumes valid grids.**
   - `_interp_val` returns endpoint values when queried outside the grid (`src/outer.jl:134-136`).
   - `_interp_scalar` does the same (`src/composite.jl:94-99`) and assumes sorted `xs` with at least two points (`src/composite.jl:100-111`).
   - This behavior is tested as desired clamping (`test/test_bead7.jl:12-22`), but it can hide failed matching or non-overlapping domains in `composite_solution` (`src/composite.jl:63-86`) and `overlap_residual` (`src/composite.jl:122-139`). A general library needs explicit extrapolation/clamping modes and validation.

10. **Risk: error handling is inconsistent and not domain-aware.**
    - Symbolic differentiation throws generic string errors for unsupported exponent/function cases (`src/expr_ops.jl:31`, `src/expr_ops.jl:40`, `src/expr_ops.jl:59`).
    - `eval_sexpr` throws `KeyError` for missing bindings through `b[e.name]` (`src/outer_hierarchy.jl:19`) and only handles `sin`, `cos`, and `sqrt` (`src/outer_hierarchy.jl:23-28`), even though differentiation knows `exp` and `log` (`src/expr_ops.jl:52-55`).
    - Numerical RHS functions suppress invalid states by setting derivatives to zero (`src/inner.jl:37-40`, `src/outer.jl:43-54`) rather than surfacing a domain violation or terminating integration with a useful reason.

11. **Risk: symbolic numeric design is too narrow for robust algebra.**
    - `Num` stores only `Rational{Int}` (`src/expr.jl:17-22`), and `Num(x::Float64)` converts floats through `rationalize(Int, x)` (`src/expr.jl:20-21`).
    - Series expansion and binomial coefficients also use `Rational{Int}` (`src/series.jl:73-95`).
    - This is good for deterministic toy algebra, but it risks overflow, surprising rationalization of measured constants, and poor interoperability with scientific numeric types. It also makes `SExpr` unsuitable as a general symbolic/numeric backend without redesign.

12. **Risk: Unicode-heavy public names improve local readability but reduce package ergonomics without ASCII aliases.**
    - Public structs expose fields such as `ξ`, `Sξ`, `Sξξ`, `ξ₀`, `S₀`, `s₁`, `s₁ξ`, and `u₁` (`src/inner.jl:18-25`, `src/outer.jl:21-26`, `src/composite.jl:13-16`, `src/outer_hierarchy.jl:90-95`).
    - Function signatures also require Unicode keywords such as `ξ₀`, `S₀`, `Sξξ₀`, `ξ_max`, and `ε` (`src/inner.jl:95-97`, `src/outer.jl:75-76`, `src/pde.jl:113-115`).
    - Julia supports this well, but a package intended as the place to come for similarity methods should provide ASCII aliases or constructor APIs for terminals, scripts, notebooks, and non-US keyboard users.

13. **Risk: source boundaries reflect bead history more than reusable package architecture.**
    - `src/SlenderConeRecoil.jl:3-11` includes all files into one flat module.
    - The README describes layers (`README.md:17-27`), but the code has no corresponding submodules, internal namespace, or problem abstractions.
    - `src/slender.jl` and `src/similarity.jl` mix derivation comments, symbolic expression constructors, and verification helpers (`src/slender.jl:27-55`, `src/slender.jl:117-149`, `src/similarity.jl:25-105`, `src/similarity.jl:157-213`). `src/outer_hierarchy.jl` mixes CAS hierarchy derivation and production nonlinear outer solves (`src/outer_hierarchy.jl:33-85`, `src/outer_hierarchy.jl:87-157`).
    - For future generalization, this likely needs separation into model definitions, symbolic transforms, discretizations/solvers, matching/composition, and examples.

### Low Severity

14. **Confirmed code smell: unresolved derivation narration remains in source comments.**
    - `src/similarity.jl:43`, `src/similarity.jl:59`, `src/similarity.jl:67`, and `src/similarity.jl:82-99` include "wait", "hmm", and self-correction commentary.
    - This is not executable behavior, but it undermines trust in the derivation boundary. Move scratch derivation to notes or replace with a clean, checked derivation.

15. **Risk: package metadata is not ready for a reusable library.**
    - `Project.toml:1-12` has no `[compat]` entry for Julia itself, no `[extras]`/`[targets]` for tests, and keeps `Plots` as a package dependency even though visible usage is in `scripts/figures.jl:5`.
    - If this grows into a methods package, plotting should likely move to examples/docs extras or a weak extension, while core numerical dependencies should be separated from optional visualization.

16. **Risk: validation of grid and problem inputs is thin.**
    - `stretched_grid` assumes sensible `N`, `z_min`, `z_max`, and `β` (`src/pde.jl:32-36`).
    - `ddz!` indexes `z[3]`, `f[3]`, `z[N-2]`, and `f[N-2]` without checking `N >= 3`, matching lengths, or monotonic grids (`src/pde.jl:43-63`).
    - Solver entry points accept physically dubious values such as nonpositive `ε`, reversed domains, or too-small grids without explicit checks (`src/inner.jl:95-97`, `src/outer.jl:75-76`, `src/pde.jl:113-115`).

## Open Questions

1. Is the near-term goal to remain a focused Decent-King reproduction, or should the next architecture explicitly target a `SimilarityMethods.jl`-style framework with this problem as one model/example?
2. Should symbolic derivation remain a custom AST, or should the package migrate derivation/manipulation to established Julia tooling once the equations are stable?
3. What should be considered public API: only high-level solves and result objects, or also symbolic constructors, RHS functions, grid tools, interpolation, and verification helpers?
4. Should Unicode names remain primary with ASCII aliases, or should public constructors/keyword APIs move to ASCII while retaining Unicode internally?
5. What failure policy is preferred for numerical methods: throw on non-success, return rich result objects with status, or both depending on a `strict` option?
6. Are threaded/ensemble PDE solves in scope? If yes, `pde_rhs!` needs a buffer ownership design before parallel execution is introduced.

## Commands Run

No Julia commands were run.

- `rg --files`
- `git status --short --branch`
- `ls -la reviews reviews/2026-06-01-full-audit 2>/dev/null`
- `nl -ba src/SlenderConeRecoil.jl`
- `nl -ba src/expr.jl`
- `nl -ba src/expr_ops.jl`
- `nl -ba src/series.jl`
- `nl -ba src/slender.jl`
- `nl -ba src/similarity.jl`
- `nl -ba src/inner.jl`
- `nl -ba src/outer.jl`
- `nl -ba src/composite.jl`
- `nl -ba src/outer_hierarchy.jl`
- `nl -ba src/pde.jl`
- `nl -ba Project.toml`
- `rg -n "^(export|using|include|struct|function)|error\\(|TODO|FIXME|wait|hmm|_interp|solve_|InnerSolution|OuterSolution|CompositeSolution|HierarchySolution|PDESolution" src test README.md docs scripts Project.toml`
- `nl -ba test/runtests.jl`
- `nl -ba README.md`
- `nl -ba docs/method.md`
- `nl -ba test/test_bead1.jl`
- `nl -ba test/test_bead2.jl`
- `nl -ba test/test_bead3.jl`
- `nl -ba test/test_bead4.jl`
- `nl -ba test/test_bead5.jl`
- `nl -ba test/test_bead6.jl`
- `nl -ba test/test_bead7.jl`
- `nl -ba test/test_bead8.jl`
- `nl -ba scripts/figures.jl`
- `nl -ba reviews/2026-06-01-full-audit/README.md`
- `rg -n "^#|bd |subagent|audit|review|architecture|AGENTS|compaction|persist" AGENTS.md HANDOFF.md SlenderConeRecoil_CLAUDE.md package.json`
- `rg -n "using SlenderConeRecoil|import SlenderConeRecoil|include\\(joinpath\\(@__DIR__, \\\"\\.\\.\\\", \\\"src\\\"|include\\(\\\"\\.\\./src|include\\(joinpath\\(@__DIR__" test scripts docs README.md`
- `rg -n "export solve_outer_matched|export solve_outer_linearised|export inner_far_field|export far_field_slope|export curvature_full|export verify|export stretched_grid|export ddz" src`
- `rg -n "Dict\\{Symbol,Float64\\}|Vector\\{Float64\\}|::Float64|::Vector\\{Float64\\}|nothing|break|return|sol\\.retcode|retcode|success|conver" src/inner.jl src/outer.jl src/outer_hierarchy.jl src/pde.jl src/composite.jl`
- `rg -n "[ξ₀₁₂₃₄₅₆₇₈₉εℓκ≤≥∂γρπθ]|s₁|u₁|σ₁|ω₁|δ" src test scripts README.md docs Project.toml`
