# Dependency And Extension Policy

Date: 2026-06-01

This policy resolves upgrade-plan item `UP-API4` for where dependencies belong
as `SlenderConeRecoil.jl` grows from a reconstructed cone example into a
source-fidelity similarity-methods package. This issue makes no
`Project.toml`, `Manifest.toml`, or environment changes unless a separate
inconsistency requires it; the current root package still has only
`DifferentialEquations` and `LinearAlgebra` as runtime dependencies, plus
`Test` as a test extra.

## Goals

- Keep `using SlenderConeRecoil` fast enough for API, provenance, symbolic,
  and lightweight numerical workflows.
- Keep heavy optional algorithms available without forcing them onto basic
  package load or unrelated tests.
- Isolate plotting, benchmark, paper-fetch, and browser tooling from the core
  library environment.
- Avoid Julia package-operation races by running package, precompile, and test
  jobs serially in each project environment.

## Dependency Categories

| Category | Location | Policy |
| --- | --- | --- |
| Core runtime | root `Project.toml` `[deps]` | Small, stable dependencies required for the default public API and current production algorithms. A core dependency must be needed by normal package use, not only by figures, benchmarks, experiments, paper acquisition, or optional solver backends. |
| Standard libraries and extras | root `Project.toml` stdlibs, `[extras]`, `[targets]` | Julia stdlibs used by the library may be listed directly. Test-only dependencies belong in `[extras]` and the `test` target if they are needed by ordinary tests without changing basic package load. |
| Scripts environment | `scripts/Project.toml` | Figure generation, metadata refresh, examples that write tracked artifacts, and plotting dependencies. Scripts may depend on the root package through `[sources]` and may carry their own manifest when reproducibility requires it. |
| Optional extensions and weak deps | future root `[weakdeps]` plus `ext/` | Heavy solver, continuation, spectral, or PDE-tool integrations that add methods when the user explicitly loads the dependency. Extension APIs must fail loudly with a clear message when the optional package is absent. |
| Benchmark environment | future `benchmark/Project.toml` | Performance measurement, reference benchmark suites, and tools such as `BenchmarkTools` or `PkgBenchmark`. Benchmarks should be reproducible but must not become root runtime dependencies. |
| Paper-fetch and browser tooling | separate local/tooling environment or documented external commands | DOI resolution, institutional-access helpers, browser automation, PDF checksums, and source acquisition tools. These are workflow tools, not package dependencies, and should not be added to the library or scripts environment unless a tracked script truly requires them. |

## Current Placements

| Package/tool | Placement | Rationale |
| --- | --- | --- |
| `DifferentialEquations.jl` | Core runtime | Current inner, outer, and PDE implementations already use SciML ODE/PDE integrators as production package functionality. |
| `LinearAlgebra` | Core runtime stdlib | Lightweight stdlib used by numerical code. |
| `Test` | Root extra and test target | Test-only stdlib; no runtime load cost. |
| `Plots.jl` and plotting backends | Scripts environment | Figures are tracked reproducibility artifacts, but plotting is not part of package load. Current placement in `scripts/Project.toml` is correct. |
| `BoundaryValueDiffEq.jl` | Future optional extension or, after evidence, core solver dependency | Use an extension first for MIRK/FIRK/Ascher collocation and independent BVP verification. Promote only if it becomes the default source-fidelity BVP solver used by the public API and fast package tests. |
| `NonlinearSolve.jl` | Future optional extension; possible core promotion for public residual solving | Use for explicit nonlinear residual APIs, shooting residuals, or collocation systems once those APIs exist. Promote only if the package's primary solve path requires it and a load-time/precompile review shows acceptable cost. |
| `BifurcationKit.jl` | Future optional extension | Continuation, fold detection, branch switching, and parameter diagrams are advanced workflows. They should live behind `SlenderConeRecoilBifurcationKitExt` or a similarly named extension. |
| `ApproxFun.jl` | Future optional extension or verification environment | Spectral residual and collocation checks are valuable independent verification paths, but they are not part of core package load. Use an extension if public methods are exposed; otherwise keep prototypes in scripts, research, or benchmark tooling. |
| `MethodOfLines.jl` | Future optional extension or verification environment | The current PDE code is hand-discretized. MethodOfLines-generated discretizations should be isolated for independent operator checks or optional public PDE workflows. |
| `BenchmarkTools.jl` | Future benchmark environment | Needed for `BenchmarkGroup` suites and performance tracking, not for library runtime or ordinary tests. |
| `PkgBenchmark.jl` and related benchmark tools | Future benchmark environment | Useful for package-level benchmark comparison and reports; keep outside core and scripts unless a benchmark issue explicitly creates that environment. |
| Paper-fetch, browser, checksum, and PDF tooling | Separate tooling only | Source acquisition and provenance workflows should write paper manifests or research notes, but their tools should not become package dependencies. |

## Movement Triggers

Move a dependency from scripts, benchmark, or research tooling to an optional
extension when all of the following are true:

- the package provides a stable user-facing method or result type for that
  capability;
- absence of the dependency can be handled by Julia extension loading or by a
  clear fail-loud error;
- tests can exercise the integration without slowing the default package load
  path; and
- the dependency has documented compatibility bounds and a validation command.

Move an optional dependency into core only when all of the following are true:

- the default public API cannot provide its documented behavior without the
  dependency;
- replacing local code with the dependency materially improves source fidelity,
  numerical robustness, or maintenance;
- `using SlenderConeRecoil` and the fast gate remain acceptable after a
  measured load/precompile review; and
- the relevant Beads issue records why extension loading is no longer the right
  boundary.

Move a core dependency out to an extension or environment when it is no longer
required by default package functionality, when it causes disproportionate
precompile/load cost, or when its use is limited to figures, benchmarks,
source acquisition, or exploratory verification.

## Command Discipline

Run Julia package, precompile, test, figure, and benchmark jobs serially within
each project environment. Do not start concurrent root package tests, slow
solver tests, figure generation, package instantiation, or manifest-changing
commands against this checkout.

Recommended validation commands:

```bash
# Documentation-only dependency-policy changes
git diff --check

# If root package code or root Project.toml changes
julia --project test/runtests.jl

# If solver, BVP, outer, composite, PDE, or numerical regression behavior changes
SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl

# If both fast and slow gates are needed, use one Julia process
SLENDER_RECOIL_TEST_GROUP=all julia --project test/runtests.jl

# If figure scripts or plotting environment changes
julia --project=scripts scripts/figures.jl --metadata-only
julia --project=scripts scripts/figures.jl
```

Future benchmark issues should define their own serial command, for example a
`julia --project=benchmark benchmark/benchmarks.jl` style command after the
benchmark environment exists.

## Implementation Notes

- Extension files should be named after the optional integration, for example
  `ext/SlenderConeRecoilBoundaryValueDiffEqExt.jl`,
  `ext/SlenderConeRecoilBifurcationKitExt.jl`,
  `ext/SlenderConeRecoilApproxFunExt.jl`, or
  `ext/SlenderConeRecoilMethodOfLinesExt.jl`.
- Optional integrations should return the same package result and diagnostic
  types as core solvers where practical, so benchmarks and provenance metadata
  do not fork by backend.
- Benchmark and figure environments may depend on the root package via
  `[sources] SlenderConeRecoil = {path = ".."}`.
- Dependency changes should be bundled with the source, tests, documentation,
  and Beads notes that justify them; avoid standalone manifest churn.
