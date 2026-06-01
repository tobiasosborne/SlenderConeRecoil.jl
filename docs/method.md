# Mathematical Method — SlenderConeRecoil.jl

This is a working implementation note for the Decent--King slender-cone problem, not yet a line-by-line reproduction certificate. The target primary source is S. P. Decent and A. C. King, "Surface-tension-driven flow in a slender cone," *IMA Journal of Applied Mathematics* **73**(1), 37--68, 2008, DOI [`10.1093/imamat/hxm043`](https://doi.org/10.1093/imamat/hxm043). Equation forms and numerical constants below reflect the current code and remain provisional until checked against that source.

## Problem

An axisymmetric body of inviscid fluid (density ρ, surface tension γ) occupies a cone of small half-angle ε. At t=0 the fluid is at rest. Surface tension drives a flow: capillary pressure is large near the tip and small far away, so fluid recoils from the tip. Capillary waves propagate backward along the cone.

## Method

### 1. Slender-body reduction (`src/slender.jl`)

For ε ≪ 1, the current implementation uses a leading-order 1D slender-body system for the cross-sectional radius R(z,t) and axial velocity u(z,t):

    ∂(R²)/∂t + ∂(R²u)/∂z = 0        (mass)
    ∂u/∂t + u ∂u/∂z = -∂/∂z(1/R)    (momentum, nondimensional γ/ρ=1)

The curvature is dominated by the azimuthal term 1/R at leading order.

### 2. Similarity reduction (`src/similarity.jl`)

With no imposed length scale, the Keller-Miksis scaling ℓ(t) = t^{2/3} applies. Setting z = ℓξ, R = ℓS(ξ), u = ℓ̇U(ξ) reduces the PDE to ODEs:

    2S + 2S'(U - ξ) + SU' = 0                    (mass)
    -(2/9)U + (4/9)(U - ξ)U' - S'/S² - S''' = 0  (momentum, with axial curvature)

### 3. Inner solution (`src/inner.jl`)

Near the tip, the current nonlinear ODEs are solved by shooting from ξ₀ with S'(ξ₀) = 0 (rounded tip). The implemented tip conditions use U(ξ₀) = (4/5)ξ₀, U'(ξ₀) = -2.

### 4. Outer solution (`src/outer.jl`)

Far from the tip, the current code linearises about the undisturbed cone S = εξ, U = 0. The perturbation (s₁, u₁) satisfies a linear ODE with ξ-dependent coefficients, integrated inward from large ξ. The correct Decent--King outer matching conditions still need primary-source verification.

The local CAS hierarchy check in `src/outer_hierarchy.jl` supports only the
package's explicit integer-power ansatz

    S = εξ + ε³σ₁ + ε⁵σ₂
    U = ε²ω₁ + ε⁴ω₂

and the series grammar implemented by `src/series.jl`: finite sums and
products, integer powers of ε including Laurent terms, and negative integer
powers of ε-dependent bases expanded around a symbolic nonzero leading
coefficient. Noninteger powers and functions may appear only as ε-free
coefficients. Half-integer powers, stretched or multiple-scale variables, and
oscillatory source-specific hierarchies such as sin(ξ/sqrt(ε)) are not
implemented and should not be described as Decent--King source-backed.

### 5. Composite solution (`src/composite.jl`)

Additive composite: S_comp = S_inner + S_outer - S_overlap, where S_overlap is the common asymptotic form in the intermediate region.

### 6. Verification (`src/pde.jl`, `src/pde_independent.jl`)

The full 1D PDE is solved by method of lines (2nd-order FD on a non-uniform stretched grid, FBDF implicit time integration). At late times, rescaling to similarity variables R/t^{2/3} vs z/t^{2/3} should collapse onto the computed similarity solution. This is an internal consistency check, not yet an independently verified reproduction of Decent and King's published numerical results.

`pde_discretization_comparison` provides a dependency-free independent operator check for cheap benchmarks. It assembles local polynomial finite-difference weights directly from moment equations on the supplied grid and compares that RHS with `pde_rhs!`, excluding a documented boundary margin by default. The default tolerances (`atol=5e-3`, `rtol=5e-2`) are implementation-level smooth-state agreement thresholds, not source-backed Decent--King validation tolerances.

## How to run

Dependency placement for core solvers, optional solver extensions, plotting,
benchmarks, and paper tooling is documented in
[`docs/dependency_policy.md`](dependency_policy.md).

```bash
julia --project=scripts scripts/figures.jl
```

Figures are saved to `figures/` as tracked PDF/PNG pairs with
`figures/metadata.toml`.

## Code ↔ Equations

| File | Equations |
|------|-----------|
| `src/expr.jl` | Symbolic expression types |
| `src/series.jl` | ε-expansion infrastructure |
| `src/slender.jl` | 1D slender model (mass + momentum) |
| `src/similarity.jl` | Keller-Miksis similarity ODEs |
| `src/inner.jl` | Nonlinear BVP near tip |
| `src/outer.jl` | Linearised capillary wave problem |
| `src/outer_hierarchy.jl` | Local integer/Laurent ε hierarchy check, not a source-backed higher-order derivation |
| `src/composite.jl` | Matched asymptotic composite |
| `src/pde.jl` | Full 1D PDE (method of lines) |
| `src/pde_independent.jl` | Dependency-free independent PDE discretization comparison |

## Source Fidelity Status

- Correct cone-paper metadata: Decent and King, *IMA Journal of Applied Mathematics* **73**(1), 37--68, DOI `10.1093/imamat/hxm043`.
- The QJMAM DOI `10.1093/qjmam/hbm028` is not the cone paper.
- The local paper manifest in `docs/papers/README.md` records the expected source PDF and the quarantined mislabeled local artifact.
