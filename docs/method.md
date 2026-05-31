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

### 5. Composite solution (`src/composite.jl`)

Additive composite: S_comp = S_inner + S_outer - S_overlap, where S_overlap is the common asymptotic form in the intermediate region.

### 6. Verification (`src/pde.jl`)

The full 1D PDE is solved by method of lines (2nd-order FD on a non-uniform stretched grid, FBDF implicit time integration). At late times, rescaling to similarity variables R/t^{2/3} vs z/t^{2/3} should collapse onto the computed similarity solution. This is an internal consistency check, not yet an independently verified reproduction of Decent and King's published numerical results.

## How to run

```bash
julia --project scripts/figures.jl
```

Figures are saved to `figures/`.

## Code ↔ Equations

| File | Equations |
|------|-----------|
| `src/expr.jl` | Symbolic expression types |
| `src/series.jl` | ε-expansion infrastructure |
| `src/slender.jl` | 1D slender model (mass + momentum) |
| `src/similarity.jl` | Keller-Miksis similarity ODEs |
| `src/inner.jl` | Nonlinear BVP near tip |
| `src/outer.jl` | Linearised capillary wave problem |
| `src/composite.jl` | Matched asymptotic composite |
| `src/pde.jl` | Full 1D PDE (method of lines) |

## Source Fidelity Status

- Correct cone-paper metadata: Decent and King, *IMA Journal of Applied Mathematics* **73**(1), 37--68, DOI `10.1093/imamat/hxm043`.
- The QJMAM DOI `10.1093/qjmam/hbm028` is not the cone paper.
- The local paper manifest in `docs/papers/README.md` records the expected source PDF and the quarantined mislabeled local artifact.
