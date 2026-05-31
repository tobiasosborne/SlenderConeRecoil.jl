# CLAUDE.md — SlenderConeRecoil.jl

> Archival note: this is the original implementation brief, retained for
> provenance. Current workflow rules live in `AGENTS.md`; current run commands
> and limitations live in `README.md` and `HANDOFF.md`; current work state lives
> in Beads. Where this file conflicts with those sources, the current sources
> win.

## Project identity

Reproduce the core results of Decent & King, "Surface-tension-driven flow in a slender cone" (*IMA Journal of Applied Mathematics* 73(1), 37--68, doi:10.1093/imamat/hxm043). The paper studies the recoil of a conical body of inviscid fluid under surface tension immediately after droplet pinch-off. The cone half-angle ε is small; a similarity transformation (Keller–Miksis t^{2/3} scaling) plus matched asymptotic expansions in ε yield a quantitative description of the recoiling tip and the capillary wave train propagating backward along the cone.

This is a computational reproduction, not a literature survey. The deliverable is working Julia code that derives, solves, and visualises the key mathematical objects in the paper.

## Design-pattern reference

**Before writing any module**, read `../TensorGR.jl/src/` for design patterns. Key things to absorb:

- **AST/IR representation of symbolic expressions**: TensorGR.jl builds expression trees from Mathematica-style input and manipulates them (substitution, expansion, collection by order, simplification). Reuse this pattern wholesale. Do NOT use Symbolics.jl — roll our own lightweight expression type, mirroring TensorGR.jl's approach.
- **Perturbation expansion infrastructure**: Study how TensorGR.jl handles series expansion in a small parameter (index contraction, order-by-order collection). We need the same for the ε-expansion of the Euler equations and free-surface BCs.
- **Code organisation**: One concern per file, ~200 LOC max per file. Types in one file, operations in another. Tests alongside.

If TensorGR.jl has a module or utility for any of the operations below, **wrap and reuse it** rather than reimplementing.

## Physics summary (for the agent)

### The physical problem

An axisymmetric body of inviscid, incompressible fluid (density ρ, surface tension γ) occupies a cone of small half-angle ε. At t = 0 the fluid is at rest. For t > 0 surface tension drives a flow: the capillary pressure p = −γκ is large near the tip (small radius, large curvature) and small far away, so fluid is pushed away from the tip. The tip recoils and rounds off; capillary waves propagate backward along the cone.

### Governing equations (full, before reduction)

Cylindrical coordinates (r, z), axisymmetric (no θ-dependence). Velocity potential φ(r, z, t).

- Laplace: φ_{rr} + (1/r)φ_r + φ_{zz} = 0 in the fluid domain
- Kinematic BC on free surface r = R(z,t): R_t + φ_z R_z = φ_r
- Dynamic BC (Bernoulli) on r = R(z,t): φ_t + ½(φ_r² + φ_z²) = −γκ/ρ
- Curvature: κ = 1/(R√(1+R_z²)) − R_{zz}/(1+R_z²)^{3/2}

### Step 1: Slender-body (lubrication) reduction

Expand for ε ≪ 1. Set R(z,t) = εf(z,t), φ = Φ_0(z,t) + ε²Φ_1(r,z,t) + ⋯. The Laplace equation at O(ε²) fixes the radial dependence: Φ_1 = −(r²/4)Φ_{0,zz} + ⋯. Collect the kinematic and dynamic BCs order by order. The result is a 1D system for the cross-sectional area A(z,t) = πR² and axial velocity u(z,t) = Φ_{0,z}:

- ∂A/∂t + ∂(Au)/∂z = 0
- ∂u/∂t + u ∂u/∂z = −(1/ρ) ∂p/∂z

with p determined by the curvature of R (dominated by the azimuthal term 1/R at leading order in ε).

### Step 2: Similarity substitution

No imposed length scale ⟹ use Keller–Miksis scaling:
ℓ(t) = (γt²/ρ)^{1/3}

Set z = ℓ(t)ξ, R = ℓ(t)S(ξ), u = ℓ̇(t)U(ξ). Time derivatives become algebraic in ξ. The PDE system becomes an ODE system in (S(ξ), U(ξ)).

BCs:
- S(ξ) ~ εξ as ξ → ∞ (undeformed cone)
- S'(ξ₀) = 0 at the tip ξ = ξ₀ (rounded tip, part of the solution)

### Step 3: Matched asymptotics

**Outer region** (ξ large, cone barely disturbed): linearise S = εξ + ε³s₁(ξ) + ⋯, get a linear ODE for s₁. This is a dispersive wave equation with ξ-dependent coefficients. Solve via transform or direct numerics.

**Inner region** (near tip, O(1) deformations): solve the full nonlinear similarity ODEs as a BVP.

**Match**: inner far-field ↔ outer near-field. Construct composite solution.

### Step 4: Verification

Solve the time-dependent 1D PDE system directly (method of lines) and verify convergence to the similarity solution.

## Implementation plan — Beads

Each bead is an atomic unit of work. Complete them in order. After each bead, run tests and commit.

### Bead 0: Project scaffold
- `Project.toml` with runtime deps: `DifferentialEquations`, `LinearAlgebra`
- `Test` belongs in the package test target; plotting belongs in `scripts/Project.toml`
- Directory structure: `src/`, `test/`, `scripts/`, `docs/`, `figures/`
- Read `../TensorGR.jl/src/` and write a short summary of which modules/types are reusable in `docs/tensorgr_patterns.md`

### Bead 1: Expression types (`src/expr.jl`, ≤200 LOC)
- Lightweight symbolic expression tree (mirror TensorGR.jl's design)
- Types: `Sym`, `Num`, `Add`, `Mul`, `Pow`, `Func` (for sin, cos, Bessel, etc.)
- Basic operations: substitute, differentiate (structural), pretty-print
- Test: differentiate r² w.r.t. r gives 2r

### Bead 2: Series expansion (`src/series.jl`, ≤200 LOC)
- `expand_in(expr, param, order)` — Taylor expand an expression tree in a small parameter to given order
- `collect_order(expr, param, n)` — extract the coefficient of param^n
- Test: expand (1 + εx)^{-1} to O(ε³)

### Bead 3: Governing equations and slender-body derivation (`src/slender.jl`, ≤200 LOC)
- Define the cylindrical Laplace equation, kinematic BC, dynamic BC, curvature as expression trees
- Apply the slender ansatz R = εf, φ = Φ₀ + ε²Φ₁ + ⋯
- Use Bead 2 to expand and collect at each order in ε
- Output: the 1D model equations as expression trees
- Test: verify the leading-order 1D system matches the known result (mass conservation + Euler-with-capillary-pressure)

### Bead 4: Similarity reduction (`src/similarity.jl`, ≤200 LOC)
- Define the Keller–Miksis similarity substitution as a substitution rule
- Apply it to the 1D model equations from Bead 3
- Simplify (all t-dependence should cancel)
- Output: the ODE system in (S(ξ), U(ξ))
- Test: verify t cancels completely; verify correct ODE structure

### Bead 5: Inner BVP solver (`src/inner.jl`, ≤200 LOC)
- Solve the nonlinear similarity ODE as a BVP
- Strategy: shooting from the tip. Parameters: tip position ξ₀, tip curvature S''(ξ₀). Integrate outward with DifferentialEquations.jl. Newton iterate on shooting parameters until far-field matches S ~ εξ.
- Output: S(ξ), U(ξ) on a grid, and the values ξ₀, S(ξ₀)
- Test: S(ξ)/ξ → ε as ξ → ∞; volume conservation (∫πS²dξ should equal the original cone volume up to the matching region)

### Bead 6: Outer linearised problem (`src/outer.jl`, ≤200 LOC)
- Linearise the similarity ODE about S = εξ
- Solve the resulting linear ODE (variable coefficients in ξ) numerically via DifferentialEquations.jl, shooting from large ξ inward
- Alternatively, if the equation admits a WKB/Airy-type solution, construct it symbolically
- Output: the capillary wave field s₁(ξ)
- Test: compare wavelength and amplitude against inner solution's far-field oscillations

### Bead 7: Matching and composite (`src/composite.jl`, ≤150 LOC)
- Extract asymptotic forms: inner solution for large ξ, outer solution for small ξ
- Verify overlap: they must agree in an intermediate region
- Construct additive composite: inner + outer − overlap
- Test: composite is smooth; residual in the original ODE is O(ε^N) for the appropriate N

### Bead 8: Time-dependent PDE verification (`src/pde.jl`, ≤200 LOC)
- Solve the full 1D time-dependent system (Bead 3 output) via method of lines
- Spatial discretisation: 4th-order finite differences on a stretched grid (cluster near tip)
- Time integration: SSPRK or Tsit5 from DifferentialEquations.jl
- Initial condition: the undisturbed cone R(z,0) = εz, u(z,0) = 0, on a truncated domain [0, L]
- Outflow BC at z = L (no reflection)
- Output: R(z,t) snapshots
- Test: at small times, the rescaled profile R(z,t)/ℓ(t) vs z/ℓ(t) should converge to the similarity solution from Beads 5–7

### Bead 9: Visualisation (`scripts/figures.jl`, ≤200 LOC)
- Figure 1: Similarity profiles S(ξ) for several values of ε, showing beads-on-string morphology
- Figure 2: Inner, outer, and composite solutions overlaid
- Figure 3: Time-dependent PDE snapshots, rescaled to similarity variables, collapsing onto the similarity curve
- Figure 4: Capillary wave amplitude and wavelength vs ξ, compared to WKB prediction
- Save tracked PDF/PNG artifacts to `figures/` and refresh `figures/metadata.toml`
  with `julia --project=scripts scripts/figures.jl`

### Bead 10: Documentation (`docs/method.md`, ≤500 words)
- Summary of the mathematical method
- Which equations correspond to which code
- How to run: `julia --project=scripts scripts/figures.jl`

## Legacy hard rules

These rules were part of the seed brief and are not active repository policy.
Use `AGENTS.md` for current workflow rules.

1. **≤200 LOC per file.** If a file approaches this, split.
2. **No Symbolics.jl, no SymPy.** Roll own expression types following TensorGR.jl.
3. **No panic.** If an approach isn't working after 30 minutes of debugging, stop, document the state in a `HANDOFF.md`, and move to the next bead.
4. **Test every bead.** Each bead gets a corresponding `test/test_beadN.jl`. Run tests before moving on.
5. **Read TensorGR.jl first.** Bead 0 exists precisely so you understand the codebase before writing anything.
6. **Numerical before analytical.** For Beads 5 and 6, get a numerical solution first, then (optionally) construct the analytical/WKB form and verify against it.
7. **Git commit after each bead.** Message format: `bead N: short description`.
8. **Physical units.** Work in non-dimensional form throughout. The only free parameter is ε. Set γ/ρ = 1 (absorbed into the similarity scaling).

## Key references

- Keller & Miksis, "Surface tension driven flows", SIAM J. Appl. Math. 43 (1983) 268–277 — the t^{2/3} similarity scaling
- Billingham, "Surface-tension-driven flow in fat fluid wedges and cones", JFM 397 (1999) 45–71 — the fat-cone (ε = O(1)) version with Hankel–Laplace transforms
- Decent & King, "The Recoil of A Broken Liquid Bridge", IUTAM Symp. Free Surface Flows, Springer (2001) — conference proceedings version of the slender cone problem
- Decent & King, "Surface-tension-driven flow in a slender cone", *IMA Journal of Applied Mathematics* 73(1), 37--68 (2008), doi:10.1093/imamat/hxm043 — the full paper
- Eggers, "Nonlinear dynamics and breakup of free-surface flows", Rev. Mod. Phys. 69 (1997) — review with context
- Keller, King & Ting, "Blob formation", Phys. Fluids 7 (1995) 226–228 — related 1D model

## Anti-panic protocol

If you hit a wall on a bead:
1. Write what you tried and why it failed in `HANDOFF.md`
2. Write what you think the next approach should be
3. Commit the partial state
4. Move to the next bead (especially: if symbolic beads 1–4 are hard, skip to numerical beads 5–6 which are independently doable given the known 1D equations)

The numerical beads (5, 6, 8) can be done with the 1D equations written by hand if the symbolic derivation (beads 1–4) stalls. The equations are:

```
# 1D slender model (nondimensional, γ/ρ = 1):
# ∂(R²)/∂t + ∂(R²u)/∂z = 0
# ∂u/∂t + u ∂u/∂z = ∂/∂z (1/R)    [leading order curvature pressure]
#
# Similarity form (ℓ = t^{2/3}):
# ξ = z/t^{2/3}, S = R/t^{2/3}, U = u·t^{1/3}·(3/2)
# ⟹ ODE system in S(ξ), U(ξ) — derive by substitution
```

So the project has a graceful degradation path: worst case you solve the ODEs numerically with hand-derived equations and skip the CAS derivation.
