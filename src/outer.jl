# Outer linearised problem: linearise the similarity ODE about the
# undisturbed cone S = εξ, U = 0 to get the capillary wave field.
#
# Set S = εξ + ε³s₁(ξ) + ..., U = ε²u₁(ξ) + ...
# (powers chosen so the perturbation enters at the right order)
#
# The linearised equations for s₁, u₁ form a system of ODEs with
# ξ-dependent coefficients. We solve numerically by shooting inward
# from large ξ.
#
# From the mass ODE: 2S + 2S'(U - ξ) + SU' = 0
# Substituting S = εξ + ε³s₁, S' = ε + ε³s₁', U = ε²u₁, U' = ε²u₁':
#   2(εξ + ε³s₁) + 2(ε + ε³s₁')(ε²u₁ - ξ) + (εξ + ε³s₁)ε²u₁' = 0
#
# At O(ε): 2εξ + 2ε(-ξ) = 0 ✓ (base state satisfies)
# At O(ε³): 2ε³s₁ + 2ε(-ξ)... need to collect carefully.
#
# Actually, let's work with the resolved ODE form directly.
# The similarity system gives:
#   S' = S²[(2/9)U + (8/9)(U-ξ)] / [1 - (8/9)S(U-ξ)²]
#   U' = -2 - 2S'(U-ξ)/S
#
# Linearising: S = εξ(1 + ε²σ(ξ)), U = ε²v(ξ)
# where σ = s₁/(εξ) is the relative perturbation.
#
# For the numerical approach, we directly linearise the full ODE system.

using DifferentialEquations

export solve_outer, OuterSolution

struct OuterSolution
    ξ::Vector{Float64}
    s₁::Vector{Float64}     # perturbation to S
    u₁::Vector{Float64}     # perturbation velocity
    ε::Float64
end

# ── Linearised ODE ─────────────────────────────────────────────────────
#
# We linearise around the base state S₀ = εξ, U₀ = 0.
# Let S = εξ + δs, U = δu (small perturbations).
# S' = ε + δs', U' = δu'.
#
# Mass ODE: 2S + 2S'(U-ξ) + SU' = 0
#   2(εξ+δs) + 2(ε+δs')(δu-ξ) + (εξ+δs)δu' = 0
#   2εξ + 2δs + 2ε(δu-ξ) + 2δs'(δu-ξ) + εξ·δu' + δs·δu' = 0
#   [2εξ - 2εξ] + 2δs + 2εδu - 2ξδs' + εξδu' + O(δ²) = 0
#   2δs + 2εδu - 2ξδs' + εξδu' = 0   ... (linearised mass)
#
# Momentum ODE: -(2/9)U + (4/9)(U-ξ)U' + S'/S² = 0
#   -(2/9)δu + (4/9)(δu-ξ)δu' + (ε+δs')/(εξ+δs)² = 0
#
#   Expand 1/(εξ+δs)² ≈ 1/(εξ)² - 2δs/(εξ)³:
#   -(2/9)δu + (4/9)(-ξ)δu' + (4/9)δu·δu' + ε/(εξ)² + δs'/(εξ)² - 2εδs/(εξ)³ = 0
#   -(2/9)δu - (4/9)ξδu' + 1/(ε·ξ²) + δs'/(ε²ξ²) - 2δs/(ε²ξ³) + O(δ²) = 0
#
# The 1/(εξ²) term is the base-state residual. In fact the base state
# (S=εξ, U=0) does NOT satisfy the momentum equation — it has a residual
# of 1/(εξ²). This is because u=0 is only the initial condition, not the
# steady state. The momentum residual drives the perturbation.
#
# So the linearised momentum is:
#   -(2/9)δu - (4/9)ξδu' + δs'/(ε²ξ²) - 2δs/(ε²ξ³) = -1/(εξ²)
#
# This is an inhomogeneous linear ODE. The homogeneous part governs
# the capillary waves; the particular solution accounts for the driving.
#
# For cleaner numerics, let's define y₁ = δs, y₂ = δu, and write:
#
#   y₁' = [2y₁ + 2εy₂ + εξy₂'] / (2ξ)       ... from mass
#   y₂' = [-1/(εξ²) + (2/9)y₂ - y₁'/(ε²ξ²) + 2y₁/(ε²ξ³)] * (9ξ/4)  ... from momentum
#
# Actually, let me resolve these properly as a first-order system.

"""
    outer_rhs!(dy, y, p, ξ)

Linearised ODE system for the outer perturbation.
y = [s₁, u₁] where S = εξ + s₁, U = u₁.

From linearised mass: 2s₁ + 2εu₁ - 2ξs₁' + εξu₁' = 0
From linearised momentum: -(2/9)u₁ - (4/9)ξu₁' + s₁'/(ε²ξ²) - 2s₁/(ε²ξ³) = -1/(εξ²)

Resolve for (s₁', u₁'):
From momentum: s₁'/(ε²ξ²) = -1/(εξ²) + (2/9)u₁ + (4/9)ξu₁' + 2s₁/(ε²ξ³)
               s₁' = -ε + (2/9)ε²ξ²u₁ + (4/9)ε²ξ³u₁' + 2s₁/ξ

Substitute into mass: 2s₁ + 2εu₁ - 2ξ[-ε + (2/9)ε²ξ²u₁ + (4/9)ε²ξ³u₁' + 2s₁/ξ] + εξu₁' = 0
  2s₁ + 2εu₁ + 2εξ - (4/9)ε²ξ³u₁ - (8/9)ε²ξ⁴u₁' - 4s₁ + εξu₁' = 0
  -2s₁ + 2εu₁ + 2εξ - (4/9)ε²ξ³u₁ + (εξ - (8/9)ε²ξ⁴)u₁' = 0
  u₁' = [2s₁ - 2εu₁ - 2εξ + (4/9)ε²ξ³u₁] / [εξ - (8/9)ε²ξ⁴]

Then s₁' from the expression above.
"""
function outer_rhs!(dy, y, p, ξ)
    ε = p[1]
    s₁, u₁ = y

    if ξ < 1e-6
        dy[1] = 0.0
        dy[2] = 0.0
        return
    end

    # u₁' from resolved mass+momentum
    numer_u = 2*s₁ - 2*ε*u₁ - 2*ε*ξ + (4/9)*ε^2*ξ^3*u₁
    denom_u = ε*ξ - (8/9)*ε^2*ξ^4

    if abs(denom_u) < 1e-14
        dy[1] = 0.0
        dy[2] = 0.0
        return
    end

    u₁ξ = numer_u / denom_u

    # s₁' from momentum
    s₁ξ = -ε + (2/9)*ε^2*ξ^2*u₁ + (4/9)*ε^2*ξ^3*u₁ξ + 2*s₁/ξ

    dy[1] = s₁ξ
    dy[2] = u₁ξ
end

# ── Solver ─────────────────────────────────────────────────────────────
"""
    solve_outer(; ε=0.1, ξ_min=1.0, ξ_max=50.0, s₁_init=0.0, u₁_init=0.0)

Solve the outer linearised problem by integrating inward from ξ_max.
Initial conditions at ξ_max: s₁ and u₁ should be small (unperturbed cone).

Returns an OuterSolution.
"""
function solve_outer(; ε::Float64=0.1, ξ_min::Float64=1.0, ξ_max::Float64=50.0,
                      s₁_init::Float64=0.0, u₁_init::Float64=0.0)
    y0 = [s₁_init, u₁_init]
    tspan = (ξ_max, ξ_min)  # integrate inward
    p = [ε]

    prob = ODEProblem(outer_rhs!, y0, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-10, abstol=1e-12,
                maxiters=100_000, dtmin=1e-15)

    ξ_vals = reverse(sol.t)
    s₁_vals = reverse([sol.u[i][1] for i in eachindex(sol.u)])
    u₁_vals = reverse([sol.u[i][2] for i in eachindex(sol.u)])

    OuterSolution(ξ_vals, s₁_vals, u₁_vals, ε)
end

"""
    solve_outer_driven(; ε=0.1, ξ_min=1.0, ξ_max=50.0)

Solve with a small seed perturbation to excite the capillary wave response.
"""
function solve_outer_driven(; ε::Float64=0.1, ξ_min::Float64=1.0, ξ_max::Float64=50.0)
    # Seed: small perturbation at large ξ
    s₁_init = ε^2  # small but nonzero to seed the response
    u₁_init = 0.0
    solve_outer(ε=ε, ξ_min=ξ_min, ξ_max=ξ_max, s₁_init=s₁_init, u₁_init=u₁_init)
end
