# Outer linearised problem: linearise the similarity ODE about the
# undisturbed cone S = εξ, U = 0 to get the capillary wave field.
#
# The similarity ODEs (correct signs):
#   2S + 2S'(U - ξ) + SU' = 0              [mass]
#   -(2/9)U + (4/9)(U - ξ)U' - S'/S² = 0   [momentum]
#
# Linearise: S = εξ + δs, U = δu (small perturbations).
# S' = ε + δs', U' = δu'.
#
# Mass at O(δ): 2δs + 2εδu - 2ξδs' + εξδu' = 0
#
# Momentum at O(δ):
#   -(2/9)δu + (4/9)(-ξ)δu' - (ε+δs')/(εξ)² + 2εδs/(εξ)³ = 0
#   base-state residual: -ε/(εξ)² = -1/(εξ²) [nonzero! drives the perturbation]
#   linearised: -(2/9)δu - (4/9)ξδu' - δs'/(ε²ξ²) + 2δs/(ε²ξ³) = 1/(εξ²)
#
# This is an inhomogeneous linear ODE system for (δs, δu).

using DifferentialEquations

export solve_outer, OuterSolution

struct OuterSolution
    ξ::Vector{Float64}
    s₁::Vector{Float64}     # perturbation to S
    u₁::Vector{Float64}     # perturbation velocity
    ε::Float64
end

# ── Linearised ODE (correct signs) ────────────────────────────────────
"""
    outer_rhs!(dy, y, p, ξ)

Linearised ODE for the outer perturbation (correct momentum sign).
y = [s₁, u₁] where S = εξ + s₁, U = u₁.

From linearised mass: 2s₁ + 2εu₁ - 2ξs₁' + εξu₁' = 0
From linearised momentum: -(2/9)u₁ - (4/9)ξu₁' - s₁'/(ε²ξ²) + 2s₁/(ε²ξ³) = 1/(εξ²)

Resolve for (s₁', u₁'):
From momentum: s₁'/(ε²ξ²) = -(2/9)u₁ - (4/9)ξu₁' + 2s₁/(ε²ξ³) - 1/(εξ²)
               s₁' = -(2/9)ε²ξ²u₁ - (4/9)ε²ξ³u₁' + 2s₁/ξ - ε

Substitute into mass to get u₁', then back-substitute for s₁'.
"""
function outer_rhs!(dy, y, p, ξ)
    ε = p[1]
    s₁, u₁ = y

    if ξ < 1e-6
        dy[1] = 0.0
        dy[2] = 0.0
        return
    end

    # From momentum resolved for s₁':
    #   s₁' = -(2/9)ε²ξ²u₁ - (4/9)ε²ξ³u₁' + 2s₁/ξ - ε
    #
    # Substitute into mass: 2s₁ + 2εu₁ - 2ξ[-(2/9)ε²ξ²u₁ - (4/9)ε²ξ³u₁' + 2s₁/ξ - ε] + εξu₁' = 0
    #   2s₁ + 2εu₁ + (4/9)ε²ξ³u₁ + (8/9)ε²ξ⁴u₁' - 4s₁ + 2εξ + εξu₁' = 0
    #   -2s₁ + 2εu₁ + 2εξ + (4/9)ε²ξ³u₁ + (εξ + (8/9)ε²ξ⁴)u₁' = 0
    #
    #   u₁' = (2s₁ - 2εu₁ - 2εξ - (4/9)ε²ξ³u₁) / (εξ + (8/9)ε²ξ⁴)

    numer_u = 2*s₁ - 2*ε*u₁ - 2*ε*ξ - (4/9)*ε^2*ξ^3*u₁
    denom_u = ε*ξ + (8/9)*ε^2*ξ^4    # always positive for ξ > 0, ε > 0

    if abs(denom_u) < 1e-14
        dy[1] = 0.0
        dy[2] = 0.0
        return
    end

    u₁ξ = numer_u / denom_u

    # s₁' from momentum
    s₁ξ = -(2/9)*ε^2*ξ^2*u₁ - (4/9)*ε^2*ξ^3*u₁ξ + 2*s₁/ξ - ε

    dy[1] = s₁ξ
    dy[2] = u₁ξ
end

# ── Solver ─────────────────────────────────────────────────────────────
"""
    solve_outer(; ε=0.1, ξ_min=1.0, ξ_max=50.0, s₁_init=0.0, u₁_init=0.0)

Solve the outer linearised problem by integrating inward from ξ_max.
"""
function solve_outer(; ε::Float64=0.1, ξ_min::Float64=1.0, ξ_max::Float64=50.0,
                      s₁_init::Float64=0.0, u₁_init::Float64=0.0)
    y0 = [s₁_init, u₁_init]
    tspan = (ξ_max, ξ_min)  # integrate inward
    p = [ε]

    prob = ODEProblem(outer_rhs!, y0, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-10, abstol=1e-12,
                maxiters=100_000)

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
    s₁_init = ε^2
    u₁_init = 0.0
    solve_outer(ε=ε, ξ_min=ξ_min, ξ_max=ξ_max, s₁_init=s₁_init, u₁_init=u₁_init)
end
