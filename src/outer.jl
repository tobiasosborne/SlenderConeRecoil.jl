# Outer linearised problem with axial curvature.
#
# The similarity ODEs (with axial curvature):
#   2S + 2S'(U - ξ) + SU' = 0                          [mass]
#   -(2/9)U + (4/9)(U - ξ)U' - S'/S² - S''' = 0       [momentum]
#
# Linearise about S = εξ, U = 0: S = εξ + δs, U = δu.
# State vector: [s₁, s₁', s₁'', u₁] — 4 components.
#
# Linearised mass (unchanged): 2s₁ + 2εu₁ - 2ξs₁' + εξu₁' = 0
# Linearised momentum: -(2/9)u₁ - (4/9)ξu₁' - s₁'/(ε²ξ²) + 2s₁/(ε²ξ³) - s₁''' = 1/(εξ²)
#
# The -s₁''' term is the dispersive axial curvature contribution
# that produces capillary waves in the outer region.

using DifferentialEquations

export solve_outer, OuterSolution

struct OuterSolution
    ξ::Vector{Float64}
    s₁::Vector{Float64}      # perturbation to S
    s₁ξ::Vector{Float64}     # s₁'
    s₁ξξ::Vector{Float64}    # s₁''
    u₁::Vector{Float64}      # perturbation velocity
    ε::Float64
end

# ── Linearised ODE (4-component with S''') ─────────────────────────────
"""
State y = [s₁, s₁', s₁'', u₁].

From mass: u₁' = (2s₁ - 2εu₁ - 2εξ - (4/9)ε²ξ³u₁) / (εξ + (8/9)ε²ξ⁴)
From momentum: s₁''' = -(2/9)u₁ - (4/9)ξu₁' - s₁'/(ε²ξ²) + 2s₁/(ε²ξ³) - 1/(εξ²)

dy = [s₁', s₁'', s₁''', u₁']
"""
function outer_rhs!(dy, y, p, ξ)
    ε = p[1]
    s₁, s₁p, s₁pp, u₁ = y

    if ξ < 1e-6
        dy .= 0.0
        return
    end

    # u₁' from mass (same resolution as before)
    numer_u = 2*s₁ - 2*ε*u₁ - 2*ε*ξ - (4/9)*ε^2*ξ^3*u₁
    denom_u = ε*ξ + (8/9)*ε^2*ξ^4

    if abs(denom_u) < 1e-14
        dy .= 0.0
        return
    end

    u₁p = numer_u / denom_u

    # s₁''' from momentum (with axial curvature term)
    s₁ppp = -(2/9)*u₁ - (4/9)*ξ*u₁p - s₁p/(ε^2*ξ^2) + 2*s₁/(ε^2*ξ^3) - 1/(ε*ξ^2)

    dy[1] = s₁p     # ds₁/dξ
    dy[2] = s₁pp    # ds₁'/dξ = s₁''
    dy[3] = s₁ppp   # ds₁''/dξ = s₁'''
    dy[4] = u₁p     # du₁/dξ
end

# ── Solver ─────────────────────────────────────────────────────────────
"""
    solve_outer(; ε=0.1, ξ_min=3.0, ξ_max=50.0)

Solve the outer linearised problem by integrating inward from ξ_max.
At ξ_max, all perturbations are zero (unperturbed cone) except a small
seed to excite the decaying modes.
"""
function solve_outer(; ε::Float64=0.1, ξ_min::Float64=3.0, ξ_max::Float64=50.0,
                      seed::Float64=0.0)
    # Initial condition at ξ_max: small perturbation
    y0 = [seed, 0.0, 0.0, 0.0]
    tspan = (ξ_max, ξ_min)  # integrate inward
    p = [ε]

    prob = ODEProblem(outer_rhs!, y0, tspan, p)
    sol = solve(prob, Rodas5P(); reltol=1e-8, abstol=1e-10,
                maxiters=500_000)

    ξ_vals = reverse(sol.t)
    s₁_vals = reverse([sol.u[i][1] for i in eachindex(sol.u)])
    s₁p_vals = reverse([sol.u[i][2] for i in eachindex(sol.u)])
    s₁pp_vals = reverse([sol.u[i][3] for i in eachindex(sol.u)])
    u₁_vals = reverse([sol.u[i][4] for i in eachindex(sol.u)])

    OuterSolution(ξ_vals, s₁_vals, s₁p_vals, s₁pp_vals, u₁_vals, ε)
end

function solve_outer_driven(; ε::Float64=0.1, ξ_min::Float64=3.0, ξ_max::Float64=50.0)
    solve_outer(ε=ε, ξ_min=ξ_min, ξ_max=ξ_max, seed=ε^2)
end
