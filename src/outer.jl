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

export solve_outer, solve_outer_driven, OuterSolution

struct OuterSolution
    ξ::Vector{Float64}
    s₁::Vector{Float64}      # perturbation to S
    s₁ξ::Vector{Float64}     # s₁'
    s₁ξξ::Vector{Float64}    # s₁''
    u₁::Vector{Float64}      # perturbation velocity
    ε::Float64
    diagnostics::NamedTuple
end

OuterSolution(ξ, s₁, s₁ξ, s₁ξξ, u₁, ε) =
    OuterSolution(ξ, s₁, s₁ξ, s₁ξξ, u₁, ε,
                  _manual_solution_diagnostics("OuterSolution constructed directly", ξ))

# ── Linearised ODE (4-component with S''') ─────────────────────────────
"""
State y = [s₁, s₁', s₁'', u₁].

From mass:  2s₁ + 2εu₁ - 2ξs₁' + εξu₁' = 0
  ⟹ u₁' = (2ξs₁' - 2s₁ - 2εu₁) / (εξ)
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

    # u₁' from linearised mass: 2s₁ + 2εu₁ - 2ξs₁' + εξu₁' = 0
    denom_u = ε * ξ

    if abs(denom_u) < 1e-14
        dy .= 0.0
        return
    end

    u₁p = (2*ξ*s₁p - 2*s₁ - 2*ε*u₁) / denom_u

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
                      seed::Float64=0.0, maxiters::Int=500_000)
    # Initial condition at ξ_max: small perturbation
    y0 = [seed, 0.0, 0.0, 0.0]
    tspan = (ξ_max, ξ_min)  # integrate inward
    p = [ε]

    prob = ODEProblem(outer_rhs!, y0, tspan, p)
    sol = solve(prob, Rodas5P(); reltol=1e-8, abstol=1e-10,
                maxiters=maxiters)
    diagnostics = _require_successful_solution(sol, ξ_min; context="outer inward solve")

    ξ_vals = reverse(sol.t)
    s₁_vals = reverse([sol.u[i][1] for i in eachindex(sol.u)])
    s₁p_vals = reverse([sol.u[i][2] for i in eachindex(sol.u)])
    s₁pp_vals = reverse([sol.u[i][3] for i in eachindex(sol.u)])
    u₁_vals = reverse([sol.u[i][4] for i in eachindex(sol.u)])

    OuterSolution(ξ_vals, s₁_vals, s₁p_vals, s₁pp_vals, u₁_vals, ε, diagnostics)
end

function solve_outer_driven(; ε::Float64=0.1, ξ_min::Float64=3.0, ξ_max::Float64=50.0,
                             maxiters::Int=500_000)
    solve_outer(ε=ε, ξ_min=ξ_min, ξ_max=ξ_max, seed=ε^2, maxiters=maxiters)
end

"""
    solve_outer_matched(inner::InnerSolution; ξ_match, ξ_max=100.0)

Solve the outer problem by extracting the perturbation state from the inner
solution at ξ_match and integrating the linearised ODE *outward* from there.

This is the physically correct matching: the inner solution's far-field
oscillatory/decaying behavior seeds the outer solution.
"""
function solve_outer_matched(inner; ξ_match::Float64=15.0, ξ_max::Float64=100.0,
                             maxiters::Int=500_000)
    ε = inner.S[end] / inner.ξ[end]  # estimate ε from far-field slope

    # Interpolate inner state at matching point
    s₁_m = _interp_val(inner.ξ, inner.S, ξ_match) - ε * ξ_match
    s₁p_m = _interp_val(inner.ξ, inner.Sξ, ξ_match) - ε
    s₁pp_m = _interp_val(inner.ξ, inner.Sξξ, ξ_match)
    u₁_m = _interp_val(inner.ξ, inner.U, ξ_match)

    y0 = [s₁_m, s₁p_m, s₁pp_m, u₁_m]
    tspan = (ξ_match, ξ_max)
    p = [ε]

    prob = ODEProblem(outer_rhs!, y0, tspan, p)
    sol = solve(prob, Rodas5P(); reltol=1e-8, abstol=1e-10, maxiters=maxiters)
    diagnostics = _require_successful_solution(sol, ξ_max; context="outer matched solve")

    ξ_vals = sol.t
    s₁_vals = [sol.u[i][1] for i in eachindex(sol.u)]
    s₁p_vals = [sol.u[i][2] for i in eachindex(sol.u)]
    s₁pp_vals = [sol.u[i][3] for i in eachindex(sol.u)]
    u₁_vals = [sol.u[i][4] for i in eachindex(sol.u)]

    OuterSolution(ξ_vals, s₁_vals, s₁p_vals, s₁pp_vals, u₁_vals, ε, diagnostics)
end

# Simple linear interpolation for a single point
function _interp_val(xs, ys, x)
    if x <= xs[1]; return ys[1]; end
    if x >= xs[end]; return ys[end]; end
    i = searchsortedlast(xs, x)
    i = clamp(i, 1, length(xs)-1)
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    ys[i] + t * (ys[i+1] - ys[i])
end
