# Inner BVP solver: solve the nonlinear similarity ODE near the tip.
#
# With axial curvature, the similarity ODE system is:
#   2S + 2S'(U - ξ) + S·U' = 0                          [mass]
#   -(2/9)U + (4/9)(U - ξ)U' - S'/S² - S''' = 0         [momentum]
#
# State vector: [S, S', S'', U] — 4 components (3rd-order in S, 1st in U).
# The mass equation gives U' algebraically; momentum gives S''' algebraically.
#
# Shooting parameters: (ξ₀, S₀, S''₀) matched to 3 far-field conditions.

using DifferentialEquations
using LinearAlgebra: norm, det

export solve_inner_bvp, InnerSolution

struct InnerSolution
    ξ::Vector{Float64}
    S::Vector{Float64}
    Sξ::Vector{Float64}     # S'
    Sξξ::Vector{Float64}    # S''
    U::Vector{Float64}
    ξ₀::Float64
    S₀::Float64
    Sξξ₀::Float64           # S''(ξ₀)
end

# ── ODE right-hand side (4-component) ──────────────────────────────────
"""
State: u = [S, Sp, Spp, U] where Sp=S', Spp=S''.

From mass:  U' = -2 - 2·S'·(U-ξ)/S
From momentum: S''' = -(2/9)U + (4/9)(U-ξ)U' - S'/S²
"""
function inner_rhs!(du, u, p, ξ)
    S, Sp, Spp, U = u
    if S ≤ 0
        du .= 0.0
        return
    end

    v = U - ξ

    # U' from mass (algebraic)
    Up = -2 - 2 * Sp * v / S

    # S''' from momentum (algebraic, once U' is known)
    Sppp = -(2/9) * U + (4/9) * v * Up - Sp / S^2

    du[1] = Sp      # dS/dξ
    du[2] = Spp     # dS'/dξ = S''
    du[3] = Sppp    # dS''/dξ = S'''
    du[4] = Up      # dU/dξ
end

# ── Tip initial conditions ─────────────────────────────────────────────
"""
At tip ξ₀: S'(ξ₀) = 0 (smooth rounded tip).
From mass with S'=0:  U' = -2.
From momentum with S'=0:  S''' = -(2/9)U + (4/9)(U-ξ₀)(-2)
  = -(2/9)U - (8/9)(U-ξ₀) = -(10/9)U + (8/9)ξ₀
With U = (4/5)ξ₀:  S''' = -(10/9)(4/5)ξ₀ + (8/9)ξ₀ = -(8/9)ξ₀ + (8/9)ξ₀ = 0.

So at the tip: S=S₀, S'=0, S''=Sξξ₀ (free), U=(4/5)ξ₀, S'''=0.
"""
function tip_initial_conditions(ξ₀::Float64, S₀::Float64, Sξξ₀::Float64)
    U₀ = 4/5 * ξ₀
    [S₀, 0.0, Sξξ₀, U₀]
end

# ── Single-shot integration ───────────────────────────────────────────
function _shoot(ξ₀, S₀, Sξξ₀, ξ_max)
    u0 = tip_initial_conditions(ξ₀, S₀, Sξξ₀)
    prob = ODEProblem(inner_rhs!, u0, (ξ₀ + 1e-6, ξ_max))
    # Use stiff solver — S''' term produces rapid oscillations
    sol = solve(prob, Rodas5P(); reltol=1e-8, abstol=1e-10, maxiters=2_000_000)
    ξ_vals = sol.t
    S_vals = [sol.u[i][1] for i in eachindex(sol.u)]
    Sp_vals = [sol.u[i][2] for i in eachindex(sol.u)]
    Spp_vals = [sol.u[i][3] for i in eachindex(sol.u)]
    U_vals = [sol.u[i][4] for i in eachindex(sol.u)]
    (ξ_vals, S_vals, Sp_vals, Spp_vals, U_vals)
end

# ── 3D Newton shooting ────────────────────────────────────────────────
"""
    solve_inner_bvp(; ξ₀=2.5, S₀=0.5, Sξξ₀=-1.0, ξ_max=100.0, ε=0.1,
                      newton_iters=40, newton_tol=1e-6)

Solve the inner BVP by 3D Newton shooting over (ξ₀, S₀, S''₀) matching:
  1. S(ξ_max)/ξ_max → ε   (far-field slope)
  2. U(ξ_max) → 0          (far-field velocity)
  3. S''(ξ_max) → 0         (curvature perturbation decays)
"""
function solve_inner_bvp(; ξ₀::Float64=2.79, S₀::Float64=0.28, Sξξ₀::Float64=0.57,
                           ξ_max::Float64=60.0, ε::Float64=0.1,
                           newton_iters::Int=30, newton_tol::Float64=1e-4)
    function residual(x)
        ξv, Sv, _, Sppv, Uv = _shoot(x[1], x[2], x[3], ξ_max)
        [Sv[end] / ξv[end] - ε, Uv[end], Sppv[end]]
    end

    x = [ξ₀, S₀, Sξξ₀]
    δ = 1e-5

    for _ in 1:newton_iters
        F = residual(x)
        norm(F) < newton_tol && break

        # 3x3 finite-difference Jacobian
        J = zeros(3, 3)
        for j in 1:3
            xp = copy(x)
            xp[j] += δ
            J[:, j] = (residual(xp) - F) / δ
        end
        abs(det(J)) < 1e-20 && break

        step = J \ F

        # Damped Newton: reduce step until residual decreases
        α = 1.0
        x_trial = x .- α .* step
        for _ in 1:10
            x_trial[1] = max(0.1, x_trial[1])
            x_trial[2] = max(0.01, x_trial[2])
            F_trial = residual(x_trial)
            norm(F_trial) < norm(F) && break
            α *= 0.5
            x_trial = x .- α .* step
        end
        x .= x_trial
        x[1] = max(0.1, x[1])
        x[2] = max(0.01, x[2])
    end

    ξv, Sv, Spv, Sppv, Uv = _shoot(x[1], x[2], x[3], ξ_max)
    InnerSolution(ξv, Sv, Spv, Sppv, Uv, x[1], x[2], x[3])
end

# ── Far-field diagnostics ──────────────────────────────────────────────
function far_field_slope(sol::InnerSolution)
    n = length(sol.ξ)
    idx = max(1, n-10):n
    [sol.S[i] / sol.ξ[i] for i in idx if sol.ξ[i] > 0]
end
