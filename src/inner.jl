# Inner BVP solver: solve the nonlinear similarity ODE near the tip.
#
# The similarity ODE system (correct signs from Bernoulli):
#   2S + 2S'(U - ξ) + S·U' = 0              [mass]
#   -(2/9)U + (4/9)(U - ξ)U' - S'/S² = 0    [momentum]
#
# Solved as an IVP by shooting from the tip ξ₀:
#   BCs at tip: S'(ξ₀) = 0 (smooth rounded tip), S(ξ₀) = S₀ > 0
#   Far-field: S(ξ)/ξ → ε, U(ξ) → 0 as ξ → ∞
#
# Two free parameters (ξ₀, S₀) determined by two far-field conditions.
# The tip position ξ₀ > 0 represents the recoil distance.

using DifferentialEquations
using LinearAlgebra: norm, det

export solve_inner_bvp, InnerSolution

struct InnerSolution
    ξ::Vector{Float64}
    S::Vector{Float64}
    U::Vector{Float64}
    ξ₀::Float64
    S₀::Float64
end

# ── ODE right-hand side ────────────────────────────────────────────────
"""
Resolved similarity ODEs:
  S' = -S²[(2/9)U + (8/9)(U-ξ)] / [1 + (8/9)S(U-ξ)²]
  U' = -2 - 2S'(U-ξ)/S

Denominator 1 + (8/9)S(U-ξ)² > 0 always (no singularity).
"""
function inner_rhs!(du, u, p, ξ)
    S, U = u
    if S ≤ 0
        du[1] = 0.0
        du[2] = 0.0
        return
    end

    v = U - ξ
    numer = -S^2 * (2/9 * U + 8/9 * v)
    denom = 1 + 8/9 * S * v^2

    Sξ = numer / denom
    Uξ = -2 - 2 * Sξ * v / S

    du[1] = Sξ
    du[2] = Uξ
end

# ── Initial conditions at the tip ──────────────────────────────────────
"""
At tip ξ₀: S'(ξ₀) = 0.
Mass → U' = -2. Momentum → U = (4/5)ξ₀.
"""
function tip_initial_conditions(ξ₀::Float64, S₀::Float64)
    U₀ = 4/5 * ξ₀
    (S₀, U₀)
end

# ── Single-shot integration ───────────────────────────────────────────
function _shoot(ξ₀, S₀, ξ_max)
    S₀_val, U₀_val = tip_initial_conditions(ξ₀, S₀)
    prob = ODEProblem(inner_rhs!, [S₀_val, U₀_val], (ξ₀ + 1e-6, ξ_max))
    sol = solve(prob, Tsit5(); reltol=1e-11, abstol=1e-13, maxiters=2_000_000)
    ξ_vals = sol.t
    S_vals = [sol.u[i][1] for i in eachindex(sol.u)]
    U_vals = [sol.u[i][2] for i in eachindex(sol.u)]
    (ξ_vals, S_vals, U_vals)
end

# ── 2D Newton shooting ────────────────────────────────────────────────
"""
    solve_inner_bvp(; ξ₀=2.5, S₀=0.5, ξ_max=100.0, ε=0.1,
                      newton_iters=30, newton_tol=1e-6)

Solve the inner BVP by 2D Newton shooting over (ξ₀, S₀) matching:
  1. S(ξ_max)/ξ_max → ε   (far-field slope)
  2. U(ξ_max) → 0          (far-field velocity)

The tip position ξ₀ > 0 represents how far the tip has recoiled.
"""
function solve_inner_bvp(; ξ₀::Float64=2.5, S₀::Float64=0.5,
                           ξ_max::Float64=100.0, ε::Float64=0.1,
                           newton_iters::Int=30, newton_tol::Float64=1e-6)
    function residual(x)
        ξv, Sv, Uv = _shoot(x[1], x[2], ξ_max)
        [Sv[end] / ξv[end] - ε, Uv[end]]
    end

    x = [ξ₀, S₀]
    δ = 1e-5

    for _ in 1:newton_iters
        F = residual(x)
        norm(F) < newton_tol && break

        # Finite-difference Jacobian
        F_d1 = residual(x .+ [δ, 0])
        F_d2 = residual(x .+ [0, δ])
        J = hcat((F_d1 - F)/δ, (F_d2 - F)/δ)
        abs(det(J)) < 1e-20 && break

        step = J \ F
        x .-= step
        x[1] = max(0.01, x[1])
        x[2] = max(0.01, x[2])
    end

    ξv, Sv, Uv = _shoot(x[1], x[2], ξ_max)
    InnerSolution(ξv, Sv, Uv, x[1], x[2])
end

# ── Far-field diagnostics ──────────────────────────────────────────────
function far_field_slope(sol::InnerSolution)
    n = length(sol.ξ)
    idx = max(1, n-10):n
    [sol.S[i] / sol.ξ[i] for i in idx if sol.ξ[i] > 0]
end
